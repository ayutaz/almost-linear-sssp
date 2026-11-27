// Sec 3.4, Lemmas 3.7-3.9. See forest_queries.hpp.
#include "paper/primitives/forest_queries.hpp"

#include <limits>
#include <stdexcept>
#include <utility>

namespace sssp::paper {

std::uint32_t paper_block_size(std::size_t k) {
  // floor(log2(k+1)) by repeated halving (Sec 2 word operations), then the shift by 3 is /8.
  std::uint64_t x = static_cast<std::uint64_t>(k) + 1;
  std::uint32_t lg = 0;
  while (x > 1) {
    x >>= 1;
    ++lg;
  }
  const std::uint32_t b = lg >> 3;
  return b < 2 ? 2 : b;
}

namespace detail {

FloorLogTable::FloorLogTable(std::size_t n) : lg_(n + 1, 0) {
  for (std::size_t i = 2; i <= n; ++i) lg_[i] = static_cast<std::uint8_t>(lg_[i >> 1] + 1);
}

}  // namespace detail

// ------------------------------------------------------------------------------ OrderedTree

OrderedTree::OrderedTree(std::span<const Vertex> parent) : parent_(parent.begin(), parent.end()) {
  const std::size_t n = parent_.size();
  if (n >= static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max() / 2)) {
    throw std::length_error("OrderedTree: too many vertices");
  }
  std::size_t roots = 0;
  child_begin_.assign(n + 2, 0);
  for (std::size_t v = 0; v < n; ++v) {
    const Vertex p = parent_[v];
    if (p == kNoVertex) {
      ++roots;
      root_ = static_cast<Vertex>(v);
    } else if (p >= n) {
      throw std::out_of_range("OrderedTree: parent index out of range");
    } else {
      ++child_begin_[p + 2];
    }
  }
  if (roots != 1) throw std::invalid_argument("OrderedTree: parent array needs exactly one root");
  // Counting placement: children of each vertex in increasing index order.
  for (std::size_t v = 0; v < n; ++v) child_begin_[v + 2] += child_begin_[v + 1];
  child_.resize(n - 1);
  for (std::size_t v = 0; v < n; ++v) {
    const Vertex p = parent_[v];
    if (p != kNoVertex) child_[child_begin_[p + 1]++] = static_cast<Vertex>(v);
  }
  child_begin_.pop_back();  // now child_begin_[v] .. child_begin_[v+1] are v's children

  // Preorder from the root with an explicit stack; a vertex not reached lies on a cycle.
  depth_.assign(n, 0);
  preorder_.reserve(n);
  std::vector<Vertex> stack{root_};
  while (!stack.empty()) {
    const Vertex v = stack.back();
    stack.pop_back();
    preorder_.push_back(v);
    const auto ch = children(v);
    for (std::size_t i = ch.size(); i-- > 0;) {
      depth_[ch[i]] = depth_[v] + 1;
      stack.push_back(ch[i]);
    }
  }
  if (preorder_.size() != n) throw std::invalid_argument("OrderedTree: parent array has a cycle");
}

// ------------------------------------------------------------------------------- LcaTable

LcaTable::LcaTable(const OrderedTree& tree) {
  build(tree, paper_block_size(tree.size()));
}

LcaTable::LcaTable(const OrderedTree& tree, std::uint32_t block_size) {
  if (block_size < 2 || block_size > 8) {
    throw std::invalid_argument("LcaTable: block size must lie in [2, 8]");
  }
  build(tree, block_size);
}

void LcaTable::build(const OrderedTree& tree, std::uint32_t block_size) {
  const std::size_t k = tree.size();
  block_ = block_size;
  const std::uint32_t bsz = block_size;

  // Depth-first tour from the root recording the reached endpoint after each edge traversal
  // (Sec 3.4); 2k - 1 entries. Block number and offset are advanced during the scan.
  tour_vertex_.clear();
  tour_depth_.clear();
  block_start_.clear();
  block_code_.clear();
  tour_vertex_.reserve(2 * k - 1);
  tour_depth_.reserve(2 * k - 1);
  first_pos_.assign(k, 0);
  first_block_.assign(k, 0);
  first_off_.assign(k, 0);
  std::uint32_t pos = 0;
  std::uint32_t blk = 0;
  std::uint32_t off = 0;
  std::uint32_t next_start = 0;
  auto emit = [&](Vertex v) {
    const std::uint32_t d = tree.depth(v);
    if (off == 0) {
      block_start_.push_back(next_start);
      block_code_.push_back(0);
      next_start += bsz;
    } else if (d > tour_depth_.back()) {
      block_code_[blk] |= 1u << (off - 1);  // sign of the change into offset `off`
    }
    tour_vertex_.push_back(v);
    tour_depth_.push_back(d);
    ++pos;
    if (++off == bsz) {
      off = 0;
      ++blk;
    }
  };
  auto first_visit = [&](Vertex v) {
    first_pos_[v] = pos;
    first_block_[v] = blk;
    first_off_[v] = static_cast<std::uint8_t>(off);
    emit(v);
  };
  {
    std::vector<std::pair<Vertex, std::uint32_t>> stack;  // (vertex, next child index)
    first_visit(tree.root());
    stack.emplace_back(tree.root(), 0);
    while (!stack.empty()) {
      auto& [v, next] = stack.back();
      const auto ch = tree.children(v);
      if (next < ch.size()) {
        const Vertex c = ch[next++];
        first_visit(c);
        stack.emplace_back(c, 0);
      } else {
        stack.pop_back();
        if (!stack.empty()) emit(stack.back().first);
      }
    }
  }
  const auto num_blocks = static_cast<std::uint32_t>(block_start_.size());
  // A short final block keeps zero signs for its missing changes ("arbitrary signs").
  const std::uint32_t last_len = off == 0 ? bsz : off;

  // Sign tables: for every code of B-1 signs and every subinterval [i, j], the first offset of
  // minimum relative depth. Valid for every block since each tour step changes depth by one.
  const std::uint32_t codes = 1u << (bsz - 1);
  in_table_.assign(static_cast<std::size_t>(codes) * bsz * bsz, 0);
  for (std::uint32_t code = 0; code < codes; ++code) {
    for (std::uint32_t i = 0; i < bsz; ++i) {
      std::int32_t rel = 0;
      std::int32_t best = 0;
      std::uint32_t best_at = i;
      const std::size_t row = (static_cast<std::size_t>(code) * bsz + i) * bsz;
      in_table_[row + i] = static_cast<std::uint8_t>(i);
      for (std::uint32_t j = i + 1; j < bsz; ++j) {
        rel += ((code >> (j - 1)) & 1u) != 0 ? 1 : -1;
        if (rel < best) {
          best = rel;
          best_at = j;
        }
        in_table_[row + j] = static_cast<std::uint8_t>(best_at);
      }
    }
  }

  // Whole-block table: level 0 is each block's first minimum; level i combines the two
  // half-length answers of level i-1 (Sec 3.4). Floor logs through 2k (Sec 2).
  lg_ = detail::FloorLogTable(2 * k);
  sparse_.clear();
  sparse_level_.clear();
  sparse_level_.push_back(0);
  for (std::uint32_t b = 0; b < num_blocks; ++b) {
    sparse_.push_back(in_block(b, 0, b + 1 == num_blocks ? last_len - 1 : bsz - 1));
  }
  for (std::uint32_t half = 1; 2 * half <= num_blocks; half <<= 1) {
    const std::size_t prev = sparse_level_.back();
    sparse_level_.push_back(sparse_.size());
    for (std::uint32_t b = 0; b + 2 * half <= num_blocks; ++b) {
      sparse_.push_back(first_min(sparse_[prev + b], sparse_[prev + b + half]));
    }
  }
}

Vertex LcaTable::lca(Vertex u, Vertex v) const {
  // Sec 3.4 query: the first minimum-depth entry between the two first occurrences.
  if (first_pos_[u] > first_pos_[v]) std::swap(u, v);
  const std::uint32_t bl = first_block_[u];
  const std::uint32_t br = first_block_[v];
  const std::uint32_t ol = first_off_[u];
  const std::uint32_t orr = first_off_[v];
  if (bl == br) return tour_vertex_[in_block(bl, ol, orr)];
  std::uint32_t best = in_block(bl, ol, block_ - 1);  // partial left block
  const std::uint32_t h = br - bl - 1;                // whole blocks in the middle
  if (h > 0) {
    // First and last intervals of 2^floor(log2 h) whole blocks cover all h blocks.
    const std::uint32_t i = lg_(h);
    const std::size_t level = sparse_level_[i];
    best = first_min(best, sparse_[level + bl + 1]);
    best = first_min(best, sparse_[level + br - (1u << i)]);
  }
  best = first_min(best, in_block(br, 0, orr));  // partial right block
  return tour_vertex_[best];
}

}  // namespace sssp::paper
