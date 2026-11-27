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

// -------------------------------------------------------------------------- LevelAncestor

LevelAncestor::LevelAncestor(const OrderedTree& tree) {
  build(tree, paper_block_size(tree.size()));
}

LevelAncestor::LevelAncestor(const OrderedTree& tree, std::uint32_t block_size) {
  if (block_size < 2 || block_size > 8) {
    throw std::invalid_argument("LevelAncestor: block size must lie in [2, 8]");
  }
  build(tree, block_size);
}

void LevelAncestor::build(const OrderedTree& tree, std::uint32_t block_size) {
  const std::size_t k = tree.size();
  const std::uint32_t bsz = block_size;
  block_ = bsz;
  const auto pre = tree.preorder();
  const Vertex root = tree.root();
  depth_.resize(k);
  for (Vertex v = 0; v < k; ++v) depth_[v] = tree.depth(v);
  lg_ = detail::FloorLogTable(2 * k);

  // (1) Marks. Each vertex returns to its parent the size of its residual subtree (itself plus
  // its children's residuals); if the sum reaches B the vertex is marked and returns zero. The
  // augmented root is marked as well. Residuals only need children before parents, so reverse
  // preorder gives the same marks as postorder.
  std::vector<std::uint32_t> residual(k, 0);
  std::vector<std::uint8_t> marked(k, 0);
  for (std::size_t idx = k; idx-- > 0;) {
    const Vertex v = pre[idx];
    std::uint32_t r = 1 + residual[v];
    if (r >= bsz || v == root) {
      marked[v] = 1;
      r = 0;
    }
    if (v != root) residual[tree.parent(v)] += r;
  }
  mark_id_.assign(k, kUnmarked);
  std::uint32_t marks = 0;
  for (const Vertex v : pre) {
    if (marked[v] != 0) mark_id_[v] = marks++;
  }

  // (2) Micro trees: components of unmarked vertices, each with fewer than B vertices (its top
  // returned its residual size, which is < B). Global preorder restricted to a component is the
  // component's own preorder, which gives the local positions.
  micro_.assign(k, kUnmarked);
  local_pos_.assign(k, 0);
  micro_top_depth_.clear();
  micro_top_parent_.clear();
  std::vector<std::uint32_t> micro_size;
  for (const Vertex v : pre) {
    if (marked[v] != 0) continue;
    const Vertex p = tree.parent(v);  // exists: the root is marked
    std::uint32_t m = 0;
    if (marked[p] != 0) {
      m = static_cast<std::uint32_t>(micro_size.size());
      micro_size.push_back(0);
      micro_top_depth_.push_back(depth_[v]);
      micro_top_parent_.push_back(p);
    } else {
      m = micro_[p];
    }
    micro_[v] = m;
    local_pos_[v] = static_cast<std::uint8_t>(micro_size[m]++);
  }
  const std::size_t num_micro = micro_size.size();
  micro_begin_.assign(num_micro + 1, 0);
  max_micro_size_ = 0;
  for (std::size_t m = 0; m < num_micro; ++m) {
    micro_begin_[m + 1] = micro_begin_[m] + micro_size[m];
    if (micro_size[m] > max_micro_size_) max_micro_size_ = micro_size[m];
  }
  micro_vertex_.assign(micro_begin_[num_micro], kNoVertex);
  for (Vertex v = 0; v < k; ++v) {
    if (marked[v] == 0) micro_vertex_[micro_begin_[micro_[v]] + local_pos_[v]] = v;
  }
  // Shape codes: a leading one, then one entering bit (1) and one leaving bit (0) per vertex
  // during the depth-first traversal of the micro tree.
  micro_shape_.assign(num_micro, 1);
  {
    std::vector<std::pair<Vertex, std::uint32_t>> stack;  // (vertex, next child index)
    stack.emplace_back(root, 0);
    while (!stack.empty()) {
      const Vertex v = stack.back().first;
      const std::uint32_t next = stack.back().second;
      const auto ch = tree.children(v);
      if (next < ch.size()) {
        ++stack.back().second;
        const Vertex c = ch[next];
        if (marked[c] == 0) micro_shape_[micro_[c]] = (micro_shape_[micro_[c]] << 1) | 1u;
        stack.emplace_back(c, 0);
      } else {
        if (marked[v] == 0) micro_shape_[micro_[v]] <<= 1;
        stack.pop_back();
      }
    }
  }

  // (3) Shape tables: for every valid traversal string of s < B vertices, the local ancestor of
  // every preorder position at every allowable relative depth, by tracing parent chains.
  const std::size_t stride = static_cast<std::size_t>(bsz) * bsz;
  shape_table_.assign((std::size_t{1} << (2 * bsz - 1)) * stride, 0);
  for (std::uint32_t s = 1; s < bsz; ++s) {
    const std::uint32_t len = 2 * s;
    for (std::uint32_t bits = 0; bits < (1u << len); ++bits) {
      std::uint8_t local_parent[8] = {};
      std::uint8_t local_depth[8] = {};
      std::uint8_t stack[8] = {};
      std::uint32_t sp = 0;
      std::uint32_t count = 0;
      bool valid = true;
      for (std::uint32_t t = 0; t < len && valid; ++t) {
        if (((bits >> (len - 1 - t)) & 1u) != 0) {
          if ((t > 0 && sp == 0) || count == s) {
            valid = false;  // a second root, or too many vertices
            break;
          }
          local_parent[count] = sp > 0 ? stack[sp - 1] : 0;
          local_depth[count] = static_cast<std::uint8_t>(sp);
          stack[sp++] = static_cast<std::uint8_t>(count++);
        } else if (sp == 0) {
          valid = false;  // leaving more vertices than entered
        } else {
          --sp;
        }
      }
      if (!valid || sp != 0 || count != s) continue;
      const std::size_t code = (std::size_t{1} << len) | bits;
      for (std::uint32_t pos = 0; pos < s; ++pos) {
        for (std::uint32_t rd = 0; rd <= local_depth[pos]; ++rd) {
          std::uint32_t a = pos;
          for (std::uint32_t up = local_depth[pos]; up > rd; --up) a = local_parent[a];
          shape_table_[code * stride + static_cast<std::size_t>(pos) * bsz + rd] =
              static_cast<std::uint8_t>(a);
        }
      }
    }
  }

  // (4) Jump pointers: at each mark p, its ancestor 2^i edges above whenever it exists, read
  // from a stack indexed by depth that holds the current root path during the traversal
  // (in preorder, the last vertex seen at a smaller depth is the ancestor at that depth).
  jump_begin_.assign(1, 0);
  jump_.clear();
  {
    std::vector<Vertex> path(k, kNoVertex);
    for (const Vertex v : pre) {
      path[depth_[v]] = v;
      if (marked[v] == 0) continue;
      for (std::uint32_t step = 1; step <= depth_[v]; step <<= 1) {
        jump_.push_back(path[depth_[v] - step]);
      }
      jump_begin_.push_back(static_cast<std::uint32_t>(jump_.size()));
    }
  }

  // (5) Ladders. Heights (in edges) in postorder order; at every nonleaf select a child of
  // maximum height (first in child order on ties). The selected edges partition the vertices
  // into paths; each path is stored in depth order and extended upward by its own number of
  // vertices or as far as the root, whichever is shorter.
  std::vector<std::uint32_t> height(k, 0);
  std::vector<Vertex> selected(k, kNoVertex);
  for (std::size_t idx = k; idx-- > 0;) {
    const Vertex v = pre[idx];
    for (const Vertex c : tree.children(v)) {
      if (selected[v] == kNoVertex || height[c] + 1 > height[v]) {
        height[v] = height[c] + 1;
        selected[v] = c;
      }
    }
  }
  ladder_.clear();
  ladder_.reserve(2 * k);
  ladder_pos_.assign(k, 0);
  std::vector<Vertex> path;
  for (const Vertex head : pre) {
    if (head != root && selected[tree.parent(head)] == head) continue;  // not a path head
    path.clear();
    for (Vertex x = head; x != kNoVertex; x = selected[x]) path.push_back(x);
    const std::size_t ext = path.size() < depth_[head] ? path.size() : depth_[head];
    const std::size_t base = ladder_.size();
    ladder_.resize(base + ext + path.size());
    Vertex a = head;
    for (std::size_t t = ext; t-- > 0;) {
      a = tree.parent(a);
      ladder_[base + t] = a;
    }
    for (std::size_t i = 0; i < path.size(); ++i) {
      ladder_[base + ext + i] = path[i];
      ladder_pos_[path[i]] = static_cast<std::uint32_t>(base + ext + i);
    }
  }
}

Vertex LevelAncestor::query(Vertex v, std::uint32_t j) const {
  if (j > depth_[v]) throw std::out_of_range("LevelAncestor: depth exceeds the vertex's depth");
  Vertex p = v;
  if (mark_id_[v] == kUnmarked) {
    // Unmarked v and the target inside its micro tree: one shape-table lookup.
    const std::uint32_t m = micro_[v];
    const std::uint32_t top_depth = micro_top_depth_[m];
    if (j >= top_depth) {
      const std::size_t cell =
          (static_cast<std::size_t>(micro_shape_[m]) * block_ + local_pos_[v]) * block_ +
          (j - top_depth);
      return micro_vertex_[micro_begin_[m] + shape_table_[cell]];
    }
    p = micro_top_parent_[m];  // marked, and on the query path
  }
  const std::uint32_t d = depth_[p] - j;
  if (d == 0) return p;
  // Jump up h = 2^floor(log2 d) edges, then d - h < h positions up x's extended path.
  const std::uint32_t i = lg_(d);
  const std::uint32_t h = 1u << i;
  const Vertex x = jump_[jump_begin_[mark_id_[p]] + i];
  return ladder_[ladder_pos_[x] - (d - h)];
}

// -------------------------------------------------------------------------- ForestQueries

ForestQueries::ForestQueries(std::span<const Vertex> parent) {
  build(parent, 0);
}

ForestQueries::ForestQueries(std::span<const Vertex> parent, std::uint32_t block_size) {
  if (block_size < 2 || block_size > 8) {
    throw std::invalid_argument("ForestQueries: block size must lie in [2, 8]");
  }
  build(parent, block_size);
}

void ForestQueries::build(std::span<const Vertex> parent, std::uint32_t block_size) {
  // Lemma 3.9 / Sec 3.4: augment the forest by a common root R = n.
  const std::size_t n = parent.size();
  if (n + 1 >= static_cast<std::size_t>(kNoVertex)) {
    throw std::length_error("ForestQueries: too many vertices");
  }
  const auto aug_root = static_cast<Vertex>(n);
  std::vector<Vertex> aug(n + 1, kNoVertex);
  for (std::size_t v = 0; v < n; ++v) {
    const Vertex p = parent[v];
    if (p != kNoVertex && p >= n) throw std::out_of_range("ForestQueries: parent out of range");
    aug[v] = p == kNoVertex ? aug_root : p;
  }
  const OrderedTree tree(aug);  // rejects cycles
  // Record the original component (the original root) of every vertex.
  component_.assign(n, kNoVertex);
  for (const Vertex v : tree.preorder()) {
    if (v == aug_root) continue;
    component_[v] = aug[v] == aug_root ? v : component_[aug[v]];
  }
  if (block_size == 0) {
    lca_ = LcaTable(tree);
    la_ = LevelAncestor(tree);
  } else {
    lca_ = LcaTable(tree, block_size);
    la_ = LevelAncestor(tree, block_size);
  }
}

std::optional<Vertex> ForestQueries::lca(Vertex u, Vertex v) const {
  if (component_[u] != component_[v]) return std::nullopt;  // different original trees
  return lca_.lca(u, v);
}

Vertex ForestQueries::level_ancestor(Vertex v, std::uint32_t j) const {
  if (j > depth(v)) throw std::out_of_range("ForestQueries: depth exceeds the vertex's depth");
  return la_.query(v, j + 1);  // augmented depth = original depth + 1
}

// ----------------------------------------------------------------------- BinaryLiftingRef

BinaryLiftingRef::BinaryLiftingRef(std::span<const Vertex> parent)
    : depth_(parent.size(), 0), root_(parent.size(), kNoVertex) {
  const std::size_t n = parent.size();
  for (std::size_t v = 0; v < n; ++v) {
    if (parent[v] != kNoVertex && parent[v] >= n) {
      throw std::out_of_range("BinaryLiftingRef: parent out of range");
    }
  }
  // Depths and roots by walking up to a resolved vertex (0 unseen, 1 on this walk, 2 resolved).
  std::vector<std::uint8_t> state(n, 0);
  std::vector<Vertex> walk;
  std::uint32_t max_depth = 0;
  for (std::size_t s = 0; s < n; ++s) {
    auto x = static_cast<Vertex>(s);
    walk.clear();
    while (state[x] == 0) {
      state[x] = 1;
      walk.push_back(x);
      if (parent[x] == kNoVertex) break;
      x = parent[x];
    }
    if (state[x] == 1 && parent[x] != kNoVertex) {
      throw std::invalid_argument("BinaryLiftingRef: parent array has a cycle");
    }
    for (std::size_t i = walk.size(); i-- > 0;) {
      const Vertex w = walk[i];
      const Vertex p = parent[w];
      depth_[w] = p == kNoVertex ? 0 : depth_[p] + 1;
      root_[w] = p == kNoVertex ? w : root_[p];
      state[w] = 2;
      if (depth_[w] > max_depth) max_depth = depth_[w];
    }
  }
  up_.emplace_back(parent.begin(), parent.end());
  for (std::uint64_t step = 2; step <= max_depth; step <<= 1) {
    const auto& prev = up_.back();
    std::vector<Vertex> next(n, kNoVertex);
    for (std::size_t v = 0; v < n; ++v) {
      if (prev[v] != kNoVertex) next[v] = prev[prev[v]];
    }
    up_.push_back(std::move(next));
  }
}

Vertex BinaryLiftingRef::level_ancestor(Vertex v, std::uint32_t j) const {
  if (j > depth_[v]) throw std::out_of_range("BinaryLiftingRef: depth exceeds the vertex's depth");
  std::uint32_t d = depth_[v] - j;
  for (std::size_t i = 0; d != 0; ++i, d >>= 1) {
    if ((d & 1u) != 0) v = up_[i][v];
  }
  return v;
}

std::optional<Vertex> BinaryLiftingRef::lca(Vertex u, Vertex v) const {
  if (root_[u] != root_[v]) return std::nullopt;
  if (depth_[u] < depth_[v]) std::swap(u, v);
  u = level_ancestor(u, depth_[v]);
  if (u == v) return u;
  for (std::size_t i = up_.size(); i-- > 0;) {
    if (up_[i][u] != up_[i][v]) {
      u = up_[i][u];
      v = up_[i][v];
    }
  }
  return up_[0][u];
}

}  // namespace sssp::paper
