#pragma once
// Sec 3.4 "Lowest Common Ancestors and Level Ancestors", Lemmas 3.7-3.9.
//
// Trees are given as parent arrays (kNoVertex marks a root). The child order the paper fixes
// ("fix an order of the children", Sec 3.4) is increasing vertex index. Only integer word
// operations are used (Sec 2): lookups, shifts, bit operations, additions and comparisons.
// Positions inside fixed-size blocks are tracked during scans, so no integer division is used.

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "sssp/graph.hpp"

namespace sssp::paper {

/// Sec 3.4: B = max(2, floor(log2(k+1)/8)) for k (augmented) vertices. Note that
/// floor(log2(k+1)/8) = floor(floor(log2(k+1))/8), so integer shifts suffice.
std::uint32_t paper_block_size(std::size_t k);

namespace detail {

/// Sec 2 / Sec 3.4: floor logarithms of 1..n stored in an array, lg[i] = lg[i >> 1] + 1.
class FloorLogTable {
 public:
  FloorLogTable() = default;
  explicit FloorLogTable(std::size_t n);
  [[nodiscard]] std::uint32_t operator()(std::size_t i) const { return lg_[i]; }

 private:
  std::vector<std::uint8_t> lg_;
};

}  // namespace detail

/// A rooted tree with a fixed child order (increasing vertex index), built from a parent array
/// with exactly one root. Throws std::out_of_range for a parent index >= n and
/// std::invalid_argument if the array does not describe exactly one rooted tree.
class OrderedTree {
 public:
  explicit OrderedTree(std::span<const Vertex> parent);

  [[nodiscard]] std::size_t size() const noexcept { return parent_.size(); }
  [[nodiscard]] Vertex root() const noexcept { return root_; }
  [[nodiscard]] Vertex parent(Vertex v) const { return parent_[v]; }
  [[nodiscard]] std::uint32_t depth(Vertex v) const { return depth_[v]; }
  [[nodiscard]] std::span<const Vertex> children(Vertex v) const {
    return {child_.data() + child_begin_[v], child_.data() + child_begin_[v + 1]};
  }
  /// Vertices in preorder (children visited in the fixed order).
  [[nodiscard]] std::span<const Vertex> preorder() const noexcept { return preorder_; }

 private:
  std::vector<Vertex> parent_;
  std::vector<std::uint32_t> child_begin_;
  std::vector<Vertex> child_;
  std::vector<std::uint32_t> depth_;
  std::vector<Vertex> preorder_;
  Vertex root_ = kNoVertex;
};

/// Lemma 3.7: Euler-tour LCA with +-1 block tables. O(k) preprocessing and space (for the
/// paper's B), O(1) query.
class LcaTable {
 public:
  LcaTable() = default;
  /// Uses the paper's block size B = paper_block_size(tree.size()).
  explicit LcaTable(const OrderedTree& tree);
  /// Test hook: explicit block size in [2, 8], so that the block tables are exercised on small
  /// trees (the paper's formula gives B = 2 for every k < 2^24 - 1). Mode A callers use the
  /// constructor above.
  LcaTable(const OrderedTree& tree, std::uint32_t block_size);

  [[nodiscard]] Vertex lca(Vertex u, Vertex v) const;

  [[nodiscard]] std::uint32_t block_size() const noexcept { return block_; }
  [[nodiscard]] std::size_t tour_length() const noexcept { return tour_vertex_.size(); }
  [[nodiscard]] std::size_t num_blocks() const noexcept { return block_start_.size(); }
  /// Entries of the whole-block sparse table, O((1 + k/B) log(k+1)) (proof of Lemma 3.7).
  [[nodiscard]] std::size_t sparse_entries() const noexcept { return sparse_.size(); }

 private:
  void build(const OrderedTree& tree, std::uint32_t block_size);
  [[nodiscard]] std::uint32_t in_block(std::uint32_t block, std::uint32_t i,
                                       std::uint32_t j) const {
    return block_start_[block] + in_table_[((block_code_[block] * block_) + i) * block_ + j];
  }
  [[nodiscard]] std::uint32_t first_min(std::uint32_t a, std::uint32_t b) const {
    return tour_depth_[b] < tour_depth_[a] ? b : a;  // a precedes b: ties keep the first
  }

  std::uint32_t block_ = 2;
  std::vector<Vertex> tour_vertex_;         ///< Euler tour, 2k - 1 entries
  std::vector<std::uint32_t> tour_depth_;   ///< depth of each tour entry
  std::vector<std::uint32_t> first_pos_;    ///< first occurrence of each vertex
  std::vector<std::uint32_t> first_block_;  ///< block of the first occurrence
  std::vector<std::uint8_t> first_off_;     ///< offset of the first occurrence in its block
  std::vector<std::uint32_t> block_start_;  ///< first tour position of each block
  std::vector<std::uint32_t> block_code_;   ///< B-1 depth-change signs (1 = step down the tree)
  std::vector<std::uint8_t> in_table_;      ///< [code][i][j]: first minimum offset in [i, j]
  std::vector<std::uint32_t> sparse_;       ///< level i: first minimum of 2^i whole blocks
  std::vector<std::size_t> sparse_level_;   ///< start of each level in sparse_
  detail::FloorLogTable lg_;
};

/// Lemma 3.8: level ancestors by macro-micro marking, shape tables for micro trees, jump
/// pointers at marks and ladders. O(k) preprocessing and space (for the paper's B), O(1) query.
class LevelAncestor {
 public:
  LevelAncestor() = default;
  /// Uses the paper's block size B = paper_block_size(tree.size()).
  explicit LevelAncestor(const OrderedTree& tree);
  /// Test hook: explicit block size in [2, 8] (see LcaTable). Mode A callers use the paper's B.
  LevelAncestor(const OrderedTree& tree, std::uint32_t block_size);

  /// The ancestor of v at depth j. Throws std::out_of_range unless 0 <= j <= depth(v).
  [[nodiscard]] Vertex query(Vertex v, std::uint32_t j) const;
  [[nodiscard]] std::uint32_t depth(Vertex v) const { return depth_[v]; }

  [[nodiscard]] std::uint32_t block_size() const noexcept { return block_; }
  [[nodiscard]] std::size_t num_marks() const noexcept { return jump_begin_.size() - 1; }
  [[nodiscard]] std::size_t num_micro_trees() const noexcept { return micro_begin_.size() - 1; }
  [[nodiscard]] std::size_t max_micro_tree_size() const noexcept { return max_micro_size_; }
  [[nodiscard]] std::size_t jump_entries() const noexcept { return jump_.size(); }
  [[nodiscard]] std::size_t ladder_entries() const noexcept { return ladder_.size(); }

 private:
  void build(const OrderedTree& tree, std::uint32_t block_size);

  static constexpr std::uint32_t kUnmarked = 0xffffffffu;

  std::uint32_t block_ = 2;
  std::vector<std::uint32_t> depth_;
  std::vector<std::uint32_t> mark_id_;  ///< index among marks, kUnmarked for unmarked vertices
  // Micro trees (components of unmarked vertices).
  std::vector<std::uint32_t> micro_;        ///< micro tree of each unmarked vertex
  std::vector<std::uint8_t> local_pos_;     ///< local preorder position in its micro tree
  std::vector<std::uint32_t> micro_begin_;  ///< start of each micro tree in micro_vertex_
  std::vector<Vertex> micro_vertex_;        ///< local preorder position -> original vertex
  std::vector<std::uint32_t> micro_shape_;  ///< leading 1, then enter = 1 / leave = 0 bits
  std::vector<std::uint32_t> micro_top_depth_;
  std::vector<Vertex> micro_top_parent_;   ///< marked parent of the micro tree's top vertex
  std::vector<std::uint8_t> shape_table_;  ///< [shape][pos][relative depth] -> local pos
  std::size_t max_micro_size_ = 0;
  // Jump pointers at marks: jump_[jump_begin_[id] + i] is the ancestor 2^i edges above.
  std::vector<std::uint32_t> jump_begin_;
  std::vector<Vertex> jump_;
  // Ladders: maximum-height paths in depth order, each extended upward by its own length.
  std::vector<Vertex> ladder_;
  std::vector<std::uint32_t> ladder_pos_;  ///< position of v in its own path's extended array
  detail::FloorLogTable lg_;
};

}  // namespace sssp::paper
