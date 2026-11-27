// Definition 4.1 / 4.2 (Sec 4, l.646-672): builders and the independent checker.
// The checker is a verification tool, not a step of the paper's algorithm; it may use a
// heap-based Dijkstra and naive ancestor walks because its cost is not part of any bound.

#include "paper/hierarchy/ordered_hierarchy.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace sssp::paper {

namespace {

template <class W>
bool valid_length(W w) noexcept {
  if constexpr (std::is_floating_point_v<W>) {
    return std::isfinite(w) && w >= W{0};
  } else {
    return w != WeightTraits<W>::infinity();
  }
}

template <class W>
bool positive_height(W h) noexcept {
  return valid_length(h) && h > W{0};
}

template <class W>
W saturating_add(W a, W b) noexcept {
  if constexpr (std::is_floating_point_v<W>) {
    return a + b;
  } else {
    return a >= WeightTraits<W>::infinity() - b ? WeightTraits<W>::infinity()
                                                : static_cast<W>(a + b);
  }
}

// c0 * h by c0 - 1 additions (Sec 2: reals are only added and compared).
template <class W>
W times_constant(W h, std::uint32_t c0) noexcept {
  if (c0 == 0) return W{0};
  W r = h;
  for (std::uint32_t i = 1; i < c0; ++i) r = saturating_add(r, h);
  return r;
}

struct SimpleDsu {
  std::vector<std::size_t> p;
  explicit SimpleDsu(std::size_t n) : p(n) { std::iota(p.begin(), p.end(), std::size_t{0}); }
  std::size_t find(std::size_t x) {
    while (p[x] != x) x = p[x] = p[p[x]];
    return x;
  }
  bool unite(std::size_t a, std::size_t b) {
    a = find(a);
    b = find(b);
    if (a == b) return false;
    p[a] = b;
    return true;
  }
};

struct QuotientEdge {
  std::size_t from;
  std::size_t to;
  ArcId arc;
};

// Kahn's algorithm on a quotient with k vertices. Returns kNoArc if acyclic, otherwise an arc
// whose endpoints both stay unsorted (it lies on or downstream of a cycle).
ArcId find_cycle_witness(std::size_t k, const std::vector<QuotientEdge>& edges) {
  std::vector<std::size_t> indeg(k, 0);
  std::vector<std::vector<std::size_t>> out(k);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    ++indeg[edges[i].to];
    out[edges[i].from].push_back(i);
  }
  std::vector<std::size_t> stack;
  for (std::size_t v = 0; v < k; ++v) {
    if (indeg[v] == 0) stack.push_back(v);
  }
  std::size_t sorted = 0;
  while (!stack.empty()) {
    const std::size_t v = stack.back();
    stack.pop_back();
    ++sorted;
    for (const std::size_t e : out[v]) {
      if (--indeg[edges[e].to] == 0) stack.push_back(edges[e].to);
    }
  }
  if (sorted == k) return kNoArc;
  for (const auto& e : edges) {
    if (indeg[e.from] > 0 && indeg[e.to] > 0) return e.arc;
  }
  return edges.front().arc;  // unreachable: an unsorted vertex has an unsorted predecessor
}

// Exact physical distances from a set of weighted sources (checker oracle).
template <class W>
std::vector<W> physical_distances(const CsrGraph<W>& g,
                                  const std::vector<std::pair<Vertex, W>>& sources) {
  const W inf = WeightTraits<W>::infinity();
  std::vector<W> dist(g.num_vertices(), inf);
  using Item = std::pair<W, Vertex>;
  std::priority_queue<Item, std::vector<Item>, std::greater<>> pq;
  for (const auto& [v, c] : sources) {
    if (c < dist[v]) {
      dist[v] = c;
      pq.emplace(c, v);
    }
  }
  while (!pq.empty()) {
    const auto [d, u] = pq.top();
    pq.pop();
    if (d > dist[u]) continue;
    for (const auto& a : g.arcs(u)) {
      const W nd = saturating_add(d, a.weight);
      if (nd < dist[a.to]) {
        dist[a.to] = nd;
        pq.emplace(nd, a.to);
      }
    }
  }
  return dist;
}

template <class W>
class Checker {
 public:
  Checker(const OrderedHierarchy<W>& h, const CsrGraph<W>& phys, ValidationReport& rep)
      : h_(h), phys_(phys), rep_(rep) {}

  bool check_structure();
  void prepare();
  void check_h1();
  void check_h2();
  void check_h3_h4();
  void check_h5_h6();

 private:
  void add(ViolationKind k, NodeId node, StateId state, ArcId arc, std::size_t link,
           std::string msg) {
    rep_.violations.push_back(Violation{k, node, state, arc, link, std::move(msg)});
  }

  // Is there a physical walk a ~> b of length <= len?
  bool walk_within(Vertex a, Vertex b, W len) {
    if (a == b) return true;
    for (const auto& arc : phys_.arcs(a)) {
      if (arc.to == b && arc.weight <= len) return true;
    }
    auto it = dist_cache_.find(a);
    if (it == dist_cache_.end()) {
      it = dist_cache_.emplace(a, physical_distances(phys_, {{a, W{0}}})).first;
    }
    const W d = it->second[b];
    return d != WeightTraits<W>::infinity() && d <= len;
  }

  NodeId lca(NodeId a, NodeId b) const {
    if (root_[a] != root_[b]) return kNoNode;
    while (depth_[a] > depth_[b]) a = h_.nodes[a].parent;
    while (depth_[b] > depth_[a]) b = h_.nodes[b].parent;
    while (a != b) {
      a = h_.nodes[a].parent;
      b = h_.nodes[b].parent;
    }
    return a;
  }

  // The child of `c` whose subtree contains state x, or kNoNode if x is not below c.
  NodeId child_below(NodeId c, StateId x) const {
    if (x >= h_.num_states()) return kNoNode;
    NodeId v = h_.leaf_of_state[x];
    if (depth_[v] <= depth_[c]) return kNoNode;
    while (depth_[v] > depth_[c] + 1) v = h_.nodes[v].parent;
    return h_.nodes[v].parent == c ? v : kNoNode;
  }

  const OrderedHierarchy<W>& h_;
  const CsrGraph<W>& phys_;
  ValidationReport& rep_;
  std::vector<std::size_t> depth_;
  std::vector<NodeId> root_;
  std::vector<std::size_t> pos_in_parent_;
  std::vector<std::size_t> root_index_;
  std::map<Vertex, std::vector<W>> dist_cache_;
};

template <class W>
bool Checker<W>::check_structure() {
  const std::size_t before = rep_.violations.size();
  const std::size_t n = h_.num_states();
  const std::size_t nn = h_.nodes.size();
  if (h_.projection.size() != n) {
    add(ViolationKind::kStateCountMismatch, kNoNode, kNoState, kNoArc, kNoIndex,
        "projection has " + std::to_string(h_.projection.size()) + " entries for " +
            std::to_string(n) + " states");
  }
  if (h_.leaf_of_state.size() != n) {
    add(ViolationKind::kStateCountMismatch, kNoNode, kNoState, kNoArc, kNoIndex,
        "leaf_of_state has " + std::to_string(h_.leaf_of_state.size()) + " entries for " +
            std::to_string(n) + " states");
  }
  for (std::size_t x = 0; x < h_.projection.size(); ++x) {
    if (h_.projection[x] >= phys_.num_vertices()) {
      add(ViolationKind::kProjectionOutOfRange, kNoNode, static_cast<StateId>(x), kNoArc, kNoIndex,
          "pi(" + std::to_string(x) + ") is not a physical vertex");
    }
  }
  for (std::size_t i = 0; i < h_.graph.arcs.size(); ++i) {
    const auto& a = h_.graph.arcs[i];
    if (a.tail >= n || a.head >= n) {
      add(ViolationKind::kArcEndpointOutOfRange, kNoNode, kNoState, static_cast<ArcId>(i), kNoIndex,
          "arc " + std::to_string(i) + " has an endpoint that is not a state");
    }
  }
  // Parent/child consistency and leaf/state correspondence.
  bool indices_ok = true;
  std::vector<std::size_t> occurrences(nn, 0);
  std::vector<std::size_t> leaves_of(n, 0);
  for (std::size_t v = 0; v < nn; ++v) {
    const auto& node = h_.nodes[v];
    const auto id = static_cast<NodeId>(v);
    if (node.parent != kNoNode && node.parent >= nn) {
      indices_ok = false;
      add(ViolationKind::kNodeIndexOutOfRange, id, kNoState, kNoArc, kNoIndex,
          "parent of node " + std::to_string(v) + " is out of range");
    }
    for (const NodeId c : node.children) {
      if (c >= nn) {
        indices_ok = false;
        add(ViolationKind::kNodeIndexOutOfRange, id, kNoState, kNoArc, kNoIndex,
            "child of node " + std::to_string(v) + " is out of range");
        continue;
      }
      ++occurrences[c];
      if (h_.nodes[c].parent != id) {
        indices_ok = false;
        add(ViolationKind::kParentChildMismatch, c, kNoState, kNoArc, kNoIndex,
            "node " + std::to_string(c) + " is listed as child of " + std::to_string(v) +
                " but its parent pointer disagrees");
      }
    }
    if (node.children.empty()) {
      if (node.state >= n) {
        add(ViolationKind::kLeafMismatch, id, node.state, kNoArc, kNoIndex,
            "leaf " + std::to_string(v) + " does not carry a valid state");
      } else {
        ++leaves_of[node.state];
      }
    } else if (node.state != kNoState) {
      add(ViolationKind::kLeafMismatch, id, node.state, kNoArc, kNoIndex,
          "internal node " + std::to_string(v) + " carries a state");
    }
  }
  for (std::size_t v = 0; v < nn; ++v) {
    const bool has_parent = h_.nodes[v].parent != kNoNode;
    if (occurrences[v] != (has_parent ? 1U : 0U)) {
      indices_ok = false;
      add(ViolationKind::kParentChildMismatch, static_cast<NodeId>(v), kNoState, kNoArc, kNoIndex,
          "node " + std::to_string(v) + " appears " + std::to_string(occurrences[v]) +
              " times in child lists");
    }
  }
  for (std::size_t x = 0; x < n; ++x) {
    const auto sx = static_cast<StateId>(x);
    if (leaves_of[x] != 1) {
      add(ViolationKind::kLeafMismatch, kNoNode, sx, kNoArc, kNoIndex,
          "state " + std::to_string(x) + " has " + std::to_string(leaves_of[x]) + " leaves");
    }
    if (x < h_.leaf_of_state.size()) {
      const NodeId l = h_.leaf_of_state[x];
      if (l >= nn || !h_.nodes[l].children.empty() || h_.nodes[l].state != sx) {
        add(ViolationKind::kLeafMismatch, l, sx, kNoArc, kNoIndex,
            "leaf_of_state[" + std::to_string(x) + "] is not the leaf of that state");
      }
    }
  }
  if (indices_ok) {
    // Parent pointers must reach a root without repeating a node.
    std::vector<std::uint8_t> colour(nn, 0);  // 0 new, 1 on current path, 2 done
    for (std::size_t s = 0; s < nn; ++s) {
      std::vector<std::size_t> path;
      std::size_t v = s;
      bool cycle = false;
      while (colour[v] == 0) {
        colour[v] = 1;
        path.push_back(v);
        if (h_.nodes[v].parent == kNoNode) break;
        v = h_.nodes[v].parent;
        if (colour[v] == 1) {
          cycle = true;
          break;
        }
      }
      for (const std::size_t p : path) colour[p] = 2;
      if (cycle) {
        indices_ok = false;
        add(ViolationKind::kParentCycle, static_cast<NodeId>(v), kNoState, kNoArc, kNoIndex,
            "parent pointers through node " + std::to_string(v) + " form a cycle");
      }
    }
    std::vector<std::uint8_t> listed(nn, 0);
    bool roots_ok = true;
    for (const NodeId r : h_.roots) {
      if (r >= nn || h_.nodes[r].parent != kNoNode || listed[r] != 0) {
        roots_ok = false;
        break;
      }
      listed[r] = 1;
    }
    for (std::size_t v = 0; v < nn && roots_ok; ++v) {
      if (h_.nodes[v].parent == kNoNode && listed[v] == 0) roots_ok = false;
    }
    if (!roots_ok) {
      add(ViolationKind::kRootsMismatch, kNoNode, kNoState, kNoArc, kNoIndex,
          "roots is not exactly the set of parentless nodes");
    }
  }
  return rep_.violations.size() == before;
}

template <class W>
void Checker<W>::prepare() {
  const std::size_t nn = h_.nodes.size();
  depth_.assign(nn, 0);
  root_.assign(nn, kNoNode);
  pos_in_parent_.assign(nn, 0);
  root_index_.assign(nn, kNoIndex);
  std::vector<NodeId> order;
  order.reserve(nn);
  for (std::size_t i = 0; i < h_.roots.size(); ++i) {
    const NodeId r = h_.roots[i];
    root_[r] = r;
    root_index_[r] = i;
    order.push_back(r);
  }
  for (std::size_t i = 0; i < order.size(); ++i) {
    const NodeId v = order[i];
    const auto& ch = h_.nodes[v].children;
    for (std::size_t j = 0; j < ch.size(); ++j) {
      depth_[ch[j]] = depth_[v] + 1;
      root_[ch[j]] = root_[v];
      pos_in_parent_[ch[j]] = j;
      order.push_back(ch[j]);
    }
  }
}

// H1: each arc projects to a physical walk no longer than the arc.
template <class W>
void Checker<W>::check_h1() {
  for (std::size_t i = 0; i < h_.graph.arcs.size(); ++i) {
    const auto& a = h_.graph.arcs[i];
    const auto id = static_cast<ArcId>(i);
    if (!valid_length(a.length)) {
      add(ViolationKind::kArcLengthInvalid, kNoNode, a.tail, id, kNoIndex,
          "arc " + std::to_string(i) + " has a negative, NaN or infinite length");
      continue;
    }
    if (!walk_within(h_.projection[a.tail], h_.projection[a.head], a.length)) {
      add(ViolationKind::kArcWalkTooLong, kNoNode, a.tail, id, kNoIndex,
          "arc " + std::to_string(i) + " is shorter than every physical walk between its images");
    }
  }
}

// H2: >= 2 children, positive height, heights nondecreasing towards the root.
template <class W>
void Checker<W>::check_h2() {
  for (std::size_t v = 0; v < h_.nodes.size(); ++v) {
    const auto& node = h_.nodes[v];
    if (node.children.empty()) continue;
    const auto id = static_cast<NodeId>(v);
    if (node.children.size() < 2) {
      add(ViolationKind::kTooFewChildren, id, kNoState, kNoArc, kNoIndex,
          "internal node " + std::to_string(v) + " has fewer than two children");
    }
    if (!positive_height(node.height)) {
      add(ViolationKind::kNonPositiveHeight, id, kNoState, kNoArc, kNoIndex,
          "internal node " + std::to_string(v) + " has a non-positive or non-finite height");
      continue;
    }
    if (node.parent != kNoNode) {
      const W ph = h_.nodes[node.parent].height;
      if (positive_height(ph) && node.height > ph) {
        add(ViolationKind::kHeightDecreases, id, kNoState, kNoArc, kNoIndex,
            "node " + std::to_string(v) + " is higher than its parent " +
                std::to_string(node.parent));
      }
    }
  }
}

// H3 (light arcs at their LCA) and H4 (all inter-root arcs).
template <class W>
void Checker<W>::check_h3_h4() {
  std::vector<std::vector<QuotientEdge>> light(h_.nodes.size());
  std::vector<QuotientEdge> roots;
  for (std::size_t i = 0; i < h_.graph.arcs.size(); ++i) {
    const auto& a = h_.graph.arcs[i];
    if (a.tail == a.head) continue;
    const NodeId lx = h_.leaf_of_state[a.tail];
    const NodeId ly = h_.leaf_of_state[a.head];
    const NodeId l = lca(lx, ly);
    const auto id = static_cast<ArcId>(i);
    if (l == kNoNode) {
      roots.push_back({root_index_[root_[lx]], root_index_[root_[ly]], id});
      continue;
    }
    // l is internal because the two leaves differ. A NaN length or height compares false and
    // is reported by H1/H2 instead.
    if (a.length < h_.nodes[l].height) {
      const NodeId ci = child_below(l, a.tail);
      const NodeId cj = child_below(l, a.head);
      light[l].push_back({pos_in_parent_[ci], pos_in_parent_[cj], id});
    }
  }
  for (std::size_t v = 0; v < h_.nodes.size(); ++v) {
    if (light[v].empty()) continue;
    const ArcId w = find_cycle_witness(h_.nodes[v].children.size(), light[v]);
    if (w != kNoArc) {
      add(ViolationKind::kLightArcCycle, static_cast<NodeId>(v), kNoState, w, kNoIndex,
          "light arcs between the children of node " + std::to_string(v) +
              " form a cycle (witness arc " + std::to_string(w) + ")");
    }
  }
  if (!roots.empty()) {
    const ArcId w = find_cycle_witness(h_.roots.size(), roots);
    if (w != kNoArc) {
      add(ViolationKind::kRootArcCycle, kNoNode, kNoState, w, kNoIndex,
          "inter-root arcs form a cycle (witness arc " + std::to_string(w) + ")");
    }
  }
}

// H5 (link tree per internal node) and H6 (claimed backing arcs are genuine).
template <class W>
void Checker<W>::check_h5_h6() {
  const std::size_t n = h_.num_states();
  for (std::size_t v = 0; v < h_.nodes.size(); ++v) {
    const auto& node = h_.nodes[v];
    if (node.children.empty()) continue;
    const auto id = static_cast<NodeId>(v);
    const std::size_t d = node.children.size();
    if (node.links.size() + 1 != d) {
      add(ViolationKind::kLinkCount, id, kNoState, kNoArc, kNoIndex,
          "node " + std::to_string(v) + " has " + std::to_string(node.links.size()) +
              " links for " + std::to_string(d) + " children");
    }
    const bool height_ok = positive_height(node.height);
    const W bound = height_ok ? times_constant(node.height, h_.c0) : W{0};
    SimpleDsu dsu(d);
    std::size_t components = d;
    for (std::size_t j = 0; j < node.links.size(); ++j) {
      const auto& lk = node.links[j];
      const bool len_ok = valid_length(lk.length);
      if (!len_ok) {
        add(ViolationKind::kLinkLengthInvalid, id, lk.x, kNoArc, j,
            "link " + std::to_string(j) + " has a negative, NaN or infinite length");
      }
      const NodeId ci = child_below(id, lk.x);
      const NodeId cj = child_below(id, lk.y);
      if (ci == kNoNode || cj == kNoNode) {
        add(ViolationKind::kLinkEndpointOutside, id, ci == kNoNode ? lk.x : lk.y, kNoArc, j,
            "link " + std::to_string(j) + " has an endpoint that is not below node " +
                std::to_string(v));
      } else if (ci == cj) {
        add(ViolationKind::kLinkSameChild, id, lk.x, kNoArc, j,
            "link " + std::to_string(j) + " joins two states below the same child");
      } else if (dsu.unite(pos_in_parent_[ci], pos_in_parent_[cj])) {
        --components;
      } else {
        add(ViolationKind::kLinkCycle, id, lk.x, kNoArc, j,
            "link " + std::to_string(j) + " closes a cycle among the children");
      }
      if (len_ok && height_ok && lk.length > bound) {
        add(ViolationKind::kLinkTooLong, id, lk.x, kNoArc, j,
            "link " + std::to_string(j) + " is longer than c0 * h(node " + std::to_string(v) + ")");
      }
      if (len_ok && lk.x < n && lk.y < n &&
          !walk_within(h_.projection[lk.x], h_.projection[lk.y], lk.length)) {
        add(ViolationKind::kLinkWalkTooLong, id, lk.x, kNoArc, j,
            "link " + std::to_string(j) + " is shorter than every physical walk");
      }
      // H6: a link that names a backing state arc must be realised by that arc.
      if (lk.arc != kNoArc) {
        if (lk.arc >= h_.graph.arcs.size()) {
          add(ViolationKind::kLinkArcOutOfRange, id, lk.x, lk.arc, j,
              "link " + std::to_string(j) + " names a non-existent backing arc");
        } else {
          const auto& a = h_.graph.arcs[lk.arc];
          const bool same_pair =
              (a.tail == lk.x && a.head == lk.y) || (a.tail == lk.y && a.head == lk.x);
          const bool shorter = len_ok && valid_length(a.length) && a.length <= lk.length;
          if (!same_pair || !shorter) {
            add(ViolationKind::kLinkArcMismatch, id, lk.x, lk.arc, j,
                "link " + std::to_string(j) + " is not realised by its backing arc " +
                    std::to_string(lk.arc));
          }
        }
      }
    }
    if (components > 1) {
      add(ViolationKind::kLinksNotSpanning, id, kNoState, kNoArc, kNoIndex,
          "links leave " + std::to_string(components) + " components among the children of " +
              std::to_string(v));
    }
  }
}

template <class W>
bool is_finite_offer(W w) noexcept {
  if constexpr (std::is_floating_point_v<W>) {
    return std::isfinite(w);
  } else {
    return w != WeightTraits<W>::infinity();
  }
}

}  // namespace

// ------------------------------------------------------------------ builders

template <class W>
OrderedHierarchy<W> make_leaf_forest(StateGraph<W> graph, std::vector<Vertex> projection) {
  if (projection.size() != graph.num_states) {
    throw std::invalid_argument("make_leaf_forest: projection size differs from state count");
  }
  OrderedHierarchy<W> h;
  const std::size_t n = graph.num_states;
  h.graph = std::move(graph);
  h.projection = std::move(projection);
  h.nodes.resize(n);
  h.leaf_of_state.resize(n);
  for (std::size_t x = 0; x < n; ++x) {
    h.nodes[x].state = static_cast<StateId>(x);
    h.leaf_of_state[x] = static_cast<NodeId>(x);
  }
  collect_roots(h);
  return h;
}

template <class W>
NodeId add_internal_node(OrderedHierarchy<W>& h, W height, std::vector<NodeId> children,
                         std::vector<PhysicalLink<W>> links) {
  const auto id = static_cast<NodeId>(h.nodes.size());
  for (const NodeId c : children) {
    if (c >= h.nodes.size()) throw std::invalid_argument("add_internal_node: child out of range");
    if (h.nodes[c].parent != kNoNode) {
      throw std::invalid_argument("add_internal_node: child already has a parent");
    }
  }
  for (const NodeId c : children) h.nodes[c].parent = id;
  HierarchyNode<W> node;
  node.children = std::move(children);
  node.height = height;
  node.links = std::move(links);
  h.nodes.push_back(std::move(node));
  return id;
}

template <class W>
void collect_roots(OrderedHierarchy<W>& h) {
  h.roots.clear();
  for (std::size_t v = 0; v < h.nodes.size(); ++v) {
    if (h.nodes[v].parent == kNoNode) h.roots.push_back(static_cast<NodeId>(v));
  }
}

// ------------------------------------------------------------------ report

Condition condition_of(ViolationKind kind) noexcept {
  switch (kind) {
    case ViolationKind::kStateCountMismatch:
    case ViolationKind::kProjectionOutOfRange:
    case ViolationKind::kArcEndpointOutOfRange:
    case ViolationKind::kNodeIndexOutOfRange:
    case ViolationKind::kLeafMismatch:
    case ViolationKind::kParentChildMismatch:
    case ViolationKind::kParentCycle:
    case ViolationKind::kRootsMismatch:
      return Condition::kStructure;
    case ViolationKind::kArcLengthInvalid:
    case ViolationKind::kArcWalkTooLong:
      return Condition::kH1;
    case ViolationKind::kTooFewChildren:
    case ViolationKind::kNonPositiveHeight:
    case ViolationKind::kHeightDecreases:
      return Condition::kH2;
    case ViolationKind::kLightArcCycle:
      return Condition::kH3;
    case ViolationKind::kRootArcCycle:
      return Condition::kH4;
    case ViolationKind::kLinkLengthInvalid:
    case ViolationKind::kLinkCount:
    case ViolationKind::kLinkEndpointOutside:
    case ViolationKind::kLinkSameChild:
    case ViolationKind::kLinkCycle:
    case ViolationKind::kLinksNotSpanning:
    case ViolationKind::kLinkTooLong:
    case ViolationKind::kLinkWalkTooLong:
      return Condition::kH5;
    case ViolationKind::kLinkArcOutOfRange:
    case ViolationKind::kLinkArcMismatch:
      return Condition::kH6;
    case ViolationKind::kSourceOutOfRange:
    case ViolationKind::kSourceOffsetInvalid:
    case ViolationKind::kOfferStateOutOfRange:
    case ViolationKind::kOfferNotFinite:
    case ViolationKind::kOfferInfeasible:
      return Condition::kSeed;
  }
  return Condition::kStructure;
}

const char* to_string(Condition c) noexcept {
  switch (c) {
    case Condition::kStructure:
      return "structure";
    case Condition::kH1:
      return "H1";
    case Condition::kH2:
      return "H2";
    case Condition::kH3:
      return "H3";
    case Condition::kH4:
      return "H4";
    case Condition::kH5:
      return "H5";
    case Condition::kH6:
      return "H6";
    case Condition::kSeed:
      return "Def4.2";
  }
  return "?";
}

const char* to_string(ViolationKind kind) noexcept {
  switch (kind) {
    case ViolationKind::kStateCountMismatch:
      return "StateCountMismatch";
    case ViolationKind::kProjectionOutOfRange:
      return "ProjectionOutOfRange";
    case ViolationKind::kArcEndpointOutOfRange:
      return "ArcEndpointOutOfRange";
    case ViolationKind::kNodeIndexOutOfRange:
      return "NodeIndexOutOfRange";
    case ViolationKind::kLeafMismatch:
      return "LeafMismatch";
    case ViolationKind::kParentChildMismatch:
      return "ParentChildMismatch";
    case ViolationKind::kParentCycle:
      return "ParentCycle";
    case ViolationKind::kRootsMismatch:
      return "RootsMismatch";
    case ViolationKind::kArcLengthInvalid:
      return "ArcLengthInvalid";
    case ViolationKind::kArcWalkTooLong:
      return "ArcWalkTooLong";
    case ViolationKind::kTooFewChildren:
      return "TooFewChildren";
    case ViolationKind::kNonPositiveHeight:
      return "NonPositiveHeight";
    case ViolationKind::kHeightDecreases:
      return "HeightDecreases";
    case ViolationKind::kLightArcCycle:
      return "LightArcCycle";
    case ViolationKind::kRootArcCycle:
      return "RootArcCycle";
    case ViolationKind::kLinkLengthInvalid:
      return "LinkLengthInvalid";
    case ViolationKind::kLinkCount:
      return "LinkCount";
    case ViolationKind::kLinkEndpointOutside:
      return "LinkEndpointOutside";
    case ViolationKind::kLinkSameChild:
      return "LinkSameChild";
    case ViolationKind::kLinkCycle:
      return "LinkCycle";
    case ViolationKind::kLinksNotSpanning:
      return "LinksNotSpanning";
    case ViolationKind::kLinkTooLong:
      return "LinkTooLong";
    case ViolationKind::kLinkWalkTooLong:
      return "LinkWalkTooLong";
    case ViolationKind::kLinkArcOutOfRange:
      return "LinkArcOutOfRange";
    case ViolationKind::kLinkArcMismatch:
      return "LinkArcMismatch";
    case ViolationKind::kSourceOutOfRange:
      return "SourceOutOfRange";
    case ViolationKind::kSourceOffsetInvalid:
      return "SourceOffsetInvalid";
    case ViolationKind::kOfferStateOutOfRange:
      return "OfferStateOutOfRange";
    case ViolationKind::kOfferNotFinite:
      return "OfferNotFinite";
    case ViolationKind::kOfferInfeasible:
      return "OfferInfeasible";
  }
  return "?";
}

bool ValidationReport::has(Condition c) const noexcept {
  return std::any_of(violations.begin(), violations.end(),
                     [c](const Violation& v) { return v.condition() == c; });
}

bool ValidationReport::has(ViolationKind k) const noexcept {
  return std::any_of(violations.begin(), violations.end(),
                     [k](const Violation& v) { return v.kind == k; });
}

std::size_t ValidationReport::count(Condition c) const noexcept {
  return static_cast<std::size_t>(
      std::count_if(violations.begin(), violations.end(),
                    [c](const Violation& v) { return v.condition() == c; }));
}

std::string ValidationReport::summary() const {
  if (violations.empty()) return "ok";
  std::ostringstream os;
  for (const auto& v : violations) {
    os << '[' << to_string(v.condition()) << ' ' << to_string(v.kind) << "] " << v.message << '\n';
  }
  return os.str();
}

// ------------------------------------------------------------------ checker entry points

template <class W>
ValidationReport validate_ordered_hierarchy(const OrderedHierarchy<W>& h,
                                            const CsrGraph<W>& physical) {
  ValidationReport rep;
  Checker<W> checker(h, physical, rep);
  if (!checker.check_structure()) return rep;
  checker.prepare();
  checker.check_h1();
  checker.check_h2();
  checker.check_h3_h4();
  checker.check_h5_h6();
  return rep;
}

template <class W>
ValidationReport validate_seeded_ordered_hierarchy(const SeededOrderedHierarchy<W>& s,
                                                   const CsrGraph<W>& physical) {
  ValidationReport rep = validate_ordered_hierarchy(s.hierarchy, physical);
  const auto& h = s.hierarchy;
  std::vector<std::pair<Vertex, W>> sources;
  for (std::size_t i = 0; i < s.sources.size(); ++i) {
    const auto& src = s.sources[i];
    if (src.vertex >= physical.num_vertices()) {
      rep.violations.push_back(Violation{ViolationKind::kSourceOutOfRange, kNoNode, kNoState,
                                         kNoArc, i,
                                         "source " + std::to_string(i) + " is out of range"});
    } else if (!valid_length(src.offset)) {
      rep.violations.push_back(
          Violation{ViolationKind::kSourceOffsetInvalid, kNoNode, kNoState, kNoArc, i,
                    "source " + std::to_string(i) + " has a negative, NaN or infinite offset"});
    } else {
      sources.emplace_back(src.vertex, src.offset);
    }
  }
  // delta(v) = min_i (c_i + d(v_i, v)) (Def 4.2).
  const std::vector<W> delta = physical_distances(physical, sources);
  for (std::size_t i = 0; i < s.offers.size(); ++i) {
    const auto& q = s.offers[i];
    if (q.state >= h.num_states()) {
      rep.violations.push_back(Violation{ViolationKind::kOfferStateOutOfRange, kNoNode, q.state,
                                         kNoArc, i,
                                         "offer " + std::to_string(i) + " names no state"});
      continue;
    }
    if (!is_finite_offer(q.value)) {
      rep.violations.push_back(Violation{ViolationKind::kOfferNotFinite, kNoNode, q.state, kNoArc,
                                         i, "offer " + std::to_string(i) + " is not finite"});
      continue;
    }
    if (q.state >= h.projection.size() || h.projection[q.state] >= physical.num_vertices()) {
      continue;  // already reported as a structural violation
    }
    if (q.value < delta[h.projection[q.state]]) {
      rep.violations.push_back(
          Violation{ViolationKind::kOfferInfeasible, kNoNode, q.state, kNoArc, i,
                    "offer " + std::to_string(i) + " at state " + std::to_string(q.state) +
                        " is below the target distance of its image"});
    }
  }
  return rep;
}

// Explicit instantiations for std::uint64_t.
template OrderedHierarchy<std::uint64_t> make_leaf_forest<std::uint64_t>(StateGraph<std::uint64_t>,
                                                                         std::vector<Vertex>);
template NodeId add_internal_node<std::uint64_t>(OrderedHierarchy<std::uint64_t>&, std::uint64_t,
                                                 std::vector<NodeId>,
                                                 std::vector<PhysicalLink<std::uint64_t>>);
template void collect_roots<std::uint64_t>(OrderedHierarchy<std::uint64_t>&);
template ValidationReport validate_ordered_hierarchy<std::uint64_t>(
    const OrderedHierarchy<std::uint64_t>&, const CsrGraph<std::uint64_t>&);
template ValidationReport validate_seeded_ordered_hierarchy<std::uint64_t>(
    const SeededOrderedHierarchy<std::uint64_t>&, const CsrGraph<std::uint64_t>&);

// Explicit instantiations for double.
template OrderedHierarchy<double> make_leaf_forest<double>(StateGraph<double>, std::vector<Vertex>);
template NodeId add_internal_node<double>(OrderedHierarchy<double>&, double, std::vector<NodeId>,
                                          std::vector<PhysicalLink<double>>);
template void collect_roots<double>(OrderedHierarchy<double>&);
template ValidationReport validate_ordered_hierarchy<double>(const OrderedHierarchy<double>&,
                                                             const CsrGraph<double>&);
template ValidationReport validate_seeded_ordered_hierarchy<double>(
    const SeededOrderedHierarchy<double>&, const CsrGraph<double>&);

}  // namespace sssp::paper
