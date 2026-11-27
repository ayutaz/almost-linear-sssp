# Dijkstra baselines

The baselines are the reference point for every performance claim about the paper's
algorithm, so they are written to be *fair* rather than clever: one shared scan loop, heaps
implemented in this repository (no `std::priority_queue`), all heap storage allocated before
the first push, and the same operation counters everywhere.

Public API: `include/sssp/baseline.hpp`. Sources: `src/baseline/`.

| Function | Heap (`src/baseline/heaps/`) | Decrease-key | Weight types | Time | Heap space |
|---|---|---|---|---|---|
| `dijkstra_binary_heap` | `BinaryHeap<K, V>` (lazy deletion) | no (re-push) | `uint64_t`, `double` | O((n + m) log m) | O(m) entries, reserved as `num_arcs + 1` |
| `dijkstra_4ary_heap` | `IndexedDaryHeap<K, 4>` | yes | `uint64_t`, `double` | O((n + m) log₄ n) | O(n) |
| `dijkstra_fibonacci_heap` | `FibonacciHeap<K>` | yes (cut + cascading cut) | `uint64_t`, `double` | O(m + n log n) amortised | O(n) node pool |
| `dijkstra_radix_heap` | `RadixHeap` (65 buckets) | yes (bucket move) | `uint64_t` only | O(m + 64 n) | O(n) node pool + 65 heads |

`n` = vertices, `m` = undirected edges; the CSR has `num_arcs` = 2m minus self-loops arcs.
Templates are explicitly instantiated in the `.cpp` files for `std::uint64_t` and `double`.

## Shared scan loop

`src/baseline/dijkstra_common.hpp`:

1. `init_result` throws `std::out_of_range` if `source >= n` (this includes the empty graph),
   then sets `dist = WeightTraits<W>::infinity()`, `parent = kNoVertex`,
   `parent_edge = kNoEdge`, `dist[source] = 0`.
2. Pop the minimum `(d_u, u)`; scan **every** arc `u -> v` of the CSR, including arcs to
   already-settled vertices and self-loops. No "settled" bitmap is consulted: with
   non-negative weights such arcs can never improve `dist[v]`, and skipping them would make the
   counters differ between implementations.
3. `nd = WeightTraits<W>::checked_add(d_u, w)`; if `nd < dist[v]` (strict), set `dist[v] = nd`,
   `parent[v] = u`, `parent_edge[v] = arc.edge` and push or decrease-key.

`dijkstra_binary_heap` uses the same loop with a lazy heap: an improved vertex is pushed again,
and a popped entry with `dist[u] < d` is stale and skipped. Because pushes happen only on
strict improvement, the final distance of a vertex is pushed exactly once, so each vertex is
scanned at most once.

## Tie-breaking

Heaps compare keys only (`a < b`); there is no secondary key. Among equal keys the pop order
is whatever the heap structure yields, which is deterministic for a given graph and source
(the CSR arc order is the edge-list order) but differs between heaps. Consequently:

- `dist` is identical across all baselines for integer weights;
- `parent` / `parent_edge` may differ between baselines when several shortest paths exist.
  Each is a valid shortest-path tree; tests check the tree structurally
  (`parent_edge` exists, joins `parent[v]` and `v`, `dist[v] == dist[parent] + w`, and the
  parent chain reaches the source) instead of comparing parents.
- The first improvement wins: an equally short path found later does not replace the parent.

## Zero weights

Zero-weight edges and zero-weight cycles are allowed. Strict improvement (`nd < dist[v]`)
guarantees that a zero-weight cycle cannot create a parent cycle: a vertex's parent is always
a vertex that was popped (settled) before it was last improved. The radix heap accepts keys
equal to the last popped key (bucket 0), which is exactly what a zero-weight arc produces.

## Integer weights and overflow

For `std::uint64_t`, `WeightTraits::checked_add` throws `std::overflow_error` when
`d_u + w >= 2^64 - 1` (the sum would reach `infinity()`). The check runs on **every** arc
scanned, so a run throws as soon as any candidate sum overflows, even if that candidate would
not have improved `dist[v]` (e.g. a huge edge back to the source). This is deliberate: the
result is never silently wrong. An edge of weight `2^64 - 1` therefore always throws once one
of its endpoints is scanned.

Which inputs throw (the same rule holds for `double`, see below):

- **Every Dijkstra baseline** throws on exactly the same inputs: iff some reachable `u` with a
  representable distance `d(u)` has an arc of weight `w` with `d(u) + w` not representable.
  The condition does not depend on the heap, because a vertex is scanned only with its final
  distance and every reachable vertex with a representable distance is scanned unless the run
  has already thrown.
- **`bellman_ford`** (the independent oracle, which does not call `checked_add`) skips an
  overflowing candidate and throws only iff some reachable vertex has a distance that is not
  representable. So whenever `bellman_ford` throws, every Dijkstra throws too, but not
  conversely: on the single edge `{0, 1, 2^64 - 2}` `bellman_ford` returns `{0, 2^64 - 2}`
  while every Dijkstra throws when scanning vertex 1 (the arc back to 0).

Both rules are pinned by `tests/differential/test_baselines_differential.cpp` (uint64, with an
independent saturating oracle) and `tests/differential/test_double_overflow_baselines.cpp`
(double, with an independent rounded oracle).

## Floating-point weights

For `double`, `checked_add` is the IEEE-754 addition (round to nearest; `-ffast-math` is never
enabled) of two finite distances, and it throws `std::overflow_error` if the rounded sum is not
finite. `+inf` is the "unreachable" sentinel, so an offer that rounds to `+inf` (e.g.
`1e308 + 1e308`, or `DBL_MAX + 2^970`, half an ulp of `DBL_MAX`) would otherwise be dropped by
the strict `nd < dist[v]` and a reachable vertex would be reported as unreachable (regression
`tests/regression/cases/20261127-double-sum-overflow.json`). A sum that is absorbed instead
(`DBL_MAX + 2^969 = DBL_MAX`) is finite and does not throw. The throw rule is the one of the
previous section with "representable" meaning "finite after rounding"; `bellman_ford` detects
it with its own `!= +inf` check. Properties relied on:

- `a + w >= a` for `w >= 0`, so arcs to settled vertices never improve them and the heap keys
  are popped in non-decreasing order;
- `dist[v]` is stored as exactly the value `dist[parent[v]] + w` that was compared, so the tree
  check `dist[v] == dist[parent] + w` holds *exactly*;
- distances are sums along the path in path order. A different algorithm (e.g. Floyd–Warshall,
  which adds sub-path sums) can round differently, so cross-algorithm comparisons of `double`
  results need a tolerance unless the weights are dyadic (the tests use multiples of 1/16 for
  exact comparison and a 1e-12 relative tolerance for general reals).

The radix heap is integer-only; it is not instantiated for `double`.

## Diagnostics (`SsspResult::diag.ops`)

All counters are incremented through `SSSP_COUNT` and compile to nothing when
`SSSP_DIAGNOSTICS=0` (release / timed builds), in which case they stay 0.

| Counter | Definition |
|---|---|
| `vertex_scans` | vertices whose adjacency list was scanned (= reachable vertices) |
| `edge_scans` | arcs inspected (sum of CSR degrees of scanned vertices; an undirected edge counts once per scanned endpoint, a self-loop once) |
| `relaxations` | successful strict improvements of `dist[v]` |
| `additions` | `d_u + w` evaluations (= `edge_scans`) |
| `comparisons` | `nd < dist[v]` (one per arc) + one staleness check per pop (binary heap only) + key comparisons made inside the heap (`heap.comparisons()`) |
| `queue_pushes` | heap inserts (binary: 1 + `relaxations`; indexed: 1 + number of first-time reached vertices) |
| `queue_pops` | extract-min calls, including stale ones (binary: = `queue_pushes`; indexed: = `vertex_scans`) |
| `decrease_keys` | decrease-key calls (binary: always 0; indexed: `relaxations + 1 - queue_pushes`) |

Heap-internal key comparisons:

- `BinaryHeap`, `IndexedDaryHeap`: every `key < key` in sift-up / sift-down
  (d-ary sift-down: up to D-1 to find the smallest child plus 1 against the moving item).
- `FibonacciHeap`: comparison with the minimum on insert (not for the first root) and on
  decrease-key, parent comparison on decrease-key, one per link during consolidation, and the
  minimum search over the consolidated roots.
- `RadixHeap`: only the minimum search over the bucket being redistributed (k − 1 comparisons
  for k items). Bucket indices are computed with `bit_width(key XOR last)`, which is not
  counted as a comparison. This makes radix-heap comparison counts not directly comparable to
  the comparison-based heaps; compare `edge_scans` / time instead.

Argument checks inside the heaps (id range, membership, key increase, radix monotonicity)
are not counted.

Hand-computed example (classic 6-vertex graph in
`tests/unit/baseline/dijkstra_test_support.hpp`, source 0): 6 vertex scans, 18 arc scans and
additions, 7 relaxations; binary heap 8 pushes / 8 pops (2 stale); indexed heaps 6 pushes,
2 decrease-keys, 6 pops.

## Allocation policy

No heap operation allocates during a run:

- `BinaryHeap` reserves `num_arcs + 1` entries (an upper bound on pushes) before the first push.
- `IndexedDaryHeap` allocates its heap array, position map and key array (size n) in the
  constructor.
- `FibonacciHeap` allocates one node per vertex and an n-entry root scratch buffer in the
  constructor; the degree table is a fixed `std::array<id, 64>` (max degree
  < 1.4405 log₂ n + 2 < 64 for n < 2³²).
- `RadixHeap` allocates one node per vertex (intrusive doubly-linked bucket lists) in the
  constructor; bucket heads are a fixed `std::array<id, 65>`.

The result vectors (`dist`, `parent`, `parent_edge`) are allocated once by `init_result`.

## Heap misuse

All heaps throw `std::out_of_range` for pop on an empty heap or an id outside `[0, n)`, and
`std::invalid_argument` for pushing a present id, decrease-key of an absent id or to a larger
key, and (radix) a key below the last popped key. Dijkstra never triggers these.

## Tests

- Heaps: `tests/unit/baseline/heaps/` (hand-computed pop orders, cascading cut, consolidation,
  bucket behaviour, counters) and `tests/property/baseline/` (random operation sequences vs a
  sorted `std::multiset` / `std::map` reference, seeds {1, 7, 42, 2026, 0xC0FFEE}, invariant
  validators after operations).
- Dijkstra: `tests/unit/baseline/test_dijkstra_*.cpp` with shared checks in
  `dijkstra_test_support.hpp`: hand-computed graphs (classic, zero weights / multi-edges /
  self-loops / disconnected / isolated, zero-weight cycle, single vertex, dyadic doubles),
  invalid sources, overflow, counters, and a test-only Floyd–Warshall on random multigraphs
  (seeds {1, 2, 3, 42, 2026}, 125 graphs with n ≤ 30 from every source plus one n = 200 graph
  per seed).
