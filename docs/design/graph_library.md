# Graph library

Graph representations, text IO and reproducible generators shared by every SSSP
implementation, benchmark and test. Graphs are **undirected** with **non-negative** weights;
weight types are `std::uint64_t` (exact) and `double` (IEEE-754).

## 1. Layout

| Component | Files | Namespace | Role |
|---|---|---|---|
| CSR graph | `include/sssp/graph.hpp` | `sssp` | `CsrGraph<W>`, `Edge`, `Arc`, `validate_edges`, `WeightTraits`. Shared by all algorithms; **not changed** by this library. |
| Formats | `include/sssp/graph_formats.hpp` (header-only) | `sssp` | `EdgeList<W>`, `AdjacencyList<W>`, conversions, `timed()` |
| Text IO | `include/sssp/graph_io.hpp`, `src/graph/io.cpp` | `sssp::io` | project edge-list format and DIMACS `.gr`, `ParseError` |
| Generators | `include/sssp/generators.hpp`, `src/graph/generators.cpp` | `sssp::gen` | PRNG, deterministic math, weight distributions, graph generators, catalog, `generate()` |
| CLI | `tools/sssp_gen.cpp` | — | `sssp_gen`: generate a graph and write it to a file |

Templates in `io.cpp` and `generators.cpp` are explicitly instantiated for `std::uint64_t` and `double`.

## 2. Representations and conversions

- `EdgeList<W> { std::size_t num_vertices; std::vector<Edge<W>> edges; }`: the format that IO and the generators exchange.
- `AdjacencyList<W>::from_edges(n, edges)`: a vector of vectors. It uses the same `Arc<W>`
  layout as `CsrGraph<W>` and the same accessors (`num_vertices`, `num_edges`, `num_arcs`,
  `degree`, `arcs(v)`, `edge(id)`, `edges()`). It validates its input the same way
  `CsrGraph::from_edges` does.
- Conversions are free functions with no caching:

  | Function | Input |
  |---|---|
  | `to_csr` | `const EdgeList&`, `EdgeList&&` (moves the edge vector), or `const AdjacencyList&` |
  | `to_edge_list` | `CsrGraph` or `AdjacencyList` |
  | `to_adjacency_list` | `EdgeList` or `CsrGraph` |

**Ordering contract.** Every representation uses the same edge ids: `EdgeId i` refers to the
same edge everywhere. `AdjacencyList::arcs(v)` lists arcs in exactly the order of
`CsrGraph::arcs(v)`: for each edge in id order, the u-side arc comes first and then the v-side
arc. A self-loop `{u,u}` is one arc. Because of this, an algorithm visits the arcs in the same
order on either layout.

**Timing conversions separately.** Benchmarks must not count conversion cost as algorithm time:

```cpp
auto [g, convert_seconds] = sssp::timed([&] { return sssp::to_csr(std::move(edge_list)); });
// run and time the algorithm on g; report convert_seconds as a separate column
```

`timed(f)` returns `Timed<T>{value, seconds}` and uses `std::chrono::steady_clock`.

## 3. File formats (`sssp::io`)

### Project edge-list format (`*.txt`)

```
n m
u v w      # exactly m lines, 0-based vertices, undirected edge {u,v}
```

- Blank lines, and lines whose first non-blank character is `#`, are ignored anywhere in the file.
- The file must contain exactly `m` edge lines. Too few or too many is an error.
- Self-loops and parallel edges are kept as written.

### DIMACS 9th Challenge `.gr`

```
c comment
p sp n m   # exactly once, before the first arc; m = number of ARC lines
a u v w    # 1-based, directed arc u->v
```

The graph is read as undirected in two steps.

1. `DimacsArcMode` decides how arcs become edges.
   - `kArcAsEdge` (default): each arc becomes one edge. A symmetric pair therefore becomes two
     parallel edges. Distances are unchanged, but `m` doubles.
   - `kMergeReversePairs`: arcs `u->v` and `v->u` with bit-identical weights are matched
     greedily in file order. A matched pair becomes **one** edge, placed at the position and in
     the orientation of the first arc of the pair. Arcs that find no partner (one-way arcs, or
     asymmetric weights) become one edge each. A self-loop arc is always one edge. Use this
     mode for the DIMACS road networks.
2. `DuplicateEdgePolicy` then handles parallel edges that remain. Two edges are parallel when
   they have the same unordered endpoint pair; two self-loops on the same vertex count as
   parallel.

   | Policy | Effect |
   |---|---|
   | `kKeepAll` (default) | keep every edge (multigraph) |
   | `kKeepMinWeight` | keep one edge, at its first position and orientation, with the minimum weight |
   | `kReject` | throw `ParseError` at the line of the first duplicate |

`write_dimacs` writes every edge as two arcs (`u->v`, then `v->u`) and a self-loop as one arc.
Reading the file back with `kMergeReversePairs` restores the original edge list exactly: same
order, same orientation, bit-identical weights.

### API

- Streams: `read_edge_list<W>(istream&)`, `write_edge_list(ostream&, el)`,
  `read_dimacs<W>(istream&, DimacsReadOptions{arc_mode, duplicates})`, `write_dimacs(ostream&, el)`.
- Files: `read_edge_list_file`, `write_edge_list_file`, `read_dimacs_file`, `write_dimacs_file`.
  They throw `std::runtime_error` if the file cannot be opened or written.

### Numbers

- Parsing and printing use `std::from_chars` / `std::to_chars`, so the result does not depend on the locale.
- `uint64_t` weights must be plain unsigned decimals in 0..2^64-1. A sign, a fraction or an exponent is rejected.
- `double` weights may be any finite non-negative decimal or exponent form. They are written
  with `chars_format::general` and precision 17 (the same as `%.17g`), which round-trips every
  finite double bit-exactly. `-0` is normalised to `+0`.

### Errors

All malformed input throws `io::ParseError` (derived from `std::runtime_error`). The message has
the form `line N: reason`, and `line()` returns the 1-based line number (0 means end of input or
the whole input). The errors detected are:

- negative, NaN, infinite or overflowing weights
- out-of-range vertices, including vertex 0 in DIMACS
- a missing, duplicate or malformed header
- too few or too many edge or arc lines
- non-numeric tokens, or the wrong number of tokens on a line
- `n >= kNoVertex` or `m >= kNoEdge`

The reader caps pre-reservation at 2^24 elements, so a huge `m` in a header cannot make it
allocate a huge buffer up front.

## 4. Reproducibility contract

Every generator is a pure function of (parameters, seed) and produces **bit-identical** output
on every platform and standard library. The rules below make that hold.

- **Own PRNGs.** We use `SplitMix64` for seeding and for deriving sub-stream seeds, and
  `Xoshiro256ss` (xoshiro256** 1.0) for the random streams. We do not use `<random>`: the C++
  standard leaves the algorithms of the `std::*_distribution` classes implementation-defined,
  so the same engine output gives different values on libstdc++, libc++ and MSVC. Both
  generators are checked against known-answer vectors from the reference implementations.
- **Exact integers.** `uniform_below(b)` rejects raw values below `2^64 mod b` and returns
  `r % b`. This has no modulo bias and no floating point. `uniform_range(lo, hi)` is inclusive
  and handles the full 64-bit range.
- **Exact reals.** `uniform01()` computes `(next() >> 11) * 2^-53`, a value in [0,1) with 53
  random bits and no rounding.
- **Deterministic math.** `det_log` and `det_exp` use only correctly rounded IEEE operations
  (`+ - * /`, `frexp`, `ldexp`, `floor`). libm's `log`, `exp` and `pow` are not correctly
  rounded and differ between platforms, so we do not use them. `det_log` uses an atanh series
  on a mantissa in [√½, √2). `det_exp` reduces the argument with a Cody–Waite split of ln 2
  (fdlibm `ln2_hi` / `ln2_lo`) and then sums a Taylor series. Both agree with libm to within
  about 1e-14 relative error. `sqrt` is correctly rounded by IEEE 754, so it is safe to use.
- **No FMA contraction.** A fused multiply-add rounds once instead of twice, so contracting
  expressions would change the results. `generators.cpp` contains
  `#pragma clang fp contract(off)`. GCC needs nothing, because in ISO mode (`-std=c++20`,
  `CMAKE_CXX_EXTENSIONS OFF`) it already uses `-ffp-contract=off`.
- **No unordered iteration.** `unordered_set` and `unordered_map` are used only for membership
  and counts. Output order never depends on their iteration order.
- **Independent sub-streams.** `substream_seed(seed, stream)` derives one stream for each
  purpose: topology = 0, weights = 1, connecting edges = 2. With the same seed, the topology is
  therefore identical whatever the weight distribution and whether `W` is `uint64_t` or `double`.

## 5. Weight distributions (`WeightSpec`)

`WeightSpec` has the fields `dist`, `lo` (default 1), `hi` (100), `mean` (10),
`zero_probability` (0.25) and `wide_log2` (50). `validate_weight_spec` throws
`std::invalid_argument` if they are inconsistent.

All distributions require `lo` and `hi` to be finite, `lo >= 0` and `hi <= 2^53`. For integer
`W`, real-valued draws are rounded to the nearest integer, with ties rounded up.

| `WeightDist` | CLI name | Draw | `uint64_t` | `double` | Extra constraint |
|---|---|---|---|---|---|
| `kUniformInt` | `uniform_int` | uniform over integers in [⌈lo⌉, ⌊hi⌋] | exact | integral doubles | the interval contains an integer |
| `kUniformReal` | `uniform_real` | `lo + (hi-lo)·U`, U ∈ [0,1) | rounded | real | lo ≤ hi |
| `kLogUniform` | `log_uniform` | `exp(uniform(log lo, log hi))`, clamped to [lo,hi] | rounded | real | 0 < lo ≤ hi |
| `kExponential` | `exponential` | `-mean·log(1-U)` | rounded | real | 0 < mean ≤ 2^47 |
| `kConstant` | `constant` | `lo` (use `lo = 0` for all-zero weights) | rounded | `lo` | — |
| `kZeroMixed` | `zero_mixed` | 0 with probability `zero_probability`, otherwise uniform in [lo,hi] | integer uniform | real uniform in [lo,hi) | zero_probability ∈ [0,1]; the interval contains an integer |
| `kWideRange` | `wide_range` | binary exponent e uniform in [0, K], then uniform within [2^e, 2^(e+1)); e = K gives exactly 2^K; covers [1, 2^K] | exact | `2^e·(1+U)` | K = `wide_log2` ∈ [0, 62] |

`to_string(WeightDist)` and `parse_weight_dist(name)` convert between the enum and the CLI
names; `parse_weight_dist` throws on an unknown name. `sample_weight<W>(spec, rng)` makes one draw.

## 6. Generators

Every generator returns `EdgeList<W>` and takes `GenOptions{weights, seed = 1, connected = false}`.

| CLI `--type` | Function | Parameters | Edges | Simple | Connected |
|---|---|---|---|---|---|
| `er` | `erdos_renyi(n, avg_degree, opt)` / `erdos_renyi_gnm(n, m, opt)` | m = round(n·avg_degree/2) ≤ n(n−1)/2 | m | yes (Floyd's sampling over pair indices) | only with `connected` |
| `regular` | `random_regular(n, d, policy, opt)` | n·d even; `kSwitchToSimple` needs d < n | n·d/2 (fewer with `kDrop`) | `kKeep`: no; `kDrop`: yes (degree ≤ d); `kSwitchToSimple`: yes (degree exactly d, by degree-preserving double-edge swaps) | only with `connected` |
| `ba` | `barabasi_albert(n, m, opt)` | 1 ≤ m < n | m(m+1)/2 + (n−m−1)·m | yes | yes |
| `grid` | `grid_2d(rows, cols, opt)`; vertex id = r·cols + c | rows, cols (default ⌊√n⌋ each) | rows(cols−1) + cols(rows−1) | yes | yes |
| `road` | `road_like(rows, cols, RoadLikeOptions, opt)` | jitter (default 0.3, < 0.5), keep_probability (0.9), shortcut_probability (0.15), shortcut_radius (2, ≥ 2), integer_scale (1000) | ≈ grid × keep_probability + shortcuts + restored streets | yes | yes |
| `tree` | `random_tree(n, opt)` | n | n−1 | yes | yes (uniform over labelled trees, from a random Prüfer sequence) |
| `chung_lu` | `chung_lu(n, avg_degree, exponent, opt)` | exponent > 2; m ≤ n(n−1)/4 | m = round(n·avg_degree/2) | yes (self-loops and duplicates are redrawn) | only with `connected` |
| `adversarial_decrease_key` | `adversarial_decrease_key(k)`, k = n−1 for the CLI | k | k + k(k−1)/2 | yes | yes |
| `adversarial_long_path` | `adversarial_long_path(n)` | n | (n−1) + (n−2) | yes | yes |

More detail on the less obvious generators:

- **Road-like.** Grid points are perturbed by U[−jitter, jitter] on each axis. Each street
  (grid edge) survives with probability `keep_probability`. With probability
  `shortcut_probability`, a vertex also gets one shortcut to a vertex within Chebyshev radius
  R that is not a grid neighbour. The graph is then made connected by restoring removed
  streets in Kruskal order, so every edge stays local. Weights are the Euclidean lengths of
  the perturbed points; integer weights are `max(1, round(len·integer_scale))`. The road
  generator ignores `opt.weights`.
- **Chung–Lu.** Vertex i has weight w_i ∝ (i+1)^(−1/(exponent−1)), computed with
  `det_exp`/`det_log`. Both endpoints of an edge are drawn proportionally to w through a
  prefix-sum binary search.
- **Barabási–Albert.** It starts from a clique on m+1 vertices. Each new vertex links to m
  distinct targets, sampled proportionally to degree from a list of edge endpoints.
- **Adversarial generators.** Their weights are fixed by the construction and do not use `opt.weights`.

### Adversarial constructions

- **Decrease-key** (k+1 vertices, source 0, v_i = i). The edges are `{0,1}` with weight 1,
  `{0,i}` with weight 2k for i > 1, and `{i,j}` with weight 2(j−i)−1 for 1 ≤ i < j.
  - The chain of weight-1 edges gives dist(v_j) ≤ j, and every edge `{i,j}` costs at least
    |j−i|, so dist(v_j) = j.
  - Dijkstra therefore settles the vertices in the order 1, 2, …, k.
  - Settling v_i lowers the key of **every** v_j with j > i, from 2j−i to 2j−i−1.
  - A heap-based Dijkstra therefore performs exactly k(k−1)/2 = Θ(m) decrease-keys. This is
    the worst case for binary heaps and the case Fibonacci heaps were designed for.
- **Long path** (n vertices). Path edges `{i,i+1}` have weight 1, and decoy skip edges
  `{i,i+2}` have weight 3 (more than the two path edges they bypass). The unique shortest path
  from 0 to n−1 has n−1 hops, and every decoy is relaxed and then rejected. This stresses
  methods whose cost grows with hop count or recursion depth.

### `connected` option and `connect_components`

`connect_components(el, weights, seed)` works as follows.

1. Components are ordered by their smallest vertex.
2. Each component i > 0 is linked to a uniformly random earlier component, by an edge between
   uniformly random members of the two. This forms a random tree over the components.
3. Exactly (number of components − 1) edges are appended. The original edges remain a prefix,
   unchanged.

The generators pass the connect sub-stream (2) as the seed. Generators that are connected by
construction skip this step. `count_components(el)` counts components, with each isolated
vertex counting as one.

### Catalog and generic entry point

- `generator_catalog()` returns `{name, description, parameters}` for each generator type.
- `generate<W>(GraphSpec)` dispatches on `GraphSpec::type`. `GraphSpec` holds `n`,
  `avg_degree`, `degree`, `attach`, `rows`, `cols`, `exponent`, `multi`, `road` and `options`.
  An unknown type throws `std::invalid_argument`.
- `parse_multi_edge_policy("keep"|"drop"|"switch")` converts the CLI name of the regular-graph policy.

## 7. Small test graphs

`template <class W> std::vector<NamedGraph<W>> small_graph_suite();` returns the boundary-case
graphs, where

```cpp
template <class W> struct NamedGraph { std::string name; EdgeList<W> graph; Vertex source = 0; std::string note; };
```

Sized builders are available in `sssp::gen::small`: `single_vertex()`, `no_edges(n)`,
`two_vertices(w)`, `path(n, w)`, `star(n, w)`, `cycle(n, w)`, `complete(n, w)` (edge {i,j}
has weight w·(1 + (i + 2j) mod 5)), `binary_tree(n, w)` and `many_equal_paths(k)`;
`long_shortest_path` is `adversarial_long_path(40)`. The remaining cases are fixed hand-built
edge lists that exist only inside the suite. The suite is header-only, so it works for any
weight type. The graphs are:

- **Structure:** `single_vertex`, `no_edges`, `two_vertices`, `path`, `star`, `cycle`,
  `complete`, `binary_tree`, `disconnected`, `multi_edges`, `self_loops`.
- **Weights:** `zero_weights`, `all_equal_weights`, `mixed_magnitudes`, `extreme_ratio`.
- **Shortest-path shape:** `long_shortest_path`, and `many_equal_paths`, a chain of k diamonds
  with 2^k equal shortest paths.

Differential and property tests should iterate over the whole suite instead of enumerating
the cases themselves.

## 8. CLI (`sssp_gen`)

```sh
sssp_gen --list                                   # generator types and weight distributions
sssp_gen --type er --n 1000 --avg-degree 4 --weights uniform_int --seed 1 --out g.txt
sssp_gen --type road --rows 100 --cols 100 --weight-type double --format dimacs --out road.gr
sssp_gen --type regular --n 1000 --degree 3 --multi switch --weights wide_range --wide-log2 50 --out r.txt
sssp_gen --type er --n 10000 --avg-degree 1 --connected --weights zero_mixed --zero-prob 0.3 --out z.txt
```

| Group | Options |
|---|---|
| Output | `--format edgelist` (default) or `dimacs`; `--weight-type uint64` (default) or `double`; `--out` |
| Weights | `--weights <dist>`, `--lo`, `--hi`, `--mean`, `--zero-prob`, `--wide-log2` |
| Shape | `--n`, `--avg-degree`, `--degree`, `--attach`, `--rows`, `--cols`, `--exponent`, `--multi keep\|drop\|switch` |
| Road | `--jitter`, `--keep-prob`, `--shortcut-prob`, `--shortcut-radius`, `--int-scale` |
| Other | `--connected`, `--seed`, `--stats` (n, m, components to stderr), `--list`, `--help` |

Exit codes: 0 on success, 2 for a usage error (unknown option, malformed number), 1 for a
generation or IO error (invalid parameters, unknown type, unwritable file). The same arguments
always produce byte-identical output.

## 9. Known limits

- Generators need n < 2^32, the 32-bit `Vertex` type. Pair keys pack two vertices into 64 bits.
- `chung_lu` needs m ≤ n(n−1)/4, so that redrawing duplicates terminates quickly. It throws
  `std::runtime_error` after 1000 + 100·m rejected draws.
- `random_regular` with `kSwitchToSimple` throws `std::invalid_argument` if the swaps do not
  converge within 1000 + 100·(n·d/2) proposals. In practice they converge quickly for d ≪ n.
- `DimacsArcMode::kMergeReversePairs` keeps a hash table of unmatched arcs, so its memory grows
  in proportion to the number of arcs on large road networks.
- `det_log` and `det_exp` are accurate to a few ulp, not correctly rounded. They are exactly
  reproducible across platforms, but they are not bit-equal to libm.
- Integer weight parameters must stay at or below 2^53 so that their conversion through
  `double` is exact.
