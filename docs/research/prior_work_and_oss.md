# Survey of Prior Work and Existing OSS

**English** | [日本語](prior_work_and_oss.ja.md) | [简体中文](prior_work_and_oss.zh-CN.md)

- Survey date: 2026-11-27
- Target paper: Iman Gholami, MohammadTaghi Hajiaghayi, "Linear Single-Source Shortest Path, Almost", arXiv:2611.34911v1 (submitted 2026-11-25) <https://arxiv.org/abs/2611.34911>
  - Problem setting: undirected multigraphs, nonnegative real edge lengths, exact distances from a single source to every vertex plus a shortest-path tree. Comparison-addition model. Expected time O(m + n log log log n). It terminates almost surely, and every terminating run is correct (paper Theorem 1, Section 1 of the main text).
  - The paper's arXiv license is `arXiv.org perpetual, non-exclusive license 1.0` (license link on the abs page <http://arxiv.org/licenses/nonexclusive-distrib/1.0/>). This license only grants distribution rights to arXiv, so the paper's text, figures and translations are not redistributed in this repository. `.cache/` is already in `.gitignore`.
- Policy: facts carry a source. Anything that could not be confirmed is marked "unverified". Only numbers that appear in a source are written; no speculative numbers.
- The cloned code was only read, not executed. Clones were made into a local temporary directory (`git clone --depth 1`).

---

## 0. Key points with direct impact on this project

1. **On any physically possible input, the paper's main algorithm (Algorithm 8 `ShortestPaths`) always enters the base case.**
   - Line 4 of Algorithm 8 says "if `N ≤ L_N^64`, retain a lightest edge per vertex pair and solve with the scan-based Dijkstra search (Lemma 3.2)". Here L_N = ⌈log2(N+2)⌉ (definition in Section 2 of the paper text and Algorithm 8; lines 283–290 and 2167–2185 of `.cache/paper/2611.34911v1.txt`).
   - The main text of Section 2 also says "The outer and local constructions assume N > L_N^64 and N > 2 L_N^4, respectively. The excluded sizes are absolutely bounded."
   - Computing the smallest L_N with N > L_N^64 gives 589. That is, the nontrivial path is entered only when N ≳ 2^589 (this computation was done in this survey: a Python scan over L = 2..2000 for the first L with 2^L − 2 > L^64). Exact values are in [thresholds.md](thresholds.md): the smallest N is 589^64 + 1 ≈ 1.94×10^177, the test is false again on [2^589 − 1, 590^64], and it holds for every N ≥ 590^64 + 1.
   - Consequently, an implementation with the paper's thresholds always measures Dijkstra at any measurable scale. To actually run and measure the hierarchical construction, a parameter that shrinks the threshold is needed (e.g. `--sep k` in the kingfree implementation described below). This is outside the scope of the paper's theorem, so benchmarks must report "paper thresholds" and "shrunk-threshold variant" separately.
2. **There is already one public OSS implementation of the target paper itself** (kingfree/sssp-almost-linear, created 2026-11-26). However, its Evaluate stage is incomplete and is not used by default (§1.4).
3. The three requested repositories (mizchi, lcs147, mjkozicki) all implement the **directed** DMM+25 (arXiv:2504.17033), whose problem setting (directed vs. undirected) differs from the target paper.

---

## 1. Existing OSS

### 1.1 Overview

| Repository | Exists | License (checked against the actual text) | Latest commit SHA | Date (committer, UTC) | Main language | Target paper | Target problem |
|---|---|---|---|---|---|---|---|
| [mizchi/check-arxiv-2504-17033](https://github.com/mizchi/check-arxiv-2504-17033) | ○ | **No LICENSE file**. Only "## License MIT" at the end of the README; the GitHub API `license` is also `null`. `Cargo.toml` has no `license` field | `b3f173ca43d4052c0872d824be0bd61456ab7dd3` | 2026-09-23T17:08:40Z | Rust | DMM+25 (arXiv:2504.17033), C-HD (vals.ai) | Directed, nonnegative weights (experiments hold integers in [1, 2^20] as f64), single source to all vertices |
| [lcs147/bmssp](https://github.com/lcs147/bmssp) | ○ | MIT (`LICENSE`: "Copyright (c) 2025 Lucas Castro and Thailsson Clementino") | `0a8eb82805bad7f64dfa1b51e17656c9af4c4681` | 2026-01-24T15:47:14Z | C++20 header-only (by GitHub's statistics, Jupyter Notebook is the largest) | DMM+25 (STOC 2025). Experimental paper: Castro, Clementino, de Freitas, arXiv:2511.03007 | Directed, nonnegative weights (integer or floating point selectable via template), single source to all vertices |
| [mjkozicki/SSSP](https://github.com/mjkozicki/SSSP) | ○ | MIT (`LICENSE`: "Copyright (c) 2026 Mike Kozicki") | `30231195ef664b7ba63bd9f1190c33ecd2735575` | 2026-10-08T00:35:22Z | Python and others, 8 languages (C#, Rust, C++, Go, Java, PHP, Python, TypeScript) | DMM+25 (arXiv:2504.17033v2) | Directed, nonnegative reals (double), single source to all vertices |
| [kingfree/sssp-almost-linear](https://github.com/kingfree/sssp-almost-linear) (found additionally) | ○ | MIT (`LICENSE`: "Copyright (c) 2026 kingfree"). The README says "the C code is MIT; rights to the Chinese translation belong to the original authors" | `151ed68a9771b4211edae7700a07bbe22c5b0c38` | 2026-11-26T10:26:16Z | Typst (translation), C11 (implementation) | **Target paper arXiv:2611.34911v1** | Undirected, nonnegative reals (IEEE double treated as an "opaque register"), single source to all vertices |

Verification method: `gh api repos/<o>/<r>` and `gh api repos/<o>/<r>/commits?per_page=1`, `git rev-parse HEAD` of each clone (matched the API SHA), and the text of `LICENSE`.

### 1.2 mizchi/check-arxiv-2504-17033

- Sources: README <https://github.com/mizchi/check-arxiv-2504-17033/blob/b3f173ca43d4052c0872d824be0bd61456ab7dd3/README.md>, replication report `docs/chd-replication.md` (same SHA)
- Implementation: Dijkstra (`dijkstra_csr` with a binary heap and lazy deletion, `dijkstra_dary` with an indexed 4-ary heap), `src/dmm25.rs`, a reimplementation of DMM+25 faithful to the paper (Algorithms 1–3, the Lemma 3.3 block structure, constant-degree transformation), and C-HD in `src/chd.rs`. The older 2025 implementations V1, V2 and core also remain.
- Regarding the older implementations, the repository itself acknowledged errors in 2026-09: V1 is identical code to Dijkstra, the Bellman–Ford phase of V2 is unreachable, and core falls back to a full Dijkstra at the top-level call. The old README's "1.1–1.6x faster" was measurement noise (single measurement, ∞ vertices were not compared).
- Correctness verification: `tests/bmssp_family_test.rs`. All vertices were compared against Dijkstra on 3000 random directed multigraphs with n < 40 (including zero weights, ties, self-loops, parallel edges and unreachable vertices), and on medium-sized random graphs, grids and d-regular graphs. They report 0 mismatches also when comparing across all benchmark cases (`docs/chd-replication.md` §3).
- Benchmark conditions (same document §4): 4-vCPU cloud container, Rust 1.94 `--release`, single-threaded, median of 3 runs. Random directed graphs with d out-edges per vertex (uniformly random heads), d ∈ {2, ⌈√lg n⌉, ⌊lg^{3/4} n⌋, 16}. Grids are 4-neighbor and bidirectional. Source is 0. Timing wraps only each algorithm call with `Instant` (`src/bin/compare.rs` lines 88–96). CSR construction is outside the timing.
- Results (slowdown relative to the fastest Dijkstra, README): DMM+25 is 14–121x slower (38–50x on grids). C-HD is 1.47–7.3x slower including preprocessing. n = 2^14–2^22.
- Fairness: 
  - (Good) C-HD is reported both with and without preprocessing. The Dijkstra baseline is the fastest of three variants.
  - (Caution) CSR construction for Dijkstra is outside the timing, but the constant-degree transformation for DMM+25 is inside it (the function takes a `Graph`). This asymmetry works against the algorithm.
  - The comparison is on directed graphs, so its numbers cannot be directly compared with this project's (undirected).

### 1.3 lcs147/bmssp

- Sources: README <https://github.com/lcs147/bmssp/blob/0a8eb82805bad7f64dfa1b51e17656c9af4c4681/README.md>, arXiv:2511.03007 <https://arxiv.org/abs/2511.03007>
- Implementation: `include/bmssp.hpp` (`spp::bmssp`, worst-case version) and `include/bmssp-expected.hpp` (`spp_expected::bmssp`, expected-time version with lower memory). `single_include/bmssp.hpp` is generated in CI (quom). Dependencies are bundled in the repository (`unordered_dense.h`, `floyd_rivest_select.h`, etc.).
- Abstract of arXiv:2511.03007: the worst-case version uses O(n log^{1/3} n) space and the expected-time version O(n) space, with the same time complexity. On sparse random graphs, grids and US road networks (up to 10 million vertices), Dijkstra is 3–4x faster in all settings. They estimate that the worst-case version would need far more than 10^67 vertices to beat Dijkstra.
- Correctness verification: doctest (`tests/tests-wc.cpp`, `tests/tests-expected.cpp`). Distances are compared against the bundled Dijkstra (`helpers/dijkstra.hpp`) on grids (5, 25, 100) and random degree-3 graphs (32–4096 vertices), and the endpoints of paths reconstructed from predecessors are also checked (`tests/common.cpp`). Both with and without the constant-degree transformation are tested. CI (`.github/workflows/ci.yaml`) runs `./tests`.
- Benchmark conditions: the README table is "random graphs, 3n edges"; BMSSP is 3.0–4.2x slower than Dijkstra (n = 2^12–2^25). The README does not state the machine or the number of repetitions (unverified; the main text of arXiv:2511.03007 was not read this time). In `experiments/experiment.cpp`, timing wraps only the `execute(s)` call. `prepare_graph(false)` (graph preparation) and the Dijkstra constructor are outside the timing.
- Fairness:
  - Preparation (`prepare_graph`) is outside the timing. However, the experiments call it with `false` (no constant-degree transformation), and since the graphs have degree 3, the effect is probably small (unverified). Whether to include the cost of the constant-degree transformation in the timing must be decided again when vendoring.
  - The distance type is integer weights as `long long` (`distT` in `experiment.cpp`). This differs from the conditions of the real-weight comparison-addition model.
  - For directed graphs. If used in this project, undirected graphs would be input as arcs in both directions.
  - It bundles `tests/graphs/USA-road-t.NY.gr` (a DIMACS road network). The redistribution terms of DIMACS data are unverified, so when vendoring, either exclude the data or check the terms.

### 1.4 kingfree/sssp-almost-linear (implementation of the target paper)

- Sources: README <https://github.com/kingfree/sssp-almost-linear/blob/151ed68a9771b4211edae7700a07bbe22c5b0c38/README.md>, `c/README.md` (same SHA)
- Created 2026-11-26, 2 commits, 0 stars (GitHub API).
- Contents: a Chinese translation of the full paper in `typst/` (Typst and PDF), and a C11 implementation in `c/` (about 4000 lines).
- Scope of the implementation (self-reported in `c/README.md`):
  - Simulates the comparison-addition model in `len.h` and counts operations.
  - MST (3-pass bounded growth), parallel merging, LCA (O(log n) doubling instead of the paper's O(1) structure), owner decomposition, route certificates (equivalent to Algorithms 2 and 3), owner coordinate J, outer state graph, threshold components, hierarchy assembly.
- **Known deviations (self-reported)**:
  1. The slot scan of the Evaluate stage (`--slot`) misses values or outputs values that are too large on some inputs, so it is **not used by default**. The default path is "exact search on the constructed state graph" (`eval.c: state_graph_dist`). This verifies whether the §8.1 construction realizes geodesics with their original lengths; it is not an implementation of the paper's evaluation procedure.
  2. The §7.3 width calendars are not implemented. J is computed by Dijkstra on the auxiliary graph.
  3. The §5.4 atom substitution is not implemented.
  4. The O(1) LCA and interval minima are substituted with O(log n).
- Threshold: by default it follows the §0 threshold N ≤ L_N^64, so it always enters the base case (scan-based Dijkstra). Shrinking the separator size to b = L_N^k with `--sep k` makes the hierarchy path run.
- Correctness verification: `--verify` compares each vertex against a heap-based Dijkstra. The README reports "12/12 match on random graphs and path graphs with n = 20–100, and n = 800 also matches". `hier_validate` checks the interval invariants of the hierarchy. No automated test suite (CI or test runner) was found.
- Benchmark: `--bench` outputs the "generation + hierarchy" and "evaluation" times for 5 cases (n = 500–8100). The README has no numbers.
- Fairness and open points (found by reading the code in this survey; not executed):
  - In `--bench` at `c/main.c` lines 323–336, `gen_rand` is called when `name[0]=='g' && name[1]=='r'` (i.e. "grid-…"), and `gen_grid` is called for "rand-…". **The case names and graph types appear to be swapped.**
  - The default path is not the paper's evaluation procedure (item 1 above), so this implementation's times are not the times of the paper's algorithm.
  - The README and c/README call the main algorithm "算法 9" (Algorithm 9) and the evaluation "算法 7、8" (Algorithms 7 and 8). In the arXiv v1 text, however, the main algorithm is Algorithm 8 (`ShortestPaths`) and the evaluation is Algorithm 3 (`Evaluate`). Which version's numbering they refer to is unverified.
  - The Chinese translation is a translation of a paper under arXiv's non-exclusive license, so permission to redistribute it is unverified. It is not vendored into this project.

### 1.5 mjkozicki/SSSP

- Sources: README <https://github.com/mjkozicki/SSSP/blob/30231195ef664b7ba63bd9f1190c33ecd2735575/README.md>, `benchmark/README.md` (same SHA)
- Implementation: each of the 8 languages has Dijkstra and `duan_mao_shu_yin`. The README itself says "the frontier data structure is a simplified version that preserves correctness; using the Lemma 3.3 block structure would improve constants".
  - The Rust `FrontierDS` (`rust/src/sssp.rs` lines 342–411) removes linearly with `Vec::retain` on every `insert` and sorts everything on every `pull`. There is no guarantee of the Lemma 3.3 complexity, and it appears not to achieve O(m log^{2/3} n) (this survey's reading; not measured).
- Correctness verification: each language has tests comparing against Dijkstra, for example Rust's `matches_dijkstra_small` and `matches_dijkstra_parallel_paths` (tolerance 1e-10) and Python's `test_matches_dijkstra`. However, they focus on small graphs.
- Benchmark conditions (`benchmark/README.md`): 50,000 vertices, a directed graph imitating a network of about 20 regions (2–8 intra-region out-edges with weights 0.5–5, 0–2 inter-region out-edges with weights 20–150, seed = 42). A Docker container is started for each language, and `wall_sec` is **the time from `docker run` to `docker wait`** (`benchmark/run_benchmarks.py` lines 209–233). The default is a single measurement.
- Results (`benchmark/results.json`): in most languages both DMSY and Dijkstra take about 1.9–2.0 seconds, and `peak_mem_mb` and `total_cpu_sec` are 0.0.
- Fairness:
  - The measurements include container startup and reading the graph file. The algorithmic difference is barely visible, so they cannot be used to compare algorithms.
  - Single measurement. Memory and CPU collection is 0 in many languages, so the measurement appears not to work.
  - `rust/src/bin/compare.rs` contains a separate comparison that takes the median over repetitions, but its results are not in the repository (unverified).

---

## 2. Prior work (focusing on 2023–2026)

"Match" means whether the setting matches this project's problem setting (undirected, nonnegative reals, single source to all vertices, comparison-addition model).

| Abbrev. | Reference | Problem setting | Complexity | Model | Match | Use as a comparison target? |
|---|---|---|---|---|---|---|
| Dij59 / FT87 | Dijkstra 1959, Fredman–Tarjan 1987 (Fibonacci heap). DOI 10.1007/BF01386390, 10.1145/28869.28874 | Directed and undirected, nonnegative reals | O(m + n log n). O((m+n) log n) with a binary heap | Comparison-addition | ○ | **Required (baseline)**. Implemented in-house (binary heap, 4-ary heap, pairing heap if possible) |
| Tho99 | Thorup, "Undirected Single-Source Shortest Paths with Positive Integer Weights in Linear Time", J. ACM 46(3), 1999. DOI 10.1145/316542.316548 | Undirected, **positive integer** weights | O(m) (worst case) | word RAM (uses the bit representation of weights) | △ (different weights and model) | Optional. Measured only as a reference on the integer-weight subset. Not included in real-weight comparisons. Public implementation unverified |
| PR05 | Pettie, Ramachandran, "A Shortest Path Algorithm for Real-Weighted Undirected Graphs", SIAM J. Comput. 34(6), 2005. DOI 10.1137/S0097539702419650 | Undirected, nonnegative reals | Preprocessing O(MST(m,n) + min{n log n, n log log r}), where r is the ratio of weights | Comparison-addition | ○ | Theoretical comparison target. Public implementation unverified. Implementing it in-house is a large undertaking |
| DMS+23 | Duan, Mao, Shu, Yin, "A Randomized Algorithm for Single-Source Shortest Path on Undirected Real-Weighted Graphs", FOCS 2023. arXiv:2307.04139 <https://arxiv.org/abs/2307.04139> | Undirected, nonnegative reals | Randomized, O(m √(log n · log log n)) | Comparison-addition | ○ | **Strong candidate**. The problem setting matches exactly. The abstract says it is "simple in structure and easy to implement". No public implementation was found in this search (`gh search repos`) |
| Yan26 | Shuyi Yan, "Lossless Derandomization for Undirected Single-Source Shortest Paths and Approximate Distance Oracles", SOSA 2026. arXiv:2510.12598 <https://arxiv.org/abs/2510.12598> | Undirected, nonnegative reals | Deterministic, O(m √(log n log log n)) | Comparison-addition | ○ | Deterministic version of DMS+23. If DMS+23 is implemented, add it if capacity allows |
| DMM+25 | Duan, Mao, Mao, Shu, Yin, "Breaking the Sorting Barrier for Directed Single-Source Shortest Paths", STOC 2025. arXiv:2504.17033 <https://arxiv.org/abs/2504.17033> | **Directed**, nonnegative reals | Deterministic, O(m log^{2/3} n) | Comparison-addition | △ (directed) | Secondary comparison target. Undirected graphs can be solved as arcs in both directions, but undirectedness is not exploited. lcs147/bmssp can be used (§3) |
| DMS+26 | Duan, Mao, Shu, Yin, "A Faster Directed Single-Source Shortest Path Algorithm", ICALP 2026 (LIPIcs 374, 81:1–81:23). arXiv:2602.07868 <https://arxiv.org/abs/2602.07868> | Directed, nonnegative reals | Deterministic, O(m √(log n) + √(mn log n log log n)). O(m √(log n log log n)) on sparse graphs | Comparison-addition | △ (directed) | Reference. Public implementation unverified |
| GZ26 | Grishin, Zabokritskiy, "Improved Bounds for Directed SSSP in Sparse and Intermediate-Density Graphs", arXiv:2611.26453 <https://arxiv.org/abs/2611.26453> (submitted 2026-11-19) | Directed, nonnegative reals | Deterministic, O(n √(log n) (log log n)^{1/4}) for m = O(n) | Comparison-addition | △ (directed) | Reference. Preprint |
| KR26 | Kadria, Roditty, "A Faster Undirected Single-Source Shortest Path Algorithm", arXiv:2609.15247 <https://arxiv.org/abs/2609.15247> (submitted 2026-09-14) | Undirected, nonnegative reals | O(m log^{1/2} n · (log log n)^{1/4} · (log log log n)^{1/4}). The target paper describes it as "with high probability" | Comparison-addition | ○ | Direct theoretical competitor (the best undirected result before the target paper). Public implementation unverified. Preprint |
| HHR+24 | Haeupler, Hladík, Rozhoň, Tarjan, Tětek, "Universal Optimality of Dijkstra via Beyond-Worst-Case Heaps", FOCS 2024 (pp. 2099–2130, Princeton bibliographic record <https://collaborate.princeton.edu/en/publications/universal-optimality-of-dijkstra-via-beyond-worst-case-heaps/>). arXiv:2311.11793 <https://arxiv.org/abs/2311.11793> | Directed and undirected, nonnegative reals. Mainly the "distance ordering problem" (ordering vertices by distance) | Dijkstra with a heap having the working-set property is universally optimal | Comparison-addition | △ (a problem that outputs the order, not only the distance values) | Important for interpretation. Closely related to the opening of Section 1 of the target paper (computing distances is different from sorting them). Public implementation of the working-set heap unverified. The FOCS 2024 Best Paper Award could only be confirmed through a secondary source (arXiv:2504.17327) |
| Castro+25 | Castro, Clementino, de Freitas, "Implementation and Brief Experimental Analysis of the Duan et al. (2025) Algorithm", arXiv:2511.03007 <https://arxiv.org/abs/2511.03007> | Implementation and experiments of DMM+25 | — | — | — | Reference for benchmark methodology (road networks, grids, random) |

Notes:

- The references cited in Section 1.1 (Related Work) of the target paper are Dij59, FT87, Tho99, PR05, DMS+23, Yan26, DMM+25, DMS+26, GZ26 and KR26 from the table above, plus FRT04 (the tree embedding of Fakcharoenphol, Rao and Talwar; DOI 10.1016/J.JCSS.2004.04.011). FRT04 is a component (the route estimator), so it is not a comparison target. HHR+24 is not in the target paper's references (added on request).
- The target paper's text notes that its bounds are stated "on connected undirected inputs or on the reachable subgraph of a directed input". Identifying the reachable vertices and writing all n output labels adds at most O(m + n) (Section 1.1). When comparing, the handling of unreachable vertices (writing ∞) must be made consistent across methods.
- The target paper's time complexity is an **expected** bound (randomized). DMS+23 is also randomized, while Yan26, DMM+25, DMS+26 and GZ26 are deterministic. Randomized methods must be run with multiple seeds and reported as distributions.

---

## 3. Vendoring existing OSS vs. reimplementing

### 3.1 Recommendation

**Implement the target paper and the baseline Dijkstra in-house. Do not vendor code from existing OSS.** The sole exception: only when a secondary comparison with directed BMSSP is needed, pin lcs147/bmssp as an external dependency.

| Repository | Recommendation | Reason |
|---|---|---|
| kingfree/sssp-almost-linear | **Do not vendor**. Use only as a read-only cross-reference | It is MIT, so vendoring is legally possible. However, the Evaluate stage is incomplete and not used by default, the benchmark case names and graphs appear to be swapped, there are no automated tests, and it has 2 commits one day after creation. It cannot serve as a reference for a correct reproduction of the paper. Permission to redistribute the translation is unverified. Once this project's implementation is complete, it may be worth using for "cross-checking independent implementations" by comparing distances on the same inputs (whether to run it is decided separately) |
| lcs147/bmssp | **Only if a secondary comparison is needed**, pin it by SHA as an external dependency | MIT, header-only C++20, with tests against Dijkstra and CI. This project is C++/CMake, so it is easy to integrate. However, its experiments assume directed graphs and integer weights, so undirected graphs must be given as arcs in both directions, and its behavior with the weight type set to `double` must be checked with our own tests before use |
| mizchi/check-arxiv-2504-17033 | **Do not vendor**. Only cite its numbers and methodology as literature | No LICENSE file (the README only says "MIT"). MIT requires including the copyright notice and permission notice, but that text is absent. It is in Rust, which also differs from this project's language. Its way of counting comparisons (`count-cmp`) and benchmark methodology are useful references |
| mjkozicki/SSSP | **Do not vendor** | The frontier is a simplified version and does not meet the DMM+25 complexity. The benchmark is a single measurement including container startup and is not an algorithm comparison |

### 3.2 How to pin the commit when vendoring (for lcs147/bmssp)

- **Write the 40-character SHA, not a tag or branch**, as in CMake's `FetchContent_Declare(bmssp GIT_REPOSITORY https://github.com/lcs147/bmssp.git GIT_TAG 0a8eb82805bad7f64dfa1b51e17656c9af4c4681)`. A git submodule is also fine (the SHA is recorded in `.gitmodules` and the tree).
- Use only `include/` (including `external/`); do not use the data in `tests/graphs/` (including DIMACS-derived data) or `experiments/`.
- As required by MIT, include its `LICENSE` in distributions and record the name, SHA and license in `THIRD_PARTY_NOTICES` or similar. The licenses of the bundled `include/external/*` (`unordered_dense.h`, `floyd_rivest_select.h`, `median_of_ninthers.h`, etc.) must be checked in each file's header before vendoring (unverified in this survey).
- This project itself does not yet have a license in the repository (no `LICENSE` found). Compatibility with MIT code will be finally confirmed once this project's license is decided (MIT is generally compatible with almost all other licenses).

### 3.3 Priority of comparison targets (proposal)

1. Dijkstra (in-house implementation; binary heap and 4-ary heap; also counts comparisons). As described in §0, an implementation with the target paper's thresholds enters the Dijkstra base case, so this effectively becomes the main comparison.
2. Implementation of the target paper (in-house). Report "paper thresholds" and "shrunk-threshold variant" separately.
3. DMS+23 (if implemented in-house). The only prior non-Dijkstra method that both matches the problem setting exactly and is easy to implement.
4. DMM+25 (using pinned lcs147/bmssp). State clearly that these are reference values from applying a directed method to undirected inputs.

---

## 4. Unverified items

- The measurement environment and number of repetitions in the main text of arXiv:2511.03007 (only the abstract was read).
- The machine and conditions on which the README table of lcs147/bmssp was measured.
- The licenses of `include/external/*` bundled with lcs147/bmssp.
- The redistribution terms of the DIMACS road network data.
- Whether public implementations of Tho99, PR05, DMS+23, Yan26, DMS+26, KR26, GZ26 and HHR+24 exist (searching paper IDs and titles with `gh search repos` found none, but the search was not exhaustive).
- Whether the swapped case names in the kingfree implementation's `--bench` are intentional (the code was only read, not executed).
- Which version of the paper the kingfree implementation's algorithm numbers correspond to.
- The FOCS 2024 Best Paper Award for HHR+24 (could only be confirmed through a secondary source).
