# almost-linear-sssp

**English** | [日本語](README.ja.md) | [简体中文](README.zh-CN.md)

An independent reproduction, verification and benchmark project for the paper
**"Linear Single-Source Shortest Path, Almost"** by Iman Gholami and MohammadTaghi Hajiaghayi
([arXiv:2611.34911v1](https://arxiv.org/abs/2611.34911)).

We implement the paper in C++20, check every result against independent oracles, and plan to
measure it against well-tuned Dijkstra baselines. This project is not affiliated with the
authors.

> **Status: early work in progress.** The graph library, the Dijkstra baselines and the
> verification tools are done. The paper's own components are being implemented.
> **Nothing has been benchmarked yet**, so this README contains no performance numbers.

## Contents

- [Status](#status)
- [The paper in one paragraph](#the-paper-in-one-paragraph)
- [Why the thresholds matter](#why-the-thresholds-matter)
- [Quick start](#quick-start)
- [Algorithms](#algorithms)
- [How correctness is checked](#how-correctness-is-checked)
- [Benchmarks](#benchmarks)
- [Repository layout](#repository-layout)
- [Documentation](#documentation)
- [Related work](#related-work)
- [Citing](#citing)
- [License](#license)
- [Disclaimer](#disclaimer)

## Status

| Area | State |
|---|---|
| Graph library: CSR, adjacency list, edge list; edge-list and DIMACS `.gr` IO | done |
| Reproducible graph generators and the `sssp_gen` command-line tool | done |
| Baselines: Bellman–Ford; Dijkstra with binary, 4-ary, Fibonacci and radix heaps | done |
| Verification: shortest-path-tree checker, exact rational reference, `Counted<W>` operation counter | done |
| Paper components, Sections 2–5 (primitives, hierarchy evaluation, local hierarchy) | in progress |
| Paper components, Sections 6–8 (pivot sampling, owner coordinates, outer hierarchy) | planned |
| Differential test harness and benchmark runner | in progress |
| Benchmark results | **not measured** |

The detailed task list is [docs/TASKS.md](docs/TASKS.md); the current summary is
[docs/PROGRESS.md](docs/PROGRESS.md).

## The paper in one paragraph

Given an undirected graph with `n` vertices, `m` edges and non-negative real edge lengths, the
single-source shortest path (SSSP) problem asks for the distance from one source vertex to every
other vertex. Dijkstra's algorithm with a Fibonacci heap takes `O(m + n log n)` time. The paper
claims an expected running time of **`O(m + n log log log n)`** in the *comparison–addition
model*, where edge lengths may only be compared and added (no rounding, no bit tricks). That is
linear whenever `m` is at least `n log log log n`. The algorithm is randomised but *Las Vegas*:
every run that terminates returns exact distances, and randomness only affects the running time.
It combines a minimum-spanning-tree based decomposition into regions with one "pivot" each,
random sampling to certify short routes between pivots, auxiliary per-vertex coordinates, and a
procedure that evaluates distances along an ordered hierarchy of strongly connected components.

## Why the thresholds matter

This is the most important finding so far, and it shapes the whole project.

The paper uses its full construction only on inputs that are large enough. With
`L_N = ⌈log2(N + 2)⌉` for an input with `N` vertices (Section 2), it requires:

- **outer construction:** `N > L_N^64`; smaller inputs go to a base case that runs Dijkstra;
- **local construction** (run on small pieces): `N > 2 · L_N^4`; smaller pieces go to a simple
  Kruskal-style base case.

We computed the exact boundaries with arbitrary-precision integers
([docs/research/thresholds.md](docs/research/thresholds.md)):

| Condition | First `N` where it holds |
|---|---|
| outer, `N > L_N^64` | `589^64 + 1` ≈ 1.9 × 10^177 (it fails again on a short interval just above; it holds for every `N ≥ 590^64 + 1` ≈ 2.2 × 10^177) |
| local, `N > 2 · L_N^4` | 209,953 (`2 · 18^4 = 209,952`) |

| `N` | `L_N` | outer construction runs? | local condition holds? |
|---|---|---|---|
| 10^3 | 10 | no (`L_N^64` ≈ 2^213) | no (needs `N > 20,000`) |
| 10^5 | 17 | no (≈ 2^262) | no (needs `N > 167,042`) |
| 10^6 | 20 | no (≈ 2^277) | yes (`N > 320,000`) |
| 10^8 | 27 | no (≈ 2^304) | yes (`N > 1,062,882`) |
| any `N < 2^64` | ≤ 65 | **no**, because `L_N^64 ≥ 2^64 > N` | — |

So a literal implementation of the paper runs **only its base case (Dijkstra) on every input
that fits in memory**. The local construction is invoked only from inside the outer one, so it
never runs either. Timing a literal implementation measures Dijkstra plus some preprocessing,
not the paper's new ideas.

To study the construction honestly, the project keeps three modes apart
([ADR 0003](docs/adr/0003-three-implementation-modes.md)):

| Mode | Directory | What it is | What it is for |
|---|---|---|---|
| **A** paper-faithful | `src/paper/` | Thresholds exactly as written, checked with overflow-safe integer arithmetic. Every entry into a base case is recorded. | Reproducing the paper as written. On real inputs this is the base case; the variant `paper-faithful-d1` swaps the scan-based base-case Dijkstra for a binary-heap one (deviation D-1, [ADR 0008](docs/adr/0008-mode-a-outer-base-case-variants.md)). |
| **B** paper-experimental | `src/experimental/paper_experimental/` | The paper's structure with smaller, configurable size thresholds. Constants that the correctness proofs depend on are kept, or re-derived from the new sizes ([thresholds.md §6](docs/research/thresholds.md)). | Actually running Sections 5–8 on realistic sizes, with every result verified. An input whose correctness cannot be established counts as a failure. |
| **C** practical-optimized | `src/experimental/practical/` | Any optimisation, as long as it is correct and labelled. | Practical speed. |

Only Mode A is ever called "the paper's algorithm". Every reported timing will show the share of
vertices and edges handled by a base case next to it.

## Quick start

Requirements: macOS or Linux, a C++20 compiler (developed and tested with Apple Clang and LLVM
Clang), CMake ≥ 3.24, Ninja, and [uv](https://docs.astral.sh/uv/) for Python ≥ 3.12. CI and
contributors also need clang-format and clang-tidy.

```bash
make bootstrap                 # check the toolchain, install what is missing, run `uv sync`, configure the debug preset
                               # (./scripts/bootstrap.sh --check only reports, installs nothing)
make build                     # configure and build the debug preset into build/debug/
make test                      # run the C++ tests (ctest) and the Python tests (pytest)

./build/debug/sssp_gen --list  # list graph generators and weight distributions
./build/debug/sssp_gen --type er --n 1000 --avg-degree 4 --weights uniform_int --seed 1 --out g.txt

make fetch-paper               # download the paper into .cache/paper/ (never committed)
```

The same seed always produces the same graph on every platform. Other build presets (`release`,
`sanitizer`, `coverage`, `ci`, …) are listed in [CMakePresets.json](CMakePresets.json); for
example `make build PRESET=release`. Before contributing, run `make verify`
(see [CONTRIBUTING.md](CONTRIBUTING.md)).

## Algorithms

| Algorithm | Function / location | Weights | State |
|---|---|---|---|
| Bellman–Ford (round-based; used as an oracle) | `sssp::baseline::bellman_ford` | `uint64_t`, `double` | done |
| Dijkstra, binary heap (lazy deletion) | `sssp::baseline::dijkstra_binary_heap` | `uint64_t`, `double` | done |
| Dijkstra, 4-ary heap (decrease-key) | `sssp::baseline::dijkstra_4ary_heap` | `uint64_t`, `double` | done |
| Dijkstra, Fibonacci heap | `sssp::baseline::dijkstra_fibonacci_heap` | `uint64_t`, `double` | done |
| Dijkstra, radix heap | `sssp::baseline::dijkstra_radix_heap` | `uint64_t` only | done |
| Paper, Sec 3: primitives (MST, ordered union–find, LCA / level ancestor, split–findmin) | `src/paper/` (Mode A) | — | in progress |
| Paper, Sec 4: evaluation of an ordered hierarchy | `src/paper/` (Mode A) | — | in progress |
| Paper, Sec 5: local hierarchy | `src/paper/` (Mode A) | — | in progress |
| Paper, Sec 6–8: pivot routes, owner coordinates, outer hierarchy, main algorithm | `src/paper/` (Mode A) | — | planned |
| Tuned thresholds / practical optimisations | Mode B / Mode C | — | planned |

All baselines share one scan loop, use heaps written in this repository and count the same
operations ([docs/design/baselines.md](docs/design/baselines.md)). The mapping from each lemma
and algorithm of the paper to code and tests is [docs/PAPER_MAPPING.md](docs/PAPER_MAPPING.md).

## How correctness is checked

Correctness comes before speed. Every implementation is checked by tools that do not share its
code:

- **Independent oracle.** Bellman–Ford, itself compared with Floyd–Warshall on *every* small
  weighted simple graph from every source (exhaustive enumeration).
- **Shortest-path-tree verifier** (`sssp::verify_sssp`). It checks `dist[source] = 0`, that no
  edge can still shorten a distance, that every tree edge is tight, that parents lead back to
  the source, and that exactly the unreachable vertices are at infinity. This proves a result
  correct without trusting the algorithm.
- **Exact rational reference.** A Python oracle using exact fractions
  (`tests/python/reference/`) for floating-point inputs.
- **Exact versus tolerant comparison.** Integer weights must match exactly. `double` results are
  compared with a rounding-error bound derived in
  [docs/design/numerics.md](docs/design/numerics.md); tolerances are never loosened to make a
  test pass, and fast-math is disabled in every build.
- **Operation counting.** `Counted<W>` wraps a weight type and counts comparisons and additions,
  so that the comparison–addition model can be checked.
- **Regression tests.** Every counterexample found is minimised and kept in `tests/regression/`.
- **No silent fallbacks.** Whenever an implementation handles part of the input by a direct
  search instead of the paper's construction, it records that in `Diagnostics.fallbacks`.

Code is written test-first, and the failing and passing test runs are recorded per task in
`docs/tdd/tasks/`.

## Benchmarks

**No results yet. Every performance figure is "not measured" until a saved run backs it.**

The protocol is fixed in advance ([docs/design/benchmark_protocol.md](docs/design/benchmark_protocol.md)):
Release build with diagnostics compiled out; the same CSR graph, inputs, machine, warm-up and
repetitions for every algorithm; graph construction excluded from timing; every result verified
before it counts; environment metadata saved with the results; median, mean, p95, standard
deviation and sample count reported; speedups relative to binary-heap Dijkstra; and the base-case
share shown next to every timing of the paper's algorithm.

The suites will be run with `make benchmark-smoke`, `make benchmark` and `make benchmark-full`,
and summarised with `make report`. The C++ benchmark runner is still being written, so these
targets do not produce results yet.

## Repository layout

| Path | Contents |
|---|---|
| `include/sssp/` | Public headers: graph types (`CsrGraph`), results (`SsspResult`, `Diagnostics`), baselines, generators, IO, verification |
| `src/baseline/` | Bellman–Ford and the Dijkstra variants, with their heaps |
| `src/graph/`, `src/verify/` | Graph IO and generators; shortest-path verification |
| `src/paper/` | Mode A, paper-faithful (in progress) |
| `src/experimental/` | Mode B (`paper_experimental/`) and Mode C (`practical/`) |
| `tools/` | Command-line tools (`sssp_gen`) |
| `tests/` | C++ tests (Catch2: unit, integration, differential, property, regression) and Python tests (pytest) |
| `benchmarks/`, `scripts/` | Benchmark drivers and scripts (setup, paper download, benchmark runs, analysis) |
| `docs/` | Research notes, design documents, decisions, task tracking |

## Documentation

Start with [docs/README.md](docs/README.md), which explains every document under `docs/`, who it
is for, and the order to read it in. Some documents are currently available only in Japanese. Where
a document exists in both languages, the Japanese version is the `.ja.md` file next to the
English one.

## Related work

The previous best bounds for this problem setting (undirected, non-negative real weights,
comparison–addition model) include Duan, Mao, Shu and Yin (FOCS 2023,
`O(m √(log n · log log n))`, randomised), its derandomisation by Yan (SOSA 2026), and a 2026
preprint by Kadria and Roditty. Classic reference points are Dijkstra with Fibonacci heaps,
Thorup's linear-time algorithm for integer weights (1999), and Pettie and Ramachandran (2005).
At least one other public implementation of the same paper exists.
[docs/research/prior_work_and_oss.md](docs/research/prior_work_and_oss.md) lists these works and
existing open-source code, with licences and notes on how comparable they are.

## Citing

If you use this software, please cite the original paper and, where relevant, this software.

```bibtex
@misc{gholami2026linear,
  title         = {Linear Single-Source Shortest Path, Almost},
  author        = {Gholami, Iman and Hajiaghayi, MohammadTaghi},
  year          = {2026},
  eprint        = {2611.34911},
  archivePrefix = {arXiv},
  note          = {Version 1}
}
```

Citation metadata for the software itself is in [CITATION.cff](CITATION.cff).

## License

The code and documentation in this repository are released under the [MIT License](LICENSE).

The paper is **not** included. arXiv distributes it under a licence that does not allow
redistribution, so `make fetch-paper` downloads it into `.cache/paper/`, which is git-ignored
([ADR 0005](docs/adr/0005-paper-not-redistributed.md)). The documents here cite the paper by
section, lemma and algorithm number instead of quoting it.

## Disclaimer

This is an independent reproduction. It is not affiliated with or endorsed by the authors of the
paper. Any mistakes in the implementation or in our reading of the paper are ours. When a result
here disagrees with the paper, that is a question about our reproduction first; we record such
issues and deviations openly in [docs/PAPER_MAPPING.md](docs/PAPER_MAPPING.md) and
[docs/adr/](docs/adr/).
