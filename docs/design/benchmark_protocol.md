# Benchmark measurement protocol

**English** | [日本語](benchmark_protocol.ja.md) | [简体中文](benchmark_protocol.zh-CN.md)

This document defines the measurement procedure for comparing the implementation of arXiv 2611.34911v1 (undirected, non-negative-weight SSSP, expected O(m + n log log log n)) with
Dijkstra-family baselines. It contains **only the procedure and the settings**, and no measured
numbers at all. Only results generated from actual measurements are placed in `results/` and `docs/reports/`.

Related documents:

- Record schema: [benchmark_schema.md](benchmark_schema.md)
- Runner CLI: [benchmark_runner_cli.md](benchmark_runner_cli.md)
- Graph generator (names and definitions of graph types and weight distributions): [graph_library.md](graph_library.md)

Implementation: `scripts/run_benchmarks.sh` (execution), `scripts/bench/plan.py` (job list),
`scripts/bench/envinfo.py` (environment), `scripts/bench/schema.py` (validation), `scripts/bench/stats.py` (statistics),
`scripts/analyze_results.py` (aggregation and figures, `make report`).

---

## 1. Purpose and scope

### 1.1 Purpose

- Measure, in a reproducible way, on which sizes and graph types the paper's algorithm is faster or slower than the Dijkstra family,
  on the same graph, the same weights and the same memory layout (the shared CSR of `include/sssp/graph.hpp`).
- Record time, memory and operation counters separately, and fix the build and the metric that serve as evidence for each claim.

### 1.2 Compared algorithms (`algorithm`)

| `algorithm` | Meaning | `mode` | Supported weight types |
|---|---|---|---|
| `dijkstra-binary` | Binary-heap Dijkstra (the reference for Speedup) | `baseline` | `uint64`, `double` |
| `dijkstra-4ary` | 4-ary-heap Dijkstra | `baseline` | `uint64`, `double` |
| `dijkstra-fibonacci` | Fibonacci-heap Dijkstra | `baseline` | `uint64`, `double` |
| `dijkstra-radix` | Radix-heap Dijkstra | `baseline` | **`uint64` only** |
| `paper-faithful` | Entry point of Mode A (Alg 8). The outer threshold N > L_N^64 does not hold at practical sizes, so the base case (l.5) is solved by scan-based Dijkstra (Lemma 3.2). Recorded as a fallback | `paper-faithful` | `uint64`, `double` |
| `paper-faithful-d1` | Same as above. The base case is solved by binary-heap Dijkstra (deviation D-1, ADR 0008) | `paper-faithful` | `uint64`, `double` |
| `paper` | The outer construction of the paper's algorithm (Alg 8 l.8-15). Reserved name, not implemented | `paper-faithful` | `uint64`, `double` |
| `paper-experimental` | The paper's framework with components and parameters changed experimentally | `paper-experimental` | `uint64`, `double` |
| `paper-optimized` | Optimised version aiming at practical speed. Reserved name, not implemented | `practical-optimized` | `uint64`, `double` |
| `practical-local` | Mode C local SSSP (`practical::local_sssp`). A paper-inspired implementation. See below | `practical-optimized` | **`uint64` only** |

Meaning of `mode`:

- `baseline`: the Dijkstra family.
- `paper-faithful`: exactly as the paper describes (constants and parameters also follow the paper). For checking the asymptotic complexity claim.
- `paper-experimental`: the paper's framework is kept while parameters or components are changed experimentally.
- `practical-optimized`: optimised aiming at practical speed. Differences from the paper must always be recorded.

Jobs of `dijkstra-radix` × `double` are not created by the plan (`scripts/bench/plan.py`). They are out of scope and are not made into failure records either.
`practical-local` × `double` is treated the same way (the runner refuses with exit code 3; it does not substitute another algorithm).

About `paper-faithful` / `paper-faithful-d1`:

- In one call, `paper_sssp` cannot separate from outside the zero-edge contraction (l.2), the base case (l.5) and the tree reconstruction (l.6).
  The whole is measured as `sssp`, and `preprocess` and `tree` are null (§4).
- Every time, one `FallbackRecord` (the vertex count N and edge count after contraction) goes as it is into the record's `fallback`
  (even in a build with diagnostics OFF). When reporting times, always put the share of this fallback next to them.
  These times are the time of "the paper's preprocessing + Dijkstra", not the time of the paper's outer construction.
- Verification of `double` uses `Tolerance::for_path_hops(2(n − 1))` (numerics.md §3.2). Integer weights must match exactly.
- The scan-based base case of `paper-faithful` is Θ(N²), so the runner refuses n > 2^17 with exit code 3
  (it does not substitute another algorithm). The plan does not create such pairs, and prints the number excluded and the reason to standard error
  (`plan.VERTEX_LIMITS`). `paper-faithful-d1` has no limit.

About `practical-local` (task P5-PRACTICAL-LOCAL-BENCH):

- It is a paper-inspired implementation and not the paper's algorithm itself. It has no outer construction (Sec 6 to 8), and
  passes the whole contracted connected component to a single local hierarchy (a component of Alg 4) for evaluation. **It does not claim the complexity of Theorem 1.**
- Only `uint64` is registered. For general `double`, tree reconstruction (Lemma 8.5) may throw `invalid_argument`, so
  no `double` solver is registered (`--weight-type double` gives exit code 3). Tree reconstruction for general `double` is not supported.
- If the vertex count N after contraction is 16 or less (default `direct_search_max_vertices = 16`; this is not a value tuned by measurement),
  it solves directly with binary-heap Dijkstra and records one `FallbackRecord` (N = 1 gives one record for a different reason).
  For larger N it builds the hierarchy, and there are 0 fallbacks (the record's `fallback.count = 0`).
- As with `paper_sssp`, one call is measured as `sssp`, and `preprocess` and `tree` are null (§4).
- The fallback share is read per reason as `fallback.records[i].vertices / n` (N is after contraction, so ≤ n) and
  `fallback.records[i].edges / m`. Do not add across reasons (Q-69). This algorithm has at most one record, so
  `fallback.vertices` / `fallback.edges` equal that one record. When reporting times, always put this share next to them.
- Performance is **not measured** (it was not measured in this task). No number appears in this table.

### 1.3 Weight types (`weight_type`)

The names are the same as `sssp_gen --weight-type`.

- `uint64`: integer weights. Distances are computed exactly (overflow is detected and raises an exception).
- `double`: floating-point weights. Subject to IEEE-754 rounding. Comparison uses a tolerance (§5).

---

## 2. Build

Two presets of `CMakePresets.json` are used for different purposes.

| Preset | `SSSP_DIAGNOSTICS` | Use |
|---|---|---|
| `release` | `OFF` | For timing. Published times are taken only from this build |
| `release-diag` | `ON` | For recording operation counters (number of heap operations, relaxations, etc.) and fallbacks |

```sh
cmake --preset release      && cmake --build --preset release
cmake --preset release-diag && cmake --build --preset release-diag
```

Rules:

- **Claims about time are made from `release` only.** `release-diag` includes counter updates, so its times are not used for published values.
- **Claims about counters are made from `release-diag` only.** In `release` the counters are not recorded (null).
- Each record stores the build preset name, the compiler and its version, and `CMAKE_BUILD_TYPE`.

---

## 3. Trial design

| Item | Default | Notes |
|---|---|---|
| Warm-up | 1 for standard, 2 for full (0 for smoke) | Recorded as `warmup=true`, but excluded from statistics |
| Repetitions (development) | 5 (standard) | For initial checks and comparisons during development |
| Repetitions (publication) | 20 or more (full) | Published numbers must satisfy this condition |
| Repetitions of `release-diag` | 1 for smoke, 3 for standard, 5 for full | For counters. Not used for time |
| Execution order | The order of jobs (one process each) is shuffled with the order seed | `BENCH_ORDER_SEED` (default 12345). When one process runs several algorithms, the runner shuffles per trial with `--order-seed` |
| Sources | 1 for smoke, 2 for standard, 4 for full (`--sources`) | Chosen deterministically from `--seed` |

- The execution order is shuffled so that temporal biases in heat, cache and frequency are not
  concentrated on a particular algorithm. Because each run uses one process per algorithm for Peak RSS (§6), the unit of shuffling is
  the job (build × graph × algorithm).
- For the same (graph, seed), the same k sources are used for all algorithms. One trial = one run of one algorithm from one source.

---

## 4. Separation of times

Timestamps are taken with `std::chrono::steady_clock`, and all units are nanoseconds [ns].

| JSONL key (`times_ns.*`) | Name | Range |
|---|---|---|
| `graph_gen` | `graph_gen_ns` | Generation or loading of the graph |
| `preprocess` | `preprocess_ns` | Algorithm-specific preprocessing. null if not applicable |
| `sssp` | `sssp_ns` | The SSSP body (until the distances are final) |
| `tree` | `tree_ns` | Construction of the shortest-path tree (parent array). null for implementations that cannot separate it from the body, with a note to that effect |
| `end_to_end` | `end_to_end_ns` | Preprocessing + body + tree. **Does not include graph generation** |
| `verify` | `verify_ns` | Verification of correctness. Must always be done outside the measured interval and is not included in `end_to_end` |

Rules:

- An unmeasured value is **null**. Never put 0 (it could not be told apart from "measured, and it was 0").
- For an implementation whose `tree` is null, the tree construction time is included in `sssp`. Note this difference when comparing.
- Verification, log output and JSON writing are not put into the measured interval.

---

## 5. Correctness verification

Verification is done outside the measured interval (recorded separately as `verify_ns`).

1. **Agreement of distances**: compare with the distances of a reference implementation at all vertices.
   - The reference implementation is `dijkstra-binary`. For small graphs `bellman-ford` is also used, to check the reference implementation itself.
   - `uint64`: exact match.
   - `double`: comparison with a tolerance (relative and absolute error together). The tolerance follows `docs/design/numerics.md`.
   - An unreachable vertex must be infinity in both.
2. **Tree consistency**: for every reachable vertex v other than the source, the parent edge (parent[v], v) must exist,
   `dist[parent[v]] + w(parent[v], v) = dist[v]` must hold (with a tolerance for `double`), and following parent must return to the source.
3. A record that fails verification gets `status=error` and is excluded from time statistics.

---

## 6. Peak RSS (memory)

- Run only one (algorithm, graph) per process (equivalent to the runner's `--isolate`).
  `scripts/run_benchmarks.sh` starts a separate process per algorithm, so that the allocations of the previous algorithm are not mixed in.
- Just before the process exits, obtain `getrusage(RUSAGE_SELF).ru_maxrss`.
- The unit differs by OS. **macOS uses bytes, Linux uses KiB.** On Linux multiply by 1024 to normalise to bytes and record it in `peak_rss_bytes`.
- The value just after graph construction is also recorded as `rss_after_graph_bytes`, and `peak_rss_bytes - rss_after_graph_bytes` is
  taken as a rough measure of the algorithm's own increment.
  - Note: `ru_maxrss` is a monotonically increasing maximum, so this difference is only a rough measure, not an exact increment
    (for example the difference is too small when temporary space of graph generation created the peak).

---

## 7. Notes on Apple Silicon (macOS)

### 7.1 CPU pinning is not possible

- Apple Silicon has P cores (performance) and E cores (efficiency). The numbers can be checked with:
  ```sh
  sysctl hw.perflevel0.physicalcpu   # number of P cores
  sysctl hw.perflevel1.physicalcpu   # number of E cores
  ```
- The macOS scheduler may place the measured process on an E core.
- **On macOS there is no way to guarantee core pinning (CPU pinning).** `taskpolicy` and the like give priority or QoS instructions and are
  not pinning to a specific core. Measured values on macOS are therefore treated as including this variability.

### 7.2 Sources of variability

- Thermal throttling (frequency drops during long continuous runs).
- Power: AC or battery, and whether low power mode is on. Check with `pmset -g` and `pmset -g batt`.
- Background processing: Spotlight indexing, Time Machine, browsers, sync clients, etc.

### 7.3 Countermeasures

- Connect to AC power. Turn low power mode off.
- Quit other applications. If possible, run at a time when Spotlight indexing and Time Machine backups are not running.
- Put a cool-down (waiting time) between trials.
- The state that can be recorded is recorded by `scripts/bench/envinfo.py` for each run (in each record and in
  `results/raw/<run_id>.<preset>.env.json`; the state at the end is in `<run_id>.<preset>.env-end.json`).
  - Power state (AC/battery, low power mode), `CPU_Speed_Limit` from `pmset -g therm`, CPU name,
    OS and version, number of P/E cores, load average (`os.getloadavg()`, the same value as `uptime`).
- A run whose `cpu_speed_limit` at the start or at the end is below 100 is invalid, as having been limited by heat.
- A run for which the environment information at the start or at the end (`env.json` / `env-end.json`) could not be obtained is
  invalid, because the above judgement cannot be made. `scripts/run_benchmarks.sh` prints the failure to obtain it to standard error, records the number of failures at the end
  in the last `done:` line as `env_end_failed=<number>`, and finishes with exit code 1 (a failure at the start exits with code 1 on the spot).
  Records already measured are left in place, but the values of that run are not used for publication.
  - If it was obtained successfully but `cpu_speed_limit` is `null` (an environment with no means of obtaining it, such as Linux), that is not a failure.
    Write in the report that the run could not be judged for heat.
- For a run whose recorded load average is high, decide individually whether to use it for published values, and keep the reason if it is not used.

### 7.4 On Linux

- `taskset -c <core>` can pin the core.
- `cpupower frequency-set -g performance` can fix the frequency governor (if you have permission).
- `envinfo.py` does not obtain these settings automatically. If they were used, write them in the run's notes (the report).

---

## 8. Timeouts and memory limits

### 8.1 Timeouts

`scripts/run_benchmarks.sh` applies an external timeout to each invocation.

| Profile | Default timeout (per job) |
|---|---|
| `smoke` | 60 seconds |
| `standard` | 600 seconds |
| `full` | 3600 seconds (`make benchmark-full` passes the Makefile's `BENCH_TIMEOUT`, default 3600) |

If the environment variable `BENCH_TIMEOUT` is set, that value is used for any profile.

macOS does not ship GNU `timeout` by default, so it falls back in this order:

1. `timeout` (GNU coreutils)
2. `gtimeout` (Homebrew's coreutils)
3. `perl -e 'alarm ...; exec ...'` (the macOS default; the exit code on timeout is 142, and 124 for GNU `timeout`)

### 8.2 Memory limit

- The limit is applied with `ulimit -v` (environment variable `BENCH_MEM_LIMIT_MB`).
- On macOS `RLIMIT_AS` cannot be set (`ulimit -v` fails with `Invalid argument`). In that case a warning is printed and the run
  continues without a limit, and after the run records whose `peak_rss_bytes` exceeded the limit are shown as a warning (the status is not changed).
- If the runner returns exit code 5 because of `std::bad_alloc`, the status is `memory_limit`.

### 8.3 Handling of failures

- Timeouts, memory overruns and abnormal termination (including verification failure) are each kept as a record with `status=timeout` / `memory_limit` / `error`.
- If the runner could not write anything, `run_benchmarks.sh` appends one failure record (`python -m bench.schema failure`).
- These are excluded from time statistics and reported separately per job as the **failure rate** (`results/summary/status.csv`).
- Timeouts and memory overruns are expected, so the exit code of `run_benchmarks.sh` is 0. If there is even one abnormal termination (`error`), it is 1.
- Values of failed pairs are not estimated or interpolated.

---

## 9. Metrics and statistics

### 9.1 Speedup

```
Speedup = T(dijkstra-binary) / T(target)
```

- Within the same combination of (`graph_type`, `n`, `avg_degree`, `weight_dist`, `weight_type`), the **median** of each algorithm is used.
- Speedup > 1 means the target is faster than `dijkstra-binary`.
- The default T is `sssp_ns`. The `end_to_end_ns` version is also output.
- No Speedup is given for a combination in which either the baseline or the target has no successful measurement (`status=ok`, not warm-up, `release` build).
- Those with different `params` are separate series, and Speedup is given for each.

### 9.2 Statistics

| Statistic | Handling |
|---|---|
| median | Primary metric |
| mean | For reference |
| p95 | Flagged as reference-only when the sample count is < 20 |
| stddev | Sample standard deviation (ddof=1). null for a single sample |
| min / max | For reference |

- p95 is linear interpolation of order statistics (numpy's default).

- Warm-up (`warmup=true`) and failure records are not included in the statistics.

### 9.3 Confidence interval of Speedup

- Obtained by bootstrap. Defaults are 10000 resamples, random seed 20261127 (fixed), and a 95% percentile interval
  (`scripts/analyze_results.py --bootstrap-resamples / --bootstrap-seed`).
- The baseline side and the target side are resampled independently, and the ratio of medians is computed each time.
- The bootstrap seed is a setting of the aggregation side and does not enter the records. If the aggregation is done with changed default values, state it in the report.

---

## 10. Graph sizes and profiles

### 10.1 Size grid

- n ∈ {1e3, 1e4, 1e5, 1e6, 1e7}
- Average degree d ∈ {2, 4, 8, 16, 32} (d = 2m/n for an undirected graph)

### 10.2 Profiles

| Profile | n | d | Graph types | Weight distributions | trials | warmup | diag trials | Sources |
|---|---|---|---|---|---|---|---|---|
| `smoke` | {1e3, 1e4} | 4 | `er`, `grid` | `uniform_int` | 1 | 0 | 1 | 1 |
| `standard` | {1e3, 1e4, 1e5, 1e6} | {4, 16} | all types | `constant`, `uniform_int`, `wide_range` | 5 | 1 | 3 | 2 |
| `full` | {1e3, …, 1e7} | {2, 4, 8, 16, 32} | all types | all distributions | 20 | 2 | 5 | 4 |

- The definition is `PROFILES` in `scripts/bench/plan.py`. The job list can be checked with `python -m bench.plan <profile> --order-seed S`.
- Narrowing down: `BENCH_ALGORITHMS`, `BENCH_GRAPH_TYPES`, `BENCH_SIZES`, `BENCH_BUILDS` (comma-separated).

- `smoke` is for checking that things work, and its results are not used for performance claims.
- Depending on the size, the time or memory limit may be reached. Pairs that reach it are recorded as failure records and are not interpolated.
- Families whose degree is fixed by the family (`grid`, `road`, `tree`, `adversarial_long_path`) are not duplicated per d (`avg_degree` is null).
- Families whose weights are fixed by construction (`road`, `adversarial_long_path`) do not vary the weight distribution (`weight_dist=builtin`).

---

## 11. Graph types (`graph_type`)

Names and definitions follow `sssp_gen --type` (docs/design/graph_library.md §6).

| `graph_type` | Meaning | How d is passed | Notes |
|---|---|---|---|
| `er` | Uniform random G(n, m), m = round(n·d/2) | `--avg-degree d` | |
| `regular` | Random d-regular | `--degree d` | |
| `ba` | Barabási–Albert (preferential attachment) | `--attach d/2` | Average degree is about 2·attach |
| `chung_lu` | Chung–Lu (power-law expected degrees) | `--avg-degree d` | Large degree skew |
| `grid` | 2-D grid ⌊√n⌋ × ⌊√n⌋ | none (degree about 4) | n is rounded to a perfect square |
| `road` | Road-network-like (perturbed grid + short-range shortcuts) | none | Weights are Euclidean lengths (`builtin`) |
| `tree` | Uniform random tree | none | m = n − 1 |
| `adversarial_long_path` | Long path + detour edges | none | Deep shortest-path tree. Weights are fixed by construction (`builtin`) |

- `adversarial_decrease_key` has an edge count on the order of n², so it is not put into the size grid. If it is used, run it separately with a small n.
- Record the n and m that were actually generated (they can differ from the specified values because of rounding and removal of duplicate edges).
- Connectivity: vertices unreachable from the source are handled as they are, and the number of reachable vertices is recorded as `reachable`.

---

## 12. Weight distributions (`weight_dist`)

Names and definitions follow `sssp_gen --weights` (docs/design/graph_library.md §5). The parameters (`lo`, `hi`, `mean`,
`zero_probability`, `wide_log2`) use the generator's defaults in the profiles.

| `weight_dist` | Meaning | Weight type in profiles |
|---|---|---|
| `constant` | All edges have the same weight (default `lo`) | `uint64` |
| `uniform_int` | Uniform integer in [lo, hi] | `uint64` |
| `log_uniform` | Log-uniform | `uint64` |
| `wide_range` | Binary exponent chosen uniformly (large dynamic range) | `uint64` |
| `zero_mixed` | 0 with probability `zero_probability`, otherwise uniform | `uint64` |
| `uniform_real` | Uniform real in [lo, hi) | `double` |
| `exponential` | Exponential distribution | `double` |
| `builtin` | Weights determined by the generator's construction (`road`, `adversarial_long_path`) | `uint64` |

- For a combination of `uint64` where the sum of distances may overflow, `status=error` is used if it is detected.

---

## 13. Reproducibility

Every record stores the following.

- Seeds: the graph seed and the source-selection seed (`seed`), `order_seed`, and the algorithm's random seed (put in `params`).
  The bootstrap seed is a default of the aggregation side (§9.3)
- The git SHA and the dirty flag (whether the working tree has uncommitted changes)
- The paper version being reproduced: `2611.34911v1`
- Build preset, compiler and environment information (§7.3)

Rule:

- **Records with `dirty=true` are not used for published values.** They are used only for checks during development.

---

## 14. Output and aggregation

```
results/raw/<run_id>.jsonl          (already gitignored, one record per line)
        │
        └─ make report  (scripts/analyze_results.py)
              ├─ results/summary/*.csv
              └─ docs/reports/figures/*.png
```

- CSVs in `results/summary/`:
  - `timing_summary.csv`: samples / median / mean / p95 / stddev / min / max /
    p95_reference_only / dirty_records / run_ids for each configuration × algorithm × params × interval
  - `speedup.csv`: Speedup and 95% CI (`sssp` and `end_to_end`)
  - `memory.csv`: median / max of Peak RSS, and median of the RSS just after graph construction
  - `status.csv`: ok / timeout / memory_limit / error per job and the failure rate
  - `fallback.csv`: fallback occurrence rate in `release-diag`, counts, vertex and edge shares, and reasons
  - `ops.csv`: median of the operation counters and `counters`
- Figures in `docs/reports/figures/`: `runtime_vs_n__<type>`, `runtime_vs_m__<type>`, `speedup_vs_n__<type>`,
  `memory_vs_n__<type>`, `fallback_rate__<type>`, `breakdown__<type>` (breakdown into preprocessing/body/tree),
  `speedup_by_graph_type`, `speedup_by_weight_dist`, `param__<algorithm>__<param>` (threshold parameter and time).
- For the field definitions of records see [benchmark_schema.md](benchmark_schema.md), and for the runner's arguments see
  [benchmark_runner_cli.md](benchmark_runner_cli.md).
- In figures and tables, unmeasured pairs and failed pairs are not drawn. If a point is missing in the middle of a series, the line is broken there (no interpolation).
- If there are records that violate the schema, `make report` fails (`--skip-invalid` excludes them and continues).
- When there is no data to aggregate, it prints `データなし` (no data) and does not make up values.

---

## 15. Notes on interpreting results

- The difference in asymptotic complexity (O(m + n log n) versus O(m + n log log log n)) can be hidden by constant factors and cache effects at practical n.
  Not being faster is not necessarily a bug in the implementation nor an error in the paper.
- The op counters (`release-diag`) are a proxy for complexity in the comparison-addition model, and need not be proportional to wall-clock time
  (memory access patterns and branch prediction are not counted).
- The expected complexity of `paper` is an expectation over the randomisation. The variation between trials includes both the randomisation and the variability of the system.
- Advantages and disadvantages can change with graph type and weight distribution. Do not generalise the result of a particular combination.
- Always state in the report that results on macOS were obtained without CPU pinning.
