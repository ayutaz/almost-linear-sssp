# CLI specification of `sssp_bench_runner`

**English** | [日本語](benchmark_runner_cli.ja.md) | [简体中文](benchmark_runner_cli.zh-CN.md)

This document defines the arguments, output and exit codes of the C++ benchmark runner `build/<preset>/sssp_bench_runner`.
The implementation is `tools/sssp_bench_runner.cpp` and `src/bench_runner/`. `scripts/run_benchmarks.sh` calls it exactly as specified here.
The format of the output records is in [benchmark_schema.md](benchmark_schema.md), and the measurement procedure is in
[benchmark_protocol.md](benchmark_protocol.md).

Graph specification uses **the same names and the same meanings** as the graph generator `sssp_gen` (docs/design/graph_library.md).
The runner is expected to link the generator library (`sssp::gen`) and use the same option parsing.

## 1. Invocation example

`run_benchmarks.sh` starts it once per job (1 build x 1 graph x 1 algorithm):

```sh
build/release/sssp_bench_runner \
  --type er --n 100000 --avg-degree 4 --weights uniform_int --weight-type uint64 --seed 1 \
  --sources 2 --algorithms paper-faithful-d1 --trials 5 --warmup 1 --order-seed 12345 \
  --timeout-sec 600 --verify reference \
  --env-json results/raw/<run_id>.release.env.json --run-id <run_id> \
  --out results/raw/<run_id>.jsonl --isolate
```

## 2. Arguments

### 2.1 Graph (same as `sssp_gen`)

| Argument | Meaning |
|---|---|
| `--type T` | `er`, `regular`, `ba`, `chung_lu`, `grid`, `road`, `tree`, `adversarial_long_path` (`adversarial_decrease_key` may also be accepted, but its edge count is on the order of n², so profiles do not use it) |
| `--n N` | Number of vertices (`grid`/`road` use `rows = cols = ⌊√N⌋`) |
| `--avg-degree D` | Average degree for `er`, `chung_lu` |
| `--degree D` | Degree for `regular` |
| `--attach M` | Number of attachments per vertex for `ba` (average degree is about 2M) |
| `--weights DIST` | `constant`, `uniform_int`, `uniform_real`, `log_uniform`, `exponential`, `zero_mixed`, `wide_range`. For families whose weights are fixed by construction (`road`, `adversarial_*`), `builtin` is passed (treated as "no weight specification") |
| `--weight-type uint64\|double` | Weight type. This name is also written to the record's `weight_type` |
| `--lo`, `--hi`, `--mean`, `--zero-prob`, `--wide-log2`, `--exponent`, ... | Same as `sssp_gen` (generator defaults when omitted). When given, they are kept in a separate run note, not in `params` (the current schema has no field for graph parameters) |
| `--seed S` | Seed for graph generation and source selection |
| `--graph-file PATH` (optional) | Reads a file made by `sssp_gen`. In that case `graph_gen` is the loading time |

The record's `avg_degree` is the value of `--avg-degree` / `--degree`, and `2 × --attach` for `ba`.
For families whose degree is fixed by the family (`grid`, `road`, `tree`, `adversarial_long_path`) it is null.

### 2.2 Measurement

| Argument | Default | Meaning |
|---|---|---|
| `--sources K` | 1 | Number of sources. K are chosen deterministically from `--seed` (the same graph and seed give the same sources for every algorithm) |
| `--source-list a,b,...` | — | Explicit sources (mutually exclusive with `--sources`) |
| `--algorithms A[,B,...]` | required | Names are the schema's `ALGORITHM_MODES`. `mode` is determined by the name |
| `--param KEY=VALUE` | — | Repeatable. Parameters of non-baseline algorithms. Put into the record's `params` as they are (numbers as numbers) |
| `--trials T` | 5 | Number of measured trials (per source) |
| `--warmup W` | 1 | Number of warm-up runs. Records are written with `warmup=true` |
| `--order-seed S` | 0 | With several algorithms, the execution order of each trial is shuffled with this seed |
| `--timeout-sec X` | none | Internal guideline. When it has elapsed, finish without starting a new trial (forced termination is by the outer `timeout`) |
| `--verify none\|reference\|bellman-ford` | `reference` | Verify outside the measured interval. `reference` checks distance agreement with `dijkstra-binary` and tree consistency |
| `--isolate` | off | One process, one algorithm. Exit code 2 unless `--algorithms` has exactly one entry. Peak RSS is recorded only in this mode |

### 2.3 Output

| Argument | Meaning |
|---|---|
| `--env-json PATH` | Output of `scripts/bench/envinfo.py`. All its keys (`schema.ENV_FIELDS`) are copied as they are into every record |
| `--run-id ID` | The record's `run_id` |
| `--out PATH` | Writes JSONL in **append** mode. One record per line, flushed per line |
| `--list` | Prints the supported `--type`, `--weights` and algorithms (`name`, `mode`, `weight_types`, `takes_params`, `max_work`, `max_vertices`) to standard output as JSON and exits. `max_work` (upper bound on n × m) and `max_vertices` (upper bound on n) are null when there is no bound |
| `--version` | Prints the version and the paper version `2611.34911v1` and exits |

## 3. Behaviour

1. Generate (or load) the graph exactly once and measure `graph_gen_ns`. Immediately afterwards take `ru_maxrss` as
   `rss_after_graph_bytes`. The same `CsrGraph` is used for all algorithms and all trials.
2. For each source: W warm-up runs, then T measured runs. One record per run.
3. Time is `std::chrono::steady_clock`, in ns. `preprocess` / `sssp` / `tree` / `end_to_end` are measured separately.
   An interval that cannot be separated is null (not 0).
4. Verification, building the JSON and writing it are outside the measured interval. If verification fails: `correctness.ok=false`, `status=error`.
5. Counters are written to `ops` / `counters` / `fallback` only in a build with `SSSP_DIAGNOSTICS=1` (`release-diag`).
   In a build with `SSSP_DIAGNOSTICS=0`, `ops` and `counters` are always null. `fallback` is written if it can be recorded.
6. With `--isolate`, the records are written all at once just before process exit, and every record's `peak_rss_bytes` gets
   `getrusage(RUSAGE_SELF).ru_maxrss` taken just before exit. **macOS reports bytes and Linux reports KiB, so multiply by 1024** to
   normalise to bytes. Without `--isolate`, `peak_rss_bytes` is null and `isolated=false`.
7. If `dijkstra-radix` receives `--weight-type double`, write nothing and exit with code 3.
   Likewise, when a size limit is exceeded, write nothing and exit with code 3 (print the reason to standard error; never substitute another algorithm):
   `bellman-ford` for n × m > 2^30, `paper-faithful` for n > 2^17 (the scan-based Dijkstra base case is Θ(N²);
   `kScanBaseCaseMaxVertices`). `paper-faithful-d1` has no limit.
8. Written records must pass `python -m bench.schema validate` (checked by the CI smoke job).

## 4. Exit codes

| Code | Meaning | Handling in `run_benchmarks.sh` |
|---|---|---|
| 0 | Normal (also 0 when records with failed verification are included; failure is expressed by the record's `status=error`) | Success |
| 1 | Other abnormal termination | Append a failure record with `status=error` and finally exit with code 1 |
| 2 | Invalid arguments | Same as above |
| 3 | Unsupported combination (e.g. radix × double), reserved name that is not implemented (such as `paper`), size limit exceeded (§3, item 7) | Same as above (the plan never produces this combination) |
| 5 | Out of memory such as `std::bad_alloc` | Failure record with `status=memory_limit` |
| 124 / 142 | Forced termination by the outer `timeout` / perl `alarm` | Failure record with `status=timeout` |
| 137 | SIGKILL (e.g. the Linux OOM killer) | `status=error` (reason noted) |

## 5. Correspondence with `run_benchmarks.sh`

- Profile → job list is produced by `scripts/bench/plan.py` (`python -m bench.plan <profile> --order-seed S` prints a TSV).
  Columns: `build, graph_type, n, avg_degree, shape_args, weight_dist, weight_type, algorithm, trials, warmup, sources, seed`.
  `shape_args` is `--avg-degree 4` / `--degree 4` / `--attach 2` / `-` (none).
- The job order is already shuffled by the order seed. One job = one process (`--isolate`).
- When the runner is missing: `run_benchmarks.sh` prints "sssp_bench_runner が見つかりません" (sssp_bench_runner not found) and exits with code 2.
  Only with `smoke --allow-missing-runner` (for CI) does it print "未実装" (not implemented) and exit with code 0.
