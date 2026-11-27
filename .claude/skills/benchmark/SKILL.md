---
name: benchmark
description: >-
  Run fair, reproducible SSSP benchmarks (Release build, warm-up, repetitions, environment capture, median/mean/p95/stddev), save results and plots, and report them without mixing in unmeasured numbers. Use whenever timing or comparing algorithms, or writing any performance statement.
when_to_use: >-
  "benchmark", "how fast is", "compare with Dijkstra", "measure scaling", "performance report", "perf regression".
argument-hint: "[smoke|standard|full]"
---
# Benchmark procedure ($ARGUMENTS, default: smoke)

## 1. Build
```bash
cmake --preset release && cmake --build --preset release -j
```
- `release` preset = `-O3`, `SSSP_DIAGNOSTICS=OFF` (counters compiled out). Operation counts
  come from a separate `release-diag` run; never time a diagnostics build.
- Never add `-ffast-math`/`-Ofast`/`-march=native` to a compared build unless every algorithm
  gets the same flags and the report says so.
- Correctness first: `ctest --preset release` must pass on the same build.

## 2. Prepare the machine (macOS)
- AC power, Low Power Mode off (`pmset -g | grep -i lowpowermode` → 0), close heavy apps.
- Apple Silicon has P- and E-cores and macOS does not allow pinning. Run in the foreground at
  default QoS (never under `taskpolicy -b`), one benchmark process at a time, and note
  `hw.perflevel0/1` core counts.
- Thermal: check `pmset -g therm` before and after; if `CPU_Speed_Limit < 100` the run is
  invalid. Insert cool-down pauses between long configurations.
- Linux CI runners are noisy: CI runs only `smoke` for regressions, never for reported numbers.

## 3. Run
```bash
make benchmark-smoke      # quick sanity, small graphs
make benchmark            # standard suite
make benchmark-full       # large graphs, BENCH_TIMEOUT applies
```
Direct Google Benchmark invocation (same flags the scripts use):
```bash
build/release/benchmarks/<bench_exe> \
  --benchmark_min_warmup_time=1 --benchmark_repetitions=15 \
  --benchmark_enable_random_interleaving=true \
  --benchmark_report_aggregates_only=false \
  --benchmark_out=results/raw/<run-id>/<bench>.json --benchmark_out_format=json
```
- Same `CsrGraph` instance, same sources (fixed seed list), same repetitions for every algorithm.
- Graph construction and I/O are excluded from the timed region; say so in the report.
- Time-out or crash = recorded as such (`status: timeout`), never dropped silently.

## 4. Environment metadata (stored next to the results)
```bash
{ date -Iseconds; git rev-parse HEAD; git status --porcelain | wc -l;
  sysctl -n machdep.cpu.brand_string hw.ncpu hw.perflevel0.physicalcpu hw.perflevel1.physicalcpu hw.memsize;
  sw_vers; c++ --version | head -1; cmake --version | head -1; pmset -g therm; } > results/raw/<run-id>/env.txt
```
A dirty working tree (`porcelain` > 0) must be flagged in the report.

## 5. Analyse and plot
```bash
uv run python scripts/analyze_results.py results/raw/<run-id>   # or: make report
```
Per (algorithm, graph family, n, m): **n reps, median, mean, p95, stddev, min**; plus the Mode
A/B fallback share from a `release-diag` run. Plots: time vs m (log-log) per family, time/m,
and speed-up relative to the binary-heap Dijkstra. Commit summaries (CSV/Markdown/PNG) under
`results/` and `docs/reports/`; raw JSON stays in `results/raw/` (git-ignored).

## 6. Reporting rules
- Every number in docs/PRs comes from a saved run; cite the run id and commit.
- Not measured → write "not measured". Never extrapolate, never copy figures from the paper as
  ours, never mix machines in one table.
- Scaling claims (e.g. "near linear") need a fit over ≥ 4 sizes spanning ≥ 2 orders of magnitude.
- Mode A results show the fallback share; a mostly-fallback run is not evidence about the
  paper's algorithm.
- Update the "Latest measurements" section of `docs/PROGRESS.md`.
