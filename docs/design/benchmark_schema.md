# Schema of benchmark results (JSONL, schema_version "1")

**English** | [日本語](benchmark_schema.ja.md) | [简体中文](benchmark_schema.zh-CN.md)

Benchmark results are stored in `results/raw/<run_id>.jsonl`, one JSON record per line.
The validator `scripts/bench/schema.py` (`validate_record`) implements the same rules as this document.
For the measurement procedure see [benchmark_protocol.md](benchmark_protocol.md), and for the runner's arguments see
[benchmark_runner_cli.md](benchmark_runner_cli.md).

```sh
PYTHONPATH=scripts uv run python -m bench.schema validate results/raw/<run_id>.jsonl   # 0: valid, 1: violations found
```

## 1. General rules

- **Every key in the tables below is required.** When there is no value, do not omit the key; write `null`. A key not in the tables is a violation (to catch typos).
- **An unmeasured value is `null`. Never put `0`.** The times `times_ns.*` are an integer >= 1 or `null`.
  `peak_rss_bytes`, `rss_after_graph_bytes` and `ram_bytes` are an integer >= 1 or `null`.
- Do not put `true`/`false` or `1000.0` in an integer field (write it as an integer even in JSON numbers).
- One record = one SSSP run (1 algorithm x 1 graph x 1 source x 1 trial).
  A failure record (§4), however, means that one runner process (job) could not output anything.
- Do not record host names, user names or local absolute paths (this is a public repository).

## 2. Fields

### 2.1 Run and reproducibility

| Key | Type | null | Meaning |
|---|---|---|---|
| `schema_version` | string | not allowed | Always `"1"` |
| `run_id` | string | not allowed | ID of one execution of `scripts/run_benchmarks.sh` (e.g. `standard-20261127T013000Z-d5a55e7`) |
| `timestamp` | string | not allowed | ISO-8601, **UTC offset required** (e.g. `2026-11-27T10:30:00+09:00`) |
| `paper_version` | string | not allowed | Paper version being reproduced. `2611.34911v1` (format `NNNN.NNNNNvN`) |
| `git_sha` | string | not allowed | Commit of the measured source (7 to 40 hexadecimal digits) |
| `git_dirty` | bool | not allowed | Whether the working tree had uncommitted changes. Results with `true` are not used for published values |

### 2.2 Algorithm

| Key | Type | null | Meaning |
|---|---|---|---|
| `algorithm` | string | not allowed | `dijkstra-binary` / `dijkstra-4ary` / `dijkstra-fibonacci` / `dijkstra-radix` / `bellman-ford` / `paper-faithful` / `paper-faithful-d1` / `paper` / `paper-experimental` / `paper-optimized` / `practical-local` (unknown names are also accepted) |
| `mode` | enum | not allowed | `baseline` / `paper-faithful` / `paper-experimental` / `practical-optimized` |
| `params` | object | not allowed | Algorithm parameters (thresholds, random seeds, etc.). `{}` if none |

The correspondence between known names and `mode` is fixed (`schema.ALGORITHM_MODES`): the Dijkstra family and `bellman-ford` are `baseline`;
`paper-faithful` / `paper-faithful-d1` (the Mode A entry point, `algorithm_name` in `include/sssp/paper.hpp`) and
`paper` are `paper-faithful`; `paper-experimental` is `paper-experimental`; `paper-optimized` is
`practical-optimized`; `practical-local` (Mode C, a paper-inspired implementation that does not claim the complexity of Theorem 1; `uint64` only,
and a record with `weight_type = double` is rejected by the validator) is `practical-optimized`. `baseline` is reserved for these 5 baseline algorithms. `paper` is a reserved name for the outer construction
(Alg 8 l.8-15) and is not implemented yet (the runner exits with code 3).

### 2.3 Environment (made by `scripts/bench/envinfo.py`; the runner copies it as is from `--env-json`)

| Key | Type | null | Meaning |
|---|---|---|---|
| `os` | string | not allowed | e.g. `macOS 26.2 (25C56); Darwin 25.2.0`, `Ubuntu 24.04 LTS; Linux 6.8.0` |
| `cpu` | object | not allowed | `{brand, logical_cores, physical_cores, perf_cores, efficiency_cores}`. Each value may be null. P/E core counts come from `hw.perflevel0/1.physicalcpu` on macOS and are null on Linux |
| `ram_bytes` | int | allowed | Amount of physical memory |
| `arch` | string | not allowed | `platform.machine()` (`arm64`, `x86_64`, etc.) |
| `compiler` | string | allowed | First line of `<CMAKE_CXX_COMPILER> --version` |
| `build_type` | string | not allowed | `CMAKE_BUILD_TYPE`. A record from an environment where CMakeCache.txt cannot be read is invalid |
| `build_flags` | string | allowed | Effective flags from `compile_commands.json` (excluding paths and inputs/outputs). If absent, `CMAKE_CXX_FLAGS[_<TYPE>]` from CMakeCache |
| `diagnostics_enabled` | bool | not allowed | `SSSP_DIAGNOSTICS` (true for `release-diag`) |
| `power` | object | allowed | `{source: "AC"/"Battery"/"UPS"/null, battery_percent, low_power_mode, cpu_speed_limit}`. `cpu_speed_limit` is `CPU_Speed_Limit` from `pmset -g therm` (null if there is no record) |
| `load_average` | [number×3] | allowed | 1/5/15-minute load average at the start of measurement |

### 2.4 Graph

| Key | Type | null | Meaning |
|---|---|---|---|
| `graph_type` | string | not allowed | Name from `sssp_gen --type` (`er`, `grid`, `ba`, ...) |
| `n` | int ≥ 1 | not allowed | Actual number of vertices |
| `m` | int ≥ 0 | only on failure | Actual number of undirected edges (after generation and removal of duplicates) |
| `avg_degree` | number | allowed | Nominal average degree d that was specified (the target of 2m/n). null for families whose degree is fixed by the family (`grid`, `road`, `tree`, `adversarial_long_path`) |
| `weight_dist` | string | not allowed | Name from `sssp_gen --weights`. `builtin` for families whose weights are fixed by construction |
| `weight_type` | enum | not allowed | `uint64` / `double` (same as `sssp_gen --weight-type`) |
| `seed` | int ≥ 0 | not allowed | Seed for graph generation and source selection |
| `source` | int | only on failure | Source. `0 ≤ source < n` |
| `reachable` | int | allowed | Number of vertices reachable from the source (`≤ n`) |

### 2.5 Trial

| Key | Type | null | Meaning |
|---|---|---|---|
| `trial` | int ≥ 0 | only on failure | Trial number (warm-up runs are also counted from 0) |
| `warmup` | bool | not allowed | Whether this is a warm-up run. true is excluded from statistics |
| `order_seed` | int | allowed | Seed of the execution-order shuffle |
| `isolated` | bool | not allowed | Whether this process ran only one algorithm (`--isolate`) |
| `status` | enum | not allowed | `ok` / `timeout` / `memory_limit` / `error` |
| `error` | string | allowed | Description of the failure |

### 2.6 Measured values

| Key | Type | null | Meaning |
|---|---|---|---|
| `times_ns` | object | not allowed | Has exactly the 6 keys below. Each value is an integer >= 1 [ns] or null |
| `peak_rss_bytes` | int | allowed | `ru_maxrss` just before exit (normalised to bytes). Always null when `isolated=false` |
| `rss_after_graph_bytes` | int | allowed | `ru_maxrss` just after graph construction (bytes) |
| `ops` | object | allowed | A subset of the 8 keys of `OpCounters` in `include/sssp/result.hpp` (`vertex_scans`, `edge_scans`, `relaxations`, `comparisons`, `additions`, `queue_pushes`, `queue_pops`, `decrease_keys`). Values are integers ≥ 0 |
| `counters` | object | allowed | `Diagnostics::counters` (name → integer ≥ 0) |
| `fallback` | object | allowed | `{count, records: [{reason, vertices, edges}], vertices, edges}`. `count = len(records)`, and `vertices`/`edges` are the sums over records. null when not recorded |
| `correctness` | object | not allowed | `{verified, ok, method}` |

Keys of `times_ns` (the definition of each interval is in benchmark_protocol.md §4):

| Key | Meaning |
|---|---|
| `graph_gen` | Generation or loading of the graph |
| `preprocess` | Algorithm-specific preprocessing (null if none) |
| `sssp` | The SSSP body |
| `tree` | Construction of the shortest-path tree (null if it cannot be separated from the body) |
| `end_to_end` | Preprocessing + body + tree (excludes graph generation and verification) |
| `verify` | Verification of correctness (outside the measured interval) |

`correctness`:

- `verified=false` (not verified): `ok` and `method` are null.
- `verified=true`: `ok` is a bool and `method` is a string (e.g. `reference:dijkstra-binary`, `bellman-ford`, `tree-check`).

## 3. Consistency rules

| Rule | Reason |
|---|---|
| If `status=ok`, `m`, `source`, `trial`, `times_ns.sssp` and `times_ns.end_to_end` must not be null | A successful measurement always has them |
| `times_ns.end_to_end ≥ preprocess, sssp, tree` (when both have values) | The parts do not exceed the total |
| If `correctness.ok=false`, then `status=error` | The time of an incorrect result is not put into statistics |
| If `diagnostics_enabled=false`, `ops` and `counters` are null | In release the counters are compiled out. Writing 0 could not be told apart from "counted, and it was 0" |
| If `isolated=false`, `peak_rss_bytes` is null | When several algorithms run in the same process, it is unclear which one a Peak RSS belongs to |
| `mode` of a known algorithm is fixed. `baseline` is for the baseline algorithms only | Do not misclassify the comparison |
| `dijkstra-radix` is `weight_type=uint64` only | The radix heap requires monotone integer keys |
| `0 ≤ source < n`, `reachable ≤ n` | Range |

## 4. Failure records

When the runner could not write anything because of a timeout, out of memory or abnormal termination, `scripts/run_benchmarks.sh` appends one record with
`python -m bench.schema failure ...` (`schema.make_failure_record`).

- `status` is `timeout` / `memory_limit` / `error`, and `error` holds the reason (e.g. `timeout after 600 s (exit 142)`).
- Unknown values are null: `m`, `source`, `reachable`, `trial`, `times_ns.*`, `peak_rss_bytes`, `rss_after_graph_bytes`,
  `ops`, `counters`, `fallback`. `correctness` is `{verified: false, ok: null, method: null}`.
- The environment fields are copied from that job's `--env-json`.

## 5. Example (to illustrate the format; the numbers are not measured values)

```json
{"schema_version": "1", "run_id": "EXAMPLE", "timestamp": "2026-11-27T10:00:00+09:00",
 "paper_version": "2611.34911v1", "git_sha": "0000000", "git_dirty": false,
 "algorithm": "paper", "mode": "paper-faithful", "params": {"algo_seed": 1},
 "os": "EXAMPLE", "cpu": {"brand": "EXAMPLE", "logical_cores": null, "physical_cores": null,
 "perf_cores": null, "efficiency_cores": null}, "ram_bytes": null, "arch": "arm64",
 "compiler": null, "build_type": "Release", "build_flags": null, "diagnostics_enabled": false,
 "power": null, "load_average": null,
 "graph_type": "er", "n": 1000, "m": 2000, "avg_degree": 4, "weight_dist": "uniform_int",
 "weight_type": "uint64", "seed": 1, "source": 17, "reachable": 1000,
 "trial": 0, "warmup": false, "order_seed": 12345, "isolated": true, "status": "ok", "error": null,
 "times_ns": {"graph_gen": 1, "preprocess": null, "sssp": 1, "tree": null, "end_to_end": 1, "verify": 1},
 "peak_rss_bytes": null, "rss_after_graph_bytes": null,
 "ops": null, "counters": null, "fallback": null,
 "correctness": {"verified": true, "ok": true, "method": "reference:dijkstra-binary"}}
```

(In an actual file, one record is written on one line.)

## 6. Procedure for changes

When changing the schema, raise `schema_version`, and fix `scripts/bench/schema.py` and
`tests/python/test_bench_schema.py` first (starting from the tests). Records of an older version are rejected by validation.
