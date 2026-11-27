# 基准测试结果的 schema(JSONL,schema_version "1")

[English](benchmark_schema.md) | [日本語](benchmark_schema.ja.md) | **简体中文**

基准测试结果保存在 `results/raw/<run_id>.jsonl` 中,每行一条 JSON 记录。
验证器 `scripts/bench/schema.py`(`validate_record`)实现了与本文档相同的规则。
测量步骤见 [benchmark_protocol.zh-CN.md](benchmark_protocol.zh-CN.md),运行器的参数见
[benchmark_runner_cli.zh-CN.md](benchmark_runner_cli.zh-CN.md)。

```sh
PYTHONPATH=scripts uv run python -m bench.schema validate results/raw/<run_id>.jsonl   # 0: 有效, 1: 存在违规
```

## 1. 总体规则

- **下表中的键全部必填。** 没有值时也不要省略键,而要写 `null`。表中没有的键视为违规(用于检测拼写错误)。
- **未测量的值为 `null`。不要填 `0`。** 时间 `times_ns.*` 为 1 以上的整数或 `null`。
  `peak_rss_bytes`、`rss_after_graph_bytes`、`ram_bytes` 为 1 以上的整数或 `null`。
- 整数字段中不要放 `true`/`false` 或 `1000.0`(即使是 JSON 数值,也要按整数书写)。
- 1 条记录 = 1 次 SSSP 运行(1 个算法 × 1 个图 × 1 个起点 × 1 次试验)。
  但失败记录(§4)表示运行器的 1 个进程(作业)没能输出任何内容。
- 不记录主机名、用户名和本地绝对路径(因为是公开仓库)。

## 2. 字段

### 2.1 运行与可复现性

| 键 | 类型 | null | 含义 |
|---|---|---|---|
| `schema_version` | string | 不可 | 始终为 `"1"` |
| `run_id` | string | 不可 | `scripts/run_benchmarks.sh` 一次执行的 ID(例如 `standard-20261127T013000Z-d5a55e7`) |
| `timestamp` | string | 不可 | ISO-8601,**必须带 UTC 偏移**(例如 `2026-11-27T10:30:00+09:00`) |
| `paper_version` | string | 不可 | 复现对象的论文版本。`2611.34911v1`(格式 `NNNN.NNNNNvN`) |
| `git_sha` | string | 不可 | 所测量源码的提交(十六进制 7 到 40 位) |
| `git_dirty` | bool | 不可 | 工作树中是否有未提交的更改。`true` 的结果不用于公开的数值 |

### 2.2 算法

| 键 | 类型 | null | 含义 |
|---|---|---|---|
| `algorithm` | string | 不可 | `dijkstra-binary` / `dijkstra-4ary` / `dijkstra-fibonacci` / `dijkstra-radix` / `bellman-ford` / `paper-faithful` / `paper-faithful-d1` / `paper` / `paper-experimental` / `paper-optimized` / `practical-local`(也允许未知名称) |
| `mode` | enum | 不可 | `baseline` / `paper-faithful` / `paper-experimental` / `practical-optimized` |
| `params` | object | 不可 | 算法的参数(阈值、随机数 seed 等)。没有则为 `{}` |

已知名称与 `mode` 的对应是固定的(`schema.ALGORITHM_MODES`):Dijkstra 系和 `bellman-ford` 为 `baseline`;
`paper-faithful` / `paper-faithful-d1`(Mode A 的入口,`include/sssp/paper.hpp` 的 `algorithm_name`)和
`paper` 为 `paper-faithful`;`paper-experimental` 为 `paper-experimental`;`paper-optimized` 为
`practical-optimized`;`practical-local`(Mode C,参考论文的实现,不主张 Theorem 1 的复杂度;仅限 `uint64`,
`weight_type = double` 的记录会被验证器拒绝)为 `practical-optimized`。`baseline` 专用于这 5 个基准算法。`paper` 是外层构造
(Alg 8 l.8-15)的保留名称,尚未实现(运行器以退出码 3 结束)。

### 2.3 环境(由 `scripts/bench/envinfo.py` 生成,运行器从 `--env-json` 原样复制)

| 键 | 类型 | null | 含义 |
|---|---|---|---|
| `os` | string | 不可 | 例如 `macOS 26.2 (25C56); Darwin 25.2.0`、`Ubuntu 24.04 LTS; Linux 6.8.0` |
| `cpu` | object | 不可 | `{brand, logical_cores, physical_cores, perf_cores, efficiency_cores}`。各值可为 null。P/E 核心数在 macOS 上取自 `hw.perflevel0/1.physicalcpu`,在 Linux 上为 null |
| `ram_bytes` | int | 可 | 物理内存量 |
| `arch` | string | 不可 | `platform.machine()`(`arm64`、`x86_64` 等) |
| `compiler` | string | 可 | `<CMAKE_CXX_COMPILER> --version` 的第 1 行 |
| `build_type` | string | 不可 | `CMAKE_BUILD_TYPE`。无法读取 CMakeCache.txt 的环境中的记录无效 |
| `build_flags` | string | 可 | `compile_commands.json` 中的实际有效标志(不含路径和输入输出)。如果没有,则取 CMakeCache 的 `CMAKE_CXX_FLAGS[_<TYPE>]` |
| `diagnostics_enabled` | bool | 不可 | `SSSP_DIAGNOSTICS`(`release-diag` 时为 true) |
| `power` | object | 可 | `{source: "AC"/"Battery"/"UPS"/null, battery_percent, low_power_mode, cpu_speed_limit}`。`cpu_speed_limit` 为 `pmset -g therm` 的 `CPU_Speed_Limit`(没有记录则为 null) |
| `load_average` | [number×3] | 可 | 测量开始时的 1/5/15 分钟 load average |

### 2.4 图

| 键 | 类型 | null | 含义 |
|---|---|---|---|
| `graph_type` | string | 不可 | `sssp_gen --type` 的名称(`er`、`grid`、`ba`、…) |
| `n` | int ≥ 1 | 不可 | 实际的顶点数 |
| `m` | int ≥ 0 | 仅失败时可 | 实际的无向边数(生成并去重之后) |
| `avg_degree` | number | 可 | 所指定的名义平均度数 d(2m/n 的目标)。度数由图族决定的(`grid`、`road`、`tree`、`adversarial_long_path`)为 null |
| `weight_dist` | string | 不可 | `sssp_gen --weights` 的名称。权重由构造决定的图族为 `builtin` |
| `weight_type` | enum | 不可 | `uint64` / `double`(与 `sssp_gen --weight-type` 相同) |
| `seed` | int ≥ 0 | 不可 | 图生成和起点选择的 seed |
| `source` | int | 仅失败时可 | 起点。`0 ≤ source < n` |
| `reachable` | int | 可 | 从起点可达的顶点数(`≤ n`) |

### 2.5 试验

| 键 | 类型 | null | 含义 |
|---|---|---|---|
| `trial` | int ≥ 0 | 仅失败时可 | 试验编号(预热也从 0 开始计数) |
| `warmup` | bool | 不可 | 是否为预热运行。true 的记录不计入统计 |
| `order_seed` | int | 可 | 执行顺序打乱的 seed |
| `isolated` | bool | 不可 | 该进程是否只运行了一个算法(`--isolate`) |
| `status` | enum | 不可 | `ok` / `timeout` / `memory_limit` / `error` |
| `error` | string | 可 | 失败的说明 |

### 2.6 测量值

| 键 | 类型 | null | 含义 |
|---|---|---|---|
| `times_ns` | object | 不可 | 恰好包含下面的 6 个键。各值为 1 以上的整数 [ns] 或 null |
| `peak_rss_bytes` | int | 可 | 退出前的 `ru_maxrss`(已换算为 bytes)。`isolated=false` 时必为 null |
| `rss_after_graph_bytes` | int | 可 | 图构建完成后立即取得的 `ru_maxrss`(bytes) |
| `ops` | object | 可 | `include/sssp/result.hpp` 中 `OpCounters` 的 8 个键(`vertex_scans`、`edge_scans`、`relaxations`、`comparisons`、`additions`、`queue_pushes`、`queue_pops`、`decrease_keys`)的子集。值为 0 以上的整数 |
| `counters` | object | 可 | `Diagnostics::counters`(名称 → 0 以上的整数) |
| `fallback` | object | 可 | `{count, records: [{reason, vertices, edges}], vertices, edges}`。`count = len(records)`,`vertices`/`edges` 为 records 的总和。未记录时为 null |
| `correctness` | object | 不可 | `{verified, ok, method}` |

`times_ns` 的键(各区间的定义见 benchmark_protocol.md §4):

| 键 | 含义 |
|---|---|
| `graph_gen` | 图的生成或读取 |
| `preprocess` | 算法特有的预处理(没有则为 null) |
| `sssp` | SSSP 主体 |
| `tree` | 最短路径树的构建(无法与主体分离时为 null) |
| `end_to_end` | 预处理 + 主体 + 树(不含图生成和验证) |
| `verify` | 正确性验证(在测量区间之外) |

`correctness`:

- `verified=false`(未验证):`ok` 和 `method` 为 null。
- `verified=true`:`ok` 为 bool,`method` 为字符串(例如 `reference:dijkstra-binary`、`bellman-ford`、`tree-check`)。

## 3. 一致性规则

| 规则 | 理由 |
|---|---|
| `status=ok` 时,`m`、`source`、`trial`、`times_ns.sssp`、`times_ns.end_to_end` 不可为 null | 成功的测量必定包含这些值 |
| `times_ns.end_to_end ≥ preprocess, sssp, tree`(两者都有值时) | 各部分之和不超过总量 |
| `correctness.ok=false` 时,`status=error` | 错误结果的时间不计入统计 |
| `diagnostics_enabled=false` 时,`ops` 和 `counters` 为 null | 在 release 中计数器会在编译时被去除。写 0 的话,无法与"统计过,结果为 0"区分 |
| `isolated=false` 时,`peak_rss_bytes` 为 null | 多个算法在同一进程中运行时,无法判断 Peak RSS 属于哪一个 |
| 已知算法的 `mode` 固定。`baseline` 专用于基准算法 | 避免对比的分类出错 |
| `dijkstra-radix` 只能是 `weight_type=uint64` | radix heap 以单调的整数键为前提 |
| `0 ≤ source < n`,`reachable ≤ n` | 范围 |

## 4. 失败记录

当运行器因超时、内存不足或异常终止而没能写出任何内容时,`scripts/run_benchmarks.sh` 会用
`python -m bench.schema failure ...` 追加 1 条记录(`schema.make_failure_record`)。

- `status` 为 `timeout` / `memory_limit` / `error`,`error` 中写原因(例如 `timeout after 600 s (exit 142)`)。
- 未知的值为 null:`m`、`source`、`reachable`、`trial`、`times_ns.*`、`peak_rss_bytes`、`rss_after_graph_bytes`、
  `ops`、`counters`、`fallback`。`correctness` 为 `{verified: false, ok: null, method: null}`。
- 环境字段从该作业的 `--env-json` 复制。

## 5. 示例(用于说明格式;数值不是测量值)

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

(在实际文件中,1 条记录写在 1 行里。)

## 6. 变更步骤

修改 schema 时,要提高 `schema_version`,并先(从测试开始)修改 `scripts/bench/schema.py` 和
`tests/python/test_bench_schema.py`。旧版本的记录会被验证拒绝。
