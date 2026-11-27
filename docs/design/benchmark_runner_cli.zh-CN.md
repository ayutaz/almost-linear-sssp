# `sssp_bench_runner` 的 CLI 规范

[English](benchmark_runner_cli.md) | [日本語](benchmark_runner_cli.ja.md) | **简体中文**

本文档规定 C++ 基准运行器 `build/<preset>/sssp_bench_runner` 的参数、输出和退出码。
实现位于 `tools/sssp_bench_runner.cpp` 和 `src/bench_runner/`。`scripts/run_benchmarks.sh` 严格按本规范调用它。
输出记录的格式见 [benchmark_schema.zh-CN.md](benchmark_schema.zh-CN.md),测量步骤见
[benchmark_protocol.zh-CN.md](benchmark_protocol.zh-CN.md)。

图的指定方式与图生成器 `sssp_gen`(docs/design/graph_library.md)**使用相同的名称、相同的含义**。
运行器应链接生成器库(`sssp::gen`),并使用相同的选项解析。

## 1. 调用示例

`run_benchmarks.sh` 对每个作业(1 个构建 × 1 个图 × 1 个算法)启动一次:

```sh
build/release/sssp_bench_runner \
  --type er --n 100000 --avg-degree 4 --weights uniform_int --weight-type uint64 --seed 1 \
  --sources 2 --algorithms paper-faithful-d1 --trials 5 --warmup 1 --order-seed 12345 \
  --timeout-sec 600 --verify reference \
  --env-json results/raw/<run_id>.release.env.json --run-id <run_id> \
  --out results/raw/<run_id>.jsonl --isolate
```

## 2. 参数

### 2.1 图(与 `sssp_gen` 相同)

| 参数 | 含义 |
|---|---|
| `--type T` | `er`、`regular`、`ba`、`chung_lu`、`grid`、`road`、`tree`、`adversarial_long_path`(也可以接受 `adversarial_decrease_key`,但其边数为 n² 量级,因此配置档不使用) |
| `--n N` | 顶点数(`grid`/`road` 取 `rows = cols = ⌊√N⌋`) |
| `--avg-degree D` | `er`、`chung_lu` 的平均度数 |
| `--degree D` | `regular` 的度数 |
| `--attach M` | `ba` 中每个顶点的连接数(平均度数约为 2M) |
| `--weights DIST` | `constant`、`uniform_int`、`uniform_real`、`log_uniform`、`exponential`、`zero_mixed`、`wide_range`。对于权重由构造决定的图族(`road`、`adversarial_*`),传入 `builtin`(视为"未指定权重") |
| `--weight-type uint64\|double` | 权重类型。记录的 `weight_type` 中也写入这个名称 |
| `--lo`、`--hi`、`--mean`、`--zero-prob`、`--wide-log2`、`--exponent`、... | 与 `sssp_gen` 相同(省略时取生成器的默认值)。指定后不放入 `params`,而是另行记在 run 的备注中(当前 schema 没有图参数字段) |
| `--seed S` | 图生成和起点选择的 seed |
| `--graph-file PATH`(可选) | 读取由 `sssp_gen` 生成的文件。此时 `graph_gen` 为读取时间 |

记录的 `avg_degree` 写入 `--avg-degree` / `--degree` 的值,`ba` 则写 `2 × --attach`。
度数由图族决定的(`grid`、`road`、`tree`、`adversarial_long_path`)为 null。

### 2.2 测量

| 参数 | 默认值 | 含义 |
|---|---|---|
| `--sources K` | 1 | 起点数量。由 `--seed` 确定性地选出 K 个(相同的图和 seed 下,所有算法使用相同的起点) |
| `--source-list a,b,...` | — | 显式指定起点(与 `--sources` 互斥) |
| `--algorithms A[,B,...]` | 必填 | 名称取自 schema 的 `ALGORITHM_MODES`。`mode` 由名称决定 |
| `--param KEY=VALUE` | — | 可重复。非基准算法的参数。原样放入记录的 `params`(数值按数值处理) |
| `--trials T` | 5 | 正式测量的试验次数(每个起点) |
| `--warmup W` | 1 | 预热次数。记录以 `warmup=true` 输出 |
| `--order-seed S` | 0 | 有多个算法时,用该 seed 打乱每次试验的执行顺序 |
| `--timeout-sec X` | 无 | 内部参考值。超时后不再开始新的试验并结束(强制终止由外层 `timeout` 负责) |
| `--verify none\|reference\|bellman-ford` | `reference` | 在测量区间之外进行验证。`reference` 检查与 `dijkstra-binary` 的距离一致性以及树的一致性 |
| `--isolate` | off | 一个进程只运行一个算法。若 `--algorithms` 不是恰好一个,则退出码为 2。仅在此时记录 Peak RSS |

### 2.3 输出

| 参数 | 含义 |
|---|---|
| `--env-json PATH` | `scripts/bench/envinfo.py` 的输出。其全部键(`schema.ENV_FIELDS`)原样复制到每条记录中 |
| `--run-id ID` | 记录的 `run_id` |
| `--out PATH` | 以**追加**方式写入 JSONL。每行一条记录,每行 flush |
| `--list` | 以 JSON 形式向标准输出打印所支持的 `--type`、`--weights` 和算法(`name`、`mode`、`weight_types`、`takes_params`、`max_work`、`max_vertices`)后退出。`max_work`(n × m 的上限)和 `max_vertices`(n 的上限)在没有上限时为 null |
| `--version` | 打印版本和论文版本 `2611.34911v1` 后退出 |

## 3. 行为

1. 图只生成(或读取)一次,并测量 `graph_gen_ns`。紧接着取 `ru_maxrss` 作为
   `rss_after_graph_bytes`。所有算法、所有试验使用同一个 `CsrGraph`。
2. 对每个起点:预热 W 次,然后正式测量 T 次。每次一条记录。
3. 时间使用 `std::chrono::steady_clock`,单位 ns。`preprocess` / `sssp` / `tree` / `end_to_end` 分别测量。
   无法分开的区间为 null(不写 0)。
4. 验证、JSON 的组装和写出都在测量区间之外。验证失败时:`correctness.ok=false`、`status=error`。
5. 仅在 `SSSP_DIAGNOSTICS=1` 的构建(`release-diag`)中才向 `ops` / `counters` / `fallback` 写入计数器。
   在 `SSSP_DIAGNOSTICS=0` 的构建中,`ops` 和 `counters` 一律为 null。`fallback` 在能够记录时写入。
6. 使用 `--isolate` 时,记录在进程退出前一次性写出,所有记录的 `peak_rss_bytes` 都填入
   退出前的 `getrusage(RUSAGE_SELF).ru_maxrss`。**macOS 为 bytes,Linux 为 KiB,因此要 ×1024** 统一换算为
   bytes。不使用 `--isolate` 时,`peak_rss_bytes` 为 null,`isolated=false`。
7. 若 `dijkstra-radix` 收到 `--weight-type double`,不写任何内容并以退出码 3 结束。
   超过规模上限时同样不写任何内容并以退出码 3 结束(将原因输出到标准错误;不替换为其他算法):
   `bellman-ford` 在 n × m > 2^30 时,`paper-faithful` 在 n > 2^17 时(扫描型 Dijkstra 的基底情形为 Θ(N²);
   `kScanBaseCaseMaxVertices`)。`paper-faithful-d1` 没有上限。
8. 写出的记录必须通过 `python -m bench.schema validate`(由 CI 的 smoke 确认)。

## 4. 退出码

| 退出码 | 含义 | `run_benchmarks.sh` 的处理 |
|---|---|---|
| 0 | 正常(包含验证失败的记录时也为 0;失败由记录的 `status=error` 表示) | 成功 |
| 1 | 其他异常 | 追加 `status=error` 的失败记录,最后以退出码 1 结束 |
| 2 | 参数错误 | 同上 |
| 3 | 不支持的组合(如 radix × double)、尚未实现的保留名称(如 `paper`)、超过规模上限(§3 的第 7 条) | 同上(计划不会生成这种组合) |
| 5 | `std::bad_alloc` 等内存不足 | `status=memory_limit` 的失败记录 |
| 124 / 142 | 被外层 `timeout` / perl `alarm` 强制终止 | `status=timeout` 的失败记录 |
| 137 | SIGKILL(Linux 的 OOM killer 等) | `status=error`(在原因中说明) |

## 5. 与 `run_benchmarks.sh` 的对应

- 配置档 → 作业列表由 `scripts/bench/plan.py` 生成(`python -m bench.plan <profile> --order-seed S` 输出 TSV)。
  列:`build, graph_type, n, avg_degree, shape_args, weight_dist, weight_type, algorithm, trials, warmup, sources, seed`。
  `shape_args` 为 `--avg-degree 4` / `--degree 4` / `--attach 2` / `-`(无)。
- 作业顺序已由 order seed 打乱。1 个作业 = 1 个进程(`--isolate`)。
- 运行器不存在时:`run_benchmarks.sh` 输出"sssp_bench_runner が見つかりません"(找不到 sssp_bench_runner)并以退出码 2 结束。
  仅在 `smoke --allow-missing-runner`(用于 CI)时输出"未実装"(尚未实现)并以退出码 0 结束。
