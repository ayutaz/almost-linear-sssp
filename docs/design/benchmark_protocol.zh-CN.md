# 基准测试测量协议

[English](benchmark_protocol.md) | [日本語](benchmark_protocol.ja.md) | **简体中文**

本文档规定将 arXiv 2611.34911v1(无向、非负权重 SSSP,期望 O(m + n log log log n))的实现与
Dijkstra 系基线进行比较时的测量步骤。这里**只写步骤和设定值**,完全不包含
实测数值。`results/` 和 `docs/reports/` 中只存放由实测生成的结果。

相关文档:

- 记录的 schema:[benchmark_schema.zh-CN.md](benchmark_schema.zh-CN.md)
- 运行器的 CLI:[benchmark_runner_cli.zh-CN.md](benchmark_runner_cli.zh-CN.md)
- 图生成器(图类型和权重分布的名称与定义):[graph_library.md](graph_library.md)

实现:`scripts/run_benchmarks.sh`(执行)、`scripts/bench/plan.py`(作业列表)、
`scripts/bench/envinfo.py`(环境)、`scripts/bench/schema.py`(验证)、`scripts/bench/stats.py`(统计)、
`scripts/analyze_results.py`(汇总与图表,`make report`)。

---

## 1. 目的与范围

### 1.1 目的

- 以可复现的方式测量:在相同的图、相同的权重、相同的内存布局(`include/sssp/graph.hpp` 的共用 CSR)上,
  论文算法在哪些规模、哪些图类型上比 Dijkstra 系更快或更慢。
- 将时间、内存和操作计数器分开记录,并为每一项主张固定作为依据的构建和指标。

### 1.2 比较对象(`algorithm`)

| `algorithm` | 含义 | `mode` | 支持的权重类型 |
|---|---|---|---|
| `dijkstra-binary` | 二叉堆 Dijkstra(Speedup 的基准) | `baseline` | `uint64`、`double` |
| `dijkstra-4ary` | 4 叉堆 Dijkstra | `baseline` | `uint64`、`double` |
| `dijkstra-fibonacci` | Fibonacci 堆 Dijkstra | `baseline` | `uint64`、`double` |
| `dijkstra-radix` | radix 堆 Dijkstra | `baseline` | **仅 `uint64`** |
| `paper-faithful` | Mode A 的入口(Alg 8)。外层阈值 N > L_N^64 在实际规模下不成立,因此基底情形(l.5)用扫描型 Dijkstra(Lemma 3.2)求解。记录为 fallback | `paper-faithful` | `uint64`、`double` |
| `paper-faithful-d1` | 同上。基底情形用二叉堆 Dijkstra 求解(偏离 D-1,ADR 0008) | `paper-faithful` | `uint64`、`double` |
| `paper` | 论文算法的外层构造(Alg 8 l.8-15)。保留名称,尚未实现 | `paper-faithful` | `uint64`、`double` |
| `paper-experimental` | 在论文框架内试验性地改变了部件和参数的版本 | `paper-experimental` | `uint64`、`double` |
| `paper-optimized` | 以实用速度为目标的优化版。保留名称,尚未实现 | `practical-optimized` | `uint64`、`double` |
| `practical-local` | Mode C 的局部 SSSP(`practical::local_sssp`)。参考论文(paper-inspired)的实现。见下文 | `practical-optimized` | **仅 `uint64`** |

`mode` 的含义:

- `baseline`:Dijkstra 系。
- `paper-faithful`:完全按论文的描述(常数和参数也与论文一致)。用于验证渐近复杂度的主张。
- `paper-experimental`:保持论文的框架,同时试验性地改变参数或部件。
- `practical-optimized`:以实用速度为目标进行了优化。与论文的差异必须留下记录。

`dijkstra-radix` × `double` 的作业不会由计划(`scripts/bench/plan.py`)生成。它不在对象范围内,也不作为失败记录。
`practical-local` × `double` 同样处理(运行器以退出码 3 拒绝;不替换为其他算法)。

关于 `paper-faithful` / `paper-faithful-d1`:

- `paper_sssp` 在一次调用中,无法从外部将零权重边的收缩(l.2)、基底情形(l.5)和树的恢复(l.6)分开。
  整体作为 `sssp` 测量,`preprocess` 和 `tree` 为 null(§4)。
- 每次都会有 1 条 `FallbackRecord`(收缩后的顶点数 N 和边数)原样放入记录的 `fallback`
  (即使在诊断 OFF 的构建中也是如此)。报告时间时,必须并列给出这个 fallback 的占比。
  这些时间是"论文的预处理 + Dijkstra"的时间,不是论文外层构造的时间。
- `double` 的验证使用 `Tolerance::for_path_hops(2(n − 1))`(numerics.md §3.2)。整数权重要求严格一致。
- `paper-faithful` 的扫描型基底情形为 Θ(N²),因此运行器以退出码 3 拒绝 n > 2^17
  (不替换为其他算法)。计划不会生成这样的组合,并将被排除的数量和原因输出到标准错误
  (`plan.VERTEX_LIMITS`)。`paper-faithful-d1` 没有上限。

关于 `practical-local`(任务 P5-PRACTICAL-LOCAL-BENCH):

- 这是参考论文(paper-inspired)的实现,不是论文算法本身。它没有外层构造(Sec 6 到 8),而是
  将收缩后的整个连通分量交给一个局部层次结构(Alg 4 的部件)进行评估。**不主张 Theorem 1 的复杂度。**
- 仅登记 `uint64`。对于一般的 `double`,树的恢复(Lemma 8.5)可能抛出 `invalid_argument`,因此
  不登记 `double` 的解法(`--weight-type double` 为退出码 3)。一般 `double` 的树恢复尚不支持。
- 若收缩后的顶点数 N 不超过 16(默认值 `direct_search_max_vertices = 16`;这不是通过实测调整的值),
  则用二叉堆 Dijkstra 直接求解,并记录 1 条 `FallbackRecord`(N = 1 时因另一个原因记录 1 条)。
  N 更大时构建层次结构,fallback 为 0 条(记录的 `fallback.count = 0`)。
- 与 `paper_sssp` 相同,将一次调用作为 `sssp` 测量,`preprocess` 和 `tree` 为 null(§4)。
- fallback 的占比按原因分别读取:`fallback.records[i].vertices / n`(N 是收缩后的,因此 ≤ n)和
  `fallback.records[i].edges / m`。不跨原因相加(Q-69)。该算法最多只有 1 条记录,因此
  `fallback.vertices` / `fallback.edges` 与这 1 条一致。报告时间时,必须并列给出这个占比。
- 性能为 **not measured**(本任务中未测量)。此表中没有任何数字。

### 1.3 权重类型(`weight_type`)

名称与 `sssp_gen --weight-type` 相同。

- `uint64`:整数权重。距离被精确计算(检测到溢出时抛出异常)。
- `double`:浮点权重。受 IEEE-754 舍入影响。比较带有容差(§5)。

---

## 2. 构建

分别使用 `CMakePresets.json` 的 2 个 preset。

| Preset | `SSSP_DIAGNOSTICS` | 用途 |
|---|---|---|
| `release` | `OFF` | 用于计时。公开的时间只取自此构建 |
| `release-diag` | `ON` | 用于记录操作计数器(堆操作次数、松弛次数等)和 fallback |

```sh
cmake --preset release      && cmake --build --preset release
cmake --preset release-diag && cmake --build --preset release-diag
```

规则:

- **关于时间的主张只来自 `release`。** `release-diag` 包含计数器的更新,因此其时间不用于公开的数值。
- **关于计数器的主张只来自 `release-diag`。** 在 `release` 中不记录计数器(null)。
- 每条记录都要记录构建 preset 名称、编译器及其版本、`CMAKE_BUILD_TYPE`。

---

## 3. 试验设计

| 项目 | 默认值 | 备注 |
|---|---|---|
| 预热 | standard 为 1 次,full 为 2 次(smoke 为 0 次) | 作为 `warmup=true` 保留记录,但不计入统计 |
| 重复次数(开发时) | 5 次(standard) | 用于初期确认和开发中的比较 |
| 重复次数(公开用) | 20 次以上(full) | 公开的数值只能采用满足此条件的 |
| `release-diag` 的重复 | smoke 1、standard 3、full 5 次 | 用于计数器。不用于时间 |
| 执行顺序 | 作业(1 个进程)的顺序用 order seed 打乱 | `BENCH_ORDER_SEED`(默认 12345)。1 个进程运行多个算法时,运行器按每次试验用 `--order-seed` 打乱 |
| 起点 | smoke 1、standard 2、full 4 个(`--sources`) | 由 `--seed` 确定性地选出 |

- 打乱执行顺序,是为了避免热量、缓存、频率在时间上的偏差
  集中在某个特定算法上。为了 Peak RSS 而让 1 个进程只运行 1 个算法(§6),因此打乱的单位是
  作业(构建 × 图 × 算法)。
- 对相同的 (graph, seed),所有算法使用相同的 k 个起点。1 次试验 = 1 个起点 × 1 个算法的 1 次运行。

---

## 4. 时间的分离

时刻使用 `std::chrono::steady_clock` 获取,单位一律为纳秒 [ns]。

| JSONL 的键(`times_ns.*`) | 名称 | 范围 |
|---|---|---|
| `graph_gen` | `graph_gen_ns` | 图的生成或读取 |
| `preprocess` | `preprocess_ns` | 算法特有的预处理。不适用时为 null |
| `sssp` | `sssp_ns` | SSSP 主体(直到距离确定) |
| `tree` | `tree_ns` | 最短路径树(parent 数组)的构建。对于无法与主体分离的实现为 null,并加以注明 |
| `end_to_end` | `end_to_end_ns` | 预处理 + 主体 + 树。**不含图生成** |
| `verify` | `verify_ns` | 正确性验证。必须在测量区间之外进行,不计入 `end_to_end` |

规则:

- 未测量的值为 **null**。不得填 0(否则无法与"测量过,结果为 0"区分)。
- 在 `tree` 为 null 的实现中,树的构建时间包含在 `sssp` 中。比较时要注明这一差异。
- 验证、日志输出和 JSON 写出不放入测量区间。

---

## 5. 正确性验证

验证在测量区间之外进行(另外记录为 `verify_ns`)。

1. **距离一致**:在所有顶点上与参考实现的距离进行比较。
   - 参考实现为 `dijkstra-binary`。对小规模的图也使用 `bellman-ford`,以检验参考实现本身。
   - `uint64`:完全一致。
   - `double`:带容差的比较(同时使用相对误差和绝对误差)。容差值遵循 `docs/design/numerics.md`。
   - 不可达的顶点在双方都必须是 infinity。
2. **树的一致性**:对于起点以外的每个可达顶点 v,parent 边 (parent[v], v) 必须真实存在,
   `dist[parent[v]] + w(parent[v], v) = dist[v]` 必须成立(`double` 带容差),并且沿 parent 追溯必须回到起点。
3. 验证失败的记录标为 `status=error`,并从时间统计中排除。

---

## 6. Peak RSS(内存)

- 1 个进程只运行 1 个 (algorithm, graph)(相当于运行器的 `--isolate`)。
  `scripts/run_benchmarks.sh` 为每个 algorithm 启动单独的进程。这样可以避免混入前一个算法分配的内存。
- 在进程退出前获取 `getrusage(RUSAGE_SELF).ru_maxrss`。
- 单位因 OS 而异。**macOS 为 bytes,Linux 为 KiB。** 在 Linux 上乘以 1024 换算为 bytes,记录到 `peak_rss_bytes`。
- 图构建完成后立即取得的值也记录为 `rss_after_graph_bytes`,并将 `peak_rss_bytes - rss_after_graph_bytes`
  作为算法自身增量的参考值。
  - 注意:`ru_maxrss` 是单调递增的最大值,因此这个差值不是严格的增量,只是参考值
    (例如图生成的临时空间造成了峰值时,差值会偏小)。

---

## 7. 关于 Apple Silicon(macOS)的注意事项

### 7.1 无法固定 CPU

- Apple Silicon 有 P 核心(性能)和 E 核心(能效)。数量可以用下面的命令确认。
  ```sh
  sysctl hw.perflevel0.physicalcpu   # P 核心数
  sysctl hw.perflevel1.physicalcpu   # E 核心数
  ```
- macOS 的调度器有时会把测量进程放到 E 核心上。
- **在 macOS 上没有能保证固定核心(CPU pinning)的手段。** `taskpolicy` 等只是优先级或 QoS 的指示,
  并不是固定到特定核心。因此,macOS 上的测量值应视为包含这种波动。

### 7.2 波动的原因

- 热节流(长时间连续运行导致频率下降)。
- 电源:AC 还是电池,是否开启低电量模式。可用 `pmset -g`、`pmset -g batt` 确认。
- 后台处理:Spotlight 的索引建立、Time Machine、浏览器、同步客户端等。

### 7.3 对策

- 连接 AC 电源。关闭低电量模式。
- 退出其他应用。如果可能,安排在 Spotlight 索引建立和 Time Machine 备份没有运行的时间段。
- 在试验之间加入冷却(等待时间)。
- 可记录的状态由 `scripts/bench/envinfo.py` 为每个 run 记录(在每条记录和
  `results/raw/<run_id>.<preset>.env.json` 中;结束时的状态在 `<run_id>.<preset>.env-end.json` 中)。
  - 电源状态(AC/电池、低电量模式)、`pmset -g therm` 的 `CPU_Speed_Limit`、CPU 名称、
    OS 及版本、P/E 核心数、load average(`os.getloadavg()`,与 `uptime` 相同的值)。
- 开始时或结束时的 `cpu_speed_limit` 低于 100 的 run,视为受到热限制而无效。
- 无法获取开始时或结束时的环境信息(`env.json` / `env-end.json`)的 run,由于无法做出上述判断,
  视为无效。`scripts/run_benchmarks.sh` 会将获取失败输出到标准错误,在最后的
  `done:` 行中以 `env_end_failed=<数>` 记录结束时的失败数,并以退出码 1 结束(开始时的失败当场以退出码 1 结束)。
  已测量的记录不删除而保留,但该 run 的值不用于公开。
  - 获取成功但 `cpu_speed_limit` 为 `null`(Linux 等没有获取手段的环境)时,不算失败。
    在报告中写明该 run 无法做热方面的判断。
- 对于记录的 load average 较高的 run,逐个判断是否用于公开的数值,不使用时要留下理由。

### 7.4 Linux 的情况

- 可以用 `taskset -c <core>` 固定核心。
- 可以用 `cpupower frequency-set -g performance` 固定频率调节器(在有权限的情况下)。
- `envinfo.py` 不会自动获取这些设置。使用了的话,要写在 run 的备注(报告)中。

---

## 8. 超时与内存上限

### 8.1 超时

`scripts/run_benchmarks.sh` 对每次调用施加外部超时。

| 配置档 | 默认超时(每个作业) |
|---|---|
| `smoke` | 60 秒 |
| `standard` | 600 秒 |
| `full` | 3600 秒(`make benchmark-full` 传入 Makefile 的 `BENCH_TIMEOUT`,默认 3600) |

设置环境变量 `BENCH_TIMEOUT` 后,任何配置档都使用该值。

macOS 默认不带 GNU `timeout`,因此按以下顺序回退。

1. `timeout`(GNU coreutils)
2. `gtimeout`(Homebrew 的 coreutils)
3. `perl -e 'alarm ...; exec ...'`(macOS 的默认方式;超时的退出码为 142,GNU `timeout` 为 124)

### 8.2 内存上限

- 用 `ulimit -v`(环境变量 `BENCH_MEM_LIMIT_MB`)设置上限。
- 在 macOS 上无法设置 `RLIMIT_AS`(`ulimit -v` 因 `Invalid argument` 而失败)。这时输出警告并
  在无上限的情况下继续,结束后将 `peak_rss_bytes` 超过上限的记录作为警告显示(不改变 status)。
- 如果运行器因 `std::bad_alloc` 返回退出码 5,则 `status=memory_limit`。

### 8.3 失败的处理

- 超时、内存超限和异常终止(包括验证失败),分别保留为 `status=timeout` / `memory_limit` / `error` 的记录。
- 运行器没能写出任何内容时,`run_benchmarks.sh` 会追加 1 条失败记录(`python -m bench.schema failure`)。
- 这些从时间统计中排除,并按作业单独报告**失败率**(`results/summary/status.csv`)。
- 超时和内存超限在预期之内,因此 `run_benchmarks.sh` 的退出码为 0。只要有 1 个异常终止(`error`),就为 1。
- 不对失败的组合的值进行推测或插值。

---

## 9. 指标与统计

### 9.1 Speedup

```
Speedup = T(dijkstra-binary) / T(对象)
```

- 在相同的 (`graph_type`, `n`, `avg_degree`, `weight_dist`, `weight_type`) 组合内,使用各算法的 **median**。
- Speedup > 1 表示对象比 `dijkstra-binary` 更快。
- 默认的 T 是 `sssp_ns`。同时也输出 `end_to_end_ns` 版本。
- 基准和对象任一方没有成功测量(`status=ok`、非 warmup、`release` 构建)的组合,不给出 Speedup。
- `params` 不同的视为不同的系列,分别给出 Speedup。

### 9.2 统计量

| 统计量 | 处理 |
|---|---|
| median | 主指标 |
| mean | 参考 |
| p95 | 样本数 < 20 时标记为仅供参考 |
| stddev | 样本标准差(ddof=1)。只有 1 个样本时为 null |
| min / max | 参考 |

- p95 是顺序统计量的线性插值(numpy 的默认方式)。

- 预热(`warmup=true`)和失败记录不计入统计。

### 9.3 Speedup 的置信区间

- 用 bootstrap 求得。默认重采样次数为 10000,随机数种子为 20261127(固定),95% 百分位区间
  (`scripts/analyze_results.py --bootstrap-resamples / --bootstrap-seed`)。
- 基准侧和对象侧分别独立重采样,每次计算 median 之比。
- bootstrap seed 是汇总侧的设置,不进入记录。若更改默认值进行汇总,要在报告中写明。

---

## 10. 图的规模与配置档

### 10.1 规模网格

- n ∈ {1e3, 1e4, 1e5, 1e6, 1e7}
- 平均度数 d ∈ {2, 4, 8, 16, 32}(无向图中 d = 2m/n)

### 10.2 配置档

| 配置档 | n | d | 图类型 | 权重分布 | trials | warmup | diag trials | 起点 |
|---|---|---|---|---|---|---|---|---|
| `smoke` | {1e3, 1e4} | 4 | `er`、`grid` | `uniform_int` | 1 | 0 | 1 | 1 |
| `standard` | {1e3, 1e4, 1e5, 1e6} | {4, 16} | 全部类型 | `constant`、`uniform_int`、`wide_range` | 5 | 1 | 3 | 2 |
| `full` | {1e3, …, 1e7} | {2, 4, 8, 16, 32} | 全部类型 | 全部分布 | 20 | 2 | 5 | 4 |

- 定义在 `scripts/bench/plan.py` 的 `PROFILES` 中。可以用 `python -m bench.plan <profile> --order-seed S` 确认作业列表。
- 筛选:`BENCH_ALGORITHMS`、`BENCH_GRAPH_TYPES`、`BENCH_SIZES`、`BENCH_BUILDS`(逗号分隔)。

- `smoke` 用于确认能否正常运行,其结果不用于性能主张。
- 根据规模的不同,可能会达到时间或内存的上限。达到上限的组合记录为失败记录,不做插值。
- 度数由图族决定的(`grid`、`road`、`tree`、`adversarial_long_path`)不按 d 复制(`avg_degree` 为 null)。
- 权重由构造决定的(`road`、`adversarial_long_path`)不改变权重分布(`weight_dist=builtin`)。

---

## 11. 图类型(`graph_type`)

名称和定义遵循 `sssp_gen --type`(docs/design/graph_library.md §6)。

| `graph_type` | 含义 | d 的传递方式 | 备注 |
|---|---|---|---|
| `er` | 均匀随机 G(n, m),m = round(n·d/2) | `--avg-degree d` | |
| `regular` | 随机 d-正则 | `--degree d` | |
| `ba` | Barabási–Albert(优先连接) | `--attach d/2` | 平均度数约为 2·attach |
| `chung_lu` | Chung–Lu(幂律的期望度数) | `--avg-degree d` | 度数偏斜较大 |
| `grid` | 二维网格 ⌊√n⌋ × ⌊√n⌋ | 无(度数约 4) | n 会被取整为平方数 |
| `road` | 类道路网(网格的扰动 + 近距离捷径) | 无 | 权重为欧氏长度(`builtin`) |
| `tree` | 均匀随机树 | 无 | m = n − 1 |
| `adversarial_long_path` | 长路径 + 绕行边 | 无 | 深的最短路径树。权重由构造固定(`builtin`) |

- `adversarial_decrease_key` 的边数为 n² 量级,因此不放入规模网格。要使用时,以较小的 n 单独运行。
- 记录实际生成的 n 和 m(由于取整和重复边的去除,可能与指定值有偏差)。
- 连通性:从起点无法到达的顶点照常处理,并将可达的顶点数记录为 `reachable`。

---

## 12. 权重分布(`weight_dist`)

名称和定义遵循 `sssp_gen --weights`(docs/design/graph_library.md §5)。参数(`lo`、`hi`、`mean`、
`zero_probability`、`wide_log2`)在配置档中使用生成器的默认值。

| `weight_dist` | 含义 | 配置档中的权重类型 |
|---|---|---|
| `constant` | 所有边权重相同(默认 `lo`) | `uint64` |
| `uniform_int` | [lo, hi] 上的均匀整数 | `uint64` |
| `log_uniform` | 对数均匀 | `uint64` |
| `wide_range` | 均匀选取二进制指数(动态范围大) | `uint64` |
| `zero_mixed` | 以概率 `zero_probability` 取 0,其余为均匀 | `uint64` |
| `uniform_real` | [lo, hi) 上的均匀实数 | `double` |
| `exponential` | 指数分布 | `double` |
| `builtin` | 由生成器的构造决定的权重(`road`、`adversarial_long_path`) | `uint64` |

- 对于 `uint64` 下距离之和可能溢出的组合,若检测到溢出,则 `status=error`。

---

## 13. 可复现性

所有记录都要记录以下内容。

- 种子:图的 seed 和起点选择的 seed(`seed`)、`order_seed`、算法的随机数 seed(放入 `params`)。
  bootstrap seed 是汇总侧的默认值(§9.3)
- git 的 SHA 和 dirty 标志(工作树中是否有未提交的更改)
- 复现对象的论文版本:`2611.34911v1`
- 构建 preset、编译器、环境信息(§7.3)

规则:

- **`dirty=true` 的记录不用于公开的数值。** 只用于开发中的确认。

---

## 14. 输出与汇总

```
results/raw/<run_id>.jsonl          (已 gitignore,每行 1 条记录)
        │
        └─ make report  (scripts/analyze_results.py)
              ├─ results/summary/*.csv
              └─ docs/reports/figures/*.png
```

- `results/summary/` 的 CSV:
  - `timing_summary.csv`:每个 配置 × 算法 × params × 区间 的 samples / median / mean / p95 / stddev / min / max /
    p95_reference_only / dirty_records / run_ids
  - `speedup.csv`:Speedup 和 95% CI(`sssp` 和 `end_to_end`)
  - `memory.csv`:Peak RSS 的 median / max,以及图构建完成后立即取得的 RSS 的 median
  - `status.csv`:按作业统计的 ok / timeout / memory_limit / error 及失败率
  - `fallback.csv`:`release-diag` 中 fallback 的发生率、次数、顶点和边的占比、原因
  - `ops.csv`:操作计数器和 `counters` 的 median
- `docs/reports/figures/` 的图:`runtime_vs_n__<type>`、`runtime_vs_m__<type>`、`speedup_vs_n__<type>`、
  `memory_vs_n__<type>`、`fallback_rate__<type>`、`breakdown__<type>`(预处理/主体/树的明细)、
  `speedup_by_graph_type`、`speedup_by_weight_dist`、`param__<algorithm>__<param>`(阈值参数与时间)。
- 记录的字段定义见 [benchmark_schema.zh-CN.md](benchmark_schema.zh-CN.md),运行器的参数见
  [benchmark_runner_cli.zh-CN.md](benchmark_runner_cli.zh-CN.md)。
- 在图表中,不绘制未测量的组合和失败的组合。如果系列中途缺少某个点,则在此处断开折线(不插值)。
- 如果存在违反 schema 的记录,`make report` 会失败(用 `--skip-invalid` 可排除它们并继续)。
- 没有可汇总的数据时,输出 `データなし`(无数据),不编造数值。

---

## 15. 解读结果时的注意事项

- 渐近复杂度的差异(O(m + n log n) 对 O(m + n log log log n))在实际的 n 下,可能被常数倍和缓存效应掩盖。
  没有变快,不一定是实现有误,也不一定是论文有误。
- op 计数器(`release-diag`)是比较加法模型中复杂度的替代指标,不一定与实际时间成正比
  (内存访问模式和分支预测没有被计入)。
- `paper` 的期望复杂度是对随机化的期望。试验间的波动既包含随机化,也包含系统的波动。
- 根据图类型和权重分布,优劣可能发生变化。不要将特定组合的结果一般化。
- 报告中必须写明:macOS 上的结果是在没有固定 CPU 的条件下取得的。
