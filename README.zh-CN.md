# almost-linear-sssp

[English](README.md) | [日本語](README.ja.md) | **简体中文**

本项目独立复现、验证并评测论文
**"Linear Single-Source Shortest Path, Almost"**(作者 Iman Gholami、MohammadTaghi Hajiaghayi,
[arXiv:2611.34911v1](https://arxiv.org/abs/2611.34911))。

我们用 C++20 实现论文中的算法,用独立的参照实现(oracle)检验每一个结果,并计划与调优过的
Dijkstra 基线进行性能比较。本项目与论文作者没有任何关系。

> **状态:早期开发中。** 图库、Dijkstra 基线和验证工具已经完成,论文本身的各个组件正在实现中。
> **目前尚未进行任何基准测试**,因此本 README 不包含任何性能数据。

## 目录

- [当前状态](#当前状态)
- [一段话介绍论文](#一段话介绍论文)
- [为什么阈值很重要](#为什么阈值很重要)
- [快速开始](#快速开始)
- [算法](#算法)
- [如何检验正确性](#如何检验正确性)
- [基准测试](#基准测试)
- [仓库结构](#仓库结构)
- [文档](#文档)
- [相关工作](#相关工作)
- [引用](#引用)
- [许可证](#许可证)
- [免责声明](#免责声明)

## 当前状态

| 模块 | 状态 |
|---|---|
| 图库:CSR、邻接表、边列表;edge-list 与 DIMACS `.gr` 格式的读写 | 已完成 |
| 可复现的图生成器与命令行工具 `sssp_gen` | 已完成 |
| 基线:Bellman–Ford;基于二叉堆、4 叉堆、Fibonacci 堆和 radix 堆的 Dijkstra | 已完成 |
| 验证:最短路径树检查器、精确有理数参照实现、操作计数器 `Counted<W>` | 已完成 |
| 论文组件,第 2–5 节(基础原语、层次结构求值、局部层次结构) | 实现中 |
| 论文组件,第 6–8 节(枢纽点采样、owner 坐标、外层层次结构) | 计划中 |
| 差分测试框架与基准测试运行器 | 实现中 |
| 基准测试结果 | **未测量** |

详细的任务列表见 [docs/TASKS.md](docs/TASKS.md),当前进展摘要见
[docs/PROGRESS.md](docs/PROGRESS.md)。

## 一段话介绍论文

给定一个有 `n` 个顶点、`m` 条边、边长为非负实数的无向图,单源最短路径(SSSP)问题要求出从一个源点到
其余所有顶点的距离。使用 Fibonacci 堆的 Dijkstra 算法耗时 `O(m + n log n)`。论文声称在
*比较–加法模型*(comparison–addition model,只允许对边长做比较和加法,不能取整,也不能利用位运算技巧)
下,算法的期望运行时间为 **`O(m + n log log log n)`**。只要 `m` 不小于 `n log log log n`,这就是线性时间。
该算法是随机化的,但属于 *Las Vegas* 型:每一次终止的运行都会返回精确的距离,随机性只影响运行时间。
算法的组成部分包括:基于最小生成树把图划分成若干区域、每个区域设一个"枢纽点"(pivot);通过随机采样
为枢纽点之间的短路径给出证明;为每个顶点计算辅助坐标;以及一个沿着强连通分量的有序层次结构求出距离的过程。

## 为什么阈值很重要

这是到目前为止最重要的发现,它决定了整个项目的结构。

论文只在输入足够大时才使用完整的构造。对于有 `N` 个顶点的输入,记 `L_N = ⌈log2(N + 2)⌉`(第 2 节),
论文要求:

- **外层构造:** `N > L_N^64`;更小的输入进入运行 Dijkstra 的基础情形(base case);
- **局部构造**(作用于小的子块):`N > 2 · L_N^4`;更小的子块进入一个简单的 Kruskal 式基础情形。

我们用任意精度整数算出了精确的边界
([docs/research/thresholds.md](docs/research/thresholds.md)):

| 条件 | 首次成立的 `N` |
|---|---|
| 外层,`N > L_N^64` | `589^64 + 1` ≈ 1.9 × 10^177(在其上方不远处的一个短区间内又不成立;对所有 `N ≥ 590^64 + 1` ≈ 2.2 × 10^177 都成立) |
| 局部,`N > 2 · L_N^4` | 209,953(`2 · 18^4 = 209,952`) |

| `N` | `L_N` | 外层构造是否运行? | 局部条件是否成立? |
|---|---|---|---|
| 10^3 | 10 | 否(`L_N^64` ≈ 2^213) | 否(需要 `N > 20,000`) |
| 10^5 | 17 | 否(≈ 2^262) | 否(需要 `N > 167,042`) |
| 10^6 | 20 | 否(≈ 2^277) | 是(`N > 320,000`) |
| 10^8 | 27 | 否(≈ 2^304) | 是(`N > 1,062,882`) |
| 任意 `N < 2^64` | ≤ 65 | **否**,因为 `L_N^64 ≥ 2^64 > N` | — |

因此,严格按论文实现的代码**在任何能放进内存的输入上都只会运行基础情形(Dijkstra)**。局部构造只会在
外层构造内部被调用,所以它同样不会运行。对严格实现计时,测到的是 Dijkstra 加上一些预处理,而不是论文
提出的新方法。

为了如实地研究论文的构造,本项目把三种模式严格分开
([ADR 0003](docs/adr/0003-three-implementation-modes.md)):

| 模式 | 目录 | 内容 | 用途 |
|---|---|---|---|
| **A** paper-faithful | `src/paper/` | 阈值与论文完全一致,用不会溢出的整数运算判定。每次进入基础情形都会被记录。 | 按原样复现论文。在实际输入上运行的就是基础情形;变体 `paper-faithful-d1` 把基础情形中基于扫描的 Dijkstra 换成二叉堆版本(偏差 D-1,[ADR 0008](docs/adr/0008-mode-a-outer-base-case-variants.md))。 |
| **B** paper-experimental | `src/experimental/paper_experimental/` | 保留论文的结构,但使用更小、可配置的规模阈值。正确性证明所依赖的常数保持不变,或根据新的规模重新推导([thresholds.md §6](docs/research/thresholds.md))。 | 在实际规模上真正运行第 5–8 节的构造,并验证每一个结果。无法确认正确性的输入一律算作失败。 |
| **C** practical-optimized | `src/experimental/practical/` | 任何优化都可以,前提是结果正确并且明确标注。 | 追求实际速度。 |

只有模式 A 才会被称为"论文的算法"。今后报告的每一个计时结果旁边,都会注明由基础情形处理的顶点和边所占的比例。

## 快速开始

环境要求:macOS 或 Linux,支持 C++20 的编译器(开发和测试使用 Apple Clang 与 LLVM Clang),
CMake ≥ 3.24,Ninja,以及用于管理 Python ≥ 3.12 的 [uv](https://docs.astral.sh/uv/)。
CI 和贡献者还需要 clang-format 与 clang-tidy。

```bash
make bootstrap                 # 检查工具链,安装缺失的工具,运行 `uv sync`
                               # (./scripts/bootstrap.sh --check 只报告,不安装任何东西)
make build                     # 配置并构建 debug preset,输出到 build/debug/
make test                      # 运行 C++ 测试(ctest)和 Python 测试(pytest)

./build/debug/sssp_gen --list  # 列出图生成器和权重分布
./build/debug/sssp_gen --type er --n 1000 --avg-degree 4 --weights uniform_int --seed 1 --out g.txt

make fetch-paper               # 把论文下载到 .cache/paper/(永不提交到仓库)
```

同一个种子在任何平台上都会生成同一张图。其他构建 preset(`release`、`sanitizer`、`coverage`、`ci` 等)
见 [CMakePresets.json](CMakePresets.json),例如 `make build PRESET=release`。参与贡献之前,请先运行
`make verify`(见 [CONTRIBUTING.zh-CN.md](CONTRIBUTING.zh-CN.md))。

## 算法

| 算法 | 函数 / 位置 | 权重类型 | 状态 |
|---|---|---|---|
| Bellman–Ford(按轮松弛;用作 oracle) | `sssp::baseline::bellman_ford` | `uint64_t`、`double` | 已完成 |
| Dijkstra,二叉堆(惰性删除) | `sssp::baseline::dijkstra_binary_heap` | `uint64_t`、`double` | 已完成 |
| Dijkstra,4 叉堆(decrease-key) | `sssp::baseline::dijkstra_4ary_heap` | `uint64_t`、`double` | 已完成 |
| Dijkstra,Fibonacci 堆 | `sssp::baseline::dijkstra_fibonacci_heap` | `uint64_t`、`double` | 已完成 |
| Dijkstra,radix 堆 | `sssp::baseline::dijkstra_radix_heap` | 仅 `uint64_t` | 已完成 |
| 论文第 3 节:基础原语(MST、有序 union–find、LCA / level ancestor、split–findmin) | `src/paper/`(模式 A) | — | 实现中 |
| 论文第 4 节:有序层次结构的求值 | `src/paper/`(模式 A) | — | 实现中 |
| 论文第 5 节:局部层次结构 | `src/paper/`(模式 A) | — | 实现中 |
| 论文第 6–8 节:枢纽点路径、owner 坐标、外层层次结构、主算法 | `src/paper/`(模式 A) | — | 计划中 |
| 调整后的阈值 / 实用优化 | 模式 B / 模式 C | — | 计划中 |

所有基线共用同一个扫描循环,使用本仓库自己实现的堆,并统计相同的操作
([docs/design/baselines.md](docs/design/baselines.md))。论文中每个引理和算法与代码、测试的对应关系见
[docs/PAPER_MAPPING.md](docs/PAPER_MAPPING.md)。

## 如何检验正确性

正确性优先于速度。每个实现都由与其不共享代码的工具来检验:

- **独立的 oracle。** Bellman–Ford;它本身又在*所有*小规模带权简单图上、从每个源点出发,与
  Floyd–Warshall 做了比对(穷举枚举)。
- **最短路径树验证器**(`sssp::verify_sssp`)。它检查 `dist[source] = 0`、没有任何边还能缩短距离、
  每条树边都是紧的(tight)、沿父指针能回到源点,以及恰好是不可达的顶点距离为无穷大。这样无需信任算法本身
  就能证明结果正确。
- **精确有理数参照实现。** 针对浮点输入的 Python oracle,使用精确分数计算(`tests/python/reference/`)。
- **精确比较与带容差比较。** 整数权重必须完全一致。`double` 结果按
  [docs/design/numerics.md](docs/design/numerics.md) 中推导的舍入误差界进行比较;绝不为了让测试通过而放宽容差,
  并且所有构建都禁用 fast-math。
- **操作计数。** `Counted<W>` 包装权重类型,统计比较和加法的次数,以便检验是否遵守比较–加法模型。
- **回归测试。** 发现的每一个反例都会被最小化,并保存在 `tests/regression/` 中。
- **不允许静默回退。** 只要某个实现用直接搜索代替论文的构造来处理输入的一部分,就会把这件事记录到
  `Diagnostics.fallbacks` 中。

代码采用测试先行的方式编写,每个任务失败与通过的测试运行记录都保存在 `docs/tdd/tasks/` 中。

## 基准测试

**目前还没有结果。在有保存下来的运行记录作为依据之前,所有性能数据一律视为"未测量"。**

测量规程已事先确定([docs/design/benchmark_protocol.md](docs/design/benchmark_protocol.md)):
使用 Release 构建并在编译时去掉诊断代码;所有算法使用相同的 CSR 图、相同的输入、相同的机器、相同的预热和
重复次数;图的构建不计入时间;每个结果都必须先通过验证才算数;运行环境的元数据与结果一同保存;报告中位数、
平均值、p95、标准差和样本数;加速比以二叉堆 Dijkstra 为基准;论文算法的每个计时结果旁都注明基础情形所占的比例。

各测试套件将通过 `make benchmark-smoke`、`make benchmark` 和 `make benchmark-full` 运行,并用
`make report` 汇总。C++ 基准测试运行器仍在编写中,因此这些目标目前还不会产生结果。

## 仓库结构

| 路径 | 内容 |
|---|---|
| `include/sssp/` | 公共头文件:图类型(`CsrGraph`)、结果(`SsspResult`、`Diagnostics`)、基线、生成器、IO、验证 |
| `src/baseline/` | Bellman–Ford 和各种 Dijkstra 变体,以及它们使用的堆 |
| `src/graph/`、`src/verify/` | 图的读写与生成器;最短路径验证 |
| `src/paper/` | 模式 A,忠实于论文(实现中) |
| `src/experimental/` | 模式 B(`paper_experimental/`)和模式 C(`practical/`) |
| `tools/` | 命令行工具(`sssp_gen`) |
| `tests/` | C++ 测试(Catch2:unit、integration、differential、property、regression)和 Python 测试(pytest) |
| `benchmarks/`、`scripts/` | 基准测试驱动程序和脚本(环境准备、论文下载、基准测试运行、结果分析) |
| `docs/` | 研究笔记、设计文档、决策记录、任务跟踪 |

## 文档

请从 [docs/README.zh-CN.md](docs/README.zh-CN.md) 开始阅读,其中说明了 `docs/` 下每份文档的内容、面向的读者
以及阅读顺序。部分文档目前只有日语版;如果某份文档同时有两种语言版本,日语版是与英文版放在同一目录下的
`.ja.md` 文件。

## 相关工作

在本问题设定(无向图、非负实数权重、比较–加法模型)下,此前的最好结果包括 Duan、Mao、Shu 和 Yin
(FOCS 2023,`O(m √(log n · log log n))`,随机化)、Yan 对其的去随机化(SOSA 2026),以及 Kadria 和
Roditty 2026 年的预印本。经典的参照点包括使用 Fibonacci 堆的 Dijkstra、Thorup 针对整数权重的线性时间算法
(1999),以及 Pettie 和 Ramachandran(2005)。同一篇论文至少还有一个其他的公开实现。
[docs/research/prior_work_and_oss.md](docs/research/prior_work_and_oss.md) 列出了这些工作和现有的开源代码,
并注明了许可证以及它们与本项目的可比程度。

## 引用

如果使用本软件,请引用原论文,并在适当时引用本软件。

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

本软件自身的引用元数据见 [CITATION.cff](CITATION.cff)。

## 许可证

本仓库中的代码和文档以 [MIT 许可证](LICENSE) 发布。

本仓库**不**包含论文。arXiv 发布该论文所用的许可证不允许再分发,因此 `make fetch-paper` 会把论文下载到
已被 git 忽略的 `.cache/paper/` 中([ADR 0005](docs/adr/0005-paper-not-redistributed.md))。本仓库中的文档
通过节号、引理编号和算法编号来引用论文,而不是摘录原文。

## 免责声明

本项目是独立的复现工作,与论文作者没有关联,也未获得其认可。实现中或我们对论文理解上的任何错误,都由我们
负责。如果这里的结果与论文不一致,首先应当怀疑的是我们的复现;我们会在
[docs/PAPER_MAPPING.md](docs/PAPER_MAPPING.md) 和 [docs/adr/](docs/adr/) 中公开记录这类问题和偏差。
