# 文档导读

[English](README.md) | [日本語](README.ja.md) | **简体中文**

本页说明 `docs/` 下有哪些内容、各部分面向哪些读者,以及建议的阅读顺序。
项目概览请见[顶层 README](../README.zh-CN.md)。

部分文档目前只有日语版;如果某份文档同时有两种语言版本,日语版是与英文版放在同一目录下的 `.ja.md` 文件。文件名、节号、命令和标识符在所有语言版本中都相同。

## 从哪里开始

| 如果你想…… | 按以下顺序阅读 |
|---|---|
| 了解项目目前的发现 | [../README.zh-CN.md](../README.zh-CN.md) → [research/thresholds.md](research/thresholds.md) §2–§3 → [ADR 0003](adr/0003-three-implementation-modes.md) → [ADR 0008](adr/0008-mode-a-outer-base-case-variants.md) |
| 查看哪些已完成、接下来做什么 | [PROGRESS.md](PROGRESS.md) → [TASKS.md](TASKS.md) |
| 实现论文的某个部分 | [research/README.md](research/README.md) → 你负责的那一节的笔记 → [PAPER_MAPPING.md](PAPER_MAPPING.md) → [research/open-questions.md](research/open-questions.md) → [tdd/README.md](tdd/README.md) |
| 开发图库、基线或验证工具 | [design/graph_library.md](design/graph_library.md)、[design/baselines.md](design/baselines.md)、[design/numerics.md](design/numerics.md) |
| 运行或评审基准测试 | [design/benchmark_protocol.md](design/benchmark_protocol.md) → [design/benchmark_schema.md](design/benchmark_schema.md) → [design/benchmark_runner_cli.md](design/benchmark_runner_cli.md) → [research/prior_work_and_oss.md](research/prior_work_and_oss.md) |
| 修改 CI 或工具 | [design/ci.md](design/ci.md),然后阅读其中引用的 ADR |
| 贡献代码 | [../CONTRIBUTING.zh-CN.md](../CONTRIBUTING.zh-CN.md) → [tdd/README.md](tdd/README.md) → [TASKS.md](TASKS.md) |

## 各部分的内容

### `research/`:解读论文

用我们自己的话逐节解读论文(arXiv 2611.34911v1)的笔记:定义、引理、算法和常数,它们之间的依赖关系,
以及实现时需要注意的地方。论文本身不存放在仓库中;运行 `make fetch-paper` 即可在本地阅读。
`l.284` 这样的行号指的是该脚本提取出的文本中的行号。

- **读者:** 实现或评审论文组件的人,以及想核实我们对论文所作论断的人。
- **索引和阅读顺序:** [research/README.md](research/README.md)。
- **最重要的笔记:** [research/thresholds.md](research/thresholds.md):论文的各个构造从多大的规模开始
  生效,以及哪些常数影响正确性、哪些只影响运行时间。
- **相关工作:** [research/prior_work_and_oss.md](research/prior_work_and_oss.md):相关算法和现有的开源实现,
  附许可证信息。

### `PAPER_MAPPING.md`:论文 ↔ 代码

论文中的每个定义、引理、算法和关键常数各占一行,列出实现它的模块、对应的测试、所属模式(A 或 B)和状态。
文件中还列出了与论文的每一处**偏差**(`D-xx`)。修改 `src/paper/` 或
`src/experimental/paper_experimental/` 之前请先阅读它,并在每次修改时同步更新。

### `design/`:代码如何构建

不直接取自论文的那些部分的设计文档。

| 文档 | 主题 |
|---|---|
| [graph_library.md](design/graph_library.md) | 图的表示、文本格式、可复现的生成器、`sssp_gen` |
| [baselines.md](design/baselines.md) | Dijkstra 各变体与堆、平局处理、操作计数器 |
| [numerics.md](design/numerics.md) | 比较–加法模型与 IEEE-754 `double` 的区别;`double` 结果使用的容差及其依据 |
| [benchmark_protocol.md](design/benchmark_protocol.md) | 基准测试的执行方式:构建、重复次数、计时范围、验证、统计、测试套件 |
| [benchmark_schema.md](design/benchmark_schema.md) | 基准测试结果记录的格式 |
| [benchmark_runner_cli.md](design/benchmark_runner_cli.md) | C++ 基准测试运行器的命令行接口 |
| [ci.md](design/ci.md) | GitHub Actions 工作流和环境准备脚本 |

- **读者:** 修改相应部分的贡献者,以及评审者。

### `adr/`:决策记录

架构决策记录(Architecture Decision Records):每个会约束后续工作的决策各写成一个简短的文件,说明其背景
和影响。决策从不被直接删改,而是由新的 ADR 取代。索引见 [adr/README.md](adr/README.md)。建议先读
[0003](adr/0003-three-implementation-modes.md)(三种模式)和
[0008](adr/0008-mode-a-outer-base-case-variants.md)(忠实于论文的基础情形)。

- **读者:** 在提出修改之前想知道"为什么要这样做"的任何人。

### `TASKS.md` 与 `PROGRESS.md`:进度跟踪

- [TASKS.md](TASKS.md) 是唯一的工作项列表:ID、优先级、状态、依赖、验收标准、测试证据和完成时的提交。
- [PROGRESS.md](PROGRESS.md) 是简短的现状快照:当前阶段、已实现的内容、最新的测试和基准测试状态、已知问题
  和接下来的任务。恢复工作时请先读它。

这两个文件由负责协调的维护者维护;其他贡献者提出修改建议,而不是并行地直接编辑。

### `tdd/` 与 `TDD_LOG.md`:测试先行的证据

每个功能或缺陷修复都有一个任务文件 `tdd/tasks/<ID>.json`,记录失败的测试运行(RED)、通过的运行(GREEN)
和重构后的运行,包括命令、退出码和输出。格式见 [tdd/README.md](tdd/README.md)。`TDD_LOG.md` 是由
`scripts/tdd.py log` 根据这些文件生成的日志,不要手工编辑(在第一次生成之前它并不存在)。

- **读者:** 想确认某个改动确实是测试先行的评审者。

### `reports/`:结果

根据已保存的运行记录生成的基准测试报告和评审报告。目前为空:**尚未进行任何测量**。这里出现的每一个数字,
都必须来自保存在 `results/` 下的运行记录。

## 约定

- 引用论文时使用节号、引理编号、算法编号和行号,绝不大段摘录原文。
- 测量值与估计值绝不混用;凡是未测量的内容一律写作"未测量"。
- 不写本地绝对路径或个人信息:这是一个公开仓库。发布修改之前请运行 `make public-safety`。
