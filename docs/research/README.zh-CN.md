# 论文分析文档(docs/research/)

[English](README.md) | [日本語](README.ja.md) | **简体中文**

研究对象是 arXiv 2611.34911v1 "Linear Single-Source Shortest Path, Almost"。论文本身不可再分发,请用 `scripts/fetch_paper.sh` 将其下载到 `.cache/paper/`(不受 git 跟踪)。本目录中的文档是对论文的转述与说明,不复制其原文。行号 `l.NNN` 指的是 `.cache/paper/2611.34911v1.txt` 中的行号。

## 索引

| 文档 | 范围 | 内容 |
|---|---|---|
| [01_preliminaries_primitives.zh-CN.md](01_preliminaries_primitives.zh-CN.md) | Sec 1–3 | 主定理、符号一览、阈值的使用方式、基本情形的约定、计算模型(长度运算、整数、随机数)、double 会出问题的地方、Sec 3 的原语 P1–P13(搜索、排序、选择、Dijkstra、Fibonacci 堆、生长路径与 MST、union-find、LCA、level ancestor、split–findmin)、与参考实现的对应、模块方案 |
| [02_ordered_hierarchy_local.zh-CN.md](02_ordered_hierarchy_local.zh-CN.md) | Sec 4–5 | 有序层次结构(Def 4.1–4.2)、规范化、子节点顺序、Scan 与 Evaluate(Alg 2–3)、可检查的不变量 I1–I8、owner 分解、radial 状态表示、阈值分量、局部层次结构(Alg 4)、原子替换、降低阈值时会失效的引理、模块方案与实现顺序 |
| [03_pivots_aux_outer.md](03_pivots_aux_outer.md) | Sec 6–8 | 枢纽路由证书(Alg 5–6、Lemma 6.2–6.8)、辅助度量与 owner 坐标、宽度类、日历与原子、OwnerDistances(Alg 7)、外层状态图与数值尺度、ShortestPaths(Alg 8)与主定理运行时间的分解、中间状态不变量一览 |
| [04_outer_scales_pipeline.md](04_outer_scales_pipeline.md) | Sec 8.2–8.3 | 外层尺度(Lemma 8.4)与外层流水线(Alg 8 l.8–17、Lemma 8.5)的实现规格:以精确的(值, 指数)形式表示的因子 F、逐阶段的 Lemma 8.4、Alg 8 的各阶段及回退的记录位置、Mode A / Mode B、现有组件所欠缺之处(M-1–M-12)、提议的问题 Q-65–Q-74、提议的验收测试 SC-1–14 与 PL-1–12 |
| [thresholds.zh-CN.md](thresholds.zh-CN.md) | 阈值与常数 | 所有常数的一览与分类(影响正确性还是只影响时间)、满足各阈值的最小 N、现实规模 N 下所走的路径、Mode B 中改动数值时会失效之处、对实现的建议 |
| [prior_work_and_oss.md](prior_work_and_oss.md) | 先行研究与现有 OSS | 现有实现(包括目标论文的一个实现)的许可证、验证方法与测量条件,先行研究,关于引入还是重新实现的建议 |
| [open-questions.md](open-questions.md) | 整篇论文 | 将上述文档中的含糊之处汇总为带 ID(Q-xx)的一张表:暂定解释、受影响的模块、解决状态(Decided by interpretation / Needs agreement / Open) |
| [../PAPER_MAPPING.md](../PAPER_MAPPING.md) | 整篇论文 | 论文中每个条目到模块、测试、Mode、状态与偏差(D-xx)的对应表。统一后的模块名记录也在这里 |

## 阅读顺序

1. **`thresholds.md` 的 §2–§3** — 先建立这样的认识:"在论文的阈值下,外层构造在任何现实输入上都不会运行,一切都落入基于扫描的 Dijkstra 的基本情形"。Mode A / B / C 划分(ADR 0003)与基本情形的两个版本(ADR 0008、D-1)的理由都在这里。
2. **`01_preliminaries_primitives.md` 的 §0 与 §2** — 计算模型与记号。编写任何模块的人都应阅读。
3. **`../PAPER_MAPPING.md`** — 查看你负责的模块所在的行、命名统一(§0)以及偏差一览(§3)。
4. 你所负责章节的分析文档:
   - Sec 3 原语 → `01` 的 §3 与 §5
   - Sec 4 求值器 → `02` 的 §2
   - Sec 5 局部层次结构 → `02` 的 §3–§4
   - Sec 6–8 → `03`(数值量级见 `thresholds.md` 的 §5–§7);Sec 8.2–8.3 的实现规格 → `04`
5. **`open-questions.md`** — 按你负责的模块所在的列筛选,在实现之前先就 "Needs agreement" 的行达成一致。新的含糊之处连同其位置(节、引理、行号)和候选解释一起添加到这里。
6. `prior_work_and_oss.md` — 在决定比较对象与测量条件时阅读(基准测试负责人)。

## 编写规则

- 将论文所述内容、由其推出的内容和你自己的解释分开(各文档中的 [Supplement]、[Interpretation] 与 [Ambiguous] 标记)。
- 不要复制长段原文。仅限引用短语,并引用节、引理与算法行号。
- 不要写绝对路径或个人信息(这是公开的 OSS 项目)。推送之前请运行 `uv run python scripts/check_public_safety.py`。
