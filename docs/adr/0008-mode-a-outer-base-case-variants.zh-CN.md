# 0008. Mode A 外层基础情形:忠实的扫描式 Dijkstra,以及作为偏离 D-1 的二叉堆变体

[English](0008-mode-a-outer-base-case-variants.md) | [日本語](0008-mode-a-outer-base-case-variants.ja.md) | **简体中文**

- 状态: Accepted
- 日期: 2026-11-27

## 背景 (Context)
Algorithm 8 (Sec 8.3,伪代码第 4–7 行)在 `N <= L_N^64` 时返回直接搜索的结果。满足外层条件 `N > L_N^64` 的 `N` 的集合,最早在 `N = 589^64 + 1`(约 1.94e177)时变为非空,而且它不是向上封闭的:在 `[2^589 - 1, 590^64]` 上它再次为假,对所有更大的 `N` 都成立则要从 `590^64 + 1` 开始(见 `docs/research/thresholds.md`)。因此,对每个能放进内存的输入(`N < 2^64`,此时 `L_N^64 >= 2^64 > N`),Mode A 只执行基础情形。

论文把基础情形规定为 Lemma 3.2 的扫描式 Dijkstra 搜索,`O(N^2 + m)`。在 `N = 10^6` 时约为 `10^12` 次操作,因此忠实的基础情形无法在评估计划所要求的规模上做基准测试。Lemma 3.2 还给出了一种二叉堆式搜索,输出(精确的距离)相同,只是运行时间不同。

## 决定 (Decision)
- `paper-faithful`(Mode A,默认)完全按原文运行扫描式基础情形,并记录原因为 `"outer: N <= L_N^64 (paper base case, scan Dijkstra)"` 的 `FallbackRecord`。
- `paper-faithful-d1` 在基础情形中运行 Lemma 3.2 的二叉堆式搜索。它在 `docs/PAPER_MAPPING.md` 中列为偏离 **D-1**,并记录原因 `"outer: N <= L_N^64 (paper base case, binary-heap Dijkstra, deviation D-1)"`。
- 阈值判断用饱和运算即时求值(用饱和方式计算 `L_N^64` 并与 `N` 比较),绝不与预先算好的分界值比较,因为满足条件的集合不是向上封闭的。
- 基准测试要标明测量的是哪一种变体。两种变体都不会被报告为对论文非平凡构造的测量;该构造只在 Mode B 中被执行(`docs/adr/0003-three-implementation-modes.md`)。

## 后果 (Consequences)
- 大规模的 "paper-faithful" 计时使用 `paper-faithful-d1`,按其构造,它测的是一次 Dijkstra 搜索加上论文的预处理(零长度边的收缩、每个顶点对中最轻的边)。报告必须在每个这样的数字旁边注明这一点。
- 扫描版本仍会在小的 `N` 上测试正确性,以验证忠实的代码路径。
