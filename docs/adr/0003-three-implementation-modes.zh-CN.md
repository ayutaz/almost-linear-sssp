# 0003. 分离论文忠实、论文实验和实用三种模式

[English](0003-three-implementation-modes.md) | [日本語](0003-three-implementation-modes.ja.md) | **简体中文**

- 状态: Accepted
- 日期: 2026-11-27

## 背景 (Context)
论文的阈值(外层 `N > L_N^64`,局部 `N > 2 L_N^4`,`L_N = ceil(log2(N+2))`,Sec 2)意味着,对任何 `N < 2^64`,外层构造都不会运行。因此,忠实的实现会把现实的输入都花在基础情形上。要研究构造本身,以及追求实用的速度,需要修改常数或使用不同的数据结构,而这绝不能与论文的算法混淆。

## 决定 (Decision)
- Mode A,`src/paper/`:论文忠实。阈值完全按原文,用防溢出的整数运算比较。进入基础情形会记录在 `Diagnostics.fallbacks` 中。
- Mode B,`src/experimental/paper_experimental/`:论文的结构,阈值参数化。每一处偏离都列在 `docs/PAPER_MAPPING.md` 中。无法确认正确性的输入视为失败;自动替换为参考算法,绝不报告为成功。
- Mode C,`src/experimental/practical/`:实用优化,每项都有标注。
- 基线算法放在 `src/baseline/`,不会被修改得像论文。

## 后果 (Consequences)
- 所谓"论文的算法",只指 Mode A。
- 报告中始终在耗时旁边显示回退的占比。
- 模式之间的代码共享通过明确的共享组件进行,而不是复制后再改。
