# 0007. 不使用 fast-math;正确性验证保持 IEEE-754 语义

[English](0007-no-fast-math.md) | [日本語](0007-no-fast-math.ja.md) | **简体中文**

- 状态: Accepted
- 日期: 2026-11-27

## 背景 (Context)
`double` 权重的正确性检查,是在有文档记录的容差下,与精确的有理数参考结果比较。fast-math 的重新结合和有限数假设,会同时破坏这些比较和对无穷大的处理 (`WeightTraits<double>::infinity()`)。

## 决定 (Decision)
`cmake/Warnings.cmake` 为每个目标添加 `-fno-fast-math`。在用于正确性验证的构建中,绝不启用 `-ffast-math`、`-Ofast` 和 `-funsafe-math-optimizations`。若基准测试构建要使用这类标志,就必须同等地应用于所有算法,通过完整的测试套件,并在报告中注明。

## 后果 (Consequences)
- 浮点结果在不同编译器之间,可在 IEEE-754 舍入的范围内复现。
- 有些向量化的机会没有被利用。
