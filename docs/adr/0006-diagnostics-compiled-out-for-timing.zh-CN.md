# 0006. 计时构建中编译掉诊断功能

[English](0006-diagnostics-compiled-out-for-timing.md) | [日本語](0006-diagnostics-compiled-out-for-timing.ja.md) | **简体中文**

- 状态: Accepted
- 日期: 2026-11-27

## 背景 (Context)
执行路径计数器和操作次数(`include/sssp/result.hpp` 中的 `Diagnostics`)用于验证论文的哪些构造运行过、回退发生得有多频繁,但它们会扰动计时。

## 决定 (Decision)
计数器通过 `SSSP_COUNT(...)` 更新,由编译定义 `SSSP_DIAGNOSTICS` 控制。`release` 预设把 `SSSP_DIAGNOSTICS` 设为 `OFF`,并且是计时唯一使用的预设。次数来自在相同输入上单独运行的 `release-diag`。`FallbackRecord` 始终记录(它们是结果的一部分,而不是计数器)。

## 后果 (Consequences)
- 计时运行不带计数器的开销。
- 报告把来自同一提交和同一输入的耗时 (release) 与次数 (release-diag) 结合起来。
