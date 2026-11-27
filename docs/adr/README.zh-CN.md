# 架构决策记录 (ADR)

[English](README.md) | [日本語](README.ja.md) | **简体中文**

每个决策一个文件,命名为 `NNNN-short-title.md`,按顺序编号,不重新编号。
要更改某个决策,只能新建一份 ADR,并把旧的标记为 "Superseded by NNNN"。

模板:

```markdown
# NNNN. Title
- Status: Proposed | Accepted | Superseded by NNNN
- Date: YYYY-MM-DD
## Context
## Decision
## Consequences
```

| ADR | 标题 | 状态 |
|---|---|---|
| [0001](0001-common-csr-representation.zh-CN.md) | 所有算法共用的 CSR 图表示 | Accepted |
| [0002](0002-per-task-tdd-evidence.zh-CN.md) | 按任务划分的 TDD 证据文件 | Accepted |
| [0003](0003-three-implementation-modes.zh-CN.md) | 分离论文忠实、论文实验和实用三种模式 | Accepted |
| [0004](0004-pinned-dependencies.zh-CN.md) | 固定依赖(CMake 用 SHA-256,Python 用 uv.lock) | Accepted |
| [0005](0005-paper-not-redistributed.zh-CN.md) | 不再分发论文 | Accepted |
| [0006](0006-diagnostics-compiled-out-for-timing.zh-CN.md) | 计时构建中编译掉诊断功能 | Accepted |
| [0007](0007-no-fast-math.zh-CN.md) | 不使用 fast-math;正确性验证保持 IEEE-754 语义 | Accepted |
| [0008](0008-mode-a-outer-base-case-variants.zh-CN.md) | Mode A 外层基础情形:忠实的扫描式 Dijkstra,以及作为偏离 D-1 的二叉堆变体 | Accepted |
