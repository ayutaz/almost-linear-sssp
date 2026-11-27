# Architecture Decision Records

**English** | [日本語](README.ja.md) | [简体中文](README.zh-CN.md)

One file per decision: `NNNN-short-title.md`, numbered sequentially, never renumbered.
A decision is changed only by a new ADR that marks the old one "Superseded by NNNN".

Template:

```markdown
# NNNN. Title
- Status: Proposed | Accepted | Superseded by NNNN
- Date: YYYY-MM-DD
## Context
## Decision
## Consequences
```

| ADR | Title | Status |
|---|---|---|
| [0001](0001-common-csr-representation.md) | Common CSR graph representation for all algorithms | Accepted |
| [0002](0002-per-task-tdd-evidence.md) | Per-task TDD evidence files | Accepted |
| [0003](0003-three-implementation-modes.md) | Separate paper-faithful, paper-experimental and practical modes | Accepted |
| [0004](0004-pinned-dependencies.md) | Pinned dependencies (CMake by SHA-256, Python by uv.lock) | Accepted |
| [0005](0005-paper-not-redistributed.md) | The paper is not redistributed | Accepted |
| [0006](0006-diagnostics-compiled-out-for-timing.md) | Diagnostics compiled out in timed builds | Accepted |
| [0007](0007-no-fast-math.md) | No fast-math; IEEE-754 semantics for correctness | Accepted |
| [0008](0008-mode-a-outer-base-case-variants.md) | Mode A outer base case: faithful scan Dijkstra, binary-heap variant as deviation D-1 | Accepted |
