# アーキテクチャ決定記録 (ADR)

[English](README.md) | **日本語** | [简体中文](README.zh-CN.md)

決定ごとに 1 ファイルとし、名前は `NNNN-short-title.md` の形にします。番号は連番で、振り直しません。
決定を変えるときは、新しい ADR を作り、古い ADR を "Superseded by NNNN" とします。

テンプレート:

```markdown
# NNNN. Title
- Status: Proposed | Accepted | Superseded by NNNN
- Date: YYYY-MM-DD
## Context
## Decision
## Consequences
```

| ADR | タイトル | 状態 |
|---|---|---|
| [0001](0001-common-csr-representation.ja.md) | すべてのアルゴリズムで共通の CSR グラフ表現 | Accepted |
| [0002](0002-per-task-tdd-evidence.ja.md) | タスクごとの TDD 証跡ファイル | Accepted |
| [0003](0003-three-implementation-modes.ja.md) | 論文忠実・論文実験・実用の 3 モードの分離 | Accepted |
| [0004](0004-pinned-dependencies.ja.md) | 依存の固定(CMake は SHA-256、Python は uv.lock) | Accepted |
| [0005](0005-paper-not-redistributed.ja.md) | 論文は再配布しない | Accepted |
| [0006](0006-diagnostics-compiled-out-for-timing.ja.md) | 計測用ビルドでは診断をコンパイルで除く | Accepted |
| [0007](0007-no-fast-math.ja.md) | fast-math を使わない。正しさの検証では IEEE-754 の意味を保つ | Accepted |
| [0008](0008-mode-a-outer-base-case-variants.ja.md) | Mode A の外側の基底ケース:忠実なスキャン版 Dijkstra と、逸脱 D-1 としての二分ヒープ版 | Accepted |
