# 0002. タスクごとの TDD 証跡ファイル

[English](0002-per-task-tdd-evidence.md) | **日本語** | [简体中文](0002-per-task-tdd-evidence.zh-CN.md)

- 状態: Accepted
- 日付: 2026-11-27

## 背景 (Context)
TDD は必須であり、あとから監査できなければなりません。複数のエージェントが並列に作業するので、共有のログファイルが 1 つだけだと、編集の衝突が絶えません。

## 決定 (Decision)
各タスクに JSON ファイル `docs/tdd/tasks/<TASK-ID>.json`(形式は `docs/tdd/README.md`)を 1 つ持たせ、`scripts/tdd.py`(`start`, `red`, `green`, `refactor`, `done`)を通してのみ書きます。RED のイベントには、意図した理由で終了コードが 0 以外であることが必要です。GREEN と REFACTOR には終了コード 0 が必要です。
`docs/TDD_LOG.md` は生成物 (`tdd.py log`) であり、手で編集しません。`.claude/hooks/` のフックは、`.tdd/config.json` が `enforce=true` のとき、RED の証跡が無ければ `src/`、`include/`、`tools/` への編集を拒否します。

## 結果 (Consequences)
- 並列の作業でも、証跡は衝突しません。
- `tdd.py check` で、CI が履歴全体を検証できます。
- 証跡はコマンドを実行して生成されるので、文章を編集して偽造することはできません。
