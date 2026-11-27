# 0002. 按任务划分的 TDD 证据文件

[English](0002-per-task-tdd-evidence.md) | [日本語](0002-per-task-tdd-evidence.ja.md) | **简体中文**

- 状态: Accepted
- 日期: 2026-11-27

## 背景 (Context)
TDD 是强制的,而且必须可审计。多个代理并行工作,只用一个共享的日志文件会不断产生编辑冲突。

## 决定 (Decision)
每个任务有一个 JSON 文件 `docs/tdd/tasks/<TASK-ID>.json`(格式见 `docs/tdd/README.md`),只通过 `scripts/tdd.py`(`start`, `red`, `green`, `refactor`, `done`)写入。RED 事件要求退出码非零,且原因符合预期;GREEN/REFACTOR 要求退出码为 0。
`docs/TDD_LOG.md` 是生成的 (`tdd.py log`),不手工编辑。当 `.tdd/config.json` 设置了 `enforce=true` 时,`.claude/hooks/` 中的钩子会在没有 RED 证据的情况下拒绝编辑 `src/`、`include/`、`tools/`。

## 后果 (Consequences)
- 并行工作不会在证据上产生冲突。
- `tdd.py check` 可以在 CI 中验证整个历史。
- 证据由运行命令生成,因此无法通过编辑文本来伪造。
