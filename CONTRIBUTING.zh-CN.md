# 贡献指南

[English](CONTRIBUTING.md) | [日本語](CONTRIBUTING.ja.md) | **简体中文**

本项目复现 arXiv 2611.34911v1("Linear Single-Source Shortest Path, Almost")。
正确性和如实报告优先于速度。具有约束力的规则写在 `CLAUDE.md` 中;本文件是面向人的摘要。
各文档的概览见 [docs/README.zh-CN.md](docs/README.zh-CN.md)。

## 环境准备

```bash
./scripts/bootstrap.sh          # 检查工具链 + uv sync + 配置
make build && make test
./scripts/fetch_paper.sh        # 把论文下载到 .cache/paper/(切勿提交)
```

环境要求:

- 支持 C++20 的 Clang(Apple Clang 或 LLVM)、CMake ≥ 3.24、Ninja
- uv、Python ≥ 3.12
- clang-format、clang-tidy
- `gh`(用于查看 CI)

Python 依赖:用 `uv add <pkg>` 添加(开发依赖:`uv add --dev <pkg>`),用 `uv run` 运行。
不要使用 pip。

## 工作流程

1. 从 `docs/TASKS.md` 中选择一个依赖项都已 DONE 的任务。
2. 按 TDD 进行并留下证据(`.claude/skills/tdd/SKILL.md`):
   `scripts/tdd.py start` → 编写会失败的测试 → `red` → 实现 → `green` → `refactor`。
3. 论文相关代码:保持 `docs/PAPER_MAPPING.md` 同步,并记录所有偏差
   (`.claude/skills/paper-reproduction/SKILL.md`)。
4. 验证:运行 `make verify`(算法改动还需运行 `make sanitizer`)。
5. 按 Conventional Commits 提交,并写明任务 ID;然后运行 `tdd.py done <ID> --commit <SHA>`。
6. 推送后观察 CI:`gh run watch <id> --exit-status`。

## 代码放在哪里

| 目录 | 内容 |
|---|---|
| `include/sssp/` | 公共头文件(共用的 `CsrGraph`、`SsspResult`、`Diagnostics`) |
| `src/baseline/` | Dijkstra 各变体、Bellman–Ford |
| `src/paper/` | 模式 A,忠实于论文 |
| `src/experimental/paper_experimental/` | 模式 B,调整阈值 |
| `src/experimental/practical/` | 模式 C,实用优化 |
| `tests/{unit,integration,differential,property,regression}/` | Catch2 测试(自动 glob 收集) |
| `tests/harness/`、`tests/python/` | pytest |
| `benchmarks/` | Google Benchmark 驱动程序 |

## 评审者执行的规则

- 测试先写,并且要因预期的原因而失败;期望值来自规格、论文或 oracle,绝不来自实现的输出。
- 不允许隐藏的回退:每一次进入基础情形或替换算法都要记录到 `Diagnostics.fallbacks` 中。
- 不允许出现未经测量的性能数据;每个数字都要注明对应的已保存运行记录。
- 用于正确性检验的构建中不允许使用 fast-math。
- 不允许 force-push,不允许改写已推送分支的历史,不允许修改远程仓库设置。
- 决策以 ADR 的形式记录在 `docs/adr/` 中。

## 代码风格

- C++:`.clang-format`、`.clang-tidy`;CI 中警告视为错误。
- Python:Ruff(行宽 100)。
- `make format` 修正格式;`make lint` 检查格式。
