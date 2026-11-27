# Contributing

**English** | [日本語](CONTRIBUTING.ja.md) | [简体中文](CONTRIBUTING.zh-CN.md)

This project reproduces arXiv 2611.34911v1 ("Linear Single-Source Shortest Path, Almost").
Correctness and honest reporting come before speed. `CLAUDE.md` holds the binding rules; this
file is the human-oriented summary. For an overview of the documents, see
[docs/README.md](docs/README.md).

## Setup

```bash
./scripts/bootstrap.sh          # toolchain check + uv sync + configure
make build && make test
./scripts/fetch_paper.sh        # paper into .cache/paper/ (never commit it)
```

Requirements:

- Clang (Apple Clang or LLVM) with C++20, CMake ≥ 3.24, Ninja
- uv, Python ≥ 3.12
- clang-format, clang-tidy
- `gh` (for CI checks)

Python dependencies: `uv add <pkg>` (dev: `uv add --dev <pkg>`), run with `uv run`.
Do not use pip.

## Workflow

1. Pick a task from `docs/TASKS.md` whose dependencies are DONE.
2. TDD with evidence (`.claude/skills/tdd/SKILL.md`):
   `scripts/tdd.py start` → failing test → `red` → implement → `green` → `refactor`.
3. Paper code: keep `docs/PAPER_MAPPING.md` in sync and record deviations
   (`.claude/skills/paper-reproduction/SKILL.md`).
4. Verify: `make verify` (and `make sanitizer` for algorithmic changes).
5. Commit with Conventional Commits, mention the task ID; then `tdd.py done <ID> --commit <SHA>`.
6. Push, then watch CI: `gh run watch <id> --exit-status`.

## Where code goes

| Directory | Content |
|---|---|
| `include/sssp/` | public headers (shared `CsrGraph`, `SsspResult`, `Diagnostics`) |
| `src/baseline/` | Dijkstra variants, Bellman–Ford |
| `src/paper/` | Mode A, paper-faithful |
| `src/experimental/paper_experimental/` | Mode B, tuned thresholds |
| `src/experimental/practical/` | Mode C, practical optimisations |
| `tests/{unit,integration,differential,property,regression}/` | Catch2 tests (globbed) |
| `tests/harness/`, `tests/python/` | pytest |
| `benchmarks/` | Google Benchmark drivers |

## Rules that reviewers enforce

- Tests are written first and fail for the intended reason; expected values come from the
  spec/paper/oracle, never from the implementation's output.
- No hidden fallbacks: every base case or substitution is recorded in `Diagnostics.fallbacks`.
- No unmeasured performance numbers; every figure cites a saved run.
- No fast-math in correctness builds.
- No force-push, no history rewriting on pushed branches, no remote changes.
- Decisions are recorded as ADRs in `docs/adr/`.

## Style

- C++: `.clang-format`, `.clang-tidy`; warnings are errors in CI.
- Python: Ruff (line length 100).
- `make format` fixes formatting; `make lint` checks it.
