---
name: git-workflow
description: >-
  Commit and push safely in this repo — quality gate, commit granularity, Conventional Commits, push, then watch GitHub Actions with gh run watch / gh run list and fix CI failures through TDD. Use when asked to commit, push, open a PR, or check CI.
when_to_use: >-
  "commit this", "push", "open a PR", "is CI green", "CI failed".
---
# Git workflow

Commit or push only when the user (or the lead acting on the user's request) asked for it.

## 1. Quality gate (before every commit)
```bash
make verify
uv run python scripts/tdd.py check
git status --short
```
If `make verify` fails, fix it (through `/tdd` for behaviour changes). Do not commit red.

## 2. Commit granularity
- One task ID per commit where practical; tests and implementation of a cycle go together.
- Stage explicit paths, never `git add -A` / `git add .` (avoids `.cache/paper/`, raw
  results, local settings): `git add src/paper/primitives/ tests/unit/paper/ docs/tdd/tasks/P3-DSU.json`
- Check nothing forbidden is staged: `git diff --cached --name-only`.
- Message (Conventional Commits):
  ```
  feat(paper): ordered merge events from MST (Sec 3.3) [P3-MERGE-EVENTS]

  <what and why; deviations; evidence file>
  ```
  Types: feat, fix, test, refactor, perf, docs, build, ci, chore. Scope = module.
- After committing: `uv run python scripts/tdd.py done <ID> --commit "$(git rev-parse --short HEAD)"`
  and commit `docs/tdd/tasks/<ID>.json` + `docs/TASKS.md` as `docs(tdd): close <ID>`.

## 3. Push
```bash
git fetch origin && git status -sb          # behind? rebase local-only commits: git pull --rebase
git push origin <branch>
```
Forbidden: `git push --force` / `--force-with-lease`, rewriting pushed commits,
`git remote add|set-url|remove`, pushing with failing local tests, `--no-verify`.

## 4. Watch CI
```bash
gh run list --branch <branch> --limit 5
gh run watch <run-id> --exit-status        # blocks until done; non-zero on failure
gh run view <run-id> --log-failed          # failing step logs
```
Record the result (run id, conclusion) in `docs/PROGRESS.md` → "CI".

## 5. CI failed
1. Read `--log-failed`; reproduce locally with the same preset (`ci` / `ci-release` /
   `sanitizer`): `cmake --preset ci && cmake --build --preset ci && ctest --preset ci`.
2. Platform-only failure (Linux vs macOS): note compiler/version from the log.
3. Fix through `/tdd`: failing test (or the CI reproduction command) as RED, then fix, GREEN.
4. New commit (`fix(ci): ...` or `fix(<scope>): ...`); never amend a pushed commit, never
   disable the job/test, never mark it `continue-on-error`.
5. Push and watch again.
