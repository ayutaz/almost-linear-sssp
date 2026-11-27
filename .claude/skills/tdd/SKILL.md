---
name: tdd
description: >-
  Run one RED→GREEN→REFACTOR cycle for a task in docs/TASKS.md and record the evidence in docs/tdd/tasks/<ID>.json with scripts/tdd.py. Use before writing or changing any code under src/, include/ or tools/, and when fixing a bug.
when_to_use: >-
  "implement a feature", "fix a bug", "add a test", "start task P3-...", or when a hook blocks an edit because no RED evidence exists.
argument-hint: "[TASK-ID]"
allowed-tools: Bash(uv run python scripts/tdd.py *) Bash(cmake --build --preset *) Bash(ctest --preset *) Bash(uv run pytest *)
---
# TDD cycle for task $ARGUMENTS

Evidence format: `docs/tdd/README.md`. One task = one file `docs/tdd/tasks/<ID>.json`.

## 0. Orient
```bash
uv run python scripts/tdd.py status          # is another task active? finish or park it first
```
Read the task row in `docs/TASKS.md` (acceptance criteria, dependencies). Do not start a task
whose dependencies are not DONE.

## 1. Start
```bash
uv run python scripts/tdd.py start <ID> \
  --spec "<observable behaviour, one sentence>" \
  --paper "<Sec x.y, Lemma a.b, Algorithm c line d>" \
  --impl "src/paper/<module>/*.cpp" \
  --test "tests/unit/<module>/test_*.cpp"
```
Use `--paper none` for non-paper work. The `--impl`/`--test` globs define the files this task
may touch; hooks compare edits against them (check `tdd.py start --help` for repeating flags).

## 2. RED — write the failing test first
- Put tests under the category directory: `tests/{unit,integration,differential,property,regression}/`.
  Files are globbed; do not edit `tests/CMakeLists.txt`.
- Assert behaviour from the spec / paper, never values copied from an implementation run.
- Run it and record the failure:
```bash
uv run python scripts/tdd.py red <ID> --reason "<why it fails: missing symbol / wrong value X vs Y>" \
  -- bash -c 'cmake --build --preset debug && ctest --preset debug -R "<test name regex>"'
```
- **Check the reason.** Read the output tail. Acceptable: missing header/symbol for the new
  feature, assertion failure on the new behaviour. Not acceptable: unrelated compile error,
  broken toolchain, typo in the test, wrong test filter (0 tests run). Fix and re-run `red`.
- `red` refuses exit code 0. If the test already passes, the test does not specify anything
  new: strengthen it.

## 3. GREEN — minimal implementation
- Write the least code that makes the test pass *for the right reason* (no special-casing
  inputs, no hard-coded outputs, no fallback to Dijkstra to cover a missing component).
```bash
uv run python scripts/tdd.py green <ID> -- bash -c 'cmake --build --preset debug && ctest --preset debug'
```
Run the whole suite, not just the new test.

## 4. REFACTOR
Clean up with tests green. Record:
```bash
uv run python scripts/tdd.py refactor <ID> -- bash -c 'cmake --build --preset debug && ctest --preset debug'
```
Repeat 2–4 for the next behaviour of the same task (several cycles per task are fine).

## 5. Finish
```bash
make verify
uv run python scripts/tdd.py check
git add <impl files> <test files> docs/tdd/tasks/<ID>.json
git commit -m "feat(<scope>): <summary> [<ID>]"
uv run python scripts/tdd.py done <ID> --commit "$(git rev-parse --short HEAD)"
uv run python scripts/tdd.py log            # regenerates docs/TDD_LOG.md
```
Then update `docs/TASKS.md` (Status, Test Evidence, Commit SHA) and commit the evidence/doc
update as `docs(tdd): close <ID>`. `docs/TASKS.md`, `docs/PROGRESS.md` and `docs/TDD_LOG.md`
are shared: a parallel agent leaves them (and `tdd.py log`) to the lead.

## Test-only tasks (verification tests that add no behaviour)
```bash
uv run python scripts/tdd.py start <ID> --test-only "<why no behaviour changes>" \
  --spec "<property checked>" --test "tests/property/test_*.cpp" --covers "src/paper/<module>/*.cpp"
```
- Exercised files go in `--covers`, never `--impl`; `tdd.py` rejects `--impl` here.
- The commit holds tests and evidence only. A comment or note fix in an implementation file is
  a separate commit (and task).
- A test that is GREEN at once is the expected result here. Prove it has teeth with a
  **mutation check**: break the covered code in one place at a time, run the test, confirm it
  fails, then revert. Never commit a mutation; check `git diff` is clean afterwards. Report the
  score, e.g. "detected 10 of 11 mutations", and name the survivors.

## Recording evidence honestly
- Never steer the exit code of a RED/GREEN command (`exit 1`, `|| true`, `; false`, `!`).
  The recorded command must be the real test command.
- Never hand-edit an event. A wrong event stays; explain it in the next event's `--reason`
  or note.
- `files.impl` lists **every** implementation file the task changes, not just the main one
  (seen: a benchmark-registration task listed only docs while it changed `benchmarks/`).
- A change to a shared header (`include/sssp/*.hpp`) is a task of its own: `start` needs
  `--spec` (the changed contract) and `--test` (a test pinning it), as for any other task.

## Rules
- Never edit a test's expected value to match output. If the expectation was wrong, record why
  in the task's `--reason` of a new RED event and cite the paper/spec that justifies the change.
- Never delete or `SKIP` a failing test to get GREEN.
- Never hand-edit evidence JSON or `docs/TDD_LOG.md`.
- Bug fix = regression test reproducing the bug (RED) before the fix.
