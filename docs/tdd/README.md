# TDD evidence format

Every feature or bug fix has one evidence file `docs/tdd/tasks/<TASK-ID>.json`
(one file per task so that parallel work never conflicts). `docs/TDD_LOG.md` is
generated from these files by `uv run python scripts/tdd.py log`.

```json
{
  "id": "P2-DIJKSTRA-BINARY",
  "spec": "what behaviour is being specified",
  "paper_ref": "Sec 3.2, Lemma 3.5 (or null for non-paper work)",
  "files": {"impl": ["src/baseline/dijkstra_binary.cpp"], "tests": ["tests/unit/baseline/test_dijkstra.cpp"]},
  "phase": "RED | GREEN | REFACTOR | DONE",
  "events": [
    {"phase": "RED", "tests": ["test case name"], "command": "ctest --preset debug -R ...",
     "exit_code": 8, "reason": "why it failed (must be the intended reason)",
     "timestamp": "2026-11-27T10:30:00+09:00", "output_tail": "last ~20 lines"},
    {"phase": "GREEN", "tests": ["..."], "command": "...", "exit_code": 0, "timestamp": "...", "output_tail": "..."},
    {"phase": "REFACTOR", "command": "ctest --preset debug", "exit_code": 0, "timestamp": "...", "output_tail": "..."}
  ],
  "commit": null
}
```

Rules:
- A RED event must have a non-zero exit code caused by the missing/incorrect behaviour,
  not by an unrelated broken build environment.
- A GREEN/REFACTOR event must have exit code 0.
- Several RED→GREEN cycles may be appended to the same task.
- Never edit a test's expected value to match the implementation's output.

## Test-only tasks

A task that only adds verification tests for behaviour that already exists (a differential
test of finished baselines, an integration test composing components that were each developed
test-first) cannot honestly produce a RED event. Never manufacture one (for example by
breaking the implementation or the test on purpose). Record the task as test-only instead:

```json
{
  "id": "P2-DIFF-BASELINES",
  "kind": "test-only",
  "test_only_reason": "why no RED is possible: which existing behaviour is verified and why it is not new",
  "spec": "...",
  "paper_ref": null,
  "files": {"impl": [], "tests": ["tests/differential/test_baselines_differential.cpp"],
            "covers": ["src/baseline/dijkstra_binary.cpp"]},
  "phase": "GREEN | REFACTOR | DONE",
  "events": [
    {"phase": "GREEN", "tests": ["..."], "command": "...", "exit_code": 0, "timestamp": "...", "output_tail": "..."}
  ],
  "commit": null
}
```

Create one with
`uv run python scripts/tdd.py start <ID> --spec "..." --test <glob> [--covers <glob>] --test-only "<reason>"`,
then record `tdd.py green <ID> -- <test cmd>` (`tdd.py red` is refused for test-only tasks: a
failing verification test means a bug, which gets a normal task with a RED event).

Rules checked by `tdd.py check`:
- `kind` is optional; the only allowed value is `"test-only"`.
- Only for `"kind": "test-only"` may the first event be GREEN. Without `kind`, a task whose
  first event is not RED stays invalid.
- `test_only_reason` must be a non-empty string.
- `files.impl` must be empty. Files the tests exercise but the task does not change go into
  `files.covers` (optional, a list of strings; allowed for every task). The only exception is
  a task whose `commit` is set: then `files.impl` may be non-empty as long as that commit
  changes no file under the implementation directories (`impl_dirs` in `.tdd/config.json`:
  `src/`, `include/`, `tools/`).
