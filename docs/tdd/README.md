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
