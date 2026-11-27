"""PreCompact hook: writes a state snapshot into docs/PROGRESS.md (AUTO block)."""

from __future__ import annotations

import json
from pathlib import Path

from harness_util import GREEN_EVENT, RED_EVENT, write_task

BEGIN = "<!-- AUTO:precompact -->"
END = "<!-- /AUTO -->"


def _payload(trigger: str = "auto") -> dict:
    return {
        "hook_event_name": "PreCompact",
        "trigger": trigger,
        "custom_instructions": "",
        "transcript_path": "/nonexistent/transcript.jsonl",
        "permission_mode": "default",
    }


def _progress(root: Path) -> str:
    return (root / "docs" / "PROGRESS.md").read_text(encoding="utf-8")


def _block(text: str) -> str:
    start = text.index(BEGIN)
    end = text.index(END, start)
    return text[start : end + len(END)]


def _write_last_check(root: Path, categories: dict) -> None:
    rt = root / ".tdd" / "runtime"
    rt.mkdir(parents=True, exist_ok=True)
    (rt / "last_check.json").write_text(
        json.dumps({"updated": "2026-01-01T00:00:00+09:00", "categories": categories})
    )


def test_creates_progress_file_when_missing(repo, run_hook):
    res = run_hook("pre_compact", _payload("manual"), repo)
    assert res.code == 0, res.stderr
    text = _progress(repo)
    assert text.startswith("# Progress\n\n")
    assert text.count(BEGIN) == 1
    assert text.count(END) == 1
    assert "manual" in _block(text)


def test_appends_block_and_keeps_human_text(repo, run_hook):
    (repo / "docs" / "PROGRESS.md").write_text("# Progress\n\nhand-written notes\n")
    res = run_hook("pre_compact", _payload(), repo)
    assert res.code == 0, res.stderr
    text = _progress(repo)
    assert text.startswith("# Progress\n\nhand-written notes\n")
    assert text.count(BEGIN) == 1


def test_replaces_only_the_block_and_is_idempotent(repo, run_hook):
    before = "# Progress\n\nintro line\n\n"
    after = "\n\n## Human section\n\nkeep me exactly\n"
    (repo / "docs" / "PROGRESS.md").write_text(before + BEGIN + "\nOLD CONTENT\n" + END + after)
    for _ in range(2):
        res = run_hook("pre_compact", _payload(), repo)
        assert res.code == 0, res.stderr
    text = _progress(repo)
    assert "OLD CONTENT" not in text
    assert text.count(BEGIN) == 1
    assert text.count(END) == 1
    assert text.startswith(before + BEGIN)
    assert text.endswith(END + after)


def test_writes_snapshot_with_same_content(repo, run_hook):
    res = run_hook("pre_compact", _payload(), repo)
    assert res.code == 0, res.stderr
    snap = repo / ".tdd" / "runtime" / "precompact_snapshot.md"
    assert snap.is_file()
    body = snap.read_text(encoding="utf-8")
    assert body.strip()
    assert body.strip() in _progress(repo)


def test_active_tasks_and_next_commands(repo, run_hook):
    write_task(repo, {"id": "T-NONE", "files": {"impl": ["src/core/x*.cpp"], "tests": []}})
    write_task(repo, {"id": "T-RED", "events": [RED_EVENT]})
    write_task(repo, {"id": "T-GREEN", "phase": "GREEN", "events": [RED_EVENT, GREEN_EVENT]})
    write_task(repo, {"id": "T-DONE", "phase": "DONE", "events": [RED_EVENT]}, active=False)
    res = run_hook("pre_compact", _payload(), repo)
    assert res.code == 0, res.stderr
    block = _block(_progress(repo))
    assert "src/core/x*.cpp" in block
    assert "tdd.py red T-NONE" in block
    assert "tdd.py green T-RED -- ctest --preset debug" in block
    assert "tdd.py refactor T-GREEN" in block
    assert "tdd.py done T-GREEN" in block
    assert "T-DONE" in block
    assert "make verify" in block
    # stdout mentions the active tasks for the summary
    assert "T-RED" in res.stdout
    assert "docs/PROGRESS.md" in res.stdout


def test_failures_and_unresolved_are_transcribed(repo, run_hook):
    _write_last_check(
        repo,
        {
            "cpp-test": {
                "status": "fail",
                "timestamp": "2026-01-01T00:00:00+09:00",
                "file": "src/core/a.cpp",
                "command": "ctest --preset debug -R widget",
                "exit_code": 8,
                "output_tail": "line1\nline2\nFAILED: widget_boundary",
            },
            "lint": {"status": "pass", "command": "ruff check", "exit_code": 0},
        },
    )
    (repo / ".tdd" / "runtime").mkdir(parents=True, exist_ok=True)
    (repo / ".tdd" / "runtime" / "stop_gate_unresolved.md").write_text("# x\n")
    (repo / "docs" / "TASKS.md").write_text("- [ ] P9-FOO BLOCKED: waiting on bar\n- [ ] ok\n")
    res = run_hook("pre_compact", _payload(), repo)
    assert res.code == 0, res.stderr
    block = _block(_progress(repo))
    assert "cpp-test" in block
    assert "ctest --preset debug -R widget" in block
    assert "FAILED: widget_boundary" in block
    assert ".tdd/runtime/stop_gate_unresolved.md" in block
    assert "P9-FOO BLOCKED: waiting on bar" in block
    assert "ruff check" not in block  # passing categories are not problems
    assert "cpp-test" in res.stdout


def test_paper_mapping_excerpt(repo, run_hook):
    lines = ["# Paper mapping", ""] + [f"row {i}" for i in range(1, 20)]
    (repo / "docs" / "PAPER_MAPPING.md").write_text("\n".join(lines) + "\n")
    res = run_hook("pre_compact", _payload(), repo)
    assert res.code == 0, res.stderr
    block = _block(_progress(repo))
    assert "row 9" in block
    assert "row 10" not in block  # only the first 10 non-empty lines


def test_paper_mapping_absent(repo, run_hook):
    res = run_hook("pre_compact", _payload(), repo)
    assert res.code == 0, res.stderr
    block = _block(_progress(repo))
    assert "PAPER_MAPPING" in block
    assert "none" in block


def test_last_commit_is_shown(repo, run_hook):
    res = run_hook("pre_compact", _payload(), repo)
    assert res.code == 0, res.stderr
    assert "init" in _block(_progress(repo))


def test_stdout_is_short_plain_text(repo, run_hook):
    for i in range(40):
        write_task(repo, {"id": f"T-{i:03d}-" + "x" * 60, "events": [RED_EVENT]})
    res = run_hook("pre_compact", _payload(), repo)
    assert res.code == 0, res.stderr
    out = res.stdout.strip()
    assert out
    assert len(res.stdout) <= 2000
    try:
        json.loads(out)
    except ValueError:
        pass
    else:  # pragma: no cover - must not be JSON
        raise AssertionError("stdout must be plain text, not JSON")


def test_exit_zero_with_broken_files(repo, run_hook):
    rt = repo / ".tdd" / "runtime"
    rt.mkdir(parents=True, exist_ok=True)
    (rt / "active.json").write_text("{not json")
    (rt / "last_check.json").write_text("[1, 2")
    (repo / "docs/tdd/tasks/BROKEN.json").write_text("{{{")
    (repo / ".tdd" / "config.json").write_text("nope")
    res = run_hook("pre_compact", _payload(), repo)
    assert res.code == 0, res.stderr
    assert BEGIN in _progress(repo)


def test_exit_zero_on_invalid_stdin(repo, run_hook):
    for raw in ("", "not json", "[1,2,3]", '{"trigger": 5}'):
        res = run_hook("pre_compact", {}, repo, raw_stdin=raw)
        assert res.code == 0, (raw, res.stderr)


def test_tasks_legend_line_is_not_a_blocked_task(repo, run_hook):
    (repo / "docs" / "TASKS.md").write_text(
        "- **Status**: TODO / IN_PROGRESS / BLOCKED / DONE.\n\n"
        "| ID | Task | Status |\n| :- | :- | :- |\n"
        "| P9-BAR | waits on upstream | BLOCKED |\n| P9-OK | fine | TODO |\n"
    )
    res = run_hook("pre_compact", _payload(), repo)
    assert res.code == 0, res.stderr
    block = _block(_progress(repo))
    assert "P9-BAR" in block
    assert "TODO / IN_PROGRESS / BLOCKED" not in block
    assert "P9-OK" not in block


def test_paper_mapping_falls_back_to_translation(repo, run_hook):
    (repo / "docs" / "PAPER_MAPPING.ja.md").write_text("# 論文対応表\n\nja row 1\n")
    res = run_hook("pre_compact", _payload(), repo)
    assert res.code == 0, res.stderr
    block = _block(_progress(repo))
    assert "ja row 1" in block


def test_blocked_is_decided_by_the_status_column_only(repo, run_hook):
    """Regression: a DONE task named P3-SFM-BLOCKED was reported as blocked."""
    (repo / "docs" / "TASKS.md").write_text(
        "| ID | Priority | Status | Dependencies |\n|---|---|---|---|\n"
        "| P3-SFM-BLOCKED | P0 | DONE | part of P3-INTERVAL-MIN |\n"
        "| P3-INTERVAL-MIN | P0 | IN_PROGRESS | see docs/tdd/tasks/P3-SFM-BLOCKED.json |\n"
        "| P9-WAIT | P1 | BLOCKED | upstream fix |\n"
        "| P9-SPACES | P1 |  BLOCKED  | padded cell |\n"
        "| P9-PARTIAL | P1 | BLOCKED-ish | not an exact match |\n"
    )
    res = run_hook("pre_compact", _payload(), repo)
    assert res.code == 0, res.stderr
    block = _block(_progress(repo))
    assert "P9-WAIT" in block
    assert "P9-SPACES" in block
    assert "P3-SFM-BLOCKED" not in block
    assert "P3-INTERVAL-MIN" not in block
    assert "P9-PARTIAL" not in block


def test_snapshot_is_written_in_english(repo, run_hook):
    """The AUTO block lands in a public English document: no Japanese boilerplate."""
    import re

    _write_last_check(
        repo,
        {
            "unit": {
                "status": "fail",
                "command": "ctest -L unit",
                "exit_code": 8,
                "output_tail": "x",
            }
        },
    )
    write_task(
        repo, {"id": "T-RED", "files": {"impl": ["src/a.cpp"], "tests": []}, "events": [RED_EVENT]}
    )
    write_task(repo, {"id": "T-NEW"})
    write_task(
        repo, {"id": "T-OLD", "phase": "DONE", "events": [RED_EVENT, GREEN_EVENT]}, active=False
    )
    (repo / ".tdd" / "runtime" / "stop_gate_unresolved.md").write_text("# x\n")
    (repo / "docs" / "TASKS.md").write_text("| ID | Status |\n|---|---|\n| P9-B | BLOCKED |\n")
    res = run_hook("pre_compact", _payload(), repo)
    assert res.code == 0, res.stderr
    block = _block(_progress(repo))
    assert not re.search(r"[\u3040-\u30ff\u3400-\u9fff\uff00-\uffef]", block), block
    for heading in (
        "## Automatic snapshot (PreCompact)",
        "### Completed tasks",
        "### Active tasks",
        "### Open problems",
        "### TDD state and next commands",
        "### Last commit",
        "### PAPER_MAPPING.md excerpt",
    ):
        assert heading in block
    assert "- Generated: " in block
    assert "Implement, then `uv run python scripts/tdd.py green T-RED" in block
    assert "Finally run `make verify`" in block
