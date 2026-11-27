"""Stop hook: keeps Claude working while tests fail / RED is left / work is uncommitted."""

from __future__ import annotations

import json
import os
import subprocess
import time
from pathlib import Path

import pytest
from harness_util import GREEN_EVENT, RED_EVENT, write_task

REPORT = Path(".tdd") / "runtime" / "stop_gate_unresolved.md"


def _git(root: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True)


def _clean(root: Path) -> None:
    """Make every check pass: PROGRESS.md committed, working tree clean."""
    (root / "docs" / "PROGRESS.md").write_text("# Progress\n")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "progress")


def _config(root: Path, **stop_gate) -> None:
    """Write .tdd/config.json and commit it so it is not an uncommitted change."""
    cfg = {"enforce": False, "stop_gate": stop_gate}
    (root / ".tdd" / "config.json").write_text(json.dumps(cfg))
    _git(root, "add", ".tdd/config.json")
    if subprocess.run(["git", "-C", str(root), "diff", "--cached", "--quiet"]).returncode:
        _git(root, "commit", "-q", "-m", "cfg")


def _payload(**kw) -> dict:
    data = {
        "hook_event_name": "Stop",
        "stop_hook_active": False,
        "transcript_path": "/nonexistent/transcript.jsonl",
        "permission_mode": "default",
        "last_assistant_message": "done",
    }
    data.update(kw)
    return data


def _state(root: Path) -> dict:
    return json.loads((root / ".tdd" / "runtime" / "stop_gate.json").read_text())


def _fail_check(root: Path, cat: str = "cpp-test", cmd: str = "ctest --preset debug") -> None:
    rt = root / ".tdd" / "runtime"
    rt.mkdir(parents=True, exist_ok=True)
    (rt / "last_check.json").write_text(
        json.dumps(
            {
                "updated": "2026-01-01T00:00:00+09:00",
                "categories": {
                    cat: {
                        "status": "fail",
                        "timestamp": "2026-01-01T00:00:00+09:00",
                        "file": "src/core/a.cpp",
                        "command": cmd,
                        "exit_code": 8,
                        "output_tail": "FAILED",
                    },
                    "lint": {"status": "pass", "command": "ruff", "exit_code": 0},
                },
            }
        )
    )


def _blocked(res) -> bool:
    data = res.json or {}
    return data.get("decision") == "block"


def _transcript(root: Path, rows: list) -> str:
    path = root.parent / "transcript.jsonl"
    with open(path, "w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row, ensure_ascii=False) + "\n")
    return str(path)


def _user(content) -> dict:
    return {"type": "user", "message": {"role": "user", "content": content}}


def _assistant(text: str) -> dict:
    return {
        "type": "assistant",
        "message": {"role": "assistant", "content": [{"type": "text", "text": text}]},
    }


# ------------------------------------------------------------------ single issues
def test_clean_repo_allows_and_resets_counter(repo, run_hook):
    _fail_check(repo)
    _clean(repo)
    assert _blocked(run_hook("stop_gate", _payload(), repo))
    assert _state(repo)["consecutive"] == 1
    (repo / ".tdd" / "runtime" / "last_check.json").unlink()
    res = run_hook("stop_gate", _payload(), repo)
    assert res.code == 0, res.stderr
    assert res.stdout.strip() == ""
    assert _state(repo)["consecutive"] == 0


def test_failed_tests_block(repo, run_hook):
    _clean(repo)
    _fail_check(repo, "cpp-test", "ctest --preset debug -R widget")
    res = run_hook("stop_gate", _payload(), repo)
    assert res.code == 0, res.stderr
    assert _blocked(res)
    reason = res.json["reason"]
    assert "cpp-test" in reason
    assert "ctest --preset debug -R widget" in reason
    assert "lint" not in reason


def test_timed_out_check_does_not_block(repo, run_hook):
    """A light check that ran out of time is inconclusive, not a failure."""
    _clean(repo)
    _fail_check(repo)
    path = repo / ".tdd" / "runtime" / "last_check.json"
    data = json.loads(path.read_text())
    data["categories"]["cpp-test"].update(status="timeout", exit_code=124)
    path.write_text(json.dumps(data))
    res = run_hook("stop_gate", _payload(), repo)
    assert res.code == 0, res.stderr
    assert not _blocked(res)


def test_red_task_blocks(repo, run_hook):
    write_task(repo, {"id": "T-RED", "events": [RED_EVENT]})
    write_task(repo, {"id": "T-OK", "phase": "GREEN", "events": [RED_EVENT, GREEN_EVENT]})
    _clean(repo)
    res = run_hook("stop_gate", _payload(), repo)
    assert _blocked(res)
    reason = res.json["reason"]
    assert "T-RED" in reason
    assert "T-OK" not in reason
    assert "tdd.py green" in reason


def test_inactive_red_task_does_not_block(repo, run_hook):
    write_task(repo, {"id": "T-OLD", "events": [RED_EVENT]}, active=False)
    _clean(repo)
    res = run_hook("stop_gate", _payload(), repo)
    assert not _blocked(res)


def test_uncommitted_changes_block(repo, run_hook):
    _clean(repo)
    (repo / "src/core/b.cpp").write_text("int b();\n")
    # keep PROGRESS.md newer so only the uncommitted check fires
    future = time.time() + 100
    os.utime(repo / "docs" / "PROGRESS.md", (future, future))
    _config(repo, check_progress_updated=False)
    res = run_hook("stop_gate", _payload(), repo)
    assert _blocked(res)
    reason = res.json["reason"]
    assert "src/core/b.cpp" in reason
    assert "コミット" in reason


def test_progress_missing_blocks(repo, run_hook):
    _config(repo, check_uncommitted=False)
    res = run_hook("stop_gate", _payload(), repo)
    assert _blocked(res)
    assert "PROGRESS.md" in res.json["reason"]


def test_progress_older_than_changes_blocks(repo, run_hook):
    _clean(repo)
    _config(repo, check_uncommitted=False)
    past = time.time() - 1000
    os.utime(repo / "docs" / "PROGRESS.md", (past, past))
    (repo / "src/core/a.cpp").write_text("int a() { return 2; }\n")
    res = run_hook("stop_gate", _payload(), repo)
    assert _blocked(res)
    assert "PROGRESS.md" in res.json["reason"]
    # updating PROGRESS.md (newer mtime) resolves it
    future = time.time() + 1000
    os.utime(repo / "docs" / "PROGRESS.md", (future, future))
    res = run_hook("stop_gate", _payload(), repo)
    assert not _blocked(res), res.stdout


def test_deleted_files_do_not_require_progress_update(repo, run_hook):
    _clean(repo)
    _config(repo, check_uncommitted=False)
    past = time.time() - 1000
    os.utime(repo / "docs" / "PROGRESS.md", (past, past))
    (repo / "src/core/a.cpp").unlink()
    res = run_hook("stop_gate", _payload(), repo)
    assert not _blocked(res), res.stdout


# ----------------------------------------------------------------- loop guards
def test_stop_hook_active_allows_and_reports(repo, run_hook):
    _clean(repo)
    _fail_check(repo)
    res = run_hook("stop_gate", _payload(stop_hook_active=True), repo)
    assert res.code == 0, res.stderr
    assert not _blocked(res)
    assert "systemMessage" in (res.json or {})
    report = (repo / REPORT).read_text(encoding="utf-8")
    assert "test-session" in report
    assert "stop_hook_active" in report
    assert "cpp-test" in report
    assert "python3 .claude/hooks/stop_gate.py" in report
    assert "git status --porcelain" in report
    assert "cat .tdd/runtime/last_check.json" in report
    assert "uv run python scripts/tdd.py status" in report


def test_same_reason_third_time_allows_and_reports(repo, run_hook):
    _clean(repo)
    _fail_check(repo)
    _config(repo, check_uncommitted=False)
    r1 = run_hook("stop_gate", _payload(), repo)
    r2 = run_hook("stop_gate", _payload(), repo)
    assert _blocked(r1) and _blocked(r2)
    assert not (repo / REPORT).exists()
    assert _state(repo)["consecutive"] == 2
    r3 = run_hook("stop_gate", _payload(), repo)
    assert r3.code == 0
    assert not _blocked(r3)
    report = (repo / REPORT).read_text(encoding="utf-8")
    assert "3" in report
    assert "cpp-test" in report


def test_report_is_not_duplicated(repo, run_hook):
    _clean(repo)
    _fail_check(repo)
    for _ in range(3):
        run_hook("stop_gate", _payload(stop_hook_active=True), repo)
    report = (repo / REPORT).read_text(encoding="utf-8")
    assert report.count("stop_hook_active") >= 1
    sections = [ln for ln in report.splitlines() if ln.startswith("## ")]
    assert len(sections) == 1


def test_changed_reason_restarts_counter(repo, run_hook):
    _clean(repo)
    _config(repo, check_uncommitted=False)
    _fail_check(repo, "cpp-test")
    run_hook("stop_gate", _payload(), repo)
    run_hook("stop_gate", _payload(), repo)
    assert _state(repo)["consecutive"] == 2
    _fail_check(repo, "py-test", "uv run pytest")
    res = run_hook("stop_gate", _payload(), repo)
    assert _blocked(res)
    assert _state(repo)["consecutive"] == 1
    res = run_hook("stop_gate", _payload(), repo)
    assert _blocked(res)
    assert _state(repo)["consecutive"] == 2


# ------------------------------------------------------ explicit user stop request
def _blocked_with_transcript(repo, run_hook, rows) -> bool:
    _clean(repo)
    _fail_check(repo)
    path = _transcript(repo, rows)
    res = run_hook("stop_gate", _payload(transcript_path=path), repo)
    assert res.code == 0, res.stderr
    return _blocked(res)


def test_user_says_stop_in_japanese(repo, run_hook):
    rows = [_user("実装して"), _assistant("ok"), _user("止めて"), _assistant("了解")]
    assert not _blocked_with_transcript(repo, run_hook, rows)


def test_user_says_stop_in_english_list_content(repo, run_hook):
    rows = [_user([{"type": "text", "text": "Please STOP now."}]), _assistant("ok")]
    assert not _blocked_with_transcript(repo, run_hook, rows)


def test_stop_inside_long_prompt_is_ignored(repo, run_hook):
    long_text = "Implement the stop condition of the bucket loop. " + "detail " * 60
    assert len(long_text) > 200
    rows = [_user(long_text), _assistant("ok")]
    assert _blocked_with_transcript(repo, run_hook, rows)


def test_tool_result_and_system_rows_are_ignored(repo, run_hook):
    rows = [
        _user("続けて"),
        _assistant("running"),
        _user([{"type": "tool_result", "tool_use_id": "x", "content": "stop"}]),
        _user("<system-reminder>stop</system-reminder>"),
        _user("<task-notification>abort</task-notification>"),
        _assistant("done"),
    ]
    assert _blocked_with_transcript(repo, run_hook, rows)


def test_older_stop_request_superseded_by_newer_prompt(repo, run_hook):
    rows = [_user("止めて"), _assistant("ok"), _user("やっぱり続けて"), _assistant("ok")]
    assert _blocked_with_transcript(repo, run_hook, rows)


def test_broken_transcript_does_not_crash(repo, run_hook):
    _clean(repo)
    _fail_check(repo)
    path = repo.parent / "transcript.jsonl"
    path.write_text('{"type":"user","message":{"content":"止めて"}}\n{broken\n')
    res = run_hook("stop_gate", _payload(transcript_path=str(path)), repo)
    assert res.code == 0, res.stderr
    assert not _blocked(res)  # the last valid human prompt still counts


# ------------------------------------------------------------------- switches
def test_disabled_in_config(repo, run_hook):
    _fail_check(repo)
    _config(repo, enabled=False)
    res = run_hook("stop_gate", _payload(), repo)
    assert res.code == 0
    assert res.stdout.strip() == ""


def test_disabled_by_env(repo, run_hook):
    _fail_check(repo)
    res = run_hook("stop_gate", _payload(), repo, env={"TDD_STOP_GATE": "off"})
    assert res.code == 0
    assert res.stdout.strip() == ""


def test_individual_checks_can_be_disabled(repo, run_hook):
    _clean(repo)
    _fail_check(repo)
    _config(repo, check_failed_tests=False)
    res = run_hook("stop_gate", _payload(), repo)
    assert not _blocked(res), res.stdout


# ------------------------------------------------------------------ robustness
def test_non_git_directory(tmp_path, run_hook):
    root = tmp_path / "plain"
    (root / "docs").mkdir(parents=True)
    (root / "docs" / "PROGRESS.md").write_text("# Progress\n")
    res = run_hook("stop_gate", _payload(), root)
    assert res.code == 0, res.stderr
    assert not _blocked(res), res.stdout


def test_invalid_stdin(repo, run_hook):
    for raw in ("", "not json", "[1]", '{"stop_hook_active": "yes"}'):
        res = run_hook("stop_gate", {}, repo, raw_stdin=raw)
        assert res.code == 0, (raw, res.stderr)
        if res.stdout.strip():
            json.loads(res.stdout)


def test_broken_runtime_files(repo, run_hook):
    _clean(repo)
    rt = repo / ".tdd" / "runtime"
    rt.mkdir(parents=True, exist_ok=True)
    (rt / "last_check.json").write_text("{oops")
    (rt / "stop_gate.json").write_text("[]")
    (rt / "active.json").write_text("null")
    res = run_hook("stop_gate", _payload(), repo)
    assert res.code == 0, res.stderr
    assert not _blocked(res), res.stdout


def test_own_report_does_not_cause_a_new_block(repo, run_hook):
    _clean(repo)
    _fail_check(repo)
    run_hook("stop_gate", _payload(stop_hook_active=True), repo)
    assert (repo / REPORT).exists()
    (repo / ".tdd" / "runtime" / "last_check.json").unlink()
    res = run_hook("stop_gate", _payload(), repo)
    assert not _blocked(res), res.stdout


def test_runtime_state_is_ignored_even_without_gitignore(repo, run_hook):
    """The hook's own .tdd/runtime/ files must not change the reason hash."""
    (repo / ".gitignore").write_text("/build/\n")
    _clean(repo)
    _fail_check(repo)
    r1 = run_hook("stop_gate", _payload(), repo)
    r2 = run_hook("stop_gate", _payload(), repo)
    assert _blocked(r1) and _blocked(r2)
    assert ".tdd/runtime" not in r2.json["reason"]
    assert _state(repo)["consecutive"] == 2
    r3 = run_hook("stop_gate", _payload(), repo)
    assert not _blocked(r3)


def test_uncommitted_reason_is_stable_while_files_change(repo, run_hook):
    """'Same reason' means the same kind of problem, not the same file list."""
    _clean(repo)
    _config(repo, check_progress_updated=False)
    results = []
    for i in range(3):
        (repo / f"wip{i}.txt").write_text("x")  # the dirty set grows every turn
        results.append(_blocked(run_hook("stop_gate", _payload(), repo)))
    assert results == [True, True, False]
    assert (repo / REPORT).exists()


def test_allowed_continuation_stops_are_not_counted_as_blocks(repo, run_hook):
    _clean(repo)
    _config(repo, check_progress_updated=False)
    (repo / "wip.txt").write_text("x")
    assert _blocked(run_hook("stop_gate", _payload(), repo))  # turn 1: block
    assert not _blocked(run_hook("stop_gate", _payload(stop_hook_active=True), repo))
    assert _state(repo)["consecutive"] == 1
    assert _blocked(run_hook("stop_gate", _payload(), repo))  # turn 2 still blocks


def test_report_does_not_leak_absolute_paths(repo, run_hook):
    _clean(repo)
    (repo / "wip.txt").write_text("x")
    run_hook("stop_gate", _payload(stop_hook_active=True), repo)
    report = (repo / REPORT).read_text()
    assert str(repo) not in report
    assert str(Path.home()) not in report


def test_report_lives_outside_tracked_paths(repo, run_hook):
    """The report holds session ids: it must stay in the git-ignored runtime dir."""
    _clean(repo)
    _fail_check(repo)
    run_hook("stop_gate", _payload(stop_hook_active=True, session_id="sess-xyz"), repo)
    assert (repo / REPORT).exists()
    assert not (repo / "docs" / "reports").exists()
    status = subprocess.run(
        ["git", "-C", str(repo), "status", "--porcelain", "--untracked-files=all"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    assert "stop_gate_unresolved" not in status
    ignored = subprocess.run(
        ["git", "-C", str(repo), "check-ignore", "-q", str(REPORT)], capture_output=True
    )
    assert ignored.returncode == 0
    for path in subprocess.run(
        ["git", "-C", str(repo), "ls-files", "--others", "--exclude-standard", "--cached"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.splitlines():
        assert "sess-xyz" not in (repo / path).read_text(errors="replace")


# ------------------------------------------- work of running helpers (background_tasks)
def _bg(kind: str, status: str = "running") -> list[dict]:
    return [{"id": "t1", "type": kind, "status": status, "description": "helper"}]


def _system_message(res) -> str:
    return (res.json or {}).get("systemMessage") or ""


@pytest.mark.parametrize("kind", ["subagent", "teammate", "workflow"])
def test_running_helper_downgrades_uncommitted_to_a_warning(repo, run_hook, kind):
    _clean(repo)
    (repo / "wip.txt").write_text("x")  # uncommitted + PROGRESS.md older than it
    res = run_hook("stop_gate", _payload(background_tasks=_bg(kind)), repo)
    assert res.code == 0, res.stderr
    assert not _blocked(res), res.stdout
    msg = _system_message(res)
    assert "未コミット" in msg
    assert "PROGRESS" in msg
    assert not (repo / REPORT).exists()  # a warning is not an unresolved block


def test_shell_background_task_does_not_downgrade(repo, run_hook):
    _clean(repo)
    (repo / "wip.txt").write_text("x")
    assert _blocked(run_hook("stop_gate", _payload(background_tasks=_bg("shell")), repo))


@pytest.mark.parametrize("status", ["completed", "failed", "killed"])
def test_finished_helper_does_not_downgrade(repo, run_hook, status):
    _clean(repo)
    (repo / "wip.txt").write_text("x")
    payload = _payload(background_tasks=_bg("subagent", status))
    assert _blocked(run_hook("stop_gate", payload, repo))


def test_failed_tests_still_block_while_helper_runs(repo, run_hook):
    _clean(repo)
    _fail_check(repo)
    (repo / "wip.txt").write_text("x")
    res = run_hook("stop_gate", _payload(background_tasks=_bg("subagent")), repo)
    assert _blocked(res)
    reason = res.json["reason"]
    assert "cpp-test" in reason
    assert "未コミット" not in reason  # only the blocking problems are listed


def test_red_task_still_blocks_while_helper_runs(repo, run_hook):
    _clean(repo)
    write_task(repo, {"id": "T-RED", "events": [RED_EVENT]})
    res = run_hook("stop_gate", _payload(background_tasks=_bg("teammate")), repo)
    assert _blocked(res)
    assert "T-RED" in res.json["reason"]


def test_warning_only_stop_does_not_count_as_block(repo, run_hook):
    _clean(repo)
    (repo / "wip.txt").write_text("x")
    for _ in range(4):
        assert not _blocked(run_hook("stop_gate", _payload(background_tasks=_bg("subagent")), repo))
    assert _blocked(run_hook("stop_gate", _payload(), repo))  # helpers gone: blocks again
