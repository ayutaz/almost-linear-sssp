"""Tests for .claude/hooks/session_start.py (SessionStart context summary)."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

from harness_util import GREEN_EVENT, RED_EVENT, write_task

HOOK = "session_start"
PAYLOAD = {"hook_event_name": "SessionStart", "source": "startup"}

TASKS_MD = """# Tasks

Intro text mentioning TODO in prose should not be picked when a table exists.

## Phase 0

| ID | Priority | Status | Dependencies | Acceptance Criteria | Test Evidence |
|---|---|---|---|---|---|
| P0-DONE-ONE | P0 | DONE | — | finished work | x.json |
| P0-TODO-ONE | P0 | TODO | — | write the tdd tool | y.json |
| P0-WIP-ONE | P1 | IN_PROGRESS | — | hooks in progress | z.json |

## Phase 1

| ID | Priority | Status | Dependencies | Acceptance Criteria | Test Evidence |
|---|---|---|---|---|---|
| P1-BLOCKED-ONE | P0 | BLOCKED | P0-TODO-ONE | waiting on paper | w.json |
| P1-DONE-TWO | P0 | DONE | — | done again | v.json |
"""


def _fake_bin(tmp_path: Path, name: str, body: str) -> Path:
    d = tmp_path / "fakebin"
    d.mkdir(exist_ok=True)
    p = d / name
    p.write_text("#!/bin/sh\n" + body)
    p.chmod(0o755)
    return d


def _env_with_path(bindir: Path) -> dict[str, str]:
    return {"PATH": f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}"}


def _add_remote(root: Path) -> None:
    subprocess.run(
        ["git", "-C", str(root), "remote", "add", "origin", "https://example.invalid/x.git"],
        check=True,
        capture_output=True,
    )


def test_emits_single_session_start_json(repo, run_hook):
    res = run_hook(HOOK, PAYLOAD, repo)
    assert res.code == 0
    data = json.loads(res.stdout)
    assert data["hookSpecificOutput"]["hookEventName"] == "SessionStart"
    assert "main" in res.context  # branch line from git status


def test_missing_progress_and_tasks(repo, run_hook):
    res = run_hook(HOOK, PAYLOAD, repo)
    assert res.code == 0
    assert "docs/PROGRESS.md: なし" in res.context
    assert "docs/TASKS.md: なし" in res.context


def test_progress_auto_block_is_shown(repo, run_hook):
    (repo / "docs/PROGRESS.md").write_text(
        "# Progress\n\nintro line\n\n<!-- AUTO:precompact -->\n"
        "auto line one\nauto line two\n<!-- /AUTO -->\n"
    )
    res = run_hook(HOOK, PAYLOAD, repo)
    assert "auto line one" in res.context
    assert "auto line two" in res.context
    assert "intro line" not in res.context


def test_progress_without_block_shows_head(repo, run_hook):
    lines = [f"progress line {i}" for i in range(1, 15)]
    (repo / "docs/PROGRESS.md").write_text("\n\n".join(lines) + "\n")
    res = run_hook(HOOK, PAYLOAD, repo)
    assert "progress line 1" in res.context
    assert "progress line 8" in res.context
    assert "progress line 9" not in res.context


def test_tasks_table_extracts_open_rows(repo, run_hook):
    (repo / "docs/TASKS.md").write_text(TASKS_MD)
    res = run_hook(HOOK, PAYLOAD, repo)
    ctx = res.context
    assert "P0-TODO-ONE" in ctx
    assert "P0-WIP-ONE" in ctx
    assert "P1-BLOCKED-ONE" in ctx
    assert "IN_PROGRESS" in ctx
    assert "write the tdd tool" in ctx
    assert "P0-DONE-ONE" not in ctx
    assert "P1-DONE-TWO" not in ctx
    assert "3件" in ctx
    assert "Intro text" not in ctx


def test_tasks_without_table_falls_back_to_keyword_lines(repo, run_hook):
    (repo / "docs/TASKS.md").write_text("# Tasks\n- TODO: alpha thing\n- DONE: beta thing\n")
    res = run_hook(HOOK, PAYLOAD, repo)
    assert "alpha thing" in res.context
    assert "beta thing" not in res.context


def test_tasks_rows_capped_at_ten(repo, run_hook):
    rows = "\n".join(f"| T-{i:02d} | TODO | item {i} |" for i in range(25))
    (repo / "docs/TASKS.md").write_text(f"| ID | Status | Summary |\n|---|---|---|\n{rows}\n")
    res = run_hook(HOOK, PAYLOAD, repo)
    assert "T-09" in res.context
    assert "T-10" not in res.context
    assert "25件" in res.context


def test_active_tdd_tasks_shown(repo, run_hook):
    write_task(
        repo,
        {
            "id": "P9-ACTIVE",
            "phase": "GREEN",
            "files": {"impl": ["src/core/**"], "tests": []},
            "events": [RED_EVENT, GREEN_EVENT],
        },
    )
    res = run_hook(HOOK, PAYLOAD, repo)
    assert "P9-ACTIVE" in res.context
    assert "GREEN" in res.context
    assert "src/core/**" in res.context


def test_no_active_tdd_tasks(repo, run_hook):
    res = run_hook(HOOK, PAYLOAD, repo)
    assert "作業中 TDD タスク: なし" in res.context


def test_missing_tools_listed(repo, run_hook, tmp_path):
    bindir = _fake_bin(tmp_path, "placeholder", "exit 0\n")
    restricted = f"{bindir}{os.pathsep}/usr/bin{os.pathsep}/bin"
    res = run_hook(HOOK, PAYLOAD, repo, env={"PATH": restricted})
    assert res.code == 0
    assert "ツール" in res.context
    for tool in ("cmake", "ninja", "uv"):
        if shutil.which(tool, path=restricted) is None:
            assert tool in res.context.split("ツール", 1)[1].splitlines()[0]


def test_gh_not_called_without_remote(repo, run_hook, tmp_path):
    log = tmp_path / "gh.log"
    bindir = _fake_bin(tmp_path, "gh", f'echo "$@" >> "{log}"\necho "[]"\n')
    res = run_hook(HOOK, PAYLOAD, repo, env=_env_with_path(bindir))
    assert res.code == 0
    assert not log.exists()
    assert "CI" in res.context
    assert "remote" in res.context


def test_gh_failure_reported(repo, run_hook, tmp_path):
    _add_remote(repo)
    log = tmp_path / "gh.log"
    run = [
        {
            "status": "completed",
            "conclusion": "failure",
            "name": "ci",
            "headBranch": "main",
            "url": "https://example.invalid/run/1",
            "createdAt": "2026-11-27T00:00:00Z",
        }
    ]
    bindir = _fake_bin(tmp_path, "gh", f'echo "$@" >> "{log}"\necho \'{json.dumps(run)}\'\n')
    res = run_hook(HOOK, PAYLOAD, repo, env=_env_with_path(bindir))
    assert res.code == 0
    assert "run list" in log.read_text()
    assert "失敗" in res.context
    assert "failure" in res.context
    assert "https://example.invalid/run/1" in res.context


def test_gh_success_one_line(repo, run_hook, tmp_path):
    _add_remote(repo)
    run = [
        {
            "status": "completed",
            "conclusion": "success",
            "name": "ci",
            "headBranch": "main",
            "url": "https://example.invalid/run/2",
            "createdAt": "2026-11-27T00:00:00Z",
        }
    ]
    bindir = _fake_bin(tmp_path, "gh", f"echo '{json.dumps(run)}'\n")
    res = run_hook(HOOK, PAYLOAD, repo, env=_env_with_path(bindir))
    ci_lines = [ln for ln in res.context.splitlines() if "success" in ln]
    assert len(ci_lines) == 1
    assert "失敗" not in res.context


def test_gh_error_does_not_break(repo, run_hook, tmp_path):
    _add_remote(repo)
    bindir = _fake_bin(tmp_path, "gh", 'echo "not logged in" >&2\nexit 4\n')
    res = run_hook(HOOK, PAYLOAD, repo, env=_env_with_path(bindir))
    assert res.code == 0
    assert "CI" in res.context


def test_gh_timeout_finishes_quickly(repo, run_hook, tmp_path):
    _add_remote(repo)
    bindir = _fake_bin(tmp_path, "gh", "sleep 30\n")
    start = time.monotonic()
    res = run_hook(HOOK, PAYLOAD, repo, env=_env_with_path(bindir), timeout=20)
    elapsed = time.monotonic() - start
    assert res.code == 0
    assert elapsed < 8
    assert "タイムアウト" in res.context


def test_total_lines_capped(repo, run_hook):
    for i in range(30):
        (repo / f"new_{i}.txt").write_text("x")
    rows = "\n".join(f"| T-{i:02d} | TODO | item {i} |" for i in range(30))
    (repo / "docs/TASKS.md").write_text(f"| ID | Status | Summary |\n|---|---|---|\n{rows}\n")
    (repo / "docs/PROGRESS.md").write_text("\n".join(f"line {i}" for i in range(50)) + "\n")
    for i in range(15):
        write_task(repo, {"id": f"P9-T{i:02d}", "events": [RED_EVENT]})
    res = run_hook(HOOK, PAYLOAD, repo)
    assert res.code == 0
    assert len(res.context.splitlines()) <= 40
    m = re.search(r"変更 (\d+) 件", res.context)
    assert m and int(m.group(1)) >= 30  # changed-file count is still reported


def test_invalid_stdin_exit_zero(repo, run_hook):
    res = run_hook(HOOK, {}, repo, raw_stdin="{not json")
    assert res.code == 0
    if res.stdout.strip():
        json.loads(res.stdout)


def test_not_a_git_repo(tmp_path, run_hook):
    plain = tmp_path / "plain"
    plain.mkdir()
    res = run_hook(HOOK, PAYLOAD, plain)
    assert res.code == 0
    if res.stdout.strip():
        json.loads(res.stdout)
