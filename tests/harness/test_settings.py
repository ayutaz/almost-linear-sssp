""".claude/settings.json and .tdd/config.json wire the hooks the way the official spec expects."""

from __future__ import annotations

import json
import re

import pytest
from harness_util import HOOKS_DIR, REPO

SETTINGS = REPO / ".claude" / "settings.json"
CONFIG = REPO / ".tdd" / "config.json"
EVENTS = {
    "SessionStart": "session_start",
    "PreToolUse": "pre_tool_use",
    "PostToolUse": "post_tool_use",
    "PreCompact": "pre_compact",
    "Stop": "stop_gate",
}


def settings() -> dict:
    return json.loads(SETTINGS.read_text())


def handlers(event: str) -> list[tuple[str | None, dict]]:
    groups = settings()["hooks"][event]
    return [(g.get("matcher"), h) for g in groups for h in g["hooks"]]


@pytest.mark.parametrize(("event", "script"), sorted(EVENTS.items()))
def test_each_event_runs_its_script_in_exec_form(event: str, script: str) -> None:
    hs = handlers(event)
    assert len(hs) == 1
    _, h = hs[0]
    assert h["type"] == "command"
    # Exec form: no shell quoting issues; python3 works without uv.
    assert h["command"] == "python3"
    assert h["args"] == [f"${{CLAUDE_PROJECT_DIR}}/.claude/hooks/{script}.py"]
    assert (HOOKS_DIR / f"{script}.py").is_file()
    assert isinstance(h.get("timeout"), int) and h["timeout"] > 0


def _matcher_tools(matcher: str) -> set[str]:
    # Only letters/digits/_/|/, keeps the matcher an exact list (not a regex).
    assert re.fullmatch(r"[A-Za-z0-9_|, ]+", matcher), matcher
    return {t.strip() for t in re.split(r"[|,]", matcher)}


def test_pre_tool_use_matcher() -> None:
    ((matcher, _),) = handlers("PreToolUse")
    assert _matcher_tools(matcher) == {"Edit", "Write", "MultiEdit", "NotebookEdit", "Bash"}


def test_post_tool_use_matcher_and_timeout() -> None:
    ((matcher, h),) = handlers("PostToolUse")
    assert _matcher_tools(matcher) == {"Edit", "Write", "MultiEdit"}
    cfg = json.loads(CONFIG.read_text())
    post_timeout = cfg.get("post_check", {}).get("timeout_sec", 120)
    # build+test budget, plus formatting and waiting for the build lock
    assert h["timeout"] >= post_timeout + 30


@pytest.mark.parametrize("event", ["SessionStart", "PreCompact", "Stop"])
def test_matcherless_or_lifecycle_events(event: str) -> None:
    ((matcher, _),) = handlers(event)
    if event == "SessionStart":
        assert set(matcher.split("|")) == {"startup", "resume", "clear", "compact"}
    else:
        assert matcher in (None, "", "*")


@pytest.mark.parametrize(
    "rule",
    [
        "Bash(git push --force*)",
        "Bash(git push -f*)",
        "Bash(git reset --hard*)",
        "Bash(git clean -f*)",
        "Bash(git branch -D *)",
        "Bash(gh repo delete*)",
    ],
)
def test_dangerous_git_is_denied_by_permissions_too(rule: str) -> None:
    assert rule in settings()["permissions"]["deny"]


def test_permissions_do_not_remove_whole_tools() -> None:
    perms = settings()["permissions"]
    for rule in perms.get("deny", []):
        assert "(" in rule, f"bare tool deny would hide the tool: {rule}"
    assert not perms.get("allow"), "project settings must not pre-approve anything"


def test_tdd_enforcement_starts_disabled() -> None:
    cfg = json.loads(CONFIG.read_text())
    assert cfg["enforce"] is False


@pytest.mark.parametrize("script", sorted(EVENTS.values()))
def test_hooks_accept_minimal_input(repo, run_hook, script: str) -> None:
    """Any hook given just {} exits 0 and prints nothing or exactly one JSON object."""
    res = run_hook(script, {}, repo)
    assert res.code == 0, res.stderr
    out = res.stdout.strip()
    if out and script != "pre_compact":  # PreCompact stdout is plain-text instructions
        assert isinstance(json.loads(out), dict)
