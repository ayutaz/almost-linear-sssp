"""Shared fixtures for harness tests (hooks + scripts/tdd.py).

Hooks are exercised exactly like Claude Code runs them: a subprocess that gets
the event JSON on stdin; we check the exit code and the stdout JSON.
Plain helpers live in harness_util.py.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from harness_util import HOOKS_DIR, HookResult


def _git(root: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """A throwaway git repository with the project's directory layout."""
    root = tmp_path / "proj"
    for d in ("src/core", "include/sssp", "tests/unit", "docs/tdd/tasks", ".tdd"):
        (root / d).mkdir(parents=True, exist_ok=True)
    (root / ".gitignore").write_text("/.tdd/runtime/\n/build/\n")
    (root / ".tdd" / "config.json").write_text(json.dumps({"enforce": False}))
    (root / "src/core/a.cpp").write_text("int a() { return 1; }\n")
    (root / "tests/unit/test_a.cpp").write_text(
        'TEST_CASE("a") {\n  REQUIRE(a() == 1);\n  REQUIRE(a() != 2);\n}\n'
    )
    _git(root, "init", "-q", "-b", "main")
    _git(root, "config", "user.email", "t@example.invalid")  # public-safety: allow
    _git(root, "config", "user.name", "t")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "init")
    return root


@pytest.fixture
def run_hook() -> Callable[..., HookResult]:
    """run_hook("pre_tool_use", payload, root) -> HookResult.

    ``payload`` gets ``cwd`` (= root) and ``hook_event_name`` filled in when
    missing. ``CLAUDE_PROJECT_DIR`` is set to ``root``.
    """

    def _run(
        name: str,
        payload: dict[str, Any],
        root: Path,
        env: dict[str, str] | None = None,
        timeout: float = 60,
        raw_stdin: str | None = None,
    ) -> HookResult:
        payload = dict(payload)
        payload.setdefault("cwd", str(root))
        payload.setdefault("session_id", "test-session")
        full_env = dict(os.environ)
        full_env["CLAUDE_PROJECT_DIR"] = str(root)
        full_env.update(env or {})
        proc = subprocess.run(
            [sys.executable, str(HOOKS_DIR / f"{name}.py")],
            input=raw_stdin if raw_stdin is not None else json.dumps(payload),
            capture_output=True,
            text=True,
            cwd=str(root),
            env=full_env,
            timeout=timeout,
        )
        return HookResult(proc.returncode, proc.stdout, proc.stderr)

    return _run
