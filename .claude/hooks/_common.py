"""Shared helpers for the project's Claude Code hooks.

Every hook is a standalone script run as ``python3 <hook>.py`` (exec form, see
docs/design/harness.md). Only the standard library is used and the code stays
compatible with Python 3.9 so that the hooks work even when ``uv`` or a newer
interpreter is missing.

Contract (official spec: https://code.claude.com/docs/en/hooks):
- The event JSON arrives on stdin.
- Exit 0 + a single JSON object on stdout = structured result.
- Exit 0 + empty stdout = no opinion.
- Internal errors must never block work (fail-open), except security guards,
  which deny on error (fail-closed). Each hook decides which applies.
"""

from __future__ import annotations

import contextlib
import datetime as _dt
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

RUNTIME_REL = Path(".tdd") / "runtime"
CONFIG_REL = Path(".tdd") / "config.json"
TASKS_REL = Path("docs") / "tdd" / "tasks"

DEFAULT_CONFIG: dict[str, Any] = {
    "enforce": False,
    "impl_dirs": ["src/", "include/", "tools/"],
    "test_dirs": ["tests/"],
    "stop_gate": {
        "enabled": True,
        "check_failed_tests": True,
        "check_red_tasks": True,
        "check_uncommitted": True,
        "check_progress_updated": True,
        "max_consecutive_blocks": 3,
    },
    "post_check": {"enabled": True, "timeout_sec": 120, "preset": "debug"},
}

LLVM_CLANG_FORMAT = "/opt/homebrew/opt/llvm/bin/clang-format"


# --------------------------------------------------------------------------- I/O
def read_input() -> dict[str, Any]:
    """Read the hook event JSON from stdin. Returns {} on empty/invalid input."""
    try:
        raw = sys.stdin.read()
    except Exception:
        return {}
    if not raw.strip():
        return {}
    try:
        data = json.loads(raw)
    except ValueError:
        return {}
    return data if isinstance(data, dict) else {}


def emit(obj: dict[str, Any]) -> None:
    """Print exactly one JSON object on stdout (nothing else may be printed)."""
    sys.stdout.write(json.dumps(obj, ensure_ascii=False))
    sys.stdout.write("\n")
    sys.stdout.flush()


def emit_context(event: str, text: str, limit: int = 9000) -> None:
    """Emit ``hookSpecificOutput.additionalContext`` (capped below the 10k limit)."""
    if not text:
        return
    if len(text) > limit:
        text = text[: limit - 20] + "\n...(truncated)"
    emit({"hookSpecificOutput": {"hookEventName": event, "additionalContext": text}})


def debug(msg: str) -> None:
    """stderr is only visible in the debug log on exit 0; safe for diagnostics."""
    with contextlib.suppress(Exception):
        sys.stderr.write(f"[hook] {msg}\n")


# ------------------------------------------------------------------------- paths
def find_root(start: str | os.PathLike[str] | None) -> Path | None:
    """Nearest ancestor (inclusive) containing ``.git`` (dir or worktree file)."""
    if not start:
        return None
    p = Path(start).expanduser()
    with contextlib.suppress(OSError):
        p = p.resolve()
    if p.is_file() or not p.exists():
        p = p.parent
    for cand in [p, *p.parents]:
        if (cand / ".git").exists():
            return cand
    return None


def project_root(data: dict[str, Any] | None = None, path_hint: str | None = None) -> Path:
    """Repository root for this hook invocation.

    Order: git root of ``path_hint`` (e.g. the edited file, so worktrees work),
    git root of the event ``cwd``, ``$CLAUDE_PROJECT_DIR``, process cwd.
    """
    data = data or {}
    for hint in (path_hint, data.get("cwd")):
        root = find_root(hint) if hint else None
        if root is not None:
            return root
    env = os.environ.get("CLAUDE_PROJECT_DIR")
    if env:
        return Path(env).resolve()
    return Path.cwd().resolve()


def rel_to_root(path: str | os.PathLike[str], root: Path) -> str | None:
    """POSIX relative path of ``path`` inside ``root``; None when outside."""
    p = Path(path).expanduser()
    if not p.is_absolute():
        p = root / p
    try:
        resolved = Path(os.path.normpath(str(p)))
        with contextlib.suppress(OSError):
            resolved = resolved.resolve()
        rel = resolved.relative_to(root.resolve())
    except ValueError:
        return None
    return rel.as_posix()


def runtime_dir(root: Path) -> Path:
    d = root / RUNTIME_REL
    d.mkdir(parents=True, exist_ok=True)
    return d


# -------------------------------------------------------------------------- json
def load_json(path: Path, default: Any = None) -> Any:
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return default


def save_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp{os.getpid()}")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=2)
        fh.write("\n")
    os.replace(tmp, path)


def _merge(base: dict[str, Any], over: dict[str, Any]) -> dict[str, Any]:
    out = dict(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out


def load_config(root: Path) -> dict[str, Any]:
    user = load_json(root / CONFIG_REL, {})
    if not isinstance(user, dict):
        user = {}
    return _merge(DEFAULT_CONFIG, user)


# ------------------------------------------------------------------------- tasks
def load_tasks(root: Path) -> dict[str, dict[str, Any]]:
    """All evidence files docs/tdd/tasks/*.json keyed by id (invalid files skipped)."""
    out: dict[str, dict[str, Any]] = {}
    d = root / TASKS_REL
    if not d.is_dir():
        return out
    for f in sorted(d.glob("*.json")):
        obj = load_json(f)
        if isinstance(obj, dict):
            out[str(obj.get("id") or f.stem)] = obj
    return out


def load_active(root: Path) -> dict[str, dict[str, Any]]:
    """Active (in-progress) tasks registered by ``scripts/tdd.py start``.

    File: .tdd/runtime/active.json = {"tasks": {"<ID>": {...}}}.
    Entry fields used by hooks: ``allow_test_reduction`` (bool), ``reason``.
    """
    obj = load_json(root / RUNTIME_REL / "active.json", {})
    tasks = obj.get("tasks") if isinstance(obj, dict) else None
    return tasks if isinstance(tasks, dict) else {}


def last_event_phase(task: dict[str, Any]) -> str | None:
    events = task.get("events") or []
    if events and isinstance(events[-1], dict):
        return events[-1].get("phase")
    return None


def has_event(task: dict[str, Any], phase: str) -> bool:
    return any(isinstance(e, dict) and e.get("phase") == phase for e in task.get("events") or [])


# ------------------------------------------------------------------------- globs
def glob_to_regex(pattern: str) -> re.Pattern[str]:
    """Translate a path glob (``*``, ``?``, ``**``, ``[..]``) to an anchored regex.

    ``*`` and ``?`` never cross ``/``; ``**`` crosses directories and ``**/``
    may match zero directories.
    """
    pattern = pattern.replace("\\", "/")
    if pattern.startswith("./"):
        pattern = pattern[2:]
    i, n, out = 0, len(pattern), []
    while i < n:
        c = pattern[i]
        if c == "*":
            if pattern[i : i + 2] == "**":
                if pattern[i : i + 3] == "**/":
                    out.append("(?:.*/)?")
                    i += 3
                else:
                    out.append(".*")
                    i += 2
                continue
            out.append("[^/]*")
        elif c == "?":
            out.append("[^/]")
        elif c == "[":
            j = pattern.find("]", i + 1)
            if j == -1:
                out.append(re.escape(c))
            else:
                body = pattern[i + 1 : j]
                if body.startswith("!"):
                    body = "^" + body[1:]
                out.append("[" + body + "]")
                i = j + 1
                continue
        else:
            out.append(re.escape(c))
        i += 1
    return re.compile("^" + "".join(out) + "$")


def glob_match(rel_path: str, pattern: str) -> bool:
    if not pattern:
        return False
    rel_path = rel_path.replace("\\", "/")
    if glob_to_regex(pattern).match(rel_path):
        return True
    # A bare directory pattern ("src/foo/") matches everything below it.
    return pattern.endswith("/") and rel_path.startswith(pattern)


def any_glob(rel_path: str, patterns: list[str] | None) -> bool:
    return any(glob_match(rel_path, p) for p in patterns or [])


# ----------------------------------------------------------------------- process
def run(cmd: list[str], cwd: Path | str | None = None, timeout: float = 10.0) -> tuple[int, str]:
    """Run a command; returns (exit_code, combined_output). Never raises.

    Exit code 124 = timeout, 127 = executable not found.
    """
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            text=True,
            errors="replace",
        )
        return proc.returncode, proc.stdout or ""
    except subprocess.TimeoutExpired as exc:
        out = exc.output if isinstance(exc.output, str) else ""
        return 124, (out or "") + f"\n[timeout after {timeout}s]"
    except FileNotFoundError:
        return 127, f"{cmd[0]}: not found"
    except Exception as exc:  # pragma: no cover - defensive
        return 1, f"{type(exc).__name__}: {exc}"


def tail(text: str, lines: int = 20) -> str:
    return "\n".join((text or "").rstrip("\n").splitlines()[-lines:])


def which(name: str) -> str | None:
    if name == "clang-format" and os.access(LLVM_CLANG_FORMAT, os.X_OK):
        return LLVM_CLANG_FORMAT
    return shutil.which(name)


def git(root: Path, *args: str, timeout: float = 5.0) -> tuple[int, str]:
    return run(["git", "-C", str(root), *args], timeout=timeout)


# -------------------------------------------------------------------------- time
def now_iso() -> str:
    return _dt.datetime.now().astimezone().isoformat(timespec="seconds")
