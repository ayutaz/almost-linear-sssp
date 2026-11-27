"""PostToolUse hook (matcher: Edit|Write|MultiEdit).

1. Format the edited file only: C++ with clang-format, Python with ruff
   (``ruff format`` then ``ruff check``; findings are reported, not fixed).
2. For C++ changes, run a light check of one test category:
   ``cmake --build build/<preset> --target sssp_<cat>_tests`` then
   ``ctest -L ^<cat>$ -LE bench``. Benchmarks are never built or run.
3. Store the result in ``.tdd/runtime/last_check.json`` (read by the Stop hook).

The tool has already run, so this hook cannot block; problems are reported
through ``additionalContext``. Any internal error is fail-open (exit 0).
Python 3.9 compatible, standard library only.
"""

from __future__ import annotations

import contextlib
import fcntl
import os
import shlex
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import _common  # noqa: E402

EVENT = "PostToolUse"
CPP_EXT = {".cpp", ".cc", ".cxx", ".hpp", ".hh", ".h", ".ipp"}
PY_EXT = {".py"}
TEST_CATEGORIES = ("unit", "integration", "differential", "property", "regression")
UNIT_DIRS = ("src/", "include/", "tools/")
FORMAT_TIMEOUT = 20.0
LOCK_WAIT_DEFAULT = 10.0
CONTEXT_TAIL_LINES = 12


# ----------------------------------------------------------------------- process
def run_group(cmd: list[str], cwd: Path, timeout: float) -> tuple[int, str, bool]:
    """Run ``cmd`` in its own process group; kill the whole group on timeout.

    Returns (exit_code, combined_output, timed_out). Never raises.
    ``subprocess.run`` alone would wait for grandchildren that keep the pipe
    open, so the group is killed explicitly.
    """
    if timeout <= 0:
        return 124, "", True
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            start_new_session=True,
        )
    except FileNotFoundError:
        return 127, f"{cmd[0]}: not found", False
    except Exception as exc:
        return 1, f"{type(exc).__name__}: {exc}", False
    try:
        out, _ = proc.communicate(timeout=timeout)
        return proc.returncode, out or "", False
    except subprocess.TimeoutExpired:
        with contextlib.suppress(Exception):
            os.killpg(proc.pid, signal.SIGKILL)
        out = ""
        with contextlib.suppress(Exception):
            out, _ = proc.communicate(timeout=2)
        return 124, (out or "") + f"\n[timeout after {timeout:.0f}s]", True


def cmd_str(cmd: list[str]) -> str:
    return " ".join(shlex.quote(c) for c in cmd)


# ---------------------------------------------------------------------- tools
def _env_tool(var: str) -> tuple[bool, str | None]:
    """(set, path) for an override env var; set-but-unusable means 'missing'."""
    val = os.environ.get(var)
    if not val:
        return False, None
    return True, (val if os.access(val, os.X_OK) else None)


def clang_format_cmd() -> list[str] | None:
    is_set, path = _env_tool("TDD_CLANG_FORMAT")
    if is_set:
        return [path] if path else None
    found = _common.which("clang-format")
    return [found] if found else None


def ruff_cmd(root: Path) -> list[str] | None:
    is_set, path = _env_tool("TDD_RUFF")
    if is_set:
        return [path] if path else None
    venv = root / ".venv" / "bin" / "ruff"
    if os.access(str(venv), os.X_OK):
        return [str(venv)]
    found = _common.which("ruff")
    if found:
        return [found]
    uv = _common.which("uv")
    if uv:
        return [uv, "run", "--quiet", "ruff"]
    return None


# ----------------------------------------------------------------- formatting
def format_cpp(path: Path, root: Path) -> list[str]:
    base = clang_format_cmd()
    if not base:
        return []
    cmd = [*base, "-i", "--style=file", str(path)]
    code, out, timed_out = run_group(cmd, root, FORMAT_TIMEOUT)
    if code == 0:
        return []
    what = "タイムアウト" if timed_out else f"exit {code}"
    return [f"clang-format が失敗しました ({what}): {cmd_str(cmd)}", _common.tail(out, 8)]


def format_python(path: Path, root: Path) -> list[str]:
    base = ruff_cmd(root)
    if not base:
        return []
    notes: list[str] = []
    fmt = [*base, "format", str(path)]
    code, out, timed_out = run_group(fmt, root, FORMAT_TIMEOUT)
    if code != 0:
        what = "タイムアウト" if timed_out else f"exit {code}"
        notes += [f"ruff format が失敗しました ({what}): {cmd_str(fmt)}", _common.tail(out, 8)]
    chk = [*base, "check", str(path)]
    code, out, timed_out = run_group(chk, root, FORMAT_TIMEOUT)
    if code != 0:
        what = "タイムアウト" if timed_out else f"exit {code}"
        notes += [f"ruff check の指摘があります ({what}): {cmd_str(chk)}", _common.tail(out, 10)]
    return notes


# ---------------------------------------------------------------- light check
def category_for(rel: str) -> str | None:
    parts = rel.split("/")
    if len(parts) >= 3 and parts[0] == "tests" and parts[1] in TEST_CATEGORIES:
        return parts[1]
    if rel.startswith(UNIT_DIRS):
        return "unit"
    return None


def has_cpp_tests(root: Path, cat: str) -> bool:
    d = root / "tests" / cat
    if not d.is_dir():
        return False
    return next(iter(d.rglob("*.cpp")), None) is not None


def total_timeout(cfg: dict[str, Any]) -> float:
    for raw in (os.environ.get("TDD_POSTCHECK_TIMEOUT"), cfg.get("timeout_sec")):
        if raw is None or raw == "":
            continue
        with contextlib.suppress(TypeError, ValueError):
            val = float(raw)
            if val > 0:
                return val
    return 120.0


def lock_wait() -> float:
    with contextlib.suppress(TypeError, ValueError):
        return max(0.0, float(os.environ.get("TDD_POSTCHECK_LOCK_WAIT", "")))
    return LOCK_WAIT_DEFAULT


def acquire_lock(path: Path, wait: float) -> Any:
    """Open + flock ``path`` exclusively, polling up to ``wait`` seconds."""
    fh = open(path, "a+")  # noqa: SIM115 - held for the duration of the check
    deadline = time.monotonic() + wait
    while True:
        try:
            fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return fh
        except OSError:
            if time.monotonic() >= deadline:
                fh.close()
                return None
            time.sleep(0.1)


def run_check(root: Path, bdir: Path, cat: str, rel: str, budget: float) -> dict[str, Any]:
    start = time.monotonic()
    cmake = _common.which("cmake") or "cmake"
    ctest = _common.which("ctest") or "ctest"
    build = [cmake, "--build", str(bdir), "--target", f"sssp_{cat}_tests"]
    test = [
        ctest,
        "--test-dir",
        str(bdir),
        "-L",
        f"^{cat}$",
        "-LE",
        "bench",
        "--output-on-failure",
        "--timeout",
        "60",
    ]
    entry: dict[str, Any] = {"file": rel}

    lock_path = _common.runtime_dir(root) / "postcheck.lock"
    fh = acquire_lock(lock_path, min(lock_wait(), budget))
    if fh is None:
        entry.update(
            status="skipped",
            command=cmd_str(build),
            exit_code=0,
            output_tail="",
            reason="busy",
        )
        return entry
    try:
        remaining = budget - (time.monotonic() - start)
        code, out, timed_out = run_group(build, root, remaining)
        if timed_out:
            entry.update(
                status="timeout",
                command=cmd_str(build),
                exit_code=code,
                output_tail=_common.tail(out, 20),
                reason=f"build exceeded {budget:.0f}s",
            )
            return entry
        if code != 0:
            entry.update(
                status="fail",
                command=cmd_str(build),
                exit_code=code,
                output_tail=_common.tail(out, 20),
                reason="build failed",
            )
            return entry
        remaining = budget - (time.monotonic() - start)
        code, out, timed_out = run_group(test, root, remaining)
        if timed_out:
            entry.update(
                status="timeout",
                command=cmd_str(test),
                exit_code=code,
                output_tail=_common.tail(out, 20),
                reason=f"ctest exceeded the remaining time of {budget:.0f}s",
            )
            return entry
        entry.update(
            status="pass" if code == 0 else "fail",
            command=f"{cmd_str(build)} && {cmd_str(test)}" if code == 0 else cmd_str(test),
            exit_code=code,
            output_tail=_common.tail(out, 20),
        )
        if code != 0:
            entry["reason"] = "tests failed"
        return entry
    finally:
        with contextlib.suppress(Exception):
            fcntl.flock(fh, fcntl.LOCK_UN)
            fh.close()


def save_result(root: Path, cat: str, entry: dict[str, Any]) -> None:
    path = _common.runtime_dir(root) / "last_check.json"
    data = _common.load_json(path, {})
    if not isinstance(data, dict):
        data = {}
    cats = data.get("categories")
    if not isinstance(cats, dict):
        cats = {}
    now = _common.now_iso()
    entry = dict(entry, timestamp=now)
    cats[cat] = entry
    data["categories"] = cats
    data["updated"] = now
    data["last"] = dict(entry, category=cat)
    _common.save_json(path, data)


def describe(cat: str, entry: dict[str, Any]) -> list[str]:
    status = entry.get("status")
    if status == "fail":
        head = (
            f"軽量チェック {cat}: 失敗 ({entry.get('reason', '')}, exit {entry.get('exit_code')})"
        )
    elif status == "timeout":
        head = f"軽量チェック {cat}: タイムアウト ({entry.get('reason', '')})"
    else:
        return []
    lines = [head, f"対象ファイル: {entry.get('file')}", f"再現コマンド: {entry.get('command')}"]
    tail = _common.tail(entry.get("output_tail") or "", CONTEXT_TAIL_LINES)
    if tail:
        lines += ["出力末尾:", tail]
    lines.append("結果の記録: .tdd/runtime/last_check.json")
    return lines


# ----------------------------------------------------------------------- main
def handle(data: dict[str, Any]) -> list[str]:
    tool_input = data.get("tool_input")
    if not isinstance(tool_input, dict):
        return []
    file_path = tool_input.get("file_path")
    if not isinstance(file_path, str) or not file_path:
        return []
    root = _common.project_root(data, path_hint=file_path)
    rel = _common.rel_to_root(file_path, root)
    if rel is None or rel.startswith(("../", ".git/")):
        return []
    path = root / rel
    if not path.is_file():
        return []
    ext = path.suffix.lower()
    notes: list[str] = []

    if ext in PY_EXT:
        return format_python(path, root)
    if ext not in CPP_EXT:
        return []
    notes += format_cpp(path, root)

    cfg = _common.load_config(root).get("post_check") or {}
    if not isinstance(cfg, dict) or not cfg.get("enabled", True):
        return notes
    if os.environ.get("TDD_POSTCHECK") == "0":
        return notes
    cat = category_for(rel)
    if cat is None or not has_cpp_tests(root, cat):
        return notes
    preset = str(cfg.get("preset") or "debug")
    bdir = root / "build" / preset
    if not (bdir / "build.ninja").is_file():
        return notes
    entry = run_check(root, bdir, cat, rel, total_timeout(cfg))
    save_result(root, cat, entry)
    return notes + describe(cat, entry)


def main() -> int:
    try:
        data = _common.read_input()
        notes = [n for n in handle(data) if n]
        if notes:
            _common.emit_context(EVENT, "PostToolUse:\n" + "\n".join(notes))
    except Exception as exc:  # fail-open
        _common.debug(f"post_tool_use: {type(exc).__name__}: {exc}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
