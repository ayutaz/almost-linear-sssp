"""SessionStart hook: a short, factual project summary injected at session start.

Sections (each a few lines, 40 lines in total at most): git status, the
docs/PROGRESS.md AUTO block (or head), open rows of docs/TASKS.md, active TDD
tasks, missing tools and the latest CI run (only when a git remote and ``gh``
exist). Runs on every start/resume/clear/compact, so every external command has
a short timeout. Fail-open: always exit 0.

Python 3.9 compatible, standard library only.
"""

from __future__ import annotations

import contextlib
import json
import os
import re
import signal
import subprocess
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import _common  # noqa: E402

MAX_LINES = 40
MAX_LINE_CHARS = 160
GIT_FILES = 10
PROGRESS_BLOCK_LINES = 10
PROGRESS_HEAD_LINES = 8
TASK_ROWS = 10
ACTIVE_TASKS = 5
GH_TIMEOUT = 4.0
OPEN_STATUSES = ("TODO", "IN_PROGRESS", "BLOCKED")
CI_BAD = {"failure", "cancelled", "timed_out", "startup_failure", "action_required"}
AUTO_START = "<!-- AUTO:precompact -->"
AUTO_END = "<!-- /AUTO -->"
SUMMARY_HEADERS = (
    "summary",
    "概要",
    "title",
    "task",
    "description",
    "acceptance criteria",
    "内容",
)


def _clip(text: str, n: int = MAX_LINE_CHARS) -> str:
    text = text.rstrip()
    return text if len(text) <= n else text[: n - 3] + "..."


def run_group(cmd: list[str], cwd: Path, timeout: float) -> tuple[int, str]:
    """Like ``_common.run`` but kills the whole process group on timeout.

    Returns (code, stdout); 124 = timeout, 127 = not found.
    """
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            errors="replace",
            start_new_session=True,
        )
    except FileNotFoundError:
        return 127, ""
    except OSError as exc:
        return 1, str(exc)
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except OSError:
            proc.kill()
        with contextlib.suppress(Exception):
            proc.communicate(timeout=1)
        return 124, ""
    if proc.returncode != 0 and not out.strip():
        out = err
    return proc.returncode, out or ""


# ----------------------------------------------------------------------- sections
def section_git(root: Path) -> list[str]:
    code, out = _common.git(root, "status", "--short", "--branch", timeout=3.0)
    if code != 0:
        return ["git: 取得できず (git リポジトリではないか git が失敗)"]
    lines = out.splitlines()
    branch = "?"
    if lines and lines[0].startswith("## "):
        branch = lines[0][3:].strip()
        lines = lines[1:]
    entries = [ln.strip() for ln in lines if ln.strip()]
    result = [_clip(f"git: {branch} / 変更 {len(entries)} 件")]
    if entries:
        shown = ", ".join(entries[:GIT_FILES])
        more = " ..." if len(entries) > GIT_FILES else ""
        result.append(_clip("  " + shown, 300) + more)
    return result


def section_progress(root: Path) -> list[str]:
    path = root / "docs" / "PROGRESS.md"
    if not path.is_file():
        return ["docs/PROGRESS.md: なし"]
    text = path.read_text(encoding="utf-8", errors="replace")
    start = text.find(AUTO_START)
    end = text.find(AUTO_END, start + 1) if start >= 0 else -1
    if start >= 0 and end > start:
        body = text[start + len(AUTO_START) : end]
        block = [ln for ln in body.splitlines() if ln.strip()]
        if block:
            head = ["docs/PROGRESS.md (AUTO ブロック):"]
            return head + ["  " + _clip(ln) for ln in block[:PROGRESS_BLOCK_LINES]]
    rows = [ln for ln in text.splitlines() if ln.strip() and not ln.strip().startswith("<!--")]
    head = ["docs/PROGRESS.md (冒頭):"]
    return head + ["  " + _clip(ln) for ln in rows[:PROGRESS_HEAD_LINES]]


_SEP_RE = re.compile(r"^\s*\|?\s*:?-{2,}:?\s*(\|\s*:?-{2,}:?\s*)*\|?\s*$")
_CELL_SPLIT = re.compile(r"(?<!\\)\|")


def _cells(line: str) -> list[str]:
    s = line.strip()
    if s.startswith("|"):
        s = s[1:]
    if s.endswith("|") and not s.endswith("\\|"):
        s = s[:-1]
    return [c.strip() for c in _CELL_SPLIT.split(s)]


def _open_rows_from_tables(lines: list[str]) -> list[tuple[str, str, str]] | None:
    """(id, status, summary) of open rows; None when no table has a Status column."""
    found_table = False
    rows: list[tuple[str, str, str]] = []
    i = 0
    while i < len(lines) - 1:
        line = lines[i]
        if line.strip().startswith("|") and _SEP_RE.match(lines[i + 1]):
            header = [h.lower() for h in _cells(line)]
            status_col = next((k for k, h in enumerate(header) if h in ("status", "状態")), None)
            i += 2
            if status_col is None:
                continue
            found_table = True
            id_col = header.index("id") if "id" in header else 0
            sum_col = next((header.index(h) for h in SUMMARY_HEADERS if h in header), None)
            if sum_col is None:
                sum_col = next(
                    (k for k in range(len(header)) if k not in (id_col, status_col)), id_col
                )
            while i < len(lines) and lines[i].strip().startswith("|"):
                cells = _cells(lines[i])
                i += 1
                if status_col >= len(cells):
                    continue
                status = cells[status_col].strip("`* ").upper()
                if status in OPEN_STATUSES:
                    tid = cells[id_col] if id_col < len(cells) else "?"
                    summ = cells[sum_col] if sum_col < len(cells) else ""
                    rows.append((tid, status, summ))
            continue
        i += 1
    return rows if found_table else None


def section_tasks(root: Path) -> list[str]:
    path = root / "docs" / "TASKS.md"
    if not path.is_file():
        return ["docs/TASKS.md: なし"]
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    rows = _open_rows_from_tables(lines)
    if rows is None:
        pat = re.compile(r"\b(" + "|".join(OPEN_STATUSES) + r")\b")
        hits = [ln.strip() for ln in lines if pat.search(ln)]
        out = [f"docs/TASKS.md: 表なし、TODO/IN_PROGRESS/BLOCKED を含む行 {len(hits)}件"]
        return out + ["  " + _clip(ln, 120) for ln in hits[:TASK_ROWS]]
    counts = {s: sum(1 for r in rows if r[1] == s) for s in OPEN_STATUSES}
    detail = ", ".join(f"{s} {counts[s]}" for s in OPEN_STATUSES)
    out = [f"docs/TASKS.md: 未完了 {len(rows)}件 ({detail})"]
    for tid, status, summ in rows[:TASK_ROWS]:
        out.append("  " + _clip(f"{tid} [{status}] {_clip(summ, 70)}"))
    if len(rows) > TASK_ROWS:
        out.append(f"  ...ほか {len(rows) - TASK_ROWS}件")
    return out


def section_active(root: Path) -> list[str]:
    active = _common.load_active(root)
    if not active:
        return ["作業中 TDD タスク: なし"]
    tasks = _common.load_tasks(root)
    out = [f"作業中 TDD タスク: {len(active)}件"]
    for tid in list(active)[:ACTIVE_TASKS]:
        task = tasks.get(tid)
        if task is None:
            out.append(f"  {tid}: 証跡ファイルなし")
            continue
        impl = (task.get("files") or {}).get("impl") or []
        impl_s = ", ".join(str(p) for p in impl) if impl else "-"
        out.append(
            "  "
            + _clip(
                f"{tid}: phase={task.get('phase')}, "
                f"最終イベント={_common.last_event_phase(task) or '-'}, impl={impl_s}"
            )
        )
    if len(active) > ACTIVE_TASKS:
        out.append(f"  ...ほか {len(active) - ACTIVE_TASKS}件")
    return out


def section_tools() -> list[str]:
    names = ("cmake", "ninja", "uv", "clang-format")
    missing = [n for n in names if not _common.which(n)]
    if not missing:
        return ["ツール: " + ", ".join(names) + " はすべて利用可"]
    return ["ツール: 見つからない " + ", ".join(missing) + " (それ以外は利用可)"]


def _format_run(run: dict[str, Any]) -> str:
    status = run.get("status") or "?"
    concl = run.get("conclusion") or ""
    where = f"{run.get('name') or '?'} on {run.get('headBranch') or '?'}"
    when = run.get("createdAt") or ""
    url = run.get("url") or ""
    if status != "completed":
        label = f"実行中 ({status})"
    elif concl == "success":
        label = "成功 (success)"
    elif concl in CI_BAD:
        label = f"失敗 ({concl})"
    else:
        label = f"終了 ({concl or 'unknown'})"
    return _clip(f"前回 CI: {label} {where} {when} {url}".rstrip(), 220)


def section_ci(root: Path) -> list[str]:
    code, remotes = _common.git(root, "remote", timeout=2.0)
    if code != 0 or not remotes.strip():
        return ["前回 CI: git remote なしのため未確認"]
    gh = _common.which("gh")
    if not gh:
        return ["前回 CI: gh が無いため未確認"]
    cmd = [
        gh,
        "run",
        "list",
        "--limit",
        "1",
        "--json",
        "status,conclusion,name,headBranch,url,createdAt",
    ]
    code, out = run_group(cmd, root, GH_TIMEOUT)
    if code == 124:
        return [f"前回 CI: gh run list がタイムアウト ({GH_TIMEOUT:g} 秒) のため未確認"]
    if code != 0:
        first = (out.strip().splitlines() or [""])[0]
        return [_clip(f"前回 CI: gh run list が失敗 (exit {code}: {first})")]
    try:
        runs = json.loads(out)
    except ValueError:
        return ["前回 CI: gh の出力を解釈できず"]
    if not isinstance(runs, list) or not runs or not isinstance(runs[0], dict):
        return ["前回 CI: 実行履歴なし"]
    return [_format_run(runs[0])]


# --------------------------------------------------------------------------- main
def build_context(data: dict[str, Any]) -> str:
    root = _common.project_root(data)
    source = data.get("source") or "unknown"
    lines = [f"プロジェクト状況の要約 (SessionStart, source={source}):"]
    sections = (
        lambda: section_git(root),
        lambda: section_progress(root),
        lambda: section_tasks(root),
        lambda: section_active(root),
        section_tools,
        lambda: section_ci(root),
    )
    for fn in sections:
        try:
            lines.extend(fn())
        except Exception as exc:  # one broken section must not hide the others
            _common.debug(f"session_start section failed: {type(exc).__name__}: {exc}")
    if len(lines) > MAX_LINES:
        lines = lines[: MAX_LINES - 1] + ["...(省略)"]
    return "\n".join(lines)


def main() -> int:
    try:
        data = _common.read_input()
        _common.emit_context("SessionStart", build_context(data))
    except Exception as exc:  # fail-open
        _common.debug(f"session_start failed: {type(exc).__name__}: {exc}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
