#!/usr/bin/env python3
"""PreCompact hook: save a state snapshot before the context is compacted.

Writes a Markdown snapshot (active tasks, TDD phase and next command, open
failures, last commit, paper-mapping excerpt) into the ``<!-- AUTO:precompact -->``
block of docs/PROGRESS.md and into .tdd/runtime/precompact_snapshot.md, then
prints a short plain-text note that Claude Code uses as extra summary
instructions.

Never blocks compaction: always exits 0 (exit 2 / ``decision: block`` would
stop the compaction), and internal errors are only logged to stderr.
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import _common  # noqa: E402

BEGIN = "<!-- AUTO:precompact -->"
END = "<!-- /AUTO -->"
STDOUT_LIMIT = 2000
DONE_SHOW = 10
TEST_CMD = "ctest --preset debug"


# ------------------------------------------------------------------- collectors
def _done_ids(tasks: dict[str, dict[str, Any]]) -> list[str]:
    return [tid for tid, t in tasks.items() if str(t.get("phase", "")).upper() == "DONE"]


def _impl_globs(task: dict[str, Any]) -> list[str]:
    files = task.get("files")
    impl = files.get("impl") if isinstance(files, dict) else None
    return [str(x) for x in impl] if isinstance(impl, list) else []


def _failures(root: Path) -> list[tuple[str, dict[str, Any]]]:
    obj = _common.load_json(root / _common.RUNTIME_REL / "last_check.json", {})
    cats = obj.get("categories") if isinstance(obj, dict) else None
    if not isinstance(cats, dict):
        return []
    out = []
    for name, info in sorted(cats.items()):
        if isinstance(info, dict) and info.get("status") in ("fail", "timeout"):
            out.append((str(name), info))
    return out


def _cells(line: str) -> list[str]:
    body = line.strip()
    if body.startswith("|"):
        body = body[1:]
    if body.endswith("|"):
        body = body[:-1]
    return [c.strip() for c in body.split("|")]


def _is_separator(line: str) -> bool:
    cells = _cells(line)
    return bool(cells) and all(re.fullmatch(r":?-{1,}:?", c) for c in cells)


def _blocked_table_rows(lines: list[str]) -> list[str] | None:
    """Rows whose Status cell is exactly BLOCKED; None when no table has a Status column."""
    found_table = False
    out: list[str] = []
    status_col: int | None = None
    for i, ln in enumerate(lines):
        if not ln.lstrip().startswith("|"):
            status_col = None
            continue
        nxt = lines[i + 1] if i + 1 < len(lines) else ""
        if status_col is None:
            # a header row is followed by a |---|---| separator
            if nxt.lstrip().startswith("|") and _is_separator(nxt):
                heads = [c.lower() for c in _cells(ln)]
                if "status" in heads:
                    status_col = heads.index("status")
                    found_table = True
            continue
        if _is_separator(ln):
            continue
        cells = _cells(ln)
        if status_col < len(cells) and cells[status_col] == "BLOCKED":
            out.append(ln.strip())
    return out if found_table else None


def _blocked_lines(root: Path, limit: int = 5) -> list[str]:
    path = root / "docs" / "TASKS.md"
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    lines = text.splitlines()
    rows = _blocked_table_rows(lines)
    if rows is not None:
        # The task table is authoritative: only its Status column decides.
        return rows[:limit]
    out = []  # no task table: fall back to free-form lines ("- [ ] P9-X BLOCKED: ...")
    for ln in lines:
        if not re.search(r"\bBLOCKED\b", ln):
            continue
        # A legend such as "Status: TODO / IN_PROGRESS / BLOCKED / DONE" is not a task.
        others = sum(bool(re.search(rf"\b{w}\b", ln)) for w in ("TODO", "IN_PROGRESS", "DONE"))
        if others >= 2:
            continue
        out.append(ln.strip())
    return out[:limit]


def _paper_mapping(root: Path, limit: int = 10) -> list[str]:
    text = None
    for name in ("PAPER_MAPPING.md", "PAPER_MAPPING.ja.md"):
        try:
            text = (root / "docs" / name).read_text(encoding="utf-8", errors="replace")
            break
        except OSError:
            continue
    if text is None:
        return []
    return [ln.rstrip() for ln in text.splitlines() if ln.strip()][:limit]


def _last_commit(root: Path) -> str:
    code, out = _common.git(root, "log", "-1", "--format=%h %s (%cr)")
    out = out.strip()
    return out if code == 0 and out else "(none)"


def _next_commands(tid: str, task: dict[str, Any] | None) -> list[str]:
    tdd = "uv run python scripts/tdd.py"
    last = _common.last_event_phase(task) if task else None
    last = str(last).upper() if last else None
    if task is None or not _common.has_event(task, "RED"):
        return [f'{tdd} red {tid} --reason "..." -- <test cmd>']
    if last == "RED":
        return [f"Implement, then `{tdd} green {tid} -- {TEST_CMD}`"]
    if last in ("GREEN", "REFACTOR"):
        return [f"{tdd} refactor {tid} -- {TEST_CMD}", f"{tdd} done {tid}"]
    return [f"{tdd} status"]


# ------------------------------------------------------------------- rendering
def build_snapshot(root: Path, data: dict[str, Any]) -> tuple[str, str]:
    """Return (markdown_block_body, short_stdout_text)."""
    tasks = _common.load_tasks(root)
    active = _common.load_active(root)
    failures = _failures(root)
    unresolved = (root / ".tdd" / "runtime" / "stop_gate_unresolved.md").is_file()
    blocked = _blocked_lines(root)
    trigger = data.get("trigger")
    trigger = trigger if isinstance(trigger, str) and trigger else "unknown"

    md: list[str] = [
        "## Automatic snapshot (PreCompact)",
        "",
        f"- Generated: {_common.now_iso()}",
        f"- trigger: {trigger}",
        "",
        "### Completed tasks",
        "",
    ]
    done = _done_ids(tasks)
    if not done:
        md.append("- none")
    elif len(done) <= DONE_SHOW:
        md.append("- " + ", ".join(done))
    else:
        md.append(f"- {len(done)} tasks (last {DONE_SHOW}: {', '.join(done[-DONE_SHOW:])})")

    md += ["", "### Active tasks", ""]
    if not active:
        md.append("- none")
    for tid in sorted(active):
        task = tasks.get(tid)
        phase = task.get("phase") if task else "(no evidence file)"
        last = _common.last_event_phase(task) if task else None
        globs = ", ".join(_impl_globs(task)) if task else ""
        md.append(
            f"- {tid}: phase={phase}, last event={last or 'none'}, impl={globs or '(not set)'}"
        )

    md += ["", "### Open problems", ""]
    if not (failures or unresolved or blocked):
        md.append("- none")
    for name, info in failures:
        md.append(
            f"- Failing check `{name}` ({info.get('status')}, exit={info.get('exit_code')}): "
            f"`{info.get('command', '')}`"
        )
        tail = _common.tail(str(info.get("output_tail") or ""), 5)
        if tail:
            md += ["", "  ```", *("  " + ln for ln in tail.splitlines()), "  ```", ""]
    if unresolved:
        md.append("- Stop gate left unresolved issues: see .tdd/runtime/stop_gate_unresolved.md")
    for ln in blocked:
        md.append(f"- TASKS.md: {ln}")

    md += ["", "### TDD state and next commands", ""]
    for tid in sorted(active):
        for cmd in _next_commands(tid, tasks.get(tid)):
            md.append(f"- {tid}: `{cmd}`" if not cmd.startswith("Implement") else f"- {tid}: {cmd}")
    for name, info in failures:
        if info.get("command"):
            md.append(f"- Reproduce the failure ({name}): `{info.get('command')}`")
    md.append("- Finally run `make verify`")

    md += ["", "### Last commit", "", f"- {_last_commit(root)}", ""]
    md += ["### PAPER_MAPPING.md excerpt", ""]
    pm = _paper_mapping(root)
    if pm:
        md += ["```", *pm, "```"]
    else:
        md.append("- none")
    body = "\n".join(md).rstrip() + "\n"

    # ---------------------------------------------- short plain text for stdout
    out: list[str] = ["コンパクション後も次の状態を要約に残すこと。"]
    if active:
        items = []
        for tid in sorted(active):
            task = tasks.get(tid)
            last = _common.last_event_phase(task) if task else None
            items.append(f"{tid}({last or 'RED未記録'})")
        out.append("作業中タスク: " + ", ".join(items))
    else:
        out.append("作業中タスク: なし")
    if failures:
        out.append(
            "未解決の失敗: " + ", ".join(f"{n}: {i.get('command', '')}" for n, i in failures)
        )
    if unresolved:
        out.append("Stop gate 未解決レポート: .tdd/runtime/stop_gate_unresolved.md")
    out.append("詳細は docs/PROGRESS.md の AUTO ブロック(<!-- AUTO:precompact -->)を参照。")
    text = "\n".join(out)
    if len(text) > STDOUT_LIMIT:
        tail_line = out[-1]
        text = text[: STDOUT_LIMIT - len(tail_line) - 10] + "...\n" + tail_line
    return body, text


def replace_block(text: str, body: str) -> str:
    block = f"{BEGIN}\n{body}{END}"
    start = text.find(BEGIN)
    if start != -1:
        end = text.find(END, start)
        if end != -1:
            return text[:start] + block + text[end + len(END) :]
    sep = "" if text.endswith("\n\n") else ("\n" if text.endswith("\n") else "\n\n")
    if not text:
        sep = ""
    return text + sep + block + "\n"


def write_progress(root: Path, body: str) -> None:
    path = root / "docs" / "PROGRESS.md"
    try:
        text = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        text = "# Progress\n\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    new = replace_block(text, body)
    tmp = path.with_name(path.name + f".tmp{os.getpid()}")
    tmp.write_text(new, encoding="utf-8")
    os.replace(tmp, path)


def main() -> int:
    data = _common.read_input()
    root = _common.project_root(data)
    body, text = build_snapshot(root, data)
    try:
        write_progress(root, body)
    except Exception as exc:
        _common.debug(f"pre_compact: PROGRESS.md not written: {exc}")
    try:
        snap = _common.runtime_dir(root) / "precompact_snapshot.md"
        snap.write_text(body, encoding="utf-8")
    except Exception as exc:
        _common.debug(f"pre_compact: snapshot not written: {exc}")
    sys.stdout.write(text + "\n")
    sys.stdout.flush()
    return 0


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # fail-open: never stop the compaction
        _common.debug(f"pre_compact: {type(exc).__name__}: {exc}")
    sys.exit(0)
