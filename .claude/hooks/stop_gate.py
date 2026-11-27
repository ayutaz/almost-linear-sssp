#!/usr/bin/env python3
"""Stop hook: keep Claude working while the job is visibly unfinished.

Issues (each switchable in .tdd/config.json ``stop_gate``):
1. check_failed_tests     last_check.json has a category with status "fail"
2. check_red_tasks        an active task whose last TDD event is RED
3. check_uncommitted      ``git status --porcelain`` is non-empty
4. check_progress_updated docs/PROGRESS.md missing or older than a changed file

When there are issues the hook answers ``{"decision": "block", "reason": ...}``.
Loop guards: an explicit stop request in the user's last prompt, the
``stop_hook_active`` flag, and ``max_consecutive_blocks`` identical blocks in a
row all let Claude stop; the last two also append a report to
.tdd/runtime/stop_gate_unresolved.md (git-ignored: it holds session ids).

Fail-open: internal errors never block (exit 0, message on stderr).
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import _common  # noqa: E402

REPORT_REL = ".tdd/runtime/stop_gate_unresolved.md"
PROGRESS_REL = "docs/PROGRESS.md"
STATE_NAME = "stop_gate.json"
RUNTIME_PREFIX = _common.RUNTIME_REL.as_posix() + "/"
STOP_RE = re.compile(
    r"止めて|止まって|やめて|中断|ストップ|終了して|\bstop\b|\bhalt\b|\babort\b", re.IGNORECASE
)
STOP_PROMPT_MAX = 200
TRANSCRIPT_TAIL_BYTES = 2_000_000
TDD = "uv run python scripts/tdd.py"


# ------------------------------------------------------------------ git helpers
def _porcelain(root: Path) -> list[tuple[str, str]] | None:
    """[(status XY, path)] of uncommitted changes, or None when git fails."""
    code, out = _common.git(root, "status", "--porcelain", "--untracked-files=all")
    if code != 0:
        return None
    rows = []
    for line in out.splitlines():
        if len(line) < 4:
            continue
        xy, path = line[:2], line[3:]
        if " -> " in path:  # rename: keep the new name
            path = path.split(" -> ", 1)[1]
        if len(path) >= 2 and path.startswith('"') and path.endswith('"'):
            path = path[1:-1]
        # our own report / runtime state must not cause (or re-hash) a block
        if path == REPORT_REL or path.startswith(RUNTIME_PREFIX):
            continue
        rows.append((xy, path))
    return rows


# ---------------------------------------------------------------------- checks
def collect_issues(root: Path, gate: dict[str, Any]) -> list[dict[str, str]]:
    """Each issue = {"kind", "key" (stable, hashed), "text" (human readable)}."""
    issues: list[dict[str, str]] = []

    if gate.get("check_failed_tests", True):
        obj = _common.load_json(root / _common.RUNTIME_REL / "last_check.json", {})
        cats = obj.get("categories") if isinstance(obj, dict) else None
        if isinstance(cats, dict):
            for name in sorted(cats):
                info = cats[name]
                if isinstance(info, dict) and info.get("status") == "fail":
                    cmd = str(info.get("command") or "")
                    issues.append(
                        {
                            "kind": "failed_tests",
                            "key": f"{name}|{cmd}",
                            "text": f"検査 `{name}` が失敗中: `{cmd}`",
                        }
                    )

    if gate.get("check_red_tasks", True):
        tasks = _common.load_tasks(root)
        for tid in sorted(_common.load_active(root)):
            task = tasks.get(tid)
            if task and str(_common.last_event_phase(task) or "").upper() == "RED":
                issues.append(
                    {
                        "kind": "red_task",
                        "key": str(tid),
                        "text": f"タスク {tid} が RED のまま放置されている",
                    }
                )

    rows = None
    if gate.get("check_uncommitted", True) or gate.get("check_progress_updated", True):
        rows = _porcelain(root)

    if gate.get("check_uncommitted", True) and rows:
        head = ", ".join(p for _, p in rows[:5])
        more = f" ほか {len(rows) - 5} 件" if len(rows) > 5 else ""
        issues.append(
            {
                "kind": "uncommitted",
                # Stable across turns: "same reason" = uncommitted work, whatever the files.
                "key": "dirty",
                "text": f"未コミットの変更が {len(rows)} 件: {head}{more}",
            }
        )

    if gate.get("check_progress_updated", True):
        progress = root / PROGRESS_REL
        if not progress.is_file():
            issues.append({"kind": "progress", "key": "missing", "text": f"{PROGRESS_REL} が無い"})
        elif rows:
            try:
                pm = progress.stat().st_mtime
            except OSError:
                pm = 0.0
            newer = []
            for xy, path in rows:
                if "D" in xy or path == PROGRESS_REL:
                    continue
                try:
                    if (root / path).stat().st_mtime > pm:
                        newer.append(path)
                except OSError:
                    continue
            if newer:
                issues.append(
                    {
                        "kind": "progress",
                        "key": "stale",
                        "text": f"{PROGRESS_REL} が変更より古い(例: {', '.join(newer[:3])})",
                    }
                )
    return issues


def issues_hash(issues: list[dict[str, str]]) -> str:
    raw = json.dumps(sorted((i["kind"], i["key"]) for i in issues), ensure_ascii=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]


# ------------------------------------------------------------------ transcript
def _prompt_text(row: Any) -> str | None:
    """Text of a human prompt row, or None for tool results / system rows."""
    if not isinstance(row, dict) or row.get("type") != "user":
        return None
    msg = row.get("message")
    content = msg.get("content") if isinstance(msg, dict) else None
    if isinstance(content, str):
        text = content
    elif isinstance(content, list):
        parts = [
            str(c.get("text") or "")
            for c in content
            if isinstance(c, dict) and c.get("type") == "text"
        ]
        if not parts:
            return None
        text = "\n".join(parts)
    else:
        return None
    text = text.strip()
    if not text or text.startswith("<"):
        return None
    return text


def last_human_prompt(transcript_path: Any) -> str | None:
    if not isinstance(transcript_path, str) or not transcript_path:
        return None
    try:
        with open(transcript_path, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            fh.seek(max(0, size - TRANSCRIPT_TAIL_BYTES))
            raw = fh.read()
    except OSError:
        return None
    for line in reversed(raw.decode("utf-8", errors="replace").splitlines()):
        line = line.strip()
        if not line:
            continue
        try:
            row = json.loads(line)
        except ValueError:
            continue
        text = _prompt_text(row)
        if text is not None:
            return text
    return None


def user_requested_stop(data: dict[str, Any]) -> bool:
    text = last_human_prompt(data.get("transcript_path"))
    return bool(text) and len(text) <= STOP_PROMPT_MAX and bool(STOP_RE.search(text))


# ----------------------------------------------------------------- state/report
def load_state(root: Path) -> dict[str, Any]:
    obj = _common.load_json(root / _common.RUNTIME_REL / STATE_NAME, {})
    return obj if isinstance(obj, dict) else {}


def save_state(root: Path, state: dict[str, Any]) -> None:
    state["updated"] = _common.now_iso()
    try:
        _common.save_json(_common.runtime_dir(root) / STATE_NAME, state)
    except Exception as exc:
        _common.debug(f"stop_gate: state not saved: {exc}")


def append_report(
    root: Path, data: dict[str, Any], issues: list[dict[str, str]], cause: str
) -> None:
    path = root / REPORT_REL
    path.parent.mkdir(parents=True, exist_ok=True)
    sample = {
        "hook_event_name": "Stop",
        "session_id": data.get("session_id") or "",
        "cwd": ".",  # run from the repository root; never record absolute paths
        "stop_hook_active": False,
        "transcript_path": "",
    }
    payload = json.dumps(sample, ensure_ascii=False).replace("'", "'\\''")
    lines = [
        f"## {_common.now_iso()} session={data.get('session_id') or '(不明)'}",
        "",
        f"- 発生条件: {cause}",
        "- 未解決の理由:",
        *(f"  - {i['text']}" for i in issues),
        "- 再現手順:",
        "",
        "```sh",
        f"echo '{payload}' | python3 .claude/hooks/stop_gate.py",
        "git status --porcelain",
        "cat .tdd/runtime/last_check.json",
        f"{TDD} status",
        "```",
        "",
    ]
    new_file = not path.exists()
    with open(path, "a", encoding="utf-8") as fh:
        if new_file:
            fh.write("# Stop gate: 未解決のまま停止を許可した記録\n\n")
        fh.write("\n".join(lines) + "\n")


def block_reason(issues: list[dict[str, str]]) -> str:
    kinds = {i["kind"] for i in issues}
    out = ["作業がまだ終わっていないため停止できません。次を解消してください:"]
    out += [f"- {i['text']}" for i in issues]
    out.append("解消方法:")
    if "failed_tests" in kinds:
        out.append("- 失敗している検査を再実行し、テストを通るよう実装を修正する")
    if "red_task" in kinds:
        out.append(f"- 実装してテストを通し、`{TDD} green <ID> -- ctest --preset debug` を記録する")
    if "uncommitted" in kinds:
        out.append(
            "- 変更をコミットする(許可が無ければ利用者にコミットを依頼するか、"
            "サブエージェントならメインに変更一覧を報告する)"
        )
    if "progress" in kinds:
        out.append(f"- {PROGRESS_REL} に今回の進捗と次の作業を書く")
    out.append("どうしても解消できない場合は、理由を明記して利用者に報告してください。")
    return "\n".join(out)


# ------------------------------------------------------------------------ main
# Background work that may still be editing the main working tree.
HELPER_TYPES = ("subagent", "teammate", "workflow")
FINISHED_STATUSES = (
    "completed", "complete", "done", "failed", "error", "killed", "stopped", "cancelled",
    "canceled",
)  # fmt: skip
# Problems a running helper can cause and the lead cannot resolve alone.
SOFT_KINDS = ("uncommitted", "progress")


def running_helpers(data: dict[str, Any]) -> list[str]:
    """Labels of in-flight subagents/teammates/workflows from the Stop input."""
    out = []
    tasks = data.get("background_tasks")
    for t in tasks if isinstance(tasks, list) else []:
        if not isinstance(t, dict) or t.get("type") not in HELPER_TYPES:
            continue
        if str(t.get("status") or "").strip().lower() in FINISHED_STATUSES:
            continue
        label = t.get("agent_type") or t.get("name") or t.get("description") or t.get("id")
        out.append(f"{t.get('type')}: {str(label or '?')[:60]}")
    return out


def decide(data: dict[str, Any]) -> dict[str, Any] | None:
    if os.environ.get("TDD_STOP_GATE", "").strip().lower() == "off":
        return None
    root = _common.project_root(data)
    gate = _common.load_config(root).get("stop_gate")
    gate = gate if isinstance(gate, dict) else {}
    if not gate.get("enabled", True):
        return None

    issues = collect_issues(root, gate)
    warnings: list[dict[str, str]] = []
    helpers = running_helpers(data)
    if helpers:
        # Uncommitted work may belong to helpers that are still running: warn, do not block.
        warnings = [i for i in issues if i["kind"] in SOFT_KINDS]
        issues = [i for i in issues if i["kind"] not in SOFT_KINDS]
    state = load_state(root)
    if not issues:
        state.update({"last_reason_hash": None, "consecutive": 0, "last_reported_hash": None})
        save_state(root, state)
        if warnings:
            lines = [
                f"Stop gate: 実行中の担当が {len(helpers)} 件あるため、次は警告のみ"
                f"({', '.join(helpers[:3])}):",
                *(f"- {w['text']}" for w in warnings),
            ]
            return {"systemMessage": "\n".join(lines)}
        return None

    if user_requested_stop(data):
        _common.debug("stop_gate: user asked to stop; allowing")
        return None

    h = issues_hash(issues)
    try:
        max_blocks = int(gate.get("max_consecutive_blocks", 3))
    except (TypeError, ValueError):
        max_blocks = 3
    prev = state.get("consecutive") if state.get("last_reason_hash") == h else 0
    prev = prev if isinstance(prev, int) and prev >= 0 else 0

    cause = None
    if data.get("stop_hook_active") is True:
        cause = "stop_hook_active=true(Stop hook による継続中に再度停止しようとした)"
    elif prev + 1 >= max_blocks:
        cause = f"同一理由での block が連続 {prev + 1} 回(上限 {max_blocks})に達した"

    if cause is not None:
        if state.get("last_reported_hash") != h:
            try:
                append_report(root, data, issues, cause)
                state["last_reported_hash"] = h
            except Exception as exc:
                _common.debug(f"stop_gate: report not written: {exc}")
        # Only real blocks are counted; an allowed stop leaves the counter as it was.
        state.update({"last_reason_hash": h, "consecutive": prev})
        save_state(root, state)
        return {
            "systemMessage": (
                f"Stop gate: 未解決の問題 {len(issues)} 件を残して停止を許可しました"
                f"(詳細は {REPORT_REL})。"
            )
        }

    state.update({"last_reason_hash": h, "consecutive": prev + 1})
    save_state(root, state)
    return {"decision": "block", "reason": block_reason(issues)}


def main() -> int:
    data = _common.read_input()
    result = decide(data)
    if result:
        _common.emit(result)
    return 0


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # fail-open: never trap Claude because of a bug here
        _common.debug(f"stop_gate: {type(exc).__name__}: {exc}")
    sys.exit(0)
