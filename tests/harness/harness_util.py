"""Helpers shared by the harness tests (import as ``from harness_util import ...``).

Kept out of conftest.py on purpose: tests/python has its own conftest.py, and
``from conftest import`` would resolve to whichever one is first on sys.path.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
HOOKS_DIR = REPO / ".claude" / "hooks"
TDD_PY = REPO / "scripts" / "tdd.py"


@dataclass
class HookResult:
    code: int
    stdout: str
    stderr: str

    @property
    def json(self) -> dict[str, Any] | None:
        text = self.stdout.strip()
        if not text:
            return None
        return json.loads(text)

    @property
    def decision(self) -> str | None:
        """PreToolUse permissionDecision, or None when the hook stayed silent."""
        data = self.json or {}
        return (data.get("hookSpecificOutput") or {}).get("permissionDecision")

    @property
    def context(self) -> str:
        data = self.json or {}
        return (data.get("hookSpecificOutput") or {}).get("additionalContext") or ""


def write_task(root: Path, task: dict[str, Any], active: bool = True, **active_fields: Any) -> None:
    """Create docs/tdd/tasks/<id>.json and (optionally) register it as active."""
    task = {
        "spec": "s",
        "paper_ref": None,
        "files": {"impl": [], "tests": []},
        "phase": "RED",
        "events": [],
        "commit": None,
        **task,
    }
    (root / "docs/tdd/tasks" / f"{task['id']}.json").write_text(json.dumps(task))
    if active:
        rt = root / ".tdd" / "runtime"
        rt.mkdir(parents=True, exist_ok=True)
        path = rt / "active.json"
        data = json.loads(path.read_text()) if path.exists() else {"tasks": {}}
        data["tasks"][task["id"]] = {"started": "2026-01-01T00:00:00+09:00", **active_fields}
        path.write_text(json.dumps(data))


RED_EVENT = {
    "phase": "RED",
    "tests": ["t"],
    "command": "false",
    "exit_code": 1,
    "reason": "feature absent",
    "timestamp": "2026-01-01T00:00:00+09:00",
    "output_tail": "",
}
GREEN_EVENT = {
    "phase": "GREEN",
    "tests": ["t"],
    "command": "true",
    "exit_code": 0,
    "timestamp": "2026-01-01T00:01:00+09:00",
    "output_tail": "",
}
