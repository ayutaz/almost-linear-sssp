#!/usr/bin/env python3
"""TDD state and evidence manager.

Evidence lives in ``docs/tdd/tasks/<ID>.json`` (format: docs/tdd/README.md, one
file per task so parallel work never conflicts). Work-in-progress state lives in
``.tdd/runtime/`` (gitignored): ``active.json`` lists the tasks being worked on,
``runs/`` keeps the full output of every recorded command and ``last_run.json``
the most recent one. The Claude Code hooks in ``.claude/hooks`` read both.

Usage (standard library only; ``uv run python scripts/tdd.py ...``)::

    tdd.py start <ID> --spec TEXT [--paper REF] --impl GLOB.. --test GLOB..
    tdd.py start <ID> --spec TEXT --test GLOB.. [--covers GLOB..] --test-only REASON
    tdd.py red <ID> --reason TEXT [--tests NAME..] -- <cmd...>   # must fail
    tdd.py green <ID> [--tests NAME..] -- <cmd...>               # must pass
    tdd.py refactor <ID> [--tests NAME..] -- <cmd...>            # must pass
    tdd.py allow-test-change <ID> --reason TEXT   # permit removing/merging tests
    tdd.py done <ID> [--commit SHA]
    tdd.py status | check [--no-git] | log

A command given as one argument (``-- "ctest --preset debug -R x"``) runs
through the shell; several arguments run without a shell.
"""

from __future__ import annotations

import argparse
import contextlib
import datetime as dt
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
from collections.abc import Iterator
from pathlib import Path
from typing import Any

ROOT = Path(os.environ.get("TDD_ROOT") or Path(__file__).resolve().parents[1]).resolve()
TASKS_DIR = ROOT / "docs" / "tdd" / "tasks"
RUNTIME = ROOT / ".tdd" / "runtime"
ACTIVE = RUNTIME / "active.json"
LOG_MD = ROOT / "docs" / "TDD_LOG.md"

PHASES = ("RED", "GREEN", "REFACTOR", "DONE")
EVENT_PHASES = ("RED", "GREEN", "REFACTOR")
REQUIRED_KEYS = ("id", "spec", "files", "phase", "events", "commit")
ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,99}$")
TAIL_LINES = 20
TEST_ONLY = "test-only"
DEFAULT_IMPL_DIRS = ("src/", "include/", "tools/")


class TddError(Exception):
    """User-facing error: printed without a traceback, exit code 1."""


# ----------------------------------------------------------------------- helpers
def now() -> str:
    return dt.datetime.now().astimezone().isoformat(timespec="seconds")


def parse_ts(value: Any) -> dt.datetime | None:
    if not isinstance(value, str) or not value:
        return None
    v = value.strip()
    if v.endswith("Z"):
        v = v[:-1] + "+00:00"
    v = re.sub(r"([+-]\d{2})(\d{2})$", r"\1:\2", v)
    try:
        return dt.datetime.fromisoformat(v)
    except ValueError:
        return None


def check_id(tid: str) -> str:
    if not ID_RE.match(tid or ""):
        raise TddError(f"invalid task id {tid!r} (letters, digits, '.', '_', '-' only)")
    return tid


def task_path(tid: str) -> Path:
    return TASKS_DIR / f"{check_id(tid)}.json"


def read_json(path: Path, default: Any = None) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return default


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.name}.tmp{os.getpid()}")
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def load_task(tid: str) -> dict[str, Any]:
    path = task_path(tid)
    task = read_json(path)
    if not isinstance(task, dict):
        raise TddError(f"task {tid} not found ({path.relative_to(ROOT)}); run `tdd.py start` first")
    return task


def save_task(task: dict[str, Any]) -> None:
    write_json(task_path(task["id"]), task)


@contextlib.contextmanager
def runtime_lock() -> Iterator[None]:
    """Serialise read-modify-write of active.json between parallel agents."""
    RUNTIME.mkdir(parents=True, exist_ok=True)
    with open(RUNTIME / "tdd.lock", "w") as fh:
        try:
            import fcntl

            fcntl.flock(fh, fcntl.LOCK_EX)
        except (ImportError, OSError):  # pragma: no cover - non-POSIX
            pass
        yield


def update_active(tid: str, remove: bool = False, **fields: Any) -> None:
    with runtime_lock():
        data = read_json(ACTIVE, {})
        if not isinstance(data, dict) or not isinstance(data.get("tasks"), dict):
            data = {"tasks": {}}
        tasks = data["tasks"]
        if remove:
            tasks.pop(tid, None)
        else:
            entry = tasks.get(tid) or {
                "started": now(),
                "allow_test_reduction": False,
                "allow_reason": None,
            }
            entry.update(fields)
            tasks[tid] = entry
        write_json(ACTIVE, data)


def active_tasks() -> dict[str, Any]:
    data = read_json(ACTIVE, {})
    tasks = data.get("tasks") if isinstance(data, dict) else None
    return tasks if isinstance(tasks, dict) else {}


def phases_of(task: dict[str, Any]) -> list[str]:
    return [e.get("phase") for e in task.get("events") or [] if isinstance(e, dict)]


def is_test_only(task: dict[str, Any]) -> bool:
    return task.get("kind") == TEST_ONLY


def impl_dirs() -> tuple[str, ...]:
    cfg = read_json(ROOT / ".tdd" / "config.json", {})
    dirs = cfg.get("impl_dirs") if isinstance(cfg, dict) else None
    if isinstance(dirs, list) and all(isinstance(d, str) for d in dirs) and dirs:
        return tuple(dirs)
    return DEFAULT_IMPL_DIRS


def _git(*args: str) -> tuple[int, str] | None:
    """Run git in the repository root; None when git cannot be started."""
    try:
        proc = subprocess.run(
            ["git", "-C", str(ROOT), *args], capture_output=True, text=True, timeout=20
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    return proc.returncode, proc.stdout.strip()


def commit_files(sha: str) -> list[str] | None:
    """Paths changed by ``sha`` (None when git or the commit is unavailable)."""
    res = _git("show", "--no-renames", "--name-only", "--format=", sha)
    if res is None or res[0] != 0:
        return None
    return [ln.strip() for ln in res[1].splitlines() if ln.strip()]


def tail(text: str, n: int = TAIL_LINES) -> str:
    return "\n".join(text.rstrip("\n").splitlines()[-n:])


def redact(text: str) -> str:
    """Strip machine-specific absolute paths before text goes into committed evidence.

    The repository root becomes ``.``, the temp dir ``$TMPDIR`` and the home ``~``.
    """
    roots = {str(ROOT), os.path.realpath(ROOT)}
    for r in sorted(roots, key=len, reverse=True):
        text = text.replace(r, ".")
    tmp = tempfile.gettempdir()
    for t in sorted({tmp, os.path.realpath(tmp)}, key=len, reverse=True):
        if len(t) > 1:
            text = text.replace(t.rstrip("/"), "$TMPDIR")
    with contextlib.suppress(Exception):
        import getpass

        user = getpass.getuser()
        if user:
            text = text.replace(f"pytest-of-{user}", "pytest-of-$USER")
    home = os.path.expanduser("~")
    if home and home != "/":
        for h in sorted({home, os.path.realpath(home)}, key=len, reverse=True):
            text = text.replace(h, "~")
    # Shape-based pass for what the exact replacements miss: other machines' paths and
    # pytest's truncated reprs (a temp prefix cut short and ended with '...'). Mirrors the
    # path rules of scripts/check_public_safety.py.
    for pat, repl in _REDACT_PATTERNS:
        text = pat.sub(repl, text)
    return text


_PATH_TAIL = r"[^\s'\"\\]*"  # up to whitespace, a quote or a backslash
_REDACT_PATTERNS = [
    (re.compile(r"(?:/private)?/var/folders/" + _PATH_TAIL), "$TMPDIR"),  # public-safety: allow
    (
        re.compile(r"/private/(?:tmp|var)(?:/" + _PATH_TAIL + r"|\.\.\." + _PATH_TAIL + ")"),
        "$TMPDIR",
    ),
    (re.compile(r"/private/(?:v|va|t|tm)?\.\.\." + _PATH_TAIL), "$TMPDIR"),
    (re.compile(r"/(?:Users|home)/[A-Za-z0-9._-]+/"), "~/"),
]


def run_command(cmd: list[str]) -> tuple[int, str, str]:
    """Run ``cmd`` in the repository root, streaming output. Returns (code, text, display)."""
    if not cmd:
        raise TddError("missing command after `--`")
    if len(cmd) == 1:
        display, popen_args, shell = cmd[0], cmd[0], True
    else:
        display, popen_args, shell = shlex.join(cmd), cmd, False
    print(f"$ {display}", flush=True)
    chunks: list[str] = []
    try:
        proc = subprocess.Popen(
            popen_args,
            shell=shell,
            cwd=str(ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
        )
    except OSError as exc:
        return 127, f"{exc}\n", display
    assert proc.stdout is not None
    for line in proc.stdout:
        sys.stdout.write(line)
        chunks.append(line)
    sys.stdout.flush()
    return proc.wait(), "".join(chunks), display


def save_run(tid: str, phase: str, code: int, display: str, output: str) -> Path:
    runs = RUNTIME / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    path = runs / f"{tid}-{stamp}-{phase}.log"
    n = 1
    while path.exists():
        n += 1
        path = runs / f"{tid}-{stamp}.{n}-{phase}.log"
    path.write_text(f"$ {display}\n# exit {code}\n{output}", encoding="utf-8")
    write_json(
        RUNTIME / "last_run.json",
        {
            "task": tid,
            "phase": phase,
            "command": display,
            "exit_code": code,
            "timestamp": now(),
            "log": str(path.relative_to(ROOT)),
            "output_tail": tail(output),
        },
    )
    return path


# ---------------------------------------------------------------------- commands
def cmd_start(args: argparse.Namespace, _cmd: list[str]) -> int:
    tid = check_id(args.id)
    path = task_path(tid)
    task = read_json(path)
    if isinstance(task, dict):
        if args.spec:
            task["spec"] = args.spec
        if args.paper is not None:
            task["paper_ref"] = args.paper
        files = task.setdefault("files", {"impl": [], "tests": []})
        for key, extra in (("impl", args.impl), ("tests", args.test), ("covers", args.covers)):
            for g in extra or []:
                if g not in files.setdefault(key, []):
                    files[key].append(g)
        if task.get("phase") == "DONE":
            task["phase"] = "GREEN" if is_test_only(task) else "RED"
        save_task(task)
        print(
            f"{tid}: re-activated (phase {task.get('phase')}, "
            f"{len(task.get('events') or [])} events)"
        )
    elif args.test_only is not None:
        if not args.test_only.strip():
            raise TddError("--test-only needs a non-empty reason")
        if args.impl:
            raise TddError("a test-only task changes no implementation: use --covers, not --impl")
        missing = [n for n, v in (("--spec", args.spec), ("--test", args.test)) if not v]
        if missing:
            raise TddError(f"new task {tid} needs {', '.join(missing)}")
        files: dict[str, list[str]] = {"impl": [], "tests": list(args.test)}
        if args.covers:
            files["covers"] = list(args.covers)
        task = {
            "id": tid,
            "kind": TEST_ONLY,
            "test_only_reason": args.test_only.strip(),
            "spec": args.spec,
            "paper_ref": args.paper,
            "files": files,
            "phase": "GREEN",
            "events": [],
            "commit": None,
        }
        save_task(task)
        print(f"{tid}: created {path.relative_to(ROOT)} (test-only)")
    else:
        missing = [
            n
            for n, v in (("--spec", args.spec), ("--impl", args.impl), ("--test", args.test))
            if not v
        ]
        if missing:
            raise TddError(f"new task {tid} needs {', '.join(missing)}")
        task = {
            "id": tid,
            "spec": args.spec,
            "paper_ref": args.paper,
            "files": {"impl": list(args.impl), "tests": list(args.test)},
            "phase": "RED",
            "events": [],
            "commit": None,
        }
        if args.covers:
            task["files"]["covers"] = list(args.covers)
        save_task(task)
        print(f"{tid}: created {path.relative_to(ROOT)}")
    update_active(tid)
    if is_test_only(task):
        print("next: write the test, then `tdd.py green <ID> -- <test cmd>`")
    elif not phases_of(task):
        print("next: write the failing test, then `tdd.py red <ID> --reason .. -- <test cmd>`")
    return 0


def _record(args: argparse.Namespace, cmd: list[str], phase: str) -> int:
    tid = check_id(args.id)
    task = load_task(tid)
    seen = phases_of(task)
    if phase == "RED" and not (args.reason or "").strip():
        raise TddError("RED needs --reason (why the test fails; it must be the intended reason)")
    if phase == "RED" and is_test_only(task):
        raise TddError(
            f"{tid} is a test-only task; a failing test there means a bug: open a normal task"
        )
    if phase == "GREEN" and "RED" not in seen and not is_test_only(task):
        raise TddError(f"{tid}: GREEN needs a recorded RED first (`tdd.py red {tid} ...`)")
    if phase == "REFACTOR" and "GREEN" not in seen:
        raise TddError(f"{tid}: REFACTOR needs a recorded GREEN first")
    code, output, display = run_command(cmd)
    save_run(tid, phase, code, display, output)
    if phase == "RED" and code == 0:
        raise TddError(
            f"{tid}: command exited 0 - a RED event needs a failing test. Not recorded. "
            "Write a test that fails for the intended reason first."
        )
    if phase != "RED" and code != 0:
        raise TddError(f"{tid}: command exited {code} - {phase} needs exit 0. Not recorded.")
    event: dict[str, Any] = {
        "phase": phase,
        "tests": list(args.tests or []),
        "command": redact(display),
        "exit_code": code,
    }
    if phase == "RED":
        event["reason"] = args.reason.strip()
    event["timestamp"] = now()
    event["output_tail"] = redact(tail(output))
    task.setdefault("events", []).append(event)
    task["phase"] = phase
    save_task(task)
    fields: dict[str, Any] = {"phase": phase}
    if phase == "RED":
        fields.update(allow_test_reduction=False, allow_reason=None)
    update_active(tid, **fields)
    print(f"{tid}: recorded {phase} (exit {code})")
    return 0


def cmd_red(args: argparse.Namespace, cmd: list[str]) -> int:
    return _record(args, cmd, "RED")


def cmd_green(args: argparse.Namespace, cmd: list[str]) -> int:
    return _record(args, cmd, "GREEN")


def cmd_refactor(args: argparse.Namespace, cmd: list[str]) -> int:
    return _record(args, cmd, "REFACTOR")


def cmd_allow_test_change(args: argparse.Namespace, _cmd: list[str]) -> int:
    tid = check_id(args.id)
    task = load_task(tid)
    seen = phases_of(task)
    if not (args.reason or "").strip():
        raise TddError("allow-test-change needs --reason")
    if "GREEN" not in seen or seen[-1] == "RED":
        raise TddError(f"{tid}: test reduction is only allowed after GREEN (refactor phase)")
    task.setdefault("test_reduction_allowed", []).append(
        {"reason": args.reason.strip(), "timestamp": now()}
    )
    save_task(task)
    update_active(tid, allow_test_reduction=True, allow_reason=args.reason.strip())
    print(f"{tid}: removing/merging tests matching {task['files'].get('tests')} is allowed")
    return 0


def cmd_done(args: argparse.Namespace, _cmd: list[str]) -> int:
    tid = check_id(args.id)
    task = load_task(tid)
    seen = phases_of(task)
    if "GREEN" not in seen or seen[-1] == "RED":
        raise TddError(f"{tid}: DONE needs a GREEN (or REFACTOR) as the last event")
    task["phase"] = "DONE"
    if args.commit:
        task["commit"] = args.commit
    save_task(task)
    update_active(tid, remove=True)
    print(f"{tid}: DONE")
    return 0


def all_tasks() -> list[tuple[Path, Any]]:
    if not TASKS_DIR.is_dir():
        return []
    return [(p, read_json(p)) for p in sorted(TASKS_DIR.glob("*.json"))]


def cmd_status(_args: argparse.Namespace, _cmd: list[str]) -> int:
    act = active_tasks()
    tasks = {p.stem: t for p, t in all_tasks() if isinstance(t, dict)}
    print(f"active tasks ({len(act)}):")
    for tid in sorted(act):
        t = tasks.get(tid, {})
        seen = phases_of(t)
        last = seen[-1] if seen else "-"
        impl = ", ".join((t.get("files") or {}).get("impl") or [])
        flag = " [test reduction allowed]" if act[tid].get("allow_test_reduction") else ""
        print(
            f"  {tid}: phase={t.get('phase', '?')} last={last} events={len(seen)} impl={impl}{flag}"
        )
    if not act:
        print("  (none)")
    counts: dict[str, int] = {}
    for t in tasks.values():
        counts[str(t.get("phase"))] = counts.get(str(t.get("phase")), 0) + 1
    summary = ", ".join(f"{k}={v}" for k, v in sorted(counts.items())) or "none"
    print(f"all tasks ({len(tasks)}): {summary}")
    last = read_json(RUNTIME / "last_check.json", {})
    cats = last.get("categories") if isinstance(last, dict) else None
    if isinstance(cats, dict) and cats:
        state = ", ".join(f"{c}={(v or {}).get('status')}" for c, v in sorted(cats.items()))
        print(f"last hook check: {state}")
    return 0


def validate(path: Path, task: Any) -> list[str]:
    """Return a list of problems for one evidence file (empty = valid)."""
    if not isinstance(task, dict):
        return ["not a JSON object (or invalid JSON)"]
    errs = [f"missing key '{k}'" for k in REQUIRED_KEYS if k not in task]
    if errs:
        return errs
    if task["id"] != path.stem:
        errs.append(f"id '{task['id']}' does not match file name '{path.stem}'")
    if not isinstance(task["spec"], str) or not task["spec"].strip():
        errs.append("spec must be a non-empty string")
    files = task["files"]
    if not isinstance(files, dict) or not all(
        isinstance(files.get(k), list) and all(isinstance(g, str) for g in files.get(k))
        for k in ("impl", "tests")
    ):
        errs.append("files must be {impl: [str], tests: [str]}")
    elif "covers" in files and not (
        isinstance(files["covers"], list) and all(isinstance(g, str) for g in files["covers"])
    ):
        errs.append("files.covers must be a list of strings")
    test_only = False
    if "kind" in task:
        if task["kind"] != TEST_ONLY:
            errs.append(f"kind {task['kind']!r} is not supported (only '{TEST_ONLY}')")
        else:
            test_only = True
            errs += _validate_test_only(task)
    phase = task["phase"]
    if phase not in PHASES:
        errs.append(f"phase '{phase}' is not one of {'|'.join(PHASES)}")
    events = task["events"]
    if not isinstance(events, list):
        return [*errs, "events must be a list"]
    seen: list[str] = []
    for i, ev in enumerate(events):
        where = f"events[{i}]"
        if not isinstance(ev, dict):
            errs.append(f"{where} is not an object")
            continue
        ph = ev.get("phase")
        if ph not in EVENT_PHASES:
            errs.append(f"{where}: phase '{ph}' is not one of {'|'.join(EVENT_PHASES)}")
            continue
        code = ev.get("exit_code")
        if not isinstance(code, int) or isinstance(code, bool):
            errs.append(f"{where}: exit_code must be an integer")
        elif ph == "RED" and code == 0:
            errs.append(f"{where}: RED exit_code must be non-zero")
        elif ph != "RED" and code != 0:
            errs.append(f"{where}: {ph} exit_code must be 0 (got {code})")
        if ph == "RED" and not (isinstance(ev.get("reason"), str) and ev["reason"].strip()):
            errs.append(f"{where}: RED needs a non-empty reason")
        if not isinstance(ev.get("command"), str) or not ev["command"].strip():
            errs.append(f"{where}: command must be a non-empty string")
        if parse_ts(ev.get("timestamp")) is None:
            errs.append(f"{where}: timestamp {ev.get('timestamp')!r} is not ISO 8601")
        if ph == "GREEN" and "RED" not in seen and not test_only:
            errs.append(f"{where}: GREEN without an earlier RED")
        if ph == "REFACTOR" and "GREEN" not in seen:
            errs.append(f"{where}: REFACTOR without an earlier GREEN")
        seen.append(ph)
    if seen and seen[0] != "RED" and not test_only:
        errs.append("the first event must be RED")
    if phase == "DONE" and ("GREEN" not in seen or seen[-1] == "RED"):
        errs.append("phase DONE needs a GREEN/REFACTOR as the last event")
    elif phase in EVENT_PHASES and seen and seen[-1] != phase:
        errs.append(f"phase '{phase}' does not match the last event '{seen[-1]}'")
    return errs


def _validate_test_only(task: dict[str, Any]) -> list[str]:
    """Rules for ``"kind": "test-only"`` (verification tests that add no behaviour)."""
    errs = []
    reason = task.get("test_only_reason")
    if not (isinstance(reason, str) and reason.strip()):
        errs.append("test-only needs a non-empty test_only_reason")
    files = task.get("files")
    impl = files.get("impl") if isinstance(files, dict) else None
    if not impl:
        return errs
    sha = task.get("commit")
    if not (isinstance(sha, str) and sha.strip()):
        errs.append(
            "test-only task lists files.impl but has no commit to prove they were not changed: "
            "move the tested files to files.covers"
        )
        return errs
    changed = commit_files(sha.strip())
    if changed is None:
        errs.append(f"commit {sha} not found (needed to check a test-only task with files.impl)")
        return errs
    touched = [f for f in changed if f.startswith(impl_dirs())]
    if touched:
        errs.append(f"test-only commit {sha} changes implementation files: {', '.join(touched)}")
    return errs


def commit_check_unavailable() -> str | None:
    """Why commit reachability cannot be judged here (None = it can)."""
    res = _git("rev-parse", "--is-inside-work-tree")
    if res is None:
        return "git not found"
    if res[0] != 0 or res[1] != "true":
        return "not a git repository"
    shallow = _git("rev-parse", "--is-shallow-repository")
    if shallow is None or shallow[0] != 0 or shallow[1] != "false":
        return "shallow clone: history is cut off (use `git fetch --unshallow` or fetch-depth: 0)"
    head = _git("rev-parse", "--verify", "--quiet", "HEAD")
    if head is None or head[0] != 0:
        return "the repository has no HEAD commit"
    return None


def commit_problem(sha: Any) -> str | None:
    """Problem with a non-empty ``commit`` field, judged against HEAD (None = fine)."""
    if not isinstance(sha, str):
        return "commit must be a commit SHA string or null"
    sha = sha.strip()
    found = _git("rev-parse", "--verify", "--quiet", f"{sha}^{{commit}}")
    if found is None or found[0] != 0:
        return f"commit {sha} does not exist in this repository"
    res = _git("merge-base", "--is-ancestor", sha, "HEAD")
    if res is None:
        return None
    if res[0] == 1:
        return (
            f"commit {sha} is not reachable from HEAD "
            "(history rewritten, or the branch was never merged)"
        )
    if res[0] != 0:
        return f"cannot check commit {sha} against HEAD (git exit {res[0]})"
    return None


def cmd_check(args: argparse.Namespace, _cmd: list[str]) -> int:
    entries = all_tasks()
    commits = [(p, t) for p, t in entries if isinstance(t, dict) and t.get("commit")]
    skip_git = getattr(args, "no_git", False)
    unavailable = None if skip_git else (commit_check_unavailable() if commits else None)
    if skip_git:
        print("NOTE commit fields not checked against git history (--no-git)")
    elif unavailable:
        print(f"WARN commit fields not checked against git history: {unavailable}")
    bad = 0
    for path, task in entries:
        problems = validate(path, task)
        if not skip_git and not unavailable and isinstance(task, dict) and task.get("commit"):
            problem = commit_problem(task["commit"])
            if problem:
                problems.append(problem)
        for msg in problems:
            print(f"ERROR {path.stem}: {msg}")
        bad += bool(problems)
    print(f"checked {len(entries)} task file(s): {bad} invalid")
    return 1 if bad else 0


def _cell(value: Any) -> str:
    text = "" if value is None else str(value)
    return text.replace("\n", " ").replace("|", "\\|").strip()


def render_log(entries: list[tuple[Path, Any]]) -> str:
    out = [
        "# TDD log",
        "",
        "Generated by `uv run python scripts/tdd.py log` from `docs/tdd/tasks/*.json`.",
        "Do not edit by hand; edit the evidence files and regenerate.",
        "",
        "| Task | Phase | RED→GREEN | Paper | Commit | Spec |",
        "| :- | :- | -: | :- | :- | :- |",
    ]
    valid = [(p, t) for p, t in entries if isinstance(t, dict)]
    for path, t in valid:
        ph = phases_of(t)
        cycles = sum(1 for i, p in enumerate(ph) if p == "GREEN" and "RED" in ph[:i])
        out.append(
            f"| {_cell(t.get('id', path.stem))} | {_cell(t.get('phase'))}"
            f"{' (test-only)' if is_test_only(t) else ''} | {cycles} | "
            f"{_cell(t.get('paper_ref')) or '-'} | {_cell(t.get('commit')) or '-'} | "
            f"{_cell(t.get('spec'))} |"
        )
    for path, t in valid:
        out += ["", f"## {t.get('id', path.stem)}", "", f"- Spec: {_cell(t.get('spec'))}"]
        if is_test_only(t):
            out.append(f"- Test-only: {_cell(t.get('test_only_reason'))}")
        files = t.get("files") or {}
        out.append(f"- Impl: {', '.join(f'`{g}`' for g in files.get('impl') or []) or '-'}")
        out.append(f"- Tests: {', '.join(f'`{g}`' for g in files.get('tests') or []) or '-'}")
        if files.get("covers"):
            out.append(f"- Covers: {', '.join(f'`{g}`' for g in files['covers'])}")
        for ev in t.get("events") or []:
            if not isinstance(ev, dict):
                continue
            line = f"- {ev.get('phase')} (exit {ev.get('exit_code')}, {ev.get('timestamp')})"
            if ev.get("reason"):
                line += f": {_cell(ev.get('reason'))}"
            line += f" — `{_cell(ev.get('command'))}`"
            out.append(line)
    broken = [p.stem for p, t in entries if not isinstance(t, dict)]
    if broken:
        out += ["", f"Unreadable evidence files: {', '.join(broken)}"]
    return "\n".join(out) + "\n"


def cmd_log(_args: argparse.Namespace, _cmd: list[str]) -> int:
    LOG_MD.parent.mkdir(parents=True, exist_ok=True)
    LOG_MD.write_text(render_log(all_tasks()), encoding="utf-8")
    print(f"wrote {LOG_MD.relative_to(ROOT)}")
    return 0


# -------------------------------------------------------------------------- main
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="tdd.py", description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("start", help="create/re-activate a task")
    p.add_argument("id")
    p.add_argument("--spec")
    p.add_argument("--paper", default=None)
    p.add_argument("--impl", nargs="+")
    p.add_argument("--test", nargs="+")
    p.add_argument("--covers", nargs="+", help="files exercised but not changed")
    p.add_argument(
        "--test-only",
        metavar="REASON",
        default=None,
        help="verification-only task (no new behaviour): GREEN may come first",
    )
    p.set_defaults(func=cmd_start)

    for name, func, helptext in (
        ("red", cmd_red, "record a failing run (needs --reason)"),
        ("green", cmd_green, "record a passing run after RED"),
        ("refactor", cmd_refactor, "record a passing run after GREEN"),
    ):
        p = sub.add_parser(name, help=helptext)
        p.add_argument("id")
        p.add_argument("--reason")
        p.add_argument("--tests", nargs="*", default=[])
        p.set_defaults(func=func, needs_cmd=True)

    p = sub.add_parser("allow-test-change", help="permit test reduction in the refactor phase")
    p.add_argument("id")
    p.add_argument("--reason")
    p.set_defaults(func=cmd_allow_test_change)

    p = sub.add_parser("done", help="mark a task DONE and deactivate it")
    p.add_argument("id")
    p.add_argument("--commit", default=None)
    p.set_defaults(func=cmd_done)

    p = sub.add_parser("check", help="validate every evidence file")
    p.add_argument(
        "--no-git",
        action="store_true",
        help="do not check that commit fields are reachable from HEAD",
    )
    p.set_defaults(func=cmd_check)

    for name, func in (("status", cmd_status), ("log", cmd_log)):
        sub.add_parser(name).set_defaults(func=func)
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    cmd: list[str] = []
    if "--" in argv:
        i = argv.index("--")
        argv, cmd = argv[:i], argv[i + 1 :]
    args = build_parser().parse_args(argv)
    if getattr(args, "needs_cmd", False) and not cmd:
        print(f"error: `{args.command}` needs `-- <command...>`", file=sys.stderr)
        return 2
    try:
        return args.func(args, cmd)
    except TddError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
