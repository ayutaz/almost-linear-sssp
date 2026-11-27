#!/usr/bin/env python3
"""Wait for the dispatched `ci` and `sanitizer` runs of one or more commits, cheaply.

usage: scripts/ci/wait.py <sha>... [--since ISO] [--interval S] [--timeout-min M]

Runs started by scripts/ci/offload.sh carry the tested SHA in their run title
(`run-name: "<workflow> <sha> os=<os>"`). Instead of one `gh run watch` per run (which
exhausted the 5000 requests/hour API limit), this polls with ONE `gh run list` per workflow
every `--interval` seconds (default 90) and decides every SHA from that shared answer.

Only runs created at or after `--since` count, so an older run of the same SHA (for example a
cancelled one) is never taken for the new dispatch. Default: 10 minutes before start.
offload.sh prints the value to pass.

Repository: $SSSP_REPO (owner/name) if set, otherwise `gh repo view`.
Exit code: 0 iff every workflow of every SHA concluded `success`; 1 otherwise;
2 on timeout (`--timeout-min`, counted in poll intervals) or when the repository is unknown.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
import time
from collections.abc import Callable

WORKFLOWS = ("ci", "sanitizer")
REPO_ENV = "SSSP_REPO"
LIST_LIMIT = "60"

Runner = Callable[..., subprocess.CompletedProcess]


def parse_time(text: str) -> dt.datetime:
    t = dt.datetime.fromisoformat(text.replace("Z", "+00:00"))
    return t if t.tzinfo else t.replace(tzinfo=dt.UTC)


def pick_run(runs: list[dict], sha: str, since: dt.datetime) -> dict | None:
    """Newest run whose title contains `sha` and that was created at or after `since`."""
    mine = [
        r for r in runs if sha in r.get("displayTitle", "") and parse_time(r["createdAt"]) >= since
    ]
    return max(mine, key=lambda r: parse_time(r["createdAt"]), default=None)


def resolve_repo(run: Runner) -> str | None:
    if os.environ.get(REPO_ENV):
        return os.environ[REPO_ENV]
    p = run(
        ["gh", "repo", "view", "--json", "nameWithOwner", "-q", ".nameWithOwner"],
        capture_output=True,
        text=True,
        check=False,
    )
    if p.returncode != 0:
        return None
    return p.stdout.strip() or None


def list_runs(run: Runner, repo: str, workflow: str) -> list[dict] | None:
    p = run(
        # -b main: offload.sh dispatches with `--ref main` (main's workflow definition)
        ["gh", "run", "list", "-R", repo, "-w", workflow, "-b", "main"]
        + ["-e", "workflow_dispatch", "-L", LIST_LIMIT]
        + ["--json", "databaseId,displayTitle,status,conclusion,createdAt"],
        capture_output=True,
        text=True,
        check=False,
    )
    if p.returncode != 0:
        log(f"gh error ({workflow}): {p.stderr.strip()[:120]}")
        return None
    return json.loads(p.stdout)


def log(msg: str) -> None:
    print(time.strftime("%H:%M"), msg, flush=True)


def parse_args(argv: list[str] | None, now: dt.datetime) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("shas", nargs="+", metavar="sha")
    ap.add_argument("--since", type=parse_time, default=now - dt.timedelta(minutes=10))
    ap.add_argument("--interval", type=float, default=90.0)
    ap.add_argument("--timeout-min", type=float, default=240.0)
    return ap.parse_args(argv)


def main(
    argv: list[str] | None = None,
    *,
    run: Runner = subprocess.run,
    sleep: Callable[[float], object] = time.sleep,
) -> int:
    args = parse_args(argv, dt.datetime.now(dt.UTC))
    repo = resolve_repo(run)
    if not repo:
        print(f"cannot determine the repository: set {REPO_ENV}=owner/name", file=sys.stderr)
        return 2
    done: dict[tuple[str, str], str] = {}
    waited = 0.0
    while True:
        for wf in WORKFLOWS:
            if all((wf, s) in done for s in args.shas):
                continue
            runs = list_runs(run, repo, wf)
            if runs is None:
                continue
            for s in args.shas:
                if (wf, s) in done:
                    continue
                r = pick_run(runs, s, args.since)
                if r and r["status"] == "completed":
                    done[(wf, s)] = r["conclusion"]
                    log(f"{s[:7]} {wf} run {r['databaseId']}: {r['conclusion']}")
        if len(done) == len(WORKFLOWS) * len(args.shas):
            break
        if waited + args.interval > args.timeout_min * 60:
            log(f"timeout after {waited:.0f} s; still waiting for {' '.join(_pending(args, done))}")
            return 2
        sleep(args.interval)
        waited += args.interval
    ok = True
    for s in args.shas:
        green = all(done[(wf, s)] == "success" for wf in WORKFLOWS)
        ok &= green
        print(s[:7], "GREEN" if green else "RED", flush=True)
    return 0 if ok else 1


def _pending(args: argparse.Namespace, done: dict[tuple[str, str], str]) -> list[str]:
    return [f"{s[:7]}/{wf}" for s in args.shas for wf in WORKFLOWS if (wf, s) not in done]


if __name__ == "__main__":
    sys.exit(main())
