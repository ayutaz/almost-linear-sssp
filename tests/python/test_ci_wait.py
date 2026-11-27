"""scripts/ci/wait.py with a FAKE `gh` (no network, no GitHub API calls).

The fake answers `gh run list` with canned JSON per poll, so the tests check what the tool
decides from that output: which run belongs to which SHA (run title), which runs are too old
(`--since`), how often it polls, and the exit code.
"""

from __future__ import annotations

import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

WAIT_PY = Path(__file__).resolve().parents[2] / "scripts" / "ci" / "wait.py"

SHA_A = "a" * 40
SHA_B = "b" * 40
SINCE = "2026-11-27T10:00:00Z"
BEFORE = "2026-11-27T09:59:00Z"
AFTER = "2026-11-27T10:05:00Z"
LATER = "2026-11-27T10:20:00Z"


def load_wait():
    spec = importlib.util.spec_from_file_location("ci_wait_under_test", WAIT_PY)
    assert spec is not None and spec.loader is not None, f"{WAIT_PY} is missing"
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run(rid, sha, created, status="completed", conclusion="success", wf="ci"):
    return {
        "databaseId": rid,
        "displayTitle": f"{wf} {sha} os=ubuntu",
        "status": status,
        "conclusion": conclusion if status == "completed" else "",
        "createdAt": created,
    }


class FakeGh:
    """Stands in for subprocess.run. polls[i][workflow] is the run list of the i-th poll."""

    def __init__(self, polls, repo="owner/repo"):
        self.polls = polls
        self.repo = repo
        self.calls: list[list[str]] = []
        self.poll = 0
        self.seen_in_poll: set[str] = set()

    def __call__(self, cmd, **_kw):
        self.calls.append(list(cmd))
        if cmd[:3] == ["gh", "repo", "view"]:
            return subprocess.CompletedProcess(cmd, 0, self.repo + "\n", "")
        assert cmd[:3] == ["gh", "run", "list"], cmd
        wf = cmd[cmd.index("-w") + 1]
        if wf in self.seen_in_poll:  # a workflow asked again: the next poll has started
            self.poll += 1
            self.seen_in_poll = set()
        self.seen_in_poll.add(wf)
        runs = self.polls[min(self.poll, len(self.polls) - 1)].get(wf, [])
        if runs == "ERROR":
            return subprocess.CompletedProcess(cmd, 1, "", "HTTP 502")
        return subprocess.CompletedProcess(cmd, 0, json.dumps(runs), "")

    def run_list_calls(self):
        return [c for c in self.calls if c[:3] == ["gh", "run", "list"]]


@pytest.fixture
def wait(monkeypatch):
    monkeypatch.setenv("SSSP_REPO", "owner/repo")
    return load_wait()


def go(wait, fake, *args):
    sleeps: list[float] = []
    rc = wait.main([*args], run=fake, sleep=sleeps.append)
    return rc, sleeps


# --- matching by run title -------------------------------------------------------------------


def test_run_is_matched_by_the_sha_in_its_title(wait):
    # SHA_B's runs failed; they must not be attributed to SHA_A.
    poll = {
        "ci": [run(2, SHA_B, AFTER, conclusion="failure"), run(1, SHA_A, AFTER)],
        "sanitizer": [run(4, SHA_B, AFTER, conclusion="failure"), run(3, SHA_A, AFTER)],
    }
    rc, _ = go(wait, FakeGh([poll]), SHA_A, "--since", SINCE)
    assert rc == 0


def test_short_sha_matches_the_full_sha_in_the_title(wait):
    poll = {"ci": [run(1, SHA_A, AFTER)], "sanitizer": [run(2, SHA_A, AFTER)]}
    rc, _ = go(wait, FakeGh([poll]), SHA_A[:7], "--since", SINCE)
    assert rc == 0


def test_waits_until_a_run_with_the_sha_appears(wait):
    other = {"ci": [run(9, SHA_B, AFTER)], "sanitizer": [run(8, SHA_B, AFTER)]}
    mine = {"ci": [run(1, SHA_A, LATER)], "sanitizer": [run(2, SHA_A, LATER)]}
    fake = FakeGh([other, mine])
    rc, sleeps = go(wait, fake, SHA_A, "--since", SINCE)
    assert rc == 0
    assert sleeps == [90]


def test_in_progress_run_is_waited_for(wait):
    running = {
        "ci": [run(1, SHA_A, AFTER, status="in_progress")],
        "sanitizer": [run(2, SHA_A, AFTER, status="queued")],
    }
    done = {"ci": [run(1, SHA_A, AFTER)], "sanitizer": [run(2, SHA_A, AFTER)]}
    rc, sleeps = go(wait, FakeGh([running, running, done]), SHA_A, "--since", SINCE)
    assert rc == 0
    assert sleeps == [90, 90]


# --- --since ---------------------------------------------------------------------------------


def test_run_created_before_since_is_ignored(wait):
    # An old cancelled run of the same SHA is listed first (newest-first order is not trusted);
    # only the new dispatch counts.
    poll = {
        "ci": [run(1, SHA_A, BEFORE, conclusion="cancelled"), run(5, SHA_A, AFTER)],
        "sanitizer": [run(2, SHA_A, BEFORE, conclusion="failure"), run(6, SHA_A, AFTER)],
    }
    rc, _ = go(wait, FakeGh([poll]), SHA_A, "--since", SINCE)
    assert rc == 0


def test_only_old_runs_means_keep_waiting(wait):
    old = {
        "ci": [run(1, SHA_A, BEFORE, conclusion="success")],
        "sanitizer": [run(2, SHA_A, BEFORE, conclusion="success")],
    }
    new = {
        "ci": [run(5, SHA_A, LATER, conclusion="failure"), *old["ci"]],
        "sanitizer": [run(6, SHA_A, LATER), *old["sanitizer"]],
    }
    rc, sleeps = go(wait, FakeGh([old, new]), SHA_A, "--since", SINCE)
    assert rc == 1, "the old green run must not be taken for the new (red) dispatch"
    assert sleeps == [90]


def test_since_accepts_an_utc_offset(wait):
    poll = {"ci": [run(1, SHA_A, AFTER)], "sanitizer": [run(2, SHA_A, AFTER)]}
    # 19:00+09:00 == 10:00Z, so AFTER (10:05Z) counts
    rc, _ = go(wait, FakeGh([poll]), SHA_A, "--since", "2026-11-27T19:00:00+09:00")
    assert rc == 0


def test_newest_matching_run_after_since_wins(wait):
    poll = {
        "ci": [run(1, SHA_A, AFTER, conclusion="failure"), run(7, SHA_A, LATER)],
        "sanitizer": [run(2, SHA_A, AFTER)],
    }
    rc, _ = go(wait, FakeGh([poll]), SHA_A, "--since", SINCE)
    assert rc == 0


# --- exit code -------------------------------------------------------------------------------


@pytest.mark.parametrize("bad", ["failure", "cancelled", "timed_out", "skipped"])
def test_exit_code_is_non_zero_unless_every_run_succeeded(wait, bad):
    poll = {
        "ci": [run(1, SHA_A, AFTER), run(3, SHA_B, AFTER)],
        "sanitizer": [run(2, SHA_A, AFTER), run(4, SHA_B, AFTER, conclusion=bad)],
    }
    rc, _ = go(wait, FakeGh([poll]), SHA_A, SHA_B, "--since", SINCE)
    assert rc == 1


def test_exit_code_zero_when_all_shas_and_workflows_succeed(wait):
    poll = {
        "ci": [run(1, SHA_A, AFTER), run(3, SHA_B, AFTER)],
        "sanitizer": [run(2, SHA_A, AFTER), run(4, SHA_B, AFTER)],
    }
    rc, _ = go(wait, FakeGh([poll]), SHA_A, SHA_B, "--since", SINCE)
    assert rc == 0


def test_timeout_exits_with_2(wait):
    nothing = {"ci": [], "sanitizer": []}
    rc, sleeps = go(wait, FakeGh([nothing]), SHA_A, "--since", SINCE, "--timeout-min", "3")
    assert rc == 2
    assert sleeps == [90, 90]  # 0 s, 90 s, 180 s polled; then the 3-minute budget is spent


def test_gh_error_is_retried_not_fatal(wait):
    err = {"ci": "ERROR", "sanitizer": "ERROR"}
    ok = {"ci": [run(1, SHA_A, AFTER)], "sanitizer": [run(2, SHA_A, AFTER)]}
    rc, sleeps = go(wait, FakeGh([err, ok]), SHA_A, "--since", SINCE)
    assert rc == 0
    assert sleeps == [90]


# --- API economy and repository --------------------------------------------------------------


def test_one_run_list_per_workflow_per_poll_for_all_shas(wait):
    shas = [c * 40 for c in "abcdef"]
    poll = {
        "ci": [run(i, s, AFTER) for i, s in enumerate(shas)],
        "sanitizer": [run(100 + i, s, AFTER) for i, s in enumerate(shas)],
    }
    fake = FakeGh([poll])
    rc, _ = go(wait, fake, *shas, "--since", SINCE)
    assert rc == 0
    calls = fake.run_list_calls()
    assert len(calls) == 2
    for c in calls:
        assert c[c.index("-R") + 1] == "owner/repo"
        assert c[c.index("-e") + 1] == "workflow_dispatch"
        assert c[c.index("-b") + 1] == "main"


def test_repo_comes_from_gh_repo_view_without_env(monkeypatch):
    monkeypatch.delenv("SSSP_REPO", raising=False)
    wait = load_wait()
    poll = {"ci": [run(1, SHA_A, AFTER)], "sanitizer": [run(2, SHA_A, AFTER)]}
    fake = FakeGh([poll], repo="someone/fork")
    rc, _ = go(wait, fake, SHA_A, "--since", SINCE)
    assert rc == 0
    assert ["gh", "repo", "view", "--json", "nameWithOwner", "-q", ".nameWithOwner"] in fake.calls
    assert all(c[c.index("-R") + 1] == "someone/fork" for c in fake.run_list_calls())
