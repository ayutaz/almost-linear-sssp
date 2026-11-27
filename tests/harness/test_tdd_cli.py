"""scripts/tdd.py: TDD evidence CLI (start/red/green/refactor/done/status/check/log)."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from harness_util import GREEN_EVENT, RED_EVENT, TDD_PY, write_task


def tdd(root: Path, *args: str) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ, TDD_ROOT=str(root))
    return subprocess.run(
        [sys.executable, str(TDD_PY), *args],
        capture_output=True,
        text=True,
        cwd=str(root),
        env=env,
        timeout=60,
    )


def task_json(root: Path, tid: str) -> dict:
    return json.loads((root / "docs/tdd/tasks" / f"{tid}.json").read_text())


def active(root: Path) -> dict:
    path = root / ".tdd/runtime/active.json"
    return json.loads(path.read_text())["tasks"] if path.exists() else {}


FAIL = [sys.executable, "-c", "import sys; print('boom'); sys.exit(3)"]
PASS = [sys.executable, "-c", "print('ok')"]


def start(root: Path, tid: str = "T-1") -> subprocess.CompletedProcess[str]:
    return tdd(
        root,
        "start",
        tid,
        "--spec",
        "adds a",
        "--paper",
        "Sec 2",
        "--impl",
        "src/core/*.cpp",
        "include/sssp/a.hpp",
        "--test",
        "tests/unit/test_a.cpp",
    )


# ------------------------------------------------------------------- start
def test_start_creates_task_file_and_registers_active(repo: Path) -> None:
    res = start(repo)
    assert res.returncode == 0, res.stderr
    t = task_json(repo, "T-1")
    assert t["id"] == "T-1"
    assert t["spec"] == "adds a"
    assert t["paper_ref"] == "Sec 2"
    assert t["files"] == {
        "impl": ["src/core/*.cpp", "include/sssp/a.hpp"],
        "tests": ["tests/unit/test_a.cpp"],
    }
    assert t["phase"] == "RED"
    assert t["events"] == []
    assert t["commit"] is None
    assert "T-1" in active(repo)
    assert active(repo)["T-1"]["allow_test_reduction"] is False


def test_start_allows_several_active_tasks(repo: Path) -> None:
    assert start(repo, "T-1").returncode == 0
    assert start(repo, "T-2").returncode == 0
    assert set(active(repo)) == {"T-1", "T-2"}


def test_start_existing_task_reactivates_without_losing_events(repo: Path) -> None:
    write_task(repo, {"id": "T-1", "events": [RED_EVENT]}, active=False)
    res = tdd(repo, "start", "T-1")
    assert res.returncode == 0, res.stderr
    assert task_json(repo, "T-1")["events"] == [RED_EVENT]
    assert "T-1" in active(repo)


def test_start_new_task_requires_spec_impl_and_test(repo: Path) -> None:
    res = tdd(repo, "start", "T-9", "--spec", "x")
    assert res.returncode != 0
    assert not (repo / "docs/tdd/tasks/T-9.json").exists()


@pytest.mark.parametrize("bad", ["../evil", "a b", "", "x/y"])
def test_start_rejects_bad_ids(repo: Path, bad: str) -> None:
    res = tdd(repo, "start", bad, "--spec", "s", "--impl", "src/a", "--test", "tests/a")
    assert res.returncode != 0


# --------------------------------------------------------------------- red
def test_red_records_failing_command(repo: Path) -> None:
    start(repo)
    res = tdd(repo, "red", "T-1", "--reason", "a() missing", "--tests", "case a", "--", *FAIL)
    assert res.returncode == 0, res.stderr
    ev = task_json(repo, "T-1")["events"][-1]
    assert ev["phase"] == "RED"
    assert ev["exit_code"] == 3
    assert ev["reason"] == "a() missing"
    assert ev["tests"] == ["case a"]
    assert "boom" in ev["output_tail"]
    assert ev["timestamp"]
    assert "import sys" in ev["command"]
    assert task_json(repo, "T-1")["phase"] == "RED"


def test_red_rejects_passing_command(repo: Path) -> None:
    start(repo)
    res = tdd(repo, "red", "T-1", "--reason", "r", "--", *PASS)
    assert res.returncode != 0
    assert "RED" in (res.stderr + res.stdout)
    assert task_json(repo, "T-1")["events"] == []


def test_red_requires_reason(repo: Path) -> None:
    start(repo)
    assert tdd(repo, "red", "T-1", "--", *FAIL).returncode != 0


def test_red_keeps_only_last_20_lines(repo: Path) -> None:
    start(repo)
    cmd = [sys.executable, "-c", "import sys\nfor i in range(50): print(f'line{i}')\nsys.exit(1)"]
    assert tdd(repo, "red", "T-1", "--reason", "r", "--", *cmd).returncode == 0
    tail = task_json(repo, "T-1")["events"][-1]["output_tail"].splitlines()
    assert len(tail) == 20
    assert tail[-1] == "line49"


def test_single_string_command_runs_through_shell(repo: Path) -> None:
    start(repo)
    res = tdd(repo, "red", "T-1", "--reason", "r", "--", "echo hi && exit 4")
    assert res.returncode == 0, res.stderr
    ev = task_json(repo, "T-1")["events"][-1]
    assert ev["exit_code"] == 4
    assert ev["command"] == "echo hi && exit 4"


def test_red_on_unknown_task_fails(repo: Path) -> None:
    assert tdd(repo, "red", "NOPE", "--reason", "r", "--", *FAIL).returncode != 0


def test_red_resets_test_reduction_permission(repo: Path) -> None:
    start(repo)
    tdd(repo, "red", "T-1", "--reason", "r", "--", *FAIL)
    tdd(repo, "green", "T-1", "--", *PASS)
    assert tdd(repo, "allow-test-change", "T-1", "--reason", "dedupe").returncode == 0
    assert active(repo)["T-1"]["allow_test_reduction"] is True
    tdd(repo, "red", "T-1", "--reason", "r2", "--", *FAIL)
    assert active(repo)["T-1"]["allow_test_reduction"] is False


# ------------------------------------------------------------ green/refactor
def test_green_requires_prior_red(repo: Path) -> None:
    start(repo)
    res = tdd(repo, "green", "T-1", "--", *PASS)
    assert res.returncode != 0
    assert task_json(repo, "T-1")["events"] == []


def test_green_records_passing_command(repo: Path) -> None:
    start(repo)
    tdd(repo, "red", "T-1", "--reason", "r", "--", *FAIL)
    res = tdd(repo, "green", "T-1", "--tests", "case a", "--", *PASS)
    assert res.returncode == 0, res.stderr
    t = task_json(repo, "T-1")
    assert [e["phase"] for e in t["events"]] == ["RED", "GREEN"]
    assert t["events"][-1]["exit_code"] == 0
    assert "ok" in t["events"][-1]["output_tail"]
    assert t["phase"] == "GREEN"


def test_green_rejects_failing_command(repo: Path) -> None:
    start(repo)
    tdd(repo, "red", "T-1", "--reason", "r", "--", *FAIL)
    res = tdd(repo, "green", "T-1", "--", *FAIL)
    assert res.returncode != 0
    assert [e["phase"] for e in task_json(repo, "T-1")["events"]] == ["RED"]


def test_refactor_requires_green_and_passing(repo: Path) -> None:
    start(repo)
    tdd(repo, "red", "T-1", "--reason", "r", "--", *FAIL)
    assert tdd(repo, "refactor", "T-1", "--", *PASS).returncode != 0
    tdd(repo, "green", "T-1", "--", *PASS)
    assert tdd(repo, "refactor", "T-1", "--", *FAIL).returncode != 0
    assert tdd(repo, "refactor", "T-1", "--", *PASS).returncode == 0
    t = task_json(repo, "T-1")
    assert t["phase"] == "REFACTOR"
    assert t["events"][-1]["phase"] == "REFACTOR"


def test_runs_are_saved_under_runtime(repo: Path) -> None:
    start(repo)
    tdd(repo, "red", "T-1", "--reason", "r", "--", *FAIL)
    runs = list((repo / ".tdd/runtime/runs").glob("T-1-*-RED.log"))
    assert len(runs) == 1
    assert "boom" in runs[0].read_text()
    last = json.loads((repo / ".tdd/runtime/last_run.json").read_text())
    assert last["task"] == "T-1"
    assert last["phase"] == "RED"
    assert last["exit_code"] == 3


# -------------------------------------------------------- allow-test-change
def test_allow_test_change_requires_green(repo: Path) -> None:
    start(repo)
    tdd(repo, "red", "T-1", "--reason", "r", "--", *FAIL)
    assert tdd(repo, "allow-test-change", "T-1", "--reason", "x").returncode != 0
    tdd(repo, "green", "T-1", "--", *PASS)
    assert tdd(repo, "allow-test-change", "T-1").returncode != 0  # reason required
    assert tdd(repo, "allow-test-change", "T-1", "--reason", "merge dup").returncode == 0
    a = active(repo)["T-1"]
    assert a["allow_test_reduction"] is True
    assert a["allow_reason"] == "merge dup"
    log = task_json(repo, "T-1")["test_reduction_allowed"]
    assert log[-1]["reason"] == "merge dup"


# --------------------------------------------------------------------- done
def test_done_requires_green(repo: Path) -> None:
    start(repo)
    tdd(repo, "red", "T-1", "--reason", "r", "--", *FAIL)
    assert tdd(repo, "done", "T-1").returncode != 0


def test_done_marks_task_and_deactivates(repo: Path) -> None:
    start(repo)
    tdd(repo, "red", "T-1", "--reason", "r", "--", *FAIL)
    tdd(repo, "green", "T-1", "--", *PASS)
    res = tdd(repo, "done", "T-1", "--commit", "abc1234")
    assert res.returncode == 0, res.stderr
    t = task_json(repo, "T-1")
    assert t["phase"] == "DONE"
    assert t["commit"] == "abc1234"
    assert "T-1" not in active(repo)


# ------------------------------------------------------------------- status
def test_status_lists_active_tasks(repo: Path) -> None:
    start(repo, "T-1")
    tdd(repo, "red", "T-1", "--reason", "r", "--", *FAIL)
    res = tdd(repo, "status")
    assert res.returncode == 0
    assert "T-1" in res.stdout
    assert "RED" in res.stdout


def test_status_without_tasks(repo: Path) -> None:
    res = tdd(repo, "status")
    assert res.returncode == 0


# -------------------------------------------------------------------- check
def test_check_accepts_valid_tasks(repo: Path) -> None:
    write_task(repo, {"id": "OK-1", "phase": "GREEN", "events": [RED_EVENT, GREEN_EVENT]})
    write_task(repo, {"id": "OK-2", "phase": "RED", "events": []})
    res = tdd(repo, "check")
    assert res.returncode == 0, res.stdout + res.stderr


@pytest.mark.parametrize("name", ["P0-VERSION", "P0-GRAPH-CSR"])
def test_check_accepts_the_reference_evidence_files(repo: Path, name: str) -> None:
    """The hand-written examples the README format was derived from must stay valid."""
    src = TDD_PY.parents[1] / "docs/tdd/tasks" / f"{name}.json"
    (repo / "docs/tdd/tasks" / f"{name}.json").write_text(src.read_text())
    res = tdd(repo, "check")
    assert res.returncode == 0, res.stdout + res.stderr


@pytest.mark.parametrize(
    ("task", "needle"),
    [
        ({"id": "B-1", "phase": "GREEN", "events": [GREEN_EVENT]}, "RED"),
        ({"id": "B-1", "phase": "GREEN", "events": [GREEN_EVENT, RED_EVENT]}, "RED"),
        ({"id": "B-1", "phase": "RED", "events": [{**RED_EVENT, "exit_code": 0}]}, "exit"),
        (
            {"id": "B-1", "phase": "GREEN", "events": [RED_EVENT, {**GREEN_EVENT, "exit_code": 2}]},
            "exit",
        ),
        ({"id": "B-1", "phase": "WEIRD", "events": []}, "phase"),
        ({"id": "B-1", "phase": "DONE", "events": [RED_EVENT]}, "DONE"),
        ({"id": "B-1", "phase": "RED", "events": [{**RED_EVENT, "reason": ""}]}, "reason"),
        (
            {"id": "B-1", "phase": "RED", "events": [{**RED_EVENT, "timestamp": "yesterday"}]},
            "time",
        ),
        ({"id": "OTHER", "phase": "RED", "events": []}, "id"),
    ],
)
def test_check_rejects_invalid(repo: Path, task: dict, needle: str) -> None:
    write_task(repo, task, active=False)
    if task["id"] == "OTHER":
        (repo / "docs/tdd/tasks/OTHER.json").rename(repo / "docs/tdd/tasks/B-1.json")
    res = tdd(repo, "check")
    assert res.returncode != 0
    assert "B-1" in res.stdout + res.stderr
    assert needle.lower() in (res.stdout + res.stderr).lower()


def test_check_rejects_missing_keys_and_bad_json(repo: Path) -> None:
    (repo / "docs/tdd/tasks/X-1.json").write_text('{"id": "X-1"}')
    (repo / "docs/tdd/tasks/X-2.json").write_text("{not json")
    res = tdd(repo, "check")
    assert res.returncode != 0
    assert "X-1" in res.stdout + res.stderr
    assert "X-2" in res.stdout + res.stderr


def test_check_accepts_both_timestamp_offsets(repo: Path) -> None:
    ev1 = {**RED_EVENT, "timestamp": "2026-11-27T10:15:13+0900"}
    ev2 = {**GREEN_EVENT, "timestamp": "2026-11-27T10:15:14+09:00"}
    write_task(repo, {"id": "TS-1", "phase": "GREEN", "events": [ev1, ev2]})
    assert tdd(repo, "check").returncode == 0


# ---------------------------------------------------------------------- log
def test_log_regenerates_tdd_log(repo: Path) -> None:
    write_task(repo, {"id": "L-1", "phase": "GREEN", "events": [RED_EVENT, GREEN_EVENT]})
    write_task(repo, {"id": "L-2", "phase": "RED", "spec": "second | pipe", "events": []})
    res = tdd(repo, "log")
    assert res.returncode == 0, res.stderr
    text = (repo / "docs/TDD_LOG.md").read_text()
    assert "L-1" in text
    assert "L-2" in text
    assert "feature absent" in text  # RED reason is shown
    assert "second \\| pipe" in text  # table cells are escaped
    first = text
    assert tdd(repo, "log").returncode == 0
    assert (repo / "docs/TDD_LOG.md").read_text() == first  # deterministic


def test_unknown_command_fails(repo: Path) -> None:
    assert tdd(repo, "frobnicate").returncode != 0


def test_recorded_output_does_not_leak_absolute_paths(repo: Path) -> None:
    start(repo)
    cmd = [
        sys.executable,
        "-c",
        f"import sys; print({str(repo)!r} + '/src/a.cpp:1: error'); sys.exit(1)",
    ]
    assert tdd(repo, "red", "T-1", "--reason", "r", "--", *cmd).returncode == 0
    ev = task_json(repo, "T-1")["events"][-1]
    assert str(repo) not in ev["output_tail"]
    assert "./src/a.cpp:1: error" in ev["output_tail"]
    assert str(Path.home()) not in ev["output_tail"] + ev["command"]


def test_recorded_output_hides_temp_dirs(repo: Path, tmp_path: Path) -> None:
    other = tmp_path / "elsewhere" / "x.log"  # outside the root, inside $TMPDIR
    start(repo)
    cmd = [sys.executable, "-c", f"import sys; print({str(other)!r}); sys.exit(1)"]
    assert tdd(repo, "red", "T-1", "--reason", "r", "--", *cmd).returncode == 0
    tail = task_json(repo, "T-1")["events"][-1]["output_tail"]
    assert str(tmp_path) not in tail
    assert os.path.realpath(tmp_path) not in tail
    assert "$TMPDIR/" in tail


def test_recorded_output_hides_the_user_name(repo: Path) -> None:
    import getpass

    user = getpass.getuser()
    start(repo)
    cmd = [sys.executable, "-c", f"import sys; print('/x/pytest-of-{user}/p-1'); sys.exit(1)"]
    assert tdd(repo, "red", "T-1", "--reason", "r", "--", *cmd).returncode == 0
    tail = task_json(repo, "T-1")["events"][-1]["output_tail"]
    assert f"pytest-of-{user}" not in tail
    assert "pytest-of-$USER" in tail


# ---------------------------------------------------------------- test-only
def _git_commit(root: Path, path: str, text: str) -> str:
    (root / path).parent.mkdir(parents=True, exist_ok=True)
    (root / path).write_text(text)
    subprocess.run(["git", "-C", str(root), "add", path], check=True)
    subprocess.run(["git", "-C", str(root), "commit", "-qm", f"add {path}"], check=True)
    return subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()


TEST_ONLY = {
    "kind": "test-only",
    "test_only_reason": "verifies existing baselines; no new behaviour",
}


def test_check_accepts_test_only_task_starting_with_green(repo: Path) -> None:
    write_task(repo, {"id": "TO-1", "phase": "GREEN", "events": [GREEN_EVENT], **TEST_ONLY})
    res = tdd(repo, "check")
    assert res.returncode == 0, res.stdout + res.stderr


def test_check_accepts_test_only_with_covers_list(repo: Path) -> None:
    task = {
        "id": "TO-1",
        "phase": "GREEN",
        "events": [GREEN_EVENT],
        "files": {"impl": [], "tests": ["tests/integration/x.cpp"], "covers": ["src/graph/io.cpp"]},
        **TEST_ONLY,
    }
    write_task(repo, task)
    assert tdd(repo, "check").returncode == 0


@pytest.mark.parametrize(
    ("extra", "needle"),
    [
        ({"test_only_reason": ""}, "test_only_reason"),
        ({"test_only_reason": None}, "test_only_reason"),
        ({"files": {"impl": ["src/core/a.cpp"], "tests": []}}, "impl"),
        ({"kind": "whatever"}, "kind"),
        ({"files": {"impl": [], "tests": [], "covers": "src/x.cpp"}}, "covers"),
    ],
)
def test_check_rejects_bad_test_only_tasks(repo: Path, extra: dict, needle: str) -> None:
    write_task(
        repo, {"id": "TO-1", "phase": "GREEN", "events": [GREEN_EVENT], **TEST_ONLY, **extra}
    )
    res = tdd(repo, "check")
    assert res.returncode != 0
    assert "TO-1" in res.stdout
    assert needle in res.stdout


def test_check_still_rejects_green_first_without_kind(repo: Path) -> None:
    write_task(repo, {"id": "G-1", "phase": "GREEN", "events": [GREEN_EVENT]})
    res = tdd(repo, "check")
    assert res.returncode != 0
    assert "first event must be RED" in res.stdout


def test_test_only_commit_touching_only_tests_is_valid(repo: Path) -> None:
    sha = _git_commit(repo, "tests/unit/test_more.cpp", 'TEST_CASE("m") {}\n')
    task = {"id": "TO-1", "phase": "DONE", "events": [GREEN_EVENT], "commit": sha, **TEST_ONLY}
    task["files"] = {"impl": ["src/core/a.cpp"], "tests": ["tests/unit/test_more.cpp"]}
    write_task(repo, task, active=False)
    res = tdd(repo, "check")
    assert res.returncode == 0, res.stdout + res.stderr


def test_test_only_commit_touching_impl_is_invalid(repo: Path) -> None:
    sha = _git_commit(repo, "src/core/b.cpp", "int b() { return 2; }\n")
    task = {"id": "TO-1", "phase": "DONE", "events": [GREEN_EVENT], "commit": sha, **TEST_ONLY}
    task["files"] = {"impl": ["src/core/b.cpp"], "tests": []}
    write_task(repo, task, active=False)
    res = tdd(repo, "check")
    assert res.returncode != 0
    assert "src/core/b.cpp" in res.stdout


def test_test_only_unknown_commit_is_invalid(repo: Path) -> None:
    task = {"id": "TO-1", "phase": "DONE", "events": [GREEN_EVENT], "commit": "0" * 40, **TEST_ONLY}
    task["files"] = {"impl": ["src/core/a.cpp"], "tests": []}
    write_task(repo, task, active=False)
    res = tdd(repo, "check")
    assert res.returncode != 0
    assert "commit" in res.stdout


def test_start_and_green_for_test_only_task(repo: Path) -> None:
    res = tdd(
        repo,
        "start",
        "TO-2",
        "--spec",
        "round-trip of existing formats",
        "--test",
        "tests/integration/test_rt.cpp",
        "--test-only",
        "composes components that were each developed test-first",
    )
    assert res.returncode == 0, res.stderr
    t = task_json(repo, "TO-2")
    assert t["kind"] == "test-only"
    assert t["test_only_reason"].startswith("composes")
    assert t["files"]["impl"] == []
    res = tdd(repo, "green", "TO-2", "--", *PASS)
    assert res.returncode == 0, res.stderr
    assert tdd(repo, "check").returncode == 0
    assert tdd(repo, "done", "TO-2").returncode == 0


def test_green_without_red_still_rejected_for_normal_task(repo: Path) -> None:
    start(repo)
    assert tdd(repo, "green", "T-1", "--", *PASS).returncode != 0


def _public_safety_rules() -> list:
    import importlib.util

    path = TDD_PY.parent / "check_public_safety.py"
    spec = importlib.util.spec_from_file_location("_check_public_safety", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return [(name, pat) for name, pat in mod.LINE_RULES if "path" in name]


LEAKY_OUTPUT = "\n".join(  # public-safety: allow
    [
        "E   assert '/private/va...cpp:1: error' not in x",  # public-safety: allow
        "E   '/private/var/folde...not_leak_abso0/proj' is contained here:",  # public-safety: allow
        "E   /var/folders/ab/cd1234/T/pytest-of-someone/pytest-9/t0/proj",  # public-safety: allow
        "E   /private/tmp/claude-1/scratch/x.log",  # public-safety: allow
        "E   /Users/someone/work/x.cpp and /home/other/y.py",  # public-safety: allow
    ]
)


def test_recorded_output_passes_public_safety_path_rules(repo: Path) -> None:
    start(repo)
    script = f"import sys; print({LEAKY_OUTPUT!r}); sys.exit(1)"
    assert (
        tdd(repo, "red", "T-1", "--reason", "r", "--", sys.executable, "-c", script).returncode == 0
    )
    raw = (repo / "docs/tdd/tasks/T-1.json").read_text()
    for line in raw.splitlines():
        for name, pat in _public_safety_rules():
            assert not pat.search(line), f"{name}: {line}"
    tail = task_json(repo, "T-1")["events"][-1]["output_tail"]
    assert "x.cpp" in tail  # only the machine-specific prefix is removed
    assert "$TMPDIR" in tail
