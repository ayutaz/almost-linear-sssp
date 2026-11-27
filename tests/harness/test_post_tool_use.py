"""PostToolUse hook: format the edited file and run a light category check.

No real C++ build happens here: fake ``cmake`` / ``ctest`` / ``clang-format`` /
``ruff`` scripts log their arguments and exit with a configurable code.
"""

from __future__ import annotations

import fcntl
import json
import os
import time
from pathlib import Path
from typing import Any

import pytest

FAKE_TOOL = """#!/bin/sh
name=$(basename "$0")
printf '%s' "$name" >> "$FAKE_LOG"
for a in "$@"; do printf '\\t%s' "$a" >> "$FAKE_LOG"; done
printf '\\n' >> "$FAKE_LOG"
upper=$(echo "$name" | tr 'a-z-' 'A-Z_')
eval "sl=\\${FAKE_${upper}_SLEEP:-0}"
eval "out=\\${FAKE_${upper}_OUT:-}"
eval "code=\\${FAKE_${upper}_EXIT:-0}"
if [ "$sl" != "0" ]; then sleep "$sl"; fi
if [ -n "$out" ]; then printf '%s\\n' "$out"; fi
exit "$code"
"""


@pytest.fixture
def fake(tmp_path: Path) -> dict[str, Any]:
    """Fake tool directory + env pointing the hook at it."""
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    for name in ("cmake", "ctest", "clang-format", "ruff"):
        p = bindir / name
        p.write_text(FAKE_TOOL)
        p.chmod(0o755)
    log = tmp_path / "calls.log"
    log.write_text("")
    env = {
        "PATH": f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}",
        "FAKE_LOG": str(log),
        "TDD_CLANG_FORMAT": str(bindir / "clang-format"),
        "TDD_RUFF": str(bindir / "ruff"),
    }
    return {"bin": bindir, "log": log, "env": env}


def calls(fake: dict[str, Any]) -> list[list[str]]:
    text = Path(fake["log"]).read_text()
    return [line.split("\t") for line in text.splitlines() if line]


def calls_of(fake: dict[str, Any], name: str) -> list[list[str]]:
    return [c[1:] for c in calls(fake) if c[0] == name]


def configured(repo: Path, preset: str = "debug") -> Path:
    d = repo / "build" / preset
    d.mkdir(parents=True, exist_ok=True)
    (d / "build.ninja").write_text("# fake\n")
    return d


def payload(path: Path, tool: str = "Edit") -> dict[str, Any]:
    return {
        "hook_event_name": "PostToolUse",
        "tool_name": tool,
        "tool_input": {"file_path": str(path)},
        "tool_response": {"success": True},
        "tool_use_id": "toolu_1",
    }


def last_check(repo: Path) -> dict[str, Any]:
    return json.loads((repo / ".tdd/runtime/last_check.json").read_text())


def go(run_hook, repo: Path, fake: dict[str, Any], rel: str, extra=None, **kw):
    env = dict(fake["env"])
    env.update(extra or {})
    return run_hook("post_tool_use", payload(repo / rel, **kw), repo, env=env)


# ------------------------------------------------------------------ formatting
def test_cpp_edit_runs_clang_format_on_that_file_only(run_hook, repo, fake):
    r = go(run_hook, repo, fake, "src/core/a.cpp")
    assert r.code == 0
    cf = calls_of(fake, "clang-format")
    assert len(cf) == 1
    assert "-i" in cf[0] and "--style=file" in cf[0]
    assert Path(cf[0][-1]).resolve() == (repo / "src/core/a.cpp").resolve()
    assert calls_of(fake, "ruff") == []


def test_python_edit_runs_ruff_format_then_check(run_hook, repo, fake):
    (repo / "scripts").mkdir()
    (repo / "scripts/x.py").write_text("x=1\n")
    r = go(run_hook, repo, fake, "scripts/x.py", tool="Write")
    assert r.code == 0
    rf = calls_of(fake, "ruff")
    assert [c[0] for c in rf] == ["format", "check"]
    assert all(Path(c[-1]).resolve() == (repo / "scripts/x.py").resolve() for c in rf)
    assert calls_of(fake, "clang-format") == []
    assert calls_of(fake, "cmake") == []
    assert r.stdout.strip() == ""


def test_ruff_check_findings_reported_as_context(run_hook, repo, fake):
    (repo / "scripts").mkdir()
    (repo / "scripts/x.py").write_text("import os\n")
    r = go(
        run_hook,
        repo,
        fake,
        "scripts/x.py",
        extra={"FAKE_RUFF_EXIT": "1", "FAKE_RUFF_OUT": "scripts/x.py:1:8: F401 unused import"},
    )
    assert r.code == 0
    data = r.json
    assert data["hookSpecificOutput"]["hookEventName"] == "PostToolUse"
    assert "F401" in r.context
    assert "decision" not in data


def test_clang_format_error_reported_as_context(run_hook, repo, fake):
    r = go(
        run_hook,
        repo,
        fake,
        "src/core/a.cpp",
        extra={"FAKE_CLANG_FORMAT_EXIT": "1", "FAKE_CLANG_FORMAT_OUT": "bad .clang-format"},
    )
    assert r.code == 0
    assert "clang-format" in r.context
    assert "bad .clang-format" in r.context


def test_missing_formatter_is_skipped_silently(run_hook, repo, fake):
    r = go(
        run_hook,
        repo,
        fake,
        "src/core/a.cpp",
        extra={"TDD_CLANG_FORMAT": str(repo / "no-such-clang-format")},
    )
    assert r.code == 0
    assert r.stdout.strip() == ""
    assert calls_of(fake, "clang-format") == []


# ---------------------------------------------------------------- light check
def test_src_cpp_builds_unit_tests_only_and_records_pass(run_hook, repo, fake):
    bdir = configured(repo)
    r = go(run_hook, repo, fake, "src/core/a.cpp")
    assert r.code == 0
    assert r.stdout.strip() == ""

    cm = calls_of(fake, "cmake")
    assert len(cm) == 1
    assert cm[0][0] == "--build"
    assert Path(cm[0][1]).resolve() == bdir.resolve()
    assert cm[0][cm[0].index("--target") + 1] == "sssp_unit_tests"
    assert "all" not in cm[0]

    ct = calls_of(fake, "ctest")
    assert len(ct) == 1
    args = ct[0]
    assert Path(args[args.index("--test-dir") + 1]).resolve() == bdir.resolve()
    assert args[args.index("-L") + 1] == "^unit$"
    assert args[args.index("-LE") + 1] == "bench"
    assert "--output-on-failure" in args

    # Benchmarks are never built or run: "bench" only ever appears as the
    # ctest exclusion label.
    for c in calls(fake):
        for i, a in enumerate(c):
            if "bench" in a:
                assert c[0] == "ctest" and c[i - 1] == "-LE" and a == "bench"

    lc = last_check(repo)
    unit = lc["categories"]["unit"]
    assert unit["status"] == "pass"
    assert unit["file"] == "src/core/a.cpp"
    assert unit["exit_code"] == 0
    assert "ctest" in unit["command"] and "sssp_unit_tests" in unit["command"]
    assert isinstance(unit["output_tail"], str)
    assert unit["timestamp"] and lc["updated"]
    assert lc["last"]["category"] == "unit"
    assert lc["last"]["status"] == "pass"


@pytest.mark.parametrize("rel", ["include/sssp/x.hpp", "tools/cli.cc"])
def test_include_and_tools_map_to_unit(run_hook, repo, fake, rel):
    configured(repo)
    (repo / rel).parent.mkdir(parents=True, exist_ok=True)
    (repo / rel).write_text("// x\n")
    r = go(run_hook, repo, fake, rel)
    assert r.code == 0
    cm = calls_of(fake, "cmake")
    assert cm and cm[0][cm[0].index("--target") + 1] == "sssp_unit_tests"


def test_tests_category_from_path(run_hook, repo, fake):
    configured(repo)
    (repo / "tests/integration").mkdir()
    (repo / "tests/integration/test_i.cpp").write_text("// i\n")
    r = go(run_hook, repo, fake, "tests/integration/test_i.cpp")
    assert r.code == 0
    cm = calls_of(fake, "cmake")
    assert cm[0][cm[0].index("--target") + 1] == "sssp_integration_tests"
    ct = calls_of(fake, "ctest")
    assert ct[0][ct[0].index("-L") + 1] == "^integration$"
    assert last_check(repo)["categories"]["integration"]["status"] == "pass"


def test_category_without_cpp_tests_is_not_checked(run_hook, repo, fake):
    configured(repo)
    (repo / "tests/property").mkdir()
    (repo / "tests/property/helpers.hpp").write_text("// h\n")
    r = go(run_hook, repo, fake, "tests/property/helpers.hpp")
    assert r.code == 0
    assert calls_of(fake, "clang-format")
    assert calls_of(fake, "cmake") == []
    assert calls_of(fake, "ctest") == []


def test_benchmarks_cpp_is_formatted_but_never_checked(run_hook, repo, fake):
    configured(repo)
    (repo / "benchmarks").mkdir()
    (repo / "benchmarks/bench_x.cpp").write_text("// b\n")
    r = go(run_hook, repo, fake, "benchmarks/bench_x.cpp")
    assert r.code == 0
    assert calls_of(fake, "clang-format")
    assert calls_of(fake, "cmake") == []
    assert calls_of(fake, "ctest") == []


def test_build_failure_records_fail_and_reports(run_hook, repo, fake):
    configured(repo)
    r = go(
        run_hook,
        repo,
        fake,
        "src/core/a.cpp",
        extra={"FAKE_CMAKE_EXIT": "1", "FAKE_CMAKE_OUT": "a.cpp:1:1: error: boom"},
    )
    assert r.code == 0
    assert calls_of(fake, "ctest") == []
    unit = last_check(repo)["categories"]["unit"]
    assert unit["status"] == "fail"
    assert unit["exit_code"] == 1
    assert "error: boom" in unit["output_tail"]
    assert "cmake --build" in unit["command"]
    ctx = r.context
    assert r.json["hookSpecificOutput"]["hookEventName"] == "PostToolUse"
    assert "unit" in ctx and "error: boom" in ctx
    assert "cmake --build" in ctx  # reproduction command
    assert "decision" not in r.json


def test_test_failure_records_fail_and_reports(run_hook, repo, fake):
    configured(repo)
    r = go(
        run_hook,
        repo,
        fake,
        "tests/unit/test_a.cpp",
        extra={"FAKE_CTEST_EXIT": "8", "FAKE_CTEST_OUT": "unit:a ... Failed"},
    )
    assert r.code == 0
    unit = last_check(repo)["categories"]["unit"]
    assert unit["status"] == "fail"
    assert unit["exit_code"] == 8
    assert "Failed" in unit["output_tail"]
    assert unit["file"] == "tests/unit/test_a.cpp"
    assert "ctest" in r.context and "Failed" in r.context


def test_output_tail_is_last_20_lines(run_hook, repo, fake):
    configured(repo)
    many = "\n".join(f"line{i}" for i in range(50))
    r = go(
        run_hook,
        repo,
        fake,
        "src/core/a.cpp",
        extra={"FAKE_CMAKE_EXIT": "1", "FAKE_CMAKE_OUT": many},
    )
    assert r.code == 0
    tail = last_check(repo)["categories"]["unit"]["output_tail"].splitlines()
    assert len(tail) <= 20
    assert tail[-1] == "line49"


def test_no_build_ninja_means_no_cmake(run_hook, repo, fake):
    r = go(run_hook, repo, fake, "src/core/a.cpp")
    assert r.code == 0
    assert calls_of(fake, "cmake") == []
    assert calls_of(fake, "ctest") == []
    assert calls_of(fake, "clang-format")  # formatting still happens


def test_postcheck_env_zero_disables_check(run_hook, repo, fake):
    configured(repo)
    r = go(run_hook, repo, fake, "src/core/a.cpp", extra={"TDD_POSTCHECK": "0"})
    assert r.code == 0
    assert calls_of(fake, "cmake") == []


def test_config_disabled_and_preset(run_hook, repo, fake):
    configured(repo)
    cfg = {"enforce": False, "post_check": {"enabled": False}}
    (repo / ".tdd/config.json").write_text(json.dumps(cfg))
    go(run_hook, repo, fake, "src/core/a.cpp")
    assert calls_of(fake, "cmake") == []

    rel = configured(repo, "release")
    cfg = {"enforce": False, "post_check": {"enabled": True, "preset": "release"}}
    (repo / ".tdd/config.json").write_text(json.dumps(cfg))
    go(run_hook, repo, fake, "src/core/a.cpp")
    cm = calls_of(fake, "cmake")
    assert len(cm) == 1 and Path(cm[0][1]).resolve() == rel.resolve()


def test_timeout_records_timeout_and_hook_returns_in_time(run_hook, repo, fake):
    configured(repo)
    t0 = time.monotonic()
    r = go(
        run_hook,
        repo,
        fake,
        "src/core/a.cpp",
        extra={"FAKE_CMAKE_SLEEP": "20", "TDD_POSTCHECK_TIMEOUT": "1"},
    )
    elapsed = time.monotonic() - t0
    assert r.code == 0
    assert elapsed < 8, elapsed
    unit = last_check(repo)["categories"]["unit"]
    assert unit["status"] == "timeout"
    assert unit.get("reason")
    assert calls_of(fake, "ctest") == []
    assert "timeout" in r.context.lower() or "タイムアウト" in r.context


def test_other_categories_are_preserved(run_hook, repo, fake):
    configured(repo)
    rt = repo / ".tdd/runtime"
    rt.mkdir(parents=True, exist_ok=True)
    old = {
        "updated": "2026-01-01T00:00:00+09:00",
        "categories": {
            "integration": {
                "status": "fail",
                "timestamp": "2026-01-01T00:00:00+09:00",
                "file": "tests/integration/x.cpp",
                "command": "c",
                "exit_code": 8,
                "output_tail": "old",
            },
            "unit": {"status": "fail", "exit_code": 1},
        },
        "last": {},
    }
    (rt / "last_check.json").write_text(json.dumps(old))
    r = go(run_hook, repo, fake, "src/core/a.cpp")
    assert r.code == 0
    lc = last_check(repo)
    assert lc["categories"]["integration"] == old["categories"]["integration"]
    assert lc["categories"]["unit"]["status"] == "pass"
    assert lc["last"]["category"] == "unit"


def test_busy_lock_records_skipped(run_hook, repo, fake):
    configured(repo)
    rt = repo / ".tdd/runtime"
    rt.mkdir(parents=True, exist_ok=True)
    with open(rt / "postcheck.lock", "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        r = go(
            run_hook,
            repo,
            fake,
            "src/core/a.cpp",
            extra={"TDD_POSTCHECK_LOCK_WAIT": "0.5"},
        )
    assert r.code == 0
    assert calls_of(fake, "cmake") == []
    unit = last_check(repo)["categories"]["unit"]
    assert unit["status"] == "skipped"
    assert unit["reason"] == "busy"


# ------------------------------------------------------------------ robustness
def test_file_outside_root_is_ignored(run_hook, repo, fake, tmp_path):
    configured(repo)
    outside = tmp_path / "elsewhere" / "z.cpp"
    outside.parent.mkdir()
    outside.write_text("// z\n")
    r = run_hook("post_tool_use", payload(outside), repo, env=fake["env"])
    assert r.code == 0
    assert r.stdout.strip() == ""
    assert calls(fake) == []


def test_invalid_stdin_exits_zero(run_hook, repo, fake):
    r = run_hook("post_tool_use", {}, repo, env=fake["env"], raw_stdin="{not json")
    assert r.code == 0
    assert r.stdout.strip() == ""
    assert calls(fake) == []


def test_missing_file_path_exits_zero(run_hook, repo, fake):
    r = run_hook(
        "post_tool_use",
        {"tool_name": "Edit", "tool_input": {}},
        repo,
        env=fake["env"],
    )
    assert r.code == 0
    assert r.stdout.strip() == ""
