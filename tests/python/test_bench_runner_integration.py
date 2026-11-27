"""The real C++ runner (build/debug/sssp_bench_runner) against the Python schema validator.

Runs tiny graphs only; nothing here is a performance measurement. Skipped (with the reason)
when the debug build has not been made.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest
from bench import schema
from bench_fixtures import ROOT, SCRIPTS

BUILD_DIR = ROOT / "build" / "debug"
RUNNER = BUILD_DIR / "sssp_bench_runner"

pytestmark = pytest.mark.skipif(
    not (RUNNER.is_file() and os.access(RUNNER, os.X_OK)),
    reason=f"{RUNNER.relative_to(ROOT)} is not built (cmake --preset debug && cmake --build "
    "--preset debug)",
)


def _py(*args, **kw):
    env = {**os.environ, "PYTHONPATH": str(SCRIPTS)}
    return subprocess.run(
        [sys.executable, *args], capture_output=True, text=True, env=env, timeout=120, **kw
    )


@pytest.fixture(scope="module")
def env_json(tmp_path_factory):
    path = tmp_path_factory.mktemp("env") / "debug.env.json"
    res = _py("-m", "bench.envinfo", "--build-dir", str(BUILD_DIR), "--out", str(path))
    assert res.returncode == 0, res.stderr
    return path


def _runner(*args):
    return subprocess.run(
        [str(RUNNER), *map(str, args)], capture_output=True, text=True, timeout=120
    )


def _validate(path):
    return _py("-m", "bench.schema", "validate", str(path))


def _records(path):
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line]


def test_baselines_on_small_graphs_pass_schema_validation(env_json, tmp_path):
    out = tmp_path / "SYNTHETIC-TEST-DATA.jsonl"
    common = ("--env-json", env_json, "--run-id", "SYNTHETIC-TEST-DATA-run", "--out", out)
    jobs = [
        (
            "--type",
            "er",
            "--n",
            "200",
            "--avg-degree",
            "4",
            "--weights",
            "uniform_int",
            "--algorithms",
            "dijkstra-binary,dijkstra-4ary,dijkstra-fibonacci,dijkstra-radix,bellman-ford",
            "--sources",
            "2",
            "--trials",
            "2",
            "--warmup",
            "1",
            "--order-seed",
            "12345",
        ),
        (
            "--type",
            "grid",
            "--n",
            "100",
            "--weights",
            "uniform_real",
            "--weight-type",
            "double",
            "--algorithms",
            "dijkstra-fibonacci",
            "--trials",
            "2",
            "--warmup",
            "0",
            "--isolate",
        ),
        (
            "--type",
            "road",
            "--n",
            "100",
            "--weights",
            "builtin",
            "--algorithms",
            "bellman-ford",
            "--verify",
            "bellman-ford",
            "--trials",
            "1",
            "--warmup",
            "0",
            "--isolate",
        ),
        (
            "--type",
            "ba",
            "--n",
            "100",
            "--attach",
            "2",
            "--algorithms",
            "dijkstra-4ary",
            "--verify",
            "none",
            "--trials",
            "1",
            "--warmup",
            "0",
        ),
    ]
    for job in jobs:
        res = _runner(*job, *common)
        assert res.returncode == 0, res.stderr
    res = _validate(out)
    assert res.returncode == 0, res.stderr + res.stdout
    records = _records(out)
    assert len(records) == 2 * 3 * 5 + 2 + 1 + 1
    env = json.loads(env_json.read_text(encoding="utf-8"))
    for rec in records:
        assert list(rec) == list(schema.REQUIRED_FIELDS)
        for key in schema.ENV_FIELDS:
            assert rec[key] == env[key]
        assert rec["status"] == "ok"
        # The debug build compiles the counters in.
        assert rec["diagnostics_enabled"] is True
        assert rec["ops"] is not None
    assert {r["algorithm"] for r in records} == set(schema.BASELINE_ALGORITHMS)
    isolated = [r for r in records if r["isolated"]]
    assert isolated and all(r["peak_rss_bytes"] for r in isolated)
    assert all(r["peak_rss_bytes"] is None for r in records if not r["isolated"])


@pytest.mark.parametrize(
    ("args", "code"),
    [
        (("--type", "er", "--algorithms", "dijkstra-binary", "--trials", "0"), 2),
        (("--type", "er", "--algorithms", "dijkstra-binary,dijkstra-4ary", "--isolate"), 2),
        (
            (
                "--type",
                "er",
                "--algorithms",
                "dijkstra-binary",
                "--sources",
                "1",
                "--source-list",
                "0",
            ),
            2,
        ),
        (
            (
                "--type",
                "er",
                "--algorithms",
                "dijkstra-radix",
                "--weights",
                "uniform_real",
                "--weight-type",
                "double",
            ),
            3,
        ),
        (("--type", "er", "--algorithms", "paper"), 3),
    ],
)
def test_refusals_exit_with_documented_codes_and_write_nothing(env_json, tmp_path, args, code):
    out = tmp_path / "out.jsonl"
    res = _runner(*args, "--n", "50", "--env-json", env_json, "--run-id", "R", "--out", out)
    assert res.returncode == code, res.stderr
    assert res.stderr.strip()
    assert not out.exists()


def test_env_json_from_another_build_is_rejected(env_json, tmp_path):
    other = json.loads(env_json.read_text(encoding="utf-8"))
    other["diagnostics_enabled"] = not other["diagnostics_enabled"]
    path = tmp_path / "other.env.json"
    path.write_text(json.dumps(other), encoding="utf-8")
    out = tmp_path / "out.jsonl"
    res = _runner(
        "--type",
        "er",
        "--n",
        "50",
        "--algorithms",
        "dijkstra-binary",
        "--env-json",
        path,
        "--run-id",
        "R",
        "--out",
        out,
    )
    assert res.returncode == 2
    assert "diagnostics_enabled" in res.stderr
    assert not out.exists()


def test_list_matches_schema_names():
    res = _runner("--list")
    assert res.returncode == 0, res.stderr
    listing = json.loads(res.stdout)
    for algo in listing["algorithms"]:
        assert schema.ALGORITHM_MODES[algo["name"]] == algo["mode"]
        if algo["name"] in schema.INTEGER_ONLY_ALGORITHMS:
            assert algo["weight_types"] == ["uint64"]
    res = _runner("--version")
    assert res.returncode == 0
    assert schema.PAPER_VERSION in res.stdout
