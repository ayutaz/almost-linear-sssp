"""scripts/run_benchmarks.sh driven with a FAKE runner (a shell stub, not the C++ binary).

The stub only logs its arguments or sleeps; it produces no measurements.
"""

from __future__ import annotations

import os
import stat
import subprocess
import sys

import pytest
from bench import plan, schema
from bench_fixtures import ROOT

SCRIPT = ROOT / "scripts" / "run_benchmarks.sh"

CMAKE_CACHE = (
    "CMAKE_BUILD_TYPE:STRING=Release\n"
    "CMAKE_CXX_COMPILER:FILEPATH=/usr/bin/c++\n"
    "CMAKE_CXX_FLAGS_RELEASE:STRING=-O3 -DNDEBUG\n"
    "SSSP_DIAGNOSTICS:BOOL={diag}\n"
)


@pytest.fixture
def build_root(tmp_path):
    root = tmp_path / "build"
    for preset, diag in (("release", "OFF"), ("release-diag", "ON")):
        (root / preset).mkdir(parents=True)
        (root / preset / "CMakeCache.txt").write_text(CMAKE_CACHE.format(diag=diag))
    return root


def _install_stub(build_root, body):
    for preset in ("release", "release-diag"):
        stub = build_root / preset / "sssp_bench_runner"
        stub.write_text("#!/bin/bash\n" + body)
        stub.chmod(stub.stat().st_mode | stat.S_IEXEC)


def _run(args, build_root, tmp_path, **env):
    full_env = {
        **os.environ,
        "BENCH_NO_BUILD": "1",
        "BENCH_BUILD_ROOT": str(build_root),
        "BENCH_RESULTS_DIR": str(tmp_path / "results"),
        "BENCH_PYTHON": sys.executable,
        "BENCH_RUN_ID": "SYNTHETIC-TEST-DATA-run",
        **env,
    }
    return subprocess.run(
        ["bash", str(SCRIPT), *args], capture_output=True, text=True, env=full_env, timeout=120
    )


def test_usage_error_without_profile(build_root, tmp_path):
    assert _run([], build_root, tmp_path).returncode == 2
    assert _run(["turbo"], build_root, tmp_path).returncode == 2


def test_missing_runner_exits_2_with_clear_message(build_root, tmp_path):
    res = _run(["smoke"], build_root, tmp_path)
    assert res.returncode == 2
    assert "sssp_bench_runner" in res.stderr
    assert "見つかりません" in res.stderr


def test_missing_runner_can_be_tolerated_in_smoke_only(build_root, tmp_path):
    res = _run(["smoke", "--allow-missing-runner"], build_root, tmp_path)
    assert res.returncode == 0
    assert "未実装" in res.stdout + res.stderr
    res = _run(["standard", "--allow-missing-runner"], build_root, tmp_path)
    assert res.returncode == 2


def test_one_isolated_process_per_job(build_root, tmp_path):
    log = tmp_path / "calls.log"
    _install_stub(build_root, f'echo "$0 $*" >> "{log}"\n')
    res = _run(["smoke"], build_root, tmp_path, BENCH_ORDER_SEED="5")
    assert res.returncode == 0, res.stderr
    calls = log.read_text().splitlines()
    jobs = plan.build_jobs("smoke", order_seed=5)
    assert len(calls) == len(jobs)
    out_dir = str(tmp_path / "results" / "raw")
    for call, job in zip(calls, jobs, strict=True):
        args = call.split()
        assert args[0].endswith(f"{job.build}/sssp_bench_runner")
        assert "--isolate" in args
        algos = args[args.index("--algorithms") + 1]
        assert algos == job.algorithm  # exactly one algorithm per process (Peak RSS)
        assert args[args.index("--out") + 1].startswith(out_dir)
        assert args[args.index("--n") + 1] == str(job.n)
        assert args[args.index("--type") + 1] == job.graph_type
        if job.shape_args != "-":
            assert job.shape_args in " ".join(args)
        for flag in ("--env-json", "--timeout-sec", "--order-seed", "--trials", "--warmup"):
            assert flag in args
    env_jsons = {c.split()[c.split().index("--env-json") + 1] for c in calls}
    assert len(env_jsons) == 2  # one per build (release / release-diag)


def test_timeout_is_recorded_as_a_failure_record(build_root, tmp_path):
    _install_stub(build_root, "exec sleep 30\n")
    res = _run(
        ["smoke"],
        build_root,
        tmp_path,
        BENCH_TIMEOUT="1",
        BENCH_ALGORITHMS="dijkstra-binary",
        BENCH_GRAPH_TYPES="er",
        BENCH_SIZES="1000",
        BENCH_BUILDS="release",
    )
    assert res.returncode == 0, res.stderr
    out = tmp_path / "results" / "raw" / "SYNTHETIC-TEST-DATA-run.jsonl"
    records, errors = schema.load_jsonl(out)
    assert errors == []
    assert [r["status"] for r in records] == ["timeout"]
    assert records[0]["algorithm"] == "dijkstra-binary"


def test_runner_crash_is_recorded_as_error(build_root, tmp_path):
    _install_stub(build_root, "exit 1\n")
    res = _run(
        ["smoke"],
        build_root,
        tmp_path,
        BENCH_ALGORITHMS="paper",
        BENCH_GRAPH_TYPES="er",
        BENCH_SIZES="1000",
        BENCH_BUILDS="release",
    )
    assert res.returncode == 1  # a crash is a bug: recorded, and the run fails
    out = tmp_path / "results" / "raw" / "SYNTHETIC-TEST-DATA-run.jsonl"
    records, errors = schema.load_jsonl(out)
    assert errors == []
    assert [r["status"] for r in records] == ["error"]
    assert "exit 1" in records[0]["error"]


def _python_failing_on(tmp_path, suffix):
    """A BENCH_PYTHON wrapper that fails only for invocations writing ``*<suffix>``."""
    wrapper = tmp_path / "python-wrapper"
    wrapper.write_text(
        "#!/bin/bash\n"
        'for a in "$@"; do\n'
        f'  case "$a" in *{suffix}) echo "simulated envinfo failure" >&2; exit 3 ;; esac\n'
        "done\n"
        f'exec "{sys.executable}" "$@"\n'
    )
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IEXEC)
    return str(wrapper)


_ONE_JOB = {
    "BENCH_ALGORITHMS": "dijkstra-binary",
    "BENCH_GRAPH_TYPES": "er",
    "BENCH_SIZES": "1000",
}


def test_end_environment_is_written_for_every_preset(build_root, tmp_path):
    _install_stub(build_root, "exit 0\n")
    res = _run(["smoke"], build_root, tmp_path, **_ONE_JOB)
    assert res.returncode == 0, res.stderr
    assert "env_end_failed=0" in res.stdout
    raw = tmp_path / "results" / "raw"
    for preset in ("release", "release-diag"):
        assert (raw / f"SYNTHETIC-TEST-DATA-run.{preset}.env-end.json").is_file()


def test_end_environment_failure_is_reported_and_fails_the_run(build_root, tmp_path):
    # Protocol 7.3: a run is invalid if CPU_Speed_Limit < 100 at the start OR the end, so a
    # missing end state means the run's validity cannot be established.
    _install_stub(build_root, "exit 0\n")
    res = _run(
        ["smoke"],
        build_root,
        tmp_path,
        BENCH_PYTHON=_python_failing_on(tmp_path, ".release-diag.env-end.json"),
        **_ONE_JOB,
    )
    assert res.returncode == 1, res.stdout + res.stderr
    assert "終了時の環境情報" in res.stderr
    assert "release-diag" in res.stderr
    assert "simulated envinfo failure" in res.stderr
    assert "env_end_failed=1" in res.stdout
    raw = tmp_path / "results" / "raw"
    assert (raw / "SYNTHETIC-TEST-DATA-run.release.env-end.json").is_file()
    assert not (raw / "SYNTHETIC-TEST-DATA-run.release-diag.env-end.json").exists()
