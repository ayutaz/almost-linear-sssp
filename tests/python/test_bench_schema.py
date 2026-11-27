"""Benchmark JSONL schema (docs/design/benchmark_schema.md). Inputs are SYNTHETIC records."""

from __future__ import annotations

import json
import subprocess
import sys

import pytest
from bench import schema
from bench_fixtures import SCRIPTS, synthetic_diag_record, synthetic_record, write_jsonl
from hypothesis import given
from hypothesis import strategies as st


def errors_of(rec):
    return schema.validate_record(rec)


def test_synthetic_valid_record_has_no_errors():
    assert errors_of(synthetic_record()) == []
    assert errors_of(synthetic_diag_record()) == []


@pytest.mark.parametrize("field", schema.REQUIRED_FIELDS)
def test_every_required_field_must_be_present(field):
    rec = synthetic_record()
    del rec[field]
    errs = errors_of(rec)
    assert any(field in e and "missing" in e for e in errs), errs


def test_unknown_field_is_rejected():
    errs = errors_of(synthetic_record(sssp_time=1))
    assert any("sssp_time" in e and "unknown" in e for e in errs)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("n", "1000"),
        ("n", 1000.0),
        ("n", True),
        ("n", 0),
        ("git_dirty", "no"),
        ("warmup", 0),
        ("mode", "fast"),
        ("status", "done"),
        ("weight_type", "u32"),
        ("schema_version", "0"),
        ("paper_version", "v1"),
        ("git_sha", "not-a-sha"),
        ("timestamp", "2026-11-27T10:00:00"),  # no timezone
        ("timestamp", "yesterday"),
        ("params", []),
        ("seed", -1),
        ("peak_rss_bytes", 0),
        ("ram_bytes", "16GB"),
        ("source", 1000),  # source must be < n
    ],
)
def test_bad_values_are_rejected(field, value):
    errs = errors_of(synthetic_record(**{field: value}))
    assert any(field in e for e in errs), errs


def test_unmeasured_time_must_be_null_not_zero():
    rec = synthetic_record()
    rec["times_ns"]["preprocess"] = 0
    errs = errors_of(rec)
    assert any("times_ns.preprocess" in e and "null" in e for e in errs), errs


def test_time_keys_are_exactly_the_documented_ones():
    rec = synthetic_record()
    del rec["times_ns"]["tree"]
    assert any("times_ns.tree" in e for e in errors_of(rec))
    rec = synthetic_record()
    rec["times_ns"]["total"] = 5
    assert any("times_ns.total" in e for e in errors_of(rec))
    rec = synthetic_record()
    rec["times_ns"]["sssp"] = 1.5
    assert any("times_ns.sssp" in e for e in errors_of(rec))


def test_ok_status_requires_core_times():
    rec = synthetic_record()
    rec["times_ns"]["sssp"] = None
    assert any("times_ns.sssp" in e for e in errors_of(rec))


def test_end_to_end_cannot_be_shorter_than_its_parts():
    rec = synthetic_record()
    rec["times_ns"]["sssp"] = 6000
    rec["times_ns"]["end_to_end"] = 5000
    assert any("end_to_end" in e for e in errors_of(rec))


def test_failure_record_may_leave_unknowns_null():
    rec = synthetic_record(
        status="timeout",
        error="killed after 60 s",
        m=None,
        source=None,
        trial=None,
        reachable=None,
        peak_rss_bytes=None,
        rss_after_graph_bytes=None,
        times_ns=dict.fromkeys(schema.TIME_KEYS),
        correctness={"verified": False, "ok": None, "method": None},
    )
    assert errors_of(rec) == []


def test_ok_record_may_not_leave_m_source_trial_null():
    for field in ("m", "source", "trial"):
        assert any(field in e for e in errors_of(synthetic_record(**{field: None})))


def test_counters_require_a_diagnostics_build():
    rec = synthetic_diag_record(diagnostics_enabled=False)
    errs = errors_of(rec)
    assert any("ops" in e and "diagnostics" in e for e in errs), errs


def test_op_counters_must_be_non_negative_integers():
    rec = synthetic_diag_record()
    rec["ops"]["relaxations"] = -1
    assert any("ops.relaxations" in e for e in errors_of(rec))


def test_peak_rss_only_from_isolated_process():
    errs = errors_of(synthetic_record(isolated=False))
    assert any("peak_rss_bytes" in e and "isolated" in e for e in errs), errs
    assert errors_of(synthetic_record(isolated=False, peak_rss_bytes=None)) == []


@pytest.mark.parametrize(
    "corr",
    [
        {"verified": False, "ok": True, "method": None},
        {"verified": True, "ok": None, "method": "reference:dijkstra-binary"},
        {"verified": True, "ok": True, "method": None},
        {"verified": "yes", "ok": True, "method": "x"},
    ],
)
def test_correctness_block_is_consistent(corr):
    assert any("correctness" in e for e in errors_of(synthetic_record(correctness=corr)))


def test_fallback_totals_must_match_records():
    fb = {
        "count": 2,
        "records": [{"reason": "outer: base case", "vertices": 10, "edges": 20}],
        "vertices": 10,
        "edges": 20,
    }
    assert any("fallback" in e for e in errors_of(synthetic_diag_record(fallback=fb)))
    fb["count"] = 1
    assert errors_of(synthetic_diag_record(fallback=fb)) == []
    fb["vertices"] = 11
    assert any("fallback.vertices" in e for e in errors_of(synthetic_diag_record(fallback=fb)))


def test_baseline_algorithms_have_baseline_mode():
    errs = errors_of(synthetic_record(mode="paper-faithful"))
    assert any("mode" in e for e in errs)
    errs = errors_of(synthetic_record(algorithm="paper", mode="baseline"))
    assert any("mode" in e for e in errs)
    assert errors_of(synthetic_record(algorithm="paper", mode="paper-faithful")) == []


@pytest.mark.parametrize("algorithm", ["paper-faithful", "paper-faithful-d1"])
def test_mode_a_entry_names_are_paper_faithful(algorithm):
    # Names of include/sssp/paper.hpp algorithm_name (registered by the C++ runner).
    assert schema.ALGORITHM_MODES[algorithm] == "paper-faithful"
    assert algorithm not in schema.BASELINE_ALGORITHMS
    assert errors_of(synthetic_record(algorithm=algorithm, mode="paper-faithful")) == []
    errs = errors_of(synthetic_record(algorithm=algorithm, mode="baseline"))
    assert any("mode" in e for e in errs)
    errs = errors_of(synthetic_record(algorithm=algorithm, mode="paper-experimental"))
    assert any("mode" in e for e in errs)


def test_radix_heap_only_with_integer_weights():
    errs = errors_of(synthetic_record(algorithm="dijkstra-radix", weight_type="double"))
    assert any("weight_type" in e for e in errs)


@given(st.dictionaries(st.sampled_from(schema.OPS_KEYS), st.integers(min_value=0)))
def test_any_non_negative_op_counts_are_accepted(ops):
    assert errors_of(synthetic_diag_record(ops=ops)) == []


@given(st.integers(min_value=1, max_value=10**9), st.data())
def test_source_in_range_is_accepted(n, data):
    source = data.draw(st.integers(min_value=0, max_value=n - 1))
    assert errors_of(synthetic_record(n=n, source=source, reachable=None)) == []


def test_load_jsonl_reports_line_numbers(tmp_path):
    path = tmp_path / "run.jsonl"
    good = json.dumps(synthetic_record())
    bad = json.dumps(synthetic_record(mode="fast"))
    path.write_text(f"{good}\n\n{{not json\n{bad}\n", encoding="utf-8")
    records, errors = schema.load_jsonl(path)
    assert len(records) == 1
    assert any(":3:" in e and "JSON" in e for e in errors), errors
    assert any(":4:" in e and "mode" in e for e in errors), errors


def test_failure_record_factory_produces_valid_record():
    env = {k: synthetic_record()[k] for k in schema.ENV_FIELDS}
    rec = schema.make_failure_record(
        env,
        run_id="SYNTHETIC-TEST-DATA-run",
        algorithm="paper",
        graph_type="er",
        n=1000,
        avg_degree=4.0,
        weight_dist="uniform_int",
        weight_type="uint64",
        seed=1,
        status="timeout",
        error="timeout after 1 s",
        order_seed=3,
        isolated=True,
    )
    assert errors_of(rec) == []
    assert rec["mode"] == "paper-faithful"
    assert all(v is None for v in rec["times_ns"].values())
    with pytest.raises(ValueError):
        schema.make_failure_record(env, **{**_failure_kwargs(), "status": "ok"})


def _failure_kwargs():
    return {
        "run_id": "SYNTHETIC-TEST-DATA-run",
        "algorithm": "paper",
        "graph_type": "er",
        "n": 1000,
        "avg_degree": 4.0,
        "weight_dist": "uniform_int",
        "weight_type": "uint64",
        "seed": 1,
        "status": "timeout",
        "error": "x",
    }


def test_validate_env_rejects_missing_and_bad_fields():
    env = {k: synthetic_record()[k] for k in schema.ENV_FIELDS}
    assert schema.validate_env(env) == []
    del env["git_sha"]
    env["ram_bytes"] = "lots"
    errs = schema.validate_env(env)
    assert any("git_sha" in e for e in errs)
    assert any("ram_bytes" in e for e in errs)


def _cli(*args):
    return subprocess.run(
        [sys.executable, "-m", "bench.schema", *args],
        cwd=SCRIPTS,
        capture_output=True,
        text=True,
    )


def test_cli_validate_exit_codes(tmp_path):
    good = write_jsonl(tmp_path / "good.jsonl", [synthetic_record()])
    bad = write_jsonl(tmp_path / "bad.jsonl", [synthetic_record(n=-5)])
    assert _cli("validate", str(good)).returncode == 0
    res = _cli("validate", str(bad))
    assert res.returncode == 1
    assert "n" in res.stderr


def test_cli_failure_appends_valid_record(tmp_path):
    env = {k: synthetic_record()[k] for k in schema.ENV_FIELDS}
    env_path = tmp_path / "env.json"
    env_path.write_text(json.dumps(env), encoding="utf-8")
    out = tmp_path / "run.jsonl"
    args = [
        "failure",
        "--env-json",
        str(env_path),
        "--out",
        str(out),
        "--run-id",
        "SYNTHETIC-TEST-DATA-run",
        "--algorithm",
        "dijkstra-binary",
        "--graph-type",
        "er",
        "--n",
        "1000",
        "--avg-degree",
        "4",
        "--weights",
        "uniform_int",
        "--weight-type",
        "uint64",
        "--seed",
        "1",
        "--status",
        "timeout",
        "--error",
        "timeout after 1 s",
    ]
    assert _cli(*args).returncode == 0
    assert _cli(*args).returncode == 0
    records, errors = schema.load_jsonl(out)
    assert errors == []
    assert len(records) == 2
    assert records[0]["status"] == "timeout"


def test_failed_verification_must_be_status_error():
    bad = {"verified": True, "ok": False, "method": "reference:dijkstra-binary"}
    errs = errors_of(synthetic_record(correctness=bad))
    assert any("status" in e and "verification" in e for e in errs), errs
    assert errors_of(synthetic_record(correctness=bad, status="error", error="dist mismatch")) == []


def test_cli_failure_params_keep_numbers(tmp_path):
    # the runner writes numeric params as numbers; failure records must group with them
    env = {k: synthetic_record()[k] for k in schema.ENV_FIELDS}
    env_path = tmp_path / "env.json"
    env_path.write_text(json.dumps(env), encoding="utf-8")
    out = tmp_path / "run.jsonl"
    res = _cli(
        "failure", "--env-json", str(env_path), "--out", str(out), "--run-id", "SYNTHETIC",
        "--algorithm", "paper", "--graph-type", "er", "--n", "10", "--weights", "uniform_int",
        "--weight-type", "uint64", "--seed", "1", "--status", "error", "--error", "x",
        "--param", "threshold=4", "--param", "variant=fast",
    )  # fmt: skip
    assert res.returncode == 0, res.stderr
    records, errors = schema.load_jsonl(out)
    assert errors == []
    assert records[0]["params"] == {"threshold": 4, "variant": "fast"}
