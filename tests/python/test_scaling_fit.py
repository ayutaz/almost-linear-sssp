"""scripts/scaling_fit.py on SYNTHETIC JSONL.

Every record here is fabricated test input (see bench_fixtures.py): the "times" and "operation
counts" are computed from a known formula so that the test can check that the fit recovers the
formula's coefficient and exponent. None of these numbers is a measurement.
"""

from __future__ import annotations

import math

import pandas as pd
import pytest
import scaling_fit
from bench_fixtures import synthetic_diag_record, synthetic_record, write_jsonl

# SYNTHETIC sizes: 4 sizes, n spans 2^7 = 128x (> 2 orders of magnitude).
SIZES = (1 << 10, 1 << 12, 1 << 14, 1 << 17)
JITTER = (0.98, 1.0, 1.02)  # three trials whose median is the formula's value


def _lll(n):
    return math.log2(math.log2(math.log2(n)))


def _times(sssp):
    return {
        "graph_gen": 1000,
        "preprocess": None,
        "sssp": sssp,
        "tree": None,
        "end_to_end": sssp,
        "verify": 10,
    }


def _mode(algorithm):
    return "baseline" if algorithm.startswith("dijkstra-") else "paper-faithful"


def _all_fallback(n, m):
    return {
        "count": 1,
        "records": [{"reason": "outer: base case", "vertices": n, "edges": m}],
        "vertices": n,
        "edges": m,
    }


def _timing(formula, *, algorithm="dijkstra-binary", sizes=SIZES, fallback=None, **kw):
    """SYNTHETIC release records: sssp time = formula(n, m) with m = 2n."""
    recs = []
    for n in sizes:
        m = 2 * n
        for trial, j in enumerate(JITTER):
            recs.append(
                synthetic_record(
                    algorithm=algorithm,
                    mode=_mode(algorithm),
                    n=n,
                    m=m,
                    reachable=n,
                    trial=trial,
                    times_ns=_times(max(1, round(formula(n, m) * j))),
                    fallback=fallback(n, m) if fallback else None,
                    **kw,
                )
            )
    return recs


def _diag(formula, *, algorithm="dijkstra-binary", sizes=SIZES, key="relaxations", **kw):
    """SYNTHETIC release-diag records: ops[key] = formula(n, m); the time is a huge decoy."""
    recs = []
    for n in sizes:
        m = 2 * n
        recs.append(
            synthetic_diag_record(
                algorithm=algorithm,
                mode=_mode(algorithm),
                n=n,
                m=m,
                reachable=n,
                times_ns=_times(10**12),
                ops={key: round(formula(n, m))},
                **kw,
            )
        )
    return recs


@pytest.fixture
def dirs(tmp_path):
    return tmp_path / "raw", tmp_path / "summary", tmp_path / "figures"


def _run(dirs, *extra):
    raw, summary, figures = dirs
    return scaling_fit.main(
        [
            "--raw-dir",
            str(raw),
            "--summary-dir",
            str(summary),
            "--figures-dir",
            str(figures),
            "--bootstrap-resamples",
            "200",
            *extra,
        ]
    )


def _fit(dirs, records, *extra):
    write_jsonl(dirs[0] / "synthetic.jsonl", records)
    assert _run(dirs, "--no-figures", *extra) == 0
    return pd.read_csv(dirs[1] / "scaling_fit.csv")


def _row(df, quantity, model, **where):
    sel = df[(df["quantity"] == quantity) & (df["model"] == model)]
    for k, v in where.items():
        sel = sel[sel[k] == v]
    assert len(sel) == 1, f"expected one row for {quantity}/{model}/{where}, got {len(sel)}"
    return sel.iloc[0]


def _files(d):
    return sorted(p.name for p in d.rglob("*") if p.is_file()) if d.exists() else []


# ---------------------------------------------------------------- model terms
def test_model_list_is_the_documented_one():
    assert list(scaling_fit.MODELS) == [
        "m",
        "m_plus_n_logloglog_n",
        "m_plus_n_log_n",
        "m_log_n",
        "power_m",
    ]


def test_logloglog_is_exact_on_towers_and_undefined_below_16():
    assert scaling_fit.logloglog2(16) == 1.0
    assert scaling_fit.logloglog2(65536) == 2.0
    assert scaling_fit.logloglog2(256) == math.log2(3)
    assert scaling_fit.LOGLOGLOG_MIN_N == 16
    with pytest.raises(ValueError, match="n >= 16"):
        scaling_fit.logloglog2(15)


def test_model_terms_use_exact_log2():
    n, m = 1 << 20, 3 << 20
    assert scaling_fit.MODELS["m"](n, m) == m
    assert scaling_fit.MODELS["m_log_n"](n, m) == m * 20
    assert scaling_fit.MODELS["m_plus_n_log_n"](n, m) == m + n * 20
    assert scaling_fit.MODELS["m_plus_n_logloglog_n"](n, m) == m + n * _lll(n)


# ---------------------------------------------------------------- recovery of known laws
def test_recovers_linear_coefficient(dirs):
    df = _fit(dirs, _timing(lambda n, m: 7 * m))
    lin = _row(df, "time.sssp", "m")
    assert lin["status"] == "ok"
    assert lin["coefficient"] == pytest.approx(7, rel=1e-6)
    assert lin["rel_rmse"] == pytest.approx(0, abs=1e-6)
    assert lin["coefficient_ci_low"] <= lin["coefficient"] <= lin["coefficient_ci_high"]
    assert lin["sizes"] == 4 and lin["n_min"] == SIZES[0] and lin["n_max"] == SIZES[-1]
    assert lin["samples"] == 12
    # the wrong models leave a visible residual
    assert _row(df, "time.sssp", "m_log_n")["rel_rmse"] > 0.05
    power = _row(df, "time.sssp", "power_m")
    assert power["alpha"] == pytest.approx(1.0, abs=1e-6)
    assert power["coefficient"] == pytest.approx(7, rel=1e-4)


def test_recovers_m_plus_n_logloglog_n(dirs):
    df = _fit(dirs, _timing(lambda n, m: 3000 * (m + n * _lll(n))))
    t = df[(df["quantity"] == "time.sssp") & (df["status"] == "ok")]
    best = t.loc[t[t["model"] != "power_m"]["rel_rmse"].idxmin()]
    assert best["model"] == "m_plus_n_logloglog_n"
    assert best["coefficient"] == pytest.approx(3000, rel=1e-6)
    assert best["rel_rmse"] < 1e-6


def test_recovers_power_law_exponent_with_bootstrap_interval(dirs):
    df = _fit(dirs, _timing(lambda n, m: 2000 * m**1.25))
    power = _row(df, "time.sssp", "power_m")
    assert power["alpha"] == pytest.approx(1.25, abs=1e-6)
    assert power["coefficient"] == pytest.approx(2000, rel=1e-4)
    assert power["alpha_ci_low"] <= 1.25 + 1e-9 and power["alpha_ci_high"] >= 1.25 - 1e-9
    assert power["alpha_ci_high"] - power["alpha_ci_low"] < 0.05  # 2% jitter only
    assert power["bootstrap_resamples"] == 200
    # one-parameter models have no exponent
    assert pd.isna(_row(df, "time.sssp", "m")["alpha"])


def test_fit_is_reproducible_with_the_fixed_seed(dirs, tmp_path):
    recs = _timing(lambda n, m: 2000 * m**1.25)
    first = _fit(dirs, recs)
    other = (tmp_path / "raw2", tmp_path / "summary2", tmp_path / "figures2")
    second = _fit(other, recs)
    pd.testing.assert_frame_equal(first, second)


# ---------------------------------------------------------------- timing vs counters
def test_counters_come_from_diag_records_and_times_only_from_release(dirs):
    recs = _timing(lambda n, m: 7 * m)
    recs += _diag(lambda n, m: 5 * m * math.log2(n))  # times in these records are 1e12 decoys
    df = _fit(dirs, recs)
    ops = _row(df, "ops.relaxations", "m_log_n")
    assert ops["status"] == "ok"
    assert ops["coefficient"] == pytest.approx(5, rel=1e-9)
    assert ops["rel_rmse"] == pytest.approx(0, abs=1e-9)
    # diagnostics-ON timings did not enter the time fit
    lin = _row(df, "time.sssp", "m")
    assert lin["coefficient"] == pytest.approx(7, rel=1e-6)
    assert lin["samples"] == 12


def test_diag_only_data_gives_no_time_fit(dirs):
    df = _fit(dirs, _diag(lambda n, m: 5 * m))
    t = df[df["quantity"] == "time.sssp"]
    assert len(t) == 1
    assert t.iloc[0]["status"] == "refused"
    assert "release" in t.iloc[0]["reason"]
    assert pd.isna(t.iloc[0]["coefficient"])
    assert _row(df, "ops.relaxations", "m")["coefficient"] == pytest.approx(5)


def test_warmup_and_failed_records_are_not_used(dirs):
    recs = _timing(lambda n, m: 7 * m)
    recs.append(synthetic_record(n=SIZES[0], m=2 * SIZES[0], warmup=True, times_ns=_times(10**12)))
    recs.append(
        synthetic_record(
            n=SIZES[1],
            m=None,
            source=None,
            trial=None,
            reachable=None,
            status="timeout",
            error="timeout after 1 s",
            times_ns=dict.fromkeys(_times(1)),
            peak_rss_bytes=None,
            rss_after_graph_bytes=None,
            correctness={"verified": False, "ok": None, "method": None},
        )
    )
    lin = _row(_fit(dirs, recs), "time.sssp", "m")
    assert lin["coefficient"] == pytest.approx(7, rel=1e-6)
    assert lin["samples"] == 12


# ---------------------------------------------------------------- refusals
def _only_refusal(df, quantity="time.sssp"):
    sel = df[df["quantity"] == quantity]
    assert len(sel) == 1
    row = sel.iloc[0]
    assert row["status"] == "refused"
    assert pd.isna(row["model"]) and pd.isna(row["coefficient"]) and pd.isna(row["alpha"])
    return row


def test_refuses_fewer_than_four_sizes(dirs, capsys):
    df = _fit(dirs, _timing(lambda n, m: 7 * m, sizes=(1 << 10, 1 << 14, 1 << 17)))
    row = _only_refusal(df)
    assert "3 sizes" in row["reason"] and "4" in row["reason"]
    assert "3 sizes" in capsys.readouterr().out  # the reason is printed, not only stored


def test_refuses_less_than_two_orders_of_magnitude(dirs):
    df = _fit(dirs, _timing(lambda n, m: 7 * m, sizes=(1000, 2000, 4000, 8000, 64000)))
    row = _only_refusal(df)
    assert "orders of magnitude" in row["reason"]
    # exactly 100x is two orders of magnitude and is accepted
    ok = (dirs[0].parent / "r2", dirs[0].parent / "s2", dirs[0].parent / "f2")
    df = _fit(ok, _timing(lambda n, m: 7 * m, sizes=(1000, 5000, 20000, 100000)))
    assert _row(df, "time.sssp", "m")["status"] == "ok"


def test_refuses_paper_mode_without_fallback_share(dirs):
    df = _fit(dirs, _timing(lambda n, m: 7 * m, algorithm="paper-faithful"))
    row = _only_refusal(df)
    assert "fallback" in row["reason"]


def test_fallback_share_is_reported_next_to_every_fit(dirs):
    recs = _timing(lambda n, m: 7 * m, algorithm="paper-faithful", fallback=_all_fallback)
    recs += _timing(lambda n, m: 9 * m)  # baseline: no fallback concept, still fitted
    df = _fit(dirs, recs)
    paper = df[(df["algorithm"] == "paper-faithful") & (df["quantity"] == "time.sssp")]
    assert len(paper) == len(scaling_fit.MODELS)
    assert (paper["status"] == "ok").all()
    assert (paper["fallback_vertex_share_min"] == 1.0).all()
    assert (paper["fallback_vertex_share_max"] == 1.0).all()
    assert (paper["fallback_edge_share_max"] == 1.0).all()
    assert (paper["fallback_source"] == "release").all()
    base = _row(df, "time.sssp", "m", algorithm="dijkstra-binary")
    assert base["coefficient"] == pytest.approx(9, rel=1e-6)
    assert pd.isna(base["fallback_vertex_share_min"])
    assert base["fallback_source"] == "not applicable (baseline)"


def test_fallback_share_may_come_from_the_diag_run(dirs):
    recs = _timing(lambda n, m: 7 * m, algorithm="paper-faithful")  # fallback null in release
    recs += _diag(
        lambda n, m: 5 * m,
        algorithm="paper-faithful",
        fallback={
            "count": 1,
            "records": [{"reason": "outer: base case", "vertices": 256, "edges": 512}],
            "vertices": 256,
            "edges": 512,
        },
    )
    df = _fit(dirs, recs)
    lin = _row(df, "time.sssp", "m")
    assert lin["status"] == "ok"
    assert lin["fallback_source"] == "release-diag"
    assert lin["fallback_vertex_share_max"] == pytest.approx(256 / SIZES[0])
    assert lin["fallback_vertex_share_min"] == pytest.approx(256 / SIZES[-1])


def test_refuses_mixed_machines(dirs):
    recs = _timing(lambda n, m: 7 * m)
    other_cpu = dict(recs[0]["cpu"], brand="SYNTHETIC-TEST-DATA other machine")
    recs += _timing(lambda n, m: 7 * m, cpu=other_cpu)
    row = _only_refusal(_fit(dirs, recs))
    assert "machine" in row["reason"]


# ---------------------------------------------------------------- boundary n < 16
def test_sizes_below_16_are_excluded_and_counted(dirs):
    df = _fit(dirs, _timing(lambda n, m: 7 * m, sizes=(8, 15, *SIZES)))
    t = df[df["quantity"] == "time.sssp"]
    assert (t["status"] == "ok").all()
    assert (t["excluded_n_lt_16"] == 2).all()
    assert (t["sizes"] == 4).all() and (t["n_min"] == SIZES[0]).all()
    assert _row(df, "time.sssp", "m")["coefficient"] == pytest.approx(7, rel=1e-6)


def test_n_equal_16_is_kept(dirs):
    df = _fit(dirs, _timing(lambda n, m: 7000 * m, sizes=(16, 256, 4096, 65536)))
    lin = _row(df, "time.sssp", "m")
    assert lin["status"] == "ok" and lin["n_min"] == 16 and lin["excluded_n_lt_16"] == 0
    assert _row(df, "time.sssp", "m_plus_n_logloglog_n")["status"] == "ok"


def test_exclusion_can_leave_too_few_sizes(dirs):
    df = _fit(dirs, _timing(lambda n, m: 7 * m, sizes=(8, 1 << 10, 1 << 14, 1 << 17)))
    row = _only_refusal(df)
    assert "3 sizes" in row["reason"]
    assert row["excluded_n_lt_16"] == 1


# ---------------------------------------------------------------- outputs
def test_no_data_writes_nothing(dirs, capsys):
    assert _run(dirs) == 0
    assert "データなし" in capsys.readouterr().out
    assert _files(dirs[1]) == [] and _files(dirs[2]) == []


def test_schema_violations_fail_unless_skipped(dirs, capsys):
    bad = synthetic_record(n=0)
    write_jsonl(dirs[0] / "synthetic.jsonl", [*_timing(lambda n, m: 7 * m), bad])
    assert _run(dirs, "--no-figures") == 1
    assert _files(dirs[1]) == []
    assert _run(dirs, "--no-figures", "--skip-invalid") == 0
    assert _files(dirs[1]) == ["scaling_fit.csv"]
    capsys.readouterr()


def test_only_measured_sizes_are_listed(dirs):
    lin = _row(_fit(dirs, _timing(lambda n, m: 7 * m)), "time.sssp", "m")
    assert [int(x) for x in str(lin["n_values"]).split(";")] == list(SIZES)


def test_dirty_and_throttled_records_are_flagged(dirs):
    recs = _timing(lambda n, m: 7 * m)
    recs[0]["git_dirty"] = True
    recs[1]["power"] = {
        "source": "AC",
        "battery_percent": None,
        "low_power_mode": False,
        "cpu_speed_limit": 80,
    }
    lin = _row(_fit(dirs, recs), "time.sssp", "m")
    assert lin["dirty_records"] == 1
    assert lin["throttled_records"] == 1


def test_figures_are_optional_and_named_scaling(dirs):
    write_jsonl(dirs[0] / "synthetic.jsonl", _timing(lambda n, m: 7 * m))
    assert _run(dirs) == 0
    names = _files(dirs[2])
    assert names and all(n.startswith("scaling_") and n.endswith(".png") for n in names)


def test_refused_groups_get_no_figure(dirs):
    write_jsonl(dirs[0] / "synthetic.jsonl", _timing(lambda n, m: 7 * m, sizes=SIZES[:3]))
    assert _run(dirs) == 0
    assert _files(dirs[2]) == []


def test_missing_diag_run_is_a_refusal_without_a_fallback_source(dirs):
    recs = _timing(lambda n, m: 7 * m, algorithm="paper-faithful", fallback=_all_fallback)
    row = _only_refusal(_fit(dirs, recs), quantity="ops")
    assert "release-diag" in row["reason"]
    assert row["fallback_source"] == "not recorded"
    assert row["sizes"] == 0 and row["samples"] == 0
