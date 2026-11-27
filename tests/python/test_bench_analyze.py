"""scripts/analyze_results.py on SYNTHETIC JSONL.

All records here are fabricated test inputs (see bench_fixtures.py). The tests check that the
pipeline aggregates and plots exactly what is in the files: nothing interpolated, nothing
invented when there is no data.
"""

from __future__ import annotations

import math

import analyze_results
import pandas as pd
import pytest
from bench import schema
from bench_fixtures import synthetic_diag_record, synthetic_record, write_jsonl


def _times(sssp, preprocess=None):
    e2e = sssp + (preprocess or 0)
    return {
        "graph_gen": 1000,
        "preprocess": preprocess,
        "sssp": sssp,
        "tree": None,
        "end_to_end": e2e,
        "verify": 10,
    }


def _rec(algorithm, n, trial, sssp, **kw):
    mode = "baseline" if algorithm.startswith("dijkstra-") else "paper-faithful"
    return synthetic_record(
        algorithm=algorithm,
        mode=mode,
        n=n,
        m=2 * n,
        reachable=n,
        trial=trial,
        times_ns=_times(sssp, kw.pop("preprocess", None)),
        **kw,
    )


@pytest.fixture
def dirs(tmp_path):
    raw = tmp_path / "raw"
    summary = tmp_path / "summary"
    figures = tmp_path / "figures"
    return raw, summary, figures


def _run(dirs, *extra):
    raw, summary, figures = dirs
    return analyze_results.main(
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


def _files(d):
    return sorted(p.name for p in d.rglob("*") if p.is_file()) if d.exists() else []


def test_no_data_produces_nothing(dirs, capsys):
    assert _run(dirs) == 0
    assert "データなし" in capsys.readouterr().out
    assert _files(dirs[1]) == []
    assert _files(dirs[2]) == []


def test_empty_files_count_as_no_data(dirs, capsys):
    raw = dirs[0]
    raw.mkdir()
    (raw / "empty.jsonl").write_text("", encoding="utf-8")
    assert _run(dirs) == 0
    assert "データなし" in capsys.readouterr().out
    assert _files(dirs[1]) == []


def _dataset():
    recs = []
    # SYNTHETIC: baseline measured at n=1000 and n=10000; paper only at n=1000.
    for t, v in enumerate([100, 110, 120]):
        recs.append(_rec("dijkstra-binary", 1000, t, v))
        recs.append(_rec("dijkstra-binary", 10000, t, v * 10))
    for t, v in enumerate([50, 55, 60]):
        recs.append(
            _rec("paper", 1000, t, v, preprocess=20, params={"threshold": 4}, peak_rss_bytes=2**24)
        )
    for t, v in enumerate([70, 80, 90]):
        recs.append(_rec("paper", 1000, t, v, preprocess=20, params={"threshold": 8}))
    # a fixed-degree family (grid): avg_degree is null, not a number
    for t, v in enumerate([300, 310]):
        recs.append(_rec("dijkstra-binary", 1024, t, v, graph_type="grid", avg_degree=None))
        recs.append(_rec("paper", 1024, t, v // 2, graph_type="grid", avg_degree=None))
    # warm-up and diagnostics-build timings must not enter the timing statistics
    recs.append(_rec("dijkstra-binary", 1000, 0, 10**9, warmup=True))
    recs.append(
        synthetic_diag_record(
            algorithm="paper",
            mode="paper-faithful",
            params={"threshold": 4},
            times_ns=_times(10**9, 10),
            fallback={
                "count": 1,
                "records": [{"reason": "outer: base case", "vertices": 100, "edges": 150}],
                "vertices": 100,
                "edges": 150,
            },
        )
    )
    # a timeout is a status, not a timing
    recs.append(
        synthetic_record(
            algorithm="paper",
            mode="paper-faithful",
            n=10000,
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
    return recs


def test_summary_csvs_aggregate_only_measured_data(dirs):
    raw, summary, figures = dirs
    write_jsonl(raw / "synthetic.jsonl", _dataset())
    assert _run(dirs) == 0

    timing = pd.read_csv(summary / "timing_summary.csv")
    sssp = timing[timing["metric"] == "sssp"]
    base = sssp[(sssp["algorithm"] == "dijkstra-binary") & (sssp["n"] == 1000)]
    assert len(base) == 1
    assert base.iloc[0]["median"] == 110  # warm-up record (1e9) excluded
    assert base.iloc[0]["samples"] == 3
    assert bool(base.iloc[0]["p95_reference_only"]) is True
    paper = sssp[(sssp["algorithm"] == "paper") & (sssp["graph_type"] == "er")]
    assert set(paper["n"]) == {1000}  # n=10000 timed out: no timing row, no interpolation
    assert paper["median"].max() < 10**9  # diagnostics-build timing excluded

    sp = pd.read_csv(summary / "speedup.csv")
    sp = sp[sp["metric"] == "sssp"]
    row = sp[(sp["algorithm"] == "paper") & (sp["params"].str.contains('"threshold": 4'))]
    assert len(row) == 1
    assert row.iloc[0]["speedup"] == pytest.approx(110 / 55)
    assert row.iloc[0]["ci_low"] <= row.iloc[0]["speedup"] <= row.iloc[0]["ci_high"]
    assert set(sp[(sp["algorithm"] == "paper") & (sp["graph_type"] == "er")]["n"]) == {1000}

    status = pd.read_csv(summary / "status.csv")
    to = status[(status["algorithm"] == "paper") & (status["n"] == 10000)]
    assert to.iloc[0]["timeout"] == 1

    fb = pd.read_csv(summary / "fallback.csv")
    assert fb.iloc[0]["fallback_rate"] == 1.0
    assert fb.iloc[0]["fallback_vertex_share"] == pytest.approx(100 / 1000)

    ops = pd.read_csv(summary / "ops.csv")
    assert "relaxations" in ops.columns

    grid = sp[sp["graph_type"] == "grid"]
    assert len(grid) == 1 and grid.iloc[0]["speedup"] == pytest.approx(305 / 152.5)
    assert grid["avg_degree"].isna().all()

    mem = pd.read_csv(summary / "memory.csv")
    assert set(mem["algorithm"]) >= {"paper", "dijkstra-binary"}


def test_figures_are_written(dirs):
    raw, summary, figures = dirs
    write_jsonl(raw / "synthetic.jsonl", _dataset())
    assert _run(dirs) == 0
    names = _files(figures)
    for prefix in (
        "runtime_vs_n",
        "speedup_vs_n",
        "runtime_vs_m",
        "memory_vs_n",
        "speedup_by_graph_type",
        "speedup_by_weight_dist",
        "breakdown",
        "fallback_rate",
        "param_",
        "runtime_vs_n__grid",
    ):
        assert any(n.startswith(prefix) and n.endswith(".png") for n in names), (prefix, names)


def test_gaps_are_not_interpolated():
    ys = analyze_results.series_with_gaps([1, 2, 3], {1: 10.0, 3: 30.0})
    assert ys[0] == 10.0 and ys[2] == 30.0
    assert math.isnan(ys[1])


def test_invalid_records_fail_unless_skipped(dirs, capsys):
    raw = dirs[0]
    write_jsonl(raw / "bad.jsonl", [synthetic_record(), synthetic_record(mode="fast")])
    assert _run(dirs) == 1
    assert "mode" in capsys.readouterr().err
    assert _run(dirs, "--skip-invalid") == 0


def test_fixed_degree_rows_land_in_a_facet():
    df = pd.DataFrame(
        {
            "weight_dist": ["uniform_int", "uniform_int"],
            "weight_type": ["uint64", "uint64"],
            "avg_degree": [float("nan"), 4.0],
            "y": [1.0, 2.0],
        }
    )
    cells = analyze_results.facet_cells(df)
    assert sorted(len(c) for c in cells.values()) == [1, 1]
    assert ("uniform_int / uint64", "fixed degree") in cells


# Names in bench.schema that get no colour of their own: the Bellman-Ford oracle, and "paper",
# which is reserved for the outer construction and not produced by any runner yet, and the
# ninth name "practical-local" (the palette has eight hues; names beyond fold into OTHER_COLOR).
_UNCOLOURED = {"bellman-ford", "paper", "practical-local"}


def test_registered_algorithms_have_distinct_colours():
    names = [a for a in schema.ALGORITHM_MODES if a not in _UNCOLOURED]
    assert {"paper-faithful", "paper-faithful-d1"} <= set(names)
    colours = {a: analyze_results.ALGO_COLORS.get(a, analyze_results.OTHER_COLOR) for a in names}
    # nobody falls back to the "other" colour, and identity is never shared
    assert analyze_results.OTHER_COLOR not in colours.values(), colours
    assert len(set(colours.values())) == len(colours), colours
    style = analyze_results._SeriesStyle([("paper-faithful", "{}"), ("paper-faithful-d1", "{}")])
    assert style("paper-faithful", "{}")["color"] != style("paper-faithful-d1", "{}")["color"]
