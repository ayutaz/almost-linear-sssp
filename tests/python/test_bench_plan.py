"""Benchmark job plan (scripts/bench/plan.py): one runner process per (build, graph, algorithm)."""

from __future__ import annotations

import subprocess
import sys

import pytest
from bench import plan, schema
from bench_fixtures import SCRIPTS


def test_profiles_match_protocol():
    smoke, std, full = plan.PROFILES["smoke"], plan.PROFILES["standard"], plan.PROFILES["full"]
    assert smoke.sizes == (1_000, 10_000)
    assert smoke.degrees == (4,)
    assert smoke.trials == 1
    assert std.trials == 5
    assert full.trials >= 20
    assert full.sizes == (1_000, 10_000, 100_000, 1_000_000, 10_000_000)
    assert full.degrees == (2, 4, 8, 16, 32)
    assert set(full.graph_types) == set(plan.GRAPH_TYPES)
    assert set(full.weight_dists) == set(plan.WEIGHT_DISTS)


@pytest.mark.parametrize("profile", ["smoke", "standard", "full"])
def test_jobs_are_single_algorithm_and_well_formed(profile):
    jobs = plan.build_jobs(profile, order_seed=1)
    assert jobs
    for j in jobs:
        assert j.build in ("release", "release-diag")
        assert j.algorithm in schema.ALGORITHM_MODES
        assert j.weight_type == plan.WEIGHT_TYPE_OF[j.weight_dist]
        if j.algorithm == "dijkstra-radix":
            assert j.weight_type == "uint64"
        # road / adversarial generators fix their own weights (docs/design/graph_library.md)
        assert (j.weight_dist == "builtin") == (j.graph_type in plan.BUILTIN_WEIGHTS)
    assert {j.build for j in jobs} == {"release", "release-diag"}
    assert "dijkstra-binary" in {j.algorithm for j in jobs}  # Speedup denominator always present


def test_graph_types_are_sssp_gen_types():
    # names must be accepted by `sssp_gen --type` (docs/design/graph_library.md §6)
    assert set(plan.GRAPH_TYPES) == {
        "er",
        "regular",
        "ba",
        "grid",
        "road",
        "tree",
        "chung_lu",
        "adversarial_long_path",
    }
    assert set(plan.WEIGHT_DISTS) == {
        "uniform_int",
        "uniform_real",
        "log_uniform",
        "exponential",
        "constant",
        "zero_mixed",
        "wide_range",
        "builtin",
    }


def test_fixed_degree_families_are_not_duplicated_per_degree():
    jobs = plan.build_jobs("full", order_seed=1)
    for gt in ("grid", "road", "tree", "adversarial_long_path"):
        fixed = {(j.avg_degree, j.shape_args) for j in jobs if j.graph_type == gt}
        assert fixed == {("-", "-")}, gt


def test_shape_arguments_follow_the_generator_cli():
    jobs = plan.build_jobs("full", order_seed=1, sizes=[1000], builds=["release"])
    shape = {(j.graph_type, j.avg_degree): j.shape_args for j in jobs}
    assert shape[("er", "8")] == "--avg-degree 8"
    assert shape[("chung_lu", "8")] == "--avg-degree 8"
    assert shape[("regular", "8")] == "--degree 8"
    assert shape[("ba", "8")] == "--attach 4"  # BA with m attachments has average degree ~2m


def test_order_is_shuffled_reproducibly_without_changing_the_set():
    a = plan.build_jobs("standard", order_seed=1)
    b = plan.build_jobs("standard", order_seed=1)
    c = plan.build_jobs("standard", order_seed=2)
    assert a == b
    assert a != c
    assert sorted(a) == sorted(c)


def test_filters_restrict_the_plan():
    jobs = plan.build_jobs(
        "smoke",
        order_seed=1,
        algorithms=["dijkstra-binary"],
        graph_types=["er"],
        sizes=[1000],
        builds=["release"],
    )
    assert len(jobs) == 1
    assert jobs[0].algorithm == "dijkstra-binary"
    with pytest.raises(ValueError):
        plan.build_jobs("smoke", order_seed=1, algorithms=["no-such-algo"])
    with pytest.raises(ValueError):
        plan.build_jobs("nope", order_seed=1)


def test_diag_build_uses_fewer_trials():
    jobs = plan.build_jobs("full", order_seed=1)
    rel = {j.trials for j in jobs if j.build == "release"}
    diag = {j.trials for j in jobs if j.build == "release-diag"}
    assert max(diag) <= min(rel)


def test_cli_prints_tab_separated_jobs():
    res = subprocess.run(
        [sys.executable, "-m", "bench.plan", "smoke", "--order-seed", "3"],
        cwd=SCRIPTS,
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stderr
    lines = res.stdout.strip().splitlines()
    assert len(lines) == len(plan.build_jobs("smoke", order_seed=3))
    fields = lines[0].split("\t")
    assert len(fields) == len(plan.Job._fields)
