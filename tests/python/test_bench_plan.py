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


PAPER = ("paper-faithful", "paper-faithful-d1")


def test_mode_a_entry_is_planned_under_its_registered_names():
    assert set(PAPER) <= set(plan.BENCH_ALGORITHMS)
    assert "paper" not in plan.BENCH_ALGORITHMS  # reserved schema name, no implementation
    jobs = plan.build_jobs("smoke", order_seed=1)
    for a in PAPER:
        mine = [j for j in jobs if j.algorithm == a]
        # smoke sizes are within every limit: both sizes, both graph types, both builds
        assert {j.n for j in mine} == set(plan.PROFILES["smoke"].sizes), a
        assert {j.graph_type for j in mine} == set(plan.PROFILES["smoke"].graph_types), a
        assert {j.build for j in mine} == set(plan.BUILDS), a


def test_vertex_limits_exclude_infeasible_jobs_with_a_reason():
    limit, reason = plan.VERTEX_LIMITS["paper-faithful"]
    assert reason  # the exclusion is explained in the plan itself
    assert set(plan.VERTEX_LIMITS) == {"paper-faithful"}  # the heap variant has no limit
    for profile in plan.PROFILES:
        jobs = plan.build_jobs(profile, order_seed=1)
        scan = {j.n for j in jobs if j.algorithm == "paper-faithful"}
        heap = {j.n for j in jobs if j.algorithm == "paper-faithful-d1"}
        assert scan == {n for n in plan.PROFILES[profile].sizes if n <= limit}, profile
        assert heap == set(plan.PROFILES[profile].sizes), profile
        excluded = plan.excluded_jobs(profile, order_seed=1)
        assert bool(excluded) == (max(plan.PROFILES[profile].sizes) > limit), profile
        for job, why in excluded:
            assert job.algorithm == "paper-faithful" and job.n > limit
            assert why == reason
    assert plan.excluded_jobs("smoke", order_seed=1) == []


def test_cli_reports_excluded_jobs_on_stderr():
    res = subprocess.run(
        [sys.executable, "-m", "bench.plan", "standard", "--order-seed", "3"],
        cwd=SCRIPTS,
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stderr
    n_excluded = len(plan.excluded_jobs("standard", order_seed=3))
    assert n_excluded > 0
    assert f"excluded {n_excluded} job(s)" in res.stderr
    assert "paper-faithful" in res.stderr
    planned = [line.split("\t") for line in res.stdout.splitlines()]
    algo, n = plan.Job._fields.index("algorithm"), plan.Job._fields.index("n")
    limit = plan.VERTEX_LIMITS["paper-faithful"][0]
    assert not [f for f in planned if f[algo] == "paper-faithful" and int(f[n]) > limit]
    assert [f for f in planned if f[algo] == "paper-faithful-d1" and int(f[n]) > limit]
