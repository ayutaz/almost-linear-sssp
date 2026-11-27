"""Benchmark statistics (scripts/bench/stats.py). Inputs are SYNTHETIC numbers."""

from __future__ import annotations

import math

import pytest
from bench import stats
from hypothesis import given
from hypothesis import strategies as st


def test_summary_of_small_sample():
    s = stats.summarize([5, 1, 4, 2, 3])
    assert s["n"] == 5
    assert s["median"] == 3
    assert s["mean"] == 3
    assert s["min"] == 1
    assert s["max"] == 5
    assert s["stddev"] == pytest.approx(math.sqrt(2.5))  # sample stddev, ddof=1
    assert s["p95"] == pytest.approx(4.8)  # linear interpolation between order statistics
    assert s["p95_reference_only"] is True


def test_p95_is_reference_only_below_twenty_samples():
    assert stats.P95_MIN_SAMPLES == 20
    assert stats.summarize(list(range(19)))["p95_reference_only"] is True
    assert stats.summarize(list(range(20)))["p95_reference_only"] is False


def test_single_sample_has_no_stddev():
    s = stats.summarize([7])
    assert s["stddev"] is None
    assert s["median"] == 7


@pytest.mark.parametrize("bad", [[], [1, None], [1, float("nan")], [1, float("inf")]])
def test_summary_rejects_empty_or_missing_values(bad):
    with pytest.raises(ValueError):
        stats.summarize(bad)


def test_speedup_is_baseline_over_target_medians():
    assert stats.speedup([10, 10, 12], [5, 6, 5]) == pytest.approx(10 / 5)
    with pytest.raises(ValueError):
        stats.speedup([], [1])


def test_bootstrap_ci_is_deterministic_for_fixed_seed():
    base = [100, 110, 120, 105, 115]
    tgt = [50, 55, 60, 52, 58]
    a = stats.bootstrap_speedup_ci(base, tgt, n_resamples=2000, seed=42)
    b = stats.bootstrap_speedup_ci(base, tgt, n_resamples=2000, seed=42)
    assert a == b
    lo, hi = a
    assert lo <= stats.speedup(base, tgt) <= hi


def test_bootstrap_ci_of_constant_data_is_degenerate():
    assert stats.bootstrap_speedup_ci([4, 4, 4], [2, 2, 2], seed=0) == (2.0, 2.0)


@pytest.mark.parametrize("conf", [0, 1, 1.5, -0.1])
def test_bootstrap_rejects_bad_confidence(conf):
    with pytest.raises(ValueError):
        stats.bootstrap_speedup_ci([1, 2], [1, 2], confidence=conf)


def test_bootstrap_default_seed_is_fixed():
    base, tgt = [3, 5, 9, 4], [2, 3, 1, 4]
    assert stats.bootstrap_speedup_ci(base, tgt) == stats.bootstrap_speedup_ci(base, tgt)
    assert stats.DEFAULT_BOOTSTRAP_SEED is not None


finite = st.floats(min_value=1, max_value=1e12, allow_nan=False, allow_infinity=False)


@given(st.lists(finite, min_size=1, max_size=60))
def test_summary_order_invariants(values):
    s = stats.summarize(values)
    assert s["min"] <= s["median"] <= s["max"]
    assert s["median"] <= s["p95"] * (1 + 1e-12) + 1e-9
    assert s["p95"] <= s["max"] * (1 + 1e-12)
    assert s["p95_reference_only"] is (len(values) < 20)


@given(st.lists(finite, min_size=1, max_size=15), st.lists(finite, min_size=1, max_size=15))
def test_bootstrap_interval_is_ordered(base, tgt):
    lo, hi = stats.bootstrap_speedup_ci(base, tgt, n_resamples=200, seed=1)
    assert lo <= hi
