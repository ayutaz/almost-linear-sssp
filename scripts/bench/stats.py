"""Statistics for benchmark samples (docs/design/benchmark_protocol.md, "指標と統計").

- ``summarize``: n, median, mean, p95, stddev (sample, ddof=1), min, max. p95 uses linear
  interpolation between order statistics (numpy's default) and is flagged
  ``p95_reference_only`` when fewer than ``P95_MIN_SAMPLES`` samples exist.
- ``speedup``: median(baseline) / median(target)  (> 1 means the target is faster).
- ``bootstrap_speedup_ci``: percentile bootstrap of that ratio, resampling the two samples
  independently, with a fixed RNG seed so that reports are reproducible.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

P95_MIN_SAMPLES = 20
DEFAULT_BOOTSTRAP_SEED = 20261127
DEFAULT_RESAMPLES = 10_000


def _as_array(values: Sequence[float], what: str = "values") -> np.ndarray:
    if len(values) == 0:
        raise ValueError(f"{what}: no samples")
    if any(v is None for v in values):
        raise ValueError(f"{what}: contains null (unmeasured) values; filter them first")
    arr = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{what}: contains NaN or infinite values")
    return arr


def summarize(values: Sequence[float]) -> dict[str, float | int | bool | None]:
    arr = _as_array(values)
    n = int(arr.size)
    return {
        "n": n,
        "median": float(np.median(arr)),
        "mean": float(arr.mean()),
        "p95": float(np.percentile(arr, 95)),
        "stddev": float(arr.std(ddof=1)) if n > 1 else None,
        "min": float(arr.min()),
        "max": float(arr.max()),
        "p95_reference_only": n < P95_MIN_SAMPLES,
    }


def speedup(baseline: Sequence[float], target: Sequence[float]) -> float:
    base = _as_array(baseline, "baseline")
    tgt = _as_array(target, "target")
    t = float(np.median(tgt))
    if t <= 0:
        raise ValueError("target median must be positive")
    return float(np.median(base)) / t


def bootstrap_speedup_ci(
    baseline: Sequence[float],
    target: Sequence[float],
    *,
    n_resamples: int = DEFAULT_RESAMPLES,
    confidence: float = 0.95,
    seed: int = DEFAULT_BOOTSTRAP_SEED,
) -> tuple[float, float]:
    if not 0 < confidence < 1:
        raise ValueError(f"confidence must be in (0, 1), got {confidence}")
    if n_resamples < 1:
        raise ValueError("n_resamples must be >= 1")
    base = _as_array(baseline, "baseline")
    tgt = _as_array(target, "target")
    rng = np.random.default_rng(seed)
    b = np.median(rng.choice(base, size=(n_resamples, base.size), replace=True), axis=1)
    t = np.median(rng.choice(tgt, size=(n_resamples, tgt.size), replace=True), axis=1)
    ratios = b / t
    alpha = (1 - confidence) / 2
    lo, hi = np.quantile(ratios, [alpha, 1 - alpha])
    lo_f, hi_f = float(lo), float(hi)
    point = speedup(base, tgt)
    # Floating-point noise must not push a degenerate interval off its point estimate.
    if math.isclose(lo_f, point, rel_tol=1e-12):
        lo_f = min(lo_f, point)
    if math.isclose(hi_f, point, rel_tol=1e-12):
        hi_f = max(hi_f, point)
    return lo_f, hi_f
