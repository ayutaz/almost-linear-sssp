#!/usr/bin/env python3
"""Fit scaling models to measured benchmark records (P6-SCALING).

Reads results/raw/*.jsonl (schema: docs/design/benchmark_schema.md), validates every record,
and writes
  results/summary/scaling_fit.csv
  docs/reports/figures/scaling__<quantity>__<graph_type>.png   (optional, --no-figures)
With no data it prints "データなし" and writes nothing.

Record selection (same rules as scripts/analyze_results.py):
  * times   : status == ok, warmup == false, diagnostics_enabled == false (release build)
  * counters: status == ok, warmup == false, diagnostics_enabled == true  (release-diag build)

One series = (graph_type, avg_degree, weight_dist, weight_type, algorithm, mode, params).
One point  = one size n of the series: x = (n, mean m), y = median over the records of that n.
Quantities: time.sssp, time.end_to_end, ops.<OpCounters key>, counter.<Diagnostics counter>.

Models, y ~ c * f(n, m) (log = exact log2, never a quotient of natural logarithms):
  m                      f = m
  m_plus_n_logloglog_n   f = m + n * log2(log2(log2(n)))     (defined for n >= 16)
  m_plus_n_log_n         f = m + n * log2(n)
  m_log_n                f = m * log2(n)
  power_m                y ~ c * m^alpha
Sizes with n < 16 are dropped from every model of the series (so that all models are compared
on the same points) and counted in ``excluded_n_lt_16``.

Fitting is least squares on log y: c = exp(mean(log y - log f)); the power law is a straight
line in (log m, log y). ``rel_rmse`` is sqrt(mean(((fit - y) / y)^2)) and ``log_rmse`` the RMS
of the natural-log residuals (the minimised quantity). Confidence intervals are a 95 %
percentile bootstrap with a fixed seed that resamples the records *within* each size and
refits; they describe measurement noise, not model error, and collapse to the point estimate
when a size has one sample (deterministic counters).

A series is NOT fitted, and the reason is printed and stored (status = refused), when
  * it has no usable records,
  * its time records come from more than one machine or build,
  * it has fewer than 4 sizes (after the n < 16 exclusion),
  * n spans less than 2 orders of magnitude (n_max < 100 * n_min),
  * it is not a baseline and the fallback share is not recorded for every size,
  * a size has m <= 0 or a non-positive median (no multiplicative model exists).
Fits are evaluated at the measured sizes only; nothing is interpolated or extrapolated.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from bench import schema, stats

ROOT = Path(__file__).resolve().parents[1]
SERIES = ("graph_type", "avg_degree", "weight_dist", "weight_type", "algorithm", "mode", "params")
TIME_QUANTITIES = ("sssp", "end_to_end")
MIN_SIZES = 4
MIN_SPAN = 100  # n_max / n_min: two orders of magnitude
LOGLOGLOG_MIN_N = 16
CONFIDENCE = 0.95
NOT_APPLICABLE = "not applicable (baseline)"

INK, INK2, GRID, AXIS = "#191919", "#666666", "#e7e7e7", "#c2c2c2"
MODEL_STYLE = {
    "m": ("#257de0", "o"),
    "m_plus_n_logloglog_n": ("#f47230", "s"),
    "m_plus_n_log_n": ("#00b894", "^"),
    "m_log_n": ("#6a4ce2", "D"),
    "power_m": ("#666666", "v"),
}


# ---------------------------------------------------------------- model terms
def logloglog2(n: float) -> float:
    """log2(log2(log2(n))), defined here for n >= 16 (where it is >= 1)."""
    if n < LOGLOGLOG_MIN_N:
        raise ValueError(f"log log log n is only used for n >= {LOGLOGLOG_MIN_N}, got {n}")
    return math.log2(math.log2(math.log2(n)))


#: name -> f(n, m), or None for the two-parameter power law c * m^alpha.
MODELS: dict[str, Callable[[float, float], float] | None] = {
    "m": lambda n, m: m,
    "m_plus_n_logloglog_n": lambda n, m: m + n * logloglog2(n),
    "m_plus_n_log_n": lambda n, m: m + n * math.log2(n),
    "m_log_n": lambda n, m: m * math.log2(n),
    "power_m": None,
}


# ---------------------------------------------------------------- loading
def load_records(raw_dir: Path) -> tuple[list[dict[str, Any]], list[str]]:
    records: list[dict[str, Any]] = []
    errors: list[str] = []
    for path in sorted(raw_dir.glob("*.jsonl")) if raw_dir.is_dir() else []:
        recs, errs = schema.load_jsonl(path)
        records += recs
        errors += errs
    return records, errors


def _series_key(rec: dict[str, Any]) -> tuple[Any, ...]:
    deg = rec["avg_degree"]
    return (
        rec["graph_type"],
        None if deg is None else float(deg),
        rec["weight_dist"],
        rec["weight_type"],
        rec["algorithm"],
        rec["mode"],
        json.dumps(rec["params"], sort_keys=True),
    )


def _sort_key(key: tuple[Any, ...]) -> tuple[Any, ...]:
    return tuple((v is None, 0.0 if v is None else v) for v in key)


def _usable(rec: dict[str, Any], *, diagnostics: bool) -> bool:
    return rec["status"] == "ok" and not rec["warmup"] and rec["diagnostics_enabled"] is diagnostics


def _quantities(recs: Sequence[dict[str, Any]], *, diagnostics: bool) -> list[str]:
    if not diagnostics:
        return [f"time.{k}" for k in TIME_QUANTITIES]
    ops = {k for r in recs for k in (r["ops"] or {})}
    counters = {k for r in recs for k in (r["counters"] or {})}
    return [f"ops.{k}" for k in schema.OPS_KEYS if k in ops] + [
        f"counter.{k}" for k in sorted(counters)
    ]


def _value(rec: dict[str, Any], quantity: str) -> float | None:
    kind, _, name = quantity.partition(".")
    block = {"time": rec["times_ns"], "ops": rec["ops"], "counter": rec["counters"]}[kind]
    v = (block or {}).get(name)
    return None if v is None else float(v)


# ---------------------------------------------------------------- fallback share
def _shares(recs: Sequence[dict[str, Any]]) -> dict[int, tuple[float, float | None]]:
    """n -> (mean vertex share, mean edge share) over the records that carry a fallback block."""
    by_n: dict[int, list[tuple[float, float | None]]] = {}
    for r in recs:
        fb = r["fallback"]
        if fb is None:
            continue
        edge = fb["edges"] / r["m"] if r["m"] > 0 else None
        by_n.setdefault(r["n"], []).append((fb["vertices"] / r["n"], edge))
    out: dict[int, tuple[float, float | None]] = {}
    for n, pairs in by_n.items():
        edges = [e for _, e in pairs if e is not None]
        out[n] = (
            float(np.mean([v for v, _ in pairs])),
            float(np.mean(edges)) if edges else None,
        )
    return out


@dataclass
class Fallback:
    source: str
    shares: dict[int, tuple[float, float | None]]
    missing: list[int]

    def summary(self) -> str:
        if not self.shares:
            return f"fallback: {self.source}"
        vs = [v for v, _ in self.shares.values()]
        return f"fallback vertex share {min(vs):.3g}..{max(vs):.3g} ({self.source})"

    def columns(self) -> dict[str, Any]:
        vs = [v for v, _ in self.shares.values()]
        es = [e for _, e in self.shares.values() if e is not None]
        return {
            "fallback_source": self.source,
            "fallback_vertex_share_min": min(vs) if vs else None,
            "fallback_vertex_share_max": max(vs) if vs else None,
            "fallback_edge_share_min": min(es) if es else None,
            "fallback_edge_share_max": max(es) if es else None,
            "fallback_vertex_share_by_n": ";".join(
                f"{n}:{self.shares[n][0]:.6g}" for n in sorted(self.shares)
            )
            or None,
        }


def fallback_for(
    mode: str,
    sizes: Sequence[int],
    used: Sequence[dict[str, Any]],
    diag: Sequence[dict[str, Any]],
    *,
    used_is_diag: bool,
) -> Fallback:
    """Fallback share of the fitted sizes: from the fitted records, else from the diag run."""
    if mode == "baseline":
        return Fallback(NOT_APPLICABLE, {}, [])
    own = _shares(used)
    candidates = [("release-diag" if used_is_diag else "release", own)]
    if not used_is_diag:
        candidates.append(("release-diag", _shares(diag)))
    for source, shares in candidates:
        if sizes and all(n in shares for n in sizes):
            return Fallback(source, {n: shares[n] for n in sizes}, [])
    # Sizes no source covers; if every size is covered by some source but no single source
    # covers them all, name the sizes the fitted records themselves lack.
    covered = set().union(*(shares for _, shares in candidates))
    missing = [n for n in sizes if n not in covered] or [n for n in sizes if n not in own]
    return Fallback("not recorded", {}, missing)


# ---------------------------------------------------------------- fitting
@dataclass
class ModelFit:
    model: str
    coefficient: float
    alpha: float | None
    rel_rmse: float
    log_rmse: float
    coefficient_ci: tuple[float, float]
    alpha_ci: tuple[float, float] | None
    fitted: np.ndarray  # model value at the measured sizes only


def _slope_intercept(lx: np.ndarray, ly: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Least-squares line ly ~ a * lx + b; ``ly`` may hold one series per row."""
    dx = lx - lx.mean()
    a = (dx * (ly - ly.mean(axis=-1, keepdims=True))).sum(axis=-1) / (dx * dx).sum()
    return a, ly.mean(axis=-1) - a * lx.mean()


def _interval(samples: np.ndarray, point: float) -> tuple[float, float]:
    tail = (1 - CONFIDENCE) / 2
    lo, hi = (float(v) for v in np.quantile(samples, [tail, 1 - tail]))
    # Floating-point noise must not push a degenerate interval off its point estimate.
    if math.isclose(lo, point, rel_tol=1e-12, abs_tol=1e-12):
        lo = min(lo, point)
    if math.isclose(hi, point, rel_tol=1e-12, abs_tol=1e-12):
        hi = max(hi, point)
    return lo, hi


def bootstrap_medians(samples: Sequence[np.ndarray], resamples: int, seed: int) -> np.ndarray:
    """(resamples, sizes) medians of records resampled with replacement within each size."""
    rng = np.random.default_rng(seed)
    cols = [
        np.median(s[rng.integers(0, s.size, size=(resamples, s.size))], axis=1) for s in samples
    ]
    return np.stack(cols, axis=1)


def fit_models(
    n: np.ndarray, m: np.ndarray, y: np.ndarray, boot: np.ndarray
) -> list[ModelFit | tuple[str, str]]:
    """Every model on the points (n, m, y); ``boot`` holds resampled y rows for the intervals.

    Returns a ModelFit, or (model, reason) when that model cannot be fitted to these points.
    """
    ly, lboot, lm = np.log(y), np.log(boot), np.log(m)
    out: list[ModelFit | tuple[str, str]] = []
    for name, term in MODELS.items():
        if term is None:
            if np.ptp(lm) == 0:
                out.append((name, "m does not vary across sizes: no exponent can be fitted"))
                continue
            a, b = _slope_intercept(lm, ly)
            a_b, b_b = _slope_intercept(lm, lboot)
            alpha, log_c = float(a), float(b)
            log_fit = alpha * lm + log_c
            alpha_ci = _interval(a_b, alpha)
            c_ci = _interval(np.exp(b_b), math.exp(log_c))
        else:
            lf = np.log([term(float(ni), float(mi)) for ni, mi in zip(n, m, strict=True)])
            alpha, alpha_ci = None, None
            log_c = float((ly - lf).mean())
            log_fit = lf + log_c
            c_ci = _interval(np.exp((lboot - lf).mean(axis=1)), math.exp(log_c))
        fitted = np.exp(log_fit)
        out.append(
            ModelFit(
                model=name,
                coefficient=math.exp(log_c),
                alpha=alpha,
                rel_rmse=float(np.sqrt(np.mean(((fitted - y) / y) ** 2))),
                log_rmse=float(np.sqrt(np.mean((log_fit - ly) ** 2))),
                coefficient_ci=c_ci,
                alpha_ci=alpha_ci,
                fitted=fitted,
            )
        )
    return out


# ---------------------------------------------------------------- one series, one quantity
def _machines(recs: Sequence[dict[str, Any]]) -> set[tuple[Any, ...]]:
    return {(r["cpu"]["brand"], r["arch"], r["ram_bytes"]) for r in recs}


def _builds(recs: Sequence[dict[str, Any]]) -> set[tuple[Any, ...]]:
    return {(r["compiler"], r["build_type"], r["build_flags"]) for r in recs}


def _throttled(rec: dict[str, Any]) -> bool:
    limit = (rec["power"] or {}).get("cpu_speed_limit")
    return limit is not None and limit < 100


@dataclass
class Series:
    """The points of one quantity of one series, after record selection."""

    n: np.ndarray
    m: np.ndarray
    y: np.ndarray
    samples: list[np.ndarray]
    fits: list[ModelFit]
    fallback: str  # one-line fallback share, shown next to the series in figures


def analyse(
    key: tuple[Any, ...],
    quantity: str,
    used: Sequence[dict[str, Any]],
    diag: Sequence[dict[str, Any]],
    *,
    resamples: int,
    seed: int,
) -> tuple[list[dict[str, Any]], Series | None]:
    """CSV rows for one quantity of one series, and the fitted points (None if refused)."""
    is_time = quantity.startswith("time.")
    mode = key[SERIES.index("mode")]
    recs = [r for r in used if _value(r, quantity) is not None]
    by_n: dict[int, list[dict[str, Any]]] = {}
    for r in recs:
        by_n.setdefault(r["n"], []).append(r)
    sizes = sorted(n for n in by_n if n >= LOGLOGLOG_MIN_N)
    kept = [r for n in sizes for r in by_n[n]]
    fallback = fallback_for(mode, sizes, kept, diag, used_is_diag=not is_time)
    common: dict[str, Any] = {
        **dict(zip(SERIES, key, strict=True)),
        "quantity": quantity,
        "sizes": len(sizes),
        "n_min": sizes[0] if sizes else None,
        "n_max": sizes[-1] if sizes else None,
        "n_values": ";".join(str(n) for n in sizes) or None,
        "excluded_n_lt_16": sum(1 for n in by_n if n < LOGLOGLOG_MIN_N),
        "samples": len(kept),
        "samples_per_size_min": min((len(by_n[n]) for n in sizes), default=None),
        **fallback.columns(),
        "bootstrap_resamples": resamples,
        "bootstrap_seed": seed,
        "dirty_records": sum(1 for r in kept if r["git_dirty"]),
        "throttled_records": sum(1 for r in kept if _throttled(r)),
        "run_ids": ";".join(sorted({r["run_id"] for r in kept})) or None,
        "git_shas": ";".join(sorted({r["git_sha"] for r in kept})) or None,
    }

    def refused(reason: str, model: str | None = None) -> dict[str, Any]:
        return {**common, "model": model, "status": "refused", "reason": reason}

    m = np.array([np.mean([r["m"] for r in by_n[n]]) for n in sizes], dtype=float)
    samples = [np.array([_value(r, quantity) for r in by_n[n]], dtype=float) for n in sizes]
    y = np.array([np.median(s) for s in samples], dtype=float)
    common["m_min"] = float(m.min()) if sizes else None
    common["m_max"] = float(m.max()) if sizes else None

    reason = None
    if not recs:
        build = "release (diagnostics off)" if is_time else "release-diag"
        reason = f"no {build} records with status ok and warmup false"
    elif is_time and len(_machines(kept)) > 1:
        reason = f"records come from {len(_machines(kept))} machines; never mix machines in a fit"
    elif is_time and len(_builds(kept)) > 1:
        reason = f"records come from {len(_builds(kept))} compilers/build flags"
    elif len(sizes) < MIN_SIZES:
        reason = f"{len(sizes)} sizes with n >= {LOGLOGLOG_MIN_N} (need >= {MIN_SIZES})"
    elif sizes[-1] < MIN_SPAN * sizes[0]:
        span = math.log10(sizes[-1] / sizes[0])
        reason = f"n spans {span:.2f} orders of magnitude, {sizes[0]}..{sizes[-1]} (need >= 2)"
    elif fallback.missing or (mode != "baseline" and not fallback.shares):
        missing = ",".join(str(n) for n in fallback.missing)
        reason = f"fallback share not recorded for mode {mode} at n={missing}"
    elif (m <= 0).any():
        reason = "a size has m = 0: no multiplicative model in m exists"
    elif (y <= 0).any():
        zero = ",".join(str(n) for n, v in zip(sizes, y, strict=True) if v <= 0)
        reason = f"median is not positive at n={zero}: no multiplicative model exists"
    if reason is not None:
        return [refused(reason)], None

    n_arr = np.array(sizes, dtype=float)
    boot = bootstrap_medians(samples, resamples, seed)
    if (boot <= 0).any():
        return [refused("a bootstrap resample has a non-positive median")], None
    rows: list[dict[str, Any]] = []
    fits: list[ModelFit] = []
    for fit in fit_models(n_arr, m, y, boot):
        if isinstance(fit, tuple):
            rows.append(refused(fit[1], model=fit[0]))
            continue
        fits.append(fit)
        rows.append(
            {
                **common,
                "model": fit.model,
                "status": "ok",
                "reason": None,
                "coefficient": fit.coefficient,
                "alpha": fit.alpha,
                "rel_rmse": fit.rel_rmse,
                "log_rmse": fit.log_rmse,
                "coefficient_ci_low": fit.coefficient_ci[0],
                "coefficient_ci_high": fit.coefficient_ci[1],
                "alpha_ci_low": fit.alpha_ci[0] if fit.alpha_ci else None,
                "alpha_ci_high": fit.alpha_ci[1] if fit.alpha_ci else None,
            }
        )
    return rows, Series(n_arr, m, y, samples, fits, fallback.summary())


COLUMNS = (
    *SERIES,
    "quantity",
    "model",
    "status",
    "reason",
    "coefficient",
    "alpha",
    "rel_rmse",
    "log_rmse",
    "coefficient_ci_low",
    "coefficient_ci_high",
    "alpha_ci_low",
    "alpha_ci_high",
    "sizes",
    "n_min",
    "n_max",
    "m_min",
    "m_max",
    "n_values",
    "excluded_n_lt_16",
    "samples",
    "samples_per_size_min",
    "fallback_source",
    "fallback_vertex_share_min",
    "fallback_vertex_share_max",
    "fallback_edge_share_min",
    "fallback_edge_share_max",
    "fallback_vertex_share_by_n",
    "bootstrap_resamples",
    "bootstrap_seed",
    "dirty_records",
    "throttled_records",
    "run_ids",
    "git_shas",
)


def fit_all(
    records: Sequence[dict[str, Any]], *, resamples: int, seed: int
) -> tuple[pd.DataFrame, dict[tuple[Any, ...], Series]]:
    """The scaling_fit table and, per (series key, quantity), the fitted points."""
    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for r in records:
        groups.setdefault(_series_key(r), []).append(r)
    rows: list[dict[str, Any]] = []
    fitted: dict[tuple[Any, ...], Series] = {}
    for key in sorted(groups, key=_sort_key):
        timing = [r for r in groups[key] if _usable(r, diagnostics=False)]
        diag = [r for r in groups[key] if _usable(r, diagnostics=True)]
        if not timing and not diag:
            continue  # only warm-up or failed records: nothing was measured for this series
        jobs = [(q, timing) for q in _quantities(timing, diagnostics=False)]
        jobs += [(q, diag) for q in _quantities(diag, diagnostics=True) or ["ops"]]
        for quantity, used in jobs:
            new, series = analyse(key, quantity, used, diag, resamples=resamples, seed=seed)
            rows += new
            if series is not None:
                fitted[(*key, quantity)] = series
    return pd.DataFrame(rows, columns=list(COLUMNS)), fitted


# ---------------------------------------------------------------- figures
def _series_label(key: tuple[Any, ...]) -> str:
    _, deg, wdist, wtype, algo, _, params = key
    parts = [algo if params == "{}" else f"{algo} {params}"]
    parts.append(f"{'fixed degree' if deg is None else f'd={deg:g}'}, {wdist} / {wtype}")
    return "\n".join(parts)


def make_figures(fitted: dict[tuple[Any, ...], Series], fig_dir: Path) -> list[Path]:
    """measured / fitted at the measured sizes, one panel per series (time.* and ops.* only)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    plain = FuncFormatter(lambda v, _pos: f"{v:g}")

    panels: dict[tuple[str, str], list[tuple[tuple[Any, ...], Series]]] = {}
    for (*key, quantity), series in fitted.items():
        if quantity.startswith(("time.", "ops.")):
            panels.setdefault((quantity, key[0]), []).append((tuple(key), series))
    written: list[Path] = []
    for (quantity, graph_type), items in sorted(panels.items()):
        cols = min(3, len(items))
        rows = math.ceil(len(items) / cols)
        fig, axes = plt.subplots(
            rows, cols, squeeze=False, figsize=(3.8 * cols + 1.6, 2.9 * rows + 0.9)
        )
        handles: dict[str, Any] = {}
        for ax, (key, s) in zip(axes.flat, items, strict=False):
            for fit in s.fits:
                color, marker = MODEL_STYLE[fit.model]
                label = fit.model if fit.alpha is None else f"{fit.model} (fitted alpha)"
                (h,) = ax.plot(s.n, s.y / fit.fitted, color=color, marker=marker, linewidth=1.5)
                handles[label] = h
            ax.axhline(1.0, color=AXIS, linewidth=1, linestyle="--")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(plain)
            ax.yaxis.set_minor_formatter(plain)
            ax.grid(True, which="major", color=GRID, linewidth=0.6)
            ax.set_axisbelow(True)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            ax.tick_params(which="both", colors=INK2, labelsize=8)
            ax.set_title(f"{_series_label(key)}\n{s.fallback}", fontsize=8, color=INK)
            ax.set_xlabel("vertices n (measured sizes only)", fontsize=8, color=INK2)
            ax.set_ylabel("measured / fitted", fontsize=8, color=INK2)
        for ax in list(axes.flat)[len(items) :]:
            ax.set_visible(False)
        fig.suptitle(f"scaling fit residuals: {quantity}, {graph_type}", fontsize=11, color=INK)
        fig.tight_layout(rect=(0, 0, 0.8, 0.95))
        fig.legend(
            handles.values(),
            handles.keys(),
            loc="center left",
            bbox_to_anchor=(0.81, 0.5),
            fontsize=7,
            frameon=False,
        )
        stem = re.sub(r"[^A-Za-z0-9_-]+", "_", f"{quantity}__{graph_type}")
        path = fig_dir / f"scaling__{stem}.png"
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        written.append(path)
    return written


# ---------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("--raw-dir", type=Path, default=ROOT / "results" / "raw")
    ap.add_argument("--summary-dir", type=Path, default=ROOT / "results" / "summary")
    ap.add_argument("--figures-dir", type=Path, default=ROOT / "docs" / "reports" / "figures")
    ap.add_argument(
        "--skip-invalid", action="store_true", help="drop invalid records instead of failing"
    )
    ap.add_argument("--bootstrap-resamples", type=int, default=stats.DEFAULT_RESAMPLES)
    ap.add_argument("--bootstrap-seed", type=int, default=stats.DEFAULT_BOOTSTRAP_SEED)
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)
    if args.bootstrap_resamples < 1:
        ap.error("--bootstrap-resamples must be >= 1")

    records, errors = load_records(args.raw_dir)
    if errors:
        for e in errors:
            print(e, file=sys.stderr)
        print(f"スキーマ違反 {len(errors)} 件", file=sys.stderr)
        if not args.skip_invalid:
            return 1
    table, fitted = fit_all(records, resamples=args.bootstrap_resamples, seed=args.bootstrap_seed)
    if table.empty:
        print(f"データなし: {args.raw_dir} に使えるレコードがありません。何も生成しません。")
        return 0

    for _, row in table[table["status"] == "refused"].iterrows():
        label = " ".join(str(row[c]) for c in ("algorithm", "params", "graph_type", "weight_type"))
        model = "" if row["model"] is None or pd.isna(row["model"]) else f" [{row['model']}]"
        print(f"当てはめない: {label} {row['quantity']}{model}: {row['reason']}")
    args.summary_dir.mkdir(parents=True, exist_ok=True)
    path = args.summary_dir / "scaling_fit.csv"
    table.to_csv(path, index=False)
    ok = int((table["status"] == "ok").sum())
    print(f"wrote {path} ({ok} fits, {len(table) - ok} refused)")
    if not args.no_figures:
        for fig in make_figures(fitted, args.figures_dir):
            print(f"wrote {fig}")
    for column, what in (("dirty_records", "dirty な作業ツリー"), ("throttled_records", "熱制限")):
        bad = table[(table["status"] == "ok") & (table[column] > 0)]
        if not bad.empty:
            print(f"注意: {what}のレコードを含む当てはめが {len(bad)} 行あります(公表値に使わない)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
