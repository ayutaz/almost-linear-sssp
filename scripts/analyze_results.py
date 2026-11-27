#!/usr/bin/env python3
"""Aggregate benchmark JSONL into summary CSVs and figures (`make report`).

Reads results/raw/*.jsonl (schema: docs/design/benchmark_schema.md), validates every record,
and writes
  results/summary/{timing_summary,speedup,memory,status,fallback,ops}.csv
  docs/reports/figures/*.png
Only measured combinations are aggregated or drawn: nothing is interpolated or extrapolated,
and a series line is broken where a point is missing. With no data it prints "データなし" and
writes nothing.

Rules (docs/design/benchmark_protocol.md):
  * timing  : status == ok, warmup == false, diagnostics_enabled == false (release build)
  * counters: diagnostics_enabled == true (release-diag build)
  * Speedup : median(dijkstra-binary) / median(target) on the same graph configuration
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
from bench import schema, stats

ROOT = Path(__file__).resolve().parents[1]
BASELINE = "dijkstra-binary"
CONFIG = ["graph_type", "n", "avg_degree", "weight_dist", "weight_type"]
GROUP = [*CONFIG, "algorithm", "mode", "params"]
TIMING_METRICS = ("graph_gen", "preprocess", "sssp", "tree", "end_to_end", "verify")
SPEEDUP_METRICS = ("sssp", "end_to_end")
STATUS_COLS = ("ok", "timeout", "memory_limit", "error")

# Fixed colour per algorithm (identity follows the entity, never its rank); markers are the
# secondary encoding so that series are distinguishable without colour.
ALGO_COLORS = {
    "dijkstra-binary": "#257de0",
    "dijkstra-4ary": "#f47230",
    "dijkstra-fibonacci": "#00b894",
    "dijkstra-radix": "#f2a600",
    "paper": "#ee7fac",
    "paper-experimental": "#008a17",
    "paper-optimized": "#6a4ce2",
}
OTHER_COLOR = "#ec4647"
MARKERS = ("o", "s", "^", "D", "v", "P", "X", "*")
INK, INK2, GRID, AXIS = "#191919", "#666666", "#e7e7e7", "#c2c2c2"


# ---------------------------------------------------------------- loading
def load_records(raw_dir: Path) -> tuple[list[dict[str, Any]], list[str]]:
    records: list[dict[str, Any]] = []
    errors: list[str] = []
    for path in sorted(raw_dir.glob("*.jsonl")) if raw_dir.is_dir() else []:
        recs, errs = schema.load_jsonl(path)
        records += recs
        errors += errs
    return records, errors


def to_frame(records: Iterable[dict[str, Any]]) -> pd.DataFrame:
    rows = []
    for r in records:
        row = {
            k: r[k]
            for k in (
                *CONFIG,
                "algorithm",
                "mode",
                "m",
                "seed",
                "trial",
                "warmup",
                "status",
                "diagnostics_enabled",
                "git_dirty",
                "run_id",
                "isolated",
                "peak_rss_bytes",
                "rss_after_graph_bytes",
            )
        }
        row["params"] = json.dumps(r["params"], sort_keys=True)
        for k in schema.TIME_KEYS:
            row[f"t_{k}"] = r["times_ns"][k]
        for k, v in (r["ops"] or {}).items():
            row[f"ops.{k}"] = v
        for k, v in (r["counters"] or {}).items():
            row[f"counter.{k}"] = v
        fb = r["fallback"]
        row["fallback_count"] = fb["count"] if fb else None
        row["fallback_vertices"] = fb["vertices"] if fb else None
        row["fallback_edges"] = fb["edges"] if fb else None
        row["fallback_reasons"] = [x["reason"] for x in fb["records"]] if fb else None
        rows.append(row)
    df = pd.DataFrame(rows)
    df["avg_degree"] = df["avg_degree"].astype(float)
    return df


def _groups(df: pd.DataFrame, keys: Sequence[str]):
    return df.groupby(list(keys), dropna=False, sort=True)


def _values(series: pd.Series) -> list[float]:
    return [float(v) for v in series if v is not None and not pd.isna(v)]


# ---------------------------------------------------------------- summaries
def timing_frame(df: pd.DataFrame) -> pd.DataFrame:
    ok = (df["status"] == "ok") & (~df["warmup"]) & (~df["diagnostics_enabled"])
    return df[ok]


def timing_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, g in _groups(timing_frame(df), GROUP):
        for metric in TIMING_METRICS:
            vals = _values(g[f"t_{metric}"])
            if not vals:
                continue  # never measured: no row (not a zero)
            s = stats.summarize(vals)
            rows.append(
                {
                    **dict(zip(GROUP, key, strict=True)),
                    "metric": metric,
                    "samples": s["n"],
                    "median": s["median"],
                    "mean": s["mean"],
                    "p95": s["p95"],
                    "stddev": s["stddev"],
                    "min": s["min"],
                    "max": s["max"],
                    "p95_reference_only": s["p95_reference_only"],
                    "m_mean": float(g["m"].mean()),
                    "dirty_records": int(g["git_dirty"].sum()),
                    "run_ids": ";".join(sorted(set(g["run_id"]))),
                }
            )
    return pd.DataFrame(rows)


def speedup_table(df: pd.DataFrame, resamples: int, seed: int) -> pd.DataFrame:
    t = timing_frame(df)
    base = t[t["algorithm"] == BASELINE]
    rows = []
    for key, g in _groups(t, GROUP):
        cfg = dict(zip(GROUP, key, strict=True))
        if cfg["algorithm"] == BASELINE:
            continue
        mask = pd.Series(True, index=base.index)
        for c in CONFIG:
            mask &= base[c].isna() if pd.isna(cfg[c]) else base[c] == cfg[c]
        b = base[mask]
        for metric in SPEEDUP_METRICS:
            bv, tv = _values(b[f"t_{metric}"]), _values(g[f"t_{metric}"])
            if not bv or not tv:
                continue  # baseline or target not measured here: no Speedup
            lo, hi = stats.bootstrap_speedup_ci(bv, tv, n_resamples=resamples, seed=seed)
            rows.append(
                {
                    **cfg,
                    "metric": metric,
                    "speedup": stats.speedup(bv, tv),
                    "ci_low": lo,
                    "ci_high": hi,
                    "baseline_samples": len(bv),
                    "target_samples": len(tv),
                }
            )
    return pd.DataFrame(rows)


def memory_table(df: pd.DataFrame) -> pd.DataFrame:
    sel = timing_frame(df)
    sel = sel[sel["isolated"] & sel["peak_rss_bytes"].notna()]
    rows = []
    for key, g in _groups(sel, GROUP):
        peak = _values(g["peak_rss_bytes"])
        after = _values(g["rss_after_graph_bytes"])
        rows.append(
            {
                **dict(zip(GROUP, key, strict=True)),
                "samples": len(peak),
                "peak_rss_median": stats.summarize(peak)["median"],
                "peak_rss_max": max(peak),
                "rss_after_graph_median": stats.summarize(after)["median"] if after else None,
            }
        )
    return pd.DataFrame(rows)


def status_table(df: pd.DataFrame) -> pd.DataFrame:
    """Job-level outcome: a job (run, build, configuration, algorithm) fails if any record does."""
    sel = df[~df["warmup"]]
    rows = []
    for key, g in _groups(sel, [*GROUP, "diagnostics_enabled"]):
        counts = dict.fromkeys(STATUS_COLS, 0)
        for _, job in _groups(g, ["run_id", "seed"]):
            failed = [s for s in job["status"] if s != "ok"]
            counts[failed[0] if failed else "ok"] += 1
        total = sum(counts.values())
        rows.append(
            {
                **dict(zip([*GROUP, "diagnostics_enabled"], key, strict=True)),
                "jobs": total,
                **counts,
                "failure_rate": (total - counts["ok"]) / total,
            }
        )
    return pd.DataFrame(rows)


def fallback_table(df: pd.DataFrame) -> pd.DataFrame:
    sel = df[(df["status"] == "ok") & (~df["warmup"]) & df["diagnostics_enabled"]]
    sel = sel[sel["fallback_count"].notna()]
    rows = []
    for key, g in _groups(sel, GROUP):
        reasons: dict[str, int] = {}
        for lst in g["fallback_reasons"]:
            for r in lst:
                reasons[r] = reasons.get(r, 0) + 1
        m = g["m"].astype(float)
        rows.append(
            {
                **dict(zip(GROUP, key, strict=True)),
                "samples": len(g),
                "fallback_rate": float((g["fallback_count"] > 0).mean()),
                "fallback_count_mean": float(g["fallback_count"].mean()),
                "fallback_vertex_share": float((g["fallback_vertices"] / g["n"]).mean()),
                "fallback_edge_share": float((g["fallback_edges"] / m.where(m > 0)).mean())
                if (m > 0).any()
                else None,
                "reasons": "; ".join(f"{r} x{c}" for r, c in sorted(reasons.items())),
            }
        )
    return pd.DataFrame(rows)


def ops_table(df: pd.DataFrame) -> pd.DataFrame:
    sel = df[(df["status"] == "ok") & (~df["warmup"]) & df["diagnostics_enabled"]]
    cols = [c for c in sel.columns if c.startswith(("ops.", "counter."))]
    rows = []
    for key, g in _groups(sel, GROUP):
        row: dict[str, Any] = {**dict(zip(GROUP, key, strict=True)), "samples": len(g)}
        for c in cols:
            vals = _values(g[c])
            if vals:
                row[c.removeprefix("ops.")] = stats.summarize(vals)["median"]
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- plotting helpers
def series_with_gaps(x_all: Sequence[float], points: dict[float, float]) -> list[float]:
    """y for every x in ``x_all``; NaN where unmeasured so matplotlib breaks the line there."""
    return [points.get(x, math.nan) for x in x_all]


def _label(algorithm: str, params: str) -> str:
    return algorithm if params in ("{}", "") else f"{algorithm} {params}"


def _style(ax) -> None:
    ax.grid(True, which="major", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(which="both", colors=INK2, labelsize=8)


def _plain_ticks(ax) -> None:
    """Ratios read better as 0.5 / 1 / 2 than as powers of ten."""
    from matplotlib.ticker import FuncFormatter

    fmt = FuncFormatter(lambda v, _pos: f"{v:g}")
    ax.yaxis.set_major_formatter(fmt)
    ax.yaxis.set_minor_formatter(fmt)
    ax.tick_params(axis="y", which="minor", labelsize=7, labelcolor=INK2)


class _SeriesStyle:
    """Colour by algorithm, marker by (algorithm, params) variant -- stable across figures."""

    def __init__(self, labels: Iterable[tuple[str, str]]):
        self._marker: dict[tuple[str, str], str] = {}
        per_algo: dict[str, int] = {}
        for algo, params in sorted(set(labels)):
            i = per_algo.get(algo, 0)
            per_algo[algo] = i + 1
            self._marker[(algo, params)] = MARKERS[i % len(MARKERS)]

    def __call__(self, algo: str, params: str) -> dict[str, Any]:
        return {"color": ALGO_COLORS.get(algo, OTHER_COLOR), "marker": self._marker[(algo, params)]}


FIXED_DEGREE_LABEL = "fixed degree"


def _deg_label(v: Any) -> str:
    return FIXED_DEGREE_LABEL if v is None or pd.isna(v) else f"d={v:g}"


def _deg_order(label: str) -> tuple[bool, float]:
    return (label == FIXED_DEGREE_LABEL, 0.0 if label == FIXED_DEGREE_LABEL else float(label[2:]))


def facet_cells(data: pd.DataFrame) -> dict[tuple[str, str], pd.DataFrame]:
    """Rows by weight distribution/type, columns by nominal degree ("fixed degree" if null).

    Grouping on labels (not on the float column) keeps null-degree families: NaN != NaN.
    """
    row = data["weight_dist"] + " / " + data["weight_type"]
    col = data["avg_degree"].map(_deg_label)
    return {(str(r), str(c)): g for (r, c), g in data.groupby([row, col], sort=True)}


def _facet_lines(
    data: pd.DataFrame,
    *,
    x: str,
    y: str,
    title: str,
    xlabel: str,
    ylabel: str,
    path: Path,
    logy: bool = True,
    err: tuple[str, str] | None = None,
    hline: float | None = None,
) -> bool:
    import matplotlib.pyplot as plt

    data = data.dropna(subset=[x, y])
    if data.empty:
        return False
    cells = facet_cells(data)
    rows = sorted({r for r, _ in cells})
    cols = sorted({c for _, c in cells}, key=_deg_order)
    style = _SeriesStyle(zip(data["algorithm"], data["params"], strict=True))
    fig, axes = plt.subplots(
        len(rows),
        len(cols),
        squeeze=False,
        sharex=True,
        figsize=(3.6 * len(cols) + 1.5, 2.8 * len(rows) + 0.8),
    )
    handles: dict[str, Any] = {}
    for i, rname in enumerate(rows):
        for j, cname in enumerate(cols):
            ax = axes[i][j]
            _style(ax)
            cell = cells.get((rname, cname))
            ax.set_title(f"{rname}, {cname}", fontsize=9, color=INK)
            if cell is None:
                ax.text(
                    0.5,
                    0.5,
                    "not measured",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    color=INK2,
                    fontsize=8,
                )
                continue
            x_all = sorted(cell[x].unique())
            for (algo, params), s in cell.groupby(["algorithm", "params"], sort=True):
                pts = dict(zip(s[x], s[y], strict=True))
                ys = series_with_gaps(x_all, pts)
                st = style(algo, params)
                (h,) = ax.plot(x_all, ys, linewidth=2, markersize=6, **st)
                if err is not None:
                    lo = series_with_gaps(x_all, dict(zip(s[x], s[err[0]], strict=True)))
                    hi = series_with_gaps(x_all, dict(zip(s[x], s[err[1]], strict=True)))
                    ax.vlines(x_all, lo, hi, color=st["color"], linewidth=1)
                handles[_label(algo, params)] = h
            ax.set_xscale("log")
            if logy:
                ax.set_yscale("log")
            if hline is not None:
                ax.axhline(hline, color=AXIS, linewidth=1, linestyle="--")
                _plain_ticks(ax)
            if i == len(rows) - 1:
                ax.set_xlabel(xlabel, fontsize=8, color=INK2)
            if j == 0:
                ax.set_ylabel(ylabel, fontsize=8, color=INK2)
    fig.suptitle(title, fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 0.78, 0.95))
    fig.legend(
        handles.values(),
        handles.keys(),
        loc="center left",
        bbox_to_anchor=(0.79, 0.5),
        fontsize=7,
        frameon=False,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return True


def _strip(sp: pd.DataFrame, by: str, path: Path) -> bool:
    import matplotlib.pyplot as plt

    sp = sp[sp["metric"] == "sssp"]
    if sp.empty:
        return False
    cats = sorted(sp[by].unique())
    series = sorted(set(zip(sp["algorithm"], sp["params"], strict=True)))
    style = _SeriesStyle(series)
    fig, ax = plt.subplots(figsize=(max(4.0, 1.2 * len(cats) + 2), 3.6))
    _style(ax)
    width = 0.7 / max(1, len(series))
    for k, (algo, params) in enumerate(series):
        s = sp[(sp["algorithm"] == algo) & (sp["params"] == params)]
        xs = [cats.index(c) - 0.35 + width * (k + 0.5) for c in s[by]]
        ax.scatter(
            xs,
            s["speedup"],
            s=36,
            label=_label(algo, params),
            edgecolors="white",
            linewidths=1,
            **style(algo, params),
        )
    ax.axhline(1.0, color=AXIS, linewidth=1, linestyle="--")
    ax.set_xticks(range(len(cats)), cats)
    ax.set_yscale("log")
    _plain_ticks(ax)
    ax.set_ylabel(f"Speedup vs {BASELINE} (sssp, median)", fontsize=8, color=INK2)
    ax.set_title(
        f"Speedup by {by} (each point = one measured configuration)", fontsize=10, color=INK
    )
    ax.legend(fontsize=7, frameon=False, loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return True


def _breakdown(ts: pd.DataFrame, gt: str, path: Path) -> bool:
    import matplotlib.pyplot as plt

    sub = ts[(ts["graph_type"] == gt) & ts["metric"].isin(["preprocess", "sssp", "tree"])]
    if sub.empty or not (sub["metric"] != "sssp").any():
        return False
    if (
        sub[["n", "avg_degree", "weight_dist", "algorithm", "params"]].drop_duplicates().shape[0]
        > 60
    ):
        sub = sub[sub["n"] == sub["n"].max()]  # keep the figure readable
    keys = sorted(
        set(
            zip(
                sub["n"],
                sub["avg_degree"],
                sub["weight_dist"],
                sub["algorithm"],
                sub["params"],
                strict=True,
            )
        )
    )
    fig, ax = plt.subplots(figsize=(max(5.0, 0.6 * len(keys) + 2), 4.4))
    _style(ax)
    colors = {"preprocess": "#257de0", "sssp": "#f47230", "tree": "#00b894"}
    bottoms = [0.0] * len(keys)
    for metric in ("preprocess", "sssp", "tree"):
        vals = []
        for n, d, w, a, p in keys:
            r = sub[
                (sub["n"] == n)
                & (sub["avg_degree"] == d)
                & (sub["weight_dist"] == w)
                & (sub["algorithm"] == a)
                & (sub["params"] == p)
                & (sub["metric"] == metric)
            ]
            vals.append(float(r["median"].iloc[0]) if len(r) else 0.0)  # absent = no segment
        ax.bar(
            range(len(keys)),
            vals,
            bottom=bottoms,
            width=0.6,
            color=colors[metric],
            label=metric,
            edgecolor="white",
            linewidth=2,
        )
        bottoms = [b + v for b, v in zip(bottoms, vals, strict=True)]
    labels = [f"{_label(a, p)}\nn={n} {_deg_label(d)} {w}" for n, d, w, a, p in keys]
    ax.set_xticks(range(len(keys)), labels, rotation=60, ha="right", fontsize=6)
    ax.set_xlim(-0.75, len(keys) - 0.25)
    ax.set_ylabel("median time [ns]", fontsize=8, color=INK2)
    ax.set_title(f"{gt}: time breakdown", fontsize=10, color=INK, loc="left", pad=14)
    ax.text(
        0,
        1.01,
        "stacked medians per part; unmeasured (null) parts are not drawn",
        transform=ax.transAxes,
        fontsize=7,
        color=INK2,
        va="bottom",
    )
    ax.legend(fontsize=7, frameon=False, loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return True


def _param_plots(ts: pd.DataFrame, fig_dir: Path) -> list[Path]:
    import matplotlib.pyplot as plt

    sub = ts[(ts["metric"] == "sssp") & (ts["params"] != "{}")].copy()
    out: list[Path] = []
    if sub.empty:
        return out
    sub["_p"] = sub["params"].map(json.loads)
    for algo, g in sub.groupby("algorithm"):
        keys = sorted(
            {
                k
                for p in g["_p"]
                for k, v in p.items()
                if isinstance(v, int | float) and not isinstance(v, bool)
            }
        )
        for key in keys:
            rows = g[g["_p"].map(lambda p, k=key: isinstance(p.get(k), int | float))]
            rows = rows.assign(_x=rows["_p"].map(lambda p, k=key: float(p[k])))
            if rows["_x"].nunique() < 2:
                continue
            fig, ax = plt.subplots(figsize=(5.5, 3.6))
            _style(ax)
            for idx, (cfg, s) in enumerate(rows.groupby(CONFIG, dropna=False, sort=True)):
                s = s.sort_values("_x")
                label = f"{cfg[0]} n={cfg[1]} {_deg_label(cfg[2])} {cfg[3]}"
                ax.plot(
                    s["_x"],
                    s["median"],
                    marker=MARKERS[idx % len(MARKERS)],
                    linewidth=2,
                    color=ALGO_COLORS.get(str(algo), OTHER_COLOR),
                    label=label,
                    alpha=1.0 if idx == 0 else 0.75,
                )
            ax.set_yscale("log")
            ax.set_xlabel(key, fontsize=8, color=INK2)
            ax.set_ylabel("median sssp time [ns]", fontsize=8, color=INK2)
            ax.set_title(f"{algo}: {key} vs time (measured values only)", fontsize=10, color=INK)
            ax.legend(fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1, 0.5))
            fig.tight_layout()
            path = fig_dir / f"param__{algo}__{key}.png"
            fig.savefig(path, dpi=120, bbox_inches="tight")
            plt.close(fig)
            out.append(path)
    return out


def make_figures(tables: dict[str, pd.DataFrame], fig_dir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    written: list[Path] = []
    ts, sp = tables["timing_summary"], tables["speedup"]
    mem, fb = tables["memory"], tables["fallback"]
    gtypes = sorted(set(ts["graph_type"]) if not ts.empty else set())
    for gt in gtypes:
        t = ts[(ts["graph_type"] == gt) & (ts["metric"] == "sssp")]
        jobs = [
            (t, "n", "median", "runtime_vs_n", "vertices n", "median sssp time [ns]", {}),
            (t, "m_mean", "median", "runtime_vs_m", "edges m", "median sssp time [ns]", {}),
        ]
        if not sp.empty:
            s = sp[(sp["graph_type"] == gt) & (sp["metric"] == "sssp")]
            jobs.append(
                (
                    s,
                    "n",
                    "speedup",
                    "speedup_vs_n",
                    "vertices n",
                    f"Speedup vs {BASELINE}",
                    {"err": ("ci_low", "ci_high"), "hline": 1.0},
                )
            )
        if not mem.empty:
            jobs.append(
                (
                    mem[mem["graph_type"] == gt],
                    "n",
                    "peak_rss_median",
                    "memory_vs_n",
                    "vertices n",
                    "median Peak RSS [bytes]",
                    {},
                )
            )
        if not fb.empty:
            jobs.append(
                (
                    fb[fb["graph_type"] == gt],
                    "n",
                    "fallback_rate",
                    "fallback_rate",
                    "vertices n",
                    "share of runs with a fallback",
                    {"logy": False},
                )
            )
        for data, x, y, stem, xl, yl, kw in jobs:
            path = fig_dir / f"{stem}__{gt}.png"
            if _facet_lines(
                data,
                x=x,
                y=y,
                title=f"{stem.replace('_', ' ')}: {gt}",
                xlabel=xl,
                ylabel=yl,
                path=path,
                **kw,
            ):
                written.append(path)
        path = fig_dir / f"breakdown__{gt}.png"
        if _breakdown(ts, gt, path):
            written.append(path)
    if not sp.empty:
        for by, stem in (
            ("graph_type", "speedup_by_graph_type"),
            ("weight_dist", "speedup_by_weight_dist"),
        ):
            path = fig_dir / f"{stem}.png"
            if _strip(sp, by, path):
                written.append(path)
    if not ts.empty:
        written += _param_plots(ts, fig_dir)
    return written


# ---------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
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

    records, errors = load_records(args.raw_dir)
    if errors:
        for e in errors:
            print(e, file=sys.stderr)
        print(f"スキーマ違反 {len(errors)} 件", file=sys.stderr)
        if not args.skip_invalid:
            return 1
    if not records:
        print(f"データなし: {args.raw_dir} に有効なレコードがありません。何も生成しません。")
        return 0

    df = to_frame(records)
    tables = {
        "timing_summary": timing_summary(df),
        "speedup": speedup_table(df, args.bootstrap_resamples, args.bootstrap_seed),
        "memory": memory_table(df),
        "status": status_table(df),
        "fallback": fallback_table(df),
        "ops": ops_table(df),
    }
    args.summary_dir.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        if table.empty:
            print(f"{name}: データなし(CSV は作成しません)")
            continue
        path = args.summary_dir / f"{name}.csv"
        table.to_csv(path, index=False)
        print(f"wrote {path} ({len(table)} rows)")
    if not args.no_figures:
        for path in make_figures(tables, args.figures_dir):
            print(f"wrote {path}")
    dirty = int(df["git_dirty"].sum())
    if dirty:
        print(f"注意: dirty な作業ツリーで得たレコードが {dirty} 件あります(公表値に使わない)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
