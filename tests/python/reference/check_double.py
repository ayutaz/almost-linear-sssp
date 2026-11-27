"""Verify a C++ double SSSP result exactly, with rational arithmetic.

Result dump (JSON)::

    {"n": n, "source": s, "dist": [...], "parent": [...], "parent_edge": [...]}

``dist`` entries are doubles printed with ``%.17g`` (which round-trips every
double) or ``null`` / ``"inf"`` / ``Infinity`` for unreachable vertices.  Each
number is read as the *exact binary value* of that double.  ``parent`` (vertex,
``-1`` = none) and ``parent_edge`` (index into the edge list, ``-1`` = none)
are optional.

Tolerance (identical to the C++ checker): u = 2^-53, k = hop bound (default
max(n-1, 0)), gamma_k = k*u / (1 - k*u) with k*u < 1.  Summing k+1
non-negative terms in floating point in any order gives |s_hat - s| <=
gamma_k * s (Higham, *Accuracy and Stability of Numerical Algorithms*, around
eq. (4.4)).  Every vertex must satisfy |dist(v) - d*(v)| <= rel*d*(v) + abs
with rel = gamma_k and abs = 0 by default.  All arithmetic is exact.

By default the graph weights are the doubles the C++ code actually used,
``Fraction(float(w))``; ``weights="exact"`` uses the exact decimal/rational
values instead.

CLI::

    uv run python tests/python/reference/check_double.py GRAPH.txt RESULT.json \
        [--rel R --abs A --hops K --weights double|exact]

exit 0 = accepted, 1 = rejected (reasons printed), 2 = unreadable input.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any

if __package__:
    from . import exact_sssp
else:  # executed as a script
    import exact_sssp

UNIT_ROUNDOFF = Fraction(1, 2**53)
WEIGHT_MODES = ("double", "exact")
_INF_STRINGS = {"inf", "+inf", "infinity", "+infinity"}


def gamma(k: int) -> Fraction:
    """gamma_k = k*u / (1 - k*u); requires 0 <= k and k*u < 1."""
    if k < 0:
        raise ValueError(f"hop bound must be non-negative, got {k}")
    ku = k * UNIT_ROUNDOFF
    if ku >= 1:
        raise ValueError(f"gamma_k needs k*u < 1, got k={k}")
    return ku / (1 - ku)


@dataclass
class Report:
    failures: list[str] = field(default_factory=list)
    rel: Fraction = Fraction(0)
    abs_tol: Fraction = Fraction(0)
    hops: int = 0

    @property
    def ok(self) -> bool:
        return not self.failures


# ---------------------------------------------------------------------------
# reading results
# ---------------------------------------------------------------------------


def _json_constant(name: str) -> float:
    if name == "Infinity":
        return math.inf
    raise ValueError(f"invalid distance constant {name}")


def _to_distance(x: Any) -> Fraction | None:
    """Exact value of one dumped distance; ``None`` means unreachable."""
    if x is None:
        return None
    if isinstance(x, bool):
        raise ValueError(f"invalid distance {x!r}")
    if isinstance(x, str):
        if x.strip().lower() in _INF_STRINGS:
            return None
        raise ValueError(f"invalid distance {x!r}")
    if isinstance(x, float):
        if math.isnan(x) or x == -math.inf:
            raise ValueError(f"invalid distance {x!r}")
        return None if x == math.inf else Fraction(x)
    if isinstance(x, int | Fraction):
        return Fraction(x)
    raise ValueError(f"invalid distance {x!r}")


def parse_result_json(text: str) -> dict:
    """Parse a dump; numbers become the exact value of the printed double."""
    result = json.loads(
        text,
        parse_float=lambda s: Fraction(float(s)),
        parse_constant=_json_constant,
    )
    if not isinstance(result, dict):
        raise ValueError("result JSON must be an object")
    dist = result.get("dist")
    if not isinstance(dist, list):
        raise ValueError("result JSON needs a 'dist' array")
    result["dist"] = [_to_distance(d) for d in dist]
    return result


def read_result(path: str | Path) -> dict:
    return parse_result_json(Path(path).read_text())


# ---------------------------------------------------------------------------
# checking
# ---------------------------------------------------------------------------


def _fmt(x: Fraction | None) -> str:
    return "unreachable" if x is None else f"{float(x):.17g}"


def _int_list(result: Mapping, key: str, n: int, lo: int, hi: int, failures: list[str]):
    """Validate an optional integer array; return it or None."""
    if key not in result or result[key] is None:
        return None
    values = result[key]
    if not isinstance(values, list) or len(values) != n:
        failures.append(f"{key}: expected a list of length {n}")
        return None
    bad = [
        i
        for i, x in enumerate(values)
        if isinstance(x, bool) or not isinstance(x, int) or not lo <= x < hi
    ]
    if bad:
        failures += [f"vertex {i}: {key} {values[i]!r} out of range [{lo}, {hi})" for i in bad]
        return None
    return values


def _double_weight(w: Fraction) -> Fraction:
    try:
        return Fraction(float(w))
    except OverflowError as exc:
        raise ValueError(f"edge weight {w} overflows a double") from exc


def check(
    n: int,
    edges: Sequence[tuple[int, int, Fraction]],
    result: Mapping,
    *,
    rel: Any = None,
    abs_tol: Any = 0,
    hops: int | None = None,
    weights: str = "double",
) -> Report:
    """Check a double SSSP result against the exact rational shortest paths."""
    if weights not in WEIGHT_MODES:
        raise ValueError(f"weights must be one of {WEIGHT_MODES}, got {weights!r}")
    k = max(n - 1, 0) if hops is None else hops
    rel_f = gamma(k) if rel is None else Fraction(rel)
    abs_f = Fraction(abs_tol)
    if rel_f < 0 or abs_f < 0:
        raise ValueError("tolerances must be non-negative")
    report = Report(rel=rel_f, abs_tol=abs_f, hops=k)
    fail = report.failures

    used = [
        (u, v, _double_weight(Fraction(w)) if weights == "double" else Fraction(w))
        for u, v, w in edges
    ]

    # --- shape --------------------------------------------------------------
    if result.get("n") != n:
        fail.append(f"n: result says {result.get('n')!r}, graph has {n}")
    source = result.get("source")
    if isinstance(source, bool) or not isinstance(source, int) or not 0 <= source < n:
        fail.append(f"source {source!r} out of range [0, {n})")
    raw = result.get("dist")
    if not isinstance(raw, list) or len(raw) != n:
        fail.append(f"dist: expected a list of length {n}")
    if fail:
        return report

    dist: list[Fraction | None] = []
    for v, x in enumerate(raw):
        try:
            dist.append(_to_distance(x))
        except ValueError as exc:
            fail.append(f"vertex {v}: {exc}")
            dist.append(None)
    exact = exact_sssp.dijkstra(n, used, source)

    # --- distances ----------------------------------------------------------
    if dist[source] != 0:
        fail.append(f"vertex {source}: source distance {_fmt(dist[source])} != 0")
    for v in range(n):
        dc, ds = dist[v], exact[v]
        if (dc is None) != (ds is None):
            fail.append(f"vertex {v}: reachability differs (result {_fmt(dc)}, exact {_fmt(ds)})")
        elif ds is not None and abs(dc - ds) > rel_f * ds + abs_f:
            fail.append(
                f"vertex {v}: distance {_fmt(dc)} vs exact {_fmt(ds)}, "
                f"error {_fmt(abs(dc - ds))} > bound {_fmt(rel_f * ds + abs_f)}"
            )

    # --- parent tree --------------------------------------------------------
    parent = _int_list(result, "parent", n, -1, n, fail)
    if parent is None:
        return report
    parent_edge = _int_list(result, "parent_edge", n, -1, len(used), fail)
    if "parent_edge" in result and result["parent_edge"] is not None and parent_edge is None:
        return report
    _check_tree(n, used, source, exact, parent, parent_edge, rel_f, abs_f, fail)
    return report


def _check_tree(n, used, source, exact, parent, parent_edge, rel_f, abs_f, fail) -> None:
    lightest: dict[tuple[int, int], Fraction] = {}
    for u, v, w in used:
        key = (min(u, v), max(u, v))
        if key not in lightest or w < lightest[key]:
            lightest[key] = w

    # weight of the edge v -> parent[v], or None if it is not a valid edge
    tree_w: list[Fraction | None] = [None] * n
    for v in range(n):
        p = parent[v]
        if v == source or exact[v] is None:
            if p != -1:
                what = "source" if v == source else "unreachable vertex"
                fail.append(f"vertex {v}: {what} has parent {p}")
            continue
        if p == -1:
            fail.append(f"vertex {v}: reachable but has no parent")
            continue
        if parent_edge is not None:
            e = parent_edge[v]
            if e == -1:
                fail.append(f"vertex {v}: parent {p} given but parent_edge is -1")
                continue
            a, b, w = used[e]
            if sorted((a, b)) != sorted((p, v)):
                fail.append(f"vertex {v}: parent_edge {e} is ({a}, {b}), not ({p}, {v})")
                continue
            tree_w[v] = w
        else:
            w = lightest.get((min(p, v), max(p, v)))
            if w is None:
                fail.append(f"vertex {v}: no edge between parent {p} and {v}")
                continue
            tree_w[v] = w

    # follow parent pointers; status True = reaches source through valid edges
    status: list[bool | None] = [None] * n
    length: list[Fraction | None] = [None] * n
    status[source] = True
    length[source] = Fraction(0)
    for v in range(n):
        if status[v] is not None or exact[v] is None:
            continue
        path: list[int] = []
        on_path: set[int] = set()
        x = v
        while True:
            if status[x] is not None:
                outcome = status[x]
                break
            if x in on_path:
                cycle = path[path.index(x) :]
                fail.append(f"vertex {v}: parent pointers form a cycle {cycle}")
                outcome = False
                break
            on_path.add(x)
            path.append(x)
            if parent[x] == -1 or tree_w[x] is None:
                outcome = False  # already reported above
                break
            x = parent[x]
        for y in reversed(path):
            status[y] = outcome
            if outcome:
                length[y] = length[parent[y]] + tree_w[y]

    # dist_c(v) is a rounded sum of the tree path, so (1 - rel) L <= dist_c(v), and the distance
    # check gives dist_c(v) <= (1 + rel) d* + abs; hence L <= ((1 + rel) d* + abs) / (1 - rel).
    if rel_f >= 1:
        fail.append(f"relative tolerance {rel_f} must be < 1 for the tree check")
        return
    for v in range(n):
        L, ds = length[v], exact[v]
        if L is None or ds is None:
            continue
        upper = (ds * (1 + rel_f) + abs_f) / (1 - rel_f)
        if not ds <= L <= upper:
            fail.append(
                f"vertex {v}: parent tree path length {_fmt(L)} outside [{_fmt(ds)}, {_fmt(upper)}]"
            )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

MAX_PRINTED = 50


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Exactly verify a double SSSP result.")
    parser.add_argument("graph", type=Path, help="edge-list text file")
    parser.add_argument("result", type=Path, help="result JSON dump")
    parser.add_argument("--rel", type=Fraction, default=None, help="relative tolerance")
    parser.add_argument("--abs", type=Fraction, default=Fraction(0), help="absolute tolerance")
    parser.add_argument("--hops", type=int, default=None, help="hop bound k (default n-1)")
    parser.add_argument("--weights", choices=WEIGHT_MODES, default="double")
    args = parser.parse_args(argv)
    try:
        n, edges = exact_sssp.read_graph(args.graph)
        result = read_result(args.result)
        report = check(
            n,
            edges,
            result,
            rel=args.rel,
            abs_tol=args.abs,
            hops=args.hops,
            weights=args.weights,
        )
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    summary = (
        f"n={n} hops={report.hops} rel={float(report.rel):.3g} abs={float(report.abs_tol):.3g}"
    )
    if report.ok:
        print(f"OK {summary}")
        return 0
    print(f"FAIL {summary}: {len(report.failures)} problem(s)")
    for line in report.failures[:MAX_PRINTED]:
        print(f"  {line}")
    if len(report.failures) > MAX_PRINTED:
        print(f"  ... {len(report.failures) - MAX_PRINTED} more")
    return 1


if __name__ == "__main__":
    sys.exit(main())
