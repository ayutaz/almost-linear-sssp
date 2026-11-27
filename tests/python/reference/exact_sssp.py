"""Exact single-source shortest paths over rationals (``fractions.Fraction``).

This is the test oracle: it never rounds, so its distances are the true
shortest-path lengths of the input graph.

Edge-list text format::

    # comments and blank lines are ignored
    n m
    u v w        (m lines; 0-based vertices, undirected edge)

``w`` may be an integer (``3``), a fraction (``p/q``) or a decimal string
(``0.1``, ``1e-3``); it is read exactly with ``Fraction(str)``.  Multi-edges
and self-loops are allowed.

CLI::

    uv run python tests/python/reference/exact_sssp.py GRAPH.txt --source S

prints one distance per line (``p/q`` or an integer, ``inf`` if unreachable).
"""

from __future__ import annotations

import argparse
import heapq
import sys
from collections.abc import Sequence
from fractions import Fraction
from pathlib import Path

Edge = tuple[int, int, Fraction]
Distances = list[Fraction | None]


# ---------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------


def parse_weight(token: str) -> Fraction:
    """Read a non-negative edge length exactly."""
    try:
        w = Fraction(token)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"invalid edge weight {token!r}") from exc
    if w < 0:
        raise ValueError(f"negative edge weight {token!r}")
    return w


def _parse_int(token: str, what: str) -> int:
    try:
        return int(token)
    except ValueError as exc:
        raise ValueError(f"invalid {what} {token!r}") from exc


def parse_edge_list(text: str) -> tuple[int, list[Edge]]:
    """Parse the edge-list format; raise ``ValueError`` on any malformed input."""
    lines = [
        (lineno, line.split())
        for lineno, raw in enumerate(text.splitlines(), start=1)
        if (line := raw.strip()) and not line.startswith("#")
    ]
    if not lines:
        raise ValueError("missing header line 'n m'")
    lineno, header = lines[0]
    if len(header) != 2:
        raise ValueError(f"line {lineno}: header must be 'n m', got {' '.join(header)!r}")
    n = _parse_int(header[0], "vertex count")
    m = _parse_int(header[1], "edge count")
    if n < 0 or m < 0:
        raise ValueError(f"line {lineno}: n and m must be non-negative")
    body = lines[1:]
    if len(body) != m:
        raise ValueError(f"header declares m={m} edges but {len(body)} edge lines follow")
    edges: list[Edge] = []
    for lineno, tokens in body:
        if len(tokens) != 3:
            raise ValueError(f"line {lineno}: expected 'u v w', got {' '.join(tokens)!r}")
        u = _parse_int(tokens[0], "vertex")
        v = _parse_int(tokens[1], "vertex")
        for x in (u, v):
            if not 0 <= x < n:
                raise ValueError(f"line {lineno}: vertex {x} out of range [0, {n})")
        try:
            w = parse_weight(tokens[2])
        except ValueError as exc:
            raise ValueError(f"line {lineno}: {exc}") from exc
        edges.append((u, v, w))
    return n, edges


def read_graph(path: str | Path) -> tuple[int, list[Edge]]:
    return parse_edge_list(Path(path).read_text())


# ---------------------------------------------------------------------------
# validation shared by both algorithms (input checks only, no SSSP logic)
# ---------------------------------------------------------------------------


def _validate(n: int, edges: Sequence[tuple[int, int, Fraction]], source: int) -> None:
    if n < 0:
        raise ValueError(f"negative vertex count {n}")
    if not 0 <= source < n:
        raise ValueError(f"source {source} out of range [0, {n})")
    for i, (u, v, w) in enumerate(edges):
        if not (0 <= u < n and 0 <= v < n):
            raise ValueError(f"edge {i} ({u}, {v}) out of range [0, {n})")
        if w < 0:
            raise ValueError(f"edge {i} has negative weight {w}")


# ---------------------------------------------------------------------------
# Dijkstra (binary heap, lazy deletion)
# ---------------------------------------------------------------------------


def dijkstra_with_parent(
    n: int, edges: Sequence[tuple[int, int, Fraction]], source: int
) -> tuple[Distances, list[int]]:
    _validate(n, edges, source)
    adj: list[list[tuple[int, Fraction]]] = [[] for _ in range(n)]
    for u, v, w in edges:
        adj[u].append((v, Fraction(w)))
        adj[v].append((u, Fraction(w)))
    dist: Distances = [None] * n
    parent = [-1] * n
    settled = [False] * n
    dist[source] = Fraction(0)
    heap: list[tuple[Fraction, int]] = [(Fraction(0), source)]
    while heap:
        d, x = heapq.heappop(heap)
        if settled[x]:
            continue
        settled[x] = True
        for y, w in adj[x]:
            nd = d + w
            dy = dist[y]
            if not settled[y] and (dy is None or nd < dy):
                dist[y] = nd
                parent[y] = x
                heapq.heappush(heap, (nd, y))
    return dist, parent


def dijkstra(n: int, edges: Sequence[tuple[int, int, Fraction]], source: int) -> Distances:
    return dijkstra_with_parent(n, edges, source)[0]


# ---------------------------------------------------------------------------
# Bellman-Ford (round based: relax every edge each round, no queue)
# ---------------------------------------------------------------------------


def bellman_ford_with_parent(
    n: int, edges: Sequence[tuple[int, int, Fraction]], source: int
) -> tuple[Distances, list[int]]:
    _validate(n, edges, source)
    arcs = [(u, v, Fraction(w)) for u, v, w in edges] + [(v, u, Fraction(w)) for u, v, w in edges]
    dist: Distances = [None] * n
    parent = [-1] * n
    dist[source] = Fraction(0)
    for _ in range(max(n - 1, 0)):
        changed = False
        for a, b, w in arcs:
            da = dist[a]
            if da is None:
                continue
            db = dist[b]
            if db is None or da + w < db:
                dist[b] = da + w
                parent[b] = a
                changed = True
        if not changed:
            break
    return dist, parent


def bellman_ford(n: int, edges: Sequence[tuple[int, int, Fraction]], source: int) -> Distances:
    return bellman_ford_with_parent(n, edges, source)[0]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def format_distance(d: Fraction | None) -> str:
    return "inf" if d is None else str(d)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Exact rational SSSP (reference oracle).")
    parser.add_argument("graph", type=Path, help="edge-list text file")
    parser.add_argument("--source", type=int, default=0, help="source vertex (default 0)")
    parser.add_argument("--algorithm", choices=("dijkstra", "bellman-ford"), default="dijkstra")
    args = parser.parse_args(argv)
    try:
        n, edges = read_graph(args.graph)
        sssp = dijkstra if args.algorithm == "dijkstra" else bellman_ford
        dist = sssp(n, edges, args.source)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    sys.stdout.write("".join(format_distance(d) + "\n" for d in dist))
    return 0


if __name__ == "__main__":
    sys.exit(main())
