"""Tests for the exact rational SSSP reference and the double-result checker.

The reference (``reference.exact_sssp``) computes shortest-path distances with
``fractions.Fraction`` so that it never rounds.  The checker
(``reference.check_double``) verifies a C++ double result against it with the
rounding-error bound gamma_k = k*u / (1 - k*u), u = 2^-53.
"""

from __future__ import annotations

import heapq
import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import networkx as nx
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from reference import check_double, exact_sssp

REFERENCE_DIR = Path(__file__).resolve().parent / "reference"
EXACT_CLI = REFERENCE_DIR / "exact_sssp.py"
CHECK_CLI = REFERENCE_DIR / "check_double.py"

U = Fraction(1, 2**53)

F = Fraction
DETERMINISTIC = settings(
    derandomize=True,
    max_examples=150,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def graph_text(n: int, edges: list[tuple[int, int, str]]) -> str:
    lines = [f"{n} {len(edges)}"]
    lines += [f"{u} {v} {w}" for u, v, w in edges]
    return "\n".join(lines) + "\n"


def float_dijkstra(n: int, edges, source: int):
    """Plain IEEE-double Dijkstra, standing in for the C++ implementation."""
    adj: list[list[tuple[int, float, int]]] = [[] for _ in range(n)]
    for idx, (u, v, w) in enumerate(edges):
        wf = float(w)
        adj[u].append((v, wf, idx))
        adj[v].append((u, wf, idx))
    inf = float("inf")
    dist = [inf] * n
    parent = [-1] * n
    parent_edge = [-1] * n
    dist[source] = 0.0
    heap = [(0.0, source)]
    done = [False] * n
    while heap:
        d, x = heapq.heappop(heap)
        if done[x]:
            continue
        done[x] = True
        for y, wf, idx in adj[x]:
            nd = d + wf
            if nd < dist[y]:
                dist[y] = nd
                parent[y] = x
                parent_edge[y] = idx
                heapq.heappush(heap, (nd, y))
    return dist, parent, parent_edge


def dump(n: int, source: int, dist, parent=None, parent_edge=None) -> dict:
    out: dict = {
        "n": n,
        "source": source,
        "dist": [None if d == float("inf") else d for d in dist],
    }
    if parent is not None:
        out["parent"] = list(parent)
    if parent_edge is not None:
        out["parent_edge"] = list(parent_edge)
    return out


def dump_json(result: dict) -> str:
    """Serialise like the C++ side: every double printed with %.17g."""
    parts = []
    for key, value in result.items():
        if key == "dist":
            items = ["null" if d is None else f"{d:.17g}" for d in value]
            parts.append(f'"dist": [{", ".join(items)}]')
        else:
            parts.append(f"{json.dumps(key)}: {json.dumps(value)}")
    return "{" + ", ".join(parts) + "}"


def float_dump(n: int, edges, source: int, *, with_parent_edge: bool = True) -> dict:
    dist, parent, parent_edge = float_dijkstra(n, edges, source)
    return dump(n, source, dist, parent, parent_edge if with_parent_edge else None)


def nx_distances(n: int, edges, source: int):
    g = nx.Graph()
    g.add_nodes_from(range(n))
    for u, v, w in edges:
        if g.has_edge(u, v):
            g[u][v]["weight"] = min(g[u][v]["weight"], w)
        else:
            g.add_edge(u, v, weight=w)
    lengths = nx.single_source_dijkstra_path_length(g, source, weight="weight")
    return [Fraction(lengths[v]) if v in lengths else None for v in range(n)]


# ---------------------------------------------------------------------------
# parser
# ---------------------------------------------------------------------------


class TestParser:
    def test_integer_fraction_and_decimal_weights_are_exact(self):
        n, edges = exact_sssp.parse_edge_list("3 3\n0 1 3\n1 2 1/3\n0 2 0.1\n")
        assert n == 3
        assert edges == [(0, 1, F(3)), (1, 2, F(1, 3)), (0, 2, F(1, 10))]
        assert all(isinstance(w, Fraction) for _, _, w in edges)

    def test_exponent_weight_is_exact(self):
        _, edges = exact_sssp.parse_edge_list("2 1\n0 1 1e-3\n")
        assert edges == [(0, 1, F(1, 1000))]

    def test_blank_lines_and_comments_are_ignored(self):
        text = "# header comment\n\n2 1\n   \n# edge follows\n0 1 5\n\n"
        assert exact_sssp.parse_edge_list(text) == (2, [(0, 1, F(5))])

    def test_multi_edges_and_self_loops_are_kept(self):
        _, edges = exact_sssp.parse_edge_list("2 3\n0 1 2\n0 1 1\n1 1 7\n")
        assert edges == [(0, 1, F(2)), (0, 1, F(1)), (1, 1, F(7))]

    def test_zero_weight_and_no_edges(self):
        assert exact_sssp.parse_edge_list("1 0\n") == (1, [])
        assert exact_sssp.parse_edge_list("2 1\n0 1 0\n") == (2, [(0, 1, F(0))])

    @pytest.mark.parametrize(
        "text",
        [
            "",  # no header
            "# only a comment\n",
            "2\n",  # header with one token
            "2 1 9\n0 1 1\n",  # header with three tokens
            "-1 0\n",  # negative n
            "2 -1\n",  # negative m
            "2 2\n0 1 1\n",  # fewer edge lines than m
            "2 1\n0 1 1\n1 0 1\n",  # more edge lines than m
            "2 1\n0 2 1\n",  # vertex out of range
            "2 1\n-1 1 1\n",  # negative vertex
            "2 1\n0 1 -1\n",  # negative weight
            "2 1\n0 1 -1/2\n",  # negative fractional weight
            "2 1\n0 1\n",  # missing weight
            "2 1\n0 1 1 1\n",  # extra token
            "2 1\n0 1 abc\n",  # not a number
            "2 1\n0 1 inf\n",  # infinity is not a length
            "2 1\n0 1 nan\n",
            "2 1\n0 1 1/0\n",  # division by zero
            "2 1\n0.5 1 1\n",  # non-integer vertex
            "a b\n",
        ],
    )
    def test_malformed_input_raises_value_error(self, text):
        with pytest.raises(ValueError):
            exact_sssp.parse_edge_list(text)

    def test_read_graph_from_file(self, tmp_path):
        path = tmp_path / "g.txt"
        path.write_text("2 1\n0 1 1/2\n")
        assert exact_sssp.read_graph(path) == (2, [(0, 1, F(1, 2))])


# ---------------------------------------------------------------------------
# hand-written SSSP cases (both algorithms)
# ---------------------------------------------------------------------------

ALGORITHMS = [
    pytest.param(exact_sssp.dijkstra, id="dijkstra"),
    pytest.param(exact_sssp.bellman_ford, id="bellman_ford"),
]
TREE_ALGORITHMS = [
    pytest.param(exact_sssp.dijkstra_with_parent, id="dijkstra"),
    pytest.param(exact_sssp.bellman_ford_with_parent, id="bellman_ford"),
]


@pytest.mark.parametrize("sssp", ALGORITHMS)
class TestHandWritten:
    def test_single_vertex(self, sssp):
        assert sssp(1, [], 0) == [F(0)]

    def test_no_edges(self, sssp):
        assert sssp(3, [], 1) == [None, F(0), None]

    def test_path(self, sssp):
        edges = [(0, 1, F(2)), (1, 2, F(3)), (2, 3, F(4))]
        assert sssp(4, edges, 0) == [F(0), F(2), F(5), F(9)]
        assert sssp(4, edges, 3) == [F(9), F(7), F(4), F(0)]

    def test_thirds_sum_to_exactly_one(self, sssp):
        _, edges = exact_sssp.parse_edge_list("4 4\n0 1 1/3\n1 2 1/3\n2 3 1/3\n0 3 1\n")
        dist = sssp(4, edges, 0)
        assert dist[3] == F(1)
        assert dist == [F(0), F(1, 3), F(2, 3), F(1)]

    def test_decimal_tenths_are_exact(self, sssp):
        # In IEEE double 0.1+0.1+0.1 != 0.3; exactly it is 3/10.
        _, edges = exact_sssp.parse_edge_list("4 3\n0 1 0.1\n1 2 0.1\n2 3 0.1\n")
        assert sssp(4, edges, 0)[3] == F(3, 10)
        assert 0.1 + 0.1 + 0.1 != 0.3

    def test_zero_weights(self, sssp):
        edges = [(0, 1, F(0)), (1, 2, F(0)), (2, 0, F(0)), (2, 3, F(1, 2))]
        assert sssp(4, edges, 0) == [F(0), F(0), F(0), F(1, 2)]

    def test_multi_edges_use_the_lightest(self, sssp):
        edges = [(0, 1, F(5)), (0, 1, F(2)), (1, 0, F(3))]
        assert sssp(2, edges, 0) == [F(0), F(2)]

    def test_self_loops_are_harmless(self, sssp):
        edges = [(0, 0, F(1)), (0, 1, F(4)), (1, 1, F(0))]
        assert sssp(2, edges, 0) == [F(0), F(4)]

    def test_disconnected(self, sssp):
        edges = [(0, 1, F(1)), (2, 3, F(1))]
        assert sssp(4, edges, 0) == [F(0), F(1), None, None]

    def test_shortcut_through_more_hops(self, sssp):
        edges = [(0, 3, F(10)), (0, 1, F(1)), (1, 2, F(1)), (2, 3, F(1))]
        assert sssp(4, edges, 0)[3] == F(3)

    def test_invalid_source_raises(self, sssp):
        with pytest.raises(ValueError):
            sssp(2, [], 2)
        with pytest.raises(ValueError):
            sssp(2, [], -1)

    def test_negative_weight_raises(self, sssp):
        with pytest.raises(ValueError):
            sssp(2, [(0, 1, F(-1))], 0)

    def test_out_of_range_edge_raises(self, sssp):
        with pytest.raises(ValueError):
            sssp(2, [(0, 5, F(1))], 0)


@pytest.mark.parametrize("sssp_tree", TREE_ALGORITHMS)
def test_parent_versions_return_shortest_path_tree(sssp_tree):
    edges = [(0, 1, F(1)), (1, 2, F(1)), (0, 2, F(5)), (3, 4, F(1))]
    dist, parent = sssp_tree(5, edges, 0)
    assert dist == [F(0), F(1), F(2), None, None]
    assert parent == [-1, 0, 1, -1, -1]


# ---------------------------------------------------------------------------
# property tests
# ---------------------------------------------------------------------------

weights_st = st.builds(Fraction, st.integers(0, 20), st.integers(1, 6))


@st.composite
def small_graphs(draw, weights=weights_st):
    n = draw(st.integers(1, 8))
    edge = st.tuples(st.integers(0, n - 1), st.integers(0, n - 1), weights)
    edges = draw(st.lists(edge, max_size=16))
    source = draw(st.integers(0, n - 1))
    return n, edges, source


@DETERMINISTIC
@given(small_graphs())
def test_dijkstra_matches_bellman_ford_and_networkx(graph):
    n, edges, source = graph
    dist = exact_sssp.dijkstra(n, edges, source)
    assert dist == exact_sssp.bellman_ford(n, edges, source)
    assert dist == nx_distances(n, edges, source)
    assert dist[source] == 0


@DETERMINISTIC
@given(small_graphs())
def test_distances_satisfy_every_arc(graph):
    n, edges, source = graph
    dist = exact_sssp.dijkstra(n, edges, source)
    for u, v, w in edges:
        for a, b in ((u, v), (v, u)):
            if dist[a] is not None:
                assert dist[b] is not None
                assert dist[b] <= dist[a] + w


@DETERMINISTIC
@given(small_graphs())
def test_parent_trees_are_tight(graph):
    n, edges, source = graph
    for sssp_tree in (exact_sssp.dijkstra_with_parent, exact_sssp.bellman_ford_with_parent):
        dist, parent = sssp_tree(n, edges, source)
        assert parent[source] == -1
        for v in range(n):
            if v == source or dist[v] is None:
                assert parent[v] == -1
                continue
            p = parent[v]
            lightest = min(w for a, b, w in edges if {a, b} == {p, v} and a != b)
            assert dist[v] == dist[p] + lightest


# ---------------------------------------------------------------------------
# check_double: tolerance definition
# ---------------------------------------------------------------------------


class TestTolerance:
    def test_unit_roundoff(self):
        assert check_double.UNIT_ROUNDOFF == U

    def test_gamma_definition(self):
        assert check_double.gamma(0) == 0
        assert check_double.gamma(1) == U / (1 - U)
        assert check_double.gamma(10) == 10 * U / (1 - 10 * U)
        assert isinstance(check_double.gamma(3), Fraction)

    def test_gamma_requires_ku_below_one(self):
        with pytest.raises(ValueError):
            check_double.gamma(2**53)
        with pytest.raises(ValueError):
            check_double.gamma(-1)


# ---------------------------------------------------------------------------
# check_double: accepting correct double results
# ---------------------------------------------------------------------------

PATH_TENTHS = [(0, 1, F(1, 10)), (1, 2, F(1, 10)), (2, 3, F(1, 10)), (0, 3, F(1))]


class TestAccept:
    def test_float_dijkstra_dict_is_accepted(self):
        report = check_double.check(4, PATH_TENTHS, float_dump(4, PATH_TENTHS, 0))
        assert report.ok, report.failures
        assert report.failures == []

    def test_float_dijkstra_json_is_accepted(self):
        text = dump_json(float_dump(4, PATH_TENTHS, 0))
        result = check_double.parse_result_json(text)
        assert check_double.check(4, PATH_TENTHS, result).ok

    def test_parent_edge_is_optional(self):
        result = float_dump(4, PATH_TENTHS, 0, with_parent_edge=False)
        assert check_double.check(4, PATH_TENTHS, result).ok

    def test_parent_is_optional(self):
        result = float_dump(4, PATH_TENTHS, 0)
        del result["parent"]
        del result["parent_edge"]
        assert check_double.check(4, PATH_TENTHS, result).ok

    def test_double_result_differs_from_exact_but_is_accepted(self):
        result = float_dump(4, PATH_TENTHS, 0)
        exact = exact_sssp.dijkstra(4, PATH_TENTHS, 0)
        # Neither the double sum nor the double weights equal 3/10 exactly.
        assert Fraction(result["dist"][3]) != exact[3]
        assert check_double.check(4, PATH_TENTHS, result).ok

    def test_unreachable_vertices(self):
        edges = [(0, 1, F(1, 3)), (2, 3, F(1))]
        result = float_dump(4, edges, 0)
        assert result["dist"][2] is None
        report = check_double.check(4, edges, result)
        assert report.ok, report.failures

    @pytest.mark.parametrize("marker", ['"inf"', "Infinity", "null"])
    def test_unreachable_markers_in_json(self, marker):
        text = (
            '{"n": 3, "source": 0, "dist": [0, 0.5, '
            + marker
            + '], "parent": [-1, 0, -1], "parent_edge": [-1, 0, -1]}'
        )
        result = check_double.parse_result_json(text)
        assert result["dist"][2] is None
        assert check_double.check(3, [(0, 1, F(1, 2))], result).ok

    def test_json_numbers_are_read_as_exact_doubles(self):
        result = check_double.parse_result_json(
            '{"n": 2, "source": 0, "dist": [0, 0.10000000000000001], "parent": [-1, 0]}'
        )
        assert result["dist"] == [F(0), Fraction(0.1)]
        assert all(isinstance(d, Fraction) for d in result["dist"])
        assert result["parent"] == [-1, 0]

    def test_float_absorption_is_within_gamma(self):
        # 1e17 + 1 == 1e17 in double, so the C++ distance is off by one exactly.
        assert 1e17 + 1 == 1e17
        edges = [(0, 1, F(10**17)), (1, 2, F(1))]
        result = float_dump(3, edges, 0)
        assert result["dist"][2] == 1e17
        assert exact_sssp.dijkstra(3, edges, 0)[2] == 10**17 + 1
        report = check_double.check(3, edges, result)
        assert report.ok, report.failures
        # ...but a zero tolerance must catch the difference.
        assert not check_double.check(3, edges, result, rel=0).ok

    def test_multi_edge_parent_edge_with_lightest_edge(self):
        edges = [(0, 1, F(5)), (0, 1, F(1))]
        result = float_dump(2, edges, 0)
        assert result["parent_edge"] == [-1, 1]
        assert check_double.check(2, edges, result).ok


@DETERMINISTIC
@given(
    small_graphs(
        weights=st.one_of(
            st.integers(0, 1000).map(str),
            st.tuples(st.integers(0, 99), st.integers(0, 999)).map(lambda t: f"{t[0]}.{t[1]}"),
            st.tuples(st.integers(1, 9), st.integers(1, 20)).map(lambda t: f"{t[0]}e-{t[1]}"),
        )
    )
)
def test_float_dijkstra_is_always_accepted(graph):
    n, str_edges, source = graph
    n, edges = exact_sssp.parse_edge_list(graph_text(n, str_edges))
    # The C++ side reads each weight string with strtod.
    double_edges = [(u, v, float(w)) for u, v, w in str_edges]
    result = dump(n, source, *float_dijkstra(n, double_edges, source))
    report = check_double.check(n, edges, check_double.parse_result_json(dump_json(result)))
    assert report.ok, report.failures


# ---------------------------------------------------------------------------
# check_double: tolerance boundary
# ---------------------------------------------------------------------------


class TestBoundary:
    EDGE = [(0, 1, F(1))]

    def test_next_double_above_needs_two_hops(self):
        # d* = 1, computed 1 + 2^-52: error 2^-52 > gamma_1 but <= gamma_2.
        result = {"n": 2, "source": 0, "dist": [0.0, 1.0 + 2.0**-52], "parent": [-1, 0]}
        assert not check_double.check(2, self.EDGE, result).ok
        assert not check_double.check(2, self.EDGE, result, hops=1).ok
        assert check_double.check(2, self.EDGE, result, hops=2).ok

    def test_next_double_below_is_within_one_hop(self):
        # error 2^-53 <= gamma_1 = u / (1 - u)
        result = {"n": 2, "source": 0, "dist": [0.0, 1.0 - 2.0**-53], "parent": [-1, 0]}
        assert check_double.check(2, self.EDGE, result).ok
        assert not check_double.check(2, self.EDGE, result, hops=0).ok

    def test_exactly_gamma_is_accepted_and_just_outside_is_rejected(self):
        k = 1
        g = check_double.gamma(k)
        d_star = F(3)
        edges = [(0, 1, d_star)]
        inside = {"n": 2, "source": 0, "dist": [F(0), d_star * (1 + g)]}
        outside = {"n": 2, "source": 0, "dist": [F(0), d_star * (1 + g) + F(1, 10**30)]}
        below = {"n": 2, "source": 0, "dist": [F(0), d_star * (1 - g)]}
        below_out = {"n": 2, "source": 0, "dist": [F(0), d_star * (1 - g) - F(1, 10**30)]}
        assert check_double.check(2, edges, inside, hops=k).ok
        assert not check_double.check(2, edges, outside, hops=k).ok
        assert check_double.check(2, edges, below, hops=k).ok
        assert not check_double.check(2, edges, below_out, hops=k).ok

    def test_absolute_tolerance(self):
        edges = [(0, 1, F(3))]
        ok = {"n": 2, "source": 0, "dist": [F(0), F(3) + F(1, 1000)]}
        bad = {"n": 2, "source": 0, "dist": [F(0), F(3) + F(1, 1000) + F(1, 10**30)]}
        assert check_double.check(2, edges, ok, rel=0, abs_tol=F(1, 1000)).ok
        assert not check_double.check(2, edges, bad, rel=0, abs_tol=F(1, 1000)).ok

    def test_default_hops_is_n_minus_one(self):
        report = check_double.check(
            5, [], {"n": 5, "source": 0, "dist": [0, None, None, None, None]}
        )
        assert report.ok
        assert report.hops == 4
        assert report.rel == check_double.gamma(4)
        assert report.abs_tol == 0

    def test_weights_double_vs_exact(self):
        # 0.1 as a double is slightly larger than 1/10.
        edges = [(0, 1, F(1, 10))]
        exact_answer = {"n": 2, "source": 0, "dist": [F(0), F(1, 10)]}
        double_answer = {"n": 2, "source": 0, "dist": [F(0), Fraction(0.1)]}
        assert check_double.check(2, edges, double_answer, rel=0).ok
        assert not check_double.check(2, edges, exact_answer, rel=0).ok
        assert check_double.check(2, edges, exact_answer, rel=0, weights="exact").ok
        assert not check_double.check(2, edges, double_answer, rel=0, weights="exact").ok

    def test_unknown_weights_mode_raises(self):
        with pytest.raises(ValueError):
            check_double.check(1, [], {"n": 1, "source": 0, "dist": [0]}, weights="single")


# ---------------------------------------------------------------------------
# check_double: rejecting broken results
# ---------------------------------------------------------------------------

TRIANGLE = [(0, 1, F(1)), (1, 2, F(1)), (0, 2, F(10))]


def triangle_dump() -> dict:
    return float_dump(3, TRIANGLE, 0)


class TestReject:
    def test_baseline_triangle_is_accepted(self):
        assert check_double.check(3, TRIANGLE, triangle_dump()).ok

    @pytest.mark.parametrize("value", [3.0, 1.5, 2.001, 1.0, 0.0])
    def test_shifted_distance(self, value):
        result = triangle_dump()
        result["dist"][2] = value
        report = check_double.check(3, TRIANGLE, result)
        assert not report.ok
        assert any("2" in f for f in report.failures)

    def test_source_distance_not_zero(self):
        result = triangle_dump()
        result["dist"][0] = 1e-300
        assert not check_double.check(3, TRIANGLE, result).ok

    def test_source_distance_not_zero_even_within_abs_tolerance(self):
        # dist[source] must be exactly 0, not merely within the tolerance.
        result = triangle_dump()
        result["dist"][0] = 1e-6
        report = check_double.check(3, TRIANGLE, result, abs_tol=F(1, 1000))
        assert not report.ok
        assert any("source" in f for f in report.failures)

    def test_reachable_but_null(self):
        result = triangle_dump()
        result["dist"][1] = None
        assert not check_double.check(3, TRIANGLE, result).ok

    def test_unreachable_but_finite(self):
        edges = [(0, 1, F(1))]
        result = {"n": 3, "source": 0, "dist": [0.0, 1.0, 5.0], "parent": [-1, 0, -1]}
        assert not check_double.check(3, edges, result).ok

    def test_cyclic_parents(self):
        result = triangle_dump()
        result["parent"] = [-1, 2, 1]
        del result["parent_edge"]
        report = check_double.check(3, TRIANGLE, result)
        assert not report.ok
        assert any("cycle" in f for f in report.failures)

    def test_self_parent(self):
        edges = [(0, 1, F(1)), (1, 1, F(0))]
        result = {"n": 2, "source": 0, "dist": [0.0, 1.0], "parent": [-1, 1]}
        assert not check_double.check(2, edges, result).ok

    def test_parent_without_edge(self):
        edges = [(0, 1, F(1)), (1, 2, F(1))]
        result = {"n": 3, "source": 0, "dist": [0.0, 1.0, 2.0], "parent": [-1, 0, 0]}
        report = check_double.check(3, edges, result)
        assert not report.ok

    def test_parent_edge_with_wrong_endpoints(self):
        result = triangle_dump()
        result["parent_edge"][2] = 0  # edge 0 is (0, 1), not (1, 2)
        assert not check_double.check(3, TRIANGLE, result).ok

    @pytest.mark.parametrize("idx", [3, -2, 99])
    def test_parent_edge_out_of_range(self, idx):
        result = triangle_dump()
        result["parent_edge"][2] = idx
        assert not check_double.check(3, TRIANGLE, result).ok

    def test_parent_tree_not_shortest(self):
        # parent[2] = 0 uses the weight-10 edge although dist[2] = 2.
        result = triangle_dump()
        result["parent"][2] = 0
        result["parent_edge"][2] = 2
        report = check_double.check(3, TRIANGLE, result)
        assert not report.ok
        assert any("tree" in f for f in report.failures)

    def test_parent_tree_not_shortest_without_parent_edge(self):
        result = triangle_dump()
        result["parent"][2] = 0
        del result["parent_edge"]
        assert not check_double.check(3, TRIANGLE, result).ok

    def test_heavy_multi_edge_as_parent_edge(self):
        edges = [(0, 1, F(5)), (0, 1, F(1))]
        result = float_dump(2, edges, 0)
        result["parent_edge"][1] = 0  # the weight-5 copy
        assert not check_double.check(2, edges, result).ok
        del result["parent_edge"]  # without it the lightest copy is assumed
        assert check_double.check(2, edges, result).ok

    def test_reachable_vertex_without_parent(self):
        result = triangle_dump()
        result["parent"][2] = -1
        result["parent_edge"][2] = -1
        assert not check_double.check(3, TRIANGLE, result).ok

    def test_source_with_parent(self):
        result = triangle_dump()
        result["parent"][0] = 1
        result["parent_edge"][0] = 0
        assert not check_double.check(3, TRIANGLE, result).ok

    def test_unreachable_vertex_with_parent(self):
        edges = [(0, 1, F(1)), (2, 3, F(1))]
        result = float_dump(4, edges, 0)
        result["parent"][3] = 2
        result["parent_edge"][3] = 1
        assert not check_double.check(4, edges, result).ok

    def test_parent_out_of_range(self):
        result = triangle_dump()
        result["parent"][2] = 7
        assert not check_double.check(3, TRIANGLE, result).ok

    @pytest.mark.parametrize("key", ["dist", "parent", "parent_edge"])
    def test_wrong_lengths(self, key):
        result = triangle_dump()
        result[key] = result[key][:-1]
        assert not check_double.check(3, TRIANGLE, result).ok

    def test_wrong_n_or_source(self):
        result = triangle_dump()
        result["n"] = 4
        assert not check_double.check(3, TRIANGLE, result).ok
        result = triangle_dump()
        result["source"] = 1
        assert not check_double.check(3, TRIANGLE, result).ok

    def test_nan_distance(self):
        result = triangle_dump()
        result["dist"][1] = float("nan")
        assert not check_double.check(3, TRIANGLE, result).ok

    def test_json_nan_distance(self):
        text = dump_json(triangle_dump()).replace('"dist": [0', '"dist": [NaN', 1)
        with pytest.raises(ValueError):
            check_double.check(3, TRIANGLE, check_double.parse_result_json(text))


# ---------------------------------------------------------------------------
# CLIs
# ---------------------------------------------------------------------------


def run(*args) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, *map(str, args)], capture_output=True, text=True, timeout=60
    )


class TestCli:
    def test_exact_sssp_cli(self, tmp_path):
        graph = tmp_path / "g.txt"
        graph.write_text("# demo\n5 4\n0 1 1/3\n1 2 1/3\n2 3 1/3\n0 3 0.1e1\n")
        proc = run(EXACT_CLI, graph, "--source", 0)
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.splitlines() == ["0", "1/3", "2/3", "1", "inf"]

    def test_exact_sssp_cli_other_source(self, tmp_path):
        graph = tmp_path / "g.txt"
        graph.write_text("3 2\n0 1 0.5\n1 2 2\n")
        proc = run(EXACT_CLI, graph, "--source", 2)
        assert proc.returncode == 0, proc.stderr
        assert proc.stdout.splitlines() == ["5/2", "2", "0"]

    def test_exact_sssp_cli_bad_input(self, tmp_path):
        graph = tmp_path / "g.txt"
        graph.write_text("2 1\n0 1 -3\n")
        proc = run(EXACT_CLI, graph, "--source", 0)
        assert proc.returncode != 0
        assert "negative" in proc.stderr.lower()

    def write_case(self, tmp_path, result: dict):
        graph = tmp_path / "g.txt"
        graph.write_text(graph_text(3, [(0, 1, "1"), (1, 2, "0.1"), (0, 2, "10")]))
        res = tmp_path / "r.json"
        res.write_text(dump_json(result))
        return graph, res

    def cli_edges(self):
        return [(0, 1, 1.0), (1, 2, 0.1), (0, 2, 10.0)]

    def test_check_double_cli_accepts(self, tmp_path):
        result = dump(3, 0, *float_dijkstra(3, self.cli_edges(), 0))
        graph, res = self.write_case(tmp_path, result)
        proc = run(CHECK_CLI, graph, res)
        assert proc.returncode == 0, proc.stdout + proc.stderr

    def test_check_double_cli_rejects(self, tmp_path):
        result = dump(3, 0, *float_dijkstra(3, self.cli_edges(), 0))
        result["dist"][2] = 1.2
        graph, res = self.write_case(tmp_path, result)
        proc = run(CHECK_CLI, graph, res)
        assert proc.returncode == 1
        assert "2" in proc.stdout

    def test_check_double_cli_options(self, tmp_path):
        result = dump(3, 0, *float_dijkstra(3, self.cli_edges(), 0))
        result["dist"][2] = 1.1 + 2**-40
        graph, res = self.write_case(tmp_path, result)
        assert run(CHECK_CLI, graph, res).returncode == 1
        assert run(CHECK_CLI, graph, res, "--rel", "1e-9").returncode == 0
        assert run(CHECK_CLI, graph, res, "--rel", "0", "--abs", "1/1000").returncode == 0
        assert run(CHECK_CLI, graph, res, "--hops", "0").returncode == 1
        exact = tmp_path / "exact.json"
        exact.write_text('{"n": 3, "source": 0, "dist": [0, 1, 1.1], "parent": [-1, 0, 1]}')
        assert run(CHECK_CLI, graph, exact, "--rel", "0", "--weights", "double").returncode == 1
        # 1.1 printed by %.17g is a double, so it is never exactly 11/10.
        assert run(CHECK_CLI, graph, exact, "--rel", "0", "--weights", "exact").returncode == 1
        assert run(CHECK_CLI, graph, exact, "--rel", "1e-15", "--weights", "exact").returncode == 0


# ---------------------------------------------------------------------------
# contract with the C++ dumper (sssp::write_result_json)
# ---------------------------------------------------------------------------

# Byte-for-byte the string asserted in tests/unit/verify/test_verify.cpp
# ("write_result_json dumps round-trippable distances", double section).
CPP_DUMP = (
    '{"n":3,"source":0,"dist":[0,0.10000000000000001,null],"parent":[-1,0,-1],'
    '"parent_edge":[-1,0,-1]}\n'
)


def test_cpp_dump_format_is_readable_and_exact():
    result = check_double.parse_result_json(CPP_DUMP)
    assert result["dist"] == [F(0), F(0.1), None]  # exact binary value of the double 0.1
    edges = [(0, 1, F(1, 10))]
    assert check_double.check(3, edges, result).ok  # weight read as the double C++ used
    report = check_double.check(3, edges, result, weights="exact", rel=0)
    assert not report.ok  # double(0.1) != 1/10 exactly


def test_tree_bound_allows_rounding_on_both_sides():
    # Regression: the tree-path bound must be ((1+rel) d* + abs) / (1-rel), not (1+rel) d* + abs.
    # d* = 1 via edge 0; the result reports 1 + rel (allowed) through the heavier edge 1 whose
    # exact length L = 1.015 is itself within rel of that value, so it is a legitimate
    # rel-accurate result: (1 - rel) L <= dist <= (1 + rel) d*.
    result = {"n": 2, "source": 0, "dist": [F(0), F(101, 100)], "parent": [-1, 0]}
    result["parent_edge"] = [-1, 1]
    edges = [(0, 1, F(1)), (0, 1, F(1015, 1000))]
    report = check_double.check(2, edges, result, rel=F(1, 100), weights="exact")
    assert report.ok, report.failures
    heavier = [(0, 1, F(1)), (0, 1, F(1025, 1000))]  # L > 1.01 / 0.99 = 1.0202...
    assert not check_double.check(2, heavier, result, rel=F(1, 100), weights="exact").ok
