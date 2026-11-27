"""Cross-language check: C++ baselines and the Mode A entry (via ``sssp_run``) against the exact
Python oracle.

Graphs are produced by the C++ generator CLI ``sssp_gen``; every baseline and both Mode A variants
(``paper-faithful``, ``paper-faithful-d1``: the names of ``include/sssp/paper.hpp``) are run by
the C++ CLI ``sssp_run``; its JSON result is checked with ``reference.exact_sssp``
(``fractions.Fraction``, never rounds) and ``reference.check_double``, which share no code with
the C++ side.

* ``uint64`` results must equal the exact distances (no tolerance) and form a valid tight tree.
* ``double`` results must lie within gamma_{n-1} of the exact distances of the doubles actually
  used (docs/design/numerics.md Sec 3.1), and all relaxation-based baselines must print
  bit-identical distances (numerics.md Sec 2). The Mode A variants are held to gamma_{n-1}
  only: they are not relaxation baselines, so bit identity is not claimed for them.

The executables are looked up in ``$SSSP_BUILD_DIR`` or ``build/debug``. If ``sssp_gen`` is not
there the C++ project has not been built, and the module is skipped with that reason (CI runs
pytest after ctest, so it is found). If ``sssp_gen`` exists but ``sssp_run`` does not, the tests
fail: the build is incomplete.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import pytest
from reference import check_double, exact_sssp

REPO = Path(__file__).resolve().parents[2]
BUILD_DIR = Path(os.environ.get("SSSP_BUILD_DIR", REPO / "build" / "debug"))
SSSP_GEN = BUILD_DIR / "sssp_gen"
SSSP_RUN = BUILD_DIR / "sssp_run"
CHECK_CLI = Path(__file__).resolve().parent / "reference" / "check_double.py"

if not SSSP_GEN.is_file():
    pytest.skip(
        f"C++ tools not built: {SSSP_GEN.name} not found in "
        f"{BUILD_DIR.relative_to(REPO) if BUILD_DIR.is_relative_to(REPO) else BUILD_DIR}; "
        "run `cmake --preset debug && cmake --build --preset debug` or set SSSP_BUILD_DIR",
        allow_module_level=True,
    )

ALGORITHMS = [
    "bellman_ford",
    "dijkstra_binary",
    "dijkstra_4ary",
    "dijkstra_fibonacci",
    "dijkstra_radix",
]
PAPER_ALGORITHMS = ["paper-faithful", "paper-faithful-d1"]
MAX_U64 = 2**64 - 1


def algorithms_for(weight_type: str) -> list[str]:
    baselines = ALGORITHMS if weight_type == "uint64" else ALGORITHMS[:-1]
    return baselines + PAPER_ALGORITHMS


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def require_sssp_run() -> None:
    assert SSSP_RUN.is_file(), f"{SSSP_RUN.name} is missing next to {SSSP_GEN.name}"


def gen(tmp: Path, name: str, *args: str, fmt: str = "edgelist") -> Path:
    out = tmp / f"{name}.{'gr' if fmt == 'dimacs' else 'txt'}"
    subprocess.run(
        [str(SSSP_GEN), *args, "--format", fmt, "--out", str(out)],
        check=True,
        capture_output=True,
        text=True,
    )
    return out


def run(*args: str) -> subprocess.CompletedProcess[str]:
    require_sssp_run()
    return subprocess.run([str(SSSP_RUN), *args], capture_output=True, text=True, timeout=60)


def run_ok(graph: Path, algorithm: str, source: int, weight_type: str, *extra: str) -> dict:
    p = run(
        "--graph",
        str(graph),
        "--algorithm",
        algorithm,
        "--source",
        str(source),
        "--weight-type",
        weight_type,
        *extra,
    )
    assert p.returncode == 0, f"{algorithm}: exit {p.returncode}\n{p.stderr}"
    return check_double.parse_result_json(p.stdout)


def write_graph(tmp: Path, name: str, n: int, edges: list[tuple[int, int, int | str]]) -> Path:
    path = tmp / f"{name}.txt"
    lines = [f"{n} {len(edges)}"] + [f"{u} {v} {w}" for u, v, w in edges]
    path.write_text("\n".join(lines) + "\n")
    return path


def exact_tree_check(n, edges, result) -> list[str]:
    """uint64: exact distances and an exact tight tree (rel = abs = 0, exact weights)."""
    return check_double.check(n, edges, result, rel=0, abs_tol=0, weights="exact").failures


# ---------------------------------------------------------------------------
# generated graphs: every baseline vs the exact oracle
# ---------------------------------------------------------------------------

GRAPHS = [
    ("er", ["--n", "50", "--avg-degree", "4"]),
    ("er_sparse", ["--type", "er", "--n", "40", "--avg-degree", "1"]),
    ("grid", ["--rows", "6", "--cols", "7"]),
    ("road", ["--rows", "6", "--cols", "6"]),
    ("tree", ["--n", "45"]),
    ("ba", ["--n", "50", "--attach", "2"]),
    ("chung_lu", ["--n", "60", "--avg-degree", "3"]),
    ("regular", ["--n", "30", "--degree", "3", "--multi", "keep"]),
    ("adversarial_decrease_key", ["--n", "20"]),
    ("adversarial_long_path", ["--n", "40"]),
]
UINT64_WEIGHTS = [
    ["--weights", "uniform_int", "--lo", "0", "--hi", "3"],
    ["--weights", "wide_range", "--wide-log2", "50"],
    ["--weights", "zero_mixed", "--lo", "1", "--hi", "9"],
]
DOUBLE_WEIGHTS = [
    ["--weights", "uniform_real", "--lo", "0", "--hi", "1"],
    ["--weights", "log_uniform", "--lo", "1e-9", "--hi", "1e9"],
    ["--weights", "exponential", "--mean", "3"],
    ["--weights", "wide_range", "--wide-log2", "60"],
]
# One seed per (family, weight distribution), varied across distributions; each graph is solved
# from two sources by every algorithm. Kept small because every run is a process start.


def graph_args(kind: str, args: list[str]) -> list[str]:
    return args if "--type" in args else ["--type", kind, *args]


def cases(weight_type: str):
    weights = UINT64_WEIGHTS if weight_type == "uint64" else DOUBLE_WEIGHTS
    for kind, args in GRAPHS:
        for wi, w in enumerate(weights):
            seed = 1 + wi
            yield pytest.param(kind, args, w, seed, id=f"{kind}-{w[1]}{wi}-seed{seed}")


@pytest.mark.parametrize(("kind", "args", "weights", "seed"), list(cases("uint64")))
def test_uint64_results_are_exact(tmp_path, kind, args, weights, seed):
    path = gen(
        tmp_path,
        kind,
        *graph_args(kind, args),
        *weights,
        "--seed",
        str(seed),
        "--weight-type",
        "uint64",
    )
    n, edges = exact_sssp.read_graph(path)
    for source in sorted({0, (7 * seed + 3) % n}):
        exact = exact_sssp.dijkstra(n, edges, source)
        for alg in algorithms_for("uint64"):
            result = run_ok(path, alg, source, "uint64")
            # Integers are printed as JSON integers and read exactly.
            assert result["dist"] == exact, f"{alg} source {source}"
            assert exact_tree_check(n, edges, result) == [], f"{alg} source {source}"


@pytest.mark.parametrize(("kind", "args", "weights", "seed"), list(cases("double")))
def test_double_results_within_gamma_and_bit_identical(tmp_path, kind, args, weights, seed):
    path = gen(
        tmp_path,
        kind,
        *graph_args(kind, args),
        *weights,
        "--seed",
        str(seed),
        "--weight-type",
        "double",
    )
    n, edges = exact_sssp.read_graph(path)
    for source in sorted({0, (7 * seed + 3) % n}):
        printed: dict[str, list] = {}
        for alg in algorithms_for("double"):
            p = run(
                "--graph",
                str(path),
                "--algorithm",
                alg,
                "--source",
                str(source),
                "--weight-type",
                "double",
            )
            assert p.returncode == 0, p.stderr
            result = check_double.parse_result_json(p.stdout)
            report = check_double.check(n, edges, result)  # gamma_{n-1}, weights as doubles
            assert report.ok, f"{alg} source {source}: {report.failures[:5]}"
            printed[alg] = json.loads(p.stdout)["dist"]
        assert set(PAPER_ALGORITHMS) <= set(printed)
        first = printed["bellman_ford"]
        for alg, dist in printed.items():
            if alg in PAPER_ALGORITHMS:
                continue  # within gamma_{n-1} (checked above); bit identity is not claimed
            assert dist == first, f"{alg} differs from bellman_ford (must be bit-identical)"


def test_check_double_cli_accepts_sssp_run_output(tmp_path):
    """End to end through files: sssp_gen -> sssp_run --out -> check_double.py CLI."""
    path = gen(
        tmp_path,
        "er",
        "--type",
        "er",
        "--n",
        "40",
        "--weights",
        "log_uniform",
        "--lo",
        "1e-6",
        "--hi",
        "1e6",
        "--weight-type",
        "double",
        "--seed",
        "5",
    )
    out = tmp_path / "result.json"
    p = run(
        "--graph",
        str(path),
        "--algorithm",
        "dijkstra_fibonacci",
        "--source",
        "3",
        "--weight-type",
        "double",
        "--out",
        str(out),
    )
    assert p.returncode == 0, p.stderr
    assert p.stdout == ""
    c = subprocess.run(
        [sys.executable, str(CHECK_CLI), str(path), str(out)], capture_output=True, text=True
    )
    assert c.returncode == 0, c.stdout + c.stderr


# ---------------------------------------------------------------------------
# DIMACS input
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("weight_type", ["uint64", "double"])
def test_dimacs_input_matches_edge_list(tmp_path, weight_type):
    """The same generated graph as .gr and .txt gives the same result (arc pairs merged)."""
    common = [
        "--type",
        "road",
        "--rows",
        "5",
        "--cols",
        "6",
        "--seed",
        "3",
        "--weight-type",
        weight_type,
    ]
    txt = gen(tmp_path, "g", *common)
    gr = gen(tmp_path, "g", *common, fmt="dimacs")
    n, edges = exact_sssp.read_graph(txt)
    for alg in algorithms_for(weight_type):
        from_txt = run_ok(txt, alg, 4, weight_type)
        from_gr = run_ok(gr, alg, 4, weight_type)  # format chosen from the .gr extension
        assert from_gr == from_txt, alg
        if weight_type == "uint64":
            assert exact_tree_check(n, edges, from_gr) == []
        else:
            assert check_double.check(n, edges, from_gr).ok


# ---------------------------------------------------------------------------
# boundary inputs and CLI contract
# ---------------------------------------------------------------------------


def test_unreachable_vertices_are_null_with_no_parent(tmp_path):
    path = write_graph(tmp_path, "disc", 5, [(0, 1, 2), (3, 4, 1)])
    for alg in ALGORITHMS + PAPER_ALGORITHMS:
        r = run_ok(path, alg, 0, "uint64")
        assert r["dist"] == [0, 2, None, None, None]
        assert r["parent"] == [-1, 0, -1, -1, -1]
        assert r["parent_edge"] == [-1, 0, -1, -1, -1]
        assert r["n"] == 5 and r["source"] == 0


def test_largest_uint64_weights_are_parsed_exactly(tmp_path):
    """2^64 - 2 is the largest finite distance; Bellman-Ford returns it exactly."""
    path = write_graph(tmp_path, "big", 3, [(0, 1, MAX_U64 - 3), (1, 2, 1)])
    r = run_ok(path, "bellman_ford", 0, "uint64")
    assert r["dist"] == [0, MAX_U64 - 3, MAX_U64 - 2]


def test_overflow_behaviour_is_reported(tmp_path):
    """Dijkstra throws on any overflowing scanned arc; Bellman-Ford only if a true distance is
    not representable (docs/design/baselines.md). Both are reported as runtime errors."""
    heavy_back = write_graph(tmp_path, "hb", 3, [(0, 1, 1), (1, 2, 1), (0, 2, MAX_U64 - 1)])
    assert run_ok(heavy_back, "bellman_ford", 0, "uint64")["dist"] == [0, 1, 2]
    for alg in ALGORITHMS[1:]:
        p = run("--graph", str(heavy_back), "--algorithm", alg, "--source", "0")
        assert p.returncode == 1, (alg, p.returncode, p.stderr)
        assert "overflow" in p.stderr
        assert p.stdout == ""
    unrep = write_graph(tmp_path, "un", 3, [(0, 1, 2**63), (1, 2, 2**63)])
    for alg in ALGORITHMS + PAPER_ALGORITHMS:
        p = run("--graph", str(unrep), "--algorithm", alg, "--source", "0")
        assert p.returncode == 1, (alg, p.returncode, p.stderr)
        assert p.stdout == ""


def test_verify_flag_reports_ok(tmp_path):
    path = gen(tmp_path, "er", "--type", "er", "--n", "30", "--seed", "9")
    for alg in ALGORITHMS + PAPER_ALGORITHMS:
        p = run("--graph", str(path), "--algorithm", alg, "--source", "2", "--verify")
        assert p.returncode == 0, p.stderr
        assert "verify: OK" in p.stderr
        assert json.loads(p.stdout)["source"] == 2


def test_default_weight_type_is_uint64_and_double_is_read_exactly(tmp_path):
    path = write_graph(tmp_path, "w", 3, [(0, 1, "0.1"), (1, 2, "0.2")])
    p = run("--graph", str(path), "--algorithm", "dijkstra_binary", "--source", "0")
    assert p.returncode == 1  # "0.1" is not a uint64
    r = run_ok(path, "dijkstra_binary", 0, "double")
    assert r["dist"] == [0, Fraction(0.1), Fraction(0.1 + 0.2)]


@pytest.mark.parametrize(
    ("args", "code"),
    [
        (["--algorithm", "dijkstra_binary", "--source", "5"], 1),  # source out of range
        (["--algorithm", "paper-faithful", "--source", "5"], 1),
        (["--algorithm", "paper-faithful-d1", "--source", "5"], 1),
        (["--algorithm", "paper", "--source", "0"], 2),  # reserved schema name, not a variant
        (["--algorithm", "nope", "--source", "0"], 2),  # unknown algorithm
        (["--algorithm", "dijkstra_radix", "--source", "0", "--weight-type", "double"], 2),
        (["--source", "0"], 2),  # missing --algorithm
        (["--algorithm", "bellman_ford"], 2),  # missing --source
        (["--algorithm", "bellman_ford", "--source", "-1"], 2),
        (["--algorithm", "bellman_ford", "--source", "0", "--bogus"], 2),
    ],
)
def test_cli_errors(tmp_path, args, code):
    path = write_graph(tmp_path, "tri", 3, [(0, 1, 1), (1, 2, 1), (0, 2, 3)])
    p = run("--graph", str(path), *args)
    assert p.returncode == code, p.stderr
    assert p.stdout == ""
    assert p.stderr.startswith("sssp_run: ")


def test_cli_rejects_missing_and_malformed_files(tmp_path):
    p = run("--graph", str(tmp_path / "absent.txt"), "--algorithm", "bellman_ford", "--source", "0")
    assert p.returncode == 1 and p.stdout == ""
    bad = tmp_path / "bad.txt"
    bad.write_text("2 1\n0 5 1\n")  # vertex out of range
    p = run("--graph", str(bad), "--algorithm", "bellman_ford", "--source", "0")
    assert p.returncode == 1 and "line 2" in p.stderr
    p = run("--algorithm", "bellman_ford", "--source", "0")
    assert p.returncode == 2  # missing --graph
