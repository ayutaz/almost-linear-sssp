"""Benchmark job plans for scripts/run_benchmarks.sh (docs/design/benchmark_protocol.md §10).

A job is one runner process: one build x one graph configuration x ONE algorithm, so that
``ru_maxrss`` (Peak RSS) belongs to a single algorithm. The job order is shuffled with a fixed
``order_seed`` so that slow drifts (thermal, background load) do not favour one algorithm.

    PYTHONPATH=scripts python -m bench.plan smoke --order-seed 1   # TSV, one job per line
"""

from __future__ import annotations

import argparse
import random
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from typing import NamedTuple

from bench.schema import ALGORITHM_MODES, INTEGER_ONLY_ALGORITHMS

# Names are those of `sssp_gen --type` / `--weights` (docs/design/graph_library.md §5-6).
GRAPH_TYPES = (
    "er",
    "regular",
    "ba",
    "chung_lu",
    "grid",
    "road",
    "tree",
    "adversarial_long_path",
)
#: How the nominal average degree d is passed to each generator; None = fixed by the family.
DEGREE_ARG = {
    "er": "--avg-degree {d}",
    "chung_lu": "--avg-degree {d}",
    "regular": "--degree {d}",
    "ba": "--attach {half}",  # m attachments per new vertex -> average degree ~2m
    "grid": None,
    "road": None,
    "tree": None,
    "adversarial_long_path": None,
}
#: Generators whose weights are fixed by the construction (they ignore --weights).
BUILTIN_WEIGHTS = ("road", "adversarial_long_path")
WEIGHT_TYPE_OF = {
    "constant": "uint64",
    "uniform_int": "uint64",
    "log_uniform": "uint64",
    "wide_range": "uint64",
    "zero_mixed": "uint64",
    "uniform_real": "double",
    "exponential": "double",
    "builtin": "uint64",
}
WEIGHT_DISTS = tuple(WEIGHT_TYPE_OF)
NONE = "-"  # TSV placeholder for "not applicable"
BENCH_ALGORITHMS = (
    "dijkstra-binary",
    "dijkstra-4ary",
    "dijkstra-fibonacci",
    "dijkstra-radix",
    "paper",
)
BUILDS = ("release", "release-diag")
GRAPH_SEED = 1


@dataclass(frozen=True)
class Profile:
    sizes: tuple[int, ...]
    degrees: tuple[int, ...]
    graph_types: tuple[str, ...]
    weight_dists: tuple[str, ...]
    trials: int
    warmup: int
    diag_trials: int
    sources: int
    timeout_sec: int


PROFILES = {
    "smoke": Profile(
        sizes=(1_000, 10_000),
        degrees=(4,),
        graph_types=("er", "grid"),
        weight_dists=("uniform_int",),
        trials=1,
        warmup=0,
        diag_trials=1,
        sources=1,
        timeout_sec=60,
    ),
    "standard": Profile(
        sizes=(1_000, 10_000, 100_000, 1_000_000),
        degrees=(4, 16),
        graph_types=GRAPH_TYPES,
        weight_dists=("constant", "uniform_int", "wide_range"),
        trials=5,
        warmup=1,
        diag_trials=3,
        sources=2,
        timeout_sec=600,
    ),
    "full": Profile(
        sizes=(1_000, 10_000, 100_000, 1_000_000, 10_000_000),
        degrees=(2, 4, 8, 16, 32),
        graph_types=GRAPH_TYPES,
        weight_dists=WEIGHT_DISTS,
        trials=20,
        warmup=2,
        diag_trials=5,
        sources=4,
        timeout_sec=3600,
    ),
}


class Job(NamedTuple):
    build: str
    graph_type: str
    n: int
    avg_degree: str  # nominal d, or "-" for fixed-degree families
    shape_args: str  # generator shape option(s), or "-"
    weight_dist: str
    weight_type: str
    algorithm: str
    trials: int
    warmup: int
    sources: int
    seed: int


def _check(kind: str, values: Sequence[str] | None, allowed: Sequence[str]) -> None:
    for v in values or ():
        if v not in allowed:
            raise ValueError(f"unknown {kind} {v!r}; allowed: {list(allowed)}")


def _shapes(graph_type: str, degrees: Sequence[int]) -> list[tuple[str, str]]:
    fmt = DEGREE_ARG[graph_type]
    if fmt is None:
        return [(NONE, NONE)]
    return [(str(d), fmt.format(d=d, half=max(1, d // 2))) for d in degrees]


def build_jobs(
    profile: str,
    *,
    order_seed: int,
    algorithms: Sequence[str] | None = None,
    graph_types: Sequence[str] | None = None,
    sizes: Sequence[int] | None = None,
    builds: Sequence[str] | None = None,
) -> list[Job]:
    if profile not in PROFILES:
        raise ValueError(f"unknown profile {profile!r}; allowed: {list(PROFILES)}")
    _check("algorithm", algorithms, tuple(ALGORITHM_MODES))
    _check("graph type", graph_types, GRAPH_TYPES)
    _check("build", builds, BUILDS)
    p = PROFILES[profile]
    algos = tuple(algorithms) if algorithms else BENCH_ALGORITHMS
    gtypes = tuple(graph_types) if graph_types else p.graph_types
    ns = tuple(sizes) if sizes else p.sizes
    jobs: set[Job] = set()
    for build in builds or BUILDS:
        trials, warmup = (p.trials, p.warmup) if build == "release" else (p.diag_trials, 0)
        for gt in gtypes:
            dists = ("builtin",) if gt in BUILTIN_WEIGHTS else p.weight_dists
            for wd in dists:
                if wd == "builtin" and gt not in BUILTIN_WEIGHTS:
                    continue
                wt = WEIGHT_TYPE_OF[wd]
                for n in ns:
                    for deg, shape in _shapes(gt, p.degrees):
                        for a in algos:
                            if a in INTEGER_ONLY_ALGORITHMS and wt != "uint64":
                                continue
                            jobs.add(
                                Job(
                                    build=build,
                                    graph_type=gt,
                                    n=n,
                                    avg_degree=deg,
                                    shape_args=shape,
                                    weight_dist=wd,
                                    weight_type=wt,
                                    algorithm=a,
                                    trials=trials,
                                    warmup=warmup,
                                    sources=p.sources,
                                    seed=GRAPH_SEED,
                                )
                            )
    ordered = sorted(jobs)
    random.Random(order_seed).shuffle(ordered)
    return ordered


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="bench.plan", description=__doc__.splitlines()[0])
    ap.add_argument("profile", choices=list(PROFILES))
    ap.add_argument("--order-seed", type=int, required=True)
    ap.add_argument("--algorithms", help="comma-separated subset")
    ap.add_argument("--graph-types", help="comma-separated subset")
    ap.add_argument("--sizes", help="comma-separated vertex counts")
    ap.add_argument("--builds", help="comma-separated subset of release,release-diag")
    ap.add_argument("--timeout-sec", action="store_true", help="print the profile timeout only")
    args = ap.parse_args(argv)
    if args.timeout_sec:
        print(PROFILES[args.profile].timeout_sec)
        return 0

    def split(s: str | None) -> list[str] | None:
        return [x for x in s.split(",") if x] if s else None

    try:
        jobs = build_jobs(
            args.profile,
            order_seed=args.order_seed,
            algorithms=split(args.algorithms),
            graph_types=split(args.graph_types),
            sizes=[int(x) for x in split(args.sizes) or []] or None,
            builds=split(args.builds),
        )
    except ValueError as exc:
        print(f"bench.plan: {exc}", file=sys.stderr)
        return 2
    for j in jobs:
        print("\t".join(str(v) for v in j))
    return 0


if __name__ == "__main__":
    sys.exit(main())
