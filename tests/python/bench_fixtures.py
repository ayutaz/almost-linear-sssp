"""SYNTHETIC benchmark records for tests ONLY.

Every value produced here is fabricated test input. None of it is a measurement, and nothing
derived from it may appear in a report. Records are tagged with ``SYNTHETIC_MARKER`` in the
run id and CPU brand so that they can never be mistaken for real data.
"""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

SYNTHETIC_MARKER = "SYNTHETIC-TEST-DATA"

_BASE: dict[str, Any] = {
    "schema_version": "1",
    "run_id": f"{SYNTHETIC_MARKER}-run",
    "timestamp": "2026-11-27T10:00:00+09:00",
    "paper_version": "2611.34911v1",
    "git_sha": "0" * 40,
    "git_dirty": False,
    "algorithm": "dijkstra-binary",
    "mode": "baseline",
    "params": {},
    "os": "SyntheticOS 0.0",
    "cpu": {
        "brand": SYNTHETIC_MARKER,
        "logical_cores": 1,
        "physical_cores": 1,
        "perf_cores": None,
        "efficiency_cores": None,
    },
    "ram_bytes": 1 << 30,
    "arch": "synthetic",
    "compiler": "synthetic-cc 0",
    "build_type": "Release",
    "build_flags": "-O3",
    "diagnostics_enabled": False,
    "power": None,
    "load_average": None,
    "graph_type": "er",
    "n": 1000,
    "m": 2000,
    "avg_degree": 4.0,
    "weight_dist": "uniform_int",
    "weight_type": "uint64",
    "seed": 1,
    "source": 0,
    "reachable": 1000,
    "trial": 0,
    "warmup": False,
    "order_seed": 7,
    "isolated": True,
    "status": "ok",
    "error": None,
    "times_ns": {
        "graph_gen": 1000,
        "preprocess": None,
        "sssp": 5000,
        "tree": None,
        "end_to_end": 5000,
        "verify": 3000,
    },
    "peak_rss_bytes": 10_000_000,
    "rss_after_graph_bytes": 8_000_000,
    "ops": None,
    "counters": None,
    "fallback": None,
    "correctness": {"verified": True, "ok": True, "method": "reference:dijkstra-binary"},
}


def synthetic_record(**overrides: Any) -> dict[str, Any]:
    """A schema-valid SYNTHETIC record; top-level keys can be overridden."""
    rec = copy.deepcopy(_BASE)
    rec.update(copy.deepcopy(overrides))
    return rec


def synthetic_diag_record(**overrides: Any) -> dict[str, Any]:
    """A SYNTHETIC record from a diagnostics build (counters present)."""
    base = {
        "diagnostics_enabled": True,
        "ops": {
            "vertex_scans": 1000,
            "edge_scans": 4000,
            "relaxations": 1500,
            "comparisons": 9000,
            "additions": 4000,
            "queue_pushes": 1500,
            "queue_pops": 1500,
            "decrease_keys": 0,
        },
        "counters": {},
        "fallback": {"count": 0, "records": [], "vertices": 0, "edges": 0},
    }
    base.update(overrides)
    return synthetic_record(**base)


def write_jsonl(path: Path, records: list[dict[str, Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        for rec in records:
            fh.write(json.dumps(rec) + "\n")
    return path
