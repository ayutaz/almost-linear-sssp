"""Benchmark JSONL record schema (docs/design/benchmark_schema.md) and its validator.

One JSON object per line. Every key below is REQUIRED (present even when its value is null).
A value that was not measured is ``null`` -- never ``0``.

CLI (run from ``scripts/``, e.g. ``PYTHONPATH=scripts python -m bench.schema ...``)::

    python -m bench.schema validate FILE...        # exit 0 valid, 1 invalid
    python -m bench.schema failure --env-json E --out OUT --status timeout ...  # append record
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import re
import sys
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "1"
PAPER_VERSION = "2611.34911v1"

MODES = ("baseline", "paper-faithful", "paper-experimental", "practical-optimized")
STATUSES = ("ok", "timeout", "memory_limit", "error")
WEIGHT_TYPES = ("uint64", "double")

#: Known algorithm names and their mode. Unknown names are allowed with a non-baseline mode.
ALGORITHM_MODES: dict[str, str] = {
    "dijkstra-binary": "baseline",
    "dijkstra-4ary": "baseline",
    "dijkstra-fibonacci": "baseline",
    "dijkstra-radix": "baseline",
    "bellman-ford": "baseline",
    "paper": "paper-faithful",
    "paper-experimental": "paper-experimental",
    "paper-optimized": "practical-optimized",
}
BASELINE_ALGORITHMS = tuple(a for a, m in ALGORITHM_MODES.items() if m == "baseline")
INTEGER_ONLY_ALGORITHMS = ("dijkstra-radix",)

TIME_KEYS = ("graph_gen", "preprocess", "sssp", "tree", "end_to_end", "verify")
OPS_KEYS = (
    "vertex_scans",
    "edge_scans",
    "relaxations",
    "comparisons",
    "additions",
    "queue_pushes",
    "queue_pops",
    "decrease_keys",
)
CPU_KEYS = ("brand", "logical_cores", "physical_cores", "perf_cores", "efficiency_cores")
POWER_KEYS = ("source", "battery_percent", "low_power_mode", "cpu_speed_limit")

#: Fields filled from scripts/bench/envinfo.py (the runner copies them from --env-json).
ENV_FIELDS = (
    "paper_version",
    "git_sha",
    "git_dirty",
    "os",
    "cpu",
    "ram_bytes",
    "arch",
    "compiler",
    "build_type",
    "build_flags",
    "diagnostics_enabled",
    "power",
    "load_average",
)

_SHA = re.compile(r"[0-9a-f]{7,40}")
_PAPER = re.compile(r"\d{4}\.\d{4,5}v\d+")

Check = Callable[[Any], str | None]  # returns an error message or None


# ---------------------------------------------------------------- primitive checks
def _is_int(v: Any) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


def _int(lo: int | None = None) -> Check:
    def check(v: Any) -> str | None:
        if not _is_int(v):
            return f"must be an integer, got {v!r}"
        if lo is not None and v < lo:
            return f"must be >= {lo}, got {v}"
        return None

    return check


def _number(lo: float | None = None) -> Check:
    def check(v: Any) -> str | None:
        if isinstance(v, bool) or not isinstance(v, int | float) or not math.isfinite(v):
            return f"must be a finite number, got {v!r}"
        if lo is not None and v < lo:
            return f"must be >= {lo}, got {v}"
        return None

    return check


def _str(pattern: re.Pattern[str] | None = None) -> Check:
    def check(v: Any) -> str | None:
        if not isinstance(v, str) or not v:
            return f"must be a non-empty string, got {v!r}"
        if pattern is not None and not pattern.fullmatch(v):
            return f"must match {pattern.pattern}, got {v!r}"
        return None

    return check


def _enum(values: Iterable[str]) -> Check:
    allowed = tuple(values)

    def check(v: Any) -> str | None:
        return None if v in allowed else f"must be one of {list(allowed)}, got {v!r}"

    return check


def _bool(v: Any) -> str | None:
    return None if isinstance(v, bool) else f"must be a boolean, got {v!r}"


def _dict(v: Any) -> str | None:
    return None if isinstance(v, dict) else f"must be an object, got {type(v).__name__}"


def _const(value: str) -> Check:
    def check(v: Any) -> str | None:
        return None if v == value else f"must be {value!r}, got {v!r}"

    return check


def _timestamp(v: Any) -> str | None:
    if not isinstance(v, str):
        return f"must be an ISO-8601 string, got {v!r}"
    try:
        ts = _dt.datetime.fromisoformat(v)
    except ValueError:
        return f"must be ISO-8601, got {v!r}"
    return None if ts.tzinfo is not None else f"must include a UTC offset, got {v!r}"


def _load_average(v: Any) -> str | None:
    if not isinstance(v, list) or len(v) != 3:
        return f"must be a list of 3 numbers, got {v!r}"
    errs = [e for x in v if (e := _number(0)(x))]
    return errs[0] if errs else None


# (check, nullable) per top-level field. Order = documentation order.
FIELDS: dict[str, tuple[Check, bool]] = {
    "schema_version": (_const(SCHEMA_VERSION), False),
    "run_id": (_str(), False),
    "timestamp": (_timestamp, False),
    "paper_version": (_str(_PAPER), False),
    "git_sha": (_str(_SHA), False),
    "git_dirty": (_bool, False),
    "algorithm": (_str(), False),
    "mode": (_enum(MODES), False),
    "params": (_dict, False),
    "os": (_str(), False),
    "cpu": (_dict, False),
    "ram_bytes": (_int(1), True),
    "arch": (_str(), False),
    "compiler": (_str(), True),
    "build_type": (_str(), False),
    "build_flags": (_str(), True),
    "diagnostics_enabled": (_bool, False),
    "power": (_dict, True),
    "load_average": (_load_average, True),
    "graph_type": (_str(), False),
    "n": (_int(1), False),
    "m": (_int(0), True),  # null only when status != ok
    "avg_degree": (_number(0), True),
    "weight_dist": (_str(), False),
    "weight_type": (_enum(WEIGHT_TYPES), False),
    "seed": (_int(0), False),
    "source": (_int(0), True),  # null only when status != ok
    "reachable": (_int(0), True),
    "trial": (_int(0), True),  # null only when status != ok
    "warmup": (_bool, False),
    "order_seed": (_int(0), True),
    "isolated": (_bool, False),
    "status": (_enum(STATUSES), False),
    "error": (_str(), True),
    "times_ns": (_dict, False),
    "peak_rss_bytes": (_int(1), True),
    "rss_after_graph_bytes": (_int(1), True),
    "ops": (_dict, True),
    "counters": (_dict, True),
    "fallback": (_dict, True),
    "correctness": (_dict, False),
}
REQUIRED_FIELDS = tuple(FIELDS)


def _check_fields(rec: dict[str, Any], names: Iterable[str]) -> list[str]:
    errs: list[str] = []
    for name in names:
        if name not in rec:
            errs.append(f"{name}: missing required field")
            continue
        check, nullable = FIELDS[name]
        value = rec[name]
        if value is None:
            if not nullable:
                errs.append(f"{name}: must not be null")
            continue
        if (msg := check(value)) is not None:
            errs.append(f"{name}: {msg}")
    return errs


def _check_nested(rec: dict[str, Any]) -> list[str]:
    errs: list[str] = []
    cpu = rec.get("cpu")
    if isinstance(cpu, dict):
        errs += _exact_keys("cpu", cpu, CPU_KEYS)
        if cpu.get("brand") is not None and (m := _str()(cpu["brand"])):
            errs.append(f"cpu.brand: {m}")
        for k in CPU_KEYS[1:]:
            if cpu.get(k) is not None and (m := _int(0)(cpu[k])):
                errs.append(f"cpu.{k}: {m}")
    power = rec.get("power")
    if isinstance(power, dict):
        errs += _exact_keys("power", power, POWER_KEYS)
        if power.get("source") not in (None, "AC", "Battery", "UPS"):
            errs.append(f"power.source: must be AC/Battery/UPS or null, got {power['source']!r}")
        for k in ("battery_percent", "cpu_speed_limit"):
            if power.get(k) is not None and (m := _int(0)(power[k])):
                errs.append(f"power.{k}: {m}")
        if power.get("low_power_mode") is not None and (m := _bool(power["low_power_mode"])):
            errs.append(f"power.low_power_mode: {m}")

    times = rec.get("times_ns")
    if isinstance(times, dict):
        errs += _exact_keys("times_ns", times, TIME_KEYS)
        for k in TIME_KEYS:
            v = times.get(k)
            if v is None:
                continue
            if not _is_int(v):
                errs.append(f"times_ns.{k}: must be an integer number of ns or null, got {v!r}")
            elif v < 1:
                errs.append(f"times_ns.{k}: unmeasured values must be null, not {v}")
        e2e = times.get("end_to_end")
        if _is_int(e2e):
            for part in ("preprocess", "sssp", "tree"):
                p = times.get(part)
                if _is_int(p) and p > e2e:
                    errs.append(f"times_ns.end_to_end ({e2e}) < times_ns.{part} ({p})")

    for name in ("ops", "counters"):
        block = rec.get(name)
        if isinstance(block, dict):
            if name == "ops":
                errs += [f"ops.{k}: unknown counter" for k in block if k not in OPS_KEYS]
            for k, v in block.items():
                if (m := _int(0)(v)) is not None:
                    errs.append(f"{name}.{k}: {m}")

    fb = rec.get("fallback")
    if isinstance(fb, dict):
        errs += _check_fallback(fb)

    corr = rec.get("correctness")
    if isinstance(corr, dict):
        errs += _check_correctness(corr)
    return errs


def _exact_keys(prefix: str, block: dict[str, Any], keys: tuple[str, ...]) -> list[str]:
    errs = [f"{prefix}.{k}: missing required key" for k in keys if k not in block]
    errs += [f"{prefix}.{k}: unknown key" for k in block if k not in keys]
    return errs


def _check_fallback(fb: dict[str, Any]) -> list[str]:
    errs = _exact_keys("fallback", fb, ("count", "records", "vertices", "edges"))
    for k in ("count", "vertices", "edges"):
        if k in fb and (m := _int(0)(fb[k])):
            errs.append(f"fallback.{k}: {m}")
    records = fb.get("records")
    if not isinstance(records, list):
        return [*errs, "fallback.records: must be a list"]
    for i, r in enumerate(records):
        if not isinstance(r, dict):
            errs.append(f"fallback.records[{i}]: must be an object")
            continue
        errs += _exact_keys(f"fallback.records[{i}]", r, ("reason", "vertices", "edges"))
        if (m := _str()(r.get("reason"))) is not None:
            errs.append(f"fallback.records[{i}].reason: {m}")
        for k in ("vertices", "edges"):
            if (m := _int(0)(r.get(k))) is not None:
                errs.append(f"fallback.records[{i}].{k}: {m}")
    if errs:
        return errs
    if fb["count"] != len(records):
        errs.append(f"fallback.count ({fb['count']}) != len(fallback.records) ({len(records)})")
    for k in ("vertices", "edges"):
        total = sum(r[k] for r in records)
        if fb[k] != total:
            errs.append(f"fallback.{k} ({fb[k]}) != sum of records ({total})")
    return errs


def _check_correctness(c: dict[str, Any]) -> list[str]:
    errs = _exact_keys("correctness", c, ("verified", "ok", "method"))
    if errs:
        return errs
    if not isinstance(c["verified"], bool):
        return [f"correctness.verified: must be a boolean, got {c['verified']!r}"]
    if c["verified"]:
        if not isinstance(c["ok"], bool):
            errs.append("correctness.ok: must be a boolean when verified")
        if (m := _str()(c["method"])) is not None:
            errs.append(f"correctness.method: {m} (required when verified)")
    elif c["ok"] is not None or c["method"] is not None:
        errs.append("correctness: ok and method must be null when not verified")
    return errs


def _check_consistency(rec: dict[str, Any]) -> list[str]:
    errs: list[str] = []
    algo, mode = rec.get("algorithm"), rec.get("mode")
    if algo in ALGORITHM_MODES and mode in MODES and ALGORITHM_MODES[algo] != mode:
        errs.append(f"mode: algorithm {algo!r} must have mode {ALGORITHM_MODES[algo]!r}")
    if mode == "baseline" and isinstance(algo, str) and algo not in BASELINE_ALGORITHMS:
        errs.append(f"mode: 'baseline' is reserved for {list(BASELINE_ALGORITHMS)}")
    if algo in INTEGER_ONLY_ALGORITHMS and rec.get("weight_type") == "double":
        errs.append(f"weight_type: {algo} supports integer weights (uint64) only")

    n, source = rec.get("n"), rec.get("source")
    if _is_int(n) and _is_int(source) and source >= n:
        errs.append(f"source: must be < n ({n}), got {source}")
    reachable = rec.get("reachable")
    if _is_int(n) and _is_int(reachable) and reachable > n:
        errs.append(f"reachable: must be <= n ({n}), got {reachable}")

    if rec.get("status") == "ok":
        errs += [
            f"{f}: must not be null when status is 'ok'"
            for f in ("m", "source", "trial")
            if f in rec and rec[f] is None
        ]
        times = rec.get("times_ns")
        if isinstance(times, dict):
            errs += [
                f"times_ns.{k}: must be measured when status is 'ok'"
                for k in ("sssp", "end_to_end")
                if k in times and times[k] is None
            ]

    corr = rec.get("correctness")
    if isinstance(corr, dict) and corr.get("ok") is False and rec.get("status") != "error":
        errs.append("status: must be 'error' when verification failed (correctness.ok = false)")

    if rec.get("diagnostics_enabled") is False:
        for name in ("ops", "counters"):
            if rec.get(name) is not None:
                errs.append(f"{name}: must be null when diagnostics are disabled (compiled out)")
    if rec.get("isolated") is False and rec.get("peak_rss_bytes") is not None:
        errs.append(
            "peak_rss_bytes: must be null unless the process was isolated "
            "(one algorithm per process)"
        )
    return errs


def validate_record(rec: Any) -> list[str]:
    """All schema violations of one record (empty list = valid)."""
    if not isinstance(rec, dict):
        return [f"record: must be a JSON object, got {type(rec).__name__}"]
    errs = _check_fields(rec, REQUIRED_FIELDS)
    errs += [f"{k}: unknown field" for k in rec if k not in FIELDS]
    errs += _check_nested(rec)
    errs += _check_consistency(rec)
    return errs


def validate_env(env: Any) -> list[str]:
    """Violations of the environment subset written by envinfo.py."""
    if not isinstance(env, dict):
        return ["env: must be a JSON object"]
    errs = _check_fields(env, ENV_FIELDS)
    probe = {k: env.get(k) for k in ("cpu", "power")}
    errs += _check_nested(probe)
    return errs


def load_jsonl(path: Path | str) -> tuple[list[dict[str, Any]], list[str]]:
    """Valid records of a JSONL file and ``path:line: message`` errors for the rest."""
    path = Path(path)
    records: list[dict[str, Any]] = []
    errors: list[str] = []
    with path.open(encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip():
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError as exc:
                errors.append(f"{path}:{lineno}: invalid JSON: {exc.msg}")
                continue
            errs = validate_record(rec)
            if errs:
                errors += [f"{path}:{lineno}: {e}" for e in errs]
            else:
                records.append(rec)
    return records, errors


def make_failure_record(
    env: dict[str, Any],
    *,
    run_id: str,
    algorithm: str,
    graph_type: str,
    n: int,
    avg_degree: float | None,
    weight_dist: str,
    weight_type: str,
    seed: int,
    status: str,
    error: str,
    order_seed: int | None = None,
    isolated: bool = True,
    params: dict[str, Any] | None = None,
    mode: str | None = None,
    timestamp: str | None = None,
) -> dict[str, Any]:
    """Record for a run that produced no measurement (timeout, memory limit, crash)."""
    if status == "ok":
        raise ValueError("a failure record cannot have status 'ok'")
    if mode is None:
        mode = ALGORITHM_MODES.get(algorithm)
        if mode is None:
            raise ValueError(f"unknown algorithm {algorithm!r}: pass mode explicitly")
    rec: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "run_id": run_id,
        "timestamp": timestamp or _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        **{k: env.get(k) for k in ENV_FIELDS},
        "algorithm": algorithm,
        "mode": mode,
        "params": params or {},
        "graph_type": graph_type,
        "n": n,
        "m": None,
        "avg_degree": avg_degree,
        "weight_dist": weight_dist,
        "weight_type": weight_type,
        "seed": seed,
        "source": None,
        "reachable": None,
        "trial": None,
        "warmup": False,
        "order_seed": order_seed,
        "isolated": isolated,
        "status": status,
        "error": error,
        "times_ns": dict.fromkeys(TIME_KEYS),
        "peak_rss_bytes": None,
        "rss_after_graph_bytes": None,
        "ops": None,
        "counters": None,
        "fallback": None,
        "correctness": {"verified": False, "ok": None, "method": None},
    }
    return {k: rec[k] for k in REQUIRED_FIELDS}


# ---------------------------------------------------------------- CLI
def _cmd_validate(args: argparse.Namespace) -> int:
    total, bad = 0, 0
    for p in args.files:
        records, errors = load_jsonl(p)
        total += len(records)
        bad += len(errors)
        for e in errors:
            print(e, file=sys.stderr)
    print(f"{total} valid record(s), {bad} error(s)")
    return 1 if bad else 0


def _parse_param(kv: str) -> tuple[str, Any]:
    """``key=value``; the value is read as JSON when possible (``4`` -> 4), else kept as text."""
    key, _, raw = kv.partition("=")
    try:
        return key, json.loads(raw)
    except json.JSONDecodeError:
        return key, raw


def _cmd_failure(args: argparse.Namespace) -> int:
    env = json.loads(Path(args.env_json).read_text(encoding="utf-8"))
    params = dict(_parse_param(kv) for kv in args.param)
    rec = make_failure_record(
        env,
        run_id=args.run_id,
        algorithm=args.algorithm,
        graph_type=args.graph_type,
        n=args.n,
        avg_degree=args.avg_degree,
        weight_dist=args.weights,
        weight_type=args.weight_type,
        seed=args.seed,
        status=args.status,
        error=args.error,
        order_seed=args.order_seed,
        params=params,
    )
    errs = validate_record(rec)
    if errs:
        print("\n".join(f"failure record invalid: {e}" for e in errs), file=sys.stderr)
        return 1
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="bench.schema", description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate", help="validate JSONL files")
    v.add_argument("files", nargs="+")
    v.set_defaults(func=_cmd_validate)
    f = sub.add_parser("failure", help="append a failure record (timeout / memory / crash)")
    f.add_argument("--env-json", required=True)
    f.add_argument("--out", required=True)
    f.add_argument("--run-id", required=True)
    f.add_argument("--algorithm", required=True)
    f.add_argument("--graph-type", required=True)
    f.add_argument("--n", type=int, required=True)
    f.add_argument("--avg-degree", type=float)
    f.add_argument("--weights", required=True)
    f.add_argument("--weight-type", required=True, choices=WEIGHT_TYPES)
    f.add_argument("--seed", type=int, required=True)
    f.add_argument("--order-seed", type=int)
    f.add_argument("--param", action="append", default=[], metavar="KEY=VALUE")
    f.add_argument("--status", required=True, choices=[s for s in STATUSES if s != "ok"])
    f.add_argument("--error", required=True)
    f.set_defaults(func=_cmd_failure)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
