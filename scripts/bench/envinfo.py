"""Capture the benchmark environment as the ``ENV_FIELDS`` subset of a benchmark record.

The runner copies these fields into every record (``sssp_bench_runner --env-json``).
Works on macOS (sysctl, sw_vers, pmset) and Linux (/proc, /sys). A value that cannot be
determined is ``null``. No hostname or user name is recorded.

    PYTHONPATH=scripts python -m bench.envinfo --build-dir build/release --out env.json
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from bench.schema import ENV_FIELDS, PAPER_VERSION

REPO_ROOT = Path(__file__).resolve().parents[2]


def _run(cmd: list[str], cwd: Path | None = None) -> str | None:
    if shutil.which(cmd[0]) is None and not Path(cmd[0]).is_file():
        return None
    try:
        res = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return res.stdout if res.returncode == 0 else None


# ---------------------------------------------------------------- parsers (pure)
def parse_pmset_batt(text: str) -> dict[str, Any]:
    source = None
    if m := re.search(r"Now drawing from '([^']+)'", text):
        source = {"AC Power": "AC", "Battery Power": "Battery", "UPS Power": "UPS"}.get(m[1])
    pct = re.search(r"(\d+)%", text)
    return {"source": source, "battery_percent": int(pct[1]) if pct else None}


def parse_pmset_low_power(text: str) -> bool | None:
    """``lowpowermode N`` (older macOS) or ``powermode N`` (1 = Low Power)."""
    if m := re.search(r"^\s*lowpowermode\s+(\d+)", text, re.M):
        return m[1] != "0"
    if m := re.search(r"^\s*powermode\s+(\d+)", text, re.M):
        return m[1] == "1"
    return None


def parse_pmset_therm(text: str) -> int | None:
    m = re.search(r"CPU_Speed_Limit\s*=\s*(\d+)", text)
    return int(m[1]) if m else None


def parse_cmake_cache(text: str) -> dict[str, Any]:
    entries: dict[str, str] = {}
    for line in text.splitlines():
        if m := re.match(r"^([A-Za-z0-9_]+):[A-Z]+=(.*)$", line):
            entries[m[1]] = m[2].strip()
    build_type = entries.get("CMAKE_BUILD_TYPE") or None
    flags = [entries.get("CMAKE_CXX_FLAGS", "")]
    if build_type:
        flags.append(entries.get(f"CMAKE_CXX_FLAGS_{build_type.upper()}", ""))
    diag = entries.get("SSSP_DIAGNOSTICS")
    return {
        "build_type": build_type,
        "cxx_compiler": entries.get("CMAKE_CXX_COMPILER") or None,
        "build_flags": " ".join(f for f in flags if f) or None,
        "diagnostics_enabled": None if diag is None else diag.upper() in ("ON", "1", "TRUE"),
    }


_PATH_FLAGS = {"-I", "-isystem", "-iquote", "-idirafter", "-o", "-MF", "-MT", "-MQ"}


def flags_from_compile_commands(entries: list[dict[str, Any]]) -> str | None:
    """Effective compile flags of one translation unit, without paths, I/O and the source."""
    if not entries:
        return None
    entry = next((e for e in entries if "/src/" in str(e.get("file", ""))), entries[0])
    args = entry.get("arguments") or shlex.split(entry.get("command", ""))
    src = str(entry.get("file", ""))
    out: list[str] = []
    skip = False
    for a in args[1:]:
        if skip:
            skip = False
            continue
        if a in _PATH_FLAGS:
            skip = True
            continue
        if a in ("-c", src, "-MD", "-MMD") or a.startswith(("-I", "-isystem", "-o", "-MF")):
            continue
        if a.endswith((".cpp", ".cc", ".cxx", ".c")):
            continue
        out.append(a)
    return " ".join(out) or None


def parse_proc_cpuinfo(text: str) -> dict[str, Any]:
    brand = None
    logical = 0
    cores: set[tuple[str, str]] = set()
    phys = core = None
    for line in text.splitlines() + [""]:
        key, _, val = line.partition(":")
        key, val = key.strip(), val.strip()
        if key == "processor":
            logical += 1
        elif key in ("model name", "Model", "Hardware") and brand is None and val:
            brand = val
        elif key == "physical id":
            phys = val
        elif key == "core id":
            core = val
        elif not line.strip() and core is not None:
            cores.add((phys or "0", core))
            phys = core = None
    return {
        "brand": brand,
        "logical_cores": logical or None,
        "physical_cores": len(cores) or None,
        "perf_cores": None,
        "efficiency_cores": None,
    }


def parse_proc_meminfo(text: str) -> int | None:
    m = re.search(r"^MemTotal:\s+(\d+)\s+kB", text, re.M)
    return int(m[1]) * 1024 if m else None


# ---------------------------------------------------------------- collectors
def _sysctl(name: str) -> str | None:
    out = _run(["sysctl", "-n", name])
    return out.strip() if out and out.strip() else None


def _sysctl_int(name: str) -> int | None:
    v = _sysctl(name)
    return int(v) if v and v.isdigit() else None


def _os_string() -> str:
    if sys.platform == "darwin":
        ver = (_run(["sw_vers", "-productVersion"]) or "").strip()
        build = (_run(["sw_vers", "-buildVersion"]) or "").strip()
        if ver:
            return f"macOS {ver} ({build}); Darwin {platform.release()}"
    if sys.platform.startswith("linux"):
        try:
            text = Path("/etc/os-release").read_text(encoding="utf-8")
            if m := re.search(r'^PRETTY_NAME="?([^"\n]+)"?', text, re.M):
                return f"{m[1]}; Linux {platform.release()}"
        except OSError:
            pass
    return f"{platform.system()} {platform.release()}"


def _cpu() -> dict[str, Any]:
    if sys.platform == "darwin":
        return {
            "brand": _sysctl("machdep.cpu.brand_string"),
            "logical_cores": _sysctl_int("hw.logicalcpu"),
            "physical_cores": _sysctl_int("hw.physicalcpu"),
            "perf_cores": _sysctl_int("hw.perflevel0.physicalcpu"),
            "efficiency_cores": _sysctl_int("hw.perflevel1.physicalcpu"),
        }
    try:
        return parse_proc_cpuinfo(Path("/proc/cpuinfo").read_text(encoding="utf-8"))
    except OSError:
        return {
            "brand": platform.processor() or None,
            "logical_cores": os.cpu_count(),
            "physical_cores": None,
            "perf_cores": None,
            "efficiency_cores": None,
        }


def _ram_bytes() -> int | None:
    if sys.platform == "darwin":
        return _sysctl_int("hw.memsize")
    try:
        return parse_proc_meminfo(Path("/proc/meminfo").read_text(encoding="utf-8"))
    except OSError:
        return None


def _power() -> dict[str, Any] | None:
    if sys.platform == "darwin":
        batt = parse_pmset_batt(_run(["pmset", "-g", "batt"]) or "")
        return {
            **batt,
            "low_power_mode": parse_pmset_low_power(_run(["pmset", "-g"]) or ""),
            "cpu_speed_limit": parse_pmset_therm(_run(["pmset", "-g", "therm"]) or ""),
        }
    supplies = Path("/sys/class/power_supply")
    if not supplies.is_dir():
        return None
    source = None
    for online in supplies.glob("*/online"):
        try:
            if online.read_text().strip() == "1":
                source = "AC"
        except OSError:
            continue
    if source is None and any(supplies.glob("BAT*")):
        source = "Battery"
    return {
        "source": source,
        "battery_percent": None,
        "low_power_mode": None,
        "cpu_speed_limit": None,
    }


def _compiler(cache: dict[str, Any], build_dir: Path | None) -> str | None:
    candidates = [cache.get("cxx_compiler")]
    cc_json = build_dir / "compile_commands.json" if build_dir else None
    if cc_json and cc_json.is_file():
        try:
            entries = json.loads(cc_json.read_text(encoding="utf-8"))
            first = entries[0]
            args = first.get("arguments") or shlex.split(first.get("command", ""))
            candidates.append(args[0])
        except (OSError, ValueError, IndexError, KeyError):
            pass
    candidates.append("c++")
    for c in candidates:
        if c and (out := _run([c, "--version"])):
            return out.splitlines()[0].strip()
    return None


def _build_flags(cache: dict[str, Any], build_dir: Path | None) -> str | None:
    cc_json = build_dir / "compile_commands.json" if build_dir else None
    if cc_json and cc_json.is_file():
        try:
            flags = flags_from_compile_commands(json.loads(cc_json.read_text(encoding="utf-8")))
            if flags:
                return flags
        except (OSError, ValueError):
            pass
    return cache.get("build_flags")


def _git(repo: Path) -> tuple[str | None, bool | None]:
    sha = _run(["git", "rev-parse", "HEAD"], cwd=repo)
    status = _run(["git", "status", "--porcelain"], cwd=repo)
    return (sha.strip() if sha else None), (bool(status.strip()) if status is not None else None)


def _load_average() -> list[float] | None:
    try:
        return [round(x, 2) for x in os.getloadavg()]
    except OSError:
        return None


def collect(build_dir: Path | str | None = None, repo: Path | str = REPO_ROOT) -> dict[str, Any]:
    build_dir = Path(build_dir) if build_dir else None
    cache: dict[str, Any] = {}
    if build_dir and (build_dir / "CMakeCache.txt").is_file():
        cache = parse_cmake_cache((build_dir / "CMakeCache.txt").read_text(encoding="utf-8"))
    sha, dirty = _git(Path(repo))
    env = {
        "paper_version": PAPER_VERSION,
        "git_sha": sha,
        "git_dirty": dirty,
        "os": _os_string(),
        "cpu": _cpu(),
        "ram_bytes": _ram_bytes(),
        "arch": platform.machine(),
        "compiler": _compiler(cache, build_dir),
        "build_type": cache.get("build_type"),
        "build_flags": _build_flags(cache, build_dir),
        "diagnostics_enabled": cache.get("diagnostics_enabled"),
        "power": _power(),
        "load_average": _load_average(),
    }
    return {k: env[k] for k in ENV_FIELDS}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="bench.envinfo", description=__doc__.splitlines()[0])
    ap.add_argument("--build-dir", help="CMake build directory (reads CMakeCache.txt)")
    ap.add_argument("--out", help="write JSON here (default: stdout)")
    args = ap.parse_args(argv)
    text = json.dumps(collect(args.build_dir), ensure_ascii=False, indent=2)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
