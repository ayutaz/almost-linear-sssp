"""Environment capture (scripts/bench/envinfo.py).

Parser tests feed SYNTHETIC command output (hand-written samples in the documented formats);
the collect() tests read the real machine and only check field presence and types.
"""

from __future__ import annotations

import json
import platform
import re
import subprocess
import sys

import pytest
from bench import envinfo, schema
from bench_fixtures import SCRIPTS

PMSET_BATT_AC = (
    "Now drawing from 'AC Power'\n"
    " -InternalBattery-0 (id=1234)\t100%; charged; 0:00 remaining present: true\n"
)
PMSET_BATT_BATTERY = (
    "Now drawing from 'Battery Power'\n"
    " -InternalBattery-0 (id=1234)\t57%; discharging; 3:10 remaining present: true\n"
)
PMSET_G = "System-wide power settings:\nCurrently in use:\n powermode            1\n sleep 0\n"
PMSET_G_OLD = "Currently in use:\n lowpowermode         0\n"
PMSET_THERM = "Note: No thermal warning level has been recorded\n"
PMSET_THERM_LIMIT = "CPU Power notify\n\tCPU_Scheduler_Limit \t= 100\n\tCPU_Speed_Limit \t= 72\n"

CMAKE_CACHE = """\
# This is the CMakeCache file.
CMAKE_BUILD_TYPE:STRING=Release
CMAKE_CXX_COMPILER:FILEPATH=/usr/bin/c++
CMAKE_CXX_FLAGS:STRING=
CMAKE_CXX_FLAGS_RELEASE:STRING=-O3 -DNDEBUG
SSSP_DIAGNOSTICS:BOOL=OFF
"""

PROC_CPUINFO = """\
processor\t: 0
model name\t: Synthetic CPU 9000
physical id\t: 0
core id\t\t: 0

processor\t: 1
model name\t: Synthetic CPU 9000
physical id\t: 0
core id\t\t: 0

processor\t: 2
model name\t: Synthetic CPU 9000
physical id\t: 0
core id\t\t: 1
"""

PROC_MEMINFO = "MemTotal:       16303428 kB\nMemFree:         1000000 kB\n"


def test_parse_pmset_batt():
    assert envinfo.parse_pmset_batt(PMSET_BATT_AC) == {"source": "AC", "battery_percent": 100}
    assert envinfo.parse_pmset_batt(PMSET_BATT_BATTERY) == {
        "source": "Battery",
        "battery_percent": 57,
    }
    assert envinfo.parse_pmset_batt("") == {"source": None, "battery_percent": None}


def test_parse_pmset_low_power_mode():
    assert envinfo.parse_pmset_low_power(PMSET_G) is True
    assert envinfo.parse_pmset_low_power(PMSET_G_OLD) is False
    assert envinfo.parse_pmset_low_power("nothing here") is None


def test_parse_pmset_therm():
    assert envinfo.parse_pmset_therm(PMSET_THERM) is None
    assert envinfo.parse_pmset_therm(PMSET_THERM_LIMIT) == 72


def test_parse_cmake_cache():
    info = envinfo.parse_cmake_cache(CMAKE_CACHE)
    assert info["build_type"] == "Release"
    assert info["cxx_compiler"] == "/usr/bin/c++"
    assert info["build_flags"] == "-O3 -DNDEBUG"
    assert info["diagnostics_enabled"] is False
    assert (
        envinfo.parse_cmake_cache(CMAKE_CACHE.replace("=OFF", "=ON"))["diagnostics_enabled"] is True
    )


def test_flags_from_compile_commands_drop_paths_and_io():
    entries = [
        {
            "directory": "/somewhere/build/release",
            "command": "/usr/bin/c++ -I/abs/include -isystem /abs/dep -DSSSP_DIAGNOSTICS=0 "
            "-O3 -DNDEBUG -std=c++20 -Wall -o src/core/x.o -c /abs/src/core/x.cpp",
            "file": "/abs/src/core/x.cpp",
        }
    ]
    flags = envinfo.flags_from_compile_commands(entries)
    assert flags == "-DSSSP_DIAGNOSTICS=0 -O3 -DNDEBUG -std=c++20 -Wall"
    args_form = [{"arguments": ["c++", "-O2", "-c", "a.cpp", "-o", "a.o"], "file": "a.cpp"}]
    assert envinfo.flags_from_compile_commands(args_form) == "-O2"
    assert envinfo.flags_from_compile_commands([]) is None


def test_parse_proc_cpuinfo_and_meminfo():
    cpu = envinfo.parse_proc_cpuinfo(PROC_CPUINFO)
    assert cpu == {
        "brand": "Synthetic CPU 9000",
        "logical_cores": 3,
        "physical_cores": 2,
        "perf_cores": None,
        "efficiency_cores": None,
    }
    assert envinfo.parse_proc_meminfo(PROC_MEMINFO) == 16303428 * 1024


# SYNTHETIC sysctl answers (name -> value); a missing name is a sysctl that does not exist.
SYSCTL_APPLE_SILICON = {
    "machdep.cpu.brand_string": "Synthetic M9",
    "hw.logicalcpu": "16",
    "hw.physicalcpu": "16",
    "hw.nperflevels": "2",
    "hw.perflevel0.physicalcpu": "12",
    "hw.perflevel1.physicalcpu": "4",
}
# One performance level and no hw.perflevel1.*: what a macOS arm64 VM (e.g. a GitHub-hosted
# macos runner) exposes. Regression for CI run 45929909889, where the old assumption
# "every arm64 Mac reports perflevel1" made efficiency_cores null.
SYSCTL_SINGLE_PERFLEVEL = {
    "machdep.cpu.brand_string": "Synthetic M9 (Virtual)",
    "hw.logicalcpu": "3",
    "hw.physicalcpu": "3",
    "hw.nperflevels": "1",
    "hw.perflevel0.physicalcpu": "3",
}
# No perflevel sysctls at all (Intel Macs): the split cannot be determined.
SYSCTL_NO_PERFLEVELS = {
    "machdep.cpu.brand_string": "Synthetic Core i9",
    "hw.logicalcpu": "16",
    "hw.physicalcpu": "8",
}


@pytest.mark.parametrize(
    ("answers", "perf", "efficiency"),
    [
        (SYSCTL_APPLE_SILICON, 12, 4),
        (SYSCTL_SINGLE_PERFLEVEL, 3, 0),
        (SYSCTL_NO_PERFLEVELS, None, None),
        # two levels reported but perflevel1 unreadable: unknown, not zero
        ({**SYSCTL_APPLE_SILICON, "hw.perflevel1.physicalcpu": None}, 12, None),
    ],
    ids=["two-levels", "single-level-vm", "no-perflevels", "perflevel1-unreadable"],
)
def test_darwin_cpu_core_split_from_sysctl(monkeypatch, answers, perf, efficiency):
    monkeypatch.setattr(envinfo.sys, "platform", "darwin")
    monkeypatch.setattr(envinfo, "_sysctl", lambda name: answers.get(name))
    cpu = envinfo._cpu()
    assert cpu["perf_cores"] == perf
    assert cpu["efficiency_cores"] == efficiency
    assert cpu["brand"] == answers["machdep.cpu.brand_string"]


@pytest.fixture
def fake_build_dir(tmp_path):
    d = tmp_path / "release"
    d.mkdir()
    (d / "CMakeCache.txt").write_text(CMAKE_CACHE, encoding="utf-8")
    return d


def test_collect_has_every_env_field_with_valid_types(fake_build_dir):
    env = envinfo.collect(build_dir=fake_build_dir)
    assert set(schema.ENV_FIELDS) <= set(env)
    assert schema.validate_env(env) == []
    assert env["paper_version"] == schema.PAPER_VERSION
    assert re.fullmatch(r"[0-9a-f]{40}", env["git_sha"])
    assert isinstance(env["git_dirty"], bool)
    assert env["build_type"] == "Release"
    assert env["diagnostics_enabled"] is False
    assert env["arch"] == platform.machine()


@pytest.mark.skipif(sys.platform != "darwin", reason="Apple Silicon P/E core query")
def test_collect_reports_perf_levels_on_macos(fake_build_dir):
    env = envinfo.collect(build_dir=fake_build_dir)
    assert env["cpu"]["brand"]
    if platform.machine() == "arm64":
        assert isinstance(env["cpu"]["perf_cores"], int)
        assert isinstance(env["cpu"]["efficiency_cores"], int)
    assert env["power"] is not None
    assert set(env["power"]) == {"source", "battery_percent", "low_power_mode", "cpu_speed_limit"}


def test_collect_without_cmake_cache_leaves_build_fields_null(tmp_path):
    env = envinfo.collect(build_dir=tmp_path / "missing")
    assert env["build_type"] is None
    # an unknown build must not validate: records must not claim a build type we did not see
    assert any("build_type" in e for e in schema.validate_env(env))


def test_cli_writes_json(fake_build_dir, tmp_path):
    out = tmp_path / "env.json"
    res = subprocess.run(
        [
            sys.executable,
            "-m",
            "bench.envinfo",
            "--build-dir",
            str(fake_build_dir),
            "--out",
            str(out),
        ],
        cwd=SCRIPTS,
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stderr
    env = json.loads(out.read_text(encoding="utf-8"))
    assert schema.validate_env(env) == []
