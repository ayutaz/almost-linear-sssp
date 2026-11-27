"""Tests for scripts/bootstrap.sh (toolchain check, install of what is missing, uv sync, configure).

The script runs against a sandbox PATH: a directory of fake tools written by each test plus a
directory of symlinks to the few POSIX utilities the script needs. Nothing real is installed:
the fake ``brew`` / ``apt-get`` only log their arguments and drop further fakes into the sandbox.
``uname`` is faked too, so the macOS and Linux branches are exercised on either host.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "bootstrap.sh"

# Host utilities the script (and the fakes below) may use. Tools under test are never linked.
SYSTEM_UTILS = (
    "awk",
    "cat",
    "chmod",
    "dirname",
    "grep",
    "head",
    "id",
    "mkdir",
    "sed",
    "sort",
)

REQUIRED_LINUX = ("clang++", "cmake", "ninja", "clang-format", "clang-tidy", "uv")


@dataclass
class Sandbox:
    tmp: Path
    os_name: str
    fake: Path
    sysbin: Path
    log: Path
    brew_prefix: Path

    def tool(self, name: str, body: str = "", directory: Path | None = None) -> Path:
        """Write a fake executable; by default it answers --version and logs other calls."""
        target = (directory or self.fake) / name
        target.parent.mkdir(parents=True, exist_ok=True)
        default = (
            'if [ "${1:-}" = --version ]; then echo "' + name + ' version 99.0.0"; exit 0; fi\n'
            'echo "' + name + ' $* @ $PWD" >> "' + str(self.log) + '"\n'
        )
        target.write_text("#!/bin/sh\n" + (body or default))
        target.chmod(0o755)
        return target

    def cmake(self, version: str = "3.28.3", preset_exit: int = 0) -> None:
        self.tool(
            "cmake",
            'if [ "${1:-}" = --version ]; then echo "cmake version ' + version + '"; exit 0; fi\n'
            'echo "cmake $* CXX=${CXX:-} @ $PWD" >> "' + str(self.log) + '"\n'
            f"exit {preset_exit}\n",
        )

    def linux_toolchain(self, *, skip: tuple[str, ...] = (), suffix: str = "") -> None:
        for name in REQUIRED_LINUX + ("gh",):
            if name in skip:
                continue
            if name == "cmake":
                self.cmake()
            elif name in ("clang++", "clang-format", "clang-tidy"):
                self.tool(name + suffix)
            else:
                self.tool(name)

    def run(self, *args: str, **extra_env: str) -> subprocess.CompletedProcess:
        env = {
            "PATH": f"{self.fake}{os.pathsep}{self.sysbin}",
            "HOME": str(self.tmp / "home"),
            "HOMEBREW_PREFIX": str(self.brew_prefix),
            "LLVM_MAJOR": "22",
            **extra_env,
        }
        return subprocess.run(
            [shutil.which("bash") or "/bin/bash", str(SCRIPT), *args],
            env=env,
            capture_output=True,
            text=True,
            timeout=60,
        )

    def calls(self) -> list[str]:
        return self.log.read_text().splitlines() if self.log.exists() else []


def make_sandbox(tmp_path: Path, os_name: str) -> Sandbox:
    fake = tmp_path / "fake"
    sysbin = tmp_path / "sysbin"
    fake.mkdir()
    sysbin.mkdir()
    for util in SYSTEM_UTILS:
        found = shutil.which(util)
        assert found, f"host utility {util} not found"
        (sysbin / util).symlink_to(found)
    box = Sandbox(tmp_path, os_name, fake, sysbin, tmp_path / "calls.log", tmp_path / "brew")
    arch = "arm64" if os_name == "Darwin" else "x86_64"
    box.tool("uname", f'case "$1" in -s) echo {os_name} ;; -m) echo {arch} ;; esac\n')
    # A package manager "install" only creates fakes, mimicking a successful install.
    box.tool("sudo", 'exec "$@"\n')
    box.tool(
        "apt-get",
        'echo "apt-get $*" >> "' + str(box.log) + '"\n'
        'for p in "$@"; do\n'
        '  case "$p" in\n'
        "    ninja-build) t=ninja ;;\n"
        "    cmake) t=cmake ;;\n"
        "    clang-format-*|clang-tidy-*) t=$p ;;\n"
        "    clang-*[0-9]) t=clang++-${p#clang-} ;;\n"
        "    *) continue ;;\n"
        "  esac\n"
        '  printf \'#!/bin/sh\\necho "%s version 99.0.0"\\n\' "$t" > "' + str(fake) + '/$t"\n'
        '  chmod 755 "' + str(fake) + '/$t"\n'
        "done\n",
    )
    llvm_bin = tmp_path / "llvm" / "bin"
    box.tool(
        "brew",
        'case "$1" in\n'
        '  --prefix) echo "' + str(tmp_path / "llvm") + '" ;;\n'
        "  install)\n"
        '    echo "brew $*" >> "' + str(box.log) + '"\n'
        "    shift\n"
        '    for f in "$@"; do\n'
        '      case "$f" in\n'
        '        llvm) d="' + str(llvm_bin) + '"; set -- clang-format clang-tidy ;;\n'
        '        *) d="' + str(fake) + '"; set -- "$f" ;;\n'
        "      esac\n"
        '      mkdir -p "$d"\n'
        '      for t in "$@"; do\n'
        '        printf \'#!/bin/sh\\necho "%s version 99.0.0"\\n\' "$t" > "$d/$t"\n'
        '        chmod 755 "$d/$t"\n'
        "      done\n"
        "    done ;;\n"
        "esac\n",
        directory=box.brew_prefix / "bin",
    )
    return box


@pytest.fixture
def linux(tmp_path: Path) -> Sandbox:
    return make_sandbox(tmp_path, "Linux")


@pytest.fixture
def macos(tmp_path: Path) -> Sandbox:
    return make_sandbox(tmp_path, "Darwin")


def macos_toolchain(box: Sandbox, *, skip: tuple[str, ...] = ()) -> None:
    for name in ("clang++", "ninja", "uv", "gh"):
        if name not in skip:
            box.tool(name)
    if "cmake" not in skip:
        box.cmake()
    llvm_bin = box.tmp / "llvm" / "bin"
    if "llvm" not in skip:
        box.tool("clang-format", directory=llvm_bin)
        box.tool("clang-tidy", directory=llvm_bin)


def out(p: subprocess.CompletedProcess) -> str:
    return p.stdout + p.stderr


# ------------------------------------------------------------------ check-only mode


def test_check_passes_when_every_tool_is_present(linux):
    linux.linux_toolchain()
    p = linux.run("--check")
    assert p.returncode == 0, out(p)
    for name in REQUIRED_LINUX + ("gh",):
        assert name in p.stdout
    assert "cmake 3.28.3" in p.stdout
    # --check reports only: no sync, no configure, no install.
    assert linux.calls() == []


def test_check_accepts_versioned_llvm_tools(linux):
    linux.linux_toolchain(suffix="-22")
    p = linux.run("--check")
    assert p.returncode == 0, out(p)
    assert "clang++-22" in p.stdout


@pytest.mark.parametrize("tool", REQUIRED_LINUX)
def test_check_fails_when_a_required_tool_is_missing(linux, tool):
    linux.linux_toolchain(skip=(tool,))
    p = linux.run("--check")
    assert p.returncode == 1, out(p)
    missing_line = [line for line in p.stdout.splitlines() if line.startswith("missing:")]
    assert missing_line, out(p)
    assert tool in missing_line[0]
    assert linux.calls() == []


def test_check_rejects_cmake_older_than_3_24(linux):
    linux.linux_toolchain()
    linux.cmake(version="3.20.0")
    p = linux.run("--check")
    assert p.returncode == 1, out(p)
    assert "3.20.0" in p.stdout
    assert "3.24" in p.stdout


def test_missing_gh_is_reported_but_optional(linux):
    linux.linux_toolchain(skip=("gh",))
    p = linux.run("--check")
    assert p.returncode == 0, out(p)
    assert "gh" in out(p)
    assert "optional" in out(p)


def test_unknown_option_is_rejected(linux):
    p = linux.run("--bogus")
    assert p.returncode == 2
    assert "unknown option" in p.stderr


# ------------------------------------------------------------------ install mode (Linux)


def test_install_mode_runs_uv_sync_then_configures_debug_preset_with_clang(linux):
    linux.linux_toolchain(suffix="-22")
    p = linux.run()
    assert p.returncode == 0, out(p)
    calls = linux.calls()
    sync = [i for i, c in enumerate(calls) if c.startswith("uv sync")]
    conf = [i for i, c in enumerate(calls) if c.startswith("cmake --preset debug")]
    assert sync and conf, calls
    assert sync[0] < conf[0], calls
    assert calls[sync[0]].endswith(f"@ {ROOT}")
    assert calls[conf[0]].endswith(f"@ {ROOT}")
    # The stack is Clang: an unset CXX must not let CMake fall back to the system c++ (GCC).
    assert f"CXX={linux.fake / 'clang++-22'} " in calls[conf[0]]
    assert not any(c.startswith("apt-get") for c in calls)


def test_install_mode_keeps_a_user_supplied_cxx(linux):
    linux.linux_toolchain()
    p = linux.run(CXX="my-clang++")
    assert p.returncode == 0, out(p)
    assert any("--preset debug CXX=my-clang++ " in c for c in linux.calls()), linux.calls()


def test_install_mode_fails_when_configure_fails(linux):
    linux.linux_toolchain()
    linux.cmake(preset_exit=3)
    p = linux.run()
    assert p.returncode != 0, out(p)
    assert "done. Next" not in p.stdout


def test_install_mode_with_old_cmake_stops_before_sync_and_configure(linux):
    # An old cmake is not replaced behind the user's back (apt's cmake can itself be < 3.24,
    # e.g. Ubuntu 22.04): it is reported with a hint and the run fails.
    linux.linux_toolchain()
    linux.cmake(version="3.20.0")
    p = linux.run()
    assert p.returncode != 0, out(p)
    assert "3.20.0" in out(p)
    assert not any(c.startswith("apt-get") for c in linux.calls())
    assert not any(c.startswith(("uv sync", "cmake --preset")) for c in linux.calls())


def test_install_mode_apt_installs_only_what_is_missing(linux):
    linux.linux_toolchain(skip=("ninja",))
    p = linux.run()
    assert p.returncode == 0, out(p)
    installs = [c for c in linux.calls() if c.startswith("apt-get") and " install " in c]
    assert len(installs) == 1, linux.calls()
    assert installs[0].split()[-1] == "ninja-build"
    assert installs[0].count("ninja-build") == 1
    assert any(c.startswith("cmake --preset debug") for c in linux.calls())


def test_install_mode_without_uv_fails_and_does_not_configure(linux):
    linux.linux_toolchain(skip=("uv",))
    p = linux.run()
    assert p.returncode == 1, out(p)
    assert "docs.astral.sh/uv" in out(p)
    assert not any(c.startswith("cmake --preset") for c in linux.calls())


# ------------------------------------------------------------------ macOS


def test_macos_check_reports_brew_hint_and_never_installs(macos):
    macos_toolchain(macos, skip=("ninja",))
    p = macos.run("--check")
    assert p.returncode == 1, out(p)
    assert "brew install ninja" in p.stdout
    assert not any(c.startswith("brew install") for c in macos.calls())


def test_macos_install_mode_brews_only_missing_formulae_then_configures(macos):
    macos_toolchain(macos, skip=("ninja", "llvm"))
    p = macos.run()
    assert p.returncode == 0, out(p)
    calls = macos.calls()
    assert [c for c in calls if c.startswith("brew install")] == ["brew install ninja llvm"]
    assert any(c.startswith("uv sync") for c in calls), calls
    assert any(c.startswith("cmake --preset debug") for c in calls), calls


def test_macos_missing_compiler_points_to_command_line_tools(macos):
    macos_toolchain(macos, skip=("clang++",))
    p = macos.run("--check")
    assert p.returncode == 1, out(p)
    assert "xcode-select --install" in out(p)
    assert "clang++" in [ln for ln in p.stdout.splitlines() if ln.startswith("missing:")][0]


def test_macos_without_homebrew_fails_with_hint(macos):
    macos_toolchain(macos)
    shutil.rmtree(macos.brew_prefix)
    p = macos.run("--check")
    assert p.returncode == 1, out(p)
    assert "brew.sh" in p.stderr
