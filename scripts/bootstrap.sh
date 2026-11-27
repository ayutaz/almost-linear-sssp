#!/usr/bin/env bash
# Verify / install the local toolchain, then `uv sync` and configure the debug preset.
# Idempotent: safe to run repeatedly.
#
#   scripts/bootstrap.sh            # install what is missing, uv sync, cmake --preset debug
#   scripts/bootstrap.sh --check    # report only, change nothing (exit 1 if something is missing)
#
# Required: clang++, cmake >= 3.24, ninja, clang-format, clang-tidy, uv. Optional: gh.
# Platforms: macOS (Homebrew; prefix $HOMEBREW_PREFIX, else /opt/homebrew on arm64 and
# /usr/local on x86_64) and apt-based Linux (LLVM packages of major $LLVM_MAJOR, default 22).
# Tools are found by command, however they were installed. Only absent tools are installed;
# an old cmake is reported, never replaced.
set -euo pipefail

CHECK_ONLY=0
case "${1:-}" in
  --check) CHECK_ONLY=1 ;;
  -h | --help) sed -n '2,12p' "$0"; exit 0 ;;
  "") ;;
  *) echo "unknown option: $1" >&2; exit 2 ;;
esac

MIN_CMAKE="3.24"
LLVM_MAJOR="${LLVM_MAJOR:-22}"
OS="$(uname -s)"
ARCH="$(uname -m)"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

MISSING=() # required tools that are absent or too old
ABSENT=()  # required tools that are absent (candidates for installation)
HINTS=()   # how to fix what is missing
FOUND=""   # path found by the last find_tool
CXX_FOUND=""
CMAKE_FOUND=""

info() { printf '\033[1;34m==>\033[0m %s\n' "$*"; }
ok() { printf '  \033[32mok\033[0m  %s\n' "$*"; }
miss() { printf '  \033[33mmissing\033[0m  %s\n' "$*"; MISSING+=("$1"); }
warn() { printf '  \033[33mwarning\033[0m  %s\n' "$*"; }

version_ge() { # version_ge A B -> true if A >= B
  [ "$(printf '%s\n%s\n' "$2" "$1" | sort -V | head -n1)" = "$2" ]
}

find_tool() { # find_tool <candidate>... -> FOUND = first usable candidate (name or path)
  local c p
  FOUND=""
  for c in "$@"; do
    if p="$(command -v "$c")" && "$p" --version >/dev/null 2>&1; then
      FOUND="$p"
      return 0
    fi
  done
  return 1
}

require() { # require <logical-name> <candidate>...
  local name="$1"
  shift
  if find_tool "$@"; then
    ok "$name ($FOUND)"
  else
    miss "$name"
    ABSENT+=("$name")
  fi
}

check_cmake() {
  CMAKE_FOUND=""
  if ! find_tool cmake; then
    miss cmake
    ABSENT+=(cmake)
    return 0
  fi
  local v
  v="$("$FOUND" --version | head -n1 | awk '{print $3}')"
  if version_ge "$v" "$MIN_CMAKE"; then
    ok "cmake $v ($FOUND)"
    CMAKE_FOUND="$FOUND"
  else
    miss cmake "(found $v at $FOUND, need >= $MIN_CMAKE)"
    HINTS+=("cmake $v is older than $MIN_CMAKE; install a newer one (e.g. 'uv tool install cmake' or https://cmake.org/download/) and put it first on PATH")
  fi
}

# Candidate commands per tool; set by the platform setup below.
CXX_CANDIDATES=()
FORMAT_CANDIDATES=()
TIDY_CANDIDATES=()

check_toolchain() {
  MISSING=()
  ABSENT=()
  HINTS=()
  if require clang++ "${CXX_CANDIDATES[@]}"; then CXX_FOUND="$FOUND"; fi
  check_cmake
  require ninja ninja
  require clang-format "${FORMAT_CANDIDATES[@]}"
  require clang-tidy "${TIDY_CANDIDATES[@]}"
  require uv uv
  if find_tool gh; then
    ok "gh ($FOUND)"
  else
    warn "gh not found (optional: only needed to watch CI with 'gh run watch'; https://cli.github.com)"
  fi
}

# ------------------------------------------------------------------ macOS
BREW=""
LLVM_PREFIX=""

setup_macos() {
  local default_prefix
  case "$ARCH" in
    arm64) default_prefix=/opt/homebrew ;;
    x86_64)
      default_prefix=/usr/local
      if [ "$(sysctl -in sysctl.proc_translated 2>/dev/null || echo 0)" = "1" ]; then
        echo "warning: this shell runs under Rosetta 2 on Apple Silicon; open a native arm64 shell" >&2
      fi
      ;;
    *) echo "unsupported macOS architecture: $ARCH" >&2; exit 1 ;;
  esac
  local prefix="${HOMEBREW_PREFIX:-$default_prefix}"
  info "macOS $ARCH (Homebrew prefix $prefix)"

  BREW="$prefix/bin/brew"
  if [ ! -x "$BREW" ]; then
    if ! BREW="$(command -v brew)"; then
      echo "Homebrew not found. Install it first: https://brew.sh" >&2
      exit 1
    fi
  fi
  # Homebrew llvm is keg-only: it is deliberately not put on PATH so that Apple clang stays the
  # default compiler. The Makefile and scripts/coverage.sh look it up by path.
  LLVM_PREFIX="$("$BREW" --prefix llvm 2>/dev/null || echo "$prefix/opt/llvm")"
  CXX_CANDIDATES=(clang++)
  FORMAT_CANDIDATES=("$LLVM_PREFIX/bin/clang-format" clang-format)
  TIDY_CANDIDATES=("$LLVM_PREFIX/bin/clang-tidy" clang-tidy)
}

macos_formulae() { # print the brew formulae that provide the absent tools, deduplicated
  local t f out=()
  for t in ${ABSENT[@]+"${ABSENT[@]}"}; do
    case "$t" in
      clang-format | clang-tidy) f=llvm ;;
      cmake | ninja | uv) f="$t" ;;
      *) continue ;;
    esac
    case " ${out[*]-} " in *" $f "*) ;; *) out+=("$f") ;; esac
  done
  echo "${out[*]-}"
}

hints_macos() {
  local formulae
  formulae="$(macos_formulae)"
  [ -n "$formulae" ] && HINTS+=("run: $BREW install $formulae")
  case " ${ABSENT[*]-} " in
    *" clang++ "*) HINTS+=("clang++ comes with the Xcode Command Line Tools: run 'xcode-select --install'") ;;
  esac
  return 0
}

install_macos() {
  local formulae
  formulae="$(macos_formulae)"
  if [ -n "$formulae" ]; then
    info "brew install $formulae"
    # shellcheck disable=SC2086 # word splitting into formula names is intended
    "$BREW" install $formulae
  fi
}

# ------------------------------------------------------------------ Linux
setup_linux() {
  info "Linux $ARCH (LLVM $LLVM_MAJOR)"
  CXX_CANDIDATES=("clang++-$LLVM_MAJOR" clang++)
  FORMAT_CANDIDATES=("clang-format-$LLVM_MAJOR" clang-format)
  TIDY_CANDIDATES=("clang-tidy-$LLVM_MAJOR" clang-tidy)
}

linux_packages() { # print the apt packages that provide the absent tools
  local t out=()
  for t in ${ABSENT[@]+"${ABSENT[@]}"}; do
    case "$t" in
      clang++) out+=("clang-$LLVM_MAJOR" "llvm-$LLVM_MAJOR") ;;
      clang-format | clang-tidy) out+=("$t-$LLVM_MAJOR") ;;
      cmake) out+=(cmake) ;;
      ninja) out+=(ninja-build) ;;
    esac
  done
  echo "${out[*]-}"
}

hints_linux() {
  local pkgs
  pkgs="$(linux_packages)"
  [ -n "$pkgs" ] && HINTS+=("run: sudo apt-get install -y $pkgs")
  case " ${ABSENT[*]-} " in
    *" uv "*) HINTS+=("install uv: https://docs.astral.sh/uv/getting-started/installation/ (e.g. 'pipx install uv')") ;;
  esac
  return 0
}

install_linux() {
  local pkgs
  pkgs="$(linux_packages)"
  [ -n "$pkgs" ] || return 0
  if ! command -v apt-get >/dev/null; then
    echo "only apt-based distributions are automated; install $pkgs (or equivalents) manually" >&2
    return 0
  fi
  local sudo=""
  [ "$(id -u)" -ne 0 ] && sudo="sudo"
  info "apt-get install $pkgs"
  $sudo apt-get update
  # shellcheck disable=SC2086 # word splitting into package names is intended
  if ! $sudo apt-get install -y --no-install-recommends $pkgs; then
    local unversioned="${pkgs//-$LLVM_MAJOR/}"
    echo "versioned LLVM $LLVM_MAJOR packages unavailable; trying: $unversioned" >&2
    # shellcheck disable=SC2086
    $sudo apt-get install -y --no-install-recommends $unversioned
  fi
}

# ------------------------------------------------------------------ main
case "$OS" in
  Darwin) setup_macos; platform=macos ;;
  Linux) setup_linux; platform=linux ;;
  *) echo "unsupported OS: $OS" >&2; exit 1 ;;
esac

check_toolchain
if [ "${#ABSENT[@]}" -gt 0 ] && [ "$CHECK_ONLY" -eq 0 ]; then
  "install_$platform"
  info "re-checking"
  check_toolchain
fi

if [ "${#MISSING[@]}" -gt 0 ]; then
  "hints_$platform"
  for h in ${HINTS[@]+"${HINTS[@]}"}; do echo "$h"; done
  echo "missing: ${MISSING[*]}"
  exit 1
fi

if [ "$platform" = linux ] && ! command -v clang-format >/dev/null; then
  info "only versioned LLVM tools found; for make, export CLANG_FORMAT=clang-format-$LLVM_MAJOR CLANG_TIDY=clang-tidy-$LLVM_MAJOR"
fi

if [ "$CHECK_ONLY" -eq 1 ]; then
  info "all required tools present"
  exit 0
fi

info "uv sync"
(cd "$ROOT" && uv sync)

# The stack is Clang; without this, CMake would pick the system c++ (GCC on most Linux hosts).
# A CXX set by the user wins. An existing build directory keeps its cached compiler.
export CXX="${CXX:-$CXX_FOUND}"
info "cmake --preset debug (CXX=$CXX)"
(cd "$ROOT" && "$CMAKE_FOUND" --preset debug)

info "done. Next: make build && make test"
