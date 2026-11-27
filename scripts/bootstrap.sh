#!/usr/bin/env bash
# Install / verify the local toolchain: cmake (>= 3.24), ninja, LLVM (clang-format, clang-tidy,
# llvm-cov, llvm-profdata) and uv, then `uv sync`. Idempotent: safe to run repeatedly.
#
#   scripts/bootstrap.sh            # install what is missing, then uv sync
#   scripts/bootstrap.sh --check    # report only, install nothing (exit 1 if something is missing)
#
# Platforms: macOS Apple Silicon (arm64, Homebrew in /opt/homebrew), macOS Intel (x86_64,
# Homebrew in /usr/local) and Linux (apt; x86_64 or aarch64). Detected with `uname -s/-m`.
set -euo pipefail

CHECK_ONLY=0
case "${1:-}" in
  --check) CHECK_ONLY=1 ;;
  -h | --help) sed -n '2,10p' "$0"; exit 0 ;;
  "") ;;
  *) echo "unknown option: $1" >&2; exit 2 ;;
esac

MIN_CMAKE="3.24"
OS="$(uname -s)"
ARCH="$(uname -m)"
MISSING=()

info() { printf '\033[1;34m==>\033[0m %s\n' "$*"; }
ok() { printf '  \033[32mok\033[0m  %s\n' "$*"; }
miss() { printf '  \033[33mmissing\033[0m  %s\n' "$*"; MISSING+=("$1"); }

version_ge() { # version_ge A B -> true if A >= B
  [ "$(printf '%s\n%s\n' "$2" "$1" | sort -V | head -n1)" = "$2" ]
}

check_cmake() {
  if command -v cmake >/dev/null; then
    local v
    v="$(cmake --version | head -n1 | awk '{print $3}')"
    if version_ge "$v" "$MIN_CMAKE"; then ok "cmake $v"; else miss cmake "(found $v, need >= $MIN_CMAKE)"; fi
  else
    miss cmake
  fi
}

check_cmd() { # check_cmd <logical-name> <command>
  if command -v "$2" >/dev/null; then ok "$1 ($(command -v "$2"))"; else miss "$1"; fi
}

# ------------------------------------------------------------------ macOS
bootstrap_macos() {
  local brew_prefix
  case "$ARCH" in
    arm64) brew_prefix=/opt/homebrew ;;
    x86_64)
      brew_prefix=/usr/local
      if [ "$(sysctl -in sysctl.proc_translated 2>/dev/null || echo 0)" = "1" ]; then
        echo "warning: this shell runs under Rosetta 2 on Apple Silicon; open a native arm64 shell" >&2
      fi
      ;;
    *) echo "unsupported macOS architecture: $ARCH" >&2; exit 1 ;;
  esac
  info "macOS $ARCH (Homebrew prefix $brew_prefix)"

  local brew="$brew_prefix/bin/brew"
  if [ ! -x "$brew" ]; then
    if command -v brew >/dev/null; then brew="$(command -v brew)"; else
      echo "Homebrew not found. Install it first: https://brew.sh" >&2
      exit 1
    fi
  fi

  # A tool counts as present if it is usable, however it was installed (e.g. uv from the
  # official installer in ~/.local/bin). Only tools that are absent are installed with brew.
  local llvm_prefix
  llvm_prefix="$("$brew" --prefix llvm 2>/dev/null || echo "$brew_prefix/opt/llvm")"
  local formulae=(cmake ninja llvm uv) f found
  local to_install=()
  # Presence probes: an absent tool (non-zero status) is expected and leaves `found` empty,
  # which is reported as missing (and installed) below.
  for f in "${formulae[@]}"; do
    case "$f" in
      llvm) found="$([ -x "$llvm_prefix/bin/clang-format" ] && echo "$llvm_prefix/bin" || true)" ;;
      *) found="$(command -v "$f" || true)" ;;
    esac
    if [ -n "$found" ]; then
      ok "$f ($found)"
    else
      miss "$f"
      to_install+=("$f")
    fi
  done
  # A non-Homebrew cmake might shadow brew's; still verify the version that is on PATH.
  if command -v cmake >/dev/null; then check_cmake; fi

  if [ "${#to_install[@]}" -gt 0 ]; then
    if [ "$CHECK_ONLY" -eq 1 ]; then
      echo "run: $brew install ${to_install[*]}"
    else
      info "brew install ${to_install[*]}"
      "$brew" install "${to_install[@]}"
      MISSING=()
    fi
  fi

  if [ -x "$llvm_prefix/bin/clang-format" ]; then
    ok "clang-format/clang-tidy in $llvm_prefix/bin (Makefile picks this path up automatically)"
  fi
  # Homebrew llvm is keg-only: it is deliberately not put on PATH so that Apple clang stays the
  # default compiler. The Makefile and scripts/coverage.sh look it up by path.
}

# ------------------------------------------------------------------ Linux
bootstrap_linux() {
  info "Linux $ARCH"
  if ! command -v apt-get >/dev/null; then
    echo "only apt-based distributions are automated; install cmake>=$MIN_CMAKE, ninja, clang," \
      "clang-format, clang-tidy, llvm and uv manually" >&2
    exit 1
  fi

  # Prefer the same LLVM major as CI (LLVM_MAJOR, default 22); fall back to unversioned packages.
  local major="${LLVM_MAJOR:-22}"
  local pkgs=() p
  for p in cmake ninja-build "clang-$major" "clang-format-$major" "clang-tidy-$major" "llvm-$major"; do
    if dpkg-query -W -f='${Status}' "$p" 2>/dev/null | grep -q 'install ok installed'; then
      ok "$p"
    else
      miss "$p"
      pkgs+=("$p")
    fi
  done

  if [ "${#pkgs[@]}" -gt 0 ]; then
    local sudo=""
    [ "$(id -u)" -ne 0 ] && sudo="sudo"
    if [ "$CHECK_ONLY" -eq 1 ]; then
      echo "run: $sudo apt-get install -y ${pkgs[*]}"
    else
      info "apt-get install ${pkgs[*]}"
      $sudo apt-get update
      if ! $sudo apt-get install -y --no-install-recommends "${pkgs[@]}"; then
        echo "versioned LLVM $major packages unavailable; trying unversioned clang/llvm" >&2
        $sudo apt-get install -y --no-install-recommends cmake ninja-build clang clang-format clang-tidy llvm
      fi
      MISSING=()
    fi
  fi
  check_cmake

  if command -v uv >/dev/null; then
    ok "uv ($(uv --version))"
  else
    miss uv
    # uv is not packaged by apt; do not pipe a remote script into a shell on the user's behalf.
    echo "install uv: https://docs.astral.sh/uv/getting-started/installation/ (e.g. 'pipx install uv')"
  fi
}

case "$OS" in
  Darwin) bootstrap_macos ;;
  Linux) bootstrap_linux ;;
  *) echo "unsupported OS: $OS" >&2; exit 1 ;;
esac

if [ "$CHECK_ONLY" -eq 1 ]; then
  if [ "${#MISSING[@]}" -gt 0 ]; then
    echo "missing: ${MISSING[*]}"
    exit 1
  fi
  info "all tools present"
  exit 0
fi

if command -v uv >/dev/null; then
  info "uv sync"
  (cd "$(dirname "$0")/.." && uv sync)
else
  echo "uv is still missing; install it and re-run this script" >&2
  exit 1
fi
info "done. Next: make build && make test"
