#!/usr/bin/env bash
# usage: scripts/ci/changed_tidy.sh <base>
#
# Runs clang-tidy on the src/ and tools/ .cpp files changed in <base>..HEAD (deleted files
# skipped), with the same settings as `make tidy`: compile database build/$PRESET (default
# debug), --quiet, and the macOS SDK sysroot on Darwin. Use it before integrating into main
# and before pushing, so the CI lint job never turns red unnoticed.
#
# Environment: PRESET (default debug), CLANG_TIDY, LLVM_BIN, JOBS (parallel processes).
# Exit code: 0 clean or nothing to check, 1 clang-tidy reported errors, 2 setup problem.
set -euo pipefail

[ $# -eq 1 ] || { echo "usage: scripts/ci/changed_tidy.sh <base>" >&2; exit 2; }
base=$1
cd "$(git rev-parse --show-toplevel)"
git rev-parse --verify -q "${base}^{commit}" >/dev/null || { echo "not a commit: $base" >&2; exit 2; }

# Same tool lookup as the Makefile.
if [ -z "${LLVM_BIN:-}" ]; then
  for d in /opt/homebrew/opt/llvm/bin /usr/local/opt/llvm/bin; do
    if [ -d "$d" ]; then LLVM_BIN=$d; break; fi
  done
fi
CLANG_TIDY=${CLANG_TIDY:-${LLVM_BIN:+$LLVM_BIN/}clang-tidy}
extra=()
if [ "$(uname -s)" = Darwin ]; then extra=("--extra-arg=-isysroot$(xcrun --show-sdk-path)"); fi
PRESET=${PRESET:-debug}
JOBS=${JOBS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}

files=$(git diff --name-only --diff-filter=d "$base" HEAD -- 'src/*.cpp' 'tools/*.cpp')
if [ -z "$files" ]; then
  echo "clang-tidy: no changed src/ or tools/ .cpp since ${base}"
  exit 0
fi
if [ ! -f "build/${PRESET}/compile_commands.json" ]; then
  echo "build/${PRESET}/compile_commands.json missing: run 'cmake --preset ${PRESET}' first" >&2
  exit 2
fi
command -v "$CLANG_TIDY" >/dev/null || { echo "clang-tidy not found: $CLANG_TIDY" >&2; exit 2; }

n=$(printf '%s\n' "$files" | wc -l | tr -d ' ')
echo "clang-tidy: ${n} changed file(s) since ${base} (build/${PRESET})"
if printf '%s\n' "$files" \
  | xargs -P "$JOBS" -n 2 "$CLANG_TIDY" -p "build/${PRESET}" --quiet ${extra[@]+"${extra[@]}"}; then
  echo "clang-tidy: clean"
else
  echo "clang-tidy: FAILED (see errors above)" >&2
  exit 1
fi
