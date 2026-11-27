#!/usr/bin/env bash
# Source-based coverage (clang): configure + build the `coverage` preset, run ctest with
# LLVM_PROFILE_FILE, merge with llvm-profdata and summarise with llvm-cov.
#
#   scripts/coverage.sh           # text summary -> build/coverage/coverage/summary.txt
#   scripts/coverage.sh --html    # additionally an HTML report in build/coverage/coverage/html
#   scripts/coverage.sh --lcov    # additionally build/coverage/coverage/coverage.lcov
#
# Tool lookup (first hit wins): $LLVM_PROFDATA/$LLVM_COV; the directory of the compiler in
# CMakeCache.txt; `xcrun` for Apple clang; versioned llvm-*-<major> on Linux; Homebrew llvm
# (/opt/homebrew or /usr/local); PATH. The tools must match the compiler's LLVM version,
# because the .profraw format is version specific.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PRESET=coverage
BUILD="$ROOT/build/$PRESET"
OUT="$BUILD/coverage"
HTML=0
LCOV=0
for arg in "$@"; do
  case "$arg" in
    --html) HTML=1 ;;
    --lcov) LCOV=1 ;;
    -h | --help) sed -n '2,13p' "$0"; exit 0 ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

cd "$ROOT"
cmake --preset "$PRESET"
cmake --build --preset "$PRESET"

cache_get() { sed -n "s/^$1:[A-Z]*=//p" "$BUILD/CMakeCache.txt" | head -n1; }
CXX_PATH="$(cache_get CMAKE_CXX_COMPILER)"
CXX_ID="$(cache_get CMAKE_CXX_COMPILER_ID)"
if [ -z "$CXX_ID" ]; then
  # CMAKE_CXX_COMPILER_ID is not a cache entry; ask the compiler instead.
  if "$CXX_PATH" --version 2>/dev/null | grep -q 'Apple clang'; then CXX_ID=AppleClang
  elif "$CXX_PATH" --version 2>/dev/null | grep -qi 'clang'; then CXX_ID=Clang
  else CXX_ID=Other; fi
fi
if [ "$CXX_ID" = "Other" ]; then
  echo "coverage preset needs clang (source-based coverage); compiler is $CXX_PATH" >&2
  exit 1
fi

find_tool() { # find_tool <llvm-profdata|llvm-cov>
  local name="$1" dir major cand
  dir="$(dirname "$(readlink -f "$CXX_PATH" 2>/dev/null || echo "$CXX_PATH")")"
  if [ "$CXX_ID" = "AppleClang" ] && command -v xcrun >/dev/null; then
    xcrun --find "$name" 2>/dev/null && return 0
  fi
  [ -x "$dir/$name" ] && { echo "$dir/$name"; return 0; }
  major="$("$CXX_PATH" --version | sed -n 's/.*clang version \([0-9]*\).*/\1/p' | head -n1)"
  if [ -n "$major" ]; then
    for cand in "$name-$major" "/usr/lib/llvm-$major/bin/$name"; do
      command -v "$cand" >/dev/null && { command -v "$cand"; return 0; }
    done
  fi
  for cand in /opt/homebrew/opt/llvm/bin/$name /usr/local/opt/llvm/bin/$name; do
    [ -x "$cand" ] && { echo "$cand"; return 0; }
  done
  command -v "$name" && return 0
  return 1
}

PROFDATA="${LLVM_PROFDATA:-$(find_tool llvm-profdata || true)}"
COV="${LLVM_COV:-$(find_tool llvm-cov || true)}"
if [ -z "$PROFDATA" ] || [ -z "$COV" ]; then
  echo "llvm-profdata / llvm-cov not found (set LLVM_PROFDATA and LLVM_COV)" >&2
  exit 1
fi
echo "compiler: $CXX_PATH ($CXX_ID)"
echo "llvm-profdata: $PROFDATA"
echo "llvm-cov: $COV"

rm -rf "$OUT"
mkdir -p "$OUT/profraw"
# %p = pid, %m = binary signature: parallel ctest processes never clobber each other.
LLVM_PROFILE_FILE="$OUT/profraw/%p-%m.profraw" ctest --preset "$PRESET" --no-tests=ignore

shopt -s nullglob
raws=("$OUT"/profraw/*.profraw)
if [ "${#raws[@]}" -eq 0 ]; then
  echo "no .profraw files were produced (no tests?)" >&2
  exit 1
fi
"$PROFDATA" merge -sparse "${raws[@]}" -o "$OUT/merged.profdata"

# Instrumented binaries: test executables and tools.
objects=()
for bin in "$BUILD"/tests/sssp_*_tests "$BUILD"/*; do
  [ -f "$bin" ] && [ -x "$bin" ] || continue
  case "$(basename "$bin")" in *.sh | *.cmake | *.txt) continue ;; esac
  objects+=("$bin")
done
if [ "${#objects[@]}" -eq 0 ]; then
  echo "no instrumented binaries found under $BUILD" >&2
  exit 1
fi
obj_args=("${objects[0]}")
for o in "${objects[@]:1}"; do obj_args+=(-object "$o"); done

IGNORE='(_deps|/tests/|/benchmarks/|/build/)'
"$COV" report "${obj_args[@]}" -instr-profile="$OUT/merged.profdata" \
  -ignore-filename-regex="$IGNORE" | tee "$OUT/summary.txt"

if [ "$HTML" -eq 1 ]; then
  "$COV" show "${obj_args[@]}" -instr-profile="$OUT/merged.profdata" \
    -ignore-filename-regex="$IGNORE" -format=html -output-dir="$OUT/html" -show-line-counts-or-regions
  echo "HTML report: $OUT/html/index.html"
fi
if [ "$LCOV" -eq 1 ]; then
  "$COV" export "${obj_args[@]}" -instr-profile="$OUT/merged.profdata" \
    -ignore-filename-regex="$IGNORE" -format=lcov > "$OUT/coverage.lcov"
  echo "lcov: $OUT/coverage.lcov"
fi
echo "summary: $OUT/summary.txt"
