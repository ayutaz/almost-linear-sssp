#!/usr/bin/env bash
# Run the SSSP benchmark suite (docs/design/benchmark_protocol.md).
#
#   scripts/run_benchmarks.sh smoke|standard|full [--allow-missing-runner] [--no-build]
#
# Builds the `release` (timing) and `release-diag` (counters) presets, then launches
# build/<preset>/sssp_bench_runner ONCE PER JOB (one build x graph x algorithm, so that
# Peak RSS belongs to one algorithm), in a seeded random order, under a timeout and an
# optional memory limit. Timeouts / crashes are appended as failure records.
# Output: results/raw/<run_id>.jsonl (+ <run_id>.<preset>.env*.json). Runner CLI:
# docs/design/benchmark_runner_cli.md.
#
# Exit codes: 0 done (timeouts / memory limits are recorded, not fatal), 1 a job crashed or
# the output failed schema validation or the build failed, 2 usage error or runner missing.
#
# Environment overrides:
#   BENCH_TIMEOUT        seconds per job (default: profile value; smoke 60, standard 600)
#   BENCH_MEM_LIMIT_MB   address-space limit per job via `ulimit -v` (not enforceable on macOS)
#   BENCH_ORDER_SEED     job/algorithm order seed (default 12345)
#   BENCH_ALGORITHMS, BENCH_GRAPH_TYPES, BENCH_SIZES, BENCH_BUILDS   comma-separated filters
#   BENCH_PARAMS         comma-separated key=value passed as --param to non-baseline algorithms
#   BENCH_COOLDOWN_SEC   pause between jobs (default 0)
#   BENCH_RUN_ID, BENCH_RESULTS_DIR, BENCH_BUILD_ROOT, BENCH_NO_BUILD=1, BENCH_PYTHON
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
log() { printf '[run_benchmarks] %s\n' "$*"; }
err() { printf '[run_benchmarks] %s\n' "$*" >&2; }
usage() {
  sed -n '2,4p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' >&2
  exit 2
}

PROFILE=""
ALLOW_MISSING=0
NO_BUILD="${BENCH_NO_BUILD:-0}"
for arg in "$@"; do
  case "$arg" in
    smoke | standard | full) PROFILE="$arg" ;;
    --allow-missing-runner) ALLOW_MISSING=1 ;;
    --no-build) NO_BUILD=1 ;;
    -h | --help) usage ;;
    *) err "不明な引数: $arg"; usage ;;
  esac
done
[ -n "$PROFILE" ] || usage
if [ "$ALLOW_MISSING" = 1 ] && [ "$PROFILE" != smoke ]; then
  err "--allow-missing-runner は smoke 専用です(CI 用)。"
  exit 2
fi

BUILD_ROOT="${BENCH_BUILD_ROOT:-$ROOT/build}"
RESULTS_DIR="${BENCH_RESULTS_DIR:-$ROOT/results}"
RAW_DIR="$RESULTS_DIR/raw"
ORDER_SEED="${BENCH_ORDER_SEED:-12345}"
BUILDS="${BENCH_BUILDS:-release,release-diag}"
COOLDOWN="${BENCH_COOLDOWN_SEC:-0}"
# shellcheck disable=SC2206  # intentional word splitting: "uv run --quiet python"
PY=(${BENCH_PYTHON:-uv run --quiet python})
export PYTHONPATH="$ROOT/scripts${PYTHONPATH:+:$PYTHONPATH}"
JOBS="$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)"

TIMEOUT="${BENCH_TIMEOUT:-$("${PY[@]}" -m bench.plan "$PROFILE" --order-seed 0 --timeout-sec)}"
case "$TIMEOUT" in '' | *[!0-9]*) err "BENCH_TIMEOUT は正の整数(秒): '$TIMEOUT'"; exit 2 ;; esac

# ---------------------------------------------------------------- build
IFS=',' read -r -a PRESETS <<<"$BUILDS"
if [ "$NO_BUILD" != 1 ]; then
  for preset in "${PRESETS[@]}"; do
    log "build: $preset"
    if ! { cmake --preset "$preset" >/dev/null && cmake --build --preset "$preset" -j "$JOBS"; }; then
      err "ビルドに失敗しました: preset=$preset"
      exit 1
    fi
  done
fi

# ---------------------------------------------------------------- runner presence
missing=0
for preset in "${PRESETS[@]}"; do
  runner="$BUILD_ROOT/$preset/sssp_bench_runner"
  if [ ! -x "$runner" ]; then
    err "sssp_bench_runner が見つかりません: $runner"
    missing=1
  fi
done
if [ "$missing" = 1 ]; then
  err "C++ ベンチランナーが未実装、または未ビルドです。CLI 仕様: docs/design/benchmark_runner_cli.md"
  if [ "$ALLOW_MISSING" = 1 ]; then
    log "smoke: ベンチランナー未実装のため計測をスキップしました(--allow-missing-runner)。結果は作成していません。"
    exit 0
  fi
  exit 2
fi

# ---------------------------------------------------------------- run setup
SHA="$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo nogit)"
RUN_ID="${BENCH_RUN_ID:-$PROFILE-$(date -u +%Y%m%dT%H%M%SZ)-$SHA}"
OUT="$RAW_DIR/$RUN_ID.jsonl"
mkdir -p "$RAW_DIR"
for preset in "${PRESETS[@]}"; do
  "${PY[@]}" -m bench.envinfo --build-dir "$BUILD_ROOT/$preset" \
    --out "$RAW_DIR/$RUN_ID.$preset.env.json" || { err "envinfo に失敗しました"; exit 1; }
done

MEM_KB=""
if [ -n "${BENCH_MEM_LIMIT_MB:-}" ]; then
  MEM_KB=$((BENCH_MEM_LIMIT_MB * 1024))
  if ! (ulimit -v "$MEM_KB") 2>/dev/null; then
    err "警告: この OS では ulimit -v (RLIMIT_AS) を設定できません。メモリ上限は強制されず、"
    err "      終了後に peak_rss_bytes と ${BENCH_MEM_LIMIT_MB} MiB を比較して警告するだけです。"
    MEM_KB=""
    MEM_POSTCHECK=1
  fi
fi

if command -v timeout >/dev/null 2>&1; then
  TIMEOUT_CMD=(timeout "$TIMEOUT")
elif command -v gtimeout >/dev/null 2>&1; then
  TIMEOUT_CMD=(gtimeout "$TIMEOUT")
else # macOS default: SIGALRM via perl (exit status 142)
  TIMEOUT_CMD=(perl -e 'alarm shift; exec @ARGV or die "exec failed: $!\n"' "$TIMEOUT")
fi

PARAM_ARGS=()
if [ -n "${BENCH_PARAMS:-}" ]; then
  IFS=',' read -r -a _params <<<"$BENCH_PARAMS"
  for kv in "${_params[@]}"; do PARAM_ARGS+=(--param "$kv"); done
fi

plan_args=("$PROFILE" --order-seed "$ORDER_SEED" --builds "$BUILDS")
[ -n "${BENCH_ALGORITHMS:-}" ] && plan_args+=(--algorithms "$BENCH_ALGORITHMS")
[ -n "${BENCH_GRAPH_TYPES:-}" ] && plan_args+=(--graph-types "$BENCH_GRAPH_TYPES")
[ -n "${BENCH_SIZES:-}" ] && plan_args+=(--sizes "$BENCH_SIZES")
PLAN_FILE="$RAW_DIR/$RUN_ID.plan.tsv"
"${PY[@]}" -m bench.plan "${plan_args[@]}" >"$PLAN_FILE" || exit 2
TOTAL="$(wc -l <"$PLAN_FILE" | tr -d ' ')"
log "profile=$PROFILE run_id=$RUN_ID jobs=$TOTAL timeout=${TIMEOUT}s order_seed=$ORDER_SEED"
log "output: $OUT"

n_ok=0 n_timeout=0 n_mem=0 n_err=0 i=0
# ---------------------------------------------------------------- jobs
while IFS=$'\t' read -r build gt n deg shape wd wt algo trials warmup sources seed; do
  i=$((i + 1))
  runner="$BUILD_ROOT/$build/sssp_bench_runner"
  env_json="$RAW_DIR/$RUN_ID.$build.env.json"
  extra=()
  case "$algo" in dijkstra-* | bellman-ford) ;; *) extra=("${PARAM_ARGS[@]+"${PARAM_ARGS[@]}"}") ;; esac
  shape_args=() deg_args=()
  if [ "$shape" != - ]; then
    read -r -a shape_args <<<"$shape" # e.g. "--avg-degree 4", "--degree 4", "--attach 2"
    deg_args=(--avg-degree "$deg")
  fi
  cmd=("$runner" --type "$gt" --n "$n" "${shape_args[@]+"${shape_args[@]}"}" --weights "$wd"
    --weight-type "$wt" --seed "$seed" --sources "$sources" --algorithms "$algo"
    --trials "$trials" --warmup "$warmup" --order-seed "$ORDER_SEED"
    --timeout-sec "$TIMEOUT" --verify reference --env-json "$env_json"
    --run-id "$RUN_ID" --out "$OUT" --isolate "${extra[@]+"${extra[@]}"}")
  log "[$i/$TOTAL] $build $algo $gt n=$n d=$deg $wd/$wt"
  (
    [ -n "$MEM_KB" ] && ulimit -v "$MEM_KB"
    exec "${TIMEOUT_CMD[@]}" "${cmd[@]}"
  ) </dev/null
  rc=$?
  status="" msg=""
  case "$rc" in
    0) n_ok=$((n_ok + 1)) ;;
    124 | 142) status=timeout msg="timeout after ${TIMEOUT} s (exit $rc)"; n_timeout=$((n_timeout + 1)) ;;
    5) status=memory_limit msg="runner reported out of memory (exit 5)"; n_mem=$((n_mem + 1)) ;;
    137) status=error msg="killed by SIGKILL (exit 137; possibly the OOM killer)"; n_err=$((n_err + 1)) ;;
    *) status=error msg="runner failed (exit $rc)"; n_err=$((n_err + 1)) ;;
  esac
  if [ -n "$status" ]; then
    err "  -> $status: $msg"
    "${PY[@]}" -m bench.schema failure --env-json "$env_json" --out "$OUT" --run-id "$RUN_ID" \
      --algorithm "$algo" --graph-type "$gt" --n "$n" "${deg_args[@]+"${deg_args[@]}"}" --weights "$wd" \
      --weight-type "$wt" --seed "$seed" --order-seed "$ORDER_SEED" \
      --status "$status" --error "$msg" "${extra[@]+"${extra[@]}"}" ||
      { err "失敗レコードを書けませんでした"; n_err=$((n_err + 1)); }
  fi
  [ "$COOLDOWN" != 0 ] && sleep "$COOLDOWN"
done <"$PLAN_FILE"

for preset in "${PRESETS[@]}"; do # state at the end (thermal limit, load) for comparison
  "${PY[@]}" -m bench.envinfo --build-dir "$BUILD_ROOT/$preset" \
    --out "$RAW_DIR/$RUN_ID.$preset.env-end.json" || true
done

rc_final=0
if [ -s "$OUT" ]; then
  "${PY[@]}" -m bench.schema validate "$OUT" || { err "スキーマ検証に失敗しました: $OUT"; rc_final=1; }
  if [ "${MEM_POSTCHECK:-0}" = 1 ]; then
    "${PY[@]}" - "$OUT" "$BENCH_MEM_LIMIT_MB" <<'EOF'
import json, sys
limit = int(sys.argv[2]) * 1024 * 1024
for line in open(sys.argv[1], encoding="utf-8"):
    r = json.loads(line)
    if (r.get("peak_rss_bytes") or 0) > limit:
        print(f"[run_benchmarks] 警告: peak_rss_bytes={r['peak_rss_bytes']} が上限を超過: "
              f"{r['algorithm']} {r['graph_type']} n={r['n']}", file=sys.stderr)
EOF
  fi
else
  log "レコードはありません: $OUT"
fi
log "done: ok=$n_ok timeout=$n_timeout memory_limit=$n_mem error=$n_err / $TOTAL"
[ "$n_err" -gt 0 ] && rc_final=1
exit "$rc_final"
