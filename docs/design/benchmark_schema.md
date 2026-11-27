# ベンチマーク結果のスキーマ(JSONL, schema_version "1")

ベンチマークの結果は `results/raw/<run_id>.jsonl` に、1 行 1 レコードの JSON で保存する。
検証器は `scripts/bench/schema.py`(`validate_record`)で、この文書と同じ規則を実装している。
計測の手順は [benchmark_protocol.md](benchmark_protocol.md)、ランナーの引数は
[benchmark_runner_cli.md](benchmark_runner_cli.md) を参照する。

```sh
PYTHONPATH=scripts uv run python -m bench.schema validate results/raw/<run_id>.jsonl   # 0: 有効, 1: 違反あり
```

## 1. 全体の規則

- **下表のキーはすべて必須**。値が無いときもキーは省略せず、`null` を書く。表に無いキーは違反とする(誤字の検出のため)。
- **未計測の値は `null`。`0` を入れない。** 時間 `times_ns.*` は 1 以上の整数か `null`。
  `peak_rss_bytes`、`rss_after_graph_bytes`、`ram_bytes` は 1 以上の整数か `null`。
- 整数のフィールドに `true`/`false` や `1000.0` を入れない(JSON の数値でも整数として書く)。
- 1 レコード = 1 回の SSSP 実行(1 アルゴリズム × 1 グラフ × 1 始点 × 1 試行)。
  ただし失敗レコード(§4)は、ランナー 1 プロセス(ジョブ)が何も出力できなかったことを表す。
- ホスト名・ユーザー名・ローカルの絶対パスは記録しない(公開リポジトリのため)。

## 2. フィールド

### 2.1 実行と再現性

| キー | 型 | null | 内容 |
|---|---|---|---|
| `schema_version` | string | 不可 | 常に `"1"` |
| `run_id` | string | 不可 | `scripts/run_benchmarks.sh` の 1 回の実行の ID(例 `standard-20261127T013000Z-d5a55e7`) |
| `timestamp` | string | 不可 | ISO-8601、**UTC オフセット必須**(例 `2026-11-27T10:30:00+09:00`) |
| `paper_version` | string | 不可 | 再現対象の論文版。`2611.34911v1`(形式 `NNNN.NNNNNvN`) |
| `git_sha` | string | 不可 | 計測したソースのコミット(16 進 7〜40 桁) |
| `git_dirty` | bool | 不可 | 作業ツリーに未コミットの変更があったか。`true` の結果は公表値に使わない |

### 2.2 アルゴリズム

| キー | 型 | null | 内容 |
|---|---|---|---|
| `algorithm` | string | 不可 | `dijkstra-binary` / `dijkstra-4ary` / `dijkstra-fibonacci` / `dijkstra-radix` / `bellman-ford` / `paper` / `paper-experimental` / `paper-optimized`(未知の名前も可) |
| `mode` | enum | 不可 | `baseline` / `paper-faithful` / `paper-experimental` / `practical-optimized` |
| `params` | object | 不可 | アルゴリズムのパラメータ(閾値、乱数 seed など)。無ければ `{}` |

既知の名前と `mode` の対応は固定(`schema.ALGORITHM_MODES`): Dijkstra 系と `bellman-ford` は `baseline`、
`paper` は `paper-faithful`、`paper-experimental` は `paper-experimental`、`paper-optimized` は
`practical-optimized`。`baseline` はこの 5 つの基準アルゴリズム専用。

### 2.3 環境(`scripts/bench/envinfo.py` が作り、ランナーが `--env-json` からそのまま写す)

| キー | 型 | null | 内容 |
|---|---|---|---|
| `os` | string | 不可 | 例 `macOS 26.2 (25C56); Darwin 25.2.0`、`Ubuntu 24.04 LTS; Linux 6.8.0` |
| `cpu` | object | 不可 | `{brand, logical_cores, physical_cores, perf_cores, efficiency_cores}`。各値は null 可。P/E コア数は macOS の `hw.perflevel0/1.physicalcpu`、Linux では null |
| `ram_bytes` | int | 可 | 物理メモリ量 |
| `arch` | string | 不可 | `platform.machine()`(`arm64`、`x86_64` など) |
| `compiler` | string | 可 | `<CMAKE_CXX_COMPILER> --version` の 1 行目 |
| `build_type` | string | 不可 | `CMAKE_BUILD_TYPE`。CMakeCache.txt が読めない環境のレコードは無効になる |
| `build_flags` | string | 可 | `compile_commands.json` の実効フラグ(パス・入出力を除く)。無ければ CMakeCache の `CMAKE_CXX_FLAGS[_<TYPE>]` |
| `diagnostics_enabled` | bool | 不可 | `SSSP_DIAGNOSTICS`(`release-diag` で true) |
| `power` | object | 可 | `{source: "AC"/"Battery"/"UPS"/null, battery_percent, low_power_mode, cpu_speed_limit}`。`cpu_speed_limit` は `pmset -g therm` の `CPU_Speed_Limit`(記録が無ければ null) |
| `load_average` | [number×3] | 可 | 計測開始時の 1/5/15 分 load average |

### 2.4 グラフ

| キー | 型 | null | 内容 |
|---|---|---|---|
| `graph_type` | string | 不可 | `sssp_gen --type` の名前(`er`、`grid`、`ba`、…) |
| `n` | int ≥ 1 | 不可 | 実際の頂点数 |
| `m` | int ≥ 0 | 失敗時のみ可 | 実際の無向辺数(生成後、重複除去後) |
| `avg_degree` | number | 可 | 指定した名目の平均次数 d(= 2m/n の目標)。次数が族で決まるもの(`grid`、`road`、`tree`、`adversarial_long_path`)は null |
| `weight_dist` | string | 不可 | `sssp_gen --weights` の名前。重みが構成で決まる族は `builtin` |
| `weight_type` | enum | 不可 | `uint64` / `double`(`sssp_gen --weight-type` と同じ) |
| `seed` | int ≥ 0 | 不可 | グラフ生成と始点選択の seed |
| `source` | int | 失敗時のみ可 | 始点。`0 ≤ source < n` |
| `reachable` | int | 可 | 始点から到達可能な頂点数(`≤ n`) |

### 2.5 試行

| キー | 型 | null | 内容 |
|---|---|---|---|
| `trial` | int ≥ 0 | 失敗時のみ可 | 試行番号(ウォームアップも 0 から数える) |
| `warmup` | bool | 不可 | ウォームアップの実行か。true は統計から除外 |
| `order_seed` | int | 可 | 実行順シャッフルの seed |
| `isolated` | bool | 不可 | このプロセスで 1 アルゴリズムだけを実行したか(`--isolate`) |
| `status` | enum | 不可 | `ok` / `timeout` / `memory_limit` / `error` |
| `error` | string | 可 | 失敗の説明 |

### 2.6 計測値

| キー | 型 | null | 内容 |
|---|---|---|---|
| `times_ns` | object | 不可 | 下の 6 キーをちょうど持つ。各値は 1 以上の整数 [ns] か null |
| `peak_rss_bytes` | int | 可 | 終了直前の `ru_maxrss`(bytes に正規化)。`isolated=false` では必ず null |
| `rss_after_graph_bytes` | int | 可 | グラフ構築直後の `ru_maxrss`(bytes) |
| `ops` | object | 可 | `include/sssp/result.hpp` の `OpCounters` の 8 キー(`vertex_scans`, `edge_scans`, `relaxations`, `comparisons`, `additions`, `queue_pushes`, `queue_pops`, `decrease_keys`)の部分集合。値は 0 以上の整数 |
| `counters` | object | 可 | `Diagnostics::counters`(名前 → 0 以上の整数) |
| `fallback` | object | 可 | `{count, records: [{reason, vertices, edges}], vertices, edges}`。`count = len(records)`、`vertices`/`edges` は records の和。記録していない場合は null |
| `correctness` | object | 不可 | `{verified, ok, method}` |

`times_ns` のキー(区間の定義は benchmark_protocol.md §4):

| キー | 内容 |
|---|---|
| `graph_gen` | グラフの生成または読み込み |
| `preprocess` | アルゴリズム固有の前処理(無ければ null) |
| `sssp` | SSSP 本体 |
| `tree` | 最短経路木の構築(本体と分離できなければ null) |
| `end_to_end` | 前処理 + 本体 + 木(グラフ生成と検証は含まない) |
| `verify` | 正解検証(計測区間の外) |

`correctness`:

- `verified=false`(検証していない): `ok` と `method` は null。
- `verified=true`: `ok` は bool、`method` は文字列(例 `reference:dijkstra-binary`、`bellman-ford`、`tree-check`)。

## 3. 整合性の規則

| 規則 | 理由 |
|---|---|
| `status=ok` なら `m`、`source`、`trial`、`times_ns.sssp`、`times_ns.end_to_end` は null 不可 | 成功した計測は必ずこれらを持つ |
| `times_ns.end_to_end ≥ preprocess, sssp, tree`(両方が値を持つとき) | 内訳が合計を超えない |
| `correctness.ok=false` なら `status=error` | 誤った結果の時間を統計に入れない |
| `diagnostics_enabled=false` なら `ops` と `counters` は null | release ではカウンタがコンパイルで消える。0 を書くと「数えて 0」と区別できない |
| `isolated=false` なら `peak_rss_bytes` は null | 複数アルゴリズムを同じプロセスで動かすと Peak RSS がどれのものか分からない |
| 既知のアルゴリズムの `mode` は固定。`baseline` は基準アルゴリズム専用 | 比較の分類を誤らない |
| `dijkstra-radix` は `weight_type=uint64` のみ | radix heap は単調な整数キーが前提 |
| `0 ≤ source < n`、`reachable ≤ n` | 範囲 |

## 4. 失敗レコード

ランナーが時間切れ・メモリ不足・異常終了で何も書けなかったとき、`scripts/run_benchmarks.sh` が
`python -m bench.schema failure ...` で 1 件追記する(`schema.make_failure_record`)。

- `status` は `timeout` / `memory_limit` / `error`、`error` に理由(例 `timeout after 600 s (exit 142)`)。
- 分からない値は null: `m`、`source`、`reachable`、`trial`、`times_ns.*`、`peak_rss_bytes`、`rss_after_graph_bytes`、
  `ops`、`counters`、`fallback`。`correctness` は `{verified: false, ok: null, method: null}`。
- 環境フィールドはそのジョブの `--env-json` から写す。

## 5. 例(形式の説明用。数値は計測値ではない)

```json
{"schema_version": "1", "run_id": "EXAMPLE", "timestamp": "2026-11-27T10:00:00+09:00",
 "paper_version": "2611.34911v1", "git_sha": "0000000", "git_dirty": false,
 "algorithm": "paper", "mode": "paper-faithful", "params": {"algo_seed": 1},
 "os": "EXAMPLE", "cpu": {"brand": "EXAMPLE", "logical_cores": null, "physical_cores": null,
 "perf_cores": null, "efficiency_cores": null}, "ram_bytes": null, "arch": "arm64",
 "compiler": null, "build_type": "Release", "build_flags": null, "diagnostics_enabled": false,
 "power": null, "load_average": null,
 "graph_type": "er", "n": 1000, "m": 2000, "avg_degree": 4, "weight_dist": "uniform_int",
 "weight_type": "uint64", "seed": 1, "source": 17, "reachable": 1000,
 "trial": 0, "warmup": false, "order_seed": 12345, "isolated": true, "status": "ok", "error": null,
 "times_ns": {"graph_gen": 1, "preprocess": null, "sssp": 1, "tree": null, "end_to_end": 1, "verify": 1},
 "peak_rss_bytes": null, "rss_after_graph_bytes": null,
 "ops": null, "counters": null, "fallback": null,
 "correctness": {"verified": true, "ok": true, "method": "reference:dijkstra-binary"}}
```

(実際のファイルでは 1 レコードを 1 行に書く。)

## 6. 変更の手順

スキーマを変えるときは `schema_version` を上げ、`scripts/bench/schema.py` と
`tests/python/test_bench_schema.py` を先に(テストから)直す。古い版のレコードは検証で弾かれる。
