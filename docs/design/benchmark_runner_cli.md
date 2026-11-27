# `sssp_bench_runner` の CLI 仕様

C++ のベンチランナー `build/<preset>/sssp_bench_runner` の引数・出力・終了コードを定める。
**ランナーは未実装**(C++ 担当が実装する)。`scripts/run_benchmarks.sh` はこの仕様どおりに呼び出す。
出力レコードの形式は [benchmark_schema.md](benchmark_schema.md)、計測手順は
[benchmark_protocol.md](benchmark_protocol.md)。

グラフの指定は、グラフ生成器 `sssp_gen`(docs/design/graph_library.md)と**同じ名前・同じ意味**にする。
ランナーは生成器ライブラリ(`sssp::gen`)をリンクし、同じオプション解析を使うことを想定する。

## 1. 呼び出し例

`run_benchmarks.sh` が 1 ジョブ(1 ビルド × 1 グラフ × 1 アルゴリズム)ごとに 1 回起動する:

```sh
build/release/sssp_bench_runner \
  --type er --n 100000 --avg-degree 4 --weights uniform_int --weight-type uint64 --seed 1 \
  --sources 2 --algorithms paper --trials 5 --warmup 1 --order-seed 12345 \
  --timeout-sec 600 --verify reference \
  --env-json results/raw/<run_id>.release.env.json --run-id <run_id> \
  --out results/raw/<run_id>.jsonl --isolate
```

## 2. 引数

### 2.1 グラフ(`sssp_gen` と同じ)

| 引数 | 内容 |
|---|---|
| `--type T` | `er`, `regular`, `ba`, `chung_lu`, `grid`, `road`, `tree`, `adversarial_long_path`(`adversarial_decrease_key` も受け付けてよいが、辺数が n² 級なのでプロファイルでは使わない) |
| `--n N` | 頂点数(`grid`/`road` は `rows = cols = ⌊√N⌋`) |
| `--avg-degree D` | `er`, `chung_lu` の平均次数 |
| `--degree D` | `regular` の次数 |
| `--attach M` | `ba` の 1 頂点あたりの接続数(平均次数は約 2M) |
| `--weights DIST` | `constant`, `uniform_int`, `uniform_real`, `log_uniform`, `exponential`, `zero_mixed`, `wide_range`。重みが構成で決まる族(`road`, `adversarial_*`)には `builtin` が渡る(重み指定なしとして扱う) |
| `--weight-type uint64\|double` | 重み型。レコードの `weight_type` にもこの名前を書く |
| `--lo`, `--hi`, `--mean`, `--zero-prob`, `--wide-log2`, `--exponent`, ... | `sssp_gen` と同じ(省略時は生成器の既定値)。指定した場合は `params` ではなく別途 run のメモに残す(現行スキーマにグラフパラメータ欄は無い) |
| `--seed S` | グラフ生成と始点選択の seed |
| `--graph-file PATH`(任意) | `sssp_gen` で作ったファイルを読む。このとき `graph_gen` は読み込み時間 |

レコードの `avg_degree` には `--avg-degree` / `--degree` の値、`ba` では `2 × --attach` を書く。
次数が族で決まるもの(`grid`, `road`, `tree`, `adversarial_long_path`)は null。

### 2.2 計測

| 引数 | 既定 | 内容 |
|---|---|---|
| `--sources K` | 1 | 始点の数。`--seed` から決定的に K 個選ぶ(同じグラフ・seed なら全アルゴリズムで同じ始点) |
| `--source-list a,b,...` | — | 始点を明示(`--sources` と排他) |
| `--algorithms A[,B,...]` | 必須 | 名前は schema の `ALGORITHM_MODES`。`mode` は名前から決まる |
| `--param KEY=VALUE` | — | 繰り返し可。基準でないアルゴリズムのパラメータ。レコードの `params` にそのまま入れる(数値は数値として) |
| `--trials T` | 5 | 本計測の試行数(始点ごと) |
| `--warmup W` | 1 | ウォームアップの回数。レコードは `warmup=true` で出力する |
| `--order-seed S` | 0 | 複数アルゴリズムのとき、試行ごとの実行順をこの seed でシャッフルする |
| `--timeout-sec X` | なし | 内部の目安。経過したら新しい試行を始めずに終える(強制終了は外側の `timeout`) |
| `--verify none\|reference\|bellman-ford` | `reference` | 計測区間の外で検証する。`reference` は `dijkstra-binary` との距離一致と木の整合性 |
| `--isolate` | off | 1 プロセス 1 アルゴリズム。`--algorithms` が 1 つでなければ終了コード 2。Peak RSS はこのときだけ記録 |

### 2.3 出力

| 引数 | 内容 |
|---|---|
| `--env-json PATH` | `scripts/bench/envinfo.py` の出力。その全キー(`schema.ENV_FIELDS`)を各レコードにそのまま写す |
| `--run-id ID` | レコードの `run_id` |
| `--out PATH` | JSONL を**追記**で書く。1 行 1 レコード、行ごとに flush |
| `--list` | 対応する `--type`、`--weights`、アルゴリズム(`name`, `mode`, 対応する重み型)を JSON で標準出力に出して終了 |
| `--version` | 版と論文版 `2611.34911v1` を出して終了 |

## 3. 振る舞い

1. グラフを 1 回だけ生成(または読み込み)し、`graph_gen_ns` を測る。直後に `ru_maxrss` を
   `rss_after_graph_bytes` として取る。全アルゴリズム・全試行で同じ `CsrGraph` を使う。
2. 始点ごとに、ウォームアップ W 回 → 本計測 T 回。各回 1 レコード。
3. 時間は `std::chrono::steady_clock`、単位 ns。`preprocess` / `sssp` / `tree` / `end_to_end` を分けて測る。
   分けられない区間は null(0 にしない)。
4. 検証・JSON の組み立て・書き出しは計測区間の外。検証に失敗したら `correctness.ok=false`、`status=error`。
5. カウンタは `SSSP_DIAGNOSTICS=1` のビルド(`release-diag`)でだけ `ops` / `counters` / `fallback` に書く。
   `SSSP_DIAGNOSTICS=0` のビルドでは `ops` と `counters` は必ず null。`fallback` は記録できるなら書く。
6. `--isolate` のとき、レコードはプロセス終了直前にまとめて書き、全レコードの `peak_rss_bytes` に
   終了直前の `getrusage(RUSAGE_SELF).ru_maxrss` を入れる。**macOS は bytes、Linux は KiB なので ×1024** して
   bytes に揃える。`--isolate` でなければ `peak_rss_bytes` は null、`isolated=false`。
7. `dijkstra-radix` に `--weight-type double` が来たら、何も書かずに終了コード 3。
8. 書いたレコードは `python -m bench.schema validate` を通ること(CI の smoke で確認する)。

## 4. 終了コード

| コード | 意味 | `run_benchmarks.sh` の扱い |
|---|---|---|
| 0 | 正常(検証失敗のレコードを含む場合も 0。失敗はレコードの `status=error` で表す) | 成功 |
| 1 | その他の異常 | `status=error` の失敗レコードを追記し、最後に終了コード 1 |
| 2 | 引数の誤り | 同上 |
| 3 | 未対応の組み合わせ(例 radix × double) | 同上(プランはこの組を作らない) |
| 5 | `std::bad_alloc` などメモリ不足 | `status=memory_limit` の失敗レコード |
| 124 / 142 | 外側の `timeout` / perl `alarm` による強制終了 | `status=timeout` の失敗レコード |
| 137 | SIGKILL(Linux の OOM killer など) | `status=error`(理由に記載) |

## 5. `run_benchmarks.sh` との対応

- プロファイル → ジョブ一覧は `scripts/bench/plan.py`(`python -m bench.plan <profile> --order-seed S` が TSV を出す)。
  列: `build, graph_type, n, avg_degree, shape_args, weight_dist, weight_type, algorithm, trials, warmup, sources, seed`。
  `shape_args` は `--avg-degree 4` / `--degree 4` / `--attach 2` / `-`(なし)。
- ジョブの順序は order seed でシャッフル済み。1 ジョブ = 1 プロセス(`--isolate`)。
- ランナーが無いとき: `run_benchmarks.sh` は「sssp_bench_runner が見つかりません」と出して終了コード 2。
  `smoke --allow-missing-runner`(CI 用)のときだけ「未実装」と出して終了コード 0。
