# ベンチマーク計測プロトコル

本書は、arXiv 2611.34911v1(無向・非負重み SSSP、期待 O(m + n log log log n))の実装と
Dijkstra 系ベースラインを比較するときの計測手順を定める。ここに書くのは**手順と設定値**だけで、
実測結果の数値は一切含まない。結果は `results/` と `docs/reports/` に実測から生成されたものだけを置く。

関連文書:

- レコードのスキーマ: [benchmark_schema.md](benchmark_schema.md)
- ランナーの CLI: [benchmark_runner_cli.md](benchmark_runner_cli.md)
- グラフ生成器(グラフ種類・重み分布の名前と定義): [graph_library.md](graph_library.md)

実装: `scripts/run_benchmarks.sh`(実行)、`scripts/bench/plan.py`(ジョブ一覧)、
`scripts/bench/envinfo.py`(環境)、`scripts/bench/schema.py`(検証)、`scripts/bench/stats.py`(統計)、
`scripts/analyze_results.py`(集計と図、`make report`)。

---

## 1. 目的と範囲

### 1.1 目的

- 論文アルゴリズムが、同じグラフ・同じ重み・同じメモリ配置(`include/sssp/graph.hpp` の共通 CSR)の上で、
  Dijkstra 系と比べてどの規模・どのグラフ種類で速い/遅いかを、再現可能な形で測る。
- 時間・メモリ・操作カウンタを分けて記録し、主張ごとに根拠となるビルドと指標を固定する。

### 1.2 比較対象(`algorithm`)

| `algorithm` | 内容 | `mode` | 対応する重み型 |
|---|---|---|---|
| `dijkstra-binary` | 二分ヒープ Dijkstra(Speedup の基準) | `baseline` | `uint64`, `double` |
| `dijkstra-4ary` | 4-ary ヒープ Dijkstra | `baseline` | `uint64`, `double` |
| `dijkstra-fibonacci` | Fibonacci ヒープ Dijkstra | `baseline` | `uint64`, `double` |
| `dijkstra-radix` | radix ヒープ Dijkstra | `baseline` | **`uint64` のみ** |
| `paper-faithful` | Mode A の入口(Alg 8)。外側の閾値 N > L_N^64 は実用規模で成り立たないため、基底ケース(l.5)を走査型 Dijkstra(Lemma 3.2)で解く。fallback として記録される | `paper-faithful` | `uint64`, `double` |
| `paper-faithful-d1` | 同上。基底ケースを二分ヒープ Dijkstra で解く(逸脱 D-1、ADR 0008) | `paper-faithful` | `uint64`, `double` |
| `paper` | 論文アルゴリズムの外側の構成(Alg 8 l.8-15)。予約名で未実装 | `paper-faithful` | `uint64`, `double` |
| `paper-experimental` | 論文の枠組みで部品・パラメータを試験的に変えたもの | `paper-experimental` | `uint64`, `double` |
| `paper-optimized` | 実用速度を狙った最適化版 | `practical-optimized` | `uint64`, `double` |

`mode` の意味:

- `baseline`: Dijkstra 系。
- `paper-faithful`: 論文の記述どおり(定数・パラメータも論文に合わせる)。漸近計算量の主張を確かめる用途。
- `paper-experimental`: 論文の枠組みを保ったまま、パラメータや部品を試験的に変えたもの。
- `practical-optimized`: 実用速度を狙って最適化したもの。論文との差分は必ず記録に残す。

`dijkstra-radix` × `double` のジョブはプラン(`scripts/bench/plan.py`)が作らない。対象外であり、失敗レコードにもしない。

`paper-faithful` / `paper-faithful-d1` について:

- `paper_sssp` は 1 回の呼び出しで、ゼロ辺の縮約(l.2)・基底ケース(l.5)・木の復元(l.6)を外から分けられない。
  全体を `sssp` として測り、`preprocess` と `tree` は null(§4)。
- 毎回 1 件の `FallbackRecord`(縮約後の頂点数 N と辺数)がレコードの `fallback` にそのまま入る
  (診断 OFF のビルドでも)。時間を報告するときは、この fallback の割合を必ず並べる。
  これらの時間は「論文の前処理 + Dijkstra」の時間であり、論文の外側の構成の時間ではない。
- `double` の検証は `Tolerance::for_path_hops(2(n − 1))`(numerics.md §3.2)。整数重みは厳密一致。
- `paper-faithful` の走査型基底ケースは Θ(N²) なので、ランナーは n > 2^17 を終了コード 3 で断る
  (別のアルゴリズムに置き換えない)。プランはその組を作らず、除外した件数と理由を標準エラーに出す
  (`plan.VERTEX_LIMITS`)。`paper-faithful-d1` に上限は無い。

### 1.3 重み型(`weight_type`)

名前は `sssp_gen --weight-type` と同じ。

- `uint64`: 整数重み。距離は厳密に計算される(オーバーフローは検出して例外)。
- `double`: 浮動小数重み。IEEE-754 の丸めを受ける。比較は許容誤差つき(§5)。

---

## 2. ビルド

`CMakePresets.json` の 2 つのプリセットを使い分ける。

| プリセット | `SSSP_DIAGNOSTICS` | 用途 |
|---|---|---|
| `release` | `OFF` | 時間計測用。公表する時間はこのビルドからのみ取る |
| `release-diag` | `ON` | 操作カウンタ(ヒープ操作数、緩和回数など)と fallback の記録用 |

```sh
cmake --preset release      && cmake --build --preset release
cmake --preset release-diag && cmake --build --preset release-diag
```

規則:

- **時間の主張は `release` のみ**から行う。`release-diag` はカウンタの更新が入るため、その時間は公表値に使わない。
- **カウンタの主張は `release-diag` のみ**から行う。`release` ではカウンタは記録されない(null)。
- 各レコードにビルドプリセット名、コンパイラとそのバージョン、`CMAKE_BUILD_TYPE` を記録する。

---

## 3. 試行設計

| 項目 | 既定 | 備考 |
|---|---|---|
| ウォームアップ | standard 1 回、full 2 回(smoke 0 回) | 記録は `warmup=true` として残すが、統計からは除外する |
| 反復回数(開発時) | 5 回(standard) | 初期確認・開発中の比較用 |
| 反復回数(公開用) | 20 回以上(full) | 公表する数値はこの条件を満たすものだけ |
| `release-diag` の反復 | smoke 1、standard 3、full 5 回 | カウンタ用。時間には使わない |
| 実行順 | ジョブ(1 プロセス)の順を order seed でシャッフル | `BENCH_ORDER_SEED`(既定 12345)。1 プロセスで複数アルゴリズムを動かす場合は、ランナーが試行ごとに `--order-seed` でシャッフル |
| 始点 | smoke 1、standard 2、full 4 本(`--sources`) | `--seed` から決定的に選ぶ |

- 実行順をシャッフルするのは、熱・キャッシュ・周波数の時間的な偏りが特定のアルゴリズムに
  集中しないようにするため。Peak RSS のために 1 プロセス 1 アルゴリズムで動かすので(§6)、シャッフルの単位は
  ジョブ(ビルド × グラフ × アルゴリズム)になる。
- 始点は同じ (graph, seed) なら全アルゴリズムで同じ k 本になる。1 試行 = 1 始点 × 1 アルゴリズムの 1 回実行。

---

## 4. 時間の分離

時刻は `std::chrono::steady_clock` で取り、単位はすべてナノ秒 [ns]。

| JSONL のキー(`times_ns.*`) | 名前 | 範囲 |
|---|---|---|
| `graph_gen` | `graph_gen_ns` | グラフの生成または読み込み |
| `preprocess` | `preprocess_ns` | アルゴリズム固有の前処理。該当なしは null |
| `sssp` | `sssp_ns` | SSSP 本体(距離の確定まで) |
| `tree` | `tree_ns` | 最短経路木(parent 配列)の構築。本体と分離できない実装では null とし、その旨を注記する |
| `end_to_end` | `end_to_end_ns` | 前処理 + 本体 + 木。**グラフ生成は含まない** |
| `verify` | `verify_ns` | 正解検証。必ず計測区間の外で行い、`end_to_end` に含めない |

規則:

- 未計測の値は **null** とする。0 を入れてはいけない(0 は「計測して 0 だった」と区別できなくなる)。
- `tree` が null の実装では、木の構築時間は `sssp` に含まれている。比較の際はこの差を注記する。
- 検証・ログ出力・JSON 書き出しは計測区間に入れない。

---

## 5. 正解検証

検証は計測区間の外で行う(`verify_ns` として別に記録)。

1. **距離の一致**: 参照実装の距離と全頂点で比べる。
   - 参照実装は `dijkstra-binary`。小規模グラフでは `bellman-ford` も使い、参照実装自体も確かめる。
   - `uint64`: 完全一致。
   - `double`: 許容誤差つき比較(相対誤差と絶対誤差の併用)。許容値は `docs/design/numerics.md` に従う。
   - 到達不能な頂点は双方とも infinity であること。
2. **木の整合性**: 始点以外の到達可能な頂点 v について、parent 辺 (parent[v], v) が実在し、
   `dist[parent[v]] + w(parent[v], v) = dist[v]`(`double` は許容誤差つき)が成り立ち、parent をたどると始点に戻ること。
3. 検証に失敗したレコードは `status=error` とし、時間の統計から除外する。

---

## 6. Peak RSS(メモリ)

- 1 プロセスで 1 つの (algorithm, graph) だけを実行する(ランナーの `--isolate` 相当)。
  `scripts/run_benchmarks.sh` は algorithm ごとに別プロセスで起動する。前のアルゴリズムの確保分が混ざらないようにするため。
- プロセス終了直前に `getrusage(RUSAGE_SELF).ru_maxrss` を取得する。
- 単位は OS で異なる。**macOS は bytes、Linux は KiB**。Linux では ×1024 して bytes に正規化し、`peak_rss_bytes` に記録する。
- グラフ構築直後の値も `rss_after_graph_bytes` として記録し、`peak_rss_bytes - rss_after_graph_bytes` を
  アルゴリズム自身の増分の目安とする。
  - 注意: `ru_maxrss` は単調増加の最大値なので、この差分は厳密な増分ではなく目安にすぎない
    (グラフ生成の一時領域がピークを作った場合など、差分が過小になる)。

---

## 7. Apple Silicon(macOS)での注意

### 7.1 CPU 固定はできない

- Apple Silicon には P コア(性能)と E コア(効率)がある。数は次で確かめられる。
  ```sh
  sysctl hw.perflevel0.physicalcpu   # P コア数
  sysctl hw.perflevel1.physicalcpu   # E コア数
  ```
- macOS のスケジューラは計測プロセスを E コアに載せることがある。
- **macOS ではコアの固定(CPU pinning)を保証する手段が無い。** `taskpolicy` などは優先度や QoS の指示であり、
  特定コアへの固定ではない。したがって macOS での計測値は、この揺らぎを含むものとして扱う。

### 7.2 揺らぎの原因

- 熱スロットリング(長時間の連続実行で周波数が下がる)。
- 電源: AC かバッテリーか、低電力モードの有無。`pmset -g`、`pmset -g batt` で確かめる。
- バックグラウンド処理: Spotlight の索引作成、Time Machine、ブラウザ、同期クライアントなど。

### 7.3 対策

- AC 電源に接続する。低電力モードを off にする。
- 他のアプリを終了する。可能なら Spotlight の索引作成や Time Machine のバックアップが走っていない時間に回す。
- 試行間にクールダウン(待ち時間)を入れる。
- 記録できる状態は `scripts/bench/envinfo.py` が各 run に記録する(各レコードと
  `results/raw/<run_id>.<preset>.env.json`。終了時の状態は `<run_id>.<preset>.env-end.json`)。
  - 電源状態(AC/バッテリー、低電力モード)、`pmset -g therm` の `CPU_Speed_Limit`、CPU 名、
    OS とバージョン、P/E コア数、load average(`os.getloadavg()`、`uptime` と同じ値)。
- 開始時または終了時の `cpu_speed_limit` が 100 未満の run は、熱による制限を受けたものとして無効とする。
- 開始時または終了時の環境情報(`env.json` / `env-end.json`)を取得できなかった run は、上の判定ができないため
  無効とする。`scripts/run_benchmarks.sh` は取得の失敗を標準エラーに出し、終了時の失敗数を最後の
  `done:` 行に `env_end_failed=<数>` として記録し、終了コード 1 で終える(開始時の失敗はその場で終了コード 1)。
  計測済みのレコードは消さずに残すが、その run の値は公表に使わない。
  - 取得には成功したが `cpu_speed_limit` が `null`(Linux など、取得手段の無い環境)の場合は失敗ではない。
    その run は熱の判定ができないものとして報告書に書く。
- 記録された load average が高い run は、公表値に使うかどうかを個別に判断し、使わない場合は理由を残す。

### 7.4 Linux の場合

- `taskset -c <core>` でコアを固定できる。
- `cpupower frequency-set -g performance` で周波数ガバナを固定できる(権限がある場合)。
- `envinfo.py` はこれらの設定を自動では取得しない。使った場合は run のメモ(報告書)に書く。

---

## 8. タイムアウトとメモリ上限

### 8.1 タイムアウト

`scripts/run_benchmarks.sh` が 1 呼び出しごとに外部タイムアウトをかける。

| プロファイル | 既定のタイムアウト(1 ジョブあたり) |
|---|---|
| `smoke` | 60 秒 |
| `standard` | 600 秒 |
| `full` | 3600 秒(`make benchmark-full` は Makefile の `BENCH_TIMEOUT`、既定 3600 を渡す) |

環境変数 `BENCH_TIMEOUT` を設定すると、どのプロファイルでもその値を使う。

macOS には GNU `timeout` が標準で入っていないため、次の順にフォールバックする。

1. `timeout`(GNU coreutils)
2. `gtimeout`(Homebrew の coreutils)
3. `perl -e 'alarm ...; exec ...'`(macOS の既定。時間切れの終了コードは 142、GNU `timeout` は 124)

### 8.2 メモリ上限

- `ulimit -v`(環境変数 `BENCH_MEM_LIMIT_MB`)で上限をかける。
- macOS では `RLIMIT_AS` を設定できない(`ulimit -v` が `Invalid argument` で失敗する)。その場合は警告を出して
  上限なしで続行し、終了後に `peak_rss_bytes` が上限を超えたレコードを警告として表示する(status は変えない)。
- ランナーが `std::bad_alloc` で終了コード 5 を返した場合は `status=memory_limit` とする。

### 8.3 失敗の扱い

- タイムアウト・メモリ超過・異常終了(検証失敗を含む)は、それぞれ `status=timeout` / `memory_limit` / `error`
  のレコードとして残す。
- ランナーが何も書けなかった場合は、`run_benchmarks.sh` が失敗レコードを 1 件追記する(`python -m bench.schema failure`)。
- これらは時間の統計から除外し、**失敗率**としてジョブ単位で別に報告する(`results/summary/status.csv`)。
- 時間切れ・メモリ超過は想定内なので `run_benchmarks.sh` の終了コードは 0。異常終了(`error`)が 1 件でもあれば 1。
- 失敗した組の値を推定・補間しない。

---

## 9. 指標と統計

### 9.1 Speedup

```
Speedup = T(dijkstra-binary) / T(対象)
```

- 同じ (`graph_type`, `n`, `avg_degree`, `weight_dist`, `weight_type`) の組の中で、各アルゴリズムの **median** を使う。
- Speedup > 1 なら対象が `dijkstra-binary` より速い。
- 既定の T は `sssp_ns`。`end_to_end_ns` 版も併せて出す。
- 成功した計測(`status=ok`、warmup でない、`release` ビルド)が基準・対象のどちらかに無い組は Speedup を出さない。
- `params` が異なるものは別の系列として、それぞれ Speedup を出す。

### 9.2 統計量

| 統計量 | 扱い |
|---|---|
| median | 主指標 |
| mean | 参考 |
| p95 | 標本数 < 20 のときは参考値としてフラグを立てる |
| stddev | 標本標準偏差(ddof=1)。標本 1 個なら null |
| min / max | 参考 |

- p95 は順序統計量の線形補間(numpy の既定)。

- ウォームアップ(`warmup=true`)と失敗レコードは統計に入れない。

### 9.3 Speedup の信頼区間

- ブートストラップで求める。既定はリサンプル数 10000、乱数シード 20261127(固定)、95% パーセンタイル区間
  (`scripts/analyze_results.py --bootstrap-resamples / --bootstrap-seed`)。
- 基準側と対象側をそれぞれ独立にリサンプルし、各回で median の比を計算する。
- bootstrap seed は集計側の設定であり、レコードには入らない。既定値を変えて集計した場合は報告に書く。

---

## 10. グラフ規模とプロファイル

### 10.1 規模の格子

- n ∈ {1e3, 1e4, 1e5, 1e6, 1e7}
- 平均次数 d ∈ {2, 4, 8, 16, 32}(無向グラフで d = 2m/n)

### 10.2 プロファイル

| プロファイル | n | d | グラフ種類 | 重み分布 | trials | warmup | diag trials | 始点 |
|---|---|---|---|---|---|---|---|---|
| `smoke` | {1e3, 1e4} | 4 | `er`, `grid` | `uniform_int` | 1 | 0 | 1 | 1 |
| `standard` | {1e3, 1e4, 1e5, 1e6} | {4, 16} | 全種類 | `constant`, `uniform_int`, `wide_range` | 5 | 1 | 3 | 2 |
| `full` | {1e3, …, 1e7} | {2, 4, 8, 16, 32} | 全種類 | 全分布 | 20 | 2 | 5 | 4 |

- 定義は `scripts/bench/plan.py` の `PROFILES`。`python -m bench.plan <profile> --order-seed S` でジョブ一覧を確認できる。
- 絞り込み: `BENCH_ALGORITHMS`、`BENCH_GRAPH_TYPES`、`BENCH_SIZES`、`BENCH_BUILDS`(カンマ区切り)。

- `smoke` は動作確認用であり、その結果を性能の主張に使わない。
- 規模によっては時間・メモリの上限に達しうる。達した組は失敗レコードとして記録し、補間しない。
- 次数が族で決まるもの(`grid`, `road`, `tree`, `adversarial_long_path`)は d ごとに複製しない(`avg_degree` は null)。
- 重みが構成で決まるもの(`road`, `adversarial_long_path`)は重み分布を変えない(`weight_dist=builtin`)。

---

## 11. グラフ種類(`graph_type`)

名前と定義は `sssp_gen --type`(docs/design/graph_library.md §6)に従う。

| `graph_type` | 内容 | d の渡し方 | 注記 |
|---|---|---|---|
| `er` | 一様ランダム G(n, m)、m = round(n·d/2) | `--avg-degree d` | |
| `regular` | ランダム d-正則 | `--degree d` | |
| `ba` | Barabási–Albert(優先的選択) | `--attach d/2` | 平均次数は約 2·attach |
| `chung_lu` | Chung–Lu(べき乗則の期待次数) | `--avg-degree d` | 次数の偏りが大きい |
| `grid` | 2 次元格子 ⌊√n⌋ × ⌊√n⌋ | なし(次数約 4) | n は平方数に丸められる |
| `road` | 道路網風(格子の摂動 + 近距離ショートカット) | なし | 重みはユークリッド長(`builtin`) |
| `tree` | 一様ランダム木 | なし | m = n − 1 |
| `adversarial_long_path` | 長い経路 + 迂回辺 | なし | 深い最短経路木。重みは構成で固定(`builtin`) |

- `adversarial_decrease_key` は辺数が n² 級なので規模の格子には入れない。使う場合は小さな n で個別に回す。
- 実際に生成された n と m を記録する(丸めや重複辺の除去で指定値とずれうるため)。
- 連結性: 始点から到達できない頂点があってもそのまま扱い、到達可能な頂点数を `reachable` として記録する。

---

## 12. 重み分布(`weight_dist`)

名前と定義は `sssp_gen --weights`(docs/design/graph_library.md §5)に従う。パラメータ(`lo`, `hi`, `mean`,
`zero_probability`, `wide_log2`)はプロファイルでは生成器の既定値を使う。

| `weight_dist` | 内容 | プロファイルでの重み型 |
|---|---|---|
| `constant` | 全辺同じ重み(既定 `lo`) | `uint64` |
| `uniform_int` | [lo, hi] の一様整数 | `uint64` |
| `log_uniform` | 対数一様 | `uint64` |
| `wide_range` | 2 進指数を一様に選ぶ(ダイナミックレンジが大きい) | `uint64` |
| `zero_mixed` | 確率 `zero_probability` で 0、それ以外は一様 | `uint64` |
| `uniform_real` | [lo, hi) の一様実数 | `double` |
| `exponential` | 指数分布 | `double` |
| `builtin` | 生成器の構成で決まる重み(`road`, `adversarial_long_path`) | `uint64` |

- `uint64` で距離の和がオーバーフローしうる組は、検出された場合 `status=error` とする。

---

## 13. 再現性

全レコードに次を記録する。

- シード: graph seed と始点選択の seed(`seed`)、`order_seed`、アルゴリズムの乱数 seed(`params` に入れる)。
  bootstrap seed は集計側の既定値(§9.3)
- git の SHA と dirty フラグ(作業ツリーに未コミットの変更があるか)
- 再現対象の論文版: `2611.34911v1`
- ビルドプリセット、コンパイラ、環境情報(§7.3)

規則:

- **`dirty=true` のレコードは公表値に使わない。** 開発中の確認にだけ使う。

---

## 14. 出力と集計

```
results/raw/<run_id>.jsonl          (gitignore 済み、1 行 1 レコード)
        │
        └─ make report  (scripts/analyze_results.py)
              ├─ results/summary/*.csv
              └─ docs/reports/figures/*.png
```

- `results/summary/` の CSV:
  - `timing_summary.csv`: 構成 × アルゴリズム × params × 区間ごとの samples / median / mean / p95 / stddev / min / max /
    p95_reference_only / dirty_records / run_ids
  - `speedup.csv`: Speedup と 95% CI(`sssp` と `end_to_end`)
  - `memory.csv`: Peak RSS の median / max、グラフ構築直後の RSS の median
  - `status.csv`: ジョブ単位の ok / timeout / memory_limit / error と失敗率
  - `fallback.csv`: `release-diag` の fallback 発生率、回数、頂点・辺の割合、理由
  - `ops.csv`: 操作カウンタと `counters` の median
- `docs/reports/figures/` の図: `runtime_vs_n__<type>`、`runtime_vs_m__<type>`、`speedup_vs_n__<type>`、
  `memory_vs_n__<type>`、`fallback_rate__<type>`、`breakdown__<type>`(前処理/本体/木の内訳)、
  `speedup_by_graph_type`、`speedup_by_weight_dist`、`param__<algorithm>__<param>`(閾値パラメータと時間)。
- レコードのフィールド定義は [benchmark_schema.md](benchmark_schema.md)、ランナーの引数は
  [benchmark_runner_cli.md](benchmark_runner_cli.md) を参照する。
- 図表では、未計測の組・失敗した組を描画しない。系列の途中に欠けた点があれば線をそこで切る(補間しない)。
- スキーマ違反のレコードがあれば `make report` は失敗する(`--skip-invalid` で除外して続行)。
- 集計対象のデータが無い場合は「データなし」と出力し、値を作らない。

---

## 15. 結果の解釈上の注意

- 漸近計算量の差(O(m + n log n) 対 O(m + n log log log n))は、実用的な n では定数倍とキャッシュ効果に隠れうる。
  速くならないことは、実装の誤りとも論文の誤りとも限らない。
- op カウンタ(`release-diag`)は比較加算モデルでの計算量の代理指標であり、実時間と比例するとは限らない
  (メモリアクセスのパターンや分岐予測は数えていない)。
- `paper` の期待計算量は乱択に対する期待値である。試行間のばらつきは乱択とシステムの揺らぎの両方を含む。
- グラフ種類・重み分布によって有利不利が変わりうる。特定の組の結果を一般化しない。
- macOS での結果は CPU 固定なしの条件で取ったものであることを、報告に必ず明記する。
