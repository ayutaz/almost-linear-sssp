# グラフライブラリ

[English](graph_library.md) | **日本語** | [简体中文](graph_library.zh-CN.md)

すべての SSSP 実装、ベンチマーク、テストが共有する、グラフの表現、テキスト IO、再現可能なジェネレータ。グラフは**無向**で重みは**非負**とする。重みの型は `std::uint64_t`（厳密）と `double`（IEEE-754）である。

## 1. 構成

| コンポーネント | ファイル | 名前空間 | 役割 |
|---|---|---|---|
| CSR グラフ | `include/sssp/graph.hpp` | `sssp` | `CsrGraph<W>`、`Edge`、`Arc`、`validate_edges`、`WeightTraits`。すべてのアルゴリズムが共有する。このライブラリは**変更しない**。 |
| 形式 | `include/sssp/graph_formats.hpp`（ヘッダのみ） | `sssp` | `EdgeList<W>`、`AdjacencyList<W>`、変換、`timed()` |
| テキスト IO | `include/sssp/graph_io.hpp`、`src/graph/io.cpp` | `sssp::io` | プロジェクト独自の辺リスト形式と DIMACS `.gr`、`ParseError` |
| ジェネレータ | `include/sssp/generators.hpp`、`src/graph/generators.cpp` | `sssp::gen` | PRNG、決定的な数学関数、重みの分布、グラフジェネレータ、カタログ、`generate()` |
| CLI | `tools/sssp_gen.cpp` | — | `sssp_gen`: グラフを生成してファイルに書き出す |

`io.cpp` と `generators.cpp` のテンプレートは、`std::uint64_t` と `double` について明示的にインスタンス化している。

## 2. 表現と変換

- `EdgeList<W> { std::size_t num_vertices; std::vector<Edge<W>> edges; }`: IO とジェネレータがやり取りする形式。
- `AdjacencyList<W>::from_edges(n, edges)`: vector の vector である。`CsrGraph<W>` と同じ `Arc<W>` のレイアウトと同じアクセサ（`num_vertices`、`num_edges`、`num_arcs`、`degree`、`arcs(v)`、`edge(id)`、`edges()`）を持つ。入力の検証も `CsrGraph::from_edges` と同じ方法で行う。
- 変換はキャッシュを持たない自由関数である。

  | 関数 | 入力 |
  |---|---|
  | `to_csr` | `const EdgeList&`、`EdgeList&&`（辺の vector をムーブする）、または `const AdjacencyList&` |
  | `to_edge_list` | `CsrGraph` または `AdjacencyList` |
  | `to_adjacency_list` | `EdgeList` または `CsrGraph` |

**順序の契約。** どの表現も同じ辺 ID を使う。`EdgeId i` はどこでも同じ辺を指す。`AdjacencyList::arcs(v)` は `CsrGraph::arcs(v)` とまったく同じ順序で弧を並べる。辺を ID の順に見て、u 側の弧を先に、次に v 側の弧を置く。自己ループ `{u,u}` は弧1本である。このため、アルゴリズムはどちらのレイアウトでも弧を同じ順序でたどる。

**変換の時間は別に測る。** ベンチマークは、変換のコストをアルゴリズムの時間に含めてはならない。

```cpp
auto [g, convert_seconds] = sssp::timed([&] { return sssp::to_csr(std::move(edge_list)); });
// run and time the algorithm on g; report convert_seconds as a separate column
```

`timed(f)` は `Timed<T>{value, seconds}` を返し、`std::chrono::steady_clock` を使う。

## 3. ファイル形式（`sssp::io`）

### プロジェクト独自の辺リスト形式（`*.txt`）

```
n m
u v w      # exactly m lines, 0-based vertices, undirected edge {u,v}
```

- 空行と、最初の空白でない文字が `#` の行は、ファイルのどこにあっても無視する。
- ファイルには辺の行がちょうど `m` 行なければならない。少なくても多くてもエラーである。
- 自己ループと平行辺は、書かれたとおりに残す。

### DIMACS 9th Challenge `.gr`

```
c comment
p sp n m   # exactly once, before the first arc; m = number of ARC lines
a u v w    # 1-based, directed arc u->v
```

グラフは2段階で無向グラフとして読み込む。

1. `DimacsArcMode` が、弧をどう辺にするかを決める。
   - `kArcAsEdge`（既定）: 各弧を1本の辺にする。そのため、対称な弧の組は2本の平行辺になる。距離は変わらないが、`m` は2倍になる。
   - `kMergeReversePairs`: 弧 `u->v` と `v->u` のうち、重みがビットまで同一のものを、ファイルの順に貪欲に対応付ける。対応付けた組は**1本**の辺になり、組の最初の弧の位置と向きに置く。相手が見つからない弧（一方向の弧、または重みが非対称な弧）は、それぞれ1本の辺になる。自己ループの弧は常に1本の辺である。DIMACS の道路網にはこのモードを使う。
2. 次に `DuplicateEdgePolicy` が、残った平行辺を扱う。2本の辺は、端点の無順序対が同じなら平行である。同じ頂点上の2つの自己ループも平行とみなす。

   | ポリシー | 効果 |
   |---|---|
   | `kKeepAll`（既定） | すべての辺を残す（多重グラフ） |
   | `kKeepMinWeight` | 最小の重みを持つ辺を1本だけ、最初の位置と向きで残す |
   | `kReject` | 最初の重複がある行で `ParseError` を投げる |

`write_dimacs` は、すべての辺を2本の弧（`u->v`、続いて `v->u`）として書き、自己ループは1本の弧として書く。このファイルを `kMergeReversePairs` で読み戻すと、元の辺リストが完全に復元される。順序も向きも同じで、重みもビットまで同一である。

### API

- ストリーム: `read_edge_list<W>(istream&)`、`write_edge_list(ostream&, el)`、`read_dimacs<W>(istream&, DimacsReadOptions{arc_mode, duplicates})`、`write_dimacs(ostream&, el)`。
- ファイル: `read_edge_list_file`、`write_edge_list_file`、`read_dimacs_file`、`write_dimacs_file`。ファイルを開けない、または書き込めないときは `std::runtime_error` を投げる。

### 数値

- 解析と出力には `std::from_chars` / `std::to_chars` を使うので、結果はロケールに依存しない。
- `uint64_t` の重みは、0..2^64-1 の符号なし 10 進数でなければならない。符号、小数部、指数は拒否する。
- `double` の重みは、有限で非負の 10 進形式か指数形式であればよい。書き出しには `chars_format::general` と精度 17（`%.17g` と同じ）を使う。これは有限の double をビットまで正確にラウンドトリップする。`-0` は `+0` に正規化する。

### エラー

不正な入力はすべて `io::ParseError`（`std::runtime_error` の派生）を投げる。メッセージの形式は `line N: reason` で、`line()` は 1 始まりの行番号を返す（0 は入力の終端、または入力全体を表す）。検出するエラーは次のとおり。

- 負、NaN、無限大、またはオーバーフローする重み
- 範囲外の頂点（DIMACS での頂点 0 を含む）
- ヘッダの欠落、重複、または不正
- 辺または弧の行が少なすぎる、または多すぎる
- 数値でないトークン、または行のトークン数が違う
- `n >= kNoVertex` または `m >= kNoEdge`

リーダは事前確保を 2^24 要素に制限しているので、ヘッダの `m` が巨大でも、巨大なバッファを先に確保することはない。

## 4. 再現性の契約

すべてのジェネレータは (パラメータ, シード) の純粋関数であり、どのプラットフォームでもどの標準ライブラリでも**ビットまで同一**の出力を返す。次の規則がそれを成り立たせている。

- **独自の PRNG。** シードの生成とサブストリームのシードの導出には `SplitMix64` を、乱数ストリームには `Xoshiro256ss`（xoshiro256** 1.0）を使う。`<random>` は使わない。C++ 標準は `std::*_distribution` クラスのアルゴリズムを実装定義としているので、同じエンジン出力でも libstdc++、libc++、MSVC で違う値になるからだ。どちらの生成器も、参照実装の既知解ベクタで検査している。
- **厳密な整数。** `uniform_below(b)` は `2^64 mod b` より小さい生の値を棄却し、`r % b` を返す。モジュロの偏りがなく、浮動小数点も使わない。`uniform_range(lo, hi)` は両端を含み、64 ビット全域を扱える。
- **厳密な実数。** `uniform01()` は `(next() >> 11) * 2^-53` を計算する。53 ビットの乱数を持つ [0,1) の値で、丸めは起きない。
- **決定的な数学関数。** `det_log` と `det_exp` は、正しく丸められる IEEE 演算（`+ - * /`、`frexp`、`ldexp`、`floor`）だけを使う。libm の `log`、`exp`、`pow` は正しく丸められず、プラットフォーム間で異なるので使わない。`det_log` は [√½, √2) の仮数に対する atanh 級数を使う。`det_exp` は ln 2 の Cody–Waite 分割（fdlibm の `ln2_hi` / `ln2_lo`）で引数を縮約し、そのあとテイラー級数を足し合わせる。どちらも libm との相対誤差は約 1e-14 以内で一致する。`sqrt` は IEEE 754 により正しく丸められるので、安全に使える。
- **FMA の縮約なし。** 融合積和は2回ではなく1回だけ丸めるので、式を縮約すると結果が変わってしまう。`generators.cpp` には `#pragma clang fp contract(off)` がある。GCC は何もしなくてよい。ISO モード（`-std=c++20`、`CMAKE_CXX_EXTENSIONS OFF`）では、すでに `-ffp-contract=off` を使うからだ。
- **順序なし反復の禁止。** `unordered_set` と `unordered_map` は、所属判定と個数の計算にだけ使う。出力の順序が、それらの反復順に依存することはない。
- **独立したサブストリーム。** `substream_seed(seed, stream)` が、用途ごとに1本のストリームを導出する。トポロジ = 0、重み = 1、連結用の辺 = 2。したがって、同じシードならば、重みの分布が何であっても、`W` が `uint64_t` でも `double` でも、トポロジは同一になる。

### バージョン間の出力の変更

- **P2-GRAPH-GEOMETRIC（新しい型 `geometric`）。** ジェネレータを追加しても、既存のどの型の出力も変わらなかった。変更の前後で `sssp_gen` を 2,526 通りの引数で比較した。既存の9種類の型 × 両方の重み型 × 7つの重み分布すべて × シード {1, 2, 3, 17, 1000} × n ∈ {50, 1000} × `--connected` の有無（regular は keep/drop/switch を順に巡回）に加えて、`--int-scale` 付きの road と DIMACS 出力である。ハッシュはすべて同一で、`--list` には新しい項目が末尾に増えただけである。

- **P2-DIFF-REGULAR-SWITCH（`random_regular`、`kSwitchToSimple`、`sssp_gen --type regular --multi switch`）。** この修正の前は、有効な引数の一部で「switching did not converge」を投げていた。多くは小さく密なケース（d = n − 1、n ≤ 12 でシードの約 3〜10%。例: n = 3、d = 2、seed = 20）である。これらの引数は、いまでは単純な d 正則グラフを返す。ジェネレータは、同じトポロジのストリームから引いた新しい構成から、やり直す（§9 を参照）。以前に成功した引数は、すべてバイト単位で同一の出力になる。最初の構成とそのスワップが、これまでとまったく同じようにストリームを使うからである。これは、古いライブラリと新しいライブラリを 74,955 通りの引数で比較して確かめた。3つのポリシーすべて、n ≤ 16 × 有効なすべての d × シード 1..200、および n ∈ {100, 1000, 5000} × d ∈ {2, 3, 4, 7, 10, 20} の `connected` の有無である。以前に成功していた 74,778 通りは、ハッシュが同一だった。残る 177 通りは、以前は例外を投げていたが、いまは成功する。`kKeep` と `kDrop` は変わっていない。

## 5. 重みの分布（`WeightSpec`）

`WeightSpec` のフィールドは `dist`、`lo`（既定 1）、`hi`（100）、`mean`（10）、`zero_probability`（0.25）、`wide_log2`（50）である。これらが矛盾していると、`validate_weight_spec` が `std::invalid_argument` を投げる。

どの分布も、`lo` と `hi` が有限で、`lo >= 0` かつ `hi <= 2^53` であることを要求する。整数の `W` では、実数値の標本を最も近い整数に丸める。同じ近さのときは切り上げる。

| `WeightDist` | CLI 名 | 標本 | `uint64_t` | `double` | 追加の制約 |
|---|---|---|---|---|---|
| `kUniformInt` | `uniform_int` | [⌈lo⌉, ⌊hi⌋] の整数上の一様分布 | 厳密 | 整数値の double | 区間が整数を含む |
| `kUniformReal` | `uniform_real` | `lo + (hi-lo)·U`、U ∈ [0,1) | 丸める | 実数 | lo ≤ hi |
| `kLogUniform` | `log_uniform` | `exp(uniform(log lo, log hi))`、[lo,hi] に収める | 丸める | 実数 | 0 < lo ≤ hi |
| `kExponential` | `exponential` | `-mean·log(1-U)` | 丸める | 実数 | 0 < mean ≤ 2^47 |
| `kConstant` | `constant` | `lo`（重みをすべて 0 にするには `lo = 0` を使う） | 丸める | `lo` | — |
| `kZeroMixed` | `zero_mixed` | 確率 `zero_probability` で 0、それ以外は [lo,hi] の一様分布 | 整数の一様分布 | [lo,hi) の実数の一様分布 | zero_probability ∈ [0,1]。区間が整数を含む |
| `kWideRange` | `wide_range` | 2進指数 e を [0, K] から一様に選び、[2^e, 2^(e+1)) の中で一様に選ぶ。e = K ならちょうど 2^K。[1, 2^K] を覆う | 厳密 | `2^e·(1+U)` | K = `wide_log2` ∈ [0, 62] |

`to_string(WeightDist)` と `parse_weight_dist(name)` は、列挙型と CLI 名を相互に変換する。`parse_weight_dist` は未知の名前に対して例外を投げる。`sample_weight<W>(spec, rng)` は標本を1つ引く。

## 6. ジェネレータ

すべてのジェネレータは `EdgeList<W>` を返し、`GenOptions{weights, seed = 1, connected = false}` を受け取る。

| CLI `--type` | 関数 | パラメータ | 辺数 | 単純 | 連結 |
|---|---|---|---|---|---|
| `er` | `erdos_renyi(n, avg_degree, opt)` / `erdos_renyi_gnm(n, m, opt)` | m = round(n·avg_degree/2) ≤ n(n−1)/2 | m | はい（対の添字に対する Floyd の標本抽出） | `connected` のときだけ |
| `regular` | `random_regular(n, d, policy, opt)` | n·d は偶数。`kSwitchToSimple` は d < n が必要 | n·d/2（`kDrop` ではそれより少ない） | `kKeep`: いいえ。`kDrop`: はい（次数 ≤ d）。`kSwitchToSimple`: はい（次数はちょうど d。次数を保つ二重辺スワップによる） | `connected` のときだけ |
| `ba` | `barabasi_albert(n, m, opt)` | 1 ≤ m < n | m(m+1)/2 + (n−m−1)·m | はい | はい |
| `grid` | `grid_2d(rows, cols, opt)`。頂点 ID = r·cols + c | rows、cols（既定はそれぞれ ⌊√n⌋） | rows(cols−1) + cols(rows−1) | はい | はい |
| `road` | `road_like(rows, cols, RoadLikeOptions, opt)` | jitter（既定 0.3、< 0.5）、keep_probability（0.9）、shortcut_probability（0.15）、shortcut_radius（2、≥ 2）、integer_scale（1000） | ≈ grid × keep_probability + ショートカット + 復元した街路 | はい | はい |
| `tree` | `random_tree(n, opt)` | n | n−1 | はい | はい（ランダムな Prüfer 列による、ラベル付き木上の一様分布） |
| `chung_lu` | `chung_lu(n, avg_degree, exponent, opt)` | exponent > 2。m ≤ n(n−1)/4 | m = round(n·avg_degree/2) | はい（自己ループと重複は引き直す） | `connected` のときだけ |
| `adversarial_decrease_key` | `adversarial_decrease_key(k)`。CLI では k = n−1 | k | k + k(k−1)/2 | はい | はい |
| `adversarial_long_path` | `adversarial_long_path(n)` | n | (n−1) + (n−2) | はい | はい |
| `geometric` | `random_geometric(n, r, GeometricOptions, opt)` | r ≥ 0 で有限（既定は `geometric_radius(n, avg_degree)` = √(avg_degree/(πn))）。重みは `kFromSpec` または `kEuclidean`。integer_scale（1000） | ≈ n·(n−1)·F(r)/2。F は下記 | はい | `connected` のときだけ |

分かりにくいジェネレータの詳細は次のとおり。

- **道路風（Road-like）。** 格子点は各軸方向に U[−jitter, jitter] だけずらす。各街路（格子の辺）は確率 `keep_probability` で残る。確率 `shortcut_probability` で、頂点はショートカットを1本持つ。つなぎ先は、チェビシェフ半径 R 以内にあって格子の隣接でない頂点である。そのあと、取り除いた街路を Kruskal の順に復元してグラフを連結にするので、どの辺も局所的なままである。重みは、ずらした点のユークリッド長である。整数の重みは `max(1, round(len·integer_scale))` とする。road ジェネレータは `opt.weights` を無視する。
- **Chung–Lu。** 頂点 i の重みは w_i ∝ (i+1)^(−1/(exponent−1)) で、`det_exp`/`det_log` で計算する。辺の両端点は、累積和に対する二分探索を通して、w に比例して引く。
- **Barabási–Albert。** m+1 頂点のクリークから始める。新しい頂点は、それぞれ m 個の相異なる対象につなぐ。対象は、辺の端点のリストから、次数に比例して標本抽出する。
- **敵対的ジェネレータ。** 重みは構成によって固定されており、`opt.weights` は使わない。
- **ランダム幾何グラフ**（`geometric`）。次の小節を参照。

### ランダム幾何グラフ

`random_geometric<W>(n, r, GeometricOptions{weights, integer_scale}, opt)` は、単位正方形上の G(n, r) を作る。

- **点。** `random_geometric_points(n, seed)` が点を返す。頂点 i は i 番目の点にあり、トポロジのサブストリーム（0）から `x = uniform01()`、続いて `y = uniform01()` として引く。したがって、点は、どの半径でも、どの重みの選択でも同じである。
- **辺。** u < v の `{u, v}` が辺であるための必要十分条件は、IEEE double で `dx*dx + dy*dy <= r*r`（dx = x_u − x_v、dy = y_u − y_v、FMA なし）となることである。この式が契約である。テストは、すべての対に対してこれを再計算できる。辺は (u, v) の辞書式順で並べるので、出力は近傍の見つけ方に依存しない。r = 0 では辺がない（2点が一致しなければならなくなる）。r ≥ √2（正方形の直径）では完全グラフになる。
- **近傍探索。** 点を k × k のセルに振り分ける。k·r ≤ 1 − 10⁻⁹ かつ k² ≤ n とする。すると、距離 r 以内のすべての対は、同じセルか隣のセルに入る（10⁻⁹ の余裕が、述語とセル座標の丸めを吸収する。証明は `geometric_cells_per_side` にある）。カウンティングソートで点をセルごとにまとめ、各頂点 v が自分の 3 × 3 ブロックを u < v について走査し、u による安定なカウンティングソートで対を順序どおりに並べる。期待時間は O(n + m)、メモリは O(n + m) である。
- **次数から半径へ。** `geometric_radius(n, d)` = √(d/(πn))（n < 2 では 0）。境界は無視する。期待される平均次数は (n−1)·F(r) である。ここで F(r) = πr² − (8/3)r³ + r⁴/2（r ≤ 1）は、正方形上の一様な2点が距離 r 以内にある確率である。実現される平均次数は d より相対的に ≈ 8r/(3π) だけ低い（n = 2·10⁴、d = 8 で約 1%）。
- **重み。** `GeometricWeights::kFromSpec`（CLI の `--geo-weights random`、既定）は、ほかのジェネレータと同じように、重みのサブストリームから、辺の順に `opt.weights` で引く。`kEuclidean`（CLI の `--geo-weights euclidean`）は辺の長さ √(dx² + dy²) を使う。`double` はそのまま保持し、整数の W は `round(length · integer_scale)`（最近傍、同じ近さは切り上げ）を使う。ごく近い点では 0 になりうるが、重み 0 は正当である。integer_scale は (0, 2⁵²] になければならない。`kEuclidean` では `opt.weights` は無視され、検証もされない。
- **`connected`。** `connect_components` は、ほかのジェネレータと同様に、連結用のサブストリーム（2）で成分をつなぐ。`kEuclidean` では、つなぐ辺の重みもその長さにするので、すべての重みがユークリッド距離になる。

### 敵対的な構成

- **Decrease-key**（k+1 頂点、始点 0、v_i = i）。辺は、重み 1 の `{0,1}`、i > 1 に対する重み 2k の `{0,i}`、1 ≤ i < j に対する重み 2(j−i)−1 の `{i,j}` である。
  - 重み 1 の辺の鎖から dist(v_j) ≤ j となり、どの辺 `{i,j}` も少なくとも |j−i| かかるので、dist(v_j) = j である。
  - したがって Dijkstra は、頂点を 1, 2, …, k の順に確定する。
  - v_i を確定すると、j > i の**すべての** v_j のキーが 2j−i から 2j−i−1 に下がる。
  - したがって、ヒープ版の Dijkstra は、ちょうど k(k−1)/2 = Θ(m) 回の decrease-key を行う。これは二分ヒープにとって最悪のケースであり、Fibonacci ヒープが想定した場面である。
- **長い経路**（n 頂点）。経路の辺 `{i,i+1}` の重みは 1、おとりのスキップ辺 `{i,i+2}` の重みは 3（迂回する2本の経路辺より大きい）。0 から n−1 への最短経路はただ1つで n−1 ホップあり、おとりはすべて緩和されたあとで棄却される。これは、コストがホップ数や再帰の深さとともに増える手法に負荷をかける。

### `connected` オプションと `connect_components`

`connect_components(el, weights, seed)` は次のように動く。

1. 成分を、最小の頂点の順に並べる。
2. i > 0 の各成分を、一様にランダムな前の成分につなぐ。つなぐ辺は、2つの成分の一様にランダムな成員同士の間に張る。これで成分上にランダムな木ができる。
3. ちょうど（成分数 − 1）本の辺を追加する。元の辺は先頭部分として、変更されずに残る。

ジェネレータは、連結用のサブストリーム（2）をシードとして渡す。構成上すでに連結なジェネレータは、この手順を省く。`count_components(el)` は成分を数える。孤立した頂点は、それぞれ1つと数える。

### カタログと汎用の入口

- `generator_catalog()` は、ジェネレータの型ごとに `{name, description, parameters}` を返す。
- `generate<W>(GraphSpec)` は `GraphSpec::type` で振り分ける。`GraphSpec` は `n`、`avg_degree`、`degree`、`attach`、`rows`、`cols`、`exponent`、`multi`、`road`、`radius`（`std::optional<double>`。未設定なら `geometric_radius(n, avg_degree)`）、`geometric`、`options` を持つ。未知の型は `std::invalid_argument` を投げる。
- `parse_multi_edge_policy("keep"|"drop"|"switch")` は、正則グラフのポリシーの CLI 名を変換する。
- `parse_geometric_weights("random"|"euclidean")` は、幾何グラフの重みの選択の CLI 名を変換する。

## 7. 小さなテスト用グラフ

`template <class W> std::vector<NamedGraph<W>> small_graph_suite();` は境界ケースのグラフを返す。ここで

```cpp
template <class W> struct NamedGraph { std::string name; EdgeList<W> graph; Vertex source = 0; std::string note; };
```

である。サイズを指定するビルダーは `sssp::gen::small` にある: `single_vertex()`、`no_edges(n)`、`two_vertices(w)`、`path(n, w)`、`star(n, w)`、`cycle(n, w)`、`complete(n, w)`（辺 {i,j} の重みは w·(1 + (i + 2j) mod 5)）、`binary_tree(n, w)`、`many_equal_paths(k)`。`long_shortest_path` は `adversarial_long_path(40)` である。残りのケースは、手で作った固定の辺リストで、スイートの中にしか存在しない。スイートはヘッダのみなので、どの重み型でも使える。グラフは次のとおり。

- **構造:** `single_vertex`、`no_edges`、`two_vertices`、`path`、`star`、`cycle`、`complete`、`binary_tree`、`disconnected`、`multi_edges`、`self_loops`。
- **重み:** `zero_weights`、`all_equal_weights`、`mixed_magnitudes`、`extreme_ratio`。
- **最短経路の形:** `long_shortest_path`、および `many_equal_paths`（k 個のダイヤモンドの鎖で、等しい最短経路が 2^k 本ある）。

差分テストと性質テストは、ケースを自分で列挙せずに、スイート全体を反復すべきである。

## 8. CLI（`sssp_gen`）

```sh
sssp_gen --list                                   # generator types and weight distributions
sssp_gen --type er --n 1000 --avg-degree 4 --weights uniform_int --seed 1 --out g.txt
sssp_gen --type road --rows 100 --cols 100 --weight-type double --format dimacs --out road.gr
sssp_gen --type regular --n 1000 --degree 3 --multi switch --weights wide_range --wide-log2 50 --out r.txt
sssp_gen --type er --n 10000 --avg-degree 1 --connected --weights zero_mixed --zero-prob 0.3 --out z.txt
sssp_gen --type geometric --n 100000 --avg-degree 8 --geo-weights euclidean --weight-type double --out geo.txt
sssp_gen --type geometric --n 5000 --radius 0.03 --weights log_uniform --lo 1 --hi 1000 --connected --out geo2.txt
```

| グループ | オプション |
|---|---|
| 出力 | `--format edgelist`（既定）または `dimacs`、`--weight-type uint64`（既定）または `double`、`--out` |
| 重み | `--weights <dist>`、`--lo`、`--hi`、`--mean`、`--zero-prob`、`--wide-log2` |
| 形 | `--n`、`--avg-degree`、`--degree`、`--attach`、`--rows`、`--cols`、`--exponent`、`--multi keep\|drop\|switch` |
| Road | `--jitter`、`--keep-prob`、`--shortcut-prob`、`--shortcut-radius`、`--int-scale` |
| Geometric | `--radius`（既定は `--avg-degree` から決まる）、`--geo-weights random\|euclidean`、`--int-scale`（road と共有。両方を設定する） |
| その他 | `--connected`、`--seed`、`--stats`（n、m、成分数を stderr に出す）、`--list`、`--help` |

終了コードは、成功が 0、使い方の誤り（未知のオプション、不正な数値）が 2、生成または IO のエラー（不正なパラメータ、未知の型、書き込めないファイル）が 1 である。同じ引数は、常にバイト単位で同一の出力を生む。

## 9. 既知の限界

- ジェネレータは n < 2^32、つまり 32 ビットの `Vertex` 型を必要とする。対のキーは2つの頂点を 64 ビットに詰める。
- `chung_lu` は m = round(n·avg_degree/2) ≤ ⌊⌊n(n−1)/2⌋ / 2⌋（整数除算）を必要とする。つまり頂点対の半分以下であり、重複の引き直しが速やかに終わるようにするためである。そうでなければ `std::invalid_argument` を投げる。棄却された引きが 1000 + 100·m 回に達すると `std::runtime_error` を投げる。
- `kSwitchToSimple` の `random_regular` は、単純な d 正則グラフが存在しない（n·d が奇数、または d > 0 で d ≥ n）ときに `std::invalid_argument` を投げる。有効な引数でも、スワップが 1000 + 100·(n·d/2) 回の提案以内に終わらないことがある。たとえば n = 3、d = 2 の自己ループ3つのように、決して修復できない構成もある。そうなったときは、ジェネレータは同じシード付きストリームから新しい構成を引いて、やり直す。100 個の構成がすべて失敗したときに限り、`std::runtime_error` を投げる。n ≤ 12 の (n, d) ごとに、300 個のシードのうち、1つの構成が失敗したのは最大で 31 個だったので、これは観測されていない。密なケース（d が n に近い）は、提案の大半が棄却されるので遅い。
- `random_geometric` の辺数は m ≈ n²·F(r)/2 なので、半径が小さくない場合（たとえば r ≥ √2 の完全グラフ）、出力は n の2乗になる。疎なグラフには `--avg-degree` を使うこと。
- `DimacsArcMode::kMergeReversePairs` は、対応の付いていない弧のハッシュ表を持つので、大きな道路網ではメモリが弧の数に比例して増える。
- `det_log` と `det_exp` の精度は数 ulp であり、正しい丸めではない。プラットフォーム間では厳密に再現できるが、libm とビットまで等しくはない。
- 整数の重みのパラメータは、`double` を通した変換が厳密になるように、2^53 以下でなければならない。
