# Dijkstra ベースライン

[English](baselines.md) | **日本語** | [简体中文](baselines.zh-CN.md)

ベースラインは、論文のアルゴリズムに関するあらゆる性能の主張の基準点になる。そのため、巧妙さではなく*公平さ*を優先して書いてある。走査ループは共有の1つだけ、ヒープはこのリポジトリ内で実装したもの（`std::priority_queue` は使わない）、ヒープの記憶領域は最初の push の前にすべて確保し、操作カウンタはどこでも同じものを使う。

公開 API: `include/sssp/baseline.hpp`。ソース: `src/baseline/`。

| 関数 | ヒープ（`src/baseline/heaps/`） | Decrease-key | 重みの型 | 時間 | ヒープの空間 |
|---|---|---|---|---|---|
| `dijkstra_binary_heap` | `BinaryHeap<K, V>`（遅延削除） | なし（再 push） | `uint64_t`、`double` | O((n + m) log m) | O(m) 個のエントリ、`num_arcs + 1` 個を予約 |
| `dijkstra_4ary_heap` | `IndexedDaryHeap<K, 4>` | あり | `uint64_t`、`double` | O((n + m) log₄ n) | O(n) |
| `dijkstra_fibonacci_heap` | `FibonacciHeap<K>` | あり（cut + cascading cut） | `uint64_t`、`double` | O(m + n log n) 償却 | O(n) のノードプール |
| `dijkstra_radix_heap` | `RadixHeap`（65 個のバケット） | あり（バケット移動） | `uint64_t` のみ | O(m + 64 n) | O(n) のノードプール + 65 個の先頭 |

`n` = 頂点数、`m` = 無向辺数。CSR の `num_arcs` は 2m から自己ループの分の弧を引いた数である。
テンプレートは、各 `.cpp` ファイルで `std::uint64_t` と `double` に対して明示的にインスタンス化している。

## 共有の走査ループ

`src/baseline/dijkstra_common.hpp`:

1. `init_result` は、`source >= n`（空グラフもここに含まれる）なら `std::out_of_range` を投げる。そうでなければ `dist = WeightTraits<W>::infinity()`、`parent = kNoVertex`、`parent_edge = kNoEdge`、`dist[source] = 0` を設定する。
2. 最小の `(d_u, u)` を取り出し、CSR の `u -> v` の弧を**すべて**走査する。確定済みの頂点への弧や自己ループも含める。「確定済み」ビットマップは参照しない。非負の重みでは、そうした弧が `dist[v]` を改善することはなく、読み飛ばすと実装間でカウンタが食い違ってしまうからである。
3. `nd = WeightTraits<W>::checked_add(d_u, w)` を計算し、`nd < dist[v]`（厳密な不等号）なら `dist[v] = nd`、`parent[v] = u`、`parent_edge[v] = arc.edge` を設定して、push または decrease-key を行う。

`dijkstra_binary_heap` は同じループを遅延ヒープで使う。改善された頂点は再度 push され、取り出したエントリが `dist[u] < d` なら古いエントリとして読み飛ばす。push は厳密な改善のときにしか起こらないので、頂点の最終的な距離はちょうど1回だけ push される。したがって、各頂点が走査されるのは多くても1回である。

## 同順位の扱い

ヒープはキーだけを比較する（`a < b`）。第2キーは無い。キーが等しいものの取り出し順は、ヒープの構造が決めたとおりになる。これは、与えられたグラフと始点に対しては決定的（CSR の弧の順は辺リストの順）だが、ヒープごとに異なる。その結果、次のようになる。

- 整数の重みでは、`dist` はすべてのベースラインで同一である。
- 最短路が複数あるとき、`parent` / `parent_edge` はベースラインごとに異なりうる。どれも正しい最短路木である。テストは親どうしを比べるのではなく、木を構造として検査する（`parent_edge` が存在し、`parent[v]` と `v` を結び、`dist[v] == dist[parent] + w` が成り立ち、親をたどると始点に着く）。
- 最初の改善が勝つ。同じ長さの経路があとから見つかっても、親は置き換えない。

## 重み 0

重み 0 の辺と重み 0 の閉路は許す。厳密な改善（`nd < dist[v]`）により、重み 0 の閉路が親の閉路を作ることはない。頂点の親は、その頂点が最後に改善される前に取り出された（確定した）頂点だからである。radix ヒープは、最後に取り出したキーと等しいキー（バケット 0）を受け付ける。これはまさに重み 0 の弧が作るキーである。

## 整数の重みとオーバーフロー

`std::uint64_t` では、`d_u + w >= 2^64 - 1`（和が `infinity()` に達する）のとき、`WeightTraits::checked_add` が `std::overflow_error` を投げる。この検査は、走査する弧の**すべて**で行う。そのため、ある候補の和がオーバーフローした時点で、その候補が `dist[v]` を改善しない場合でも（たとえば始点へ戻る巨大な辺）、実行は例外を投げる。これは意図的であり、結果が黙って誤ることはない。したがって、重み `2^64 - 1` の辺は、その端点のどちらかが走査されると必ず例外を投げる。

どの入力で例外が出るか（`double` でも同じ規則が成り立つ。後述）。

- **Dijkstra のベースラインはすべて**、まったく同じ入力で例外を投げる。到達可能な頂点 `u` のうち、距離 `d(u)` が表現可能なものに、重み `w` の弧があって `d(u) + w` が表現できないとき、かつそのときに限る。この条件はヒープに依存しない。頂点は最終距離でのみ走査され、距離が表現可能な到達可能頂点は、実行がすでに例外を投げていない限りすべて走査されるからである。
- **`bellman_ford`**（独立したオラクルで、`checked_add` を呼ばない）は、オーバーフローする候補を読み飛ばし、到達可能な頂点のうち距離が表現できないものがあるとき、かつそのときに限り例外を投げる。したがって、`bellman_ford` が例外を投げるなら Dijkstra もすべて例外を投げるが、逆は成り立たない。1本の辺 `{0, 1, 2^64 - 2}` では、`bellman_ford` は `{0, 2^64 - 2}` を返すが、Dijkstra はすべて、頂点 1 を走査するとき（0 へ戻る弧）に例外を投げる。

この2つの規則は、`tests/differential/test_baselines_differential.cpp`（uint64、独立した飽和演算のオラクル付き）と `tests/differential/test_double_overflow_baselines.cpp`（double、独立した丸めのオラクル付き）で固定している。

## 浮動小数点の重み

`double` では、`checked_add` は2つの有限な距離の IEEE-754 の加算（最近接への丸め。`-ffast-math` は決して有効にしない）であり、丸めた和が有限でなければ `std::overflow_error` を投げる。`+inf` は「到達不能」の目印なので、`+inf` に丸められる候補（たとえば `1e308 + 1e308` や、`DBL_MAX` の半 ulp である `DBL_MAX + 2^970`）を放置すると、厳密な `nd < dist[v]` で捨てられ、到達可能な頂点が到達不能と報告されてしまう（回帰テスト `tests/regression/cases/20261127-double-sum-overflow.json`）。代わりに吸収される和（`DBL_MAX + 2^969 = DBL_MAX`）は有限なので、例外は投げない。例外を投げる規則は前の節のものと同じで、「表現可能」を「丸めたあとで有限」と読み替える。`bellman_ford` は、独自の `!= +inf` の検査でこれを検出する。前提にしている性質は次のとおり。

- `w >= 0` なら `a + w >= a` なので、確定済みの頂点への弧がそれを改善することはなく、ヒープのキーは非減少の順に取り出される。
- `dist[v]` には、比較した `dist[parent[v]] + w` とまったく同じ値が格納されるので、木の検査 `dist[v] == dist[parent] + w` は*厳密に*成り立つ。
- 距離は、経路の順に沿った和である。別のアルゴリズム（たとえば部分経路の和を足す Floyd–Warshall）は丸めが異なりうるので、`double` の結果をアルゴリズム間で比較するには、重みが 2 進有理数でない限り許容誤差が要る（テストでは、厳密な比較には 1/16 の倍数を、一般の実数には 1e-12 の相対許容誤差を使う）。

radix ヒープは整数専用であり、`double` ではインスタンス化しない。

## 診断（`SsspResult::diag.ops`）

カウンタはすべて `SSSP_COUNT` を通して増やし、`SSSP_DIAGNOSTICS=0`（リリース・計測用ビルド）ではコードが消えて 0 のままになる。

| カウンタ | 定義 |
|---|---|
| `vertex_scans` | 隣接リストを走査した頂点（= 到達可能な頂点） |
| `edge_scans` | 検査した弧（走査した頂点の CSR の次数の和。無向辺は走査した端点ごとに1回、自己ループは1回数える） |
| `relaxations` | `dist[v]` の厳密な改善に成功した回数 |
| `additions` | `d_u + w` の評価（= `edge_scans`） |
| `comparisons` | `nd < dist[v]`（弧ごとに1回）+ pop ごとの古さの検査1回（二分ヒープのみ）+ ヒープ内部で行ったキー比較（`heap.comparisons()`） |
| `queue_pushes` | ヒープへの挿入（二分: 1 + `relaxations`、インデックス付き: 1 + 初めて到達した頂点の数） |
| `queue_pops` | extract-min の呼び出し（古いものを含む）（二分: = `queue_pushes`、インデックス付き: = `vertex_scans`） |
| `decrease_keys` | decrease-key の呼び出し（二分: 常に 0、インデックス付き: `relaxations + 1 - queue_pushes`） |

ヒープ内部のキー比較:

- `BinaryHeap`、`IndexedDaryHeap`: sift-up / sift-down のすべての `key < key`（d 分木の sift-down: 最小の子を探すのに最大 D-1 回、動かす要素との比較に1回）。
- `FibonacciHeap`: 挿入時（最初の根は除く）と decrease-key 時の最小値との比較、decrease-key 時の親との比較、consolidation 中の link ごとに1回、そして consolidation 後の根の中での最小値探索。
- `RadixHeap`: 再配置するバケットの最小値探索だけ（k 個の要素に対して k − 1 回）。バケットの添字は `bit_width(key XOR last)` で計算し、これは比較として数えない。このため、radix ヒープの比較回数は、比較ベースのヒープのものと直接は比べられない。`edge_scans` や時間で比べること。

ヒープ内部の引数検査（id の範囲、所属、キーの増加、radix の単調性）は数えない。

手計算の例（`tests/unit/baseline/dijkstra_test_support.hpp` にある古典的な 6 頂点のグラフ、始点 0）: 頂点走査 6 回、弧の走査と加算が 18 回、relaxation が 7 回。二分ヒープは push 8 回 / pop 8 回（うち古いもの 2 回）、インデックス付きヒープは push 6 回、decrease-key 2 回、pop 6 回。

## 割り当ての方針

実行中、ヒープ操作はメモリを確保しない。

- `BinaryHeap` は、最初の push の前に `num_arcs + 1` 個のエントリ（push 回数の上限）を予約する。
- `IndexedDaryHeap` は、ヒープ配列、位置の対応表、キー配列（大きさ n）をコンストラクタで確保する。
- `FibonacciHeap` は、頂点ごとに1つのノードと、n 個分の根の作業用バッファをコンストラクタで確保する。次数の表は固定の `std::array<id, 64>` である（n < 2³² なら最大次数 < 1.4405 log₂ n + 2 < 64）。
- `RadixHeap` は、頂点ごとに1つのノード（侵入型の双方向連結バケットリスト）をコンストラクタで確保する。バケットの先頭は固定の `std::array<id, 65>` である。

結果のベクタ（`dist`、`parent`、`parent_edge`）は、`init_result` が1回だけ確保する。

## ヒープの誤用

すべてのヒープは、空のヒープからの pop や `[0, n)` の外の id に対して `std::out_of_range` を投げ、すでにある id の push、無い id の decrease-key、より大きなキーへの decrease-key、および（radix では）最後に取り出したキーより小さいキーに対して `std::invalid_argument` を投げる。Dijkstra がこれらを起こすことはない。

## テスト

- ヒープ: `tests/unit/baseline/heaps/`（手計算した取り出し順、cascading cut、consolidation、バケットの挙動、カウンタ）と `tests/property/baseline/`（ランダムな操作列を、整列済みの `std::multiset` / `std::map` の参照と比較。シード {1, 7, 42, 2026, 0xC0FFEE}、操作後に不変条件の検証器を実行）。
- Dijkstra: `tests/unit/baseline/test_dijkstra_*.cpp`。共有の検査は `dijkstra_test_support.hpp` にある。手計算したグラフ（古典的なもの、重み 0 / 多重辺 / 自己ループ / 非連結 / 孤立頂点、重み 0 の閉路、単一頂点、2 進有理数の double）、不正な始点、オーバーフロー、カウンタ、そしてテスト専用の Floyd–Warshall によるランダムな多重グラフの検査（シード {1, 2, 3, 42, 2026}、n ≤ 30 の 125 グラフをすべての始点から、加えてシードごとに n = 200 のグラフ1つ）。
