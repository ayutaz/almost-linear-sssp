# 論文と実装の対応表(arXiv 2611.34911v1)

[English](PAPER_MAPPING.md) | 日本語

論文の各項目(定義・補題・定理・アルゴリズム・主要な定数)を、モジュール、テスト、Mode、状態、逸脱に対応づける。
表の形式は `.claude/skills/paper-reproduction/SKILL.md` に従う。行番号 `l.NNN` は `.cache/paper/2611.34911v1.txt`(pdftotext)のもの。論文本文は転載せず、言い換えで書く。

- 根拠の文書: `docs/research/01_preliminaries_primitives.md`(以下 01)、`02_ordered_hierarchy_local.md`(02)、`03_pivots_aux_outer.md`(03)、`thresholds.md`。
- 曖昧な点は `docs/research/open-questions.md`(Q-xx)に集めた。この表の「逸脱」は D-xx(表の下の一覧)。
- Status: TODO / RED / GREEN / VERIFIED / DEVIATES(SKILL の定義)。これに加えて、着手済みで RED の証跡がまだ揃っていない行を **IN_PROGRESS** とする(暫定の値。RED になったら置き換える)。
- Mode: A = `src/paper/`、B = `src/experimental/paper_experimental/`。「A, B」は、同じ部品を Mode B が引数を変えて使うことを表す(ADR 0003: 部品の共有は明示的な共通部品で行う)。

## 0. 命名の統一

3 つの解析文書のモジュール案で、名前や置き場所がずれていた箇所は、次のように統一した。

| 対象 | 01 の案 | 02 の案 | 03 の案 | 統一後 |
|---|---|---|---|---|
| split–findmin | `src/paper/primitives/split_findmin.hpp`(本実装と `SplitFindminRef`) | `src/paper/hierarchy/split_findmin.hpp`(`SegTreeSplitFindmin`) | — | **`src/paper/primitives/split_findmin.hpp`** に一本化。簡易版の名前は **`SplitFindminRef`**(境界の集合 + セグメント木)。`src/paper/hierarchy/` には置かない |
| 局所の基本の場合 | `src/paper/base_cases.cpp: local_base_hierarchy` | `src/paper/local/base_case.cpp: local_hierarchy_base_case` | — | **`src/paper/local/base_case.cpp: local_hierarchy_base_case`**(実装中の名前) |
| 外側の基本の場合 | `src/paper/base_cases.cpp: outer_base_case` | — | (名前の案なし) | **`src/paper/outer/base_case.cpp: outer_base_case`**。D-1 の版は `outer_base_case_d1`(ADR 0008) |
| 頂点対ごとの最軽量辺 | `mst.hpp: lightest_edge_per_pair` | (Alg 4 l.4 の中で使う) | (Alg 5 l.4、Alg 8 l.5 で使う) | 3 か所で使うので MST から切り離し、**`src/paper/primitives/edge_pairs.hpp: lightest_edge_per_pair`**。実装では `src/paper/primitives/mst.hpp` に置いた |
| 等重みの同時併合(Kruskal の併合木) | `base_cases` の中 | `base_case` と `assemble`(Y の木)の中 | Alg 6 の併合森 | **`src/paper/primitives/merge_events.hpp: build_merge_events`**(SKILL の例の行と同じ名前) |
| 零辺の縮約と Lemma 8.5 の復元 | — | — | (名前の案なし) | **`src/paper/outer/contract_restore.cpp: contract_zero_edges, restore_predecessor_tree`** |
| Sec 6 | — | — | (名前の案なし) | **`src/paper/pivots/name_routes.cpp`**(Alg 5)と **`src/paper/pivots/estimate.cpp`**(Alg 6) |
| Sec 7 | — | — | (名前の案なし) | **`src/paper/owner/`**: `aux_metric.cpp`(7.1–7.2)、`calendars.cpp`(7.3)、`owner_distances.cpp`(7.4、Alg 7)。`aux` という名前のディレクトリは Windows で予約名なので使わない |
| Sec 8 | — | — | (名前の案なし) | **`src/paper/outer/`**: `outer_states.cpp`(8.1)、`scales.cpp`(8.2)、`shortest_paths.cpp`(Alg 8 の入口) |
| Lemma 7.7 のセルで切られた走査 | — | `evaluate.cpp` に `PortionedEvaluator` を見越す | Alg 7 の中で使う | **`src/paper/hierarchy/evaluate.hpp: PortionedEvaluator`**(Sec 4 の評価器の公開 API。Alg 7 が使う) |
| テストの置き場所 | `tests/unit/paper/<subdir>/test_*.cpp` | `tests/unit/paper/test_*.cpp`、`tests/differential/test_*.cpp`、`tests/property/prop_*.cpp` | — | **`tests/unit/paper/<src と同じ subdir>/test_<module>.cpp`**、**`tests/differential/paper/test_*_vs_*.cpp`**、**`tests/property/paper/test_*_property.cpp`**。Mode B は `tests/*/paper_experimental/` |
| Dijkstra のテスト名 | `test_dijkstra.cpp` | — | — | 実装中の担当が使っている **`test_paper_dijkstra.cpp`**(baseline の Dijkstra のテストと区別する) |
| Mode B の引数 | `PaperParams{outer_exponent, local_exponent, local_coeff}` | `LocalParams{size_threshold, b, h, K, ...}` | (thresholds.md §7.2 の候補表) | **`src/experimental/paper_experimental/params.hpp`** に `PaperParams`(閾値)、`LocalParams`、`OuterParams` をまとめる。Mode A は `src/paper/model/params.hpp` の定数に固定 |
| 多重定義の記号 | `growth_limit`, `lca_block`, `sfm_block` | — | — | 01 の名前を採用(Q-09)。h は `hop_bound`、高さは `height`、乱数の範囲は `range` と書き分ける |

## 1. 定理・定義・補題・アルゴリズム

### Sec 1–2 主定理と前提

| Paper item | Statement (paraphrase) | Module / symbol | Tests | Mode | Status | Deviations |
|---|---|---|---|---|---|---|
| Sec 1, Theorem 1(l.91、証明 l.2160) | 無向多重グラフと有限の非負実数長、始点から、全頂点の正確な距離と到達頂点の最短路木を返す。期待時間 O(m + n log log log n)。ほとんど確実に停止し、停止した実行は正しい | `src/paper/outer/shortest_paths.cpp: paper_shortest_paths` | `tests/differential/paper/test_paper_sssp_vs_oracle.cpp` | A | TODO | D-1(変種)。計算量は測定で確かめる(SKILL §4) |
| Sec 2, l.284: L_k と τ(k) | L_k = ⌈log2(k+2)⌉、τ(k) = 1 + log log log(k+16)。τ は解析専用 | `src/paper/model/params.hpp: L` | `tests/unit/paper/model/test_params.cpp` | A | GREEN | none(解釈: L はすべての 64 ビットの k で定義し、L_0 = 1、2^64−1 での値は 65) |
| Sec 2, l.286–290: 大きさの前提 | 外側の構成は N > L_N^64、局所の構成は N > 2L_N^4 を仮定し、除いた大きさは絶対的に有界とする。局所の小さい場合の規約(頂点ごとに 1 状態、等重みの MST 辺を同時併合) | `params.hpp: outer_nontrivial, local_nontrivial, exceeds_scaled_pow, capped_pow, saturating_pow, kModeAParams` | `test_params.cpp`(局所は 209,952 で偽・209,953 で真、外側は N < 2^64 で常に偽、飽和で桁あふれしない、Mode B の引数がそのまま効く) | A, B | GREEN | Mode B は D-7 |
| Sec 2, l.292–300: 実数の計算モデル | 長さはコピー・比較・加算だけ。無限大は別のフラグ。差は非負の和 2 つの組。2^j 倍は j 回の倍加。床・log・ビット・アドレスとしての利用は禁止 | `src/paper/model/length.hpp: Len, SignedPair, EdgeKey, LenOps` | `tests/unit/paper/model/test_length.cpp` | A, B | GREEN | D-12(候補。実際に効いている: `Len::finite` が −0.0 を +0.0 にする)。ほかは none(解釈: 和や倍加が桁あふれしたら、巻き戻しや飽和ではなく `std::overflow_error` を投げる)。D-13 はここでは使っていない |
| Sec 2, l.301–306: 整数の語演算 | 添字とカウンタは語の演算 O(1)。⌊log2 i⌋ は表で作る。大きさの判定の積は上限で打ち切ってよい | `src/paper/model/ilog.hpp: IlogTable, floor_log2, ceil_log2` | `tests/unit/paper/model/test_ilog.cpp`(表を 1..2^20 で `bit_width − 1` と比べ、定義の不等式とも比べる) | A, B | GREEN | none。D-2 は使っていない: 語の版の `floor_log2` / `ceil_log2` はシフトと比較だけを使う(最大 64 回)。操作のたびに走る場所(`uniform_below` など)で、O(1) の表の代わりにこのループを使ってよいかは Q-51 |
| Sec 2, l.307–310: 乱数 | 一様な乱数ワード。0..h−1 は ⌈log2 h⌉ ビットの棄却法。乱数は Sec 6 だけで使う | `src/paper/model/random.hpp: WordRng, uniform_below, random_permutation, random_shuffle` | `tests/unit/paper/model/test_random.cpp`(生成器の参照列、h の剰余ではなく棄却、カイ二乗と順列の一様性) | A, B | GREEN | none(1 回の抽選ごとに使う対数は Q-51 を参照) |

### Sec 3 プリミティブ

| Paper item | Statement (paraphrase) | Module / symbol | Tests | Mode | Status | Deviations |
|---|---|---|---|---|---|---|
| Sec 3.1, Lemma 3.1(l.336) | DFS、BFS、トポロジカル整列は O(k+l)。安定な計数整列と定数座標の基数整列。BFPRT 選択は最悪線形 | `src/paper/primitives/search.hpp: Digraph, dfs_order, bfs_order, connected_components, topological_sort, counting_sort_stable, radix_sort_stable, select_bfprt` | `tests/unit/paper/primitives/test_search.cpp`(再帰の先行順とホップ距離と比較), `test_sort.cpp`(`std::stable_sort` と比較), `test_select.cpp`(`std::nth_element` と比較、比較回数が線形) | A, B | GREEN | none |
| Sec 3.1, Lemma 3.2(l.351) | Dijkstra は到達頂点を真の距離で確定する。半径 R で打ち切ると距離 ≤ R の頂点をちょうど確定する。二分ヒープ版と走査版 O(k²+l) | `src/paper/primitives/dijkstra.hpp: IndexedBinaryHeap, SearchArena, dijkstra(半径の引数), dijkstra_scan` | `tests/unit/paper/primitives/test_paper_dijkstra.cpp`(両版を uint64 と double で Floyd–Warshall と比較、半径の境界、Lemma 3.2 の操作回数) | A, B | GREEN | none |
| Sec 3.2, Lemma 3.3(l.409) | Fibonacci ヒープ: 挿入と減少は償却 O(1)、大きさ s からの取り出しは償却 O(1+log(s+1)) | `src/paper/primitives/fib_heap.hpp: FibHeap` | `tests/unit/paper/primitives/test_fib_heap.cpp` | A, B | GREEN | none(解釈: Alg 1 l.5 は成長ごとに新しいヒープを求め、ベースラインの `FibonacciHeap` は O(1) で捨てられないので、`discard()` が O(1) の専用ヒープにした。次数配列は証明どおりの大きさ 3 + 2⌈log2(s+1)⌉ で、s にはどのヒープの大きさも超えないプールの容量を使う。統合では使った枠だけに触れるので、上界も結果も変わらない)。D-6 は使っていない |
| Sec 3.2, Algorithm 1(l.376)+ Lemma 3.4(l.437) | 上限 B の成長パス。残した辺は共通の MST に属する。1 パスの費用と、残る商頂点 ≤ 1+2l'/B | `src/paper/primitives/mst.hpp: growth_pass` | `tests/unit/paper/primitives/test_growth_pass.cpp` | A, B | GREEN | none(解釈: 成長は頂点番号の昇順に始める。商頂点はその成分の最小の元頂点で名付ける。商の並行辺は残す) |
| Sec 3.2, Lemma 3.5(l.457) | 3 パス(B1 = max(2,⌈log2 L_k⌉)、B2 = L_k、∞)で MST を O(l + kτ(k))。木のボトルネック b_T は最小ミニマックス路の値 | `mst.hpp: mst_three_pass, three_pass_limits`、テスト用のオラクル `mst_kruskal_ref` | `tests/unit/paper/primitives/test_mst.cpp`; `tests/property/paper/test_mst_property.cpp`(Kruskal と辺集合まで完全に一致。総重みの一致より強い。ボトルネックの性質) | A, B | VERIFIED | none(解釈: 商に辺が無くなったら残りのパスを飛ばす。連結な入力では論文の「1 頂点ならパスは要らない」と同じ。出力の辺番号は、選択の印を O(m) で 1 回なめて昇順に集める) |
| (補完、01-A12。論文には無い)非連結な入力の全域森 | 非連結な入力に同じ 3 パスを当てると最小全域森になる。残る頂点数の上界は、成分数 c に対して c + 2l'/B と読む | `mst.hpp: spanning_forest_three_pass` | `test_mst.cpp`, `test_mst_property.cpp` | A(補完) | VERIFIED | none(論文には無い) |
| Sec 3.3 / Lemma 3.5 の利用: 順序付きの併合事象 | 等重みの MST 辺を同時に併合する Kruskal の併合木(子 ≥ 2、高さ = 併合の重み) | `src/paper/primitives/merge_events.hpp: build_merge_events` | `tests/unit/paper/primitives/test_merge_events.cpp` | A, B | TODO | none |
| Alg 4 l.4、Alg 5 l.4、Alg 8 l.5: 頂点対ごとの最軽量辺 | 同じ頂点対の辺のうち (長さ, 添字) が最小の 1 本を残し、ループを捨てる。計数整列で O(m+N)(01 §2.4 の補完) | `src/paper/primitives/mst.hpp: lightest_edge_per_pair`(安定な計数整列 2 回、O(n+m)。予定の `edge_pairs.hpp` ではなく `mst.hpp` に置いた。§0 を参照) | `tests/unit/paper/primitives/test_mst.cpp` | A, B | GREEN | none。Alg 5 l.4 と Alg 8 l.5 の呼び出し側はまだ main に無い |
| Sec 3.3, Lemma 3.6(l.481) | union-find の初期化と Q 回の操作で O(Q + kτ(k)) | `src/paper/primitives/dsu.hpp: DisjointSets` | `tests/unit/paper/primitives/test_dsu.cpp`(素朴なラベルとの差分), `tests/property/paper/test_dsu_property.cpp` | A, B | VERIFIED | none(解釈: ランクが等しいときは、第 1 引数の根が残る) |
| Sec 3.4, Lemma 3.7(l.521) | オイラー巡回の表で LCA を O(1)、前処理 O(k) | `src/paper/primitives/forest_queries.hpp: LcaTable, OrderedTree, paper_block_size` | `tests/unit/paper/primitives/test_lca.cpp`(`BinaryLiftingRef` との差分), `tests/property/paper/test_forest_queries_property.cpp` | A, B | VERIFIED | D-16(テスト専用のブロック長)。ほかは none(解釈: 子は頂点番号の昇順に並べる)。D-5 は使っていない |
| Sec 3.4, Lemma 3.8(l.550) | 印付けで B 未満の小木に分け、形の表・跳躍・ladder で level ancestor を O(1) | `forest_queries.hpp: LevelAncestor` | `tests/unit/paper/primitives/test_level_ancestor.cpp`(`BinaryLiftingRef` との差分), `test_forest_queries_property.cpp` | A, B | VERIFIED | D-16(テスト専用のブロック長)。ほかは none(解釈: 印付けを後行順ではなく逆の先行順で行う。付く印は同じ)。D-5 は使っていない |
| Sec 3.4, Lemma 3.9(l.570) | 根付き森で LCA と LA を O(1)、前処理 O(k)。違う木は成分で判定 | `forest_queries.hpp: ForestQueries`、テスト用の参照 `BinaryLiftingRef` | `test_lca.cpp`, `test_level_ancestor.cpp`, `tests/property/paper/test_forest_queries_property.cpp` | A, B | VERIFIED | D-16(テスト専用のブロック長)。注: 論文の B は k < 2^24−1 ではすべて 2 なので、現実の大きさでは疎表と跳躍表は O(k log k) のようにふるまう。線形性は漸近的な主張で、計測していない。D-5 は使っていない |
| Sec 3.5, Lemma 3.10(l.592) | dyadic な表現: 初期化と全 split で O(k(1+log k))、decrease・query・記述子の参照は O(1) | `src/paper/primitives/split_findmin.hpp: DyadicSplitFindmin` | `tests/unit/paper/primitives/test_split_findmin.cpp`, `tests/property/paper/test_split_findmin_property.cpp` | A, B | VERIFIED | none(解釈。論文が決めていない API の選択: 区間 ID は初め 0、s 回目に成功した split で左が 2s−1、右が 2s(Q-14: 両側とも新しい ID)。既存の境界での split は false を返す。キーを上げる decrease は release では無視し、debug では assert で止める。NaN のキーは拒む。k ≤ 2^31) |
| Sec 3.5, Lemma 3.11(l.625) | 2 段のブロック化で、初期化・全 split・q 回の decrease・q' 回の query が O(q + q' + kτ(k))。同順位は元の添字 | `split_findmin.hpp: BlockedSplitFindmin, OneLevelBlockedSplitFindmin, SplitFindmin` | `test_split_findmin.cpp`; `tests/property/paper/test_split_findmin_property.cpp` | A, B | VERIFIED | none(解釈: 2 段目のブロック長は B' = max(2, ⌈log2(B+1)⌉)。短い最後のブロックにも同じ値を使う。ほかに Lemma 3.10 の API の選択)。D-18: ブロックの番号を整数の除算で求める |
| (テスト支援。論文には無い)参照用の split–findmin | 境界の集合とセグメント木、操作ごとに O(log k)。差分テストの参照にだけ使う | `split_findmin.hpp: SplitFindminRef` | `test_split_findmin.cpp` | A(テスト支援) | GREEN | none(論文の手順には使わない。論文のモジュールはどれも使っていないので、D-3 は効いていない) |

### Sec 4 順序付き階層の評価

| Paper item | Statement (paraphrase) | Module / symbol | Tests | Mode | Status | Deviations |
|---|---|---|---|---|---|---|
| Sec 4, Definition 4.1(l.646) | 有向状態グラフ、物理頂点への射影、葉が状態の森。高さは正で根に向かって非減少、子 ≥ 2、軽い弧の商と根の商は非巡回、子は長さ ≤ c0·h のリンクの木でつながる(02 §2.1 の H1–H6) | `src/paper/hierarchy/ordered_hierarchy.hpp: OrderedHierarchy, StateGraph, PhysicalLink, validate_ordered_hierarchy` | `tests/unit/paper/hierarchy/test_ordered_hierarchy_validate.cpp`(H1–H6 の違反を 1 つずつ検出) | A, B | GREEN | none(解釈: H6 はリンクと弧の整合として検査する。その核心である「証明だけに使うリンクは緩和しない」は、リンクと弧を別の型にすることで守る) |
| Sec 4, Definition 4.2(l.659) | 重み付きの物理的な源と実行可能な初期申し出。タイトな証拠を持つ状態が保護され、要求された目標ごとに保護された状態がある | `ordered_hierarchy.hpp: SeededOrderedHierarchy, PhysicalSource, StateOffer, validate_seeded_ordered_hierarchy` | `test_ordered_hierarchy_validate.cpp` | A, B | GREEN(検証) | none。予定していたテスト用のオラクル `protected_states` と `tests/unit/paper/hierarchy/test_protected_states.cpp` はまだ main に無い(TODO) |
| (補完)参照階層: 閾値ごとの SCC(02 §2.11) | 任意の有向状態グラフに Def 4.1 を満たす階層(c0 = 1)を作る。論文の構成ではなく、評価器を Sec 5 と切り離して試すためのテスト支援 | `src/paper/hierarchy/reference_hierarchy.hpp: build_threshold_scc_hierarchy, physical_graph_of_arcs` | `tests/unit/paper/hierarchy/test_reference_hierarchy.cpp` | A(テスト支援) | GREEN | none(論文の手順には使わない) |
| Sec 4.1, Lemma 4.3(l.674, l.690) | 正規化は O(M+N)。単位 ν は親で 2 倍以上、Σ S/ν = O(N)。親子は倍数で整列するか 2S(C) < ν(P) | `src/paper/hierarchy/normalize.hpp: normalize, NormalizedHierarchy, LightArc` | `tests/unit/paper/hierarchy/test_normalize.cpp` | A, B | GREEN | none(解釈: Q-23 では、縮約するのは正規化のリンクだけで、縮約した節点は最上位の節点の質量を持つ。Q-24 では、リンクの無い節点は ν = h。Q-50(要合意)では、縮約後の森は高さ = ν、接続定数 2c0 とする。補完: 重みの型を超える倍加は飽和させる) |
| Sec 4.2, Lemma 4.4(l.731) | 子の代表を逆の終了順に巡る物理的な歩道の長さ ≤ 2ΣS(i) + 2dν(X) + 2S(X) | `src/paper/hierarchy/child_order.hpp: order_children, build_link_forest, LinkForest` | `tests/unit/paper/hierarchy/test_child_order.cpp` | A, B | GREEN | D-15(候補)。ほかは none(解釈: Q-25 では、仮の代表を F の先行順で最初の状態にする)。巡回の開始位置は l.773–777 の O(1) の表で引き(D-4 は使っていない)、森の距離は Lemma 3.9 の O(1) の LCA/LA で求める(D-5 は使っていない) |
| Sec 4.2, Lemma 4.5(l.758) | 子のトポロジカル順と間隔 g。Σ g/ν = O(N)、全体 O(M+N) | `child_order.hpp: order_children` | `test_child_order.cpp` | A, B | GREEN | Lemma 4.4 と同じ(D-15 は候補。D-4 と D-5 は使っていない) |
| Sec 4.3(l.781–832): スロット、期限、run、補足の申し出 | 窓 [a,b) と閉じた期限 c、開く時刻 D0 と満了 D0+S(C)、原点の整列規則、未登録の run の分割、補足の申し出 max{b, z+S(i)+g(i,j)} | `src/paper/hierarchy/evaluate.cpp`(内部) | `tests/unit/paper/hierarchy/test_evaluate_basic.cpp` | A, B | TODO | none |
| Sec 4.3, Algorithm 2(l.834) Scan | live な節点の子をトポロジカル順に走査し、due な run を分割・登録し、子のスロットを進め、満了で引退させる | `evaluate.cpp: scan`(内部)、`evaluate.hpp: PortionedEvaluator` | `test_evaluate_basic.cpp`, `tests/unit/paper/hierarchy/test_evaluate_invariants.cpp` | A, B | TODO | D-14(テスト専用の切り替え) |
| Sec 4.3, Algorithm 3(l.866) Evaluate | 正規化・子の順序の後、根の商の順に根を開き、スロットごとに Scan を呼ぶ | `evaluate.hpp: evaluate, EvalOptions` | `tests/differential/paper/test_evaluate_vs_dijkstra.cpp` | A, B | TODO | none |
| Sec 4.3, Lemma 4.6(l.885) | 保護された状態は δ を出力し、ほかの出力はすべて δ 以上 | `evaluate.cpp`(デバッグ時の検査 I1–I8、02 §2.8) | `test_evaluate_invariants.cpp`, `test_evaluate_vs_dijkstra.cpp` | A, B | TODO | none |
| Sec 4.3, Lemma 4.7(l.922) | 評価は O(M + Nτ(N) + Q) | `evaluate.cpp`(`Diagnostics.counters` に訪問回数) | `tests/property/paper/test_evaluate_counts_property.cpp` | A, B | TODO | D-3, D-4, D-5(簡易版を使う間) |

### Sec 5 局所階層

| Paper item | Statement (paraphrase) | Module / symbol | Tests | Mode | Status | Deviations |
|---|---|---|---|---|---|---|
| Sec 5.1(l.960), Definition 5.1(l.975) | MST で支配される辺を消し、印付けとピボット集合 A、回廊の最大辺の切断で owner セルを作る。各セルにピボット 1 つ、T−A の成分 < b、h = b+1、R_v ≤ h·r_v | `src/paper/local/owner_decomposition.cpp: owner_decomposition(G, mst, b)` | `tests/unit/paper/local/test_owner_decomposition.cpp` | A, B(外側は b = L^64、Sec 7.1) | TODO | none |
| Sec 5.1, Lemma 5.2(l.993) | セル間の辺は w ≥ max{r_x, r_y, σ}、全辺で w ≤ h(r_x+r_y) + Nσ | `owner_decomposition.cpp` | `test_owner_decomposition.cpp` | A, B | TODO | none |
| Sec 5.1, Lemma 5.3(l.1011) | 任意のピボット q について b_T(v,q) ≥ r_v | `owner_decomposition.cpp` | `test_owner_decomposition.cpp` | A, B | TODO | none |
| Sec 5.2(l.1020–1055) 辺の分類・shell・窓 | macro は w ≤ N³σ のセル間辺、残りは radial。K = 8h 以上の最小の 2 冪、near は比 ≤ K²。inward / outward の 2 殻の窓ごとにコピーを作る | `src/paper/local/radial_states.cpp: build_radial_states(owner, RadialParams)` | `tests/unit/paper/local/test_radial_states.cpp` | A, B(外側は座標 J・K = 16、Sec 8.1) | TODO | D-13(候補) |
| Sec 5.2, Definition 5.4(l.1056) | radial 状態表現 = 上の有向グラフ + 射影 + 各頂点の native inward 入口。座標・比 K・分類が引数 | `radial_states.hpp: RadialStateGraph` | `test_radial_states.cpp` | A, B | TODO | none |
| Sec 5.2, Lemma 5.5(l.1061) | 最短路の radial 区間で、両端の半径がともに M/K 以下にはならない | `radial_states.cpp`(テストで確認) | `tests/property/paper/test_radial_states_property.cpp` | A, B | TODO | none |
| Sec 5.2, Lemma 5.6(l.1074) | 極大な radial run を最小半径で分けると、前半・後半の振れ幅は K 倍以内 | 同上 | `test_radial_states_property.cpp` | A, B | TODO | none |
| Sec 5.2, Lemma 5.7(l.1084) | 状態 O(N)・弧 O(m+N)。native inward から全測地線を実現。macro を除き block を縮めると DAG | `radial_states.cpp` | `test_radial_states.cpp`(コピーの最小 = 物理 Dijkstra、block DAG) | A, B | TODO | none |
| Sec 5.3(l.1109), Lemma 5.8(l.1131) | I0・O0・I*・O* の伝播と θ。B0 有限の状態で r ≤ θ、B0 ≤ θ ≤ 3hK²B0 | `src/paper/local/threshold_joins.cpp: compute_joins` | `tests/unit/paper/local/test_threshold_joins.cpp` | A, B | TODO | none |
| Sec 5.3, Lemma 5.9(l.1152) | 重み ≤ W の結合レコードの無向成分を状態に制限したものは、弧長 ≤ W の SCC と一致(< W 版も) | `threshold_joins.cpp` | `tests/property/paper/test_threshold_joins_property.cpp`(全 W で Tarjan と照合) | A, B | TODO | none |
| Sec 5.4(l.1182–1243) 組み立て | s_U・t_U、aligned 判定 s_U ≥ σ_U/N^c(c = 10)、支え [·, N^16σ] の大域倍加格子、Y の正確な Kruskal 木、階級ごとの同時併合(c0 = 2) | `src/paper/local/assemble.cpp: assemble_local_hierarchy` | `tests/unit/paper/local/test_local_assemble.cpp` | A, B | TODO | D-13(候補)。Mode B は D-8, D-11 |
| Sec 5.4, Algorithm 4 l.2–5(l.1244、Sec 2 l.286–290) 基本の場合 | k ≤ 2L_k^4 なら、頂点ごとに 1 状態・全非ループ辺を両向きの弧・最軽量辺の MST の等重み辺を同時併合した Kruskal 階層(c0 = 1) | `src/paper/local/base_case.hpp: local_hierarchy_base_case` | `tests/unit/paper/local/test_local_base_case.cpp` | A, B | GREEN | none(解釈: 非連結な入力では森になる(補完)。論文はこの定数の大きさの場合の MST の求め方を決めていないので、比較整列の Kruskal を使い、同順位は辺番号で決める。閾値の分岐とその FallbackRecord は、次の行の Alg 4 の入口が受け持つ) |
| Sec 5.4, Algorithm 4 l.6–12 + Lemma 5.10(l.1265) | 局所階層の入口。閾値で分岐し、Def 4.1 を満たす階層(c0 ≤ 2)、全頂点の native entry、O(m + N log log(N+4)) | `src/paper/local/local_hierarchy.cpp: local_hierarchy` | `tests/unit/paper/local/test_local_hierarchy.cpp`, `tests/differential/paper/test_local_hierarchy_vs_dijkstra.cpp` | A, B | TODO | none |
| Sec 5.5, Definition 5.11(l.1303) | 置き換え可能な原子(cutoff c): 内部が連結、全域木の辺 < 2c、親の高さ ≥ c | `src/paper/local/atom_replacement.hpp: ReplaceableAtom` | `tests/unit/paper/local/test_atom_replacement.cpp` | A, B | TODO | none |
| Sec 5.5, Definition 5.12(l.1308) | 分割階層: 葉が状態の分割(1 状態か原子)で、部分の上で Def 4.1 を満たす | `atom_replacement.hpp: PartitionHierarchy` | `test_atom_replacement.cpp` | A, B | TODO | none |
| Sec 5.5, Lemma 5.13(l.1320) | 原子の局所階層を高さ c で切ると、根 1 つ・高さ ≤ c・根のリンク < 2c の順序付き階層 | `atom_replacement.cpp: cap_local_hierarchy` | `test_atom_replacement.cpp` | A, B | TODO | none |
| Sec 5.5, Lemma 5.14(l.1342) | 置き換えてもタイトな証拠は保たれる(native 入口、任意のコピーからの出口) | `atom_replacement.cpp: substitute_atoms` | `test_atom_replacement.cpp`(置き換えの前後で評価結果が一致) | A, B | TODO | none |
| Sec 5.5, Lemma 5.15(l.1355) | 原子が b 状態未満なら、全ての置き換えで順序付き階層が得られる。O(m + N log log(b+4)) | `atom_replacement.cpp: substitute_atoms` | `test_atom_replacement.cpp` | A, B | TODO | none |

### Sec 6 ピボット経路の証明書

| Paper item | Statement (paraphrase) | Module / symbol | Tests | Mode | Status | Deviations |
|---|---|---|---|---|---|---|
| Sec 6, Definition 6.1(l.1368) | t 個のピボットのループの無い多重グラフ。証明書は λ_e ≤ a_e で、長さ ≤ Pλ_e の物理 walk がある。辞書は整列済みの重みの一覧 | `src/paper/pivots/name_routes.hpp: PivotRoutingInstance, RouteNames, check_route_names` | `tests/unit/paper/pivots/test_name_routes.cpp`(G_P の Floyd–Warshall), `tests/property/paper/test_name_routes_property.cpp`(Floyd–Warshall とベースラインの Bellman–Ford) | A, B | VERIFIED | none(解釈: ループは `std::invalid_argument` で拒み、呼び出し側が取り除く。証明書の walk はすべてインスタンスの辺をつないだものなので、証明書は dist_{G_P}(x,y) ≤ Pλ_e で確かめる) |
| Sec 6.1(l.1382–1392) 前処理と中央値の絞り込み | 平行辺を最軽量の代表に置き換える。群の中で P·a* ≥ Z の最小の a* を BFPRT と倍加で探す(中央値と等しいキーは捨てる) | `name_routes.hpp: median_prune`、代表の選択は `name_routes.cpp` の中 | `test_name_routes.cpp`(乱数のキーで総当たりと比較) | A, B | GREEN | none(解釈: 等重みの平行辺は最小の辺番号で代表する。端点の対ごとの最軽量辺は、§0 で予定した共通部品ではなく `name_routes.cpp` の中で選ぶ) |
| Sec 6.1, Algorithm 5(l.1394) NameRoutes | 残りが 32t より多い間、16t 本を標本にして推定し、半数以上が合格なら受理、そうでなければ乱数を引き直す。最後は自分の重みを証明書にする | `name_routes.hpp: name_routes, NameRoutesParams, NameRoutesConstants` | `test_name_routes.cpp`; `tests/property/paper/test_name_routes_property.cpp`(多数の乱数種で証明書を Floyd–Warshall と Bellman–Ford で確かめ、受理率を計測) | A, B | VERIFIED | Mode B は D-9(`EstimateParams::h_log2`)、D-10(`retry_cap`。既定では切ってあり、打ち切ったら FallbackRecord を残す)、D-17(`sample_factor`、`residual_factor`、`p_log2`)。Mode A では none(解釈: 「半数以上」は 2·\|合格\| ≥ \|U\|(Q-44)。辞書は (重み, レコード番号) のマージ整列で並べる。辞書の大きさは数えるだけ(Q-40)) |
| Sec 6.1, Lemma 6.2(l.1424) | g = 4⌈log2(t+2)⌉+3 のとき、critical な辺は ≤ 2t 本 | (テストだけ) | `tests/property/paper/test_name_routes_property.cpp`(総当たり) | A | GREEN | none |
| Sec 6.1, Lemma 6.3(l.1439) | 一様な 16t 本の標本で d_S > g·a_e となる辺の期待数 ≤ h/8 | `name_routes.hpp: partial_shuffle`(一様な k 個の部分集合) | `test_name_routes.cpp`(`partial_shuffle` の一様性), `test_name_routes_property.cpp`(証拠の無い辺の割合の計測平均を 1/8 と比べる) | A | GREEN | none |
| Sec 6.2, Algorithm 6(l.1454) Estimate | 標本の Kruskal 併合森と帯、商 Q_a、尺度ごとの刈り込み付き多中心探索、分割木と表の値 4Δ+a | `src/paper/pivots/estimate.hpp: estimate, EstimateParams, EstimateResult, TableBound, EstimateTrace, check_estimate` | `tests/unit/paper/pivots/test_estimate.cpp`, `tests/property/paper/test_estimate_property.cpp`(どちらも Floyd–Warshall と比較) | A, B | VERIFIED | none(`estimate.cpp` の冒頭の解釈 I-1〜I-7: I-1(Q-06、Q-43)2 冪の倍数は正確に扱い、整数の重みは最大 512 ビットの多語整数で持つ。double では倍加は正確で加算は丸める。I-2 ヒープは尺度を掛けない距離で並べ、H·d は半径の判定でだけ作る(比較も加算の順も同じ)。I-3(Q-41)最上位の尺度の部分は商の成分(検査し、違えば `std::logic_error`)。I-4 標本の等しい重みは辺番号で順位を付ける。I-5 計数整列の配列は帯ごとに \|V(Q_a)\| の大きさ。I-6 乱数は粗い尺度から細かい尺度へ引く。I-7 2 分の併合森は、予定した共通の併合事象の部品ではなく、`estimate.cpp` の中で `DisjointSets` を使って作る) |
| Sec 6.2, Lemma 6.4(l.1509) | 商の路の持ち上げ: a ≤ β ≤ d ≤ d_S ≤ d + a | `estimate.cpp` | `test_estimate_property.cpp`(標本とすべての商で Floyd–Warshall、ミニマックスのボトルネック) | A, B | VERIFIED | none |
| Sec 6.2, Lemma 6.5(l.1530) | 1 つの尺度で距離 d の 2 点が分かれる確率 ≤ H_t(d/Δ + 1/H) | (テストだけ、統計) | まだ書いていない(P4-RANDOMIZED で予定) | A | TODO | none |
| Sec 6.2, Lemma 6.6(l.1542) | 刈り込み付きの探索は各頂点の最初に覆う中心を返す。期待時間 O((v+e) log²(t+2)) | `estimate.cpp`(内部の `first_covering_centers`) | `test_estimate_property.cpp`(最初に覆う中心を、Floyd–Warshall の商の距離の上で総当たりと比較。走査の回数を計測し、H_t で重みを付けた上界と比べる) | A, B | VERIFIED | none(解釈 I-2) |
| Sec 6.2, Lemma 6.7(l.1560) | d_S ≤ Z は決定的、E[Z \| S] ≤ 256u²·d_S、表は O(t log(t+2)) 項 | `estimate.hpp: estimate, check_estimate` | `test_estimate.cpp`(Z ≥ d_S を毎回)、`test_estimate_property.cpp`(乱数の標本で Z ≥ d_S。Z/d_S の計測平均を 256u² と比べる。表の大きさ) | A, B | VERIFIED | Mode B は D-9。注: double の重みでは Z = 4Δ + a は丸めた和で、正確な値をわずかに下回りうる。正確さを確かめるのは整数の重みと 2 進の有理数だけ |
| Sec 6.3, Lemma 6.8(l.1605) | t = O(1+N/L^64) で全辺に証明書、辞書 O(t log²(N+2)) 項、ほとんど確実に停止、期待 O(M+N) | `name_routes.hpp: name_routes` | `test_name_routes_property.cpp`(全辺で λ_e ≤ a_e と dist ≤ Pλ_e を Floyd–Warshall と Bellman–Ford で確かめる、多数の乱数種。受理確率 ≥ 1/2 を計測) | A, B | GREEN | none。期待 O(M+N) の時間は計測していない |

### Sec 7 owner ピボットへの補助距離

| Paper item | Statement (paraphrase) | Module / symbol | Tests | Mode | Status | Deviations |
|---|---|---|---|---|---|---|
| Sec 7.1(l.1642–1651) 辺の分類 | 外側は b = L_N^64 の owner 分解。射影長 A_e = R_x + w + R_y を Lemma 6.8 に渡し、d_e = Pλ_e。2d_e ≤ w なら radial、それ以外は macro。セル内の辺は d_e = 0 | `src/paper/owner/aux_metric.cpp: classify_interowner_edges` | `tests/unit/paper/owner/test_aux_metric.cpp` | A, B | TODO | none |
| Sec 7.1, Lemma 7.1(l.1653) | 切断された回廊辺は macro。macro は λ/(2h+1) ≤ w < 2Pλ を満たす名前を持つ | `aux_metric.cpp` | `test_aux_metric.cpp` | A, B | TODO | none |
| Sec 7.2(l.1662), Lemma 7.2(l.1668) | ピボットを 1 点に同一視し、radial 辺を w+d_e で残したグラフ G' と木 T'。J_v は根までの距離で r_v ≤ J_v ≤ R_v ≤ h·r_v、radial 辺で \|ΔJ\| ≤ 3w/2 | `aux_metric.cpp: build_aux_graph`、テスト用の `owner_coordinates_ref`(Dijkstra) | `test_aux_metric.cpp` | A, B | TODO | none |
| Sec 7.3(l.1692–1710) 片・格子・カレンダー・原子 | 重なり [r/4, 4R] で木片とグラフ片を作り、原点の 7 オクターブ下から格子、32q ≤ r < 64q、U+ > R、カレンダーのセル、片をまたぐマーカー弧 U+_u → 16q_v、原子は w' < 2max(q_u,q_v) の成分 | `src/paper/owner/calendars.cpp: build_calendars, find_atoms` | `tests/unit/paper/owner/test_calendars.cpp` | A, B | TODO | none |
| Sec 7.3, Lemma 7.3(l.1711) | 格子点の総数 O(N)、各頂点の幅と上側マーカーは O(1 + log log(b+4)) | `calendars.cpp` | `test_calendars.cpp` | A, B | TODO | none |
| Sec 7.3, Lemma 7.4(l.1736) | カレンダーのセルの総数とマーカーの併合は O(N) | `calendars.cpp` | `tests/property/paper/test_calendars_property.cpp`(セル数 / N) | A, B | TODO | none |
| Sec 7.3, Lemma 7.5(l.1751) | 原子は 1 つの片と 1 つの幅に属し、< b 頂点。定義辺は 1 辺 < 2q の全域木 | `calendars.cpp: find_atoms` | `test_calendars.cpp` | A, B | TODO | none |
| Sec 7.3, Lemma 7.6(l.1762) | 項数 O(N)、構築 O(m + N log log(b+4))、マーカー弧は数値時間を増やす | `calendars.cpp` | `test_calendars.cpp`(マーカー DAG が非巡回で増加) | A, B | TODO | none |
| Sec 7.4(l.1780–1841) 族・run・原子の訪問 | 族の地平線と起動リスト、族の中の run 分割と補足の申し出 max{b0, z+L_i+g(i,j)}、原子の局所表現と出口の抑制、根スロットの Scan | `src/paper/owner/owner_distances.cpp` | `tests/unit/paper/owner/test_owner_distances.cpp` | A, B | TODO | none |
| Sec 7.4, Algorithm 7(l.1843) OwnerDistances | カレンダーのマーカーをトポロジカル順に処理し、族を登録・訪問して全 J_v を出す | `owner_distances.cpp: owner_distances(T', G')` | `tests/differential/paper/test_owner_distances_vs_dijkstra.cpp` | A, B | TODO | none |
| Sec 7.4, Lemma 7.7(l.1871) | 外部の格子で切られたスロットでも、保護された接頭辞は自分のセルで出力される。各局所走査は高々 2 回 | `src/paper/hierarchy/evaluate.hpp: PortionedEvaluator` | `tests/unit/paper/hierarchy/test_portioned_scan.cpp` | A, B | TODO | none |
| Sec 7.4, Lemma 7.8(l.1919) | 許された向きと z ≤ R_u、z+w' ≤ R_v を満たす原子間の申し出は、行き先セルの開始マーカーより前に届く | `owner_distances.cpp`(到着ログ) | `test_owner_distances.cpp` | A, B | TODO | none |
| Sec 7.4, Lemma 7.9(l.1944) | 名前が与えられれば、Alg 7 は全 J_v を決定的に O(m + Nτ(N)) で求める | `owner_distances.cpp` | `test_owner_distances_vs_dijkstra.cpp` | A, B | TODO | none |

### Sec 8 外側の階層と単一始点アルゴリズム

| Paper item | Statement (paraphrase) | Module / symbol | Tests | Mode | Status | Deviations |
|---|---|---|---|---|---|---|
| Sec 8.1(l.1966–2001) 外側の状態グラフ | w > J_x + J_y + d_e の radial 辺を削除し、座標 J・K = 16 で 2 殻の構成。ピボット以外の各頂点に特別状態 v_in, v_out。ブロック順の伝播で θ、結合レコード | `src/paper/outer/outer_states.cpp: build_outer_states`(`radial_states` と `threshold_joins` を使う) | `tests/unit/paper/outer/test_outer_states.cpp` | A, B | TODO | none |
| Sec 8.1, Lemma 8.1(l.2002) | 残った radial 辺は w ≤ 2(J_x+J_y)。振れ幅は 16 倍以内。全測地線を実現し、macro を除いた商は非巡回 | `outer_states.cpp` | `test_outer_states.cpp`, `tests/property/paper/test_outer_states_property.cpp` | A, B | TODO | none |
| Sec 8.1, Lemma 8.2(l.2027) | B0 と θ は同時に有限で、r ≤ θ、B0 ≤ θ ≤ 4K²h·B0 | `outer_states.cpp` | `test_outer_states_property.cpp` | A, B | TODO | none |
| Sec 8.1, Lemma 8.3(l.2053) | 重み ≤ W の結合で、弧長 ≤ W の強連結な状態成分を正確に記述(< W 版も) | `outer_states.cpp` | `test_outer_states_property.cpp`(Tarjan と照合) | A, B | TODO | none |
| Sec 8.2(l.2076–2104) 数値尺度と原子の置換 | 名前を最小の macro 重みに置き換え、F = 2P(2h+1)(8K²h)² より大きい最小の 2 冪、支え [a, F³a]、aligned は F²s_U ≥ a_U、t_U の決定、階級ごとの同時併合、原子に Alg 4 を当てて置換 | `src/paper/outer/scales.cpp: classify_outer_events, build_outer_hierarchy` | `tests/unit/paper/outer/test_outer_scales.cpp` | A, B | TODO | D-13(候補) |
| Sec 8.2, Lemma 8.4(l.2105) | 意味のある局所事象と有限な結合を O(m + Nτ(N)) で分類できる。閾値未満の原子は < b 頂点 | `scales.cpp` | `test_outer_scales.cpp` | A, B | TODO | none |
| Sec 8.3, Algorithm 8 l.2 零辺の縮約 | 零辺の成分を縮約し、ループを捨て、始点の成分だけを残し、元の辺の ID を保つ | `src/paper/outer/contract_restore.hpp: contract_zero_edges, ZeroContraction` | `tests/unit/paper/outer/test_contract_restore.cpp`, `tests/differential/paper/test_contract_restore_differential.cpp` | A, B | GREEN | none(解釈: 重複の除去は別の手順 Alg 8 l.5 なので、並行辺は残す。成分は union-find ではなくグラフ探索(BFS)で求める) |
| Sec 8.3, Lemma 8.5(l.2145) | 縮約後の距離から、元の全距離と最短路木を O(m+n) で復元する。有限の tight な弧を始点からたどり、最初に発見した弧で木を作る | `contract_restore.hpp: restore_predecessor_tree` | `test_contract_restore.cpp`(木の辺が tight、到達頂点を張る、零閉路), `test_contract_restore_differential.cpp` | A, B | GREEN | none(解釈: Evaluate は距離だけを出すので、復元は距離だけを受け取り tight な弧をたどる。正確でない入力の距離(丸めた double など、Q-06)は黙って受け入れず `std::invalid_argument` を投げる。d(u)+w の整数の桁あふれは tight でないとみなす) |
| Sec 8.3, Algorithm 8 l.4–7 基本の場合 | N ≤ L_N^64 なら、頂点対ごとに最軽量辺を残し、走査型 Dijkstra(Lemma 3.2)で距離を求め、Lemma 8.5 で復元する | `src/paper/outer/base_case.cpp: outer_base_case`(走査版)、`outer_base_case_d1`(二分ヒープ版) | `tests/unit/paper/outer/test_outer_base_case.cpp` | A | TODO | D-1(`outer_base_case_d1` だけ) |
| Sec 8.3, Algorithm 8 l.8–17 外側の構成 | MST と owner セル → NameRoutes → 分類と G' → OwnerDistances → radial 辺の削除 → 外側の状態と cap 付き局所階層 → 始点の種付け → Evaluate → 復元 | `src/paper/outer/shortest_paths.cpp: outer_construction` | Mode A では N < 2^64 で到達しない。`tests/differential/paper_experimental/test_outer_forced_vs_dijkstra.cpp`(Mode B で外側を強制) | A(到達しない), B | TODO | Mode B は D-7 |
| Sec 8.3, Algorithm 8(l.2167) ShortestPaths の入口 | 縮約 → 閾値の判定をその場で飽和演算で行い分岐 → 結果と FallbackRecord を返す | `shortest_paths.cpp: paper_shortest_paths` | `tests/differential/paper/test_paper_sssp_vs_oracle.cpp` | A | TODO | D-1(変種 `paper-faithful-d1`) |

## 2. 主要な定数・閾値

種別: 正 = 変えると距離の正しさが壊れうる、時 = 計算量と停止だけに効く(`thresholds.md` §1, §6)。Mode B で値を変えたら、その値を D-7〜D-10 に書く。

| Paper item | Statement (paraphrase) | Module / symbol | Tests | Mode | Status | Deviations |
|---|---|---|---|---|---|---|
| Sec 2 l.286、Alg 8 l.4: 外側 N > L_N^64(時) | 最小の N は 589^64+1。満たす集合は上に閉じず、[2^589−1, 590^64] で再び偽。その場で飽和演算で判定する | `src/paper/model/params.hpp: outer_nontrivial, kOuterExponent = 64` | `test_params.cpp`(N < 2^64 で常に偽)。予定していた多倍長の参照値(589^64, 589^64+1, 2^589−1, 590^64, 590^64+1)との照合は main に無い | A | GREEN | Mode B は D-7 |
| Sec 2 l.286、Alg 4 l.2: 局所 k > 2L_k^4(時) | 最小の k は 209,953(L = 18)。上に閉じている | `params.hpp: local_nontrivial, kLocalExponent = 4, kLocalCoeff = 2` | `test_params.cpp`(209,952 で偽、209,953 で真)。予定していた `test_local_hierarchy.cpp` での確認(209,952 は基本の場合、209,953 は radial)は main に無い | A | GREEN | Mode B は D-7 |
| Sec 5.1 l.966、Sec 7.1 l.1644: 分離の大きさ b(構) | 局所 b = L_N^4、外側 b = L_N^64、h = b+1 | `owner_decomposition(G, mst, b)` の引数 | `test_owner_decomposition.cpp`(b ∈ {1,2,3,N+1}) | A, B | TODO | Mode B は D-7 |
| Sec 5.2: 局所 K = 8h 以上の最小の 2 冪、near ≤ K²、macro w ≤ N³σ(正) | Lemma 5.5 に K ≥ 8h と N³ が要る | `RadialParams` | `test_radial_states.cpp` | A, B | TODO | Mode B でも変更不可 |
| Sec 5.4 l.1206: c = 10、N ≥ 32、h ≤ N、K < 16h、支え N^16σ(正) | aligned 判定と支えの範囲。N > 2L^4 から導かれる | `assemble.cpp` | `test_local_assemble.cpp` | A | TODO | Mode B は D-8(Q-03) |
| Sec 6, Alg 5 l.2–3: u = ⌈log2(t+2)⌉、g = 4u+3、P = 8·256·g·u² 以上の最小の 2 冪(時) | 合格判定 P·a* ≥ Z と d_e = Pλ_e。F は P から導く | `name_routes.hpp: name_routes_constants` | `test_name_routes.cpp` | A, B | GREEN | Mode B は D-17(`p_log2`) |
| Alg 5 l.5–6, l.13: 32t、16t、半数以上で受理(時) | 受理確率 ≥ 1/2 の根拠 | `name_routes.hpp: NameRoutesParams` | `test_name_routes_property.cpp` | A, B | GREEN | Mode B は D-17 |
| Sec 6.2 l.1473, l.1490: T = 4t 以上の最小の 2 冪、尺度の上端 2^ℓ ≥ tT²、窓 [a/T², T²a](正) | Lemma 6.4 に T ≥ t、Z ≥ d_S に最上位尺度が要る | `estimate.cpp` | `test_estimate.cpp` | A, B | GREEN | Mode B でも変更不可 |
| Sec 6.2 l.1491: H = (16t)^5 以上の最小の 2 冪、j ∈ [0,H)(時) | E[Z] の評価だけに効く。多語の整数になりうる(Q-43) | `estimate.hpp: EstimateParams::h_log2` | `test_estimate.cpp` | A | GREEN | Mode B は D-9 |
| Sec 8.1 l.1977: 外側 K = 16、near ≤ 256(正) | Lemma 8.1 の振れ幅 16 に必須 | `outer_states.cpp` | `test_outer_states_property.cpp` | A, B | TODO | Mode B でも変更不可 |
| Sec 7.3: 格子は原点の 7 オクターブ下〜最大 R の 2 オクターブ上、32q ≤ r < 64q、重なり [r/4, 4R]、片またぎ r_v > 16R_u、原子 w' < 2max(q)(正) | 時刻の順序(Lemma 7.6, 7.8)と単位 ≤ q(Lemma 7.7)の前提 | `calendars.cpp` | `test_calendars.cpp` | A, B | TODO | Mode B でも変更不可 |
| Sec 8.2 l.2080: F = 2P(2h+1)(8K²h)² より大きい最小の 2 冪、支え [a, F³a](正) | 支えが意味のある重みを覆う。F³ は double の範囲を超えうる(Q-06) | `scales.cpp` | `test_outer_scales.cpp` | A, B | TODO | D-13(候補)。Mode B でも式は変更不可 |
| Sec 4.1: 正規化の連結 4S(C) ≥ h(P)、ν ∈ (h/2, h](正) | 単位の整列の二者択一 | `normalize.hpp` | `test_normalize.cpp` | A, B | GREEN | Mode B でも変更不可 |
| Sec 3: B1 = max(2,⌈log2 L_k⌉)、B2 = L_k、LCA の B = max(2,⌊log2(k+1)/8⌋)、split–findmin の B = max(2,⌈log2(k+1)⌉)、Fibonacci の次数配列 3 + 2⌈log2(s+1)⌉(時) | ブロック長と上限。次数配列は小さすぎると範囲外アクセス | `mst.hpp`, `forest_queries.hpp`, `split_findmin.hpp`, `fib_heap.hpp` | 各モジュールのテスト | A, B | GREEN(fib_heap)/ VERIFIED(mst, forest_queries, split_findmin) | none |

## 3. 逸脱の一覧

形式は SKILL §3: `D-<n>: 何を — なぜ — 計算量と正しさへの影響`。「状態」は 採用(ADR あり)/ 候補(使うなら採用の手続きが要る)/ Mode B(Mode A では使わない)/ テスト専用。

| ID | 内容 | 状態 | 関係する問い |
|---|---|---|---|
| D-1 | Alg 8 l.5 の走査型 Dijkstra を、Lemma 3.2 の二分ヒープ版に替える(`paper-faithful-d1`、`outer_base_case_d1`) — 走査版は N = 10^6 で約 10^12 回の操作になり、計画した規模で測れない — 基本の場合の費用が O(N²+m) から O((N+m) log N) に変わる。距離は同じ。既定の `paper-faithful` は走査版のまま | 採用(ADR 0008) | Q-01 |
| D-2 | ⌊log2 i⌋ を Sec 2 の表ではなく `std::bit_width(i) − 1`(clz)で求める — 速く単純 — 論文のワード演算に clz が明記されていない。結果は同じ | 候補(ADR が要る) | Q-05 |
| D-3 | 上位の節(Sec 4, 7)を、Lemma 3.11 の 2 段 split–findmin の代わりに `SplitFindminRef`(セグメント木、操作 O(log k))で組む — 本実装より先に上位を組むため — Lemma 4.7 と 7.9 の計算量が O(log N) 倍になる。正しさは同じ | 候補(暫定。本実装に差し替えたら外す) | Q-30 |
| D-4 | Lemma 4.5 の巡回の開始位置を、(親, 子, preorder) の計数整列で作る O(1) の表(l.773–777)ではなく、二分探索で引く — 実装が単純 — Lemma 4.5 の O(M+N) に log 倍がかかる。正しさは同じ | 候補(使っていない: `child_order` は O(1) の表を作る) | Q-30 |
| D-5 | LCA と LA を Lemma 3.7–3.8 の O(1) 構造ではなく `BinaryLiftingRef`(前処理 O(k log k)、query O(log k))で引く。Sec 5.4 の「LCA の子」だけならオフラインの二分探索でもよい — 上位を先に組むため — 計算量が log 倍。正しさは同じ | 候補(暫定。Lemma 3.7–3.9 を実装した `forest_queries` と `child_order` は使っていない) | Q-30 |
| D-6 | GrowthPass の Fibonacci ヒープを添字つき二分ヒープ(または 4 分ヒープ)で置き換える — 実装と性能の都合 — decrease が O(log B) になり、パスの費用が O(l log B) に増える。重複要素のある遅延削除ヒープは Alg 1 l.9 の大きさの判定を変えるので不可 | 候補(使っていない: `FibHeap` が Lemma 3.3 を実装) | Q-17 |
| D-7 | Mode B: 閾値の指数・係数(外側 64、局所 4 と 2)と分離の大きさ b(外側 L^64、局所 L^4)を引数にする。h = b+1、K、P、F は変えた b から導き直す — 外側と局所の構成を実際の大きさで動かすため — 正しさは保たれる(`thresholds.md` §6)。線形性の主張は外れる。使った値は Mode B の結果ごとに記録する | Mode B | Q-03 |
| D-8 | Mode B: Sec 5.4 の N^c(c = 10)を Λ = 2·N·3hK² に、支えの上端 N^16σ を 2(3hK²)²N³σ に置き換える(02 §4.3) — 局所の閾値を下げると h ≤ N が崩れ、論文の不等式が成り立たない — 置き換えた式なら N ≥ 2、b ≥ 1 で証明が通る(02 の主張。検証が要る) | Mode B(要合意) | Q-03 |
| D-9 | Mode B: H を (16t)^5 以上の最小の 2 冪から 4tT² 以上の最小の 2 冪に下げる — 証明に必要なのは H ≥ 4tT² だけで、論文の値は 64 ビットに収まらない — E[Z] の定数が変わるだけで、Z ≥ d_S は決定的に保たれる | Mode B | Q-42 |
| D-10 | Mode B: NameRoutes の棄却が r 回(例: 64)続いたら、残りの辺に自分の重みを証明書として与えて終える — P を下げたときの停止を保証する — 論文に無い手順。証明書は有効なまま、辞書が大きくなる | Mode B(要合意) | Q-45 |
| D-11 | Mode B: Sec 5.4 の格子を使わず、結合レコードを重みで整列して等重みを同時併合する「正確な高さの組み立て」を正しさの基準に使う — 格子版の差分テストの参照 — O(N log N) の整列が入る。Def 4.1(c0 = 1)を満たす | Mode B(テストの参照) | Q-39 |
| D-12 | 入力の −0.0 を +0.0 に正規化する — ビット列が違う値を同じ長さとして扱うため — 論文に無い前処理。距離は変わらない | 候補 | Q-07 |
| D-13 | 実数の多項式倍(N³σ、N^cσ、N^16σ、F³a、3hK²)を、飽和演算か「実際の最大重み × 2」で打ち切って表す — uint64 と double で桁あふれする — 分類と見張りにしか使わないので、打ち切りの上限がすべての実在の重みを超えていれば判定は同じ(要証明) | 候補(要合意) | Q-06, Q-38 |
| D-14 | 評価器の `EvalOptions::supplemental_offers = false`(補足の申し出を切る) — Lemma 4.6 の正しさが補足の申し出に依存しないことを差分で確かめる — 計算量だけが悪くなる。報告する計測では使わない | テスト専用 | Q-49 |
| D-15 | Lemma 4.4/4.5 の間隔 g(i,j) を、符号付きの対の 2 成分を 1 回引いて 1 つの実数にする(Sec 2 l.292–296 は符号付きの値を対で持ち、実数の引き算を持たない) — 評価器が z + S(i) + g(i,j)(Sec 4.3 l.826)を普通の長さとして扱えるようにする — 整数の重みでは正確。double では間隔ごとに丸めが 1 回増え、根からの距離の丸めと合わせて、g が証明された値を全長の O(N²ε) 倍だけ下回りうる(g ≥ 0 は保つ)。そのため補足の申し出が数 ulp だけ実現不能になりうる | 候補 | Q-06 |
| D-16 | `LcaTable`、`LevelAncestor`、`ForestQueries` のテスト専用のコンストラクタが、論文の B = max(2,⌊log2(k+1)/8⌋) の代わりに [2, 8] のブロック長 B を明示的に受け取る — 論文の B は k < 2^24−1 ではすべて 2 で、ほかの大きさのブロック表が試されない — Mode A の呼び出し側が使う既定の経路は論文の B のまま。答えはどの B でも同じ | テスト専用 | — |
| D-17 | Mode B: `NameRoutesParams::sample_factor`、`residual_factor`、`p_log2` を、論文の 16t、32t(Alg 5 l.5–6)と P(Alg 5 l.3)から変えてよい — 小さい t で標本の回を何度も回すため(テストでは標本 1t、残りの閾値 2t、打ち切りのテストでは P = 1) — 証明書の正しさはこれらの値によらない(Lemma 6.8(1))。変わるのは受理確率と時間だけで、Lemma 6.8(2) は成り立たなくなる | Mode B | Q-03, Q-45 |
| D-18 | ブロック化した split–findmin(Lemma 3.11)で、レコードのブロック番号を、要素ごとに持つブロック番号とオフセットではなく、整数の除算 `i / B` で求める(Sec 2 l.305–306 はこれらを 1 回の走査で作り、整数の除算は要らないと書く。l.302–303 が定数時間とする語演算に除算は入っていない) — コードが単純 — ブロックも結果も同じ。対象のハードウェアでは除算は 1 命令なので計測には影響しないが、論文の語モデルでは列挙された演算の外になる | 候補 | Q-51 |
