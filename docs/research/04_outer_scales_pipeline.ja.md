# 04. 外側の尺度と外側の流れ（論文 Sec 8.2〜8.3）：実装のための仕様

[English](04_outer_scales_pipeline.md) | **日本語** | [简体中文](04_outer_scales_pipeline.zh-CN.md)

対象: arXiv 2611.34911v1 の Sec 8.2（Lemma 8.4、l.2076〜2142）と Sec 8.3（Lemma 8.5 l.2145、Algorithm 8 l.2167〜2196、Theorem 1 の証明 l.2158〜2217）。行番号は `.cache/paper/2611.34911v1.txt` のもの。
対応する課題: **P3-OUTER-SCALES**（`src/paper/outer/scales`）と **P3-OUTER-PIPELINE**（`src/paper/outer/shortest_paths` の `outer_construction`）。
この文書は `03_pivots_aux_outer.md` の §3〜§4 の続きで、実装の担当がそのまま使える仕様にしたもの。

記法: `【論文】` は論文が述べていること。`【導出】` は論文の記述から短い議論で出ること（議論はここに書く）。`【解釈】` は曖昧な箇所の読み方（§5 に Q の案がある）。`【補完】` は論文に無い補い。`【実測】` はこの文書のために計算または実行したもの（多倍長のスクリプト、または `main` にある部品をつないだ作業用のプログラム。どちらもコミットしていない）。

この文書の範囲の限界:
- P3-OUTER-STATES（Sec 8.1）の出力は、まだ `main` に無い。尺度の段がそこから必要とするものを、§2.0 に**前提 A-1** として並べた。その課題が入ったら、名前を合わせること。
- 端から端までは、何も実行していない。外側の構成がまだ無いためである。実行したのは、半径状態までの段だけ（§3.5）。

---

## 0. 要点

1. `F = 2^f`、`f = p + 23 + bit_width((2h+1)·h²)`。ここで `P = 2^p`、`h = b + 1`、`K = 16`。持つのは指数 `f` だけにする。`F²s` と `F³a` は、`2f` 回と `3f` 回の飽和する倍加で作る（D-13）。`h` は**正確な** `b + 1` でなければならず、飽和した語の値を使ってはいけない（§1.3）。
2. `uint64` の重みで論文の `P` を使う（すると `f ≥ 45`）と、`F²·s` はどの `s ≥ 1` でも飽和する。つまり、**付着のあるブロックはすべて aligned になり、support はすべて1つの格子に併合される**。`uint64` で unaligned の枝に届くのは `f ≤ 31`（Mode B、`P ≤ 8`、`b ≤ 3`）のときだけ。それ以外では、`2f` オクターブを超える幅を持つ `double` でしか届かない（§1.5）。
3. Sec 8.2 が使い回す Sec 5.4 の格子・併合・リンクの仕組みは、`src/paper/local/assemble.cpp` の無名名前空間の中にある。しかも小さな成分には正確な Kruskal 木を作るが、Sec 8.2 が要るのは**分割階層（Def 5.12）の atom の葉**である。そのままでは呼べない（§4、M-6、M-7）。
4. Mode A は `N < 2^64` では Alg 8 l.8〜17 に届かない（Q-04）。同じコードを Mode B のパラメータで動かして確かめる。外側を強制する最小の集合は、外側の閾値 `0·L^0` と定数の分離子サイズ `b`。ほかはすべて論文の値のまま（§3.4）。
5. 外側の枝が自明でなくなる最小の（縮約後の）`N`【実測、§3.5】: `N = 2`（ピボット2個、`b = 1`）、`N = 4`（付着のある、状態2個以上のブロック、`b = 2`）、**`N = 6`（状態2個以上の atom。Sec 8.2 のすべての段が実際に働く。`b = 3`）**。
6. いま宣言されている `outer_construction(c, diag)` には、乱数生成器もパラメータも無い。`PaperOptions` には seed が無い。Alg 8 l.9 を呼ぶ前に、どちらも要る（§4、M-3、Q-70）。

---

## 1. 係数 F（Sec 8.2 l.2080〜2081）

### 1.1 定義

【論文】`F` は、`X = 2·P·(2h+1)·(8K²h)²` より真に大きい最小の2冪。整数のパラメータから、グラフごとに1回だけ決める。`P` は NameRoutes の共通係数（Alg 5 l.3）、`h = b + 1` は外側の owner 分解のホップ上界（Sec 7.1 l.1644）、`K = 16`（Sec 8.1）。

【導出】各因子の由来: `2P(2h+1)` は、同じ名前を持つ2本のマクロ辺の重みの比の上界（Lemma 7.1: `λ/(2h+1) ≤ w < 2Pλ`）。`4K²h` は Lemma 8.2 の比 `θ/B0`。新しい名前 `a`（§2.1）を持つ有限の join では `a ≤ B0 ≤ θ < 8P(2h+1)K²h·a = X·a/(8K²h)` なので、`θ < F·a` が `8K²h` 倍の余裕つきで成り立つ。

### 1.2 厳密な（値, 指数）の形（Q-06）

`F` は（値 1、指数 `f`）の組。`F²` は（1, `2f`）、`F³` は（1, `3f`）。`F` を長さや機械整数として実体化することはしない。

【導出】`2^f > X` となる最小の `f` は `bit_width(X)`。`K = 16`、`P = 2^p` のとき:

```
X = 2 · 2^p · (2h+1) · (2^11 · h)²  =  2^(p+23) · (2h+1) · h²
f = p + 23 + bit_width((2h+1) · h²)
```

【実測】この閉じた形は、§1.4 のすべての行で `bit_width(X)` に一致した。

### 1.3 あふれない f の計算

`f` は大きさについての整数の語演算（Sec 2 l.302〜304）であって、長さの演算ではない。

```cpp
// 入力: p_log2（unsigned）、h は正確な自然数（多倍長）
// Y = (2h + 1) * h * h      -- 多倍長の演算3回。Y < 2^(3·bits(h) + 2)
// f = p_log2 + 23 + bit_width(Y)
```

- `h` は正確でなければならない。Mode A では `b = L_N^64` はどの `N` でも1語に収まらない（`L_N ≥ 2`）。`owner::outer_separator_size` と `OwnerDecomposition::hop_bound` が持っているのは、**飽和した**値 `2^64 − 1` である。飽和した `h` から `f` を出すと誤る。`N = 10^6` では、871 ではなく `f = 233` になる【実測】。`L^64` と `h` は多倍長で作る（語と多倍長の乗算を64回。`L ≤ 65` なので `h < 2^386`）。
- `src/experimental/paper_experimental/params.cpp` は、すでにこれをしている（`OuterDerived::f_log2`、無名名前空間の `Nat` クラス）。`src/paper/` は `src/experimental/` を include できない。M-2 を参照。
- `f`、`2f`、`3f` は、表現できるどの入力でも `unsigned` に収まる（`f ≤ 29 + 23 + 3·386 + 2`）。

長さへの `F` の適用（D-13、Q-38、`src/paper/local/saturating.hpp`）:

```cpp
// aligned(s, a)  :=  !(sat::mul_pow2(s, 2 * f) < a)        // F² s ≥ a、Sec 8.2 l.2086
// support_end(a) :=  sat::mul_pow2(a, 3 * f)               // F³ a、   Sec 8.2 l.2081
```

【導出】どちらも、判定としては厳密である。`sat::mul_pow2` が `WeightTraits<W>::infinity()` を返すのは、真の積が型に収まらないときに限る。そのとき真の積は `infinity()` 以上で、その型のどの有限の長さよりも大きい。したがって `F²s ≥ a` の真偽は、どの有限の `a` についても、正確な算術と同じになる。飽和した support の上端が後ろの support と重なるのは、真の上端が重なるときと一致する（後ろの support は有限の値から始まる）。`sat::mul_pow2` は無限大で止まるので、ループは整数で高々 `min(2f, 語のビット数)` 回、`double` で高々 2098 回。
ここで `LenOps::mul_pow2` を使ってはいけない。あふれると、飽和せずに例外を投げる。

### 1.4 指数の表【実測、Python の多倍長整数】

`u = ⌈log2(t+2)⌉`、`g = 4u+3`、`P = 2^p` は `2048·g·u²` 以上の最小の2冪（Alg 5 l.2〜3）。

論文の大きさと Mode A の行（`b = L_N^64`）:

| N | b | t | p | f | 2f | 3f |
|---|---|---|---|---|---|---|
| 589^64 + 1（外側が動く最小の N） | 589^64 ≈ 2^588.94 | 1（これしかありえない） | 17 | **1808** | 3616 | 5424 |
| 589^64 + 1（式だけ） | 589^64 | 66 | 22 | 1813 | 3626 | 5439 |
| 589^64 + 1（式だけ） | 589^64 | 2^20 | 27 | 1818 | 3636 | 5454 |
| 590^64 + 1 | 590^64 ≈ 2^589.09 | 1 | 17 | 1809 | 3618 | 5427 |
| 10^3（Mode A では届かない） | 10^64 | 1 | 17 | 679 | 1358 | 2037 |
| 10^6（Mode A では届かない） | 20^64 | 1 | 17 | 871 | 1742 | 2613 |
| 2^63（Mode A では届かない） | 64^64 = 2^384 | 1 | 17 | 1194 | 2388 | 3582 |

Mode B の行（`t` の列は、`thresholds.md` §5.2 と同じ上界 `2⌊N/b⌋ − 1`）:

| N | b | t | p | f | 2f | 3f |
|---|---|---|---|---|---|---|
| 10^6 | 20 | 99999 | 26 | 64 | 128 | 192 |
| 10^6 | 256 | 7811 | 25 | 74 | 148 | 222 |
| 10^6 | 400 | 4999 | 25 | 75 | 150 | 225 |
| 10^6 | 4096 | 487 | 23 | 84 | 168 | 252 |
| 10^6 | 8000 | 249 | 23 | 86 | 172 | 258 |
| 10^6 | 160000 | 11 | 20 | 96 | 192 | 288 |
| 10^8 | 27 | 7407405 | 27 | 66 | 132 | 198 |
| 10^8 | 729 | 274347 | 26 | 79 | 158 | 237 |
| 10^8 | 4096 | 48827 | 26 | 87 | 174 | 261 |
| 10^8 | 19683 | 10159 | 25 | 92 | 184 | 276 |
| 10^8 | 531441 | 375 | 23 | 105 | 210 | 315 |

小さな分離子サイズ（外側を強制するテスト用）。論文の `P` での `f`（`t = 1〜2` は `p = 17`、`t = 6` は `p = 19`、`t = 66` は `p = 22`）と、D-17 の値 `P = 1`、`P = 2` での `f`:

| b | h | bit_width((2h+1)h²) | f (t ≤ 2) | f (t = 6) | f (t = 66) | f (P = 1) | f (P = 2) |
|---|---|---|---|---|---|---|---|
| 1 | 2 | 5 | 45 | 47 | 50 | 28 | 29 |
| 2 | 3 | 6 | 46 | 48 | 51 | 29 | 30 |
| 3 | 4 | 8 | 48 | 50 | 53 | 31 | 32 |
| 4 | 5 | 9 | 49 | 51 | 54 | 32 | 33 |
| 8 | 9 | 11 | 51 | 53 | 56 | 34 | 35 |
| 16 | 17 | 14 | 54 | 56 | 59 | 37 | 38 |
| 33 | 34 | 17 | 57 | 59 | 62 | 40 | 41 |
| 64 | 65 | 20 | 60 | 62 | 65 | 43 | 44 |

これらは、`tests/unit/paper_experimental/test_experimental_params.cpp` ですでに確かめている値（871。`b = 160000`、`t = 11`、`P = 2` で 77）と、`thresholds.md` §5.1〜5.2 に一致する。

### 1.5 重みの型ごとの、飽和の意味【導出】

| 型 | `F²·s ≥ a` | support `[l, F³a]` |
|---|---|---|
| `uint64`、`2f ≥ 64`（論文の `P` では常に） | どの `s ≥ 1` でも真。有限の付着を持つブロックはすべて **aligned** | どのパラメータでも `3f ≥ 84 > 64`。support の上端はすべて飽和するので、**すべての support が1つの格子に併合される** |
| `uint64`、`2f ≤ 63`、つまり `f ≤ 31` | 偽になりうる。`p + bit_width((2h+1)h²) ≤ 8` が要る。つまり `b = 1, p ≤ 3`、`b = 2, p ≤ 2`、`b = 3, p = 0` | 上と同じ |
| `double`、`f ≤ 1022` | `a > 2^(2f)·s` のとき偽（`2f` オクターブを超える幅を持つ2進有理数の重み） | support どうしが `3f` オクターブを超えて離れると、格子の成分が分かれる |
| `double`、論文の大きさ（`f ≥ 1808`） | `2^3616` は `double` の全範囲（最小の非正規化数と最大の有限値の間で 2^2098）を超える。常に aligned | 常に1つの格子 |

テストへの帰結: unaligned の枝と、成分が複数ある格子を確かめるには、`uint64` で Mode B の `p_log2 = 0` と `b ≤ 3` を使うか、2進有理数の重みの `double` を使う（§6）。

---

## 2. Lemma 8.4 の各段（Sec 8.2 l.2076〜2142）

### 2.0 入力

`main` にある部品から:

| 入力 | 型 | 出どころ |
|---|---|---|
| 縮約後のグラフ | `CsrGraph<W>` | `outer::ZeroContraction<W>::graph` |
| owner 分解（外側の `b`） | `OwnerDecomposition<W>`: `radius`（r_v）、`owner`、`pivots`、`corridors`、`hop_bound` | `owner_decomposition(g, b)` |
| 分類 | `owner::OuterEdgeClassification<W>`: `edge_class`、`name`、`dictionary`、`p_log2` | `owner::classify_interowner_edges` |
| 補助グラフ | `owner::AuxGraph<W>`: `aux_vertex`、`radius`、`owner_distance` | `owner::build_aux_graph` |
| calendar（別の格子のため） | `owner::Calendars<W>`: `graph_piece`、`radius_index`、`graph_pieces[].origin`、`.grid` | `owner::OwnerDistancesResult<W>::calendars` |

**前提 A-1**（P3-OUTER-STATES。`main` に無い。フィールド名は仮）。尺度の段が必要とするもの:

| 必要なもの | 意味 |
|---|---|
| `StateGraph<W> graph`、`projection` | すべての状態弧（通常の状態と特別な状態）、状態ごとの物理的な像 |
| 状態ごとの `block`、`num_blocks`、ブロックの種類 | 種類 ∈ {正の通常ブロック、ピボットのコピー、特別な状態}。特別な状態とピボットのコピーは、状態1個のブロック（Sec 8.1 l.1996） |
| `arc_kind`、`arc_edge` | 窓の中の弧と、その入力辺を見分けるため |
| 頂点ごとの `native_inward` | Alg 8 l.14 の種。ピボットでは内向きのピボットコピー |
| 状態ごとの `theta`、`b0`、`b0_edge` | Lemma 8.2 の値と、それを達成するマクロ辺（「有限の B0 は、達成するマクロへのポインタを持つ」l.2078〜2079） |
| join の記録 | ブロックの MST 辺（その重みで）、状態から owner への join（θ で）、corridor の join（corridor の重みで）（Lemma 8.3）。形は `ThresholdJoins<W>::joins` と同じ |

P3-OUTER-STATES が `build_radial_states` を使い回すなら、`EdgeClassification::macro_key[e] = cls.name[e]` を渡せる。すると `ThresholdJoins::b0_key` は、達成するマクロの名前の辞書の添字になり、§2.1 が必要とするものにちょうど一致する。
単位: Alg 7 は `scaled_distance = 2^7·J` も返す（D-21）。この節のものはすべて**入力の単位**（r、θ、MST の重み、名前）。

### 2.1 段 S0：新しい名前（l.2077〜2079）

【論文】各マクロ辺の名前を、その名前を持つマクロ辺の物理的な重みの最小値で置き換える。

```cpp
// 入力: g、cls（edge_class、name、dictionary）
// 出力: std::vector<W> name_weight   // 辞書の添字ごと。その名前のマクロ辺が無ければ infinity()
for (EdgeId e : cls.edge_class[e] == kMacro である辺)
  name_weight[cls.name[e]] = min(name_weight[cls.name[e]], g.edge(e).weight);
```

- 【解釈、Q-65】「物理的な重み」は入力の重み `w_e` であり（Lemma 7.1 は `w` について述べている）、射影した長さ `A_e` ではない。最小値はマクロ辺だけについて取る。
- すでにある関数: 無い。線形の走査で済む。長さは比較するだけ。
- 不変条件（テストにできる）: どのマクロ辺 `e` についても、`a = name_weight[name[e]]` として `a ≤ w_e`、`w_e < 2P(2h+1)·a`（飽和する積で）、`a` はその名前を持つあるマクロ辺の重み。θ が有限のどの状態 `z` についても、`a_z = name_weight[name[b0_edge[z]]]` として `a_z ≤ b0[z] ≤ theta[z]` かつ `theta[z] < F·a_z`。

### 2.2 段 S1：F

`f = outer_scale_factor_log2(p_log2, 正確な h)`（§1.3）。グラフごとに1回。【補完】マクロ辺が無い（`t = 1`）ときは、`f` を計算はするが、使うことは無い。

### 2.3 段 S2：ブロックの尺度と付着（Sec 8.1 l.1999、Sec 8.2 l.2083〜2086）

```cpp
// 入力: 状態（A-1）、od.radius、name_weight
// 出力: std::vector<OuterBlockScale<W>> blocks      // ブロックごと
struct OuterBlockScale {
  Kind kind;            // kSingleton、kAligned、kUnaligned、kUnattached
  std::size_t size;
  W s;                  // s_U = U の状態 x についての od.radius[projection[x]] の最小値（Q-35）
  Vertex s_vertex;      // s_U を達成する頂点（番号が最小のもの）。別の格子で使う
  StateId z;            // U の中で theta を最小にする状態（番号が最小のもの）。すべて無限なら kNoState
  std::uint32_t name;   // b0_edge[z] の名前の辞書の添字。a_U = name_weight[name]
  W t;                  // t_U。S5 で決める
};
```

- 【論文】状態1個のブロック（ピボットのコピー、特別な状態、状態1個の通常ブロック）は、自分の有限の θ の join だけを持つ。ブロックの尺度の判定、格子、MST の分類、atom の置き換えは飛ばす（l.2083〜2085）。
- 【論文】有限の付着を持つ、状態2個以上の正の通常ブロックでは、最小にする状態を選び、その名前を `a_U` とする（l.2085〜2086）。
- 【解釈、Q-35】`s_U` は元の半径 `r_v = b_T(v, p(v))` で取る。`J` ではない。
- 【補完】同点（θ を最小にする状態が複数、半径 `s_U` の頂点が複数）は、番号が最小のものにする。どれを選んでもよい。最小にするどの `z` も `a_z ≤ θ(z) < F·a_z` を満たす。
- 不変条件: `U` のどの状態 `x` でも `0 < s_U ≤ r_x`。どの `z ∈ U` でも `s_U ≤ θ(z)`（Lemma 8.2 の `r ≤ θ`）。`θ(z_U) = min_U θ`。1つのブロックの中では、物理的な頂点1個につきコピーは1個（ブロックの状態の射影は単射）。

### 2.4 段 S3：aligned か unaligned か（l.2086〜2087、Lemma 8.4 の証明 (1)、(2)）

```cpp
// 入力: blocks（kind はまだ分けていない）、name_weight、f
// 出力: 付着のある状態2個以上のブロックについて、blocks[U].kind ∈ {kAligned, kUnaligned}
// 付着のある状態2個以上のブロックを blocks[U].name でまとめる。名前の重みが a のまとまりで:
//   s_star = pivots::median_prune(keys = {s_U}, feasible(s) = !(sat::mul_pow2(s, 2f) < a))
//   U が aligned  ⇔  s_star があり、かつ !(s_U < *s_star)
```

- 【論文】`F²·s_U ≥ a_U` のとき aligned。論文は判定を名前ごとにまとめ、中央値による枝刈りを使う。調べる中央値1個につき `O(log F)` 回の倍加（l.2130〜2131）。
- すでにある関数: `pivots::median_prune`（`src/paper/pivots/name_routes.hpp`。述語は `s` について単調）、`sat::mul_pow2`。名前でまとめるのは、辞書の添字での計数ソート。
- 【補完】ブロックごとに判定しても分類は同じで、ブロック1個につき最大 `2f` 回の倍加になる。変わるのは仕事量の上界だけ。そうするなら、時間だけに効く逸脱として記録する。
- 不変条件:
  - unaligned のブロック: ブロックのどの MST 辺（実際には `U` の窓の中のどの弧）も重みが `< a_U/2`。`sat::add(w, w) < a_U` で確かめる。【論文、l.2086〜2087】理由: `w ≤ 4K²·J_v ≤ 4K²h·s_U < a_U/2`。
  - unaligned のブロック: `U` のすべての状態で θ が等しい（Lemma 8.4 の証明 (1)、「付着の閾値が等しい」）。
  - aligned のブロック: `s_U ≤ θ(z_U) < F·a_U` かつ `a_U ≤ F²·s_U`。
  - 差分: まとめて判定した結果が、ブロックごとの判定 `!(sat::mul_pow2(s_U, 2f) < a_U)` に一致する。
  - 比較の回数: 大きさ `n_a` の名前のまとまりごとに、述語の判定は高々 `⌊log2 n_a⌋ + 1` 回。

### 2.5 段 S4：support と大域の格子（l.2081〜2083、l.2087〜2091、Lemma 8.4 の証明 (2)）

キー: ある有限の join の名前になっている辞書の添字ごとに1つ、`T_A` の corridor（`od.corridors`）ごとに1つ。

```cpp
// 入力: name_weight、od.corridors、aligned のブロック、theta の join、f
// 出力: std::vector<W> grid;                       // すべての成分の格子点。増加順
//       std::vector<std::uint32_t> component_of;   // 格子点ごと
//       std::vector<GridSupport<W>> support;       // キーごと: lower、upper、lower_index、upper_index
// 1. lower[key] = a（名前の重み、または corridor の重み）。upper[key] = sat::mul_pow2(a, 3f)
// 2. aligned のどの U についても: lower[name_U] = min(lower[name_U], s_U)
// 3. support を lower で整列する。隣どうしが重なるとき、または
//    next.lower < 2 * （いまの成分の upper）のとき（比が 2 未満の隙間）に併合する
// 4. 併合した成分ごとに: 格子点 lo, 2lo, 4lo, ... を、成分の上端より真に大きい最初の点まで作る
//    （最後の倍加1回は番兵）。infinity() で止める
// 5. キーごとに: lower_index = lower 以下で最大の格子点。upper_index = upper より大きい最初の格子点
```

- 【論文】support は `[a, F³a]`。aligned のブロックについて始点を `s_U` まで下げる。重なりと、2倍未満の隙間を併合する。添字は `l` の格子上の切り下げから、`u` より真に大きい最初の点まで。名前が `k = O(t log²(N+2))` 個のとき、オクターブは `O(k log F)` 個。
- 【解釈、Q-65】corridor の join は、マクロとしての名前ではなく、自分の重みの support を使う。
- 【解釈、Q-68】「2倍より小さい隙間」は、厳密な不等号の判定 `next.lower < 2·upper` とする。この併合のあとは、どの類も成分の始点をまたげない。したがってどの類もちょうど `[g, 2g)` になる（Sec 5.4 では類が次の成分で止まりうる、l.1279。それとは違う）。
- 【補完、`assemble.cpp` と同じ】有限の join の重みの最大値より上の格子点は、作らなくてよい。そこに分類されるものは無い。
- すでにある関数: 呼べるものは無い。同じ手順が `src/paper/local/assemble.cpp`（無名名前空間）にあり、`F³a` の代わりに `N^16σ` を使っている。M-6 を参照。
- 不変条件:
  - 有限のどの `θ(z)` も、自分のキーの `[lower, upper]` に入る。どの corridor の重みも、自分の support に入る。
  - キー `κ` を持つ aligned の `U` について: `lower[κ] ≤ s_U`。`U` のブロック MST 辺で `w ≥ t_U` のものはすべて `grid[lower_index[κ]] ≤ w < grid[upper_index[κ]]` を満たす（Lemma 8.4 の証明 (1): `w ≤ 4K²h·s_U < F²·a_U`）。
  - aligned の `U` と、θ が有限のどの `z ∈ U` についても: `name(z)` の support と `name_U` の support は、同じ併合成分に入る。【導出】`s_U ≤ θ(z) ≤ max(θ(z_U), 4K²h·s_U) < F²·a_U` なので、`θ(z)` は両方の support に入る。
  - 格子: 真に増加。成分の中では `grid[i+1] = 2·grid[i]`。成分の間では `grid[次の成分の最初] ≥ 2·（前の成分の上端）`。
  - 大きさ: 併合の前で、入力の support 1個につき高々 `3f + 2 + log2(下げ幅)` 点。下げ幅は高々 `F²`（aligned では `s_U ≥ a_U/F²`）なので、高々 `5f + 2` 点。

### 2.6 段 S5：切り捨て値 t_U（l.2091〜2093、l.2100〜2104）

```cpp
// aligned の U（キー k）:  t_U = s_U 以下で最大の格子点。[lower_index[k], upper_index[k]] の中で探す
// それ以外の、状態2個以上の正の通常ブロック（unaligned、または有限の付着が無い）:
//   v  = blocks[U].s_vertex;  x = aux.aux_vertex[v];  C = cal.graph_piece[x]
//   i  = cal.radius_index[x]              // cal の格子は 2^7 倍: grid[i] <= 2^7 r_v < 2 grid[i]
//   t_U = sat::mul_pow2(cal.graph_pieces[C].origin, i - 7)     // origin * 2^(i-7)。入力の単位
```

- 【論文】どちらの選び方でも `s_U/2 < t_U ≤ s_U`（Q-36: 左は厳密な不等号。Sec 5.4 とは違う）。
- 【論文】別の格子について: `s_U` を達成する `v` を選ぶ。ブロックの中のどの `J` も `K²·J_v` 以下。局所の辺の最大は `4K²·R_v` 以下で、意味のある辺の最小は `r_v/2` より大きい。`v` の補助格子へのポインタを使い、各 piece の表を、決まった数のオクターブだけ、1回だけ伸ばす。
- 【解釈、Q-67】「補助格子へのポインタ」は、Sec 7.3 の、`v` の graph piece の格子（`Calendars`）とする。その格子は `2^7` 倍して持っており、`grid[7] = 2^7·origin` である。したがって、入力の単位で `r_v` の下にある格子点は `origin·2^(i−7)` で、`i ≥ 7`（`r_v ≥ origin` だから）。倍加だけで届き、半分にする演算は要らない。
- 【解釈、Q-66】「それ以外の、状態2個以上の正の通常ブロック」には、有限の付着が無いブロックも入る（Sec 5.4 l.1213 と同じ）。
- 【補完】`U` の別の格子の類は `[t_U·2^j, t_U·2^(j+1))`、`j ≥ 0`。ブロックの MST 辺の最大（`≤ 4K²·R_v`、つまり `j ≤ 10 + ⌈log2 h⌉`）の類まで。これらは piece の格子点である。いまの piece の格子は `4·max R` までしか保証されていない。論文が述べるように1つの表を共有するなら、上に倍加9回ぶん（`K²` のために8回、番兵に1回）伸ばす必要がある。ブロックごとに `t_U` から倍加で点を作っても類は同じで、ブロック1個につき `O(log(K²h))` 回の倍加になる（時間だけに効く。M-8）。
- 不変条件: `t_U ≤ s_U < 2·t_U`（`!(s_U < t_U) && s_U < sat::add(t_U, t_U)` で確かめる）。aligned: `t_U` は `name_U` の成分の格子点で、`grid[lower_index] ≤ t_U`。別の格子: `t_U = origin·2^j`（`s_vertex` の piece の origin）。

### 2.7 段 S6：atom（l.2000〜2001、l.2133〜2142、Def 5.11）

```cpp
// 入力: ブロック MST の join の記録、blocks[U].t
// 出力: 状態ごとの部分の番号。std::vector<ReplaceableAtom<W>> atoms
// 状態の上の DisjointSets。重みが < t_U のブロック MST 辺の両端を併合する。
// 状態1個の成分は、単独の部分。状態2個以上の成分 Y は atom:
//   states        = Y の状態
//   edges (H)     = 両端が Y にある、窓の中の弧の組。入力辺1本につき1回。
//                   物理的な重みで（重い弦も入れる）
//   spanning_tree = Y の中の、重みが < t_U のブロック MST 辺
//   cutoff        = t_U
```

- 【論文】`t_U` より軽い辺の成分を、局所の表現に置き換える。各 atom は `T − A` の1つの成分に入り、頂点は `b` 個未満。atom の中では、owner と元の半径は一定（Lemma 8.4 の証明 (3)）。
- 【解釈、Q-53 の拡張、Q-74】内部グラフには、両端が atom にある、窓の中の辺をすべて入れる。`substitute_atoms` は、両端が1つの atom にある `D` の弧を拒むので、これらはすべて `H` に入れなければならない。
- 【導出】ほかの種類の弧は、両端が1つの atom に入ることが無い。transfer と switch は窓か相を出る。ピボットの弧と特別な状態の弧は、状態1個のブロックに触れる。マクロ辺は2つの owner セルをつなぐが、atom の owner は1つである。これを検査し、破れていれば `std::logic_error` を投げる。
- 【解釈、Q-74】状態1個の成分は、単独の部分とする（Def 5.12、Lemma 5.15 の「単独でないすべての部分」）。Alg 4 は呼ばない。
- すでにある関数: `DisjointSets`（`src/paper/primitives/dsu.hpp`）。型 `ReplaceableAtom<W>`、`kNoAtom`（`src/paper/local/atom_replacement.hpp`）。
- 不変条件: `|Y| < b`。`Y` の上で owner は1つ、`r` は1つ、`r ≥ s_U ≥ t_U`。`Y` は1つのブロックと、`T − A` の1つの成分に入る。`Y` の射影は単射。全域木のどの辺も `< t_U`（したがって `< 2·cutoff`）。`H` は連結。`z ∈ Y` の有限の `θ(z)` はすべて `≥ t_U`。

### 2.8 段 S7：事象の分類（l.2091、l.2094、Lemma 8.4 の証明 (2)）

| 事象 | 類 |
|---|---|
| ブロック MST 辺、`w < t_U` | 無し（atom の中） |
| ブロック MST 辺、`w ≥ t_U`、`U` が aligned | 大域の類 `i`: `grid[i] ≤ w < grid[i+1]`。`name_U` の添字の範囲で二分探索 |
| ブロック MST 辺、`w ≥ t_U`、`U` が unaligned または付着無し | `U` の別の格子の類 `j`: `t_U·2^j ≤ w < t_U·2^(j+1)` |
| 状態から owner への join（有限の θ(z)） | `name(z)` の添字の範囲での大域の類 |
| corridor の join | その corridor のキーの添字の範囲での大域の類 |

- 【論文】事象は `O(N)` 個。それぞれ `O(1 + log log F) = O(τ(N))` 回の比較。
- 不変条件: 見つけた類が重みを挟んでいる。事象1個あたりの二分探索の歩数は `≤ ⌈log2(upper_index − lower_index + 1)⌉`。unaligned の `U` について、`p_top` を別の格子の最上位の類の下端、`g` を `θ(z_U)` の大域の類とすると `2·p_top ≤ g`【導出: `2·p_top ≤ 2·w_max < 8K²h·a_U/F² < a_U/2 < g`】。つまり、最初の大域の事象より前に、ブロック全体が1つの節点になる（l.2093〜2094）。

### 2.9 段 S8：分割階層（l.2094〜2099、Def 5.12、Lemma 5.10 の後処理）

```cpp
// 入力: 類のついた事象、部分（単独と atom）、状態グラフ
// 出力: PartitionHierarchy<W> ph   // graph = 外側の弧だけ。nodes、leaf_atom、atoms、c0 = 2
// 1. 葉: 単独の部分ごとに1つ、atom ごとに1つ
// 2. join の頂点: 状態と、ピボットごとの補助頂点1つ。DisjointSets（Lemma 3.6）
// 3. 別の格子を持つブロックごとに、類を j の増加順に。そのあと、大域の類を添字の増加順に。
//    1つの類の事象はすべて同時に併合する。節点の高さは、類の下端。
//    状態を持たない補助節点は捨て、子が1つの節点は詰める。
// 4. リンク: 外側の弧はそれぞれ、両端の部分の LCA で、それがつなぐ子の組に割り当てる。
//    類 [g, 2g) の節点では、割り当てた弧のうち重みが < 2g のものが子の隣接関係を作る。
//    DFS で、発見辺1本につき弧を1本取る（端点は実際の状態）。
// 5. 順序: 各節点で、割り当てた弧のうち重みが < g のものによる、子の位相順。
//    根は、根の間の弧による凝縮の順。
```

- 【論文】整数の格子の類ごとに併合を処理し、補助節点を除き、各節点より軽い弧を位相順に並べると、分割階層ができる。連結木は Sec 5.4 の子の商の探索で取り出す。類 `[a, 2a)` については、Lemma 8.3 の厳密版が、`2a` より真に軽い弧による連結性を与える。
- 【解釈、Q-68】高さは類の下端、リンクは上端より軽い弧。したがって `c0 = 2`。Lemma 5.10 が Sec 5.4 についてしているのと同じ。「別の格子が先、そのあと大域の類」という順序は、こちらで決めたもの。各根への道に沿った高さと矛盾しない順序なら、どれも同じ結果になる。
- すでにある関数: `DisjointSets`、`ForestQueries`（LCA と level ancestor）、`add_internal_node`、`collect_roots`。型 `PartitionHierarchy<W>`。`assemble.cpp` の併合の組み立てと、リンク・順序の処理は呼べない（M-6、M-7）。
- 不変条件:
  - `substitute_atoms` が `ph` を受け入れる（部分が状態を分割している。`D` の弧で atom の中にあるものは無い。どの atom の cutoff も親の高さ以下）。
  - 高さ `g`、類 `[g, g')` のどの内部節点 `X` についても: `X` の下の状態は、長さ `< g'` の状態弧の1つの SCC。`X` の子は、その中での、長さ `< g` の弧の SCC（Lemma 8.3 の厳密版。小さな入力で Tarjan と照合する）。
  - 根に向かって高さは減らない。どのリンクも長さが `< 2·高さ`。
  - どの atom の親の高さも `≥ t_U`。

### 2.10 段 S9：atom の置き換え（l.2099、Lemma 5.15）

```cpp
AtomSubstitution<W> sub = substitute_atoms(ph, &diag, params.local);   // atom ごとに Alg 4 を1回
```

- すでにある関数: `substitute_atoms`（`src/paper/local/atom_replacement.hpp`）。atom ごとに `local_hierarchy` を呼ぶ。Alg 4 が基本ケースを取ると、`local_hierarchy` が `FallbackRecord{"local: k <= 2*L_k^4 (paper base case)", k, m}` を積む。
- 不変条件: `validate_ordered_hierarchy(sub.hierarchy, 縮約後のグラフ)` が、`c0 ≤ 2` で違反を報告しない（`docs/TASKS.md` の受け入れ条件）。`sub.old_state` と `projection` を合成すると、新しいどの状態についても物理的な像が得られる。atom 1個につき `local_hierarchy` の呼び出しは1回（`atom_replace_atoms`）。

### 2.11 主張されている計算量

【論文、Lemma 8.4】意味のある局所の事象と有限の join の分類は、すべての大域の格子を含めて `O(m + N·τ(N))`。置き換えには、Lemma 5.15 が `O(m + N log log(b+4))` を加える。

---

## 3. Algorithm 8 l.8〜17（Sec 8.3）

### 3.1 段、入力、出力

`c` は Alg 8 l.2 の `outer::ZeroContraction<W>`（重みは正、ループ無し、連結、多重辺は残る）。`g = c.graph`、`N = g.num_vertices()`。

| Alg 8 | 段 | 呼び出し（印が無ければ `main` にある） | 入力 → 出力 |
|---|---|---|---|
| l.8 | MST（Lemma 3.5）、owner セル、木に厳密に支配される非木辺の削除 | `owner_decomposition(g, b, &diag)` | `g`、`b` → `OwnerDecomposition<W> od`（`tree_edges`、`reduced`、`pivots`、`owner`、`radius`、`owner_distance`、`corridors`） |
| l.9 | 射影した owner 間の辺への NameRoutes | `owner::classify_interowner_edges(g, od, rng, name_routes_params, &diag)` | → `OuterEdgeClassification<W> cls`（`name`、`dictionary`、`p_log2`、`detour` = d_e） |
| l.10 | マクロと半径辺、`G'`、`T'` | （分類は上の呼び出しに入っている）`owner::build_aux_graph(g, od, cls, &diag)` | → `AuxGraph<W> aux` |
| l.11 | OwnerDistances（Alg 7） | `owner::owner_distances(aux, &diag, {.local = local_params})` | → `J`（`distance`、`scaled_distance`）、`calendars` |
| l.12 | `w_e > J_x + J_y + d_e` の半径辺 `e = xy` を削除 | **P3-OUTER-STATES** | → 残る辺の集合 |
| l.13a | 外側の状態、θ、join（Sec 8.1） | **P3-OUTER-STATES** | → A-1 |
| l.13b | 名前、F、尺度、格子、分割階層（Sec 8.2） | **P3-OUTER-SCALES**（§2、S0〜S8） | → `PartitionHierarchy<W> ph` |
| l.13c | 上限つきの局所階層 | `substitute_atoms(ph, &diag, local_params)` | → `AtomSubstitution<W> sub` |
| l.14 | 種 | `redirect_offers(sub, {{native_inward[c.source], 0}})` | → `s` の native な内向きコピーへの `StateOffer<W>` 1個。`s` が atom の中にあれば、その native な局所の副コピーへ |
| l.15 | Evaluate、頂点ごとの最小値 | `evaluate(SeededOrderedHierarchy<W>{sub.hierarchy, {{c.source, 0}}, offers})` | → 新しい状態ごとの `EvalResult<W>::value`。`dist[v]` は、`projection[sub.old_state[x]] = v` である新しい状態 `x` についての最小値 |
| l.16〜17 | 復元、tight な弧、最初に発見した弧の木 | `outer::restore_predecessor_tree(元のグラフ, s, c, dist)`（`paper_sssp` がすでに呼んでいる） | → `SsspResult<W>` |

注:
- l.8: 論文の「頂点対ごとに最も軽い辺を残す」は、基本ケース（l.5）だけのもの。外側の枝は多重辺を残す。NameRoutes は、自分のインスタンスの中で重複を除く（Alg 5 l.4）。
- l.9 が、乱数を使う唯一の段（Lemma 6.8）。生成器は `model::WordRng&` として明示的に渡す。
- l.15【解釈、Q-73】: 最小値は、`v` に射影されるすべての状態について取る。特別な状態と、すべての局所の副コピーを含む。Lemma 4.6 により、過小評価は無く、正確なコピーが少なくとも1つある。有限の値を出した状態が1つも無い頂点は、不変条件が破れていることを意味する。`std::logic_error` を投げる。始点の成分の頂点に無限大を返してはいけない。
- `outer_construction` が返すのは、縮約後の距離だけ。木は Evaluate からは出ない。

### 3.2 Lemma 8.5（l.2145〜2157）

【論文】零辺の成分を縮約したあと、縮約した始点成分の正確な距離から、元のすべての距離と最短路木が `O(m + n)` で決まる。`d(u) + w(uv) = d(v) < ∞` のとき `uv` を `u` から `v` に向け、`s` からたどり、最初に発見した弧を木とし、木を下りながら元の重みを足した和を出力する。

これは実装済みで、基本ケースの枝で VERIFIED になっている（`restore_predecessor_tree`、PAPER_MAPPING の「Sec 8.3, Lemma 8.5」の行）。距離だけを受け取るので、Alg 8 l.16 に新しいコードは要らない。外側の枝が守ること:
- 入力のラベルは正確でなければならない。有限のラベルが tight な元の弧で達成されないと、`restore_predecessor_tree` は `std::invalid_argument` を投げる。`double` では、外側の距離は木の和とは別の順序で足した和なので、2進有理数でない入力はここで失敗しうる（Q-06、Q-55 の案 A、Q-63）。Mode B ではそれは**失敗**であり、別の探索で置き換える理由にはならない。
- 到達できない元の頂点は `c.vertex_map` の外にあり、復元で無限大になる。外側の構成はそれらを見ない。

### 3.3 fallback を記録する場所

| 場所 | 記録 | 積むもの |
|---|---|---|
| Alg 8 l.4〜5（`N ≤ L_N^64`） | `kOuterBaseCaseScanReason` / `kOuterBaseCaseD1Reason`、`(N, m')` | `outer_base_case*`（ある） |
| Alg 7 l.6: Sec 7.3 の各 atom への Alg 4 | `"local: k <= <c>*L_k^<e> (paper base case)"`、`(k, m_Y)` | `cap_local_hierarchy` 経由の `local_hierarchy`（ある） |
| Alg 8 l.13 / Lemma 5.15: Sec 8.2 の各 atom への Alg 4 | 同じ理由、`(k, m_Y)` | `substitute_atoms` 経由の `local_hierarchy`（ある） |
| Alg 5 の打ち切り（D-10、Mode B だけ） | NameRoutes の記録 | `name_routes`（ある） |
| `outer_construction` そのもの | **無し** | 段を直接の探索で置き換えることは無い。扱えない入力では例外を投げる |

- 【実測】Alg 7 は、頂点1個の atom にも記録を積む（`k = 1`、`m = 0`）。§3.5 の6頂点の例では、そのような記録が4個できた。
- 1つの頂点が、複数の記録に入りうる。Alg 7 で1回、そのコピーを含む Sec 8.2 の atom ごとに1回（コピーは最大4個: `in_i`、`in_{i−1}`、`out_i`、`out_{i+1}`）。6頂点の例では、組 `{x, y}` が3つのブロックで atom になる。したがって `FallbackRecord::vertices` の和は `N` を超えうる。fallback の割合には定義が要る（Q-69）。
- 【補完】カウンタの案（fallback ではない）: `outer_constructions`、`outer.pivots`、`outer.macro_edges`、`outer.radial_edges_deleted`、`outer.states`、`outer.blocks_{singleton,aligned,unaligned,unattached}`、`outer.atoms`、`outer.atom_states`、`outer.grid_points`、`outer.grid_components`、`outer.global_classes`、`outer.separate_classes`、`outer.scale_factor_log2`。
- 退化した場合 `t = 1`（マクロ無し、どのブロックも付着無し、`J` がピボットからの全距離）は論文の構成であって、fallback ではない。`outer.pivots = 1` で報告する（Q-72）。

### 3.4 Mode A と Mode B

Mode A: `takes_outer_construction(N)` は、どの `N < 2^64` でも偽（`L_N ≥ 2` なので `L_N^64 ≥ 2^64 > N`）。Alg 8 l.8〜17 には `paper_sssp` からは届かない。届く最小の `N` は `589^64 + 1`（Q-04、`thresholds.md` §2）。したがってこのコードを動かせるのは、(i) Mode B のパラメータと、(ii) `outer_construction` を直接呼ぶ場合だけである。

【導出】小さなグラフで `outer_construction` を論文の定数のまま直接呼ぶことは、きちんと定義されている。`b` は `N` より上で飽和し、分解のピボットは1個、マクロ辺は無く、状態2個以上のブロックはすべて付着無しになり、Alg 7 が `J` をピボットからの全距離として計算する。これは Mode A の実体化のテストとして使える（§6、PL-2）。ただし Sec 8.2 は、付着無しの経路しか通らない。

外側の枝を、実質のある仕事つきで強制する、Mode B の最小のパラメータ集合（`paper_experimental::PaperParams`）:

| パラメータ | 値 | 論文 | 逸脱 |
|---|---|---|---|
| `outer.threshold` | `SizeThreshold{0, 0}`（どの `N ≥ 1` にも適用される） | `1·L^64` | D-7 |
| `outer.separator` | `SeparatorRule::constant(b)`、`b ∈ {1, 2, 3, …}` | `L^64` | D-7 |
| ほかすべて | 論文の値: 標本 `16t`、残り `32t`、論文の `P`、論文の `H`、打ち切り無し、局所の閾値 `2·L^4`、局所の `b = L^4` | — | 無し |

導かれる値: `h = b + 1`、`F` は §1 から（たとえば `b = 3`、`t = 2` で `f = 48`）。この値では、外側の atom はどれも頂点が `b` 個未満で、Alg 4 の基本ケースを取る（記録される）。NameRoutes は標本を取らない（`t ≤ 65` では `|U| ≤ 32t`）。したがって owner 間のどの辺も自分の重みを証明書にし、`d_e = P·A_e > w_e` となって、**owner 間の辺はすべてマクロ**になる。半径辺は、ちょうど owner の中の辺である。

通る範囲を広げる、任意のつまみ（どれも一覧にある逸脱）:

| 目的 | つまみ | 逸脱 |
|---|---|---|
| NameRoutes の標本の回。それによって、owner 間の半径辺と Alg 8 l.12 の削除 | `sample_factor = 1`、`residual_factor = 2`、小さな `p_log2`、`retry_cap` | D-17、D-10 |
| `uint64` での unaligned のブロック | `p_log2 = 0` と `b ≤ 3`（`f ≤ 31`） | D-17 |
| 外側の atom の中での、Alg 4 の半径の枝 | 外側の `b ≥ 33`。ある atom について `k ≥ max(32, b_local + 1)` となる、局所の閾値と局所の `b`（Q-03） | D-7 |

### 3.5 外側の構成が動く最小の N

「動く」には、いくつかの意味がある。`N` は縮約後の大きさ。外側の閾値は無効にする。NameRoutes と局所のパラメータは論文の値。

| 水準 | 意味 | 最小の N | b | 例 |
|---|---|---|---|---|
| L0 | 枝に入る | 1 | 任意 | 頂点1個。ピボット1個。ほかは何も無い |
| L1 | ピボット2個、corridor 1本、マクロ辺1本、辺のある NameRoutes のインスタンス | **2** | 1 | 辺1本。両端がピボット |
| L2 | L1 に加えて、ピボットでない頂点（正の通常の状態、自明でない `J`） | **4** | 2 | 道 0–1–2–3 |
| L3 | L2 に加えて、有限の付着を持つ、状態2個以上のブロック（aligned の判定が走る） | **4** | 2 | 道 0–1–2–3（重み 10、3、3）と、辺 1–3（重み 5） |
| L4 | L3 に加えて、状態2個以上の atom（Alg 4、上限づけ、置き換えが実際に働く）。**Sec 8.2 のすべての段が自明でない**。【訂正、Q-78】この例では唯一のマクロがピボット同士を結ぶので、状態2個以上のブロックはすべて unattached で、alignment の判定は動かない（L3 を満たさない）。SC-8 の入力（辺 2–5 を足したもの）がそれを補う | **6** | 3 | 辺 0–1 (4)、1–2 (4)、0–3 (100)、3–4 (8)、4–5 (1) |

下界【導出】: 状態2個以上の atom は `T − A` の1つの成分に入り、その成分の頂点は `b` 個未満なので、`b ≥ 3`。ピボットが2個あるには、印のついた頂点が2個要る（印が1個では分岐頂点が無い）。印1個は、自分の頂点を少なくとも `b` 個使う（Sec 5.1 の印のつけ方）。したがって `N ≥ 2b`。L4 では `N ≥ 6`。`b = 1` ではすべての頂点に印がつくので、ピボットでない頂点には `b ≥ 2` が要る。したがって L2 と L3 では `N ≥ 4`。

【実測】`main` にある部品（`owner_decomposition`、`classify_interowner_edges`、`build_aux_graph`、`owner_distances`、および l.12 の削除のあとの、座標 `J`・`K = 16` での `build_radial_states`）を、上の例と、`(N, b)` ごとに 60,000 個の無作為な連結多重グラフ（`N ≤ 7`、`b ≤ 4`）で実行した:
- L4 の例: ピボットは {0, 3}。マクロは1本（0–3）。頂点 4 と 5 で `r = 8`、`J = 8, 9`。ブロック `{4, 5}` が、内向きの窓 0 と、外向きの窓 0 と 1 にできる。どれも `s_U = 8` で、重み 1 の辺は `s_U/2` より軽い。ブロック `{1, 2}` は `s_U = 4` で、辺の重みは 4（意味のある辺。atom にならない）。【訂正、Q-78】これらのブロック6個はすべて unattached。唯一のマクロがピボット同士を結ぶので、有限の付着を持つブロックが無い。
- L3 の例: ピボットは {0, 2}。マクロは 0–1。頂点 1 と 3 で `r = J = 3`。ブロック `{1, 3}` が3つでき、`s_U = 3`、辺の重みは 5。
- 無作為な探索: `N ≤ 3` では、ピボットが2個で状態2個以上のブロックを持つ例は無かった。`N ≤ 5` では、ピボットが2個で、窓の中に `s_U` より軽い辺を持つ例は無かった。どちらも `N = 4, b = 2` と `N = 6, b = 3` で現れる。
Sec 8.1 の特別な状態、θ、それより後ろは実行していない（`main` に無い）。

L4 より先（測っていない。上界と下界だけ）:
- 論文の定数で NameRoutes が標本を取るには、相異なるピボットの対が `32t` 個より多く要る。したがって `t ≥ 66`（`C(66, 2) = 2145 > 2112`）。`b = 1` なら `N = 66`（完全グラフ。ピボットでない頂点は無い）。L4（`b = 3`）と合わせると `102 ≤ N ≤ 198`（下界は `t ≤ 2|A0| − 1` と `|A0| ≤ N/b` から。上界は、頂点3個の印を66個）。正確な最小値は計算していない。
- D-17（標本 `1t`、残り `2t`）なら、標本を取るには `t ≥ 6`。L4 と合わせて `12 ≤ N ≤ 18`。
- 外側の atom の中で Alg 4 の半径の枝を動かすには、頂点が32個以上の atom が要る。したがって外側の `b ≥ 33`、ピボット2個で `N ≥ 66`。

---

## 4. すでにある部品で足りないもの

| ID | 足りないもの | 必要とする側 | 詳細 |
|---|---|---|---|
| M-1 | Sec 8.1 の出力の型 | 尺度、流れ | P3-OUTER-STATES は作業中。必要なフィールドは前提 A-1（§2.0）。とくに、状態ごとの `b0_edge`（または辞書の添字である `b0_key`）と、ブロックの種類 |
| M-2 | `src/paper/` での正確な `f` | 尺度 | 多倍長の自然数（`experimental/paper_experimental/params.cpp` の `Nat` はファイルの中だけ）と `outer_scale_factor_log2(p_log2, 正確な h)`。正確な `b = L_N^64`（いまの `outer_separator_size` は飽和する）。そのうえで、Mode B の `OuterDerived::f_log2` も同じ関数を呼ぶ（Q-04: 実装は1つ） |
| M-3 | 外側の枝のパラメータと生成器 | 流れ | `paper::outer_construction(const ZeroContraction<W>&, Diagnostics&)` にはどちらも無い。必要なもの: `src/paper/outer/` のパラメータの構造体（分離子の決め方、`pivots::NameRoutesParams`、`LocalHierarchyParams`）と `model::WordRng&`。`PaperOptions` の seed。`PaperParams` をそれに写す Mode B の入口（P5-EXP-THRESHOLDS） |
| M-4 | 新しい名前 | 尺度 | 辞書の添字ごとの `name_weight`（§2.1） |
| M-5 | 外側の分類を `EdgeClassification` に | 外側の状態 | `OuterEdgeClass` → `EdgeClass`。l.12 で削除した辺は除き、`macro_key = name` にする（おそらく P3-OUTER-STATES の中） |
| M-6 | 使い回せる、格子・support・同時併合・リンク・子の順序の仕組み | 尺度 | Sec 5.4 のコードは `local/assemble.cpp` の無名名前空間にあり（`MergeBuilder`、`floor_index`、`bucket`）、`OwnerDecomposition` の corridor と `N^c`、`N^16` に結びついている。Mode B には、すでに2つめの写しがある（`local_assembly.cpp`、Q-59、P5-LOCAL-UNIFY）。共有の仕組みを切り出すか、3つめの写しを明示的に認めるか（lead の判断） |
| M-7 | 分割階層の出力 | 尺度 | `assemble_local_hierarchy` が返すのは、小さな成分に正確な Kruskal 木を持つ `OrderedHierarchy`。Sec 8.2 に要るのは、atom の葉を持つ `PartitionHierarchy<W>`（`leaf_atom`、`ReplaceableAtom`） |
| M-8 | 入力の単位での別の格子 | 尺度 | `Calendars::graph_pieces[].grid` は `2^7` 倍で、`4·max R` までしか無い。Sec 8.2 に要るのは、入力の単位で `4K²·R_v` までの点（§2.6）。また、参照のための、頂点 → 補助頂点の対応（`aux.aux_vertex`） |
| M-9 | 検査器 | テスト | `check_outer_scales(...)`。§2 の不変条件の違反ごとに、メッセージを1つ返す。`check_aux_metric` と `check_calendars` と同じ流儀 |
| M-10 | 外側の atom 用の局所の変種 | Mode B | `substitute_atoms` と `owner_distances` が受け取るのは `LocalHierarchyParams`（大きさの判定と分離子）だけ。atom に D-8 と D-11 の仕組みは選べない。atom が基本ケースを取るあいだは関係しない |
| M-11 | Diagnostics の受け渡し | 流れ | `absorb` は `sssp.cpp` の中だけの関数。`evaluate` は自分の `EvalResult::diag` を返すので、合わせる必要がある |
| M-12 | コピーについての最小値 | 流れ | 新しい状態 → `sub.old_state` → `projection` → 頂点。線形の走査。加えて、Q-73 の「値を出した状態が無い」ことの検査 |

そのまま使える部品: `owner_decomposition`、`classify_interowner_edges`、`build_aux_graph`、`owner_distances`、`sat::add` / `sat::mul_pow2`、`pivots::median_prune`、`DisjointSets`、`ForestQueries`、`PartitionHierarchy` / `ReplaceableAtom` / `substitute_atoms` / `redirect_offers`、`validate_ordered_hierarchy`、`evaluate`、`contract_zero_edges`、`restore_predecessor_tree`。

---

## 5. 曖昧な点：未解決の問いの案（Q-65 から）

これらは `open-questions.md` への案である（そのファイルは lead が編集する）。どれも決まっていない。

| ID | 場所 | 問い | 候補 | 案 |
|---|---|---|---|---|
| Q-65 | Sec 8.2 l.2077〜2083 | (a)「その名前を持つマクロの、最小の物理的な重み」は、入力の重み `w_e` か、射影した長さ `A_e` か。どの辺について取るか。(b) corridor の join は、自分の重みの support を使うか、マクロとしての名前の support を使うか | (a1) マクロ辺についての `w_e`。(a2) `A_e`。(b1) 自分の重み。corridor ごとにキー1つ。(b2) 名前 | (a1): Lemma 7.1 は `w` の上下界を述べている。(b1): 本文は、名前と並べて「すべてのピボット corridor の重み a」を挙げている |
| Q-66 | Sec 8.2 l.2085〜2093 | 「それ以外の、状態2個以上の正の通常ブロック」は、どのブロックを指すか | (1) unaligned のブロックと、有限の付着が無いブロック。(2) unaligned だけ。付着の無いブロックは定められていない | (1)。Sec 5.4 l.1213 が局所の場合について述べているのと同じ |
| Q-67 | Sec 8.2 l.2100〜2104 | 「v の補助格子へのポインタ」は何か。単位は何か。オクターブをいくつ足すか | (1) Sec 7.3 の graph piece の格子。origin からの倍加で入力の単位に直し、上に倍加9回ぶん伸ばす。(2) 同じ格子を `2^7` 倍の単位のまま使う（ブロックの重みをすべて拡大する。D-20 の範囲の制限がかかる）。(3) Sec 5.4 のように、`s_U` から始まるブロックごとの格子で `t_U = s_U`（Sec 8.2 の記述とは違う。逸脱） | (1)。オクターブの数はこちらで決めたもの。論文は「決まった数」としか言っていない |
| Q-68 | Sec 8.2 l.2089〜2099、Lemma 8.4 の証明 (2) | (a) 隙間の判定は厳密な不等号（`next.lower < 2·upper`）か。(b) 類 `[a, 2a)` の節点の高さと、リンクの上界。(c) 別の格子の類と大域の類を処理する順序 | (a) 厳密 / 等号を含む。(b) 高さ `a`、リンク `< 2a`、`c0 = 2`（Lemma 5.10 と同じ） / 高さ `2a`、リンク `< 2a`。(c) 別の格子が先 / 値の順に混ぜる | (a) 厳密（「2倍より小さい」）。(b) 高さ `a`。(c) 別の格子が先。高さと矛盾しないどの順序とも同じ結果 |
| Q-69 | CLAUDE.md の fallback の方針、Alg 7 l.6、Lemma 5.15 | 1つの頂点が複数の局所の基本ケースに入るとき（Alg 7 の atom と、コピーを通じた最大4個の Sec 8.2 の atom）、fallback の割合をどう定義するか | (1) 理由ごとの生の和。「N を超えうる」と明記する。(2) 相異なる物理的な頂点の数（追加の記録が要る）。(3) 基本ケースの atom に入る状態の数 ÷ 全状態の数。段ごとに報告する | (3) を段ごとに。加えて (1) の生の値。lead の合意が要る |
| Q-70 | Alg 8 l.9、`include/sssp/paper.hpp` | 外側の枝の生成器とパラメータは、Mode A の API のどこから入るか | (1) 固定の既定値を持つ `PaperOptions::seed` と、内部の `WordRng`。(2) `WordRng&` を取る別の関数。(3) 公開の入口は作らず、テストが `outer::outer_construction` を直接呼ぶ | (1) と (3) の両方。公開ヘッダの変更なので、lead が決める |
| Q-71 | Sec 8.2 l.2080〜2081、Sec 7.1 l.1644 | Mode A では `b = L_N^64` を飽和した値で持っている。`F` はどの `h` から出すか | (1) 多倍長の正確な `h`（M-2）。(2) 飽和した `b ≥ N` のときは `F` を飛ばす（そのとき `t = 1` で、`F` は使われない） | (1)。(2) は、論文に無い分岐を足すことになる |
| Q-72 | `t = 1` のときの Alg 8 l.8〜13 | ピボットが1個の実行（マクロ無し、`J` は Alg 7 による全距離）は fallback か | (1) 違う。論文の構成である。カウンタだけ。(2) FallbackRecord を積む | (1) |
| Q-73 | Alg 8 l.15 | 「頂点ごとに、出た値の最小値」は、どの状態について取るか。値を出した状態が無い頂点はどうするか | (1) その頂点に射影されるすべての状態（通常、特別、局所の副コピー）。無ければ `std::logic_error`。(2) 通常のコピーだけ | (1): 特別な状態は物理的な歩道に射影されるので、過小評価できない |
| Q-74 | Sec 8.1 l.2000〜2001、Def 5.11〜5.12 | (a) Sec 8.2 の atom の内部グラフは、どの辺から成るか。(b) 状態1個の成分は atom か | (a1) 両端が atom にある、窓の中のすべての辺（重い弦も含む。Q-53 と同じ）。(a2) `< t_U` の辺だけ。(b1) 単独の部分。Alg 4 は呼ばない。(b2) 状態1個の atom（いまの Alg 7 と同じ） | (a1)。`substitute_atoms` が要求する。(b1)。Def 5.12 に従う |

確かめて、矛盾が無かった箇所（問いにしない）: Lemma 8.4 の証明 (1) の aligned の上下界（`a_U/F² ≤ s_U ≤ F·a_U`）は、Lemma 8.2 と Lemma 7.1 から、`8K²h` 倍の余裕つきで出る。`s_U/2 < t_U ≤ s_U` は、どちらの格子でも成り立つ。atom の cutoff `t_U` は Def 5.11 を満たす（全域木の辺は `< t_U < 2t_U`、親の高さは `≥ t_U`）。

---

## 6. 受け入れテストの案

重み: 断らないかぎり `uint64`。`double` は2進有理数の値だけ（Q-55 の案 A）。NameRoutes が標本を取るところでは、seed を固定し、かつ変える。オラクルは独立したものを使う。元のグラフでの Bellman–Ford、状態弧での Tarjan の SCC、`f` には Python の多倍長整数。

### 6.1 P3-OUTER-SCALES

ファイルの案: `tests/unit/paper/outer/test_outer_scales.cpp`、`tests/property/paper/test_outer_scales_property.cpp`。

| テスト | 入力 | 確かめること |
|---|---|---|
| SC-1 `scale_factor_log2 matches big-integer values` | §1.4 の `(p, b)`。`b = 589^64`（1808）、`20^64`（871）、`b = 3, p = 17`（48）、`b = 160000, p = 1`（77）を含む | `f` が表に等しい。テストの中で、別の多倍長の手続きで `2^f > X ≥ 2^(f−1)` を計算し直す |
| SC-2 `scale factor is not derived from a saturated h` | `N = 10^6`、論文の `b` | 871 であり、233 ではない |
| SC-3 `F^2 s >= a is decided exactly` | `uint64` と2進有理数の `double` での無作為な `s`、`a`、`f`。飽和の境界の近く | 正確な答え（`unsigned __int128` のシフト、または2冪の指数の計算）に等しい |
| SC-4 `new names are the least macro weight` | 名前2つ、マクロ4本の手作りの分類 | `name_weight`。`a ≤ w < 2P(2h+1)a`。有限のどの join でも `a ≤ B0 ≤ θ < F·a` |
| SC-5 `hand example L3` | §3.5 の L3、`b = 2` | 状態2個以上のブロック `{1, 3}` が3つ。`s_U = 3`、aligned、`t_U ∈ (1.5, 3]`、atom 無し。重み 5 の辺が、それを挟む大域の類に分類される。【訂正、Q-78】3つのブロックのうち aligned は2つだけ（内向きの窓 0 と外向きの窓 0、`t_U = 3`）。外向きの窓 1 のブロックには入る弧が無く、I0 が無限大なので unattached（別の格子で `t_U = 3`） |
| SC-6 `hand example L4` | §3.5 の L4、`b = 3` | atom `{4, 5}` が3つのブロックに。cutoff は `t_U ∈ (4, 8]`。`H` は重み 1 の辺1本。ブロック `{1, 2}` に atom は無い。`ph` を `substitute_atoms` が受け入れる。`validate_ordered_hierarchy` が `c0 ≤ 2` で違反無し。【訂正、Q-78】唯一のマクロがピボット同士を結ぶので、状態2個以上のブロック6個はすべて unattached（別の格子。`t_U = 8` の atom `{4, 5}` はそのまま） |
| SC-7 `unaligned block on uint64` | L3 のグラフで corridor の重みを `2^60` に。`p_log2 = 0`、`b = 2`（`f = 29`） | ブロック `{1, 3}` が unaligned。窓の中の辺はすべて `< a_U/2`。ブロックの上で θ が等しい。`t_U` は piece の格子の上。`2·p_top ≤ g` |
| SC-8 `unaligned block with an atom` | L4 のグラフで、重みを 0–3: `2^63 + 1`、3–4: 2、4–5: 1 に。`p_log2 = 0`、`b = 3`（`f = 31`）。【訂正、Q-78】この入力のままでは、ブロック `{4, 5}` は unattached で、unaligned にならない。テストでは、非木辺 2–5（重み `2^63 + 2`）を足してブロックにマクロの端点を持たせ、名前を手で与えた | ブロック `{4, 5}` は unaligned（`F²·2 = 2^63 < a_U`）。それぞれが cutoff 2 の atom 1つ。ブロック `{1, 2}` は aligned。【訂正、Q-78】辺 2–5 を足すと、`{4, 5}` の内向きのブロックと外向きの窓 0 は unaligned、外向きの窓 1 は unattached で、それぞれ cutoff 2 の atom 1つ。ブロック `{1, 2}` は、付着のある2つが aligned、外向きの窓 1 は unattached |
| SC-9 `several grid components (double)` | corridor の重み2つが `3f` オクターブより離れた、2進有理数の重み | 成分が2つ。各成分の中で `grid[i+1] = 2·grid[i]`。成分の始点をまたぐ類は無い。隙間が2倍未満の support は併合される |
| SC-10 `unattached blocks` | ピボット1個（`b > N`） | 状態2個以上のブロックはすべて付着無し。別の格子だけ。根が複数あり、凝縮の順に並ぶ |
| SC-11 性質: Lemma 8.4 | 無作為な連結グラフ、`N ≤ 200`、`b ∈ {1, 2, 3, 4, 8, ⌈log2 N⌉}`、ハーネスの生成器の系統、複数の seed | `check_outer_scales` を通じた §2 のすべての不変条件。unaligned の MST 辺は `< a_U/2`。aligned の意味のある重みは support の中。`s_U` は元の半径の最小値。`s_U/2 < t_U ≤ s_U`。atom は頂点が `b` 個未満、owner が1つ、`r` が1つ。階層が Def 4.1 を `c0 ≤ 2` で満たす |
| SC-12 性質: SCC の一致 | 小さな `N`（≤ 40） | どの内部節点についても、下の状態は、類の上端より軽い弧の SCC 1つ。子は、高さより軽い弧の SCC（Tarjan） |
| SC-13 `alignment by median pruning equals the per-block test` | 無作為なまとまり | 分類が同じ。まとまりごとの述語の判定は `≤ ⌊log2 n⌋ + 1` 回 |
| SC-14 `distances on the hierarchy` | SC-11 の入力 | 置き換え後の階層で、Alg 8 l.14 の種を使って `evaluate` する。頂点ごとの最小値が Bellman–Ford に等しい |

### 6.2 P3-OUTER-PIPELINE

ファイルの案: `tests/unit/paper/outer/test_shortest_paths.cpp`、`tests/differential/paper_experimental/test_outer_forced_vs_dijkstra.cpp`（後者の名前は、すでに `PAPER_MAPPING.md` にある）。

| テスト | 入力 | 確かめること |
|---|---|---|
| PL-1 `Mode A never takes the outer branch` | 差分テストの系統 | `takes_outer_construction` が偽。基本ケースの FallbackRecord がちょうど1個（いまの振る舞いを保つ） |
| PL-2 `outer_construction with the paper's constants` | 小さなグラフ（`N ≤ 300`）での直接の呼び出し | 距離が Bellman–Ford に等しい。`outer.pivots = 1`。マクロ無し。理由が `outer:` で始まる FallbackRecord は無い |
| PL-3 `forced outer, levels L0–L4` | §3.5 の5つの例 | 距離と木がオラクルに等しい。`verify_sssp`。その水準の特徴がカウンタに現れる（ピボット、マクロ辺、状態2個以上のブロック、atom） |
| PL-4 `forced outer vs oracle` | 生成器のすべての系統、`N ≤ 2000`、`b ∈ {1, 2, 3, 4, 8, 16, 64}`、seed を固定し、かつ変える | 距離が Bellman–Ford と baseline の Dijkstra に等しい。`verify_sssp`（根は 0、どの弧も実行可能、木の弧は tight、親をたどると `s` に着く、到達できない頂点は無限大） |
| PL-5 `stage checkers pass on the intermediate states` | PL-4 の入力 | `check_route_names`、`check_aux_metric`、`check_calendars`、`check_atoms`、P3-OUTER-STATES の検査器、`check_outer_scales`、`validate_seeded_ordered_hierarchy` がすべて空 |
| PL-6 `fallback records` | PL-4 の入力 | 外側の基本ケースの記録は無い。どの記録も局所の基本ケースの理由（打ち切りを設定したときは D-10 の理由も）。基本ケースを取った Alg 4 の呼び出し1回につき記録1個。パラメータの集合が記録されている（`PaperParams::record`） |
| PL-7 `zero edges, unreachable vertices, loops, parallel edges` | 零の閉路と、複数の成分を持つ入力 | 復元が正しい。外側の枝が見るのは縮約後の `N` |
| PL-8 `seed independence` | NameRoutes が標本を取る入力（D-17 のパラメータ、`t ≥ 6`）で seed 16個 | どの seed でも距離が同じ。`aux_metric.radial_interowner_edges > 0` かつ `outer.radial_edges_deleted > 0` となる例が少なくとも1つある |
| PL-9 `source position` | 始点がピボット、atom の外のピボットでない頂点、atom の中 | すべてオラクルに等しい（Alg 8 l.14 の付け替え） |
| PL-10 `double weights` | 2進有理数の `double`、2進有理数でない `double` | 2進有理数: 正確な基準と厳密に等しい。そうでない場合: 文書化した許容誤差の中で等しいか、例外になるかのどちらか。黙って誤った結果を返すことは無い（Mode B の失敗の決まり、Q-63） |
| PL-11 `unaligned end to end` | SC-7 と SC-8 の入力 | 距離がオラクルに等しい |
| PL-12 `sanitizer` | PL-3 と、大きさを減らした PL-4 | `make sanitizer` で問題無し（P4-SANITIZER） |

`docs/TASKS.md` の受け入れ条件に足す文の案（lead 向け）:
- P3-OUTER-SCALES:「`f` は文書 04 §1.4 の多倍長の表に等しく、飽和した `h` から出すことは無い。unaligned の枝は、`p_log2 = 0` の `uint64` のテストと、2進有理数の `double` のテストで通す。性質テストのすべての入力で `check_outer_scales` が空」。
- P3-OUTER-PIPELINE:「文書 04 §3.5 の水準 L0〜L4 を通す。論文の定数での直接の呼び出し（ピボット1個）がオラクルに等しい。局所の基本ケース以外の FallbackRecord は無い。値を出した状態が無い頂点では例外を投げる」。
