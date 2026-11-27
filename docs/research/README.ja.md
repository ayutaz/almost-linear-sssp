# 論文の解析文書(docs/research/)

[English](README.md) | 日本語

対象は arXiv 2611.34911v1「Linear Single-Source Shortest Path, Almost」。論文の本体は再配布できないので、`.cache/paper/`(git の管理外)に `scripts/fetch_paper.sh` で取得する。ここにある文書は論文を言い換えて解説したもので、本文の転載はしていない。行番号 `l.NNN` は `.cache/paper/2611.34911v1.txt` のもの。

## 索引

| 文書 | 範囲 | 中身 |
|---|---|---|
| [01_preliminaries_primitives.md](01_preliminaries_primitives.md) / [01_preliminaries_primitives.ja.md](01_preliminaries_primitives.ja.md) | Sec 1–3 | 主定理、記号の一覧、閾値の使われ方、基本の場合の規約、計算モデル(長さの演算、整数、乱数)、double で問題になる箇所、Sec 3 のプリミティブ P1–P13(探索、整列、選択、Dijkstra、Fibonacci ヒープ、成長パスと MST、union-find、LCA、level ancestor、split–findmin)、参照実装との対応、モジュール案 |
| [02_ordered_hierarchy_local.md](02_ordered_hierarchy_local.md) / [02_ordered_hierarchy_local.ja.md](02_ordered_hierarchy_local.ja.md) | Sec 4–5 | 順序付き階層(Def 4.1–4.2)、正規化、子の順序、Scan と Evaluate(Alg 2–3)、検査できる不変条件 I1–I8、owner 分解、radial 状態表現、閾値成分、局所階層(Alg 4)、原子の置き換え、閾値を下げたときに壊れる補題、モジュール案と実装の順序 |
| [03_pivots_aux_outer.md](03_pivots_aux_outer.md) / [03_pivots_aux_outer.ja.md](03_pivots_aux_outer.ja.md) | Sec 6–8 | ピボット経路の証明書(Alg 5–6、Lemma 6.2–6.8)、補助計量と owner 座標、幅クラス・カレンダー・原子、OwnerDistances(Alg 7)、外側の状態グラフと数値尺度、ShortestPaths(Alg 8)と主定理の計算量の内訳、中間状態の不変条件の一覧 |
| [04_outer_scales_pipeline.md](04_outer_scales_pipeline.md) / [04_outer_scales_pipeline.ja.md](04_outer_scales_pipeline.ja.md) | Sec 8.2–8.3 | 外側の尺度(Lemma 8.4)と外側の流れ(Alg 8 l.8–17、Lemma 8.5)の実装のための仕様: 正確な (値, 指数) の形の F、Lemma 8.4 の段ごとの説明、Alg 8 の段と fallback を記録する場所、Mode A / Mode B、既存の部品に足りないもの(M-1〜M-12)、問いの案 Q-65〜Q-74、受け入れテストの案 SC-1〜14 と PL-1〜12 |
| [thresholds.md](thresholds.md) / [thresholds.ja.md](thresholds.ja.md) | 閾値と定数 | 全定数の一覧と種別(正しさに効くか、時間だけか)、各閾値を満たす最小の N、現実の N で通る経路、Mode B で値を変えたときに壊れるもの、実装への推奨 |
| [prior_work_and_oss.md](prior_work_and_oss.md) / [prior_work_and_oss.ja.md](prior_work_and_oss.ja.md) | 先行研究と既存 OSS | 既存の実装(対象論文の実装を含む)のライセンス・検証方法・計測条件、先行研究、取り込むか再実装するかの推奨 |
| [open-questions.md](open-questions.md) / [open-questions.ja.md](open-questions.ja.md) | 全体 | 上の文書の曖昧点を ID(Q-xx)で 1 つの表にまとめたもの。暫定の解釈、影響するモジュール、決着の状態(解釈で決定 / 要合意 / 未解決) |
| [../PAPER_MAPPING.md](../PAPER_MAPPING.md) / [../PAPER_MAPPING.ja.md](../PAPER_MAPPING.ja.md) | 全体 | 論文の全項目とモジュール・テスト・Mode・状態・逸脱(D-xx)の対応表。モジュール名の統一の記録もここにある |

## 読む順番

1. **`thresholds.md` の §2–§3** — まず「論文どおりの閾値では、現実のどの入力でも外側の構成が動かず、走査型 Dijkstra の基本の場合になる」ことを押さえる。Mode A / B / C の分け方(ADR 0003)と、基本の場合の 2 つの版(ADR 0008、D-1)の理由がここにある。
2. **`01_preliminaries_primitives.md` の §0、§2** — 計算モデルと記号。どのモジュールを書く人も読む。
3. **`../PAPER_MAPPING.md`** — 自分が担当するモジュールの行と、命名の統一(§0)、逸脱の一覧(§3)を見る。
4. 担当する節の解析文書:
   - Sec 3 のプリミティブ → `01` の §3 と §5
   - Sec 4 の評価器 → `02` の §2
   - Sec 5 の局所階層 → `02` の §3–§4
   - Sec 6–8 → `03`(数値の大きさは `thresholds.md` の §5–§7)。Sec 8.2–8.3 の実装の仕様 → `04`
5. **`open-questions.md`** — 担当するモジュールの列で絞り込み、「要合意」の行は実装の前に合意を取る。新しい曖昧点は、場所(節・補題・行番号)と候補の解釈を書いてここに足す。
6. `prior_work_and_oss.md` — 比較対象と計測条件を決めるとき(ベンチマークの担当)に読む。

## 書くときの決まり

- 論文に書かれていることと、そこから導けることと、自分の解釈を分けて書く(各文書の【補完】【解釈】【曖昧】の印)。
- 長い転載をしない。短い語句の引用と、節・補題・アルゴリズムの行番号での参照にとどめる。
- 絶対パスや個人の情報を書かない(公開する OSS)。push の前に `uv run python scripts/check_public_safety.py` を通す。
