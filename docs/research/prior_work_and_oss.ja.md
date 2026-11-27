# 先行研究と既存 OSS の調査

[English](prior_work_and_oss.md) | 日本語

- 調査日: 2026-11-27
- 対象論文: Iman Gholami, MohammadTaghi Hajiaghayi, "Linear Single-Source Shortest Path, Almost", arXiv:2611.34911v1(2026-11-25 投稿) <https://arxiv.org/abs/2611.34911>
  - 問題設定: 無向マルチグラフ、非負実数辺長、単一始点から全頂点への厳密距離と最短路木。比較・加算モデル。期待時間 O(m + n log log log n)。ほとんど確実に停止し、停止した実行は必ず正しい(論文 Theorem 1、本文 1 章)。
  - 論文の arXiv ライセンスは `arXiv.org perpetual, non-exclusive license 1.0`(abs ページのライセンスリンク <http://arxiv.org/licenses/nonexclusive-distrib/1.0/>)。このライセンスは arXiv に配布権を与えるだけなので、本文・図・翻訳を本リポジトリで再配布しない。`.cache/` は `.gitignore` 済み。
- 方針: 事実には出典を付ける。確認できなかったことは「未確認」と書く。数値は出典にあるものだけを書き、推測の数値は書かない。
- クローンしたコードは読んだだけで、実行していない。クローン先はローカルの一時ディレクトリ(`git clone --depth 1`)。

---

## 0. 本プロジェクトに直結する重要事項

1. **論文の主アルゴリズム(Algorithm 8 `ShortestPaths`)は、物理的にありうる入力では必ず基底ケースに入る。**
   - Algorithm 8 の 4 行目は「`if N ≤ L_N^64` なら、頂点対ごとに最軽量の辺を残して scan-based Dijkstra(Lemma 3.2)で解く」。ここで L_N = ⌈log2(N+2)⌉(論文テキスト 2.章の定義と Algorithm 8、`.cache/paper/2611.34911v1.txt` の 283–290 行と 2167–2185 行)。
   - 2 章の本文も「外側の構成は N > L_N^64、局所の構成は N > 2 L_N^4 を仮定する。除外されるサイズは絶対的に有界」と書いている。
   - N > L_N^64 を満たす最小の L_N を計算すると 589 になる。つまり非自明な経路に入るのは N ≳ 2^589 のときだけ(この計算は本調査で行った。Python で L = 2..2000 を走査し、2^L − 2 > L^64 となる最初の L を求めた)。正確な値は [thresholds.md](thresholds.md) にある。最小の N は 589^64 + 1 ≈ 1.94×10^177 で、[2^589 − 1, 590^64] では再び偽になり、590^64 + 1 以上ではすべて真である。
   - したがって、論文どおりの閾値で実装すると、実測できる規模では常に Dijkstra を測ることになる。層別の構成を実際に動かして測るには、閾値を縮めたパラメータ(例: 後述の kingfree 実装の `--sep k`)が必要になる。これは論文の定理の適用外なので、ベンチマークでは「論文の閾値どおり」と「閾値を縮めた変種」を区別して報告する必要がある。
2. **対象論文そのものの実装は、既に公開されている OSS が 1 件ある**(kingfree/sssp-almost-linear、2026-11-26 作成)。ただし、求値(Evaluate)段は未完成で、既定では使われていない(§1.4)。
3. 依頼された 3 件(mizchi、lcs147、mjkozicki)は、どれも**有向グラフ版**の DMM+25(arXiv:2504.17033)の実装で、対象論文とは問題設定(有向か無向か)が違う。

---

## 1. 既存 OSS

### 1.1 一覧

| リポジトリ | 実在 | ライセンス(実文で確認) | 最新コミット SHA | 日付(committer, UTC) | 主言語 | 対象論文 | 対象問題 |
|---|---|---|---|---|---|---|---|
| [mizchi/check-arxiv-2504-17033](https://github.com/mizchi/check-arxiv-2504-17033) | ○ | **LICENSE ファイル無し**。README 末尾に「## License MIT」とあるだけで、GitHub API の `license` も `null`。`Cargo.toml` に `license` フィールドは無い | `b3f173ca43d4052c0872d824be0bd61456ab7dd3` | 2026-09-23T17:08:40Z | Rust | DMM+25(arXiv:2504.17033)、C-HD(vals.ai) | 有向、非負重み(実験は整数 [1, 2^20] を f64 で保持)、単一始点・全頂点 |
| [lcs147/bmssp](https://github.com/lcs147/bmssp) | ○ | MIT(`LICENSE`: "Copyright (c) 2025 Lucas Castro and Thailsson Clementino") | `0a8eb82805bad7f64dfa1b51e17656c9af4c4681` | 2026-01-24T15:47:14Z | C++20 ヘッダオンリー(GitHub の集計上は Jupyter Notebook が最大) | DMM+25(STOC 2025)。実験の論文は Castro, Clementino, de Freitas, arXiv:2511.03007 | 有向、非負重み(テンプレートで整数・浮動小数を選べる)、単一始点・全頂点 |
| [mjkozicki/SSSP](https://github.com/mjkozicki/SSSP) | ○ | MIT(`LICENSE`: "Copyright (c) 2026 Mike Kozicki") | `30231195ef664b7ba63bd9f1190c33ecd2735575` | 2026-10-08T00:35:22Z | Python ほか 8 言語(C#, Rust, C++, Go, Java, PHP, Python, TypeScript) | DMM+25(arXiv:2504.17033v2) | 有向、非負実数(double)、単一始点・全頂点 |
| [kingfree/sssp-almost-linear](https://github.com/kingfree/sssp-almost-linear)(追加で見つけたもの) | ○ | MIT(`LICENSE`: "Copyright (c) 2026 kingfree")。README は「C コードは MIT。中文訳稿の権利は原著者に帰属」と書いている | `151ed68a9771b4211edae7700a07bbe22c5b0c38` | 2026-11-26T10:26:16Z | Typst(訳稿)、C11(実装) | **対象論文 arXiv:2611.34911v1** | 無向、非負実数(IEEE double を「不透明レジスタ」として扱う)、単一始点・全頂点 |

確認方法: `gh api repos/<o>/<r>` と `gh api repos/<o>/<r>/commits?per_page=1`、各クローンの `git rev-parse HEAD`(API の SHA と一致した)、`LICENSE` の本文。

### 1.2 mizchi/check-arxiv-2504-17033

- 出典: README <https://github.com/mizchi/check-arxiv-2504-17033/blob/b3f173ca43d4052c0872d824be0bd61456ab7dd3/README.md>、追試報告 `docs/chd-replication.md`(同 SHA)
- 実装: Dijkstra(二分ヒープ・遅延削除の `dijkstra_csr`、添字付き 4-ary ヒープの `dijkstra_dary`)、DMM+25 を論文どおりに実装し直した `src/dmm25.rs`(Algorithm 1–3、Lemma 3.3 のブロック構造、定次数化)、C-HD の `src/chd.rs`。2025 年の旧実装 V1、V2、core も残っている。
- 旧実装について、リポジトリ自身が 2026-09 に誤りを認めている。V1 は Dijkstra と同一コード、V2 の Bellman–Ford 相は到達不能、core はトップ呼び出しで完全な Dijkstra に落ちる。旧 README の「1.1〜1.6 倍速」は計測ノイズ(1 回計測、∞ の頂点を照合していなかった)。
- 正しさの検証: `tests/bmssp_family_test.rs`。n < 40 のランダム有向マルチグラフ 3000 個(重み 0、同値、自己ループ、多重辺、到達不能を含む)と、中規模のランダムグラフ、グリッド、d 正則グラフで、Dijkstra と全頂点を照合した。ベンチ全件でも照合して、不一致は 0 件と報告している(`docs/chd-replication.md` §3)。
- ベンチの条件(同 §4): 4 vCPU のクラウドコンテナ、Rust 1.94 `--release`、シングルスレッド、3 回の中央値。各頂点に出辺 d 本(終点は一様ランダム)のランダム有向グラフで、d ∈ {2, ⌈√lg n⌉, ⌊lg^{3/4} n⌋, 16}。グリッドは 4 近傍・双方向。始点は 0。計時は `Instant` で各アルゴリズムの呼び出しだけを囲んでいる(`src/bin/compare.rs` 88–96 行)。CSR の構築は計時の外。
- 結果(最速の Dijkstra に対する倍率、README): DMM+25 は 14〜121 倍遅い(グリッドでは 38〜50 倍)。C-HD は前処理込みで 1.47〜7.3 倍遅い。n = 2^14〜2^22。
- 公平性: 
  - (良い点)C-HD は前処理込みと前処理抜きを分けて報告している。比較対象の Dijkstra は 3 種類の最速値。
  - (注意)Dijkstra 用の CSR 構築は計時の外だが、DMM+25 の定次数化は計時の内側(関数が `Graph` を受け取る)。この非対称はアルゴリズムに不利な側に働いている。
  - 有向グラフでの比較なので、本プロジェクト(無向)の数値とは直接比べられない。

### 1.3 lcs147/bmssp

- 出典: README <https://github.com/lcs147/bmssp/blob/0a8eb82805bad7f64dfa1b51e17656c9af4c4681/README.md>、arXiv:2511.03007 <https://arxiv.org/abs/2511.03007>
- 実装: `include/bmssp.hpp`(`spp::bmssp`、最悪ケース版)と `include/bmssp-expected.hpp`(`spp_expected::bmssp`、期待値版で省メモリ)。`single_include/bmssp.hpp` は CI(quom)で生成する。依存はリポジトリに同梱されている(`unordered_dense.h`、`floyd_rivest_select.h` など)。
- arXiv:2511.03007 の要旨: 最悪ケース版は O(n log^{1/3} n) 空間、期待値版は O(n) 空間で同じ時間計算量。ランダム疎グラフ、グリッド、米国道路網(最大 1000 万頂点)で、Dijkstra が全条件で 3〜4 倍速い。最悪ケース版が Dijkstra に勝つには頂点数が 10^67 を大きく超える必要があると見積もっている。
- 正しさの検証: doctest(`tests/tests-wc.cpp`、`tests/tests-expected.cpp`)。グリッド(5、25、100)とランダム次数 3 グラフ(32〜4096 頂点)で、同梱の Dijkstra(`helpers/dijkstra.hpp`)と距離を照合し、さらに先行頂点から復元した経路の端点を確認している(`tests/common.cpp`)。定次数変換の有無の両方を試験する。CI(`.github/workflows/ci.yaml`)で `./tests` を実行している。
- ベンチの条件: README の表は「ランダムグラフ、辺数 3n」で、Dijkstra に対して BMSSP は 3.0〜4.2 倍遅い(n = 2^12〜2^25)。計測した機械や反復回数は README に書かれていない(未確認。arXiv:2511.03007 の本文は今回読んでいない)。`experiments/experiment.cpp` を見ると、計時は `execute(s)` の呼び出しだけを囲んでいる。`prepare_graph(false)`(グラフの準備)と Dijkstra のコンストラクタは計時の外。
- 公平性:
  - 準備(`prepare_graph`)は計時の外。ただし実験では `false`(定次数変換をしない)で呼んでおり、次数 3 のグラフなのでこの影響は小さいと思われる(未確認)。定次数変換をする場合のコストを計時に含めるかどうかは、取り込むときに決め直す必要がある。
  - 距離の型は `long long` の整数重み(`experiment.cpp` の `distT`)。実数重みの比較・加算モデルとは条件が違う。
  - 有向グラフ用。本プロジェクトで使うなら、無向グラフを両方向の弧にして入力することになる。
  - `tests/graphs/USA-road-t.NY.gr`(DIMACS の道路網)を同梱している。DIMACS データの再配布条件は未確認なので、取り込むときはデータを除くか、条件を確かめる。

### 1.4 kingfree/sssp-almost-linear(対象論文の実装)

- 出典: README <https://github.com/kingfree/sssp-almost-linear/blob/151ed68a9771b4211edae7700a07bbe22c5b0c38/README.md>、`c/README.md`(同 SHA)
- 作成は 2026-11-26、コミットは 2 件、スターは 0(GitHub API)。
- 内容: `typst/` に論文全文の中国語訳(Typst と PDF)、`c/` に C11 の実装(約 4000 行)。
- 実装の範囲(`c/README.md` の自己申告):
  - 比較・加算モデルを `len.h` で模擬し、演算回数を数える。
  - MST(3 パスの有界成長)、並列併合、LCA(論文の O(1) 構造の代わりに O(log n) の倍増法)、所有者分解、経路証明(Algorithm 2、3 相当)、所有者座標 J、外側状態グラフ、閾値成分、階層の組み立て。
- **既知の差(自己申告)**:
  1. 求値段のスロット走査(`--slot`)は、一部の入力で値の出し漏れや大きすぎる値の出力があるので、**既定では使っていない**。既定の経路は「構成した状態グラフ上の厳密探索」(`eval.c: state_graph_dist`)。これは §8.1 の構成が測地線を元の長さで実現しているかの検証であって、論文の求値手続きの実装ではない。
  2. §7.3 の幅カレンダーは未実装。J は補助グラフ上の Dijkstra で求めている。
  3. §5.4 の原子置換は未実装。
  4. O(1) の LCA と区間最小値は O(log n) で代用。
- 閾値: 既定では §0 の閾値 N ≤ L_N^64 に従うので、常に基底ケース(scan-based Dijkstra)に入る。`--sep k` で分離サイズを b = L_N^k に縮めると、階層の経路が動く。
- 正しさの検証: `--verify` でヒープ版 Dijkstra と頂点ごとに照合する。README は「n = 20〜100 のランダムグラフとパスグラフで 12/12 一致、n = 800 も一致」と報告している。`hier_validate` で階層の区間の不変量を確かめている。自動テストのスイート(CI やテストランナー)は見当たらない。
- ベンチ: `--bench` は 5 ケース(n = 500〜8100)で「生成+階層」と「求値」の時間を出力する。README に数値は無い。
- 公平性と疑問点(本調査でコードを読んで見つけたもの。実行はしていない):
  - `c/main.c` 323–336 行の `--bench` では、`name[0]=='g' && name[1]=='r'`(つまり "grid-…")のときに `gen_rand` を呼び、"rand-…" のときに `gen_grid` を呼んでいる。**ケース名とグラフの種類が入れ替わっているように読める。**
  - 既定の経路は論文の求値手続きではない(上記 1)ので、この実装の時間は論文のアルゴリズムの時間ではない。
  - README と c/README は主アルゴリズムを「算法 9」、求値を「算法 7、8」と呼んでいる。一方、arXiv v1 のテキストでは主アルゴリズムが Algorithm 8(`ShortestPaths`)、求値が Algorithm 3(`Evaluate`)になっている。どの版を参照した番号付けかは未確認。
  - 中国語訳稿は arXiv の non-exclusive ライセンスの論文の翻訳なので、再配布の許諾は未確認。本プロジェクトには取り込まない。

### 1.5 mjkozicki/SSSP

- 出典: README <https://github.com/mjkozicki/SSSP/blob/30231195ef664b7ba63bd9f1190c33ecd2735575/README.md>、`benchmark/README.md`(同 SHA)
- 実装: 8 言語それぞれに Dijkstra と `duan_mao_shu_yin` がある。README 自身が「frontier のデータ構造は簡略版で、正しさは保つ。Lemma 3.3 のブロック構造を使えば定数は改善する」と書いている。
  - Rust 版の `FrontierDS`(`rust/src/sssp.rs` 342–411 行)は、`insert` のたびに `Vec::retain` で線形に削除し、`pull` のたびに全体をソートしている。Lemma 3.3 の計算量の保証は無く、O(m log^{2/3} n) は満たさないと読める(本調査の読解。計測はしていない)。
- 正しさの検証: 各言語に Dijkstra との照合テストがある。例えば Rust の `matches_dijkstra_small` と `matches_dijkstra_parallel_paths`(許容誤差 1e-10)、Python の `test_matches_dijkstra`。ただし小さいグラフが中心。
- ベンチの条件(`benchmark/README.md`): 50,000 頂点、約 20 地域のネットワークを模した有向グラフ(地域内の出辺 2〜8 本で重み 0.5〜5、地域間の出辺 0〜2 本で重み 20〜150、seed = 42)。各言語の Docker コンテナを起動し、`wall_sec` は **`docker run` から `docker wait` までの時間**(`benchmark/run_benchmarks.py` 209–233 行)。既定は 1 回計測。
- 結果(`benchmark/results.json`): ほとんどの言語で DMSY と Dijkstra がどちらも約 1.9〜2.0 秒で、`peak_mem_mb` と `total_cpu_sec` が 0.0 になっている。
- 公平性:
  - 計測値にコンテナの起動とグラフファイルの読み込みが含まれる。アルゴリズムの差がほとんど見えないので、アルゴリズムの比較には使えない。
  - 1 回計測。メモリと CPU の採取が多くの言語で 0 になっていて、計測が機能していないように見える。
  - `rust/src/bin/compare.rs` には反復の中央値を取る別の比較もあるが、その結果はリポジトリに無い(未確認)。

---

## 2. 先行研究(2023〜2026 年を中心に)

「一致」は本プロジェクトの問題設定(無向、非負実数、単一始点・全頂点、比較・加算モデル)と一致するかどうか。

| 略号 | 文献 | 問題設定 | 計算量 | モデル | 一致 | 比較対象にするか |
|---|---|---|---|---|---|---|
| Dij59 / FT87 | Dijkstra 1959、Fredman–Tarjan 1987(Fibonacci ヒープ)。DOI 10.1007/BF01386390、10.1145/28869.28874 | 有向・無向、非負実数 | O(m + n log n)。二分ヒープでは O((m+n) log n) | 比較・加算 | ○ | **必須(基準)**。自前で実装する(二分ヒープ、4-ary ヒープ、可能なら pairing ヒープ) |
| Tho99 | Thorup, "Undirected Single-Source Shortest Paths with Positive Integer Weights in Linear Time", J. ACM 46(3), 1999。DOI 10.1145/316542.316548 | 無向、**正の整数**重み | O(m)(最悪) | word RAM(重みのビット表現を使う) | △(重みとモデルが違う) | 任意。整数重みの部分集合で参考として測るだけ。実数重みの比較には入れない。公開実装は未確認 |
| PR05 | Pettie, Ramachandran, "A Shortest Path Algorithm for Real-Weighted Undirected Graphs", SIAM J. Comput. 34(6), 2005。DOI 10.1137/S0097539702419650 | 無向、非負実数 | 前処理 O(MST(m,n) + min{n log n, n log log r})、r は重みの比 | 比較・加算 | ○ | 理論上の比較対象。公開実装は未確認。自前で実装するのは大がかり |
| DMS+23 | Duan, Mao, Shu, Yin, "A Randomized Algorithm for Single-Source Shortest Path on Undirected Real-Weighted Graphs", FOCS 2023。arXiv:2307.04139 <https://arxiv.org/abs/2307.04139> | 無向、非負実数 | 乱択、O(m √(log n · log log n)) | 比較・加算 | ○ | **有力**。問題設定が完全に一致する。要旨は「構造が単純で実装しやすい」と述べている。公開実装は今回の検索(`gh search repos`)では見つからなかった |
| Yan26 | Shuyi Yan, "Lossless Derandomization for Undirected Single-Source Shortest Paths and Approximate Distance Oracles", SOSA 2026。arXiv:2510.12598 <https://arxiv.org/abs/2510.12598> | 無向、非負実数 | 決定的、O(m √(log n log log n)) | 比較・加算 | ○ | DMS+23 の決定的版。DMS+23 を実装するなら、余力があれば追加 |
| DMM+25 | Duan, Mao, Mao, Shu, Yin, "Breaking the Sorting Barrier for Directed Single-Source Shortest Paths", STOC 2025。arXiv:2504.17033 <https://arxiv.org/abs/2504.17033> | **有向**、非負実数 | 決定的、O(m log^{2/3} n) | 比較・加算 | △(有向) | 副次的な比較対象。無向グラフは両方向の弧にすれば解けるが、無向性は使わない。lcs147/bmssp を使える(§3) |
| DMS+26 | Duan, Mao, Shu, Yin, "A Faster Directed Single-Source Shortest Path Algorithm", ICALP 2026(LIPIcs 374, 81:1–81:23)。arXiv:2602.07868 <https://arxiv.org/abs/2602.07868> | 有向、非負実数 | 決定的、O(m √(log n) + √(mn log n log log n))。疎グラフでは O(m √(log n log log n)) | 比較・加算 | △(有向) | 参考。公開実装は未確認 |
| GZ26 | Grishin, Zabokritskiy, "Improved Bounds for Directed SSSP in Sparse and Intermediate-Density Graphs", arXiv:2611.26453 <https://arxiv.org/abs/2611.26453>(2026-11-19 投稿) | 有向、非負実数 | 決定的、m = O(n) で O(n √(log n) (log log n)^{1/4}) | 比較・加算 | △(有向) | 参考。プレプリント |
| KR26 | Kadria, Roditty, "A Faster Undirected Single-Source Shortest Path Algorithm", arXiv:2609.15247 <https://arxiv.org/abs/2609.15247>(2026-09-14 投稿) | 無向、非負実数 | O(m log^{1/2} n · (log log n)^{1/4} · (log log log n)^{1/4})。対象論文は「高確率で」と紹介している | 比較・加算 | ○ | 理論上の直接の競合(対象論文以前の無向の最良)。公開実装は未確認。プレプリント |
| HHR+24 | Haeupler, Hladík, Rozhoň, Tarjan, Tětek, "Universal Optimality of Dijkstra via Beyond-Worst-Case Heaps", FOCS 2024(pp. 2099–2130、Princeton の書誌 <https://collaborate.princeton.edu/en/publications/universal-optimality-of-dijkstra-via-beyond-worst-case-heaps/>)。arXiv:2311.11793 <https://arxiv.org/abs/2311.11793> | 有向・無向、非負実数。ただし主に「距離順序問題」(頂点を距離順に並べる) | working-set 性を持つヒープを使った Dijkstra が universally optimal | 比較・加算 | △(距離の値だけでなく順序も出す問題) | 解釈上の重要文献。対象論文の 1 章の冒頭(距離を計算するのと並べるのは別)と関係が深い。working-set ヒープの公開実装は未確認。FOCS 2024 の最優秀論文賞は、二次資料(arXiv:2504.17327)でしか確認できなかった |
| Castro+25 | Castro, Clementino, de Freitas, "Implementation and Brief Experimental Analysis of the Duan et al. (2025) Algorithm", arXiv:2511.03007 <https://arxiv.org/abs/2511.03007> | DMM+25 の実装と実験 | — | — | — | ベンチの方法の参考(道路網、グリッド、ランダム) |

補足:

- 対象論文の 1.1 節(Related Work)が挙げる文献は、上の表の Dij59、FT87、Tho99、PR05、DMS+23、Yan26、DMM+25、DMS+26、GZ26、KR26 と、FRT04(Fakcharoenphol, Rao, Talwar の木埋め込み。DOI 10.1016/J.JCSS.2004.04.011)。FRT04 は部品(経路推定器)なので、比較対象ではない。HHR+24 は対象論文の参考文献には無い(依頼により追加した)。
- 対象論文のテキストは、上界を「連結な無向入力、または有向入力の到達可能部分」について書くと断っている。到達可能な頂点の特定と n 個の出力ラベルの書き込みで、高々 O(m + n) が加わる(1.1 節)。比較するときは、各手法で到達不能な頂点の扱い(∞ を書く)をそろえる必要がある。
- 対象論文の時間計算量は**期待値**(乱択)。DMS+23 も乱択で、Yan26、DMM+25、DMS+26、GZ26 は決定的。乱択の手法は、シードを複数にして分布で報告する必要がある。

---

## 3. 既存 OSS を取り込むか、独自に再実装するか

### 3.1 推奨

**対象論文の実装と、基準の Dijkstra は独自に実装する。既存 OSS のコードは取り込まない。** 唯一の例外として、有向版 BMSSP の副次比較が要るときだけ lcs147/bmssp を外部依存として固定する。

| リポジトリ | 推奨 | 理由 |
|---|---|---|
| kingfree/sssp-almost-linear | **取り込まない**。読むだけの相互参照にする | MIT なので法的には取り込める。しかし、求値段が未完成で既定では使われない、ベンチのケース名とグラフが入れ替わっているように読める、自動テストが無い、作成から 1 日でコミットが 2 件。論文の正しい再現の基準にはできない。訳稿は再配布の許諾が未確認。本プロジェクトの実装が完成したら、同じ入力で距離を照合する「独立実装どうしの突き合わせ」に使う価値はある(実行は別途判断) |
| lcs147/bmssp | **副次比較が必要な場合だけ**、外部依存として SHA で固定する | MIT で、ヘッダオンリー C++20、Dijkstra との照合テストと CI がある。本プロジェクトは C++/CMake なので組み込みやすい。ただし有向・整数重みでの実験が前提なので、無向グラフは両方向の弧にし、重みの型を `double` にしたときの動作を自前のテストで確かめてから使う |
| mizchi/check-arxiv-2504-17033 | **取り込まない**。数値と方法を文献として引用するだけにする | LICENSE ファイルが無い(README に「MIT」と書いてあるだけ)。MIT は著作権表示と許諾文の同梱を条件にしているが、その文面が無い。Rust なので本プロジェクトの言語とも違う。比較回数の数え方(`count-cmp`)やベンチの方法は参考になる |
| mjkozicki/SSSP | **取り込まない** | frontier が簡略版で、DMM+25 の計算量を満たさない。ベンチはコンテナの起動込みの 1 回計測で、アルゴリズムの比較になっていない |

### 3.2 取り込むときのコミットの固定方法(lcs147/bmssp の場合)

- CMake の `FetchContent_Declare(bmssp GIT_REPOSITORY https://github.com/lcs147/bmssp.git GIT_TAG 0a8eb82805bad7f64dfa1b51e17656c9af4c4681)` のように、**タグやブランチではなく 40 桁の SHA を書く**。git submodule でもよい(SHA が `.gitmodules` とツリーに記録される)。
- 使うのは `include/`(`external/` を含む)だけにして、`tests/graphs/` のデータ(DIMACS 由来を含む)と `experiments/` は使わない。
- MIT の条件として、配布物にその `LICENSE` を同梱し、`THIRD_PARTY_NOTICES` などに名前、SHA、ライセンスを記録する。同梱の `include/external/*`(`unordered_dense.h`、`floyd_rivest_select.h`、`median_of_ninthers.h` など)のライセンスは、取り込む前に各ファイルの見出しで確かめる(本調査では未確認)。
- 本プロジェクト自体のライセンスはリポジトリにまだ無い(`LICENSE` が見当たらない)。MIT のコードと組み合わせられるかは、本プロジェクトのライセンスが決まってから最終確認する(MIT は一般に他のほとんどのライセンスと組み合わせられる)。

### 3.3 比較対象の優先順位(提案)

1. Dijkstra(自前実装。二分ヒープ、4-ary ヒープ。比較回数も数える)。対象論文の閾値どおりの実装は §0 のとおり Dijkstra の基底ケースに入るので、これが事実上の主な比較になる。
2. 対象論文の実装(自前)。「論文の閾値どおり」と「閾値を縮めた変種」を分けて報告する。
3. DMS+23(自前で実装する場合)。問題設定が完全に一致する唯一の、実装しやすい先行の非 Dijkstra 手法。
4. DMM+25(lcs147/bmssp を固定して使う)。有向の手法を無向入力に適用した参考値と明記する。

---

## 4. 未確認事項

- arXiv:2511.03007 の本文にある計測環境と反復回数(要旨だけ読んだ)。
- lcs147/bmssp の README の表を計測した機械と条件。
- lcs147/bmssp に同梱された `include/external/*` のライセンス。
- DIMACS 道路網データの再配布条件。
- Tho99、PR05、DMS+23、Yan26、DMS+26、KR26、GZ26、HHR+24 の公開実装の有無(`gh search repos` で論文番号と題名を検索して見つからなかったが、網羅的ではない)。
- kingfree 実装の `--bench` のケース名の入れ替わりが意図的かどうか(コードを読んだだけで、実行していない)。
- kingfree 実装のアルゴリズム番号が、どの版の論文に対応するか。
- HHR+24 の FOCS 2024 最優秀論文賞(二次資料でしか確認できなかった)。
