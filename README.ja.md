# almost-linear-sssp

[English](README.md) | **日本語** | [简体中文](README.zh-CN.md)

Iman Gholami と MohammadTaghi Hajiaghayi の論文 **「Linear Single-Source Shortest Path, Almost」**
([arXiv:2611.34911v1](https://arxiv.org/abs/2611.34911))を、独立に再現・検証し、性能を測るプロジェクトです。

論文のアルゴリズムを C++20 で実装し、すべての結果を独立した参照実装で確かめます。そのうえで、よく調整した Dijkstra 法と比べて性能を測る予定です。著者とは関係のない、独立したプロジェクトです。

> **状態: 作業の初期段階です。** グラフのライブラリ、Dijkstra 系のベースライン、検証の道具はできています。論文そのものの構成要素は実装中です。
> **まだ何も計測していません。** そのため、この README には性能の数値を一切載せていません。

## 目次

- [現状](#現状)
- [論文の要旨](#論文の要旨)
- [閾値がなぜ重要か](#閾値がなぜ重要か)
- [はじめかた](#はじめかた)
- [アルゴリズム](#アルゴリズム)
- [正しさの確かめかた](#正しさの確かめかた)
- [ベンチマーク](#ベンチマーク)
- [リポジトリの構成](#リポジトリの構成)
- [文書](#文書)
- [関連研究](#関連研究)
- [引用](#引用)
- [ライセンス](#ライセンス)
- [免責](#免責)

## 現状

| 項目 | 状態 |
|---|---|
| グラフのライブラリ: CSR、隣接リスト、辺リスト。辺リスト形式と DIMACS `.gr` の入出力 | 完了 |
| 再現可能なグラフ生成器と、コマンドラインツール `sssp_gen` | 完了 |
| ベースライン: Bellman–Ford。二分ヒープ、4-ary ヒープ、Fibonacci ヒープ、radix ヒープの Dijkstra | 完了 |
| 検証: 最短路木の検証器、厳密な有理数の参照実装、演算回数を数える `Counted<W>` | 完了 |
| 論文の構成要素、Sec 2〜5(基本部品、階層の評価、局所階層) | 実装中 |
| 論文の構成要素、Sec 6〜8(ピボットの標本抽出、owner 座標、外側の階層) | 予定 |
| 差分テストの仕組みとベンチマークのランナー | 実装中 |
| ベンチマークの結果 | **未計測** |

作業項目の一覧は [docs/TASKS.md](docs/TASKS.md)、現在の要約は [docs/PROGRESS.md](docs/PROGRESS.md) にあります。

## 論文の要旨

`n` 頂点・`m` 辺で、辺の長さが非負の実数の無向グラフが与えられたとき、1 つの始点から他のすべての頂点までの距離を求めるのが単一始点最短路(SSSP)問題です。Fibonacci ヒープを使う Dijkstra 法は `O(m + n log n)` 時間かかります。この論文は、*比較・加算モデル*(辺の長さには比較と加算しか許されず、丸めやビット演算は使えない計算モデル)で、期待実行時間 **`O(m + n log log log n)`** を主張しています。`m` が `n log log log n` 以上なら線形時間です。アルゴリズムは乱択ですが *Las Vegas 型* で、停止した実行は必ず正確な距離を返し、乱数が影響するのは実行時間だけです。主な道具は、最小全域木をもとに頂点を領域に分けて各領域に 1 つの「ピボット」を置く分解、ピボット間の短い経路を乱択の標本抽出で保証する仕組み、頂点ごとの補助座標、そして強連結成分の順序付き階層に沿って距離を評価する手続きです。

## 閾値がなぜ重要か

これがここまでで最も重要な発見で、プロジェクト全体の進め方を決めています。

論文は、入力が十分に大きいときだけ本来の構成を使います。`N` 頂点の入力に対して `L_N = ⌈log2(N + 2)⌉` とおくと(Sec 2)、条件は次のとおりです。

- **外側の構成:** `N > L_N^64`。これより小さい入力は、Dijkstra 法を実行する基本ケースに回ります。
- **局所の構成**(小さな部分に対して使う): `N > 2 · L_N^4`。これより小さい部分は、Kruskal 法に似た単純な基本ケースに回ります。

多倍長整数で境界を厳密に計算しました([docs/research/thresholds.md](docs/research/thresholds.md))。

| 条件 | 初めて成り立つ `N` |
|---|---|
| 外側 `N > L_N^64` | `589^64 + 1` ≈ 1.9 × 10^177(そのすぐ上の短い区間で再び偽になる。`N ≥ 590^64 + 1` ≈ 2.2 × 10^177 からはすべて成り立つ) |
| 局所 `N > 2 · L_N^4` | 209,953(`2 · 18^4 = 209,952`) |

| `N` | `L_N` | 外側の構成は動くか | 局所の条件は成り立つか |
|---|---|---|---|
| 10^3 | 10 | いいえ(`L_N^64` ≈ 2^213) | いいえ(`N > 20,000` が必要) |
| 10^5 | 17 | いいえ(≈ 2^262) | いいえ(`N > 167,042` が必要) |
| 10^6 | 20 | いいえ(≈ 2^277) | はい(`N > 320,000`) |
| 10^8 | 27 | いいえ(≈ 2^304) | はい(`N > 1,062,882`) |
| `N < 2^64` のすべて | ≤ 65 | **いいえ**。`L_N^64 ≥ 2^64 > N` のため | — |

つまり、論文をそのまま実装すると、**メモリに載るどの入力でも基本ケース(Dijkstra 法)しか実行されません**。局所の構成は外側の構成の中からしか呼ばれないので、これも動きません。論文どおりの実装の時間を測っても、それは Dijkstra 法といくつかの前処理の時間であって、論文の新しい考え方の時間ではありません。

構成そのものを誠実に調べるために、このプロジェクトでは 3 つのモードを分けています([ADR 0003](docs/adr/0003-three-implementation-modes.md))。

| モード | ディレクトリ | 中身 | 目的 |
|---|---|---|---|
| **A** paper-faithful(論文に忠実) | `src/paper/` | 閾値は論文に書かれたとおり。桁あふれしない整数演算で判定する。基本ケースに入るたびに記録する。 | 論文を書かれたとおりに再現する。現実の入力では基本ケースになる。変種 `paper-faithful-d1` は、基本ケースの走査型 Dijkstra を二分ヒープ版に替える(逸脱 D-1、[ADR 0008](docs/adr/0008-mode-a-outer-base-case-variants.md))。 |
| **B** paper-experimental(実験版) | `src/experimental/paper_experimental/` | 論文の構造のまま、サイズの閾値を小さく、設定できるようにしたもの。正しさの証明が依存する定数は保つか、新しいサイズから導き直す([thresholds.md §6](docs/research/thresholds.md))。 | Sec 5〜8 を現実的な規模で実際に動かし、すべての結果を検証する。正しさを確かめられない入力は失敗として扱う。 |
| **C** practical-optimized(実用最適化) | `src/experimental/practical/` | 正しく、かつ明示されていれば、どんな最適化でもよい。 | 実用上の速さ。 |

「論文のアルゴリズム」と呼ぶのはモード A だけです。報告する時間にはすべて、基本ケースで処理した頂点と辺の割合を並べて示します。

## はじめかた

必要なもの: macOS または Linux、C++20 のコンパイラ(開発とテストは Apple Clang と LLVM Clang で行っています)、CMake ≥ 3.24、Ninja、Python ≥ 3.12 を管理する [uv](https://docs.astral.sh/uv/)。CI と貢献者には clang-format と clang-tidy も必要です。

```bash
make bootstrap                 # ツールを確かめ、足りないものを入れ、`uv sync` を実行
                               # (./scripts/bootstrap.sh --check は報告だけで、何も入れない)
make build                     # debug プリセットを build/debug/ に設定・ビルド
make test                      # C++ のテスト(ctest)と Python のテスト(pytest)を実行

./build/debug/sssp_gen --list  # グラフ生成器と重み分布の一覧
./build/debug/sssp_gen --type er --n 1000 --avg-degree 4 --weights uniform_int --seed 1 --out g.txt

make fetch-paper               # 論文を .cache/paper/ に取得(コミットはしない)
```

同じシードからは、どのプラットフォームでも同じグラフが生成されます。ほかのビルドプリセット(`release`、`sanitizer`、`coverage`、`ci` など)は [CMakePresets.json](CMakePresets.json) にあります。たとえば `make build PRESET=release` です。貢献する前には `make verify` を通してください([CONTRIBUTING.ja.md](CONTRIBUTING.ja.md))。

## アルゴリズム

| アルゴリズム | 関数・場所 | 重み | 状態 |
|---|---|---|---|
| Bellman–Ford(ラウンド単位。参照実装として使う) | `sssp::baseline::bellman_ford` | `uint64_t`, `double` | 完了 |
| Dijkstra、二分ヒープ(遅延削除) | `sssp::baseline::dijkstra_binary_heap` | `uint64_t`, `double` | 完了 |
| Dijkstra、4-ary ヒープ(decrease-key) | `sssp::baseline::dijkstra_4ary_heap` | `uint64_t`, `double` | 完了 |
| Dijkstra、Fibonacci ヒープ | `sssp::baseline::dijkstra_fibonacci_heap` | `uint64_t`, `double` | 完了 |
| Dijkstra、radix ヒープ | `sssp::baseline::dijkstra_radix_heap` | `uint64_t` のみ | 完了 |
| 論文 Sec 3: 基本部品(MST、順序付き union–find、LCA / level ancestor、split–findmin) | `src/paper/`(モード A) | — | 実装中 |
| 論文 Sec 4: 順序付き階層の評価 | `src/paper/`(モード A) | — | 実装中 |
| 論文 Sec 5: 局所階層 | `src/paper/`(モード A) | — | 実装中 |
| 論文 Sec 6〜8: ピボット経路、owner 座標、外側の階層、主アルゴリズム | `src/paper/`(モード A) | — | 予定 |
| 閾値の調整・実用的な最適化 | モード B / モード C | — | 予定 |

ベースラインはすべて 1 つの走査ループを共有し、このリポジトリで書いたヒープを使い、同じ演算を数えます([docs/design/baselines.md](docs/design/baselines.md))。論文の各補題・各アルゴリズムとコード・テストの対応は [docs/PAPER_MAPPING.md](docs/PAPER_MAPPING.md) にあります。

## 正しさの確かめかた

速さより正しさを優先します。どの実装も、その実装とコードを共有しない道具で確かめます。

- **独立した参照実装。** Bellman–Ford を使います。Bellman–Ford 自身も、頂点数の少ない重み付き単純グラフを*すべて*列挙し、すべての始点で Floyd–Warshall と照合しています。
- **最短路木の検証器**(`sssp::verify_sssp`)。`dist[source] = 0` であること、距離をまだ縮められる辺が無いこと、木の辺がすべてぴったり(tight)であること、親をたどると始点に戻ること、到達不能な頂点がちょうど無限大になっていることを確かめます。これにより、アルゴリズムを信用せずに結果の正しさを証明できます。
- **厳密な有理数の参照実装。** 浮動小数の入力に対しては、厳密な分数で計算する Python の参照実装(`tests/python/reference/`)と照合します。
- **厳密な比較と、許容誤差つきの比較。** 整数の重みでは完全一致を求めます。`double` の結果は、[docs/design/numerics.md](docs/design/numerics.md) で導いた丸め誤差の上界で比べます。テストを通すために許容誤差を緩めることはせず、どのビルドでも fast-math を無効にしています。
- **演算回数の計数。** `Counted<W>` は重みの型を包んで比較と加算の回数を数え、比較・加算モデルに従っているかを確かめられるようにします。
- **回帰テスト。** 見つけた反例はすべて最小化して `tests/regression/` に残します。
- **黙った代替をしない。** 論文の構成の代わりに直接の探索で入力の一部を処理したときは、必ず `Diagnostics.fallbacks` に記録します。

コードはテストを先に書き、失敗したテストの実行と通ったテストの実行を、タスクごとに `docs/tdd/tasks/` に記録しています。

## ベンチマーク

**結果はまだありません。保存された計測で裏づけられるまで、性能の数値はすべて「未計測」です。**

手順は前もって決めてあります([docs/design/benchmark_protocol.md](docs/design/benchmark_protocol.md))。診断機能を外した Release ビルドを使うこと。すべてのアルゴリズムで、同じ CSR グラフ、同じ入力、同じマシン、同じウォームアップと反復回数にすること。グラフの構築は計時に含めないこと。結果は検証を通ったものだけを数えること。環境の情報を結果と一緒に保存すること。中央値、平均、p95、標準偏差、標本数を報告すること。速度比は二分ヒープの Dijkstra を基準にすること。論文のアルゴリズムの時間には、基本ケースの割合を並べて示すこと。

計測は `make benchmark-smoke`、`make benchmark`、`make benchmark-full` で行い、`make report` で集計する予定です。C++ のベンチマーク用ランナーをまだ書いているところなので、これらのターゲットはまだ結果を出しません。

## リポジトリの構成

| パス | 中身 |
|---|---|
| `include/sssp/` | 公開ヘッダ: グラフの型(`CsrGraph`)、結果(`SsspResult`、`Diagnostics`)、ベースライン、生成器、入出力、検証 |
| `src/baseline/` | Bellman–Ford と Dijkstra の各版、およびそのヒープ |
| `src/graph/`, `src/verify/` | グラフの入出力と生成器、最短路の検証 |
| `src/paper/` | モード A、論文に忠実な実装(実装中) |
| `src/experimental/` | モード B(`paper_experimental/`)とモード C(`practical/`) |
| `tools/` | コマンドラインツール(`sssp_gen`) |
| `tests/` | C++ のテスト(Catch2: unit、integration、differential、property、regression)と Python のテスト(pytest) |
| `benchmarks/`, `scripts/` | ベンチマークのドライバとスクリプト(環境構築、論文の取得、ベンチマークの実行、集計) |
| `docs/` | 研究ノート、設計文書、決定の記録、タスクの管理 |

## 文書

まず [docs/README.ja.md](docs/README.ja.md) を読んでください。`docs/` の下の各文書が何で、誰のためのもので、どの順に読むとよいかを説明しています。一部の文書は、いまのところ日本語版しかありません。両方の言語がある文書では、日本語版は英語版と同じ場所にある `.ja.md` のファイルです。

## 関連研究

この問題設定(無向、非負の実数重み、比較・加算モデル)でのこれまでの最良の上界には、Duan, Mao, Shu, Yin(FOCS 2023、乱択で `O(m √(log n · log log n))`)、その脱乱択化である Yan(SOSA 2026)、Kadria と Roditty の 2026 年のプレプリントがあります。古典的な比較の基準は、Fibonacci ヒープの Dijkstra 法、整数重みに対する Thorup の線形時間アルゴリズム(1999)、Pettie と Ramachandran(2005)です。同じ論文の公開実装はほかにも少なくとも 1 つあります。[docs/research/prior_work_and_oss.md](docs/research/prior_work_and_oss.md) に、これらの研究と既存のオープンソースのコードを、ライセンスと比較可能性の注記つきでまとめています。

## 引用

このソフトウェアを使う場合は、元の論文と、必要に応じてこのソフトウェアを引用してください。

```bibtex
@misc{gholami2026linear,
  title         = {Linear Single-Source Shortest Path, Almost},
  author        = {Gholami, Iman and Hajiaghayi, MohammadTaghi},
  year          = {2026},
  eprint        = {2611.34911},
  archivePrefix = {arXiv},
  note          = {Version 1}
}
```

ソフトウェア自体の引用情報は [CITATION.cff](CITATION.cff) にあります。

## ライセンス

このリポジトリのコードと文書は [MIT License](LICENSE) で公開しています。

論文は**含めていません**。arXiv は再配布を認めないライセンスで論文を配布しているため、`make fetch-paper` で git の管理外の `.cache/paper/` に取得します([ADR 0005](docs/adr/0005-paper-not-redistributed.md))。このリポジトリの文書は、論文を引用・転載する代わりに、節・補題・アルゴリズムの番号で参照しています。

## 免責

これは独立した再現です。論文の著者とは関係がなく、著者の承認も受けていません。実装や論文の読み方の誤りは、すべて私たちの責任です。ここでの結果が論文と食い違うときは、まずこちらの再現を疑います。そうした問題や論文からの逸脱は、[docs/PAPER_MAPPING.md](docs/PAPER_MAPPING.md) と [docs/adr/](docs/adr/) で公開して記録します。
