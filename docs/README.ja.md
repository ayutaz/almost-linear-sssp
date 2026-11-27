# 文書の案内

[English](README.md) | **日本語** | [简体中文](README.zh-CN.md)

このページでは、`docs/` の下に何があり、それぞれ誰のためのもので、どの順に読むとよいかを説明します。プロジェクトの概要は[トップの README](../README.ja.md) を見てください。

一部の文書は、いまのところ日本語版しかありません。両方の言語がある文書では、日本語版は英語版と同じ場所にある `.ja.md` のファイルです。ファイル名、節の番号、コマンド、識別子はどの言語でも同じです。

## どこから読むか

| やりたいこと | 読む順番 |
|---|---|
| ここまでに何が分かったかを知る | [../README.ja.md](../README.ja.md) → [research/thresholds.md](research/thresholds.md) §2–§3 → [ADR 0003](adr/0003-three-implementation-modes.md) → [ADR 0008](adr/0008-mode-a-outer-base-case-variants.md) |
| 何が終わり、次に何をするかを知る | [PROGRESS.md](PROGRESS.md) → [TASKS.md](TASKS.md) |
| 論文の一部を実装する | [research/README.md](research/README.md) → 担当する節のノート → [PAPER_MAPPING.md](PAPER_MAPPING.md) → [research/open-questions.md](research/open-questions.md) → [tdd/README.md](tdd/README.md) |
| グラフ、ベースライン、検証に手を入れる | [design/graph_library.md](design/graph_library.md)、[design/baselines.md](design/baselines.md)、[design/numerics.md](design/numerics.md) |
| ベンチマークを実行する・確かめる | [design/benchmark_protocol.md](design/benchmark_protocol.md) → [design/benchmark_schema.md](design/benchmark_schema.md) → [design/benchmark_runner_cli.md](design/benchmark_runner_cli.md) → [research/prior_work_and_oss.md](research/prior_work_and_oss.md) |
| CI やツールを変える | [design/ci.md](design/ci.md)、続けてそこから参照されている ADR |
| コードで貢献する | [../CONTRIBUTING.ja.md](../CONTRIBUTING.ja.md) → [tdd/README.md](tdd/README.md) → [TASKS.md](TASKS.md) |

## それぞれの中身

### `research/` — 論文を読む

論文(arXiv 2611.34911v1)を節ごとに、自分たちの言葉で解説したノートです。定義、補題、アルゴリズム、定数、それらの依存関係、実装で気をつけるべき点を書いています。論文そのものはリポジトリに置いていません。手元で読むには `make fetch-paper` を実行します。`l.284` のような行番号は、このスクリプトが取り出したテキストの行番号です。

- **読む人:** 論文の構成要素を実装・レビューする人。論文についての私たちの主張を確かめたい人。
- **索引と読む順番:** [research/README.md](research/README.md)。
- **最も重要なノート:** [research/thresholds.md](research/thresholds.md)。論文の構成が動き始める正確なサイズと、どの定数が正しさに効き、どの定数が実行時間だけに効くか。
- **先行研究:** [research/prior_work_and_oss.md](research/prior_work_and_oss.md)。関連するアルゴリズムと、既存のオープンソース実装(ライセンスつき)。

### `PAPER_MAPPING.md` — 論文とコードの対応

論文の定義・補題・アルゴリズム・主要な定数ごとに 1 行ずつ、実装するモジュール、テスト、モード(A か B)、状態を書いた表です。論文からの**逸脱**(`D-xx`)もすべてここに載せます。`src/paper/` や `src/experimental/paper_experimental/` に手を入れる前に読み、そこを変えるたびに更新します。

### `design/` — コードの作り

論文から取ったものではない部分の設計文書です。

| 文書 | 内容 |
|---|---|
| [graph_library.md](design/graph_library.md) | グラフの表現、テキスト形式、再現可能な生成器、`sssp_gen` |
| [baselines.md](design/baselines.md) | Dijkstra の各版とヒープ、同点の扱い、演算のカウンタ |
| [numerics.md](design/numerics.md) | 比較・加算モデルと IEEE-754 の `double` の違い。`double` の結果に使う許容誤差とその根拠 |
| [benchmark_protocol.md](design/benchmark_protocol.md) | ベンチマークの実行方法: ビルド、反復、計時の範囲、検証、統計、計測の組 |
| [benchmark_schema.md](design/benchmark_schema.md) | ベンチマーク結果のレコードの形式 |
| [benchmark_runner_cli.md](design/benchmark_runner_cli.md) | C++ のベンチマーク用ランナーのコマンドライン |
| [ci.md](design/ci.md) | GitHub Actions のワークフローと、環境構築のスクリプト |

- **読む人:** その部分を変える貢献者と、レビューする人。

### `adr/` — 決定の記録

Architecture Decision Record(設計上の決定の記録)です。今後の作業を縛る決定ごとに短いファイルを 1 つ作り、背景と帰結を書きます。決定は消さずに、新しい ADR で置き換えます。索引は [adr/README.md](adr/README.md) です。まず [0003](adr/0003-three-implementation-modes.md)(3 つのモード)と [0008](adr/0008-mode-a-outer-base-case-variants.md)(論文に忠実な基本ケース)を読んでください。

- **読む人:** 変更を提案する前に「なぜこうなっているのか」を知りたい人。

### `TASKS.md` と `PROGRESS.md` — 進捗の管理

- [TASKS.md](TASKS.md) は作業項目の唯一の一覧です。ID、優先度、状態、依存関係、受け入れ条件、テストの証跡、完了したコミットを書きます。
- [PROGRESS.md](PROGRESS.md) は短い現況です。いまの段階、実装済みのもの、最新のテストとベンチマークの状態、既知の問題、次のタスクを書きます。作業を再開するときに最初に読みます。

どちらも取りまとめ役のメンテナが管理します。ほかの貢献者は、並行して直接編集せずに、変更を提案してください。

### `tdd/` と `TDD_LOG.md` — テストを先に書いた証跡

機能の追加やバグの修正には、それぞれタスクのファイル `tdd/tasks/<ID>.json` があります。失敗したテストの実行(RED)、通った実行(GREEN)、リファクタリング後の実行を、コマンド、終了コード、出力とともに記録します。形式は [tdd/README.md](tdd/README.md) にあります。`TDD_LOG.md` はこれらのファイルから `scripts/tdd.py log` で生成するログで、手で編集しません(初めて生成するまでは存在しません)。

- **読む人:** 変更が本当にテストを先に書いて作られたかを確かめるレビュー担当。

### `reports/` — 結果

保存された計測から生成した、ベンチマークとレビューの報告を置きます。いまは空です。**まだ何も計測していません。** ここに載る数値は、すべて `results/` に保存された計測から取ったものでなければなりません。

## 決まりごと

- 論文は節・補題・アルゴリズム・行番号で参照し、長く引用しません。
- 計測した数値と推定した数値を混ぜません。計測していないものは「未計測」と書きます。
- 絶対パスや個人の情報を書きません(公開リポジトリのため)。変更を公開する前に `make public-safety` を実行します。
