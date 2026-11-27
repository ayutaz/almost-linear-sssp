# 貢献の手引き

[English](CONTRIBUTING.md) | **日本語** | [简体中文](CONTRIBUTING.zh-CN.md)

このプロジェクトは arXiv 2611.34911v1(「Linear Single-Source Shortest Path, Almost」)を再現します。速さよりも、正しさと誠実な報告を優先します。拘束力のある規則は `CLAUDE.md` にあり、このファイルはそれを人間向けにまとめたものです。文書の全体像は [docs/README.ja.md](docs/README.ja.md) を見てください。

## 準備

```bash
./scripts/bootstrap.sh          # ツールの確認 + uv sync + 設定
make build && make test
./scripts/fetch_paper.sh        # 論文を .cache/paper/ に取得(コミットしないこと)
```

必要なもの:

- C++20 に対応した Clang(Apple Clang または LLVM)、CMake ≥ 3.24、Ninja
- uv、Python ≥ 3.12
- clang-format、clang-tidy
- `gh`(CI の確認に使う)

Python の依存関係は `uv add <pkg>`(開発用は `uv add --dev <pkg>`)で追加し、`uv run` で実行します。pip は使わないでください。

## 作業の流れ

1. `docs/TASKS.md` から、依存するタスクがすべて DONE になっているタスクを選ぶ。
2. 証跡を残しながら TDD で進める(`.claude/skills/tdd/SKILL.md`):
   `scripts/tdd.py start` → 失敗するテスト → `red` → 実装 → `green` → `refactor`。
3. 論文のコードでは、`docs/PAPER_MAPPING.md` を最新に保ち、逸脱を記録する(`.claude/skills/paper-reproduction/SKILL.md`)。
4. 確認: `make verify`(アルゴリズムを変えたときは `make sanitizer` も)。
5. Conventional Commits の形式で、タスク ID を書いてコミットする。続けて `tdd.py done <ID> --commit <SHA>`。
6. push したら CI を見守る: `gh run watch <id> --exit-status`。

## コードの置き場所

| ディレクトリ | 中身 |
|---|---|
| `include/sssp/` | 公開ヘッダ(共通の `CsrGraph`、`SsspResult`、`Diagnostics`) |
| `src/baseline/` | Dijkstra の各版、Bellman–Ford |
| `src/paper/` | モード A、論文に忠実な実装 |
| `src/experimental/paper_experimental/` | モード B、閾値を調整した実装 |
| `src/experimental/practical/` | モード C、実用的な最適化 |
| `tests/{unit,integration,differential,property,regression}/` | Catch2 のテスト(glob で自動登録) |
| `tests/harness/`, `tests/python/` | pytest |
| `benchmarks/` | Google Benchmark のドライバ |

## レビューで守らせる規則

- テストは先に書き、意図した理由で失敗させる。期待値は仕様・論文・参照実装から取り、実装の出力から取らない。
- 隠れた代替をしない。基本ケースや置き換えは、すべて `Diagnostics.fallbacks` に記録する。
- 計測していない性能の数値を書かない。どの数値も、保存された計測を出典にする。
- 正しさを確かめるビルドで fast-math を使わない。
- force-push をしない。push 済みのブランチの履歴を書き換えない。リモートを変更しない。
- 決定は `docs/adr/` に ADR として記録する。

## スタイル

- C++: `.clang-format`、`.clang-tidy`。CI では警告をエラーとして扱う。
- Python: Ruff(行の長さは 100)。
- `make format` で整形を直し、`make lint` で確かめる。
