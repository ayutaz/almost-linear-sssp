# 0005. 論文は再配布しない

[English](0005-paper-not-redistributed.md) | **日本語** | [简体中文](0005-paper-not-redistributed.zh-CN.md)

- 状態: Accepted
- 日付: 2026-11-27

## 背景 (Context)
arXiv 2611.34911v1 は arXiv の非独占的配布ライセンスで配布されており、第三者による再配布は許されていません。

## 決定 (Decision)
- PDF と抽出したテキストは `.cache/paper/`(git の管理外)にだけ置き、`scripts/fetch_paper.sh` で arXiv から取得します。
- リポジトリの文書は、論文を節、補題、アルゴリズム、行番号で引用し、言い換えて書きます。引用は短い語句に限ります。
- 対象のバージョンは v1 です。後のバージョンに従うには、新しい ADR が必要です。

## 結果 (Consequences)
- 貢献者は、論文を自分で取得する必要があります。
- `docs/PAPER_MAPPING.md` と `docs/research/` は、言い換えたノートとして書きます。
