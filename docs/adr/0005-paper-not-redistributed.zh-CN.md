# 0005. 不再分发论文

[English](0005-paper-not-redistributed.md) | [日本語](0005-paper-not-redistributed.ja.md) | **简体中文**

- 状态: Accepted
- 日期: 2026-11-27

## 背景 (Context)
arXiv 2611.34911v1 按 arXiv 的非独占分发许可证发布,该许可证不允许第三方再分发。

## 决定 (Decision)
- PDF 和提取出的文本只放在 `.cache/paper/`(被 git 忽略),用 `scripts/fetch_paper.sh` 从 arXiv 获取。
- 仓库中的文档按章节、引理、算法和行号引用论文,并用自己的话转述;引文仅限于简短的词句。
- 目标版本是 v1。要跟进后来的版本,需要一份新的 ADR。

## 后果 (Consequences)
- 贡献者必须自己获取论文。
- `docs/PAPER_MAPPING.md` 和 `docs/research/` 以转述的笔记形式写成。
