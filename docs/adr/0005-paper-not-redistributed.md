# 0005. The paper is not redistributed
- Status: Accepted
- Date: 2026-11-27

## Context
arXiv 2611.34911v1 is distributed under the arXiv non-exclusive distribution licence, which
does not permit redistribution by third parties.

## Decision
- The PDF and extracted text live only in `.cache/paper/` (git-ignored), obtained with
  `scripts/fetch_paper.sh` from arXiv.
- Repository documents cite the paper by section, lemma, algorithm and line number and
  paraphrase it; quotations are limited to short phrases.
- The target version is v1. Following a later version requires a new ADR.

## Consequences
- Contributors must fetch the paper themselves.
- `docs/PAPER_MAPPING.md` and `docs/research/` are written as paraphrased notes.
