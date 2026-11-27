---
name: docs-i18n
description: >-
  Write and keep in sync the three language versions of a document (English default, Japanese `.ja.md`, Simplified Chinese `.zh-CN.md`) with the language-switch line, untranslated identifiers and machine-matched strings. Use when adding, translating or updating any user-facing document under docs/ or the README.
when_to_use: >-
  "translate this doc", "add the Japanese version", "update the README", "docs in three languages", "sync the zh-CN page", any edit to a `*.md` that has `.ja.md` / `.zh-CN.md` siblings.
---
# Documentation languages

## 1. Files
- English is the default: `<name>.md`. Siblings: `<name>.ja.md` (Japanese) and
  `<name>.zh-CN.md` (Simplified Chinese), in the same directory.
- Line 3 of each version (after the `# Title`) is the language switch; the current language is
  bold, the others link to the siblings. Model: `docs/README*.md`.
  ```markdown
  **English** | [日本語](README.ja.md) | [简体中文](README.zh-CN.md)
  ```
- Links inside a translated page point to the same-language sibling when one exists.

## 2. Keep the versions in step
- Same sections, tables, rows, numbers and links in all three. A change to one version is
  applied to all three in the same commit; if a translation cannot be done now, say so in the
  report instead of leaving the versions silently different.
- English is the reference: if versions disagree, fix the translations to match it.

## 3. Never translate
- Code identifiers, commands, paths, flags, numbers, task IDs, `D-<n>` / `Q-<n>` labels,
  paper locations (`Sec 5.3, Lemma 5.7`).
- Strings a program matches literally, e.g. the CI / benchmark markers `計測をスキップしました`,
  `未実装`, `データなし`. Keep them verbatim in every version and add a translation in
  parentheses where helpful: `` `データなし` (no data) ``.

## 4. English-only working files
- `docs/TASKS.md`, `docs/PROGRESS.md`, `docs/TDD_LOG.md` are working documents: English only,
  no `.ja.md` / `.zh-CN.md` siblings.
- `docs/TDD_LOG.md` is generated: `uv run python scripts/tdd.py log`. Never edit it by hand.

## 5. Check
```bash
ls <dir>/<name>.md <dir>/<name>.ja.md <dir>/<name>.zh-CN.md   # all three exist
uv run python scripts/check_public_safety.py <changed files>
```
