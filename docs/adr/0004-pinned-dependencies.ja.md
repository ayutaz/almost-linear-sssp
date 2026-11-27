# 0004. 依存の固定

[English](0004-pinned-dependencies.md) | **日本語** | [简体中文](0004-pinned-dependencies.zh-CN.md)

- 状態: Accepted
- 日付: 2026-11-27

## 背景 (Context)
再現性のためには、テストとベンチマークが、どのマシンでも、CI でも、同じバージョンのライブラリで動く必要があります。

## 決定 (Decision)
- C++:依存は、CMake の `FetchContent` で、リリースタグと SHA-256 の**両方**を指定して取得します (`cmake/Dependencies.cmake`)。Catch2 v3.8.1、Google Benchmark v1.9.4 です。
- Python:uv で管理します。依存は `uv add`(ツール類は `uv add --dev`)で追加し、コミットされた `uv.lock` に固定します。pip は使いません。
- 固定したバージョンを上げるには、理由を述べた新しい ADR(または修正の ADR)が必要です。

## 結果 (Consequences)
- FetchContent のキャッシュがいったんできれば、ビルドはオフラインで再現できます。
- バージョンの更新は、意図的で、レビューされたものになります。
