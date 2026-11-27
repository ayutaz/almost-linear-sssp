# 0006. 計測用ビルドでは診断をコンパイルで除く

[English](0006-diagnostics-compiled-out-for-timing.md) | **日本語** | [简体中文](0006-diagnostics-compiled-out-for-timing.zh-CN.md)

- 状態: Accepted
- 日付: 2026-11-27

## 背景 (Context)
実行経路のカウンタと演算回数(`include/sssp/result.hpp` の `Diagnostics`)は、論文のどの構成が走ったか、フォールバックがどれだけ起きるかを確かめるのに必要ですが、時間の計測を乱します。

## 決定 (Decision)
カウンタは `SSSP_COUNT(...)` で更新し、コンパイル定義 `SSSP_DIAGNOSTICS` で制御します。`release` プリセットは `SSSP_DIAGNOSTICS=OFF` に設定し、時間の計測に使う唯一のプリセットです。回数は、同じ入力に対する別の `release-diag` の実行から得ます。`FallbackRecord` は常に記録します(カウンタではなく、結果の一部であるため)。

## 結果 (Consequences)
- 計測した実行には、カウンタのオーバーヘッドがありません。
- 報告では、同じコミットと入力から、時間 (release) と回数 (release-diag) を組み合わせます。
