# 0007. fast-math を使わない。正しさの検証では IEEE-754 の意味を保つ

[English](0007-no-fast-math.md) | **日本語** | [简体中文](0007-no-fast-math.zh-CN.md)

- 状態: Accepted
- 日付: 2026-11-27

## 背景 (Context)
`double` の重みの正しさの検査は、文書化した許容誤差のもとで、厳密な有理数の参照と比べます。fast-math による演算の並べ替えや、有限値だけを仮定する最適化は、この比較と、無限大の扱い (`WeightTraits<double>::infinity()`) の両方を損ないます。

## 決定 (Decision)
`cmake/Warnings.cmake` は、すべてのターゲットに `-fno-fast-math` を加えます。`-ffast-math`、`-Ofast`、`-funsafe-math-optimizations` は、正しさの検証に使うビルドでは決して有効にしません。そうしたフラグを使うベンチマークのビルドは、すべてのアルゴリズムに等しく適用し、テスト一式をすべて通し、報告にその旨を明記しなければなりません。

## 結果 (Consequences)
- 浮動小数点の結果は、IEEE-754 の丸めの範囲で、コンパイラが違っても再現できます。
- 一部のベクトル化の機会は使われません。
