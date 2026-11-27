# 0007. No fast-math; IEEE-754 semantics for correctness

**English** | [日本語](0007-no-fast-math.ja.md) | [简体中文](0007-no-fast-math.zh-CN.md)

- Status: Accepted
- Date: 2026-11-27

## Context
Correctness checks for `double` weights compare against an exact rational reference with a
documented tolerance. Fast-math reassociation and finite-math assumptions would invalidate
both the comparisons and infinity handling (`WeightTraits<double>::infinity()`).

## Decision
`cmake/Warnings.cmake` adds `-fno-fast-math` to every target. `-ffast-math`, `-Ofast` and
`-funsafe-math-optimizations` are never enabled in builds used for correctness. A benchmark
build with such flags would have to apply them to all algorithms equally, pass the full test
suite, and be labelled in the report.

## Consequences
- Floating-point results are reproducible across compilers up to IEEE-754 rounding.
- Some vectorisation opportunities are not taken.
