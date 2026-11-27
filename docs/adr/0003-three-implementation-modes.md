# 0003. Separate paper-faithful, paper-experimental and practical modes

**English** | [日本語](0003-three-implementation-modes.ja.md) | [简体中文](0003-three-implementation-modes.zh-CN.md)

- Status: Accepted
- Date: 2026-11-27

## Context
The paper's thresholds (outer `N > L_N^64`, local `N > 2 L_N^4`, `L_N = ceil(log2(N+2))`,
Sec 2) mean that the outer construction never runs for any `N < 2^64`. A faithful
implementation therefore spends realistic inputs in its base case. Studying the construction
itself, and practical speed, needs changed constants or different data structures — which must
never be confused with the paper's algorithm.

## Decision
- Mode A, `src/paper/`: paper-faithful. Thresholds exactly as written, compared with
  overflow-safe integer arithmetic. Entering a base case is recorded in
  `Diagnostics.fallbacks`.
- Mode B, `src/experimental/paper_experimental/`: paper structure, parameterised thresholds.
  Every deviation is listed in `docs/PAPER_MAPPING.md`. An input whose correctness cannot be
  established is a failure; automatic substitution of a reference algorithm is never reported
  as success.
- Mode C, `src/experimental/practical/`: practical optimisations, each labelled.
- Baselines live in `src/baseline/` and are not modified to resemble the paper.

## Consequences
- Claims of "the paper's algorithm" refer only to Mode A.
- Reports always show the fallback share next to timings.
- Code sharing between modes goes through explicit shared components, not copy-and-tweak.
