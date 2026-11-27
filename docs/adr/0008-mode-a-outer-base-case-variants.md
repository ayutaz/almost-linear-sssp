# 0008. Mode A outer base case: faithful scan Dijkstra, binary-heap variant as deviation D-1

**English** | [日本語](0008-mode-a-outer-base-case-variants.ja.md) | [简体中文](0008-mode-a-outer-base-case-variants.zh-CN.md)

- Status: Accepted
- Date: 2026-11-27

## Context
Algorithm 8 (Sec 8.3, lines 4–7 of the pseudocode) returns a direct search when
`N <= L_N^64`. The set of `N` satisfying the outer condition `N > L_N^64` first becomes
non-empty at `N = 589^64 + 1` (about 1.94e177) and is not upward-closed: it is false again on
`[2^589 - 1, 590^64]` and holds for every larger `N` only from `590^64 + 1`
(see `docs/research/thresholds.md`). Hence for every input that fits in memory
(`N < 2^64`, where `L_N^64 >= 2^64 > N`) Mode A executes only the base case.

The paper specifies the base case as the scan-based Dijkstra search of Lemma 3.2,
`O(N^2 + m)`. At `N = 10^6` that is about `10^12` operations, so the faithful base case cannot be
benchmarked at the sizes required by the evaluation plan. Lemma 3.2 also provides a
binary-heap search with identical output (exact distances), differing only in running time.

## Decision
- `paper-faithful` (Mode A, default) runs the scan-based base case exactly as written and records
  a `FallbackRecord` with reason `"outer: N <= L_N^64 (paper base case, scan Dijkstra)"`.
- `paper-faithful-d1` runs the binary-heap search of Lemma 3.2 in the base case. It is listed as
  deviation **D-1** in `docs/PAPER_MAPPING.md` and records reason
  `"outer: N <= L_N^64 (paper base case, binary-heap Dijkstra, deviation D-1)"`.
- The threshold test is evaluated on the fly with saturating arithmetic (`L_N^64` computed with
  saturation and compared with `N`), never against a precomputed cut-off, because the
  satisfying set is not upward-closed.
- Benchmarks label which variant was measured. Neither variant is reported as a measurement of
  the paper's non-trivial construction; that construction is exercised only in Mode B
  (`docs/adr/0003-three-implementation-modes.md`).

## Consequences
- Large-scale "paper-faithful" timings use `paper-faithful-d1` and are, by construction, timings
  of a Dijkstra search plus the paper's preprocessing (zero-edge contraction, lightest edge per
  vertex pair). Reports must say so next to every such number.
- The scan variant is still tested for correctness at small `N` so that the faithful code path
  is verified.
