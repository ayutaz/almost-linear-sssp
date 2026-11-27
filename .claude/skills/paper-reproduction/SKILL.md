---
name: paper-reproduction
description: >-
  Map arXiv 2611.34911v1 (theorems, lemmas, algorithms, line numbers, thresholds) to modules and tests, maintain docs/PAPER_MAPPING.md, and record every deviation. Use when implementing or reviewing anything under src/paper/ or src/experimental/paper_experimental/, or when reading the paper.
when_to_use: >-
  "implement Lemma ...", "which section covers ...", "is this faithful to the paper", "add the hierarchy evaluation", any Phase 3 task (P3-*).
---
# Paper reproduction

Source: `.cache/paper/2611.34911v1.txt` (and `.pdf`). Fetch with `scripts/fetch_paper.sh` if
missing. Never commit it, never paste more than a short phrase into the repo; cite
`Sec 5.3, Lemma 5.7, Alg 4 l.9`.

## 1. Before coding a component
1. Locate it: `grep -n "Lemma 5.7\|Threshold Components" .cache/paper/2611.34911v1.txt`.
2. Read the research notes in `docs/research/` for that section.
3. Write/extend its row in `docs/PAPER_MAPPING.md` **before** the code:

| Paper item | Statement (paraphrase) | Module / symbol | Tests | Mode | Status | Deviations |
|---|---|---|---|---|---|---|
| Sec 3.3, Lemma 3.5 | equal-weight MST edges merge simultaneously; no lighter edge between merge children | `src/paper/primitives/merge_events.cpp: build_merge_events` | `tests/unit/paper/test_merge_events.cpp` | A | GREEN | none |

   Status values: TODO / RED / GREEN / VERIFIED (differential + property tests pass) / DEVIATES.
4. Turn each lemma's claim into at least one test (invariant, bound, or exact output on a hand
   example). Pseudocode lines become named steps; cite line numbers in comments
   (`// Alg 4 l.9: open the child with its least due offer`).

## 2. Implementation rules
- Mode A (`src/paper/`): only operations the paper describes. Thresholds exactly as Sec 2:
  `L_N = ceil(log2(N+2))`, outer `N > L_N^64`, local `N > 2*L_N^4`, computed with integer
  saturating arithmetic. Note: for every `N < 2^64`, `L_N^64 >= 2^64 > N`, so the **outer
  construction never triggers on feasible inputs**; Mode A then runs the paper's base case and
  MUST record it (`FallbackRecord{"outer: N <= L_N^64 (paper base case, scan Dijkstra)", n, m}`, or
  `"outer: N <= L_N^64 (paper base case, binary-heap Dijkstra, deviation D-1)"` for the D-1
  variant; see `docs/adr/0008-mode-a-outer-base-case-variants.md`). Report this
  honestly; do not lower the threshold in Mode A.
- Mode B (`src/experimental/paper_experimental/`): same structure, thresholds are parameters.
  Every changed constant is listed in `docs/PAPER_MAPPING.md` (Deviations) with the value used.
- Model of computation (Sec 2): lengths are opaque reals (compare/add/copy only; signed values
  as two non-negative sums; multiply by `2^j` via doublings). If the code uses a real as an
  index, takes floor/log/bits of a weight, or similar, that is a deviation: record it.
- Randomness (Sec 6 pivot sampling): seedable RNG passed in explicitly; tests fix the seed and
  also run several seeds. Expected-time claims are checked by measurement, never assumed.

## 3. Recording deviations
Any of the following is a deviation and needs a row entry (and an ADR if structural):
changed constant, replaced data structure with different bounds, skipped step, extra step not
in the paper, numeric representation change, ambiguity resolved by our own choice.
Format in the Deviations column: `D-<n>: <what> — <why> — <impact on bound/correctness>`.
Ambiguities in the paper go to `docs/research/open-questions.md` with the exact location.

## 4. Claims
- "Paper-faithful" may only be said of Mode A code whose mapping rows are VERIFIED with no
  deviations. Anything else is "paper-inspired" or "experimental".
- Complexity claims about our code come from measured scaling (`/benchmark`) or counted
  operations (`Diagnostics.ops`), not from the paper's theorem.
