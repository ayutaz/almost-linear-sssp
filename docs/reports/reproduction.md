# Reproduction report (DRAFT)

**Status: draft skeleton for task P6-REPORT.** The structure and the settled content are written; every
measured result is `not measured`. Nothing in this report is an estimate. The Japanese and Simplified
Chinese versions are made after the content is fixed (see `.claude/skills/docs-i18n/SKILL.md`).

Sources: `docs/PAPER_MAPPING.md` (state at the time of this draft), `docs/reviews/`,
`docs/design/benchmark_protocol.md`, `docs/research/open-questions.md`, `docs/TASKS.md`.
Statuses below are copied from `PAPER_MAPPING.md`; if they change, regenerate this section.

## 1. Subject

- Paper: "Linear Single-Source Shortest Path, Almost" (Gholami, Hajiaghayi), arXiv:2611.34911, **v1 only**.
  Scope: SSSP on undirected graphs with finite non-negative real weights; claimed expected time
  O(m + n log log log n). A later version, if one appears, needs an ADR before it is followed.
- The paper is not redistributable (arXiv non-exclusive licence; ADR 0005). It is not in the repository.
  This report cites it only by section, lemma, algorithm and line number (`l.NNN` of the local text
  extraction). It quotes no passage.
- This project reproduces the algorithm in three modes (ADR 0003): **A** paper-faithful
  (`src/paper/`), **B** paper-experimental (`src/experimental/paper_experimental/`, paper structure with
  tuned thresholds), **C** practical-optimized (`src/experimental/practical/`). Baselines (Dijkstra
  family, Bellman–Ford oracle) are in `src/baseline/`.
- **The outer construction of the paper (Algorithm 8 l.8–17) is not implemented** (row in 2.3,
  status TODO). The public entry point of Mode A reaches only the base case (Algorithm 8 l.4–7), which
  the paper's own size threshold (outer: N > L_N^64) selects at every practical size. Every such call
  records a `FallbackRecord`. A time measured on `paper-faithful` is therefore the time of "the paper's
  preprocessing + Dijkstra", not of the paper's outer construction (see 5).

## 2. What was reproduced

### 2.1 Components

Status vocabulary (from `PAPER_MAPPING.md`): VERIFIED = differential test (against an independent
oracle; for a full SSSP, also shortest-path-tree verification) plus property tests for components;
GREEN = unit tests pass, but the VERIFIED condition is not yet met; DEVIATES = implemented with a
recorded deviation; TODO = not implemented. "Mode" is the mode that contains the component ("A, B" =
shared by Mode A and Mode B, ADR 0003). Mode C has one entry in this table's scope, `practical-local`
(`practical::local_sssp`), a paper-inspired local SSSP with no outer construction. It does not claim
the complexity of Theorem 1, and its performance is `not measured`.

Totals over the 92 rows of 2.3: VERIFIED 53 (including 2 "base branch only"), GREEN 37, the rest
DEVIATES / TODO as listed. Rows of the threshold table (`PAPER_MAPPING.md` section 2) are not repeated
here; they are GREEN except those marked VERIFIED there.

Summary by section (number of rows by status is in the tables):

| Section | Content | Review report | Review verdict |
|---|---|---|---|
| Sec 1–2 | Theorem 1, model, word operations | none | n/a |
| Sec 3 | Primitives (Lemmas 3.1–3.11) | `docs/reviews/P4-REVIEW-owner-model.md` | FAIL (overflow contract F-1, fixed in P3-CALENDAR-SCALE-OVERFLOW; no wrong value found) |
| Sec 4 | Evaluation of ordered hierarchies | `docs/reviews/P4-REVIEW-hierarchy.md` | FAIL (F-1 double, F-4 uint64 origin rule; F-2, F-3 invalid input only) |
| Sec 5 | Local hierarchy | `docs/reviews/P4-REVIEW-local.md` | FAIL (counterexample 1: `overflow_error`, decided as Q-83; counterexample 2: double Evaluate) |
| Sec 6–7 | Pivot routes, owner distances | `docs/reviews/P4-REVIEW-owner-model.md` | FAIL (same report as Sec 3) |
| Sec 8 | Outer states, scales, entry point | `docs/reviews/P4-REVIEW-outer.md` | PASS for the public entry and the base-case path |

The verdicts are the reviewers' own and are reproduced as recorded. Which counterexamples are fixed
at the time of reading is in `docs/TASKS.md` (P3-EVAL-ARITH-FIX is DONE; P3-EVAL-RANGE-COPIES,
P3-EVAL-DOUBLE-ROUNDING, P3-OUTER-DELETE-DOUBLE, P3-CALENDAR-SCALE-OVERFLOW,
P3-ESTIMATE-RADIUS-DOUBLE were IN_PROGRESS when this draft was written). Re-check before release.

### 2.2 Components still GREEN

The GREEN rows with a known cause are those with a task in the last column of 2.3:

| Component | Reason | Task |
|---|---|---|
| Sec 4.3 Scan, Algorithm 2 | suspected source of the double error | P3-EVAL-DOUBLE-ROUNDING |
| Sec 4.3 Evaluate, Algorithm 3 | wrong value for full-precision double; range of state copies | P3-EVAL-DOUBLE-ROUNDING, P3-EVAL-RANGE-COPIES |
| Sec 8.1 outer state graph, Lemma 8.1 | double deletion test keeps edges the paper deletes | P3-OUTER-DELETE-DOUBLE |
| Sec 8.2 scales, Lemma 8.4 | unaligned inputs not checked end to end (Q-79, Open) | P3-OUTER-SCALES |

The other GREEN rows have no task ID in `PAPER_MAPPING.md`; their reason is in the row there.
Review findings not yet reflected in a status (hierarchy F-1 for double gaps, D-15 / D-27) are listed
in the review report; this draft does not assign them a task.

### 2.3 Component table

Generated from the rows of `PAPER_MAPPING.md` sections 1 (Theorems, definitions, lemmas, algorithms).
Evidence per row is in the "Tests" and "Review" lines of each subsection and in the full row of
`docs/PAPER_MAPPING.md`.

#### Sec 1–2 Main theorem and preliminaries

Tests: unit, property. Review: none (Theorem 1 is covered by the base-branch differential test only).

| Paper item | Mode | Status | Why GREEN, or task |
|---|---|---|---|
| Sec 1, Theorem 1 (l.91, proof l.2160) | A | VERIFIED (base branch only) |  |
| Sec 2, l.284: L_k and τ(k) | A | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 2, l.286–290: size assumptions | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 2, l.292–300: real computation model | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 2, l.301–306: integer word operations | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 2, l.307–310: randomness | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |

#### Sec 3 Primitives

Tests: unit, differential, property (reference structures written in the tests). Review: `reviews/P4-REVIEW-owner-model.md` (primitives, model).

| Paper item | Mode | Status | Why GREEN, or task |
|---|---|---|---|
| Sec 3.1, Lemma 3.1 (l.336) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 3.1, Lemma 3.2 (l.351) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 3.2, Lemma 3.3 (l.409) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 3.2, Algorithm 1 (l.376) + Lemma 3.4 (l.437) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 3.2, Lemma 3.5 (l.457) | A, B | VERIFIED |  |
| (Supplement, 01-A12; not in the paper) spanning forest of a disconnected input | A (supplement) | VERIFIED |  |
| Sec 3.3 / use of Lemma 3.5: ordered merge events | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Alg 4 l.4, Alg 5 l.4, Alg 8 l.5: lightest edge per vertex pair | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 3.3, Lemma 3.6 (l.481) | A, B | VERIFIED |  |
| Sec 3.4, Lemma 3.7 (l.521) | A, B | VERIFIED |  |
| Sec 3.4, Lemma 3.8 (l.550) | A, B | VERIFIED |  |
| Sec 3.4, Lemma 3.9 (l.570) | A, B | VERIFIED |  |
| Sec 3.5, Lemma 3.10 (l.592) | A, B | VERIFIED |  |
| Sec 3.5, Lemma 3.11 (l.625) | A, B | VERIFIED |  |
| (Test support, not in the paper) reference split–findmin | A (test support) | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |

#### Sec 4 Evaluation of ordered hierarchies

Tests: unit (hand traces), differential against Dijkstra / Floyd–Warshall, property (I1–I8), regression. Review: `reviews/P4-REVIEW-hierarchy.md` (FAIL, see 2.2).

| Paper item | Mode | Status | Why GREEN, or task |
|---|---|---|---|
| Sec 4, Definition 4.1 (l.646) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 4, Definition 4.2 (l.659) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| (Supplement) Reference hierarchy: SCCs per threshold (02 §2.11) | A (test support) | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 4.1, Lemma 4.3 (l.674, l.690) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 4.2, Lemma 4.4 (l.731) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 4.2, Lemma 4.5 (l.758) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 4.3 (l.781–832): slots, deadlines, runs, supplemental offers | A, B | VERIFIED |  |
| Sec 4.3, Algorithm 2 (l.834) Scan | A, B | GREEN | Suspected cause of the double-weight error of P4-REVIEW-local counterexample 2 (slot / origin arithmetic); root cause to be confirmed. Task: P3-EVAL-DOUBLE-ROUNDING. |
| Sec 4.3, Algorithm 3 (l.866) Evaluate | A, B | GREEN | Known wrong value (below the true distance) for full-precision double (P4-REVIEW-local counterexample 2); a state copy can raise `overflow_error` although every d[u]+w fits (counterexample 1; Q-83 decided, D-30 adopted). Task: P3-EVAL-DOUBLE-ROUNDING, P3-EVAL-RANGE-COPIES. |
| Sec 4.3, Lemma 4.6 (l.885) | A, B | VERIFIED |  |
| Sec 4.3, Lemma 4.7 (l.922) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |

#### Sec 5 Local hierarchy

Tests: unit, differential, property, regression. Review: `reviews/P4-REVIEW-local.md` (FAIL, see 2.2).

| Paper item | Mode | Status | Why GREEN, or task |
|---|---|---|---|
| Sec 5.1 (l.960), Definition 5.1 (l.975) | A, B (outer uses b = L^64, Sec 7.1) | VERIFIED |  |
| Sec 5.1, Lemma 5.2 (l.993) | A, B | VERIFIED |  |
| Sec 5.1, Lemma 5.3 (l.1011) | A, B | VERIFIED |  |
| Sec 5.2 (l.1020–1055) edge classification, shells, windows | A, B (outer uses coordinate J and K = 16, Sec 8.1) | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 5.2, Definition 5.4 (l.1056) | A, B | VERIFIED |  |
| Sec 5.2, Lemma 5.5 (l.1061) | A, B | VERIFIED |  |
| Sec 5.2, Lemma 5.6 (l.1074) | A, B | VERIFIED |  |
| Sec 5.2, Lemma 5.7 (l.1084) | A, B | VERIFIED |  |
| Sec 5.3 (l.1109), Lemma 5.8 (l.1131) | A, B | VERIFIED |  |
| Sec 5.3, Lemma 5.9 (l.1152) | A, B | VERIFIED |  |
| Sec 5.4 (l.1182–1243) assembly | A, B | VERIFIED |  |
| Lemma 5.10 proof, post-processing (l.1279–1293) | A, B | VERIFIED |  |
| Sec 5.4 l.1206–1243; Lemma 5.10 proof l.1265–1300 (shared core) | ? | VERIFIED |  |
| Sec 5.4, Algorithm 4 l.2–5 (l.1244, Sec 2 l.286–290) base case | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 5.4, Algorithm 4 l.2 and l.6–12 + Lemma 5.10 (l.1265) | A, B | VERIFIED |  |
| Sec 5.4 l.1206–1221, Mode B (D-8) | B | VERIFIED |  |
| Sec 5.4, Mode B (D-11): exact-height assembly | B | VERIFIED |  |
| Alg 4 l.2–12, Mode B entry point | B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 5.5, Definition 5.11 (l.1303) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 5.5, Definition 5.12 (l.1308) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 5.5, Lemma 5.13 (l.1320) | A, B | VERIFIED |  |
| Sec 5.5, Lemma 5.14 (l.1342) | A, B | VERIFIED |  |
| Sec 5.5, Lemma 5.15 (l.1355) | A, B | VERIFIED |  |

#### Sec 6 Pivot route certificates

Tests: unit, property (Las Vegas), differential. Review: `reviews/P4-REVIEW-owner-model.md`.

| Paper item | Mode | Status | Why GREEN, or task |
|---|---|---|---|
| Sec 6, Definition 6.1 (l.1368) | A, B | VERIFIED |  |
| Sec 6.1 (l.1382–1392) preprocessing and median narrowing | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 6.1, Algorithm 5 (l.1394) NameRoutes | A, B | VERIFIED |  |
| Sec 6.1 Algorithm 5, Sec 6.2 l.1491, Lemma 6.8: Mode B entry point | B | VERIFIED |  |
| Sec 6.1, Lemma 6.2 (l.1424) | A | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 6.1, Lemma 6.3 (l.1439) | A | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 6.2, Algorithm 6 (l.1454) Estimate | A, B | VERIFIED |  |
| Sec 6.2, Lemma 6.4 (l.1509) | A, B | VERIFIED |  |
| Sec 6.2, Lemma 6.5 (l.1530) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 6.2, Lemma 6.6 (l.1542) | A, B | VERIFIED |  |
| Sec 6.2, Lemma 6.7 (l.1560) | A, B | VERIFIED |  |
| Sec 6.3, Lemma 6.8 (l.1605) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |

#### Sec 7 Auxiliary distances to owner pivots

Tests: unit, property, differential. Review: `reviews/P4-REVIEW-owner-model.md` (FAIL on the overflow contract, see 2.2).

| Paper item | Mode | Status | Why GREEN, or task |
|---|---|---|---|
| Sec 7.1 (l.1642–1652) edge classification | A, B | VERIFIED |  |
| Sec 7.1, Lemma 7.1 (l.1653) | A, B | VERIFIED |  |
| Sec 7.2 (l.1662–1667), Lemma 7.2 (l.1668) | A, B | VERIFIED |  |
| Sec 7.3 (l.1692–1710) pieces, grids, calendars, atoms | A, B | VERIFIED |  |
| Sec 7.3, Lemma 7.3 (l.1711) | A, B | VERIFIED |  |
| Sec 7.3, Lemma 7.4 (l.1736) | A, B | VERIFIED |  |
| Sec 7.3, Lemma 7.5 (l.1751) | A, B | VERIFIED |  |
| Sec 7.3, Lemma 7.6 (l.1762) | A, B | VERIFIED |  |
| Sec 7.4 (l.1780–1841) families, runs, atom visits | A, B | VERIFIED |  |
| Sec 7.4, Algorithm 7 (l.1843) OwnerDistances | A, B | VERIFIED |  |
| Sec 7.4, Lemma 7.7 (l.1871) | A, B | VERIFIED |  |
| Sec 7.4, Lemma 7.8 (l.1919) | A, B | VERIFIED |  |
| Sec 7.4, Lemma 7.9 (l.1944) | A, B | VERIFIED |  |

#### Sec 8 Outer hierarchy and single-source algorithm

Tests: unit, differential, integration, regression. Review: `reviews/P4-REVIEW-outer.md` (PASS for the public entry and the base-case path).

| Paper item | Mode | Status | Why GREEN, or task |
|---|---|---|---|
| Sec 8.1 (l.1966–2001) outer state graph | A, B | GREEN | The double deletion test of l.1967 keeps edges the paper deletes (P4-REVIEW-outer item 10). Task: P3-OUTER-DELETE-DOUBLE. |
| Sec 8.1, Lemma 8.1 (l.2002) | A, B | GREEN | Same cause: the kept edges can break the weight bound of Lemma 8.1 by a large factor (safe for distances). Task: P3-OUTER-DELETE-DOUBLE. |
| Sec 8.1, Lemma 8.2 (l.2027) | A, B | VERIFIED |  |
| Sec 8.1, Lemma 8.3 (l.2053) | A, B | VERIFIED |  |
| Sec 8.2 (l.2076–2104) numerical scales and atom substitution | A | GREEN | On unaligned inputs the distances are not checked end to end (Q-79, Open). Task: P3-OUTER-SCALES. |
| Sec 8.2, Lemma 8.4 (l.2105) | A | GREEN | Same as the Sec 8.2 row (Q-79, Open). Task: P3-OUTER-SCALES. |
| Sec 8.3, Algorithm 8 l.2 contraction of zero-length edges | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 8.3, Lemma 8.5 (l.2145) | A, B | GREEN | No test at the granularity required for VERIFIED yet (see the row in `PAPER_MAPPING.md`); no task recorded. |
| Sec 8.3, Algorithm 8 l.4–7 base case (Sec 2 l.286–288, Lemma 3.2) | A | VERIFIED |  |
| Sec 8.3, Algorithm 8 l.4–7 base case, binary-heap variant | A (D-1) | DEVIATES | Deviation D-1 (ADR 0008): binary heap instead of the scan. |
| Sec 8.3, Algorithm 8 l.8–17 outer construction | A (not reached), B | TODO | Not implemented. In Mode A `paper_sssp` throws `std::logic_error` on the outer branch; it is not hidden behind Dijkstra. |
| Sec 8.3, Algorithm 8 l.2–7 (l.2167) ShortestPaths entry point | A | VERIFIED (base branch only) |  |

## 3. Deviations from the paper

Full list with descriptions, status and related questions: `docs/PAPER_MAPPING.md` section 3. Decisions
that constrain future work are ADRs in `docs/adr/`; the only ADR on a structural deviation so far is
[ADR 0008](adr/0008-mode-a-outer-base-case-variants.md) (D-1). Related: [ADR 0003](adr/0003-three-implementation-modes.md)
(modes), [ADR 0007](adr/0007-no-fast-math.md).

| ID | Summary | Status | Applies to |
|---|---|---|---|
| D-1 | Base case of Alg 8 l.5 by binary-heap Dijkstra (`paper-faithful-d1`) instead of the scan version | Adopted (ADR 0008) | A |
| D-2 | `std::bit_width` instead of the paper's logarithm table | Candidate (ADR required) | A, B |
| D-3 – D-6 | Simpler structures (segment-tree split–findmin, binary-lifting LCA/LA, binary search for tour start, binary/4-ary heap for the Fibonacci heap) | Candidates; D-4 and D-6 are not used, D-3 and D-5 are provisional | A, B |
| D-7 – D-11, D-17 | Mode B parameters: thresholds, separator size, Sec 5.4 constants, H, NameRoutes cap and factors, exact-height assembly | Mode B (D-8 implemented in P5-EXP-LOCAL) | B |
| D-9, D-10 | H lowered to the least power of two at least 4tT²; NameRoutes retry cap | Mode B (D-10: agreement required) | B |
| D-12 | −0.0 normalised to +0.0 | Candidate | A, B |
| D-13 | Polynomial multiples of reals saturate instead of overflowing | Adopted (Q-38) | A, B |
| D-14, D-16 | Test-only switches (supplemental offers off; explicit block size) | Test-only | tests |
| D-15, D-19, D-22 | Signed pairs of sums materialised as one real | D-15, D-19 candidates; D-22 adopted | A, B |
| D-18 | Block index by integer division | Resolved in P3-WORD-OPS | A, B |
| D-20, D-21 | uint64 calendar range condition; Alg 7 runs in 2^7-scaled units | Adopted | A, B |
| D-23 – D-27 | Saturation, 128-bit exact forest lengths, double rejection rules for Sec 4 | Candidates | A, B |
| D-28 – D-30 | Same ideas for Sec 5.1 (128-bit root distances, saturated R_v, double deletion test) | D-30 adopted; D-28, D-29 candidates | A, B |
| D-31 | Class of a block MST edge by binary search over the piece table | Candidate | A, B |

"Candidate" means recorded in `PAPER_MAPPING.md` and not yet accepted by an ADR. Where a deviation is in
effect it is listed in the Deviations column of the component table there.

## 4. How correctness was checked

Rules: `CLAUDE.md` (correctness rules, TDD), `.claude/skills/correctness-review/SKILL.md`.

- **Differential tests** (`tests/differential/`). Each SSSP implementation is compared with an
  independent oracle (Bellman–Ford, Floyd–Warshall for small hierarchies, exact reference).
  Integer weights: exact equality. `double`: a documented tolerance
  (`Tolerance::for_path_hops`, `docs/design/numerics.md`); no tolerance was loosened to pass.
- **`verify_sssp`** (`include/sssp/verify.hpp`, `src/verify/verify.cpp`). Shortest-path-tree
  verification: `dist[s] = 0`; every arc satisfies `dist[v] <= dist[u] + w`; every tree arc is tight;
  parents reach `s`; unreachable vertices are infinite. It is also applied to every benchmark trial
  (`docs/design/benchmark_protocol.md` section 5).
- **Property tests** (`tests/property/`). Used for components: invariants stated by the lemmas
  (e.g. I1–I8 of Evaluate, Lemma 4.7 bounds, Las Vegas behaviour of Estimate, calendar properties).
  Work bounds are counted with diagnostics on; running times are not measured here.
- **Regression tests** (`tests/regression/`). Every minimised counterexample is added before its fix,
  with a record in `tests/regression/cases/`.
- **Sanitizers.** `make sanitizer` (preset `sanitizer`, `ctest --preset sanitizer`) and the `sanitizer`
  GitHub Actions workflow. Results of runs are not summarised in this draft: not measured / to be
  filled from the final CI run id.
- **Independent review.** Four reviews by a correctness reviewer who did not write the code
  (`docs/reviews/`): hierarchy (Sec 4), local (Sec 5), owner / pivots / primitives / model, and the
  outer components with the public entry. The review of the outer construction and of the Mode B
  differential check will be added once those exist. Review fuzz counts and verdicts are in the reports;
  they are not repeated here.
- **Overflow contract** (`include/sssp/paper.hpp`; uint64 weights). A result is either exact or the call
  throws; a wrong value is never returned.
  - **Q-81** (decided): in the Estimate sample search, a path length d + w reaching 2^64 throws
    `std::overflow_error`; it is the length of a real path, so the contract allows it.
  - **Q-83** (decided, D-30 adopted): the d[u] + w of the contract is a tentative path length to the
    vertex itself. The constrained distance of a state copy is an auxiliary quantity; a copy that received
    only offers beyond the range stays without a value, and `overflow_error` is thrown only when no copy
    of a vertex has a value.
- **No fast-math.** `cmake/Warnings.cmake` forces `-fno-fast-math` (ADR 0007).

## 5. Measurements

### 5.1 Method

The full procedure is `docs/design/benchmark_protocol.md` (with the schema in
`docs/design/benchmark_schema.md` and the runner in `docs/design/benchmark_runner_cli.md`). Summary:

- `release` preset only (`SSSP_DIAGNOSTICS=OFF`); counters come from a separate `release-diag` run and
  are never mixed with times.
- Same `CsrGraph`, same sources (chosen from a seed), same warm-up and repetitions for every algorithm;
  job order shuffled with a fixed seed; one process per job (peak RSS).
- Times are separated into preprocess / sssp / tree / end-to-end; graph generation and I/O are excluded.
- Warm-up: 1 (standard), 2 (full). Repetitions: 5 for development, 20 or more for published numbers.
- Every trial is verified (`verify_sssp`); a failed verification is recorded as a failure, not dropped.
- Statistics: median, mean, p95, standard deviation and n.
- Machine conditions (AC power, thermal check, load recorded before and after) and the environment
  metadata saved next to the results: protocol section 7 and the `benchmark` skill.
- GitHub-hosted runners are not used for numbers in this report.
- Compared algorithms: `dijkstra-binary`, `dijkstra-4ary`, `dijkstra-fibonacci`, `dijkstra-radix`
  (uint64 only), `paper-faithful`, `paper-faithful-d1`, `paper-experimental`, `practical-local`
  (uint64 only). `paper` (outer construction) and `paper-optimized` are reserved names and are not
  implemented.

### 5.2 Fallback share

Every `paper-faithful` call enters the base case and records one `FallbackRecord{reason, vertices,
edges}`; `practical-local` records one when the contracted vertex count is 16 or less. **A time is
always reported together with the fallback share, in the column next to it.**

**The denominator of the share is not decided (Q-82, Open).** The record stores the contracted
(N, m'), but `scripts/analyze_results.py` divides by the original n and m, so an input with many
zero-length edges or unreachable parts shows a share far below 1 although no part ran the outer
construction. The fix is task P6-FALLBACK-SHARE-DENOM. Until then the columns below are shown with
both the contracted and the original sizes, and shares are not compared across inputs.

### 5.3 Results (frame only)

Environment (machine, OS, compiler, git SHA, run id): `not measured`.

Times in nanoseconds, median [p95] over n trials. Cells marked `not measured` wait for task
P6-BENCH-STANDARD (and P6-SCALING for the fits).

| Graph type / weights | n | m | Algorithm | n trials | median | mean | p95 | stddev | Fallback share, contracted (N/m') | Fallback share, original (N/n, m'/m) |
|---|---|---|---|---|---|---|---|---|---|---|
| not measured | not measured | not measured | `dijkstra-binary` | not measured | not measured | not measured | not measured | not measured | not applicable | not applicable |
| not measured | not measured | not measured | `paper-faithful` | not measured | not measured | not measured | not measured | not measured | not measured | not measured |
| not measured | not measured | not measured | `paper-faithful-d1` | not measured | not measured | not measured | not measured | not measured | not measured | not measured |
| not measured | not measured | not measured | `paper-experimental` | not measured | not measured | not measured | not measured | not measured | not measured | not measured |
| not measured | not measured | not measured | `practical-local` | not measured | not measured | not measured | not measured | not measured | not measured | not measured |

Peak RSS: `not measured`.

Scaling fit (task P6-SCALING, tool `scripts/scaling_fit.py`, `make scaling`): at least 4 sizes over at
least 2 orders of magnitude; fitted exponent and the comparison of operation counts with
m + n log log log n: `not measured`. A fit over `paper-faithful` would describe the base case (scan or
heap Dijkstra), not the paper's bound; this must be stated beside any such result.

Speedup relative to `dijkstra-binary`: `not measured`.

No statement about the speed of the paper's algorithm is made in this draft.

## 6. Known limitations and open questions

Limitations that are facts of the repository today:

- The outer construction (Algorithm 8 l.8–17) is not implemented, so the claimed O(m + n log log log n)
  is **not reproduced**; it is not measured either.
- At every practical size N the paper's own threshold sends the input to the base case; this is legal
  and is recorded as a fallback (Mode A), but it means the paper's construction is exercised only
  through its components (Sec 4–8.2) in tests, not end to end.
- Reviews FAIL for Sec 4, Sec 5 and Sec 6–7 (see 2.2 and 2.3); some fixes are IN_PROGRESS.
- Tree reconstruction for general `double` weights is not supported in `practical-local`
  (uint64 only).
- The unaligned block of Sec 8.2 cannot be run end to end (Q-79).

Open questions (`docs/research/open-questions.md`; status as read at the time of this draft):

| ID | Subject | State |
|---|---|---|
| Q-79 | Sec 8.2 unaligned blocks cannot be run end to end (uint64 and double both blocked) | Open |
| Q-82 | Denominator of the fallback share | Open (P6-FALLBACK-SHARE-DENOM) |
| Q-05 | Integer log table or `bit_width` (D-2) | Needs agreement |
| Q-43 | Multi-word integers H and j in Estimate | Needs agreement |
| Q-45 | Capping NameRoutes in Mode B (D-10) | Needs agreement |
| Q-50 | Height of a contracted node after normalisation | Needs agreement |

Also not decided in the repository when this was written: the destination of the final report (see
`docs/TASKS.md` P6-REPORT, which names `docs/reports/reproduction.md`; this draft is `docs/REPORT.md`).

## 7. How to reproduce

Prerequisites: C++20 compiler, CMake, Ninja, `uv`. The paper is fetched locally (not committed).

```sh
make bootstrap                 # toolchain / dependencies check
make fetch-paper               # paper text into .cache/paper/ (git-ignored)

# build and test (preset ci: -Werror)
make PRESET=ci build
ctest --preset ci
make test-differential PRESET=ci
make test-property PRESET=ci
make test-regression PRESET=ci
make sanitizer                 # ASan/UBSan preset

# full local gate
make verify

# benchmarks (release preset, idle machine; see the protocol)
cmake --preset release && cmake --build --preset release -j
make benchmark-smoke           # quick sanity
make benchmark                 # standard suite
make benchmark-full            # large graphs
make report                    # aggregation and figures
make scaling                   # scaling fit
```

A benchmark run is valid only under the conditions of `docs/design/benchmark_protocol.md` section 7
(idle machine, AC power, no thermal limit). Run ids and the commit SHA of the measured build go with
the results.
