# P4-REVIEW: independent correctness review of the Sec 8 components on main

- Reviewed commit: `552b3de` (main).
- Scope: `src/paper/outer/{outer_states,contract_restore,base_case,scales_numeric}.*`,
  `src/paper/sssp.{hpp,cpp}`, `include/sssp/paper.hpp`.
- Out of scope (not on main): Alg 8 l.8-17 outer pipeline (`src/paper/outer/shortest_path*`,
  P3-OUTER-PIPELINE) and the Sec 8.2 scales assembly (P3-OUTER-SCALES).
- Overall verdict: **PASS** for the public entry point and the base-case path. I found no wrong
  distance, no bad tree, no missing or wrong `FallbackRecord`, and no `overflow_error` that the
  rule does not allow. **Needs agreement** on four points about the Sec 8.1 numerics and the
  documentation (items 10, 12, 15 and 16). No counterexample to a correctness claim on the
  public path was found, so no regression test was added. Per the launch instructions, this
  review writes only this file.

## Method and independence

Expectations come from the paper (Sec 2 l.284-304, Sec 8.1 l.1966-2001, Sec 8.2 l.2076-2104,
Sec 8.3 Lemma 8.5 l.2145-2157 and Alg 8 l.1-17), `docs/research/04_outer_scales_pipeline.md`,
`docs/research/thresholds.md`, ADR 0008 and the D-n rows of `docs/PAPER_MAPPING.md`. They do not
come from the implementer's notes.

The oracles were written for this review and share no code with `src/`:

- a 128-bit Bellman-Ford with no cap (exact uint64 distances, even when they exceed `2^64 - 1`);
- a round-based double Bellman-Ford (the minimum over paths of left-to-right rounded sums,
  which is the model `docs/design/numerics.md` uses);
- a Python `fractions.Fraction` Dijkstra (exact rational distances of double inputs);
- a union-find count of the contracted size `(N, m')`;
- a shortest-path-tree checker (source, tight parent edges in 128-bit or fp arithmetic, parent
  chains reach `s`, unreachable iff infinity);
- Python big-integer and `Fraction` checks of every numeric helper.

`verify_sssp` and `baseline::bellman_ford` were used as additional oracles. The harnesses live
in `.cache/review/`, which git ignores, and are not committed. Build: debug preset, `-j4`,
`-fno-fast-math`. The diagnostics-off and ASan/UBSan runs were compiled directly from the
listed sources.

## Verdicts

| # | Item | Verdict |
|---|---|---|
| 1 | Alg 8 l.2 zero-edge contraction (`contract_zero_edges`) | PASS |
| 2 | Alg 8 l.3-4 threshold `N > L_N^64` on the contracted N | PASS |
| 3 | Alg 8 l.5 base case: lightest edge per pair + scan Dijkstra (Lemma 3.2) | PASS |
| 4 | D-1 / ADR 0008: binary-heap variant `paper-faithful-d1` | PASS |
| 5 | Alg 8 l.6 / Lemma 8.5 restoration (`restore_predecessor_tree`) | PASS |
| 6 | Public entry `paper_sssp` vs independent oracles (uint64, double) | PASS |
| 7 | `FallbackRecord`: exactly one per run, correct reason / N / m', also with `SSSP_DIAGNOSTICS=0` | PASS |
| 8 | Boundary values (n = 0, 1, isolated source, zero weights, weights near 2^63 and 2^64 - 1, double +inf / NaN / -0.0 / DBL_MAX / subnormal) | PASS (see notes) |
| 9 | Overflow rule on the public path | PASS |
| 10 | Sec 8.1 deletion test `outer_radial_deleted<double>` (D-21) | Needs agreement |
| 11 | Sec 8.1 deletion test `outer_radial_deleted<uint64_t>` (D-21, Q-58, 128-bit) | PASS |
| 12 | Sec 8.1 `std::overflow_error` for a non-finite `2^7 w` or `2^7 d_e` (double) | Needs agreement |
| 13 | Sec 8.1 special states and arcs (l.1984-1998), owner-tree chain, realization and lower bound | PASS |
| 14 | Sec 8.2 numerics (`scales_numeric`): F, `scale_pow2`, alignment, range, window bound, half test, cutoff, global grid (D-13) | PASS |
| 15 | Header and mapping text vs behaviour | Needs agreement |
| 16 | Meaning of the fallback share in reports | Needs agreement |

## Evidence per item

**1. Contraction (PASS).** For every fuzz case (item 6), the record's `vertices` and `edges`
equal an independent union-find count: zero-edge classes in the component of `s`, and edges
inside that component that join two different classes. Loops are discarded and parallel edges
are kept, as l.2 says; deduplication is l.5. `-0.0` is treated as a zero weight and contracted.
Probes: a zero cycle 0-1-2-0 plus 2-3:7 gives N=2, m'=1; all-zero edges give N=1, m'=0.

**2. Threshold (PASS).** `model::L(k)` agrees with `bitlen(k + 1) = ceil(log2(k + 2))` on
72,446 values: 0..70,000, `2^b + {-3..3}` for b = 1..63, 2,000 random values, `2^64 - 2` and
`2^64 - 1`. `exceeds_scaled_pow(n, c, L(n), e)` agrees with the big-integer test `n > c L^e`
on 394,020 records with e in {1,2,3,4,5,8,12,16,32,63,64} and c in {1,2,3}, so the capped
power never wraps. `takes_outer_construction(n)` was false for every tested n, including
`2^64 - 1`. `paper_sssp` evaluates it on the contracted N (Alg 8 l.3). The ADR 0008 claims were
re-derived with big integers: the predicate is false at `589^64` and true at `589^64 + 1`
(about 1.94e177); true at `2^589 - 2`; false at `2^589 - 1` and at `590^64`; true at
`590^64 + 1`; and it is false for every N below `589^64 + 1`.

**3, 4. Base case and D-1 (PASS).** `outer_base_case` uses `dijkstra_scan_search`, a
minimum scan over the discovered vertices, as Lemma 3.2 describes. `outer_base_case_d1` uses
the indexed binary heap. Both run on `lightest_edge_per_pair` of the contracted graph. Both
variants gave identical distances in every run of item 6, and these equal the oracles.

**5. Lemma 8.5 (PASS).** Every returned tree passed the independent checker: parent edges
join the right endpoints and are tight (exact for uint64, fp equality for double), and parent
chains reach `s`. A traversal of tight arcs from `s` with first-discovery parents is exactly
the proof of Lemma 8.5. An overflow in `d(u) + w` is treated as "not tight", which is correct
because such a sum exceeds every finite label.

**6. Public entry differential (PASS).** Each run is one graph, one weight type and one
variant. Graph families: random multigraphs with loops and parallel edges, paths, stars,
complete graphs, trees with chords, disconnected halves, grids; n from 1 to 28. Weight modes
for uint64: all zero; {0..3} (ties, zero cycles); mostly zero; 32-bit; sets near the sentinel
(`2^62`, `2^63 ± 1`, `2^64 - 2`, `2^64 - 1`, ...); mixed. For double: zero; small integers;
dyadic; uniform non-dyadic; DBL_MAX-scale; subnormal and -0.0; values around 2^53.

| Run set | Seeds | Runs | Failures | Legitimate overflow throws |
|---|---|---|---|---|
| debug lib | 1..3,000 | 12,000 | 0 | 1,436 |
| debug lib | 1,000,000..1,199,999 | 800,000 | 0 | 103,472 |
| debug lib, double dump | 500..20,499 | 80,000 | 0 | 10,250 |
| `SSSP_DIAGNOSTICS=0` build | 7,000,000..7,049,999 | 200,000 | 0 | 25,894 |
| ASan + UBSan build | 9,000,000..9,019,999 | 80,000 | 0 | 10,398 |
| medium size, n in [500, 4000], m = 4n | 24 graphs | 48 | 0 | n/a |

How results were compared:

- uint64: exact equality with the 128-bit Bellman-Ford and with `baseline::bellman_ford`. The
  medium-size runs also compared with `baseline::dijkstra_radix_heap`.
- double: bitwise-equal distances to the fp Bellman-Ford, and `verify_sssp` passed.

The exact rational check covered the 18,304 double graphs of the dump run. The fp result
differs from the true distance at 24,033 vertices, all inside `|d - d*| <= (n-1) eps d*`. The
largest error was 0.525 of that bound.

With k taken as the hop count of a minimum-hop exact shortest path, one vertex reaches
2.0 k eps d* (seed 19402, 13 edges). It is an underestimate: the true path is a single edge,
but the double minimum is attained by a longer path whose sums round down. This is inherent to
the documented double model, and the baselines behave the same way. The tolerance in the
correctness skill ("k = hop count of the reference path") should therefore name the path that
attains the fp minimum, or use n - 1.

**7. FallbackRecord (PASS).** Every successful run, in all builds, had exactly one record. Its
reason string matches ADR 0008 for the variant, and `vertices` and `edges` equal the
independent contracted size. With diagnostics compiled out, the boundary probe output matched
the diagnostics-on build byte for byte (32 records), and the 200,000-run fuzz above found no
record failure.

**8. Boundary values (PASS).**

- n = 0 and `s = n` throw `std::out_of_range`.
- n = 1, a self-loop of weight `2^64 - 1`, an isolated source, and `0-1:0` all return correctly
  with N = 1.
- Double +inf, NaN and negative weights are rejected by `CsrGraph::from_edges` with
  `std::invalid_argument`.
- -0.0 is contracted, and the output distance is +0.
- Subnormal path sums are exact.
- `0-1:2^63 - 1` returns. `0-1:2^63`, `0-1:2^64 - 2` and `0-1:2^64 - 1` throw
  `overflow_error`. The first two are allowed (item 9); the third is required, because the
  true distance reaches the sentinel.
- A heavy parallel edge `2^64 - 1` next to a weight-1 edge returns, because l.5 removes it.
- Double `0-1:DBL_MAX` and `0-1:2^1023` throw. `0-1:2^1022` and `0-1:DBL_MAX/2` return.

**9. Overflow rule (PASS).** Of the 1,172,000 runs, 151,450 threw `std::overflow_error`. Every
one was legitimate: an independent check found a reachable u and an arc u->v with
`d(u) + w >= 2^64 - 1` (double: `fl(d(u) + w) = +inf`). No other exception type occurred.
Whenever some true distance was not representable, the run threw and did not return. The
oracle `baseline::bellman_ford` threw exactly in those cases. So on this path,
`overflow_error` arises only from a tentative path length; the auxiliary quantities here
(contraction, l.5 filter, restoration) never throw.

**10. `outer_radial_deleted<double>` (Needs agreement).** `outer_states.hpp` says "exactness is
established for integer and dyadic weights (Q-55 plan A)". The following integer-valued
inputs contradict that:

- Minimal input: `w = 2^46 + 1`, `d_e = 2^46`, `J'_x = 63`, `J'_y = 64`. Exact: `2^7 w =
  2^53 + 128 > J'_x + J'_y + 2^7 d_e = 2^53 + 127`, so l.1967 deletes the edge. The
  implementation returns `false` (kept), because `127 + 2^53` rounds to `2^53 + 128`. The
  uint64 overload returns `true`. Probe: `.cache/review/deletion_min.cpp`.
- Random search: 200,000 records gave 13 wrong decisions, 3 of them on integer-valued inputs,
  for example `(2^53 + 6, 2^53, 199, 538)`.
- Every wrong decision was in the conservative direction (the paper deletes, the
  implementation keeps). A further 3·10^6 targeted samples found no unsafe deletion.
- Consequence: a kept edge violates the survivor bound `w <= 2(J_x + J_y)` of Lemma 8.1 (here
  `w ≈ 7·10^13` against `2(J_x + J_y) < 2`). Lemma 8.1's swing bound and the realization
  argument rely on that bound.
- Suspected paper step: Sec 8.1 l.1967-1972.
- Proposal: either make the comparison exact (TwoSum-based error-free sum of three doubles), or
  narrow the claim to "exact while every partial sum is below 2^53 in units of the least
  input ulp", and record it as a deviation. The component is not yet on the public path.

**11. `outer_radial_deleted<uint64_t>` (PASS).** It agrees with the big-integer test
`128 w > J'_x + J'_y + 128 d_e` on 200,000 records drawn from small values, full-range values,
`2^64 - 1 - {0..3}`, powers of two ± 1 and `2^57 ± 2`.

**12. Sec 8.1 throw on a non-finite `2^7 w` or `2^7 d_e` for double (Needs agreement).** The
implementation throws `std::overflow_error` whenever either scaled value is not finite, even
when the decision is obvious. For example, `w = 1, d_e = DBL_MAX` gives "not deleted" in exact
arithmetic, yet it throws. These are auxiliary calendar-scaled quantities, not path lengths.
The overflow rule in `include/sssp/paper.hpp` says such quantities must saturate, not throw.
Q-58 / D-20 decide to throw for the calendars of Sec 7.3-7.4, but there is no D-n entry for
Sec 8.1. Upstream Alg 7 would already reject most such inputs. A radial edge whose ends both
map to the root of G' is a loop there and is never seen by Alg 7, so for it the Sec 8.1 throw
would be the first one. Proposal: saturate as D-30 does for Sec 5.1 (if the right side is +inf,
keep the edge), or record the throw as an explicit deviation.

**13. Sec 8.1 structure (PASS).** A Catch2 harness (`.cache/review/outer_struct.cpp`) builds
the Alg 8 l.8-11 inputs with the existing test support. It ran 4,012 cases per weight type:
300 seeds of random connected graphs (n = 6..39, weights up to 2^16; quarter-integers for
double), every separator size from `separator_sizes(n)`, and both name modes. That is
64,076 special-state pairs, 1,296 deleted radial edges and 14,847,151 assertions, all passed.
Checks:

- the deletion set equals the exact test;
- pivots have no special states, and every nonpivot has exactly one `v_in` (inward) and one
  `v_out` (outward), each a singleton block (l.1998);
- entry arcs come only from ordinary inward copies of v, with length 0 (l.1985);
- each `v_in` has exactly one outgoing tree arc, of the physical owner-edge length, to the
  parent's `v_in` or to the owner pivot's ordinary inward copy (l.1986-1987);
- each `v_out` has exactly one incoming arc, from the parent's `v_out` or from the pivot's
  ordinary outward copy (l.1987-1989, "the only ways to enter");
- the exit arc goes to the native outward copy, with length 0;
- macro arcs and all Def 5.4 arcs have ordinary endpoints (l.1990);
- the owner-parent chain reaches `p(v)` in exactly `owner_hops[v]` edges and sums to `R_v`;
- no state distance is below the physical distance of its projection, and the minimum over
  the copies equals the physical distance (Lemma 8.1 realization), by a separate 128-bit
  Dijkstra.

The existing Lemma 8.2 / 8.3 tests (Tarjan SCCs, brute-force theta) passed in ctest; I did not
repeat them.

**14. Sec 8.2 numerics (PASS).** Each helper was compared with Python big integers or
`Fraction`s:

| Helper | Records | Result |
|---|---|---|
| `outer_scale_factor` | 1,287 (p = 0..30, k = 0..8, h including 2, 3, 2^32 ± 1, 2^63, 2^64 - 1, 2^64 - 2) | exponent equals `bitlen(2P(2h+1)(8K^2h)^2)`; value is `min(2^f, 2^64-1)` |
| `scale_pow2` | 100,000 per type (k up to 2^64 - 1 for uint64, up to 5,000 for double) | equal to the saturated exact product |
| `outer_block_aligned`, `outer_aligned_scale_in_range` | 100,000 per type | exact, also for double |
| `outer_window_edge_bound` | 100,000 per type | uint64 exact; double within 64 eps |
| `outer_cutoff_valid`, `first_edge_not_below_half` | 100,000 each | exact |

`outer_scale_factor` throws for h in {0, 1}. `build_outer_global_grid` was tested on 80,000
random instances (uint64 and double, small F, saturating supports, aligned blocks lowering
starts). Checked invariants:

- the points strictly increase;
- each support's bounds equal `[min(a, s_U), sat(F^3 a)]`;
- `lower_index` is the floor of the start;
- `upper_index` is the first point strictly above the end, or the saturated sentinel;
- every class inside a support is exactly `[g, 2g)`;
- `outer_grid_class` agrees with a linear scan on boundary and random probes;
- `outer_aligned_cutoff` satisfies `t <= s < 2t` (Q-36).

All 80,000 passed. The 3,871 throws occurred only on inputs that are not aligned, or have
`s > F a`, or have an infinite `s`; these inputs are invalid by the header. The ASan + UBSan
build of `scales_numeric.cpp` and `outer_states.cpp` reported nothing.
`outer_separate_cutoff` was reviewed in code only.

**15. Header and mapping text (Needs agreement).**

- (a) `paper.hpp` and `base_case.hpp` document `std::overflow_error` "for integer weights"
  only, but double inputs throw it too (item 8: `0-1:DBL_MAX`). The rule paragraph below that
  line is generic.
- (b) `paper.hpp` gives "one edge of weight 2^63" as an example of an input that must not
  throw (the mass S(C)). Yet `paper_sssp` on exactly that input does throw, legitimately,
  through the tentative length `d(1) + w(1->0) = 2^64` (item 8). The example should say that
  it concerns Evaluate's auxiliary mass, not this function's result. Otherwise a reader takes
  the throw for a bug. A related effect: any surviving incident edge of weight `>= 2^64 - 1 -
  d(u)` makes the run throw. The usable range is therefore `max d(u) + w < 2^64 - 1` over
  reachable arcs, not `max d < 2^64 - 1`.
- (c) Already known, the lead will fix it: the stale sentences in `paper.hpp` about Sec 4.2 /
  Sec 5.1 still overflowing, and about Bellman-Ford throwing. The fuzz runs confirm that
  `baseline::bellman_ford` throws only when a true distance is not representable.
- (d) Minor: row 27 of `PAPER_MAPPING.md` names `scales.cpp` for Sec 8.2, but main has
  `scales_numeric.cpp`; `scales.cpp` is still planned under P3-OUTER-SCALES.

**16. Fallback share (Needs agreement).** The record stores the contracted `(N, m')`, which
matches ADR 0008. But `scripts/analyze_results.py` divides `N` by the original `n` and `m'` by
the original `m`. On an input with many zero edges, or with unreachable parts, the reported
"fallback share" is then far below 1, even though no part of the input ran the paper's outer
construction. With all-zero weights the share is `1/n`. Proposal: report the share relative to
`(N, m')`, or report the original size and the contracted size next to each other.

## Tests run

- `cmake --preset debug && cmake --build --preset debug -j4`: 341/341 steps.
- `ctest --preset debug -R "outer|paper_sssp|sssp|base|contract|scales" -j4`: 301/301 passed,
  46.95 s.
- `uv run pytest tests/python/test_exact_reference.py tests/python/test_cross_language.py`:
  207 passed. The cross-language test runs `paper-faithful` and `paper-faithful-d1`.
- Review harnesses (`.cache/review/`, not committed): `fuzz_entry`, `fuzz_large`, `boundary`,
  `numeric` with `numeric_check.py`, `exact_check.py`, `outer_struct`, `deletion_min`. Results
  are in the tables above.

## Not tested

- Alg 8 l.8-17 and Sec 8.2 assembly: out of scope and not on main. The outer branch of
  `paper_sssp` (`outer_construction`) is unreachable for N < 2^64. It throws
  `std::logic_error`, as the documentation says; this was reviewed in code, not executed.
- Release preset binaries. The diagnostics-off behaviour was tested by compiling the path's
  sources with `SSSP_DIAGNOSTICS=0`, not through the CMake release preset.
- Graphs with n above 4,000, and the operation-count (O(N^2 + m)) claim of the scan variant.
  The existing unit test covers the count.
- `outer_separate_cutoff` (code review only), and Lemma 8.2 / 8.3 beyond the existing tests.
