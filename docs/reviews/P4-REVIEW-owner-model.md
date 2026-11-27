# P4-REVIEW: independent correctness review of `owner/`, `pivots/`, `primitives/`, `model/`

Scope: arXiv 2611.34911v1, Sec 2–3 and Sec 6–7, as implemented in `src/paper/model/`,
`src/paper/primitives/`, `src/paper/pivots/` and `src/paper/owner/` (base commit `552b3de`).
Out of scope (another reviewer): `src/paper/outer/`, `src/paper/hierarchy/`, `src/paper/local/`.
They were treated as black boxes where the reviewed code calls them (Alg 4 per atom,
`PortionedEvaluator`).

Method: expectations were taken from the paper (cited by section, lemma and line of the text
dump), `docs/PAPER_MAPPING.md` (rows of Sec 2, 3, 6, 7, the constants table and the D table)
and `docs/research/`. Every oracle below was written for this review and shares no code with
`src/` (own Dijkstra over exact 128-bit values, own Floyd–Warshall, own Kruskal, naive LCA and
level-ancestor walks, naive interval minima, Python big integers). The throw-away harnesses live
in `.cache/review/` (git-ignored, not committed).

**Overall verdict: FAIL (one contract violation, no wrong value found).** No wrong distance,
certificate, threshold or primitive answer was found in any run. One counterexample violates the
overflow rule of `include/sssp/paper.hpp` and the stated range of D-20 (F-1). Four points need
agreement (F-2 to F-5).

## 1. Verdicts by section and lemma

| Item (paper) | Code | Verdict | Evidence |
|---|---|---|---|
| Sec 2 l.284: L_k = ⌈log2(k+2)⌉ | `model/params.hpp: L` | PASS | E1: equal to a big-integer definition (least c with 2^c ≥ k+2) on N = 0..3000, 2^k−3..2^k+3 for k = 1..63, 2^64−4..2^64−1 and 3000 random words |
| Sec 2 l.286: local N > 2L_N^4, outer N > L_N^64 | `params.hpp: local_nontrivial, outer_nontrivial` | PASS | E1: equal to the big-integer truth value on every tested N; first local N is 209,953 and the set is upward closed up to 2^22; outer is false for every 64-bit N |
| Sec 2 l.304: capped products in size tests | `capped_pow, capped_mul, exceeds_scaled_pow` | PASS | E1: `capped_pow(b, e, cap) = min(b^e, cap)` for b ∈ {0, 1, 2, 3, 5, 18, 65, 2^20, 2^64−1}, e = 0..70, 7 caps including 0 and 2^64−1; Mode B truth values N > c·L^e for e = 0..70 and c up to 2^64−1. Products are capped before they can wrap |
| Sec 2 l.300: no real log/pow/floor | whole scope | PASS | `grep` finds no `std::log*`, `pow`, `floor`, `ceil`, `ldexp`, `frexp`, `bit_width`, `countl_zero` in the four directories. The only real division is the D-21 halving (`owner_distances.cpp` `calendar_unscale`) |
| Sec 2 l.301–306: integer logarithm table, block lookups | `model/ilog.hpp: IlogTable`, `model/word_ops.hpp: BlockLayout` | PASS | table equals the word versions on 1..2^16; `BlockLayout` equals `/` and `%` for every k ≤ 300 and block size 1..k+2 (D-18 resolved as documented) |
| Sec 2 l.307–310: randomness | `model/random.hpp` | PASS | E2: `uniform_below(h)` requests exactly ⌈log2 h⌉ bits per attempt (h = 1 draws nothing; h > 2^63 uses 64), all outputs < h, fewer than 2 attempts per draw on average; χ² for h = 3, 5, 7 (7·10^5 draws), all 24 permutations of 4 (4.8·10^5), all 20 ordered 2-prefixes of `partial_shuffle(5, 2)` (4·10^5) are below the 0.999 quantile; same seed gives the same words; SplitMix64 reference value for seed 0 |
| Sec 2 l.292–296: lengths, infinity flag, signed pairs | `model/length.hpp` | PASS | E3: −0.0 becomes +0.0 (D-12); subnormal lengths are accepted and added exactly; NaN, ±inf, negative and the integer sentinel are rejected; 2^63 + (2^63−1) (= sentinel) throws, (2^63−1)·2 does not; `mul_pow2` throws exactly at overflow; `SignedPair` order, equality and subtraction equal 128-bit integer arithmetic on 2·10^5 random pairs |
| Sec 3.1 Lemma 3.1 (selection) | `primitives/search.hpp: select_bfprt` | PASS | E6: 30,000 seeds, sizes up to 5000, many ties, against `std::sort` |
| Sec 3.1 Lemma 3.1 (DFS/BFS/toposort, counting sort) | `search.hpp` | PASS (existing tests only) | not fuzzed independently; covered indirectly by the Sec 7 pipeline (`topological_sort` of the marker DAG, counting sorts) |
| Sec 3.1 Lemma 3.2 (Dijkstra, radius cutoff) | `primitives/dijkstra.hpp` | PASS | E6: heap and scan versions against Bellman–Ford on 30,000 seeds with zero weights, loops and multi-edges, radius in {∞, 0..79}: settled set is exactly {d ≤ R}, every parent edge is tight |
| Sec 3.2 Alg 1, Lemma 3.4–3.5 (three passes) | `primitives/mst.hpp` | PASS | E6: exact edge-set equality with an own Kruskal under (length, index) on 30,000 seeds (n up to 1500, ties, loops, multi-edges); B1 = max(2, ⌈log2 L_k⌉), B2 = L_k checked in `three_pass_limits` against l.397 |
| Sec 3.2 Lemma 3.3 (Fibonacci heap) | `fib_heap.hpp` | PASS (existing tests only) | degree array size 3 + 2⌈log2(s+1)⌉ as in the proof (l.433); exercised through the MST fuzz |
| Sec 3.3 Lemma 3.6 (disjoint sets) | `dsu.hpp` | PASS (existing tests only) | exercised through MST, calendars and atoms |
| Sec 3.4 Lemma 3.7–3.9 (LCA, level ancestor) | `forest_queries.hpp` | PASS | E6: LCA, depth, level ancestor and the out-of-range rejection against naive parent walks on random, path-like and caterpillar forests with shuffled labels (n up to 3000, 30,000 seeds); B = max(2, ⌊log2(k+1)/8⌋) as in l.505 |
| Sec 3.5 Lemma 3.10–3.11 (split–findmin) | `split_findmin.hpp: SplitFindmin` | PASS | E6: splits (including rejected repeated boundaries), decreases and queries against a naive array with ties broken by index, uint64 and double keys with infinite entries, k up to 700, 30,000 seeds |
| Sec 6 Def 6.1, Alg 5 l.2–3 (u, g, P) | `pivots/name_routes.cpp` | PASS | E7: P recomputed from l.2–3 in the harness equals `p_log2` on every Mode A run |
| Sec 6.1 Alg 5, Lemma 6.8 (certificates) | `name_routes` | PASS | E7: Mode A, 60 instances with t = 66..105 and density 50–100 % (31 of them sample, so Alg 5 l.5–14 runs), uint64 and non-dyadic double: every λ_e is an input weight, λ_e ≤ a_e, the dictionary is strictly increasing, and d_GP(x, y) ≤ Pλ_e against an own Floyd–Warshall; no rejected trial, no fallback record. Mode B (D-9, D-10, D-17: sample factor 1–2, residual factor up to +2, P = 2^0..2^3, H = 2^4..2^13, retry cap 16), 2000 instances, 21,358 trials, 20,763 rejections, 1236 cap hits: every certificate valid, and exactly one `FallbackRecord` per cap hit |
| Sec 6.2 Alg 6, Lemma 6.4–6.7 | `pivots/estimate.cpp` | PASS (indirect) | code review of T, ℓ, H (l.1473, l.1490–1491: `tau = ⌈log2 4t⌉`, `ell = 2 tau + ⌈log2 t⌉`, H by big-integer comparison) and of the closed radius test H d ≤ (H + j)Δ (l.1493–1494). Z ≥ d_S was checked only through certificate validity with small P in Mode B (a Z below d_S would have failed there); not checked against the table directly |
| Sec 6 numeric range | `estimate.cpp` | Needs agreement | F-2 (uint64 sample path ≥ 2^64), F-3 (double H·d beyond the double range) |
| Sec 7.1 l.1646–1652, Lemma 7.1 | `owner/aux_metric.cpp` | PASS | code matches A_e = R_x + w + R_y, d_e = Pλ by doublings (saturating), radial iff 2d_e ≤ w (equality radial), intraowner radial with d_e = 0; E4 pipeline runs below |
| Sec 7.2 Lemma 7.2 (G', T', J) | `aux_metric.cpp: build_aux_graph` | PASS | E4: J of Alg 7 equals the own exact Dijkstra on G' on every run |
| Sec 7.3 Lemma 7.3–7.6 (pieces, grids, calendars, marker arcs, atoms) | `owner/calendars.cpp` | PASS | code review: overlap of closed intervals as r_u ≤ 16R_v ∧ r_v ≤ 16R_u (l.1693–1695); search windows 4 + ⌈log2 h⌉ and 1 + ⌈log2 h⌉ follow from r_v ≤ 16 h r_u and R ≤ h r (l.1730–1734); q = grid[i − 5] gives 32q ≤ r < 64q; U⁺ is the first grid value > R; marker arc U⁺_u → 16 q_v and its existence in both calendars (cells ≥ 32); defining edges 2^7 w' < 2 max(q) on nonroot edges. Values validated by E4 |
| Sec 7.4 Alg 7, Lemma 7.7–7.9 | `owner/owner_distances.cpp` | PASS for values; FAIL for the overflow contract (F-1) | E4: uint64, 24,892 runs (n ≤ 400, 5 graph families, 5 weight ranges including powers of two and [1, 2^40], b ∈ {1, 2, 3, 4, ⌈log2 N⌉, N, N+1}, exact names with P = 1 and NameRoutes with Mode A and Mode B parameters, Alg 4 Mode A and Mode B `{1,2,2}`, `{1,2,5}`; 5833 multi-vertex atoms, 56 radial-branch atoms): 0 wrong J. Double (non-dyadic 52-bit mantissas), 22,676 completed runs: 0 values outside (hops + 1) ulp of the exact rational distance (max 4.04 ε·J); 253 runs end in the documented `std::logic_error` of Q-55 option A, none with power-of-two weights |
| D-20 (uint64 range 2^10·max R < 2^64) | `calendars.cpp`, `owner_distances.cpp`, `aux_metric.cpp` | FAIL | F-1: an `overflow_error` inside the stated range; F-5: wording |
| D-21 (2^7 scaling and halvings) | `owner_distances.cpp` | PASS | subnormal weights (units 4.9e−324 .. 2.2e−308) return exact J; integer J is checked for zero low bits as documented |
| D-22 (gap as (R_i − R_l) + (R_j − R_l)) | `owner_distances.cpp: build_families` | PASS | implemented as documented; R_l ≤ R_i because R is nondecreasing down T' |
| D-23 (saturating bookkeeping) as used by Sec 7 | `owner_distances.cpp` | PASS (Sec 7 side) | the atom deadline z + L_Y, the proposal max{b0, z + L_i + g} and the exit offers use `sat::add` and are then capped by R_v, as D-23 requires of the evaluator; evaluator internals are out of scope |
| D-24 (double slot rounding → `domain_error`) as reached from Sec 7 | `PortionedEvaluator` via `owner_distances` | Needs agreement (documentation) | F-4. Behaviour agrees with the D-24 text: the 44 rejected inputs (spokes 2^53 to 2^200 under unit-size atom edges) all have sums that are not exact in double; the 87 accepted inputs returned the same J as a double Dijkstra on G' |
| D-9, D-10, D-17 (Mode B NameRoutes) | `name_routes.cpp` | PASS | E7 Mode B above |
| D-2, D-18 | `model/`, `split_findmin.hpp` | PASS | D-2 not used (no clz); D-18 resolved (`BlockLayout`) |

## 2. Counterexample

### F-1 (FAIL): Alg 7 throws `std::overflow_error` for a G' edge that cannot be relaxed, inside the D-20 range

- Input (uint64, 14 vertices, 15 edges, separator size b = 2, route names exact with P = 1, which
  satisfy Definition 6.1):
  - path 0–1–…–11 of 11 edges of weight 2^53;
  - edge 0–11 of weight 11·2^52;
  - edges 0–12 and 11–13 of weight 1;
  - edge 12–13 of weight 11·2^53.
- Pipeline: `owner_decomposition(g, 2)`, `projected_interowner_edges`, exact names (least instance
  weight λ ≥ d_GP, λ ≤ a_e), `classify_outer_edges`, `build_aux_graph`, `owner_distances`.
  No randomness is involved.
- Observed: 7 pivots; G' has 8 vertices; max R = 2^53, so 2^10·max R = 2^63 < 2^64 (inside the
  D-20 range); the chord 12–13 is radial with w' = 11·2^53 + 11·2^52 ≈ 2^57.04. `owner_distances`
  throws `std::overflow_error("calendars: 2^7 x does not fit the weight type")` from
  `calendar_scale(ed.weight)` in the edge scan (Alg 7 l.6).
- Expected: J is returned. J of vertices 12 and 13 is 1, every other J is at most 2^53, and even
  2^7·J fits. The chord has w' > R of both endpoints, so Alg 7 discards every offer over it
  (l.1821–1823); dropping it, or saturating its scaled length, would be exact.
- Rule violated: the overflow rule of `include/sssp/paper.hpp` (an auxiliary quantity, here a
  scaled edge length, must not raise) and the D-20 statement that inputs with 2^10·max R < 2^64
  are accepted. D-20 lists L_Y and A_e as further members of the limit but not 2^7 w'. The
  function comment of `owner_distances` does mention 2^7 w', so the code and that comment agree;
  the documents above do not.
- Oracle: own exact Dijkstra on G' (gives the J above). Suspected step: Alg 7 l.6 together with the
  D-21 scaling (Lemma 7.3 l.1723–1724).
- Size sweep: K = 10 (max w' ≈ 2^56.9) is accepted with exact J; K = 11 and K = 12 throw; b = 1 and
  b = 3 do not produce the radial chord. With the paper's NameRoutes parameters the chord is macro
  (P is large), so the case needs a small P (exact names, or Mode B D-17).
- Not on the `paper_sssp` path (Mode A never takes the outer branch). Suggested regression name:
  `regression: P4-REVIEW owner_distances scaled chord overflow`.

## 3. Points that need agreement

- **F-2 (Sec 6.2, uint64).** `estimate.cpp` `first_covering_centers` forms d + w with
  `checked_add` and throws when a sample path reaches 2^64. Reproduction: complete instance on
  t = 70 pivots, weights 2^63 + (x mod 1000) from the xorshift sequence seeded with 99, `WordRng(3)`,
  Mode A parameters → `std::overflow_error("integer distance overflow")`; with base 2^62 the same
  instance succeeds and every certificate is valid. The value is a path length, but in the sample
  of the pivot-routing graph, not in G. `estimate.hpp` l.122 documents it. Is d_S a "tentative path
  length" (may throw) or an auxiliary quantity (must not throw) under the overflow rule? In Sec 7,
  A_e ≥ 2^63 needs heavy interowner edges, which D-20 does not exclude.
- **F-3 (Sec 6.2, double).** For double weights the radius test forms H·d (and T²w) by doublings in
  `double` and throws `std::overflow_error("scaled length overflows double")` once it leaves the
  range: complete instance on t = 70, weights 1e300·(1 + U·2^−52) with U from the xorshift
  sequence seeded with 12345, `WordRng(7)`, Mode A parameters → throws, although every distance
  is at most 2e300; 1e280 succeeds. This contradicts the Q-06 decision ("power-of-two factors are
  held as (value, exponent) pairs") for double, and `estimate.hpp` documents the throw for integer
  weights only. Near DBL_MAX, the same scope remark as D-30 may apply; agree whether to document it
  or to compare in (value, exponent) form.
- **F-4 (Sec 7 documentation).** D-24's `std::domain_error` reaches `owner_distances` for double
  inputs whose atom slots round back to their start. Minimal example: single pivot (b ≥ N), path
  0–1 of weight 2^53, then 1–2–3 of weight 1, chord 1–3 of weight 2. Neither the throw list of
  `owner_distances` (header l.180–183) nor the Sec 7 rows of `PAPER_MAPPING.md` mention D-24. No
  wrong value results.
- **F-5 (D-20 wording).** D-20 says that inputs with R ≥ 2^54 are not accepted. Single edges of
  weight 2^54 and 2^54 + 1 (single pivot) are accepted with exact J, because the actual test is
  whether the top grid point origin·2^j ≥ 2^9 max R fits. The path with weights 1 and 2^54 throws.
  Suggest stating 2^10·max R < 2^64 as sufficient, not necessary. This is harmless.
- Observation for the `local/` reviewer (not in this scope): Mode B `LocalHierarchyParams{1, 2, 2}`
  and `{1, 2, 5}` admit atoms of 26–31 vertices (k > L_k^2), and `local_hierarchy` then throws
  `std::invalid_argument` (Sec 5.4 needs N ≥ 32). `PAPER_MAPPING.md` documents this; it showed up
  in 92 of 20,000 E4 runs.

## 4. Checks run

Build: `cmake --preset debug && cmake --build --preset debug -j4` (Apple Clang, `-fno-fast-math`,
`SSSP_DIAGNOSTICS=1`). Harnesses: `-O2 -std=c++20 -fno-fast-math` linked against
`build/debug/libsssp.a`.

| Id | Command (in `.cache/review/`) | Size / seeds | Result |
|---|---|---|---|
| — | `ctest --test-dir build/debug -j4 -R "owner\|aux\|calendar\|pivot\|ilog\|length\|param\|random\|word\|primitive\|name_routes\|estimate\|mst\|dsu\|lca\|level_ancestor\|split_findmin\|fib\|paper_dijkstra\|search\|sort\|select\|merge_events\|edge_pairs\|growth\|las_vegas"` | 327 tests | 100 % passed |
| E1 | `./thresholds > thr.out && python3 -I thresholds_check.py thr.out` | 453,232 checks | 0 mismatches |
| E2/E3 | `./model_probe` | see §1 | 0 failures |
| E4 | `./owner_fuzz u64 20000 1000`; `NMAX=400 ./owner_fuzz u64 5000 50000` | 25,000 seeds | 24,892 runs, 0 wrong J, 0 certificate violations; 108 runs stop in Mode B Alg 4 (`invalid_argument`, see observation) |
| E4 | `./owner_fuzz dbl 20000 1000`; `NMAX=400 ./owner_fuzz dbl 3000 90000` | 23,000 seeds | 22,676 runs, 0 values outside (hops + 1) ulp; 253 documented `logic_error` (Q-55) |
| E5 | `./scaled_edge_overflow`; `./range_probe` | K = 10..12, b = 1..3; 22 range cases | F-1, F-5 |
| E6 | `./prim_fuzz 30000` | 30,000 seeds × 6 primitives | 0 failures |
| E7 | `./names_fuzz a 60 1`; `./names_fuzz b 2000 1` | 2060 instances | 0 violations |
| E8 | `./dbl_probe` (spokes 2^20..2^200 and 2^48..2^52) | 132 cases | 87 equal to a double Dijkstra on G', 44 `domain_error` (D-24), 1 `logic_error` (Q-55), 0 wrong |
| E9 | `./edge_probe`; `./u64_big_names` | subnormal and huge weights | subnormal exact; F-2, F-3 |

## 5. Not tested

- No sanitizer build (machine load); the harnesses ran on the Debug library only.
- Alg 6 internals (bands, quotients, partition trees, the E[Z] bound of Lemma 6.7) were not checked
  directly; only through certificate validity.
- Fibonacci heap, disjoint sets, merge events, edge pairs, DFS/BFS/toposort and counting sort were
  not fuzzed with own oracles; the existing tests pass and the Sec 7 pipeline exercises them.
- The radial branch of Alg 4 inside atoms ran on 93 atoms only (Mode B parameters); Mode A
  parameters always take the base case (Q-04).
- Running time and the work bounds of Lemmas 3.x, 6.6–6.8 and 7.3–7.9 were not measured.
- `PortionedEvaluator` and Alg 4 were black boxes (out of scope).
