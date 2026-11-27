# P4-REVIEW: independent correctness review of `src/paper/hierarchy/` (Sec 4)

- Paper: arXiv:2611.34911 v1, Sec 2 (length model, l.292-296), Sec 4 (Def 4.1-4.2, Sec 4.1
  Lemma 4.3, Sec 4.2 Lemma 4.4-4.5, Sec 4.3 Alg 2, Alg 3, Lemma 4.6). Cited by section, lemma
  and line only.
- Code under review: `src/paper/hierarchy/{ordered_hierarchy,normalize,child_order,evaluate,
  reference_hierarchy}.{hpp,cpp}` at `552b3de`.
- Out of scope: `PortionedEvaluator` (Sec 7.4, Lemma 7.7), the I1-I8 property tests of Evaluate
  and the Lemma 4.7 work bounds (other reviewers), `src/paper/outer/` and Sec 8.2.
- Reviewer: correctness reviewer (independent of the implementer). No implementation file and
  no test was changed. The throw-away checkers live in `.cache/review/` (git-ignored, not
  committed).

## Verdict

**FAIL.** Sec 4.1, Sec 4.2 for integers, the Def 4.1 / 4.2 checker and the reference
hierarchy pass every check below. Two major failures were found on valid input:

- F-1 (`double`): Evaluate returns wrong distances for inputs with integer weights. These inputs
  are inside the decided `double` scope (Q-55 option A: "integer and dyadic weights"). The
  cause is the materialised forest gaps of Sec 4.2 (D-15 / D-27). The minimal counterexample is
  a 4-vertex path, and it also fails through the Mode A Alg 4 base case (`local_hierarchy` +
  `evaluate`).
- F-4 (`std::uint64_t`): Evaluate throws `std::logic_error` ("origin rule left the parent's slot
  (Def 4.1 violated)") on a valid 3-state hierarchy whose distances are all small. The cause is
  that the origin rule of Sec 4.3 l.794 compares two saturated values (D-23).

Two further failures need invalid input only: F-2 (lengths that are accepted without an error)
and F-3 (a hang).

## Verdict per item

| Item | Paper | Verdict | Evidence (details below) |
|---|---|---|---|
| Ordered hierarchy, checker `validate_ordered_hierarchy` / `validate_seeded_ordered_hierarchy` | Def 4.1 (l.646-658), Def 4.2 (l.659-672) | PASS | No false positive on about 0.87 M generated inputs (C1, C4; every generator input that failed validation was traced to a generator overflow in u128 -> u64). Each of 6 one-defect mutations (H1, H2, H3, H4, H5, infeasible offer) was reported under its condition on 5 000 seeds x 2 types (C3) |
| Builders `make_leaf_forest`, `add_internal_node`, `collect_roots` | Def 4.1 | PASS | Used by every generated input (C1) |
| Normalisation: masses S(C) | Sec 4.1 l.681-684 | PASS (integer, exact doubles) | S of every surviving node equals an independent 128-bit recomputation, saturated as D-23 says (C1, `u64-tall`). For wide-range doubles S is a rounded sum (132 differences of rounding size in C4), as Q-06 allows |
| Normalisation: links 4S(C) >= h(P), unit tables, nu in (h/2, h] | Sec 4.1 l.684-689, Lemma 4.3 proof l.697-709 | PASS | nu(C) of every input node equals an independent recomputation from the paper's rule (link test, least height of the linked set, doubling table, round down) (C1) |
| Normalisation: contraction, units at least double, dyadic or 2S(C) < nu(P) | Lemma 4.3 (l.690-695), l.710-718 | PASS | Checked on every contracted forest. The contracted forest passes the Def 4.1 checker with heights nu and c0' = 2c0 (Q-50) (C1) |
| Normalisation: self-arcs, root quotient, light arcs at the new LCA | l.719-723 | PASS for valid input; FAIL (minor) on invalid lengths | See F-2: a root arc with a negative or NaN length is not rejected, although `normalize.hpp` says that it raises `std::invalid_argument` |
| Child order: topological order, representatives | Sec 4.2 l.726-730, Lemma 4.5 | PASS | Every light arc (naive LCA, length < nu) goes from a lower to a higher rank. Every rep(i) is a state below i (C1) |
| Child order: gaps g(i, j), integer | Lemma 4.4 (l.731-753), Lemma 4.5 (l.758-779); D-25, D-26 | PASS | g >= exact physical distance between the representatives (APSP), and sum of g <= 2 sum S(i) + 2 d nu(X) + 2 S(X) at every node. This also holds with root distances up to about 2^63 (`wide` mode, 20 000 seeds), and with 128-bit charges and saturation at infinity (C1, C2) |
| Child order: gaps g(i, j), double | same; D-15, D-27 | **FAIL** | F-1: rd(b) - rd(l) cancels catastrophically once root distances exceed 2^53. g falls below the physical distance by up to half an ulp of rd, which is unbounded relative to g |
| Evaluate / Scan, integer: Lemma 4.6 | Sec 4.3 l.781-864, Alg 2, Alg 3, Lemma 4.6 | PASS | Every protected state emits delta exactly, every emission is >= delta and each state emits at most once. Tested on about 0.37 M diverse integer hierarchies (0.87 M with double), with supplemental offers on and off (D-14), and also near the top of uint64 (C1, C2) |
| Evaluate, integer: origin rule with saturated values | Sec 4.3 l.794-800; D-23 | **FAIL** | F-4: `std::logic_error` on valid input when both D0 + S(C) and the parent's slot end b exceed the range. Found in 6 of 100 000 `u64-tall` runs and 5 of 2 829 `big` `u64-tall` runs, and reduced to 3 states |
| Evaluate / Scan, double | same | **FAIL** | F-1. Lemma 4.6 fails at protected states: in the minimal case 4 is returned where delta = 5. The inputs are integer-valued and every distance is below or equal to 2^53 |
| Overflow rule (`include/sssp/paper.hpp`) | Sec 2; D-23 | PASS for `std::overflow_error`; FAIL through F-4 | Every `std::overflow_error` (1 372 in `u64-big`, 69 on the reference hierarchy with 2^60 weights) came after a relaxation beyond the range. In the 69 reference cases the exact distance was >= 2^64 - 1. Saturated masses, slot ends and gaps never threw on their own (`dbl-tall`: 0 exceptions). In 6 of 100 000 `u64-tall` runs and 5 of 2 829 `big` `u64-tall` runs, a saturated expiry together with a saturated slot end led to the `std::logic_error` of F-4, which the rule ("MUST NOT throw") excludes (C1, C2) |
| Fallback records | CLAUDE.md fallback policy | PASS | Sec 4 has no direct search; `diag.fallbacks` stays empty. The D-19 origin choice for double is a counter only (needs agreement, A-3) |
| Reference hierarchy `build_threshold_scc_hierarchy` | supplement (02 sec 2.11), Def 4.1 | PASS | 60 000 outputs (3 weight scales x 20 000 seeds) pass the Def 4.1 / 4.2 checker with c0 = 1, and heights strictly increase towards the root (C2). Zero lengths are rejected as documented |
| Boundaries | - | PASS except F-2, F-3 | N = 0, N = 1 with a self-loop and two offers, zero-length arcs and links, copies of one vertex, ties (`ties` mode), lengths near 2^63, heights up to 2^60 (uint64), subnormal doubles (x 2^-1070), doubles x 2^1000, masses beyond DBL_MAX |
| Deviations D-14, D-19, D-24, D-25, D-26 | - | PASS (as documented) | Details under "Deviations" |
| Deviation D-23 | - | **FAIL** | The saturations are implemented as listed. The claim "with a saturated expiry ... the choice of the rule is the paper's" (`evaluate.cpp`, `open`) fails when b is saturated too (F-4) |
| Deviations D-15, D-27 | - | **FAIL** | Both deviations are implemented as documented. Their stated impact ("a few ulp", "no wrong result is returned silently") is false (F-1) |

## Counterexamples

### F-1 (major): Sec 4.2 gaps in double give wrong protected distances

Minimal input (found by `scc_wide`, minimised by edge deletion and weight reduction;
`seed 327`, generator `.cache/review/scc_wide.cpp`, big-weight shift 60):

- physical graph: path with the edges 0-1 weight 2^53 = 9007199254740992, 1-2 weight 2 and
  2-3 weight 3, one state per vertex;
- state graph: both directions of every edge; hierarchy `build_threshold_scc_hierarchy` (c0 = 1);
- one initial offer: state 1 with value 0 (source vertex 1, offset 0).

| vertex | exact (Dijkstra in 128 bits) | `evaluate<uint64_t>` | `evaluate<double>` | `evaluate<double>`, `supplemental_offers = false` (D-14) |
|---|---|---|---|---|
| 0 | 9007199254740992 | 9007199254740992 | 9007199254740992 | 9007199254740992 |
| 1 | 0 | 0 | 0 | 0 |
| 2 | 2 | 2 | 2 | 2 |
| 3 | **5** | 5 | **4** | 5 |

Every state is protected (its tight witness is the shortest path). All four distances and every
sum Evaluate forms along the shortest paths are exactly representable.

The same graph through the Mode A producer: `local_hierarchy(g)` with `kPaperLocalParams` (Alg 4
base case, `FallbackRecord` present, c0 = 1, the hierarchy passes the Def 4.1 / 4.2 checker),
then `evaluate`. The minimum per vertex from source 1 is `2^53 0 2 4` for double and `2^53 0 2 5`
for uint64 (`.cache/review/local_path.cpp`).

Trace (expected from Sec 4.1-4.3, confirmed by the output):

1. Sec 4.1: S(A{1,2}) = 4 is linked to B (16 >= 3), and A and B have the same unit 2, so A is
   contracted. The result is a node X = {1, 2, 3} with nu = 2 and S = 10. Its parent, the root
   with h = 2^53, is not linked to it. X has no light arcs (no length < 2).
2. Sec 4.2: F is the path 0-1-2-3, rooted at state 0. rd(2) = 2^53 + 2 is exact, but
   rd(3) = fl(2^53 + 5) = 2^53 + 4 (round half to even). So g(3, 2) = rd(3) - rd(2) = 2, while
   the certified walk is 3. The ranks at X are 3, 2, 1.
3. Sec 4.3, slot [2, 4) of X: the run [3, 2] has its minimum at state 2 (z = 2). Registering 2
   posts the supplemental offer max{b, z + S(2) + g(3, 2)} = max{4, 2 + 0 + 2} = 4 at state 3.
   delta(3) = 5, so the offer is infeasible. In slot [4, 6) state 3 emits 4.

Systematic evidence (no threshold-SCC special case):

- Reference hierarchies (c0 = 1, no copies), random graphs with n in 3..10 and 20 % of the
  weights in (1..8) x 2^56 (or x 2^60), the rest in 1..300; 20 000 seeds each. uint64: 0 wrong.
  double: 405 (2^56) and 300 (2^60) seeds have a wrong distance below 2^53. Examples: exact 73
  and got 40; exact 118 and got 84. The same generator with weights x 2^20 (all sums below
  2^53): 0 wrong.
- Generated hierarchies (`hier_fuzz ... wide`, root heights and root links x 2^56): 2 000 seeds.
  double: 85 protected states wrong, all with supplemental offers on; with D-14 off, 0 wrong.
  There are 213 emissions below delta and 1 716 gaps below the physical distance. uint64 on the
  same generator (20 000 seeds): 0 wrong.
- Also exact with c0 = 1 and copies (`.cache/review/probe.cpp gapround-c01`): the protected
  state 3 is expected to emit 327 but emits 254 (double); uint64 gives 327.

Suspected step: Sec 4.2, Lemma 4.5 proof l.767-769 ("One lowest-common-ancestor query and root
sums give each F-distance"). The paper keeps these differences as signed pairs of sums (Sec 2
l.292-296), exactly. D-15 instead materialises (rd(a) - rd(l)) + (rd(b) - rd(l)) from rounded
root distances. Once rd exceeds 2^53 the difference loses every bit below ulp(rd), and the
padding to b in Sec 4.3 l.826-831 no longer covers the loss.

Why the existing tests miss it: `test_evaluate_vs_dijkstra.cpp` restricts double to integer
lengths whose sums stay below 2^53. The double cases of `test_child_order_range.cpp` scale
every length by one power of two, which keeps them exact, or check +inf and never NaN. None of
them mixes lengths of very different size on one root path of F. The near-2^53 regression case
`20261127-evaluate-double-slot-rounding` concerns slot ends (D-24), not gaps.

Note for the lead (not a fix): the loss comes from measuring rd from the root of the whole F
tree. Any representation whose error is relative to the distances inside F[X] would remove it,
for example root distances measured inside each node, exact signed pairs, or a forest distance
summed along the F path.

### F-2 (minor, invalid input only): root arcs with negative or NaN lengths are not rejected

`normalize.hpp` says: "Lengths that are negative or NaN raise std::invalid_argument from
model::Len". This holds for arcs below an LCA, which are compared with nu, but not for arcs
between two roots, which go to the root quotient without a check (l.720). Self-arcs are not
checked either.

| Input (double; two singleton roots, arc 0 -> 1 of length L, offer 5 at state 0) | Result |
|---|---|
| L = -1.0 | returns value[1] = 4 (no error) |
| L = NaN | returns value[1] = +inf (no error; the arc is ignored) |
| L = +inf | `std::overflow_error` (not `std::invalid_argument`) |
| the same arc below a common internal node | `std::invalid_argument` (as documented) |

`validate_ordered_hierarchy` reports all three (kArcLengthInvalid). The defect is only that
`evaluate` accepts them without an error. The uint64 sentinel 2^64 - 1 as an arc length gives
`std::invalid_argument` below an LCA but `std::overflow_error` between roots (A-4).

### F-3 (minor, invalid input only): evaluate does not return on non-finite link lengths

Input: one root with height 1 and two leaf children, a link (0, 1) of length L, no arcs, and an
offer 0 at state 0, so state 1 is never reached.

| L | Result |
|---|---|
| NaN (double) | does not return within 10 s |
| +inf (double) | does not return within 10 s |
| 2^64 - 1 (uint64 sentinel) | does not return within 10 s |
| -3.0 (double) | `std::logic_error` "root distance decreases along F" (not `std::invalid_argument`) |

Cause, read from the code: the mass sums use `sat::add` without validating link lengths, so S
becomes NaN or +inf. The root loop of Alg 3 l.10 then steps through width-1 slots until a slot
end saturates (uint64: about 2^64 slots) or rounds back to its start (double: about 2^53 slots,
then D-24). The checker rejects these links (kLinkLengthInvalid); the defect is that the
documented `std::invalid_argument` is not raised. A finite but too long link (> c0 h, an H5
violation) leads to the same long loop: with uint64, h = 1 and a link of 2^40, evaluate does not
return within 10 s. For valid input the loop is bounded:
S(root) / nu(root) < 2 (N - 1)(1 + c0).

### F-4 (major): uint64 origin rule with a saturated expiry and a saturated slot end

Found by C1 in `u64-tall` (seeds 17863, 36023, 50471, 68715, 83047, 88991 of `hier_fuzz 100000 1 u64-tall`) and in `big` `u64-tall` (seeds 300457, 300935, 303949, 304394 and 304876 of
`hier_fuzz 5000 300001 u64-tall big`; seed 300457 has 51 states and 86 input nodes). It was reduced
by hand to the following input (`.cache/review/probe.cpp satorigin`). The input passes
`validate_seeded_ordered_hierarchy`:

- uint64, c0 = 1; physical vertices 0, 1, 2 with the edges 0-1 (1) and 1-2 (5 * 2^57); states
  0, 1, 2 with pi = identity; no arcs;
- C = {1, 2}: h = 5 * 2^57, link (1, 2) of length 5 * 2^57, so S(C) = 5 * 2^58;
- R = {0, C} (the root): h = H = 2^64 - 2^59, link (0, 1) of length 1;
- source (0, 0); initial offers: state 0 at 0 and state 1 at H (feasible: delta(1) = 1).

Expected (Sec 4.1, Sec 4.3): 4 S(C) = 5 * 2^60 < H, so C is not linked, nu(C) = 5 * 2^57 and
nu(R) = H. The root slots are [0, H) and [H, 2H). C opens in the second slot at D0 = H, and the
exact expiry H + 5 * 2^58 is below b = 2H, so l.794 takes the first rule (origin D0). The
protected state 0 emits 0 and no exception is thrown.

Actual: `std::logic_error`: "evaluate: origin rule left the parent's slot (Def 4.1 violated)".

Cause, read from `evaluate.cpp` (`open`): the second root slot end H + H saturates, so
b = infinity(). The expiry D0 + S(C) = H + 5 * 2^58 also saturates. The test `!less(expiry, b)`
then compares two saturated values and takes the second rule. In the saturated branch the
subtraction steps r = nu(R) - j nu(C) leave the parent's slot, because nu(C) does not divide
nu(R) (an unlinked boundary). The comment that a saturated expiry "exceeds every finite b" holds
only for a finite b. The order of D0 + S(C) and b = a + nu(P) is still decidable with terms that
fit W. Because D0 >= a, the condition D0 + S(C) >= b holds exactly when S(C) >= nu(P), or
D0 - a >= nu(P) - S(C) (a saturated S(C) counts as >= nu(P)). This is a note for the lead, not
a fix.

Suspected step: Sec 4.3 l.794-796 (choice between the two origin rules) under D-23.

## Checks run

Build: `cmake --preset debug -DSSSP_BUILD_BENCHMARKS=OFF`, then `cmake --build --preset debug
--target sssp sssp_unit_tests sssp_regression_tests sssp_differential_tests -j4` (Apple clang 17,
Debug, diagnostics ON, no fast-math). The checkers were compiled with `clang++ -std=c++20 -O1/-O2`
against `build/debug/libsssp.a`.

Existing tests: `ctest --test-dir build/debug -j4 -R "hierarch|evaluate|Evaluate|normali|
Normali|child|Child|order|Order|20261127-(evaluate|child)"`. Result: 178 / 178 passed
(unit 129, differential 27, differential-slow 11, regression 11; 157 s).

New, throw-away checkers (`.cache/review/`, not committed; the oracles share no code with
`src/` except the Def 4.1 checker, which is used only to confirm that the generator's inputs are
valid):

- C1 `hier_fuzz.cpp`. The generator builds random Def 4.1 hierarchies top-down. It is not
  threshold-SCC: heights are equal to the parent (contraction), slightly lower or much lower
  (unlinked boundaries); c0 is in 1..3; ranks are random; links have lengths in [0, c0 h] and are
  backed by physical edges no longer than the link (zero-length and copy certificates included);
  arcs are heavy in both directions; light arcs go forward and root arcs follow a random root
  order; 1-2 weighted sources; feasible offers, many of them tight. The oracles are
  Floyd-Warshall APSP and delta in 128 bits, the protected set as a tight-arc fixpoint,
  independent masses and units (Lemma 4.3 rule), the Lemma 4.3 relations, the light-arc order,
  rep below the child, gap >= APSP, the Lemma 4.4 sum bound, and Lemma 4.6 with supplemental
  offers on and off. Default sizes: <= 8 physical vertices and <= 18 states.
- C2 `scc_wide.cpp`: threshold-SCC reference hierarchies on random graphs, evaluated from one
  source and compared with an exact Dijkstra in 128 bits.
- C3 `hier_fuzz ... mut`: one-defect mutations of valid inputs, each of which must be reported
  by the checker under its condition.
- C4 `hier_fuzz ... wide`: root heights and root links x 2^56; arc lengths and offers are
  rounded up to 53 significant bits so that the double input equals the exact one.
- Probes `probe.cpp` (invalid lengths, -0.0, N = 0 and 1, F-1 constructions) and
  `local_path.cpp` (F-1 through `local_hierarchy`).

| Run | Seeds | Valid inputs run | Failures |
|---|---|---|---|
| C1 `u64` | 1..100 000 | 100 000 | 0 |
| C1 `u64-big` (lengths and heights x 2^52..2^60) | 1..100 000 | 82 333 (17 667 skipped: a scaled length does not fit) | 0 (1 372 overflow_error, each after a relaxation beyond the range) |
| C1 `u64-tall` (heights x 2^55..2^58, masses saturate) | 1..100 000 | 100 000 | 6: F-4 (seeds 17863, 36023, 50471, 68715, 83047, 88991) |
| C1 `dbl` | 1..100 000 | 100 000 | 0 |
| C1 `dbl-subnormal` (x 2^-1070) | 1..100 000 | 100 000 | 0 |
| C1 `dbl-huge` (x 2^1000) | 1..100 000 | 100 000 | 0 |
| C1 `dbl-tall` (x 2^985, heights x 2^30: masses beyond DBL_MAX) | 1..100 000 | 100 000 | 0 |
| C1 `ties` (<= 30 states, heights <= 4, weights <= 2), all 7 modes | 200 001..220 000 | 139 963 (7 modes x 20 000, 37 skipped because a length does not fit) | 0 |
| C1 `big` (<= 25 vertices, <= 60 states), all 7 modes | 300 001..305 000 | 30 552 (the generator wrapped heights in 2 171 `u64-tall` seeds; those inputs failed validation and were not run) | 5: F-4 (`u64-tall`) |
| C3 mutations, uint64 and double | 1..5 000 | 10 000 | 0 missed |
| C4 `wide`, uint64 | 1..20 000 | 19 999 | 0 (1 generator-side invalid input, 5 legitimate overflow_error) |
| C4 `wide`, double | 1..2 000 | 2 000 | F-1 (85 protected states wrong, all with supplemental offers on) |
| C2, weights x 2^20 | 1..20 000 | 20 000 | 0 (uint64 and double) |
| C2, weights x 2^56 | 1..20 000 | 20 000 | uint64 0; double 405 seeds (F-1); 889 `domain_error` (D-24) |
| C2, weights x 2^60 | 1..20 000 | 20 000 | uint64 0 (69 overflow_error, all with an exact distance >= 2^64 - 1); double 300 seeds (F-1) |

The counts of the generated structure (C1, 100 000 seeds per mode) show that the paths are
exercised: a contracted node in about 84 % of the inputs, the second origin rule in about 44 %
of the Evaluate runs (two runs per seed), and a supplemental offer that decreased a key in about
94 % of the runs with supplemental offers on.

## Deviations (PAPER_MAPPING D table) in this area

| D | Documented | Observed | Verdict |
|---|---|---|---|
| D-3, D-4, D-5 | not used | `SplitFindmin` (Lemma 3.11), the O(1) cyclic-start table and `ForestQueries` (Lemma 3.9) are used | PASS |
| D-12 (-0.0 -> +0.0, candidate) | candidate | an offer -0.0 is accepted and emitted as -0.0 (sign kept); arcs of length -0.0 propagate -0.0 | needs agreement (A-5) |
| D-14 | test-only switch; protected emissions unchanged | protected emissions identical on every exact-mode run; with the switch off, F-1 disappears | PASS |
| D-15 | one nonnegative gap; "a supplemental offer may be infeasible by a few ulp" | the gap can be smaller than the walk by up to ulp(rd)/2, and rd can be much larger than g, so the offer can be infeasible by much more than an ulp of delta | **FAIL** (F-1) |
| D-19 | b - j nu by successive subtractions; double falls back to D0, integer throws logic_error | no logic_error outside the saturated case of F-4 (about 0.87 M inputs); the second rule is used in about 44 % of the Evaluate runs; the integer logic_error is reached through the saturated values of D-23 (F-4) | PASS (the double fallback path itself is not exercised by these exact inputs) |
| D-23 | saturate S, expiry, slot ends, offers and relaxations; throw only for a state that received only beyond-range values; "with a saturated expiry ... the choice of the rule is the paper's" | as documented in `dbl-tall`, `u64-big` and C2, and in all but 6 of the 100 000 `u64-tall` runs. With a saturated expiry AND a saturated parent slot end, the rule choice is wrong and `std::logic_error` is thrown (F-4) | **FAIL** (F-4); see also F-3 for invalid links |
| D-24 | double slot or origin step that rounds back: domain_error | thrown only outside the exact regime (C2 2^56, C4); never in the exact modes, including subnormal and x 2^1000 | PASS |
| D-25 | 128-bit root distances and charges for integers | exact gaps with root distances up to about 2^63 (C4 uint64) | PASS |
| D-26 | gap beyond W saturates to infinity, the offer is not posted | no exception from `order_children` in any run | PASS |
| D-27 | double: +inf beyond the range; "no wrong result is returned silently"; g may fall about one ulp below the charge | the claim "no wrong result is returned silently" is false: F-1 returns wrong values silently, with finite root distances well inside the range | **FAIL** (F-1) |

## Points that need agreement

- A-1 (F-1 scope). Q-55 option A says correctness is established for integer and dyadic weights.
  F-1 has integer weights. The lead needs to decide whether the scope is narrowed further (for
  example, all root distances of F below 2^53) or whether D-15 is replaced, and in either case
  how double inputs outside the scope are detected. A detection would also need to reject or
  flag the producers' output (Alg 4 base case) on the same graphs.
- A-2 (expiry D0 + S(C) in double). S is a rounded sum (132 rounding-size differences in C4).
  Sec 4.3 l.790-791 needs S >= the physical diameter. The heights included in S give a margin of
  at least h(C), so a rounded S could fall below the diameter only when the rounding error of S
  exceeds the heights in it. No failure was observed (every C4 failure disappears with D-14 off),
  and this was not proven either.
- A-3 (D-19 counter). The double origin fallback is recorded only in the diagnostics counter
  `eval_origin_rounding_fallbacks` (compiled out in release). It is not a direct search, so the
  fallback policy does not require a `FallbackRecord`. The lead should confirm that a counter is
  enough.
- A-4 (exception types on invalid lengths). The uint64 sentinel 2^64 - 1 as an arc length gives
  `std::invalid_argument` below an LCA and `std::overflow_error` between roots. +inf (double)
  between roots gives `std::overflow_error`. Def 4.1 excludes these lengths, and the documented
  type is `std::invalid_argument`.
- A-5 (D-12). Emissions keep the sign of -0.0 offers. They are numerically equal to 0, so no
  distance is wrong, but D-12 is still a candidate.
- A-6 (regression tests). As instructed, no test was added and no `tdd.py red` event was
  recorded. F-1 (4-vertex path) and F-4 (3 states) are minimal and ready to become
  `tests/regression/` cases. The lead decides whether to add them.
- A-7 (F-1 and the mapping rows). The Sec 4.3 Lemma 4.6 row of `PAPER_MAPPING.md` says that
  integer-valued doubles are exact. The D-15 / D-27 rows say "a few ulp" and "no wrong result is
  returned silently". These rows need correcting once A-1 is decided. The D-19 row says that
  the code comment "still calls this proposed D-16". The code now says D-19, so that sentence is
  stale.

## Not tested

- Non-dyadic doubles (for example, 0.1-based weights) with a rounding tolerance, and so the
  D-19 double fallback path. These are covered only by the existing
  `test_evaluate_vs_dijkstra.cpp`.
- `PortionedEvaluator` (Sec 7.4, Lemma 7.7), the Lemma 4.7 work bounds and the I1-I8
  invariants (other reviewers).
- Hierarchies with more than 60 states, and the Sec 5.2-5.4 producers above the base case
  (F-1 was confirmed only through the Alg 4 base case; Mode B `local_hierarchy` needs N >= 32
  and was not run on the counterexample).
- Sanitizer build (`make sanitizer`) was not run (the local machine is busy).
