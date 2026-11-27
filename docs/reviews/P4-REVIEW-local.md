# P4-REVIEW: independent correctness review of `src/paper/local/` (Sec 5)

- Paper: arXiv 2611.34911v1, Sec 5.1–5.5 (cited by section, lemma and text line of
  `.cache/paper/2611.34911v1.txt`; no passage is quoted).
- Code reviewed: commit `552b3de` (`main` at the time of the review), files
  `src/paper/local/{owner_decomposition,radial_states,threshold_joins,assemble,assembly_core,local_hierarchy,base_case,atom_replacement}.{hpp,cpp}`
  and `saturating.hpp`; the Mode B side `src/experimental/paper_experimental/local_*` only for
  the D-11 leak check. `src/paper/outer/` is out of scope.
- Reviewer: correctness reviewer agent. The expectations come from the paper,
  `docs/research/`, the Sec 5 rows and the D table of `docs/PAPER_MAPPING.md`, and the overflow
  rule in `include/sssp/paper.hpp`. They do not come from the implementer's summaries or from
  outputs of the implementation.

## 1. Verdict

**FAIL**, for one reason. With weights near 2^62, the radial Alg 4 output fed to Evaluate throws
`std::overflow_error`, even though every true d[u] + w is below the sentinel. Section 4.1 has the
counterexample. Whether this breaks the project's overflow rule depends on how the rule is read,
and which component should change is open; both need agreement (section 5).

Every Sec 5 construction step checked against the paper passed. This covers the lemmas'
invariants, path realization and the threshold-SCC characterisation, on 16,800 fuzz instances,
4 instances at the paper's sizes, and the boundary inputs. No overflow_error, FallbackRecord gap
or D-11 leak was found in `src/paper/local/` itself.

I also found one problem outside the scope: Evaluate (Sec 4) returns wrong double values, below
the true distance. It is reported in section 4.2 for the Sec 4 owner.

## 2. Verdict per item

| Item (paper) | Verdict | Evidence (section 6 lists the commands) |
|---|---|---|
| Sec 5.1 l.961–966: MST, reduced graph (drop nontree edges with w > D_T) | PASS | The tree spans the graph, and the cut property holds for every nontree edge. `reduced[e]` equals the exact 128-bit test w ≤ D_T(u,v) on all integer instances, including weights up to 2^64 − 2. Loops are omitted. |
| Sec 5.1 l.966–974: A0 marking, \|A0\| ≤ max(1, ⌊N/b⌋), A = A0 ∪ branch vertices | PASS | A brute-force Steiner tree (leaf pruning) gives A exactly. Components of T − A have fewer than b vertices and at most 2 pivot neighbours. b ∈ {1, 2, 3, 5, 8, n/3, n − 1}, plus the paper's b = L^4 at k = 209953 and 262143. |
| Def 5.1 l.975–992: corridor cuts, one pivot per cell, r_v, R_v ≤ h·r_v, hops ≤ h, σ = b_T on T_A | PASS | Brute-force tree paths. Every cut edge is a maximum of its pivot-to-pivot path. σ(p,q) is checked for all pivot pairs. |
| Lemma 5.2 l.993 | PASS | Intercell: w ≥ max(r_x, r_y, σ). Intracell with r_x ≠ r_y: both radii ≤ w. Every reduced edge: w ≤ h(r_x + r_y) + Nσ, in exact 128-bit arithmetic. |
| Lemma 5.3 l.1011 | PASS | For every vertex and every pivot q, the bottleneck of T[v, q] is at least r_v. |
| Sec 5.2 l.1021–1029: macro iff intercell and w ≤ N³σ; K = least power of two ≥ 8h | PASS | Compared with the exact 128-bit product on integers. With full-mantissa doubles, the rounded `sat::mul_int_pow` gave the same macro decision as an exact mantissa comparison on all 5,198 intercell edges tested. |
| Sec 5.2 l.1030–1037: near components, anchors, native shells a·K^i ≤ r < a·K^(i+1) | PASS | Union-find over near edges (ratio ≤ K², both radii positive). The anchor is the minimum radius. Shells are checked exactly. |
| Def 5.4 l.1038–1058: copies, windows, arc rules | PASS | Every arc projects to its own physical edge with the same length; every switch is a zero-length inward → native outward arc on one vertex. Blocks are exactly the components of the within-window arcs. Sizes: states ≤ 4N and arcs ≤ 8m + 2N. |
| Lemma 5.7 l.1084: path realization, DAG after removing macros | PASS | The block quotient without macro arcs is acyclic (Kahn). An independent 128-bit state-graph Dijkstra from weighted offers at native inward entries equals the physical Dijkstra at every vertex, minimised over the copies, and no copy is below it. This holds for uint64 (including weights near 2^63 and 2^64 − 2), integer doubles, subnormal doubles, and the paper's constants at k = 209953 and 262143. |
| Sec 5.3 l.1110–1124: I0, O0, I*, O*, θ | PASS | Recomputed by label-correcting min and bottleneck propagation over the non-macro arcs, seeded at macro heads and tails. All six per-state values agree exactly with `compute_joins`. |
| Lemma 5.8 l.1131 | PASS | When B0 is finite: r ≤ θ and B0 ≤ θ ≤ 3hK²B0, checked exactly. B0 and θ are finite together. |
| Lemma 5.9 l.1152 | PASS | For up to about 120 thresholds W per instance (all distinct arc and join weights, sampled), and for both ≤ W and < W: the components of the joins, restricted to states, equal the SCCs of the arcs (independent Tarjan). |
| Sec 5.4 l.1183–1236: s_U, t_U, Y components, alignment, supports, grid | PASS | The class bounds satisfy g < upper ≤ 2g and are sorted. Every global-grid node height is a grid point. Every Y node lies in one owner cell, with one radius and fewer than b vertices (l.1189–1194). Every node's states are exactly one SCC of the arcs below its class bound (≤ for an exact merge), and each child is exactly one SCC of the arcs below the node's height (l.1271–1275). Coverage of the branches: aligned 1,360, unaligned 28, unattached 1,717, Y ≥ 2 states 1,464 (fuzz with 4 seeds). |
| Lemma 5.10 l.1265 / Def 4.1 on the Alg 4 output | PASS | My own Def 4.1 checker, which does not use `validate_ordered_hierarchy`, checks H1–H6. It caught each of 7 mutations in a self-test. Evaluate from native entries equals the physical Dijkstra exactly on every integer and integer-double instance. The exceptions are the overflow case of 4.1 and full-mantissa doubles (4.2, a Sec 4 issue). |
| Alg 4 l.2–5: base case (Sec 2 l.286–290) | PASS | One state per vertex, and both orientations of every nonloop edge (n = 2 with 3 parallel edges and 2 loops gives 6 arcs). The Kruskal hierarchy has height equal to the weight, and its links are MST edges of that weight (c0 = 1). n = 0 and n = 1 (with or without a loop) work. A disconnected input gives a forest with one root per component (the supplement in `base_case.hpp`). |
| Alg 4 l.2: threshold N > 2·L_N^4, b = L_N^4 | PASS | `local_construction_applies` equals an exact 128-bit evaluation for every k in [0, 2·10^6] and at k ∈ {2^40 − 2, 2^40 − 1, 2^63 − 1, 2^63, 2^64 − 2, 2^64 − 1}. In [1, 2·10^6] the predicate changes exactly once, at k = 209953. `local_separator_size` is L^4, and the l.1206 preconditions hold at 209953, 262142, 262143, 524286, 524287 and 10^6. |
| Fallback recording (CLAUDE.md policy) | PASS | The base case pushes exactly one `FallbackRecord{"local: k <= 2*L_k^4 (paper base case)", k, m}` (checked on every Mode A fuzz instance and at k = 209952). The radial branch pushes none (k = 209953, 262143, and every Mode B instance). `substitute_atoms` pushes one record per atom that takes the base case. Note: with `diag == nullptr` nothing can be recorded; this is by design (section 5, point 5). |
| Input contract | PASS | A zero weight (edge or loop) and the integer sentinel 2^64 − 1 throw `std::invalid_argument` in both branches; 2^64 − 2 is accepted and handled exactly. A disconnected input to the radial branch, N = 31 with Mode B parameters, and b = N all throw `std::invalid_argument`; b = N − 1 is accepted. `CsrGraph` already rejects double +inf and NaN. |
| Overflow promise: no `overflow_error` from Sec 5 itself | PASS | No call into `src/paper/local/` threw on 200 radial instances with weights in [2^61, 2^63) and 25% of the edges within 1000 of 2^64 − 2, nor on any fuzz instance. |
| Overflow promise: auxiliary quantities exact or saturated (D-13, D-28, D-29) | PASS | Root distances ≥ 2^64 − 1 occurred at 4,016 vertices, and the reduced-graph decisions were still exact (D-28). R_v saturated at 242 vertices, always to exactly `infinity()` (D-29). N³σ, K²r and a·K^j saturated without changing a classification (exact checks). Saturated N^16σ supports merged, and the SCC invariant of Sec 5.4 still held (D-13). |
| **Overflow promise on the Alg 4 → Evaluate path** | **FAIL (needs agreement)** | Section 4.1. |
| D-30 (double: keep an edge when the right side rd(u) + rd(v) is +inf) | PASS (as documented) | In 24 instances with root distances +inf (a path of 2^1022 edges), nontree edges heavier than their tree path were kept, as D-30 states, and no lighter edge was dropped. The local construction did not throw, the hierarchy satisfies Def 4.1, and state-graph and Evaluate distances are exact on the finite part. Evaluate threw only where a reachable d[u] + w overflows, which the rule allows: 23 overflow throws and 1 D-24 `domain_error`. |
| D-11 must not leak into `src/paper/` | PASS | Statically: nothing in `src/paper/` includes `src/experimental/`, and `kExactHeight` appears only in Mode B. Dynamically: Mode A always declares c0 = 2, every global node height is a grid point, and exact-height nodes occur only as Y trees inside one owner cell with fewer than b vertices (a paper step, l.1198–1200). |
| Def 5.11 / 5.12, Lemma 5.13 (cap) | PASS | Every node under an atom's capped root has height ≤ c, and the links of an added root are < 2c. The states below the capped root are exactly the copies of the atom's states. |
| Lemma 5.14 / 5.15 (substitution) | PASS | Two atoms (sizes 32–51 and 1–41) plus singletons, cutoff c = 100, internal weights in [1, 199], so capping removes nodes, and outside arcs in [100, 200]. Mode A (base case) and Mode B (radial) atoms, uint64 and double, 15 seeds. My Def 4.1 checker passes. The state-graph distance from the entry, minimised over the copies, equals the physical distance at every vertex. Evaluate after `redirect_offers` equals Dijkstra. |

## 3. Paper steps versus implementation (notes)

- Interpretations were taken as recorded, not as errors: Q-31 (root 0; when nothing is marked,
  the root is the single pivot; ties on corridor maxima go to the smallest id), Q-32 (pivot arcs
  from every inward copy), Q-33 (outward window 0 holds shell 0 only, so s_U = a), Q-36
  (t_U ∈ [s_U/2, s_U]), Q-37 (supports per corridor key) and Q-38 / D-13 (saturation). Each is
  consistent with the cited lines. I found no input on which another admissible reading would
  change a distance.
- The Sec 5.4 consistency checks that throw `std::logic_error` (θ < s_U, θ outside its support, a
  separate-grid class reaching the first attachment class, a light-arc cycle, children not
  connected) never fired on any input of this review.
- Coverage at the paper's sizes is small, as the mapping already says. At k = 209953 the random
  tree and the clustered graph each had one pivot, so there was no macro edge and the grid was
  empty. At k = 262143 the caterpillar had 2 pivots and 35 grid classes. The unaligned-block path
  ran only with Mode B separator sizes.

## 4. Counterexamples

### 4.1 Alg 4 (radial) + Evaluate throws `std::overflow_error` although every d[u] + w fits

Minimal input after delta debugging (edge, vertex, separator and weight reductions), uint64, 36
vertices. The edges are `u v w`:

```
0 1 1            1 2 29734300731892      2 3 1   3 4 1   4 5 1   5 6 1   6 7 1   7 8 1   8 9 1
9 10 2452652296178016055     10 11 1     11 12 2916424925676143528
12 13 1  13 14 1  14 15 1  15 16 1  16 17 1  17 18 1  18 19 1
19 20 7712470701050284560    20 21 1  21 22 1  22 23 1  23 24 1  24 25 1  25 26 1
18 27 1  26 28 1  15 29 1  9 30 1  24 31 1  14 32 1  13 33 1  18 34 1  17 35 1
```

- How the input was found: fuzz family "caterpillar", weights "logwide" (2^k + rand(2^k),
  k ≤ 62). Fuzz seed 2, generator seed 2000717, n = 62, b = 61. It was then reduced to n = 36,
  b = 35, source 12. Mode B parameters `LocalHierarchyParams{1, 1, 35}`; the Sec 5.4 rule is the
  paper's (c = 10, N^16).
- Source 12, offer 0 at its native entry. Exact distances (128-bit oracle; baseline
  `dijkstra_binary_heap` and `bellman_ford` agree and do not throw): the maximum d(v) is
  7,712,470,701,050,284,574 (about 2^62.7), and the maximum d[u] + w over all edges is
  15,424,941,402,100,569,127. Both are below the sentinel 18,446,744,073,709,551,615.
- Expected under the overflow rule of `include/sssp/paper.hpp`: no throw, and the least emission
  at every vertex equals d(v).
- Actual: `evaluate` throws `std::overflow_error` with the message "state 70 received only offers
  beyond the range of the length type and never emitted". State 70 is the `lower_inward` copy
  (inward window 0) of vertex 20. Vertices 10..28 lie in shell 1 and vertices 2..9 in shell 0, so
  this copy can only be reached by going in to shell 0 and coming back out. Its exact constrained
  distance is 18,450,625,144,758,603,735 > 2^64 − 1. The native copies of vertex 20 have the
  correct distance 7,712,470,701,050,284,567.
- **With the paper's constants (Mode A).** The same gadget with w(1,2) = 2^41 (so that
  K = 2^20 gives the same two-shell layout), plus a chain of weight-1 edges from vertex 0 up to
  k = 209953 vertices (edges 0–36, 36–37, …, 209951–209952). This takes the radial branch with no
  FallbackRecord and reproduces the throw at state 70, with the same distances.
- Suspected step: Sec 5.2 / Lemma 5.7 (l.1105–1107). The paper allows a particular phase copy to
  have a larger constrained distance than the physical distance. In reals that is harmless; in
  uint64 such a copy can exceed the range while every physical distance fits. Evaluate (Sec 4.3,
  D-23) then sees a state whose only offers are beyond the range and throws. The Sec 5 output is
  paper-faithful: the minimum over copies is correct. The throw comes from how the two
  components are composed. This is not reachable from `paper_sssp` today, because Alg 8 takes its
  base case for every representable N. It becomes reachable once Alg 4's radial output is fed to
  Evaluate (outer construction or atoms) with integer weights of about 2^62.
- No regression test was added: the launch instructions allow only this report under `docs/`.
  The edge list above is the input for a regression case `regression: P4-REVIEW local copy
  overflow`.

### 4.2 (Outside the scope, Sec 4) Evaluate on double returns a value below the true distance

This is found through Alg 4's base case, which is part of the scope, but the faulty values come
from Evaluate.

- Minimal input (Mode A base case, 5 vertices, source 3). The edges are `u v w` with w in hex
  floating point:
  `0 1 0x1.cbd40160a6eap+30` (1928659032.163002), `1 2 0x1.46b72a0c9e398p-2`,
  `2 3 0x1.3333afd471754p-3`, `3 4 0x1.4ce38b1935768p-3` (0.1625433795040252).
- Expected: d(4) = w(3,4) = 0.1625433795040252, a single edge, so the sum is exact.
- Actual: Evaluate's least emission at vertex 4 is 0.16254329681396484, which is below the true
  distance (relative error 5.1e-7). The bound k·eps·d* of the correctness rules allows about
  2e-16 here.
- Isolation: on the same instances, the state-graph Dijkstra of the Sec 5 output has relative
  error 0, and Evaluate shows the same error on the base case and on the radial output. The error
  is at the ulp scale of the 1.9e9 edge, which suggests a large slot or origin quantity being
  added to and then subtracted from a small value (D-24 neighbourhood). Over full-mantissa random
  doubles (8 seeds × 5 families × 4 separator sizes), 186 vertices exceeded n·eps, and the
  largest relative error was 1.8e-7.
- Integer-valued and quarter-integer doubles (the inputs of the existing differential tests) are
  exact, which is why the existing tests do not see this.

## 5. Points needing agreement

1. **Reading of the overflow rule for state copies (4.1).** The rule says that `overflow_error`
   MAY be thrown when a tentative path length d[u] + w of a reachable u reaches the sentinel.
   With u a physical vertex and d[u] its (tentative) Dijkstra value, no such u exists in 4.1, so
   the throw breaks the rule. With "tentative path length" also covering constrained distances
   of state copies, the throw is allowed. The lead should decide, and record the decision in the
   rule text.
2. **Where to fix 4.1, if it is a violation.** (a) Evaluate (Sec 4): a state that received only
   out-of-range offers does not need to emit; throw only if some physical target (or requested
   state) ends with no emitted copy. (b) Sec 5: no change, because the state graph is
   paper-faithful and Lemma 5.7 l.1107 explicitly allows such copies. Option (a) keeps Mode A
   faithful; it is a candidate extension of D-23.
3. **4.2 belongs to the Sec 4 review.** It should get its own regression case and owner. It is
   listed here only because Alg 4's base case is the path on which it was found.
4. **D-30.** It holds as written. Its consequence, that Lemma 5.2 is not claimed for a kept edge
   when the right side overflows, had no effect on distances on the inputs tested. Because it is
   still a candidate, it needs the lead's adoption.
5. **Fallback records with `diag == nullptr`.** `local_hierarchy(g)` without a `Diagnostics`
   takes the base case without any record. This follows the API (the caller chose not to
   collect diagnostics), but every Mode A caller that reports timings must pass a `Diagnostics`.
   Consider asserting this at the top-level entry points.

## 6. Checks run

The reviewer program is in the worktree under `.cache/review/` (git-ignored, not committed). Its
oracles are a 128-bit binary-heap Dijkstra (physical and state graph), iterative Tarjan, Kahn,
union-find and a brute-force tree-path walker. It shares no code with `src/paper/`. It was built
with `c++ -O1 -g -std=c++20 -fno-fast-math` against `build/debug/libsssp.a` (debug preset, `-j4`).

| Command (`.cache/review/review_local <mode> <arg>`) | Inputs | Result |
|---|---|---|
| `selftest` | 7 mutations (link too long, height decrease, arc shorter than its edge, light cycle, θ changed, radius changed, arcs dropped) | all 7 caught |
| `thresholds` | k ∈ [0, 2·10^6] plus 6 values near 2^40, 2^63, 2^64 | 2,000,020 checks, 0 failures |
| `boundaries` | n = 0, n = 1 (with and without a loop), n = 2 multi-edges and loops, zero weights, sentinel, 2^64 − 2, disconnected, N = 31, b = N, b = N − 1; uint64 and double | 5,479 checks, 0 failures |
| `fuzz 25` | 8 families (path, star, tree + extra, grid, complete, clustered, multi-edge with loops, caterpillar) × 6 weight distributions (ties 1–3, uniform up to 10^6, 2^0..2^45, [2^61, 2^62), log-wide up to 2^63, two-scale) × b ∈ {1, 2, 3, 5, 8, n/3, n − 1}, n ∈ [32, 72), seeds 1–25, uint64 and integer double; plus Mode A (base case) on each graph | 16,800 radial instances, 120,741,199 checks, 20 failures, all `eval.throw` = 4.1 (seed 2, caterpillar, log-wide) |
| `fuzz 4` (with coverage counters) | as above, seeds 1–4 | 18,748,622 checks, 2 failures (4.1) |
| `huge 10` | 5 families, weights in [2^61, 2^63) and 25% within 1000 of 2^64 − 2, b ∈ {1, 2, 3, 7} | 200 instances, 1,172,291 checks, 0 failures |
| `special` | D-30 (2^1022 path, 24 instances), subnormal multiples of 2^-1074 (54), mixed 1e±300 (12) | 111,496 checks, 0 failures |
| `atoms 15` | section 2, Lemma 5.14 / 5.15 row; uint64 and double | 152,558 checks; 2 `invalid_argument` throws caused by my parameter choice (Mode B {1, 1, b} on a 5-state atom, rejected by the l.1206 precondition as documented); 0 other failures |
| `large 1` | paper constants: k = 209952 (base case), 209953 (tree + extra edges, clustered), 262143 (caterpillar, log-wide) | 72,301,071 checks, 0 failures |
| `minimize`, `gadgetA 41`, `gadgetB` | 4.1 | reproduced with Mode B and Mode A constants; baselines do not throw |
| `dblreal 8`, `isolate 2`, `mindbl 2` | 4.2 | 0 macro-test mismatches; 186 Evaluate tolerance failures, isolated to Evaluate |

The project tests: `ctest --preset debug -j4 -R "local|owner|radial|join|atom|assembl"` passed
152 of 152 tests (unit 91, regression 5, plus the differential and property tests; 261 s).
`make verify` was not run, as instructed.

## 7. Not tested

- The time bounds of Lemmas 5.7, 5.10 and 5.15 (not measured).
- Sanitizer builds (`make sanitizer` was not run; the machine is busy).
- At the paper's sizes, only 4 instances, with 1–2 pivots. Unaligned blocks, many corridor keys
  and merged supports were exercised only with Mode B separator sizes.
- Full-mantissa doubles through Sec 5 were checked only against the documented tolerance and for
  the macro decision. The exact-rational comparison was not done beyond the mantissa check of
  N³σ.
- Atom partitions with more than one level above the atoms, and atoms whose cutoff is below the
  parent height (only cutoff = parent height = 100 was used).
- The Mode B entry point `paper_experimental::local_hierarchy` beyond the D-11 leak check (D-8
  and D-11 correctness are not part of this review).
