# Quantitative Analysis of Thresholds and Constants (arXiv 2611.34911v1)

English | [日本語](thresholds.ja.md)

`L_N = ⌈log2(N+2)⌉` (Sec 2, l.284). In integers, `L_N = bit_length(N+1)` (verified to agree with `ceil(log2(N+2))` for N = 1..4999).
All computations were done exactly with Python arbitrary-precision integers (binary search per band; the script lives in the scratchpad and is not kept in the repository).
For the constructions of Sec 6–8 themselves, see `03_pivots_aux_outer.md`.

---

## 1. List of Threshold Conditions and Constants in the Paper

The "Kind" column: **Corr** = changing it can break distance correctness / **Time** = affects only complexity (and termination) / **Struct** = structural definition (derived from other constants).

| # | Condition / constant | Source | Use | Kind |
|---|---|---|---|---|
| 1 | Outer `N > L_N^64` | Sec 2 l.286, Alg 8 l.4 | If false, scan-based Dijkstra base case | Time |
| 2 | Local `k > 2L_k^4` | Sec 2 l.286, Alg 4 l.2 | If false, Kruskal hierarchy base case | Time (but also serves as a premise of #14) |
| 3 | Outer separation size `b = L_N^64`, `h = b+1` | Sec 7.1 l.1644, Lemma 8.4 | Owner decomposition. Atoms < b, `R ≤ hr` | Struct (h enters F, Lemma 7.2, and 8.2) |
| 4 | Local separation size `b = L_N^4`, `h = b+1` | Sec 5.1 l.966, Alg 4 l.6 | Local owner decomposition | Struct |
| 5 | Number of pivots `t = O(1 + N/L^64)`, `\|A0\| ≤ max(1, ⌊N/b⌋)` | Sec 5.1 l.970, Sec 6 l.1374 | Size of NameRoutes | Time |
| 6 | `u = ⌈log2(t+2)⌉`, `g = 4u+3` | Alg 5 l.2, Lemma 6.2 | Path length of critical edges | Time |
| 7 | `P` = the smallest power of 2 at least `8·256·g·u²` | Alg 5 l.3 | Qualification test `P·a* ≥ Z`, `d_e = Pλ_e` | Time (P enters F, so F must be derived from P) |
| 8 | Loop condition `\|U\| > 32t`, `16t` samples | Alg 5 l.5–6, Lemma 6.3 | Acceptance probability ≥ 1/2 | Time |
| 9 | Acceptance condition "at least half of U qualify" | Alg 5 l.13 | O(log t) rounds | Time |
| 10 | `T` = the smallest power of 2 at least `4t` | Sec 6.2 l.1473 | Band `[a, Ta)`, contraction `< a/T²`, window `[a/T², T²a]` | **Corr** (Lemma 6.4: requires T ≥ t) |
| 11 | Scales `a·2^i`, `2^{ℓ_a} ≥ tT²` | l.1490 | The top scale covers the whole component | **Corr** (Z ≥ d_S and existence of the LCA) |
| 12 | `H` = the smallest power of 2 at least `(16t)^5`, `j ∈ [0, H)`, `R = Δ(1+j/H)` | l.1491 | Bound on E[Z] | Time (R < 2Δ holds for any H) |
| 13 | Table value `4Δ + a`, `E[Z\|S] ≤ 256u²d_S` | Alg 6 l.9, Lemma 6.7 | Upper bound Z | **Corr** (4Δ+a is the deterministic upper bound) |
| 14 | Local `c = 10`, `N ≥ 32`, `h ≤ N`, `K < 16h`, `3hK² ≤ 768N³ ≤ N⁵`, `N·3hK² ≤ N^{c−2}` | Sec 5.4 l.1206 | Aligned test `s_U ≥ σ_U/N^c`, support `[·, N^16σ]` | **Corr** (if the inequalities fail, the support does not cover the weights) |
| 15 | Local macro test `w ≤ N³σ` | Sec 5.2 l.1021 | Lemma 5.5 | **Corr** |
| 16 | Local `K` = the smallest power of 2 at least `8h`, near means `≤ K²` | Sec 5.2 l.1029 | Shells | **Corr** (Lemma 5.5 requires K ≥ 8h) |
| 17 | Outer `K = 16`, near means `≤ K² = 256` | Sec 8.1 l.1977, Lemma 8.1 | Two-shell window | **Corr** (K ≥ 16 is required) |
| 18 | Radial test `2d_e ≤ w`, deletion test `w > J_x+J_y+d_e` | Sec 7.1 l.1651, Sec 8.1 l.1967 | Classification and deletion | **Corr** (assumes d_e is a true upper bound) |
| 19 | `\|J_x−J_y\| ≤ 3w/2`, `5M/8`, `M/16` | Lemma 7.2, 8.1 | Swing of 16 | Value inside a proof |
| 20 | Grid: from 7 octaves below the origin to 2 octaves above the maximum R, `32q ≤ r < 64q`, `U+ > R` | Sec 7.3 l.1697–1699 | Widths and calendars | **Corr** |
| 21 | Piece overlap `[r/4, 4R]`, piece crossing `r_v > 16R_u`, arc `U+_u → 16q_v` | Sec 7.3, Lemma 7.6, 7.8 | Time order | **Corr** |
| 22 | Atom edges `w' < 2max(q_u, q_v)`, `L_Y ≤ 2q\|Y\|` | Sec 7.3, Lemma 7.5 | Atoms | **Corr** (premise unit ≤ q, Lemma 7.7) |
| 23 | Tree piece root `R ≤ 17r/16`, `R_v ≤ 17R_u` | Lemma 7.3 | O(N) cells | Value inside a proof |
| 24 | `θ ≤ 4K²h·B0` (outer), `θ ≤ 3hK²B0` (local) | Lemma 8.2, 5.8 | Width of scales | Value inside a proof |
| 25 | `F` = the smallest power of 2 strictly greater than `2P(2h+1)(8K²h)²`, aligned means `F²s_U ≥ a_U`, support `[a, F³a]` | Sec 8.2 l.2080–2087 | Global grid | **Corr** (if F is too small, the support does not cover the weights) |
| 26 | `t_U ∈ (s_U/2, s_U]`, intervals are merged across gaps of ratio < 2 | Sec 8.2 | Atom cutoff | **Corr** |
| 27 | Normalization connectivity `4S(C) ≥ h(P)`, `ν ∈ (h/2, h]`, `c0 ≤ 2` | Sec 4.1, Lemma 5.10 | Unit alignment | **Corr** (within the scope of the Sec 4 owner) |
| 28 | MST stages `B1 = max(2, ⌈log2 L_k⌉)`, `B2 = L_k`, `B3 = ∞` | Sec 3.2 l.397 | Fibonacci heap caps | Time |
| 29 | DSU bands `b3 = 1+L_k`, `b2 = ⌈log2(b3+1)⌉`, `b1 = ⌈log2(b2+1)⌉` | Lemma 3.6 | Analysis only (the implementation is plain rank + path compression) | Time |
| 30 | LCA block `B = max(2, ⌊log2(k+1)/8⌋)` | Sec 3.4 l.505 | Table size O(2^B·B³) | Time |
| 31 | split–findmin block `B = max(2, ⌈log2(k+1)⌉)` | Sec 3.5 l.602 | Two-level blocking | Time |
| 32 | Fibonacci heap degree array `3 + 2⌈log2(s+1)⌉` | Lemma 3.3 | Array size | Time (if too small, out-of-bounds array access) |
| 33 | Random integers: draw `⌈log2 h⌉` bits, reject if at least h | Sec 2 l.307 | Uniform integers | Time |

---

## 2. Smallest N Satisfying Each Condition (Exact Values)

The band of N with `L = ℓ` is `[2^{ℓ−1} − 1, 2^ℓ − 2]`. For each band we found the smallest solution of the condition `N > f(ℓ)`.

### 2.1 Outer `N > L_N^64`

- **Smallest N = 589^64 + 1** (L = 589, 589 bits, ≈ 1.938×10^177, log2 ≈ 588.936)
  ```
  1938111981431118598628903515753047043696668555260980322369489207124043350454271623211173187479582684856639938987287147281633027619968763786364928512965370630489835656691917635842
  ```
- **Important: the set of N satisfying the condition is not upward closed.** At the start of band L = 590, an interval of N that fail the condition appears once more.
  - For `N ∈ [2^589 − 1, 590^64]` the condition is false (because L rises to 590 and `590^64 ≈ 2^589.09`). The interval has length ≈ 1.34×10^176, which is 6.6 % of the band.
  - **Every larger N satisfies it: N ≥ 590^64 + 1** (≈ 2.160×10^177, log2 ≈ 589.093)
    ```
    2160373020882152640073530204217550203038369458649719437938848235440878340387193555035212269915534933515580375897610000000000000000000000000000000000000000000000000000000000000001
    ```
    (Since 590 = 59·10, 590^64 ends in 64 zeros.)
  - Every band L ≥ 591 satisfies the condition throughout (verified).
- The fixed point of the estimate "L ≈ 64·log2 L" is L ≈ 589–590, consistent with the result above.
- Implementation consequence: precomputing "N0 = the smallest value" and testing `N ≥ N0` is wrong on the interval `[2^589−1, 590^64]`. **Evaluate `N > L_N^64` on the spot, exactly as the paper does.** (For 64-bit N the test is always false either way, so there is no practical harm.)
- For reference: with the non-strict `N ≥ L^64`, the smallest value is 589^64 and the exceptional interval is `[2^589−1, 590^64 − 1]`.

### 2.2 Local `k > 2L_k^4`

- **Smallest k = 209953** (L = 18, `2·18^4 = 209952`). Upward closed (every k ≥ 209953 satisfies it).
- Boundary: k = 209952 is the base case, k = 209953 is the radial construction.

### 2.3 For Reference: Thresholds with Other Exponents (Mode B Candidates)

| Condition | Smallest N (start of the upward-closed range) | L |
|---|---|---|
| `N > L^4` | 83522 | 17 |
| `N > 2L^4` | 209953 | 18 |
| `N > L^8` | 14048223625217 (≈ 1.4×10^13) | 44 |
| `N > L^16` | ≈ 3.97×10^32 | 109 |
| `N > L^32` | ≈ 1.31×10^77 | 257 |
| `N > L^64` | ≈ 1.94×10^177 (with an exceptional interval; all N satisfy it from ≈ 2.16×10^177) | 589/590 |

---

## 3. What Runs on Realistic Inputs (N = 10^3–10^8)

| N | L_N | L^4 | 2L^4 | log2(L^64) | Alg 8 l.4 | τ(N) | log log(L^64+4) | log log(L^4+4) |
|---|---|---|---|---|---|---|---|---|
| 10^3 | 10 | 10000 | 20000 | 212.6 | **Base case** | 2.731 | 7.73 | 3.73 |
| 10^4 | 14 | 38416 | 76832 | 243.7 | **Base case** | 2.900 | 7.93 | 3.93 |
| 10^5 | 17 | 83521 | 167042 | 261.6 | **Base case** | 3.019 | 8.03 | 4.03 |
| 10^6 | 20 | 160000 | 320000 | 276.6 | **Base case** | 3.110 | 8.11 | 4.11 |
| 10^7 | 24 | 331776 | 663552 | 293.4 | **Base case** | 3.182 | 8.20 | 4.20 |
| 10^8 | 27 | 531441 | 1062882 | 304.3 | **Base case** | 3.242 | 8.25 | 4.25 |
| 2^63 | 64 | 16777216 | 33554432 | 384.0 | **Base case** | — | — | — |

- **For N < 2^64, L ≥ 2, so L^64 ≥ 2^64 > N. The outer construction never runs.** Only Alg 8 l.5–7 executes (lightest edge per vertex pair → **scan-based** Dijkstra O(N² + m) → Lemma 8.5).
- Hence what never runs in Mode A: Sec 5 (local hierarchy), Sec 6 (NameRoutes / Estimate), Sec 7 (Alg 7), Sec 8.1–8.2, the Evaluate of Sec 4, and the MST, LCA, and split–findmin of Sec 3. **No randomness is used at all.**
- If the local test (Alg 4 l.2) were applied to such a k, k ≤ 209952 would be the Kruskal base case and k ≥ 209953 the radial construction. However, Alg 4 is called only from the outer level, so in Mode A it is never called within the range of this table.
- Scan-based Dijkstra takes about 10^12 operations at N = 10^6, which is impractical (see the recommendations in §7).

---

## 4. Recursion Depth and Where the Base Case Is Entered

**The construction does not recurse.** Alg 8 does not call itself. The call depth is fixed, with only the following 3 levels.

```
Alg 8 (N)                                   N ≤ L_N^64 → base case (Dijkstra)
 ├─ Alg 7 l.6: Alg 4 on each atom Y (|Y| < b = L_N^64)
 └─ Sec 8.2 / Lemma 5.15: Alg 4 on each atom (< b)
      Alg 4 (k)                             k ≤ 2L_k^4 → Kruskal hierarchy base case
       └─ Sec 5.4: components of in-block edges < t_U (< L_k^4 vertices) → exact Kruskal merge tree (sorting)
```

- At the first N where the outer construction runs (≈ 589^64), an outer atom has at most `b − 1 ≈ 2^588.9` vertices. Locally, `L_k ≤ 589`, the local `b = L_k^4 ≤ 589^4 ≈ 1.2×10^11`, and the local `K = 2^40`. A local atom has fewer than 1.2×10^11 vertices and is handled by a sorting-based Kruskal tree.
- If an atom has size k ≤ 209952, Alg 4 is the base case. Outer atoms have no lower bound on size, so every small atom enters the base case.
- For inputs with N close to `589^64 + 1`, there can be a single pivot (`⌊N/b⌋ = 1`). Then t = 1, NameRoutes is empty, all edges are radial, and J is the full SSSP from the pivot (Alg 7 handles everything). The outer machinery pays off only when `N/b` is large.

---

## 5. Parameter Values When the Construction Runs

### 5.1 Values at the Smallest Outer N (≈ 589^64)

| Quantity | Value |
|---|---|
| b = L^64, h | ≈ 2^588.94 |
| t | 1 (can be a single pivot) |
| u, g, P | 2, 11, 2^17 |
| T, H, ℓ_a | 4, 2^20, 4 (not actually used since t = 1) |
| F | 2^1808 |
| F³ (multiplier at the top of the support) | **2^5424** (far beyond the IEEE double exponent limit 2^1024) |

### 5.2 Running the Outer Construction in Mode B at N = 10^6 and 10^8 (by Choice of b)

t uses `2⌊N/b⌋ − 1` as an upper bound.

| N | b | t | u | g | P | T | H | ℓ_a | F | F³ | t·log⁴t / N |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 10^6 | L = 20 | 99999 | 17 | 71 | 2^26 | 2^19 | 2^104 | 55 | 2^64 | 2^192 | 7611 |
| 10^6 | L² = 400 | 4999 | 13 | 55 | 2^25 | 2^15 | 2^82 | 43 | 2^75 | 2^225 | 114 |
| 10^6 | L³ = 8000 | 249 | 8 | 35 | 2^23 | 2^10 | 2^60 | 28 | 2^86 | 2^258 | 1.01 |
| 10^6 | L⁴ = 160000 | 11 | 4 | 19 | 2^20 | 2^6 | 2^38 | 16 | 2^96 | 2^288 | ≈0 |
| 10^6 | 256 | 7811 | 13 | 55 | 2^25 | 2^15 | 2^85 | 43 | 2^74 | 2^222 | 218 |
| 10^6 | 4096 | 487 | 9 | 39 | 2^23 | 2^11 | 2^65 | 31 | 2^84 | 2^252 | 3.1 |
| 10^8 | L = 27 | 7407405 | 23 | 95 | 2^27 | 2^25 | 2^135 | 73 | 2^66 | 2^198 | 20090 |
| 10^8 | L² = 729 | 274347 | 19 | 79 | 2^26 | 2^21 | 2^111 | 61 | 2^79 | 2^237 | 292 |
| 10^8 | L³ = 19683 | 10159 | 14 | 59 | 2^25 | 2^16 | 2^87 | 46 | 2^92 | 2^276 | 3.2 |
| 10^8 | L⁴ = 531441 | 375 | 9 | 39 | 2^23 | 2^11 | 2^63 | 31 | 2^105 | 2^315 | 0.02 |
| 10^8 | 4096 | 48827 | 16 | 67 | 2^26 | 2^18 | 2^98 | 52 | 2^87 | 2^261 | 28.7 |

How to read it:
- `t·log⁴t ≤ N` (the basis for linearity in Lemma 6.8) first holds at b ≳ L³–L⁴. When b is small, the sorting in NameRoutes and the estimator `O(t log³ t)` dominate.
- With small b, H becomes 2^60–2^135 and does not fit in a 64-bit integer. Since H affects only time and `H ≥ 4tT²` suffices (03 §5 item 7), lowering it to `H = pow2_ge(4tT²)` in Mode B is reasonable.
- `F³` is 2^190–2^315. Computing `a·F³` in double overflows around `a > 2^700`. Held as an integer, it does not even fit in 128 bits.

### 5.3 Values at the First Local Radial Size (k = 209953)

L = 18, b = 104976, h = 104977, K = 2^20. `3hK² ≤ 768N³ ≤ N⁵`, `N·3hK² ≤ N^8`, `h ≤ N`, and `K < 16h` all hold (verified).

---

## 6. What Breaks When Thresholds Are Lowered (Mode B): Classification by Lemma

Premise: the structural constants (h = b+1, P, F, local K) are re-derived each time from the modified b. If not, the "Correctness" column changes.

| Lemma / construction | Parameter changed | Correctness | Complexity / termination | Notes |
|---|---|---|---|---|
| Alg 8 l.4 outer threshold | Relax `N > L^64` | **Preserved** | Breaks (linearity) | If the threshold is lowered but b stays L^64, it degenerates to a single pivot and J = full SSSP. b must be lowered together |
| Def 5.1 outer b | `b = L^64` → L^c etc. | Preserved (if h = b+1 is used consistently) | `t log³ t ≤ N` fails. `log log b` gets smaller | Atoms < b holds by construction |
| Lemma 6.2, 6.3 | 16t samples, 32t loop | Preserved | The guarantee of acceptance probability 1/2 can fail | The acceptance test is deterministic, so no false certificate is produced |
| Lemma 6.7 (Z ≥ d_S) | T, ℓ_a, window | **Breaks** | — | If T < t, the lifting error exceeds a and Z < d_S can occur → d_e is no longer an upper bound → the deletion test (#18) can delete shortest-path edges |
| Lemma 6.7 (E[Z]) | H, 256u² | Preserved | Acceptance rate only | Any H works as long as R < 2Δ is kept |
| Lemma 6.8 | Smaller P | **Preserved** (by definition, any value of P is valid as the certificate multiplier; F is derived from P) | **May not terminate** (acceptance probability can become 0) | With a cutoff (§7.2), it can terminate correctly by giving the remaining edges their own weights |
| Lemma 7.1, 7.2 | h | Preserved (if h is an actual hop bound) | — | Using a value smaller than the actual h makes F and Lemma 8.2, which use `R ≤ hr`, wrong |
| Lemma 7.3–7.6, 7.8 | 32/64, [r/4, 4R], 16R, 7 octaves | **Breaks** | — | The time order (marker arcs) breaks, and external offers arrive after their cells |
| Lemma 7.5, 7.7 | Atom threshold `2max(q)` | **Breaks** | — | The premise unit ≤ q (Lemma 7.7) fails |
| Lemma 8.1 | K = 16 | **Breaks** (when K < 16) | If K > 16, time only | A swing of 16 must be realized |
| Lemma 8.2, 8.4, Sec 8.2 | F | **Breaks** (F smaller than the formula) | — | The support does not cover the meaningful weights, the order of hierarchy heights (Def 4.1) breaks → Lemma 4.6 does not hold |
| Sec 5.2, Lemma 5.5 | Local K ≥ 8h, macro `N³σ` | **Breaks** | — | The forbidden-detour argument fails |
| Sec 5.4 (lowering the local threshold) | Relax `k > 2L^4` | **Conditional**: if `k ≥ 28` and `h = b_local + 1 ≤ k`, the paper's inequalities hold as is. Otherwise c = 10 and the N^16 support fail, and it **Breaks** | — | Remedy: use the actual value of h, and compute the support multiplier as the actual product, such as `(3hK²)²N³` |
| Alg 4 l.2 **always** use the local base case | Threshold to ∞ | **Preserved** (the Kruskal hierarchy satisfies Def 4.1 for any k) | Becomes `O(k log k)` sorting | The simplest safe option for Mode B. Also suited to testing the Sec 4 evaluator in isolation |
| Sec 3 block lengths (B1, B2, LCA, split–findmin) | Anything | Preserved | Time only | The Fibonacci heap degree array (#32) breaks if too small |
| Sec 4.1 normalization 4S ≥ h | 4 | **Breaks** (Sec 4 owner) | — | The dichotomy of unit alignment fails |

### Is Any Part's Correctness Probabilistic (Monte Carlo)?

**No.** Reasons:
1. Randomness is used only for the sample S in Alg 5 and for the permutation and j in Alg 6.
2. Certificate correctness `d_S ≤ Z ≤ P·a*` follows from the **deterministic** upper bound of Lemma 6.7(1) (distance to a common center < 2Δ, and Lemma 6.4). It holds however the random draws come out.
3. The acceptance test (at least half qualify) is deterministic; randomness only affects "whether it is accepted". Proposals from rejected rounds are discarded.
4. Everything from Alg 7 onward is deterministic (Lemma 7.9 explicitly says "deterministic").

Therefore **every terminating run returns exact distances (Las Vegas)**. This property is preserved in Mode B as long as the correctness constants (the "**Corr**" rows of the table) are not touched. What can fail in Mode B is only termination (when P is made smaller) and the complexity guarantees.

---

## 7. Implementation Recommendations

### 7.1 Mode A (Faithful)

1. **Threshold tests** use 64-bit saturating arithmetic. `pow_sat(L, e)` stops once the product reaches `UINT64_MAX` (Sec 2 l.304: "Products used solely in integer size tests may be capped once they reach the tested limit"). For N < 2^64, `N > L^64` is always false, so always record `FallbackRecord{"outer: N <= L_N^64 (paper base case)"}`.
   - When comparing by logarithms, do `64·log2(L) < log2(N)` in integers: fall back to saturating multiplication only when the lower/upper estimates via `bit_length` cannot decide.
   - Comparing against a precomputed N0 is not allowed (because of the exceptional interval in §2.1).
2. **Base case**: the paper specifies scan-based Dijkstra (O(N² + m)). Since it is unusable at realistic N, we recommend using the binary-heap version from the same Lemma 3.2 (identical output) and recording the deviation `D-n: base-case Dijkstra uses binary heap — scan O(N²) infeasible — correctness unchanged, complexity O(m log N)`. Which one is the default is for the lead to decide.
3. **Weight representation**: the correctness tests (`2d_e ≤ w`, `w > J_x+J_y+d_e`, `H·dist ≤ (H+j)Δ`, `F²s ≥ a`) must be done **exactly**. Recommended approach:
   - Restrict input weights to integers (uint64), with sums in `__int128`.
   - Hold power-of-2 multipliers (P, T², H, F, F³) as "(value, exponent)" pairs, and implement `cmp_shift(x, k, y)` (comparison of x·2^k with y) as a bit-length comparison → an exact shift only when needed (saturating shift comparison).
   - Floating point (double) rounds, so do not use it for correctness tests. If it is used, record it as a deviation.
4. **Randomness**: explicitly pass an RNG that accepts a seed. Draw j for H ≥ 2^64 with multiple words (the rejection method of Sec 2).

### 7.2 Candidates to Tune in Mode B (Experimental)

| Parameter | Paper value | Candidate | Rationale |
|---|---|---|---|
| Outer threshold | `N > L^64` | `N > b_out·c` (e.g. N ≥ 4b_out) or "always use the outer construction" | Time only |
| Outer b | L^64 | **L³–L⁴** (t log³ t ≤ N), or a constant for tests (64–4096) | Time only (re-derive h = b+1 and F) |
| Local threshold | 2L^4 | If lowered, keep `k ≥ max(28, b_local+1)`. Or always use the base case | §6 |
| Local b | L^4 | Around L² | Keep h ≤ k |
| H | (16t)^5 | `pow2_ge(4tT²)` | This is all the lemma's derivation needs |
| P | 2048gu² | As is (or measure with ×1/2–1/4) | Making it smaller lowers the acceptance rate |
| NameRoutes cutoff | None | After r consecutive rejections (e.g. 64), finish by certifying the remaining edges with their own weights | Correctness preserved (same treatment as l.15). Only the dictionary gets larger |
| Must not change | — | T, ℓ_a, window, K=16, the formula for F, the constants of Sec 7.3, local K ≥ 8h, `N³σ` | The "**Corr**" entries of §6 |

### 7.3 For Testing

- To exercise paths through the outer construction, Mode B needs a parameter that "forces the outer construction and makes b small". Mode A code never enters the outer branch at realistic N, so it is best to structure both as the same code switched by parameters, and to check the invariants (03 §4) with Mode B parameters.
- Threshold boundary tests: band boundaries of `L(N)` (N = 2^ℓ − 2, 2^ℓ − 1). Local k = 209952 / 209953. Check the outer test against an arbitrary-precision reference (Python-generated values inside the tests): 589^64, 589^64 + 1, 2^589 − 2, 2^589 − 1, 590^64, 590^64 + 1.
