# 04. Outer Scales and the Outer Pipeline (Paper Sec 8.2–8.3): Implementation Spec

English | [日本語](04_outer_scales_pipeline.ja.md)

Target: arXiv 2611.34911v1, Sec 8.2 (Lemma 8.4, l.2076–2142) and Sec 8.3 (Lemma 8.5 l.2145, Algorithm 8 l.2167–2196, proof of Theorem 1 l.2158–2217). Line numbers refer to `.cache/paper/2611.34911v1.txt`.
Tasks served: **P3-OUTER-SCALES** (`src/paper/outer/scales`) and **P3-OUTER-PIPELINE** (`outer_construction` in `src/paper/outer/shortest_paths`).
This note continues `03_pivots_aux_outer.md` §3–§4 and turns it into a spec the implementers can follow directly.

Notation: `[Paper]` is stated by the paper. `[Implied]` follows from the paper's statements by a short argument given here. `[Interpretation]` is our reading of an ambiguous place (each has a proposed Q in §5). `[Supplement]` is an addition the paper does not state. `[Measured]` was computed or run for this note (big-integer script, or the components on `main` linked into a scratch program; nothing of it is committed).

Scope limits of this note:
- The output of P3-OUTER-STATES (Sec 8.1) is not on `main` yet. §2.0 lists what the scales stage needs from it as **assumption A-1**; reconcile the names when that task lands.
- Nothing here was run end to end: the outer construction does not exist yet. Only the stages up to the radial states were run (§3.5).

---

## 0. Key Points

1. `F = 2^f` with `f = p + 23 + bit_width((2h+1)·h²)`, where `P = 2^p`, `h = b + 1`, `K = 16`. Keep only the exponent `f`; `F²s` and `F³a` are formed by `2f` and `3f` saturating doublings (D-13). `h` must be the **exact** `b + 1`, not the saturated word (§1.3).
2. With `uint64` weights and the paper's `P` (so `f ≥ 45`), `F²·s` saturates for every `s ≥ 1`: **every attached block is aligned and all supports merge into one grid**. The unaligned branch can be reached in `uint64` only with `f ≤ 31` (Mode B, `P ≤ 8`, `b ≤ 3`), otherwise only with `double` and a dynamic range above `2f` octaves (§1.5).
3. The Sec 5.4 grid / merge / link machinery that Sec 8.2 reuses lives in an anonymous namespace of `src/paper/local/assemble.cpp`, and it builds exact Kruskal trees for small components, where Sec 8.2 needs **atom leaves of a partition hierarchy** (Def 5.12). It cannot be called as is (§4, M-6, M-7).
4. Mode A cannot reach Alg 8 l.8–17 for `N < 2^64` (Q-04). The same code is exercised with Mode B parameters. The minimal forcing set is: outer threshold `0·L^0` and a constant separator `b`; everything else at the paper's values (§3.4).
5. Smallest contracted `N` at which the outer branch is nontrivial [Measured, §3.5]: `N = 2` (two pivots, `b = 1`), `N = 4` (a nonsingleton attached block, `b = 2`), **`N = 6` (a nonsingleton atom, so every stage of Sec 8.2 does real work, `b = 3`)**.
6. `outer_construction(c, diag)` as declared today has no generator and no parameters; `PaperOptions` has no seed. Both are needed before Alg 8 l.9 can be called (§4, M-3; Q-70).

---

## 1. The Factor F (Sec 8.2 l.2080–2081)

### 1.1 Definition

[Paper] `F` is the least power of two strictly greater than `X = 2·P·(2h+1)·(8K²h)²`, chosen once per graph from integer parameters. `P` is the common factor of NameRoutes (Alg 5 l.3), `h = b + 1` the hop bound of the outer owner decomposition (Sec 7.1 l.1644), `K = 16` (Sec 8.1).

[Implied] Where the factors come from: `2P(2h+1)` bounds the ratio of two macro weights sharing a name (Lemma 7.1: `λ/(2h+1) ≤ w < 2Pλ`), and `4K²h` is the ratio `θ/B0` of Lemma 8.2. For a finite join with new name `a` (§2.1): `a ≤ B0 ≤ θ < 8P(2h+1)K²h·a = X·a/(8K²h)`, so `θ < F·a` with a spare factor `8K²h`.

### 1.2 Exact (value, exponent) form (Q-06)

`F` is the pair (value 1, exponent `f`). `F²` is (1, `2f`), `F³` is (1, `3f`). `F` is never materialised as a length or as a machine integer.

[Implied] The least `f` with `2^f > X` is `bit_width(X)`. With `K = 16` and `P = 2^p`:

```
X = 2 · 2^p · (2h+1) · (2^11 · h)²  =  2^(p+23) · (2h+1) · h²
f = p + 23 + bit_width((2h+1) · h²)
```

[Measured] The closed form equals `bit_width(X)` on every row of §1.4.

### 1.3 Computing f without overflow

`f` is integer word work on sizes (Sec 2 l.302–304), never on lengths.

```cpp
// Inputs: p_log2 (unsigned), h as an exact natural number (multiword).
// Y = (2h + 1) * h * h      -- three multiword operations; Y < 2^(3·bits(h) + 2)
// f = p_log2 + 23 + bit_width(Y)
```

- `h` must be exact. In Mode A `b = L_N^64` does not fit a word for any `N` (`L_N ≥ 2`), and `owner::outer_separator_size` / `OwnerDecomposition::hop_bound` hold the **saturated** value `2^64 − 1`. Deriving `f` from the saturated `h` is wrong: at `N = 10^6` it gives `f = 233` instead of 871 [Measured]. Form `L^64` and `h` with multiword arithmetic (64 word-by-multiword multiplications; `L ≤ 65`, so `h < 2^386`).
- `src/experimental/paper_experimental/params.cpp` already does this (`OuterDerived::f_log2`, class `Nat` in an anonymous namespace). `src/paper/` cannot include `src/experimental/`; see M-2.
- `f`, `2f`, `3f` fit `unsigned` for every representable input (`f ≤ 29 + 23 + 3·386 + 2`).

Applying `F` to lengths (D-13, Q-38; `src/paper/local/saturating.hpp`):

```cpp
// aligned(s, a)  :=  !(sat::mul_pow2(s, 2 * f) < a)        // F² s ≥ a, Sec 8.2 l.2086
// support_end(a) :=  sat::mul_pow2(a, 3 * f)               // F³ a,    Sec 8.2 l.2081
```

[Implied] Both are exact as decisions. `sat::mul_pow2` returns `WeightTraits<W>::infinity()` exactly when the true product does not fit, and then the true product is at least `infinity()`, which exceeds every finite length of the type. So `F²s ≥ a` has the same truth value as in exact arithmetic for every finite `a`, and a saturated support end overlaps every later support exactly when the true end does (a later support starts at a finite value). `sat::mul_pow2` stops at infinity, so the loop runs at most `min(2f, word bits)` times for integers and at most 2098 times for `double`.
Do not use `LenOps::mul_pow2` here: it throws on overflow instead of saturating.

### 1.4 Exponents [Measured, Python big integers]

`u = ⌈log2(t+2)⌉`, `g = 4u+3`, `P = 2^p` = least power of two ≥ `2048·g·u²` (Alg 5 l.2–3).

Paper-size and Mode A rows (`b = L_N^64`):

| N | b | t | p | f | 2f | 3f |
|---|---|---|---|---|---|---|
| 589^64 + 1 (least outer N) | 589^64 ≈ 2^588.94 | 1 (the only possible value) | 17 | **1808** | 3616 | 5424 |
| 589^64 + 1 (formula only) | 589^64 | 66 | 22 | 1813 | 3626 | 5439 |
| 589^64 + 1 (formula only) | 589^64 | 2^20 | 27 | 1818 | 3636 | 5454 |
| 590^64 + 1 | 590^64 ≈ 2^589.09 | 1 | 17 | 1809 | 3618 | 5427 |
| 10^3 (unreachable in Mode A) | 10^64 | 1 | 17 | 679 | 1358 | 2037 |
| 10^6 (unreachable in Mode A) | 20^64 | 1 | 17 | 871 | 1742 | 2613 |
| 2^63 (unreachable in Mode A) | 64^64 = 2^384 | 1 | 17 | 1194 | 2388 | 3582 |

Mode B rows (the `t` column is the upper bound `2⌊N/b⌋ − 1` used in `thresholds.md` §5.2):

| N | b | t | p | f | 2f | 3f |
|---|---|---|---|---|---|---|
| 10^6 | 20 | 99999 | 26 | 64 | 128 | 192 |
| 10^6 | 256 | 7811 | 25 | 74 | 148 | 222 |
| 10^6 | 400 | 4999 | 25 | 75 | 150 | 225 |
| 10^6 | 4096 | 487 | 23 | 84 | 168 | 252 |
| 10^6 | 8000 | 249 | 23 | 86 | 172 | 258 |
| 10^6 | 160000 | 11 | 20 | 96 | 192 | 288 |
| 10^8 | 27 | 7407405 | 27 | 66 | 132 | 198 |
| 10^8 | 729 | 274347 | 26 | 79 | 158 | 237 |
| 10^8 | 4096 | 48827 | 26 | 87 | 174 | 261 |
| 10^8 | 19683 | 10159 | 25 | 92 | 184 | 276 |
| 10^8 | 531441 | 375 | 23 | 105 | 210 | 315 |

Small separators (forced-outer tests). `f` with the paper's `P` at `t = 1..2` (`p = 17`), `t = 6` (`p = 19`), `t = 66` (`p = 22`), and with D-17 values `P = 1`, `P = 2`:

| b | h | bit_width((2h+1)h²) | f (t ≤ 2) | f (t = 6) | f (t = 66) | f (P = 1) | f (P = 2) |
|---|---|---|---|---|---|---|---|
| 1 | 2 | 5 | 45 | 47 | 50 | 28 | 29 |
| 2 | 3 | 6 | 46 | 48 | 51 | 29 | 30 |
| 3 | 4 | 8 | 48 | 50 | 53 | 31 | 32 |
| 4 | 5 | 9 | 49 | 51 | 54 | 32 | 33 |
| 8 | 9 | 11 | 51 | 53 | 56 | 34 | 35 |
| 16 | 17 | 14 | 54 | 56 | 59 | 37 | 38 |
| 33 | 34 | 17 | 57 | 59 | 62 | 40 | 41 |
| 64 | 65 | 20 | 60 | 62 | 65 | 43 | 44 |

These agree with the values already tested in `tests/unit/paper_experimental/test_experimental_params.cpp` (871; 77 for `b = 160000`, `t = 11`, `P = 2`) and with `thresholds.md` §5.1–5.2.

### 1.5 What saturation means for each weight type [Implied]

| Type | `F²·s ≥ a` | Supports `[l, F³a]` |
|---|---|---|
| `uint64`, `2f ≥ 64` (always with the paper's `P`) | true for every `s ≥ 1`: every block with a finite attachment is **aligned** | `3f ≥ 84 > 64` for every parameter set: every support end saturates, so **all supports merge into one grid** |
| `uint64`, `2f ≤ 63`, i.e. `f ≤ 31` | can be false: needs `p + bit_width((2h+1)h²) ≤ 8`, i.e. `b = 1, p ≤ 3`; `b = 2, p ≤ 2`; `b = 3, p = 0` | as above |
| `double`, `f ≤ 1022` | false when `a > 2^(2f)·s` (dyadic weights spanning more than `2f` octaves) | separate grid components when supports are more than `3f` octaves apart |
| `double`, paper size (`f ≥ 1808`) | `2^3616` exceeds the whole range of `double` (2^2098 between the least subnormal and the largest finite value): always aligned | always one grid |

Consequence for tests: the unaligned branch and multi-component grids need either Mode B with `p_log2 = 0` and `b ≤ 3` on `uint64`, or `double` with dyadic weights (§6).

---

## 2. Lemma 8.4, Stage by Stage (Sec 8.2 l.2076–2142)

### 2.0 Inputs

From existing components on `main`:

| Input | Type | Source |
|---|---|---|
| contracted graph | `CsrGraph<W>` | `outer::ZeroContraction<W>::graph` |
| owner decomposition (outer `b`) | `OwnerDecomposition<W>`: `radius` (r_v), `owner`, `pivots`, `corridors`, `hop_bound` | `owner_decomposition(g, b)` |
| classification | `owner::OuterEdgeClassification<W>`: `edge_class`, `name`, `dictionary`, `p_log2` | `owner::classify_interowner_edges` |
| auxiliary graph | `owner::AuxGraph<W>`: `aux_vertex`, `radius`, `owner_distance` | `owner::build_aux_graph` |
| calendars (for separate grids) | `owner::Calendars<W>`: `graph_piece`, `radius_index`, `graph_pieces[].origin`, `.grid` | `owner::OwnerDistancesResult<W>::calendars` |

**Assumption A-1** (P3-OUTER-STATES, not on `main`; field names provisional). The scales stage needs:

| Needed | Meaning |
|---|---|
| `StateGraph<W> graph`, `projection` | every state arc (ordinary and special), physical image per state |
| `block` per state, `num_blocks`, block kind | kind ∈ {positive ordinary, pivot copy, special}; special states and pivot copies are singleton blocks (Sec 8.1 l.1996) |
| `arc_kind`, `arc_edge` | to recognise within-window arcs and their input edge |
| `native_inward` per vertex | Alg 8 l.14 seed; for a pivot its inward pivot copy |
| per state `theta`, `b0`, `b0_edge` | Lemma 8.2 values and the attaining macro edge ("a finite B0 stores the pointer of an attaining macro", l.2078–2079) |
| join records | block MST edges at their weight, state-to-owner joins at θ, corridor joins at the corridor weights (Lemma 8.3), in the shape of `ThresholdJoins<W>::joins` |

If P3-OUTER-STATES reuses `build_radial_states`, it can pass `EdgeClassification::macro_key[e] = cls.name[e]`; `ThresholdJoins::b0_key` is then the dictionary index of the attaining macro's name, which is exactly what §2.1 needs.
Units: Alg 7 also returns `scaled_distance = 2^7·J` (D-21). Everything in this section is in **input units** (r, θ, MST weights, names).

### 2.1 Stage S0 — new names (l.2077–2079)

[Paper] Replace the name of each macro edge by the least physical weight of any macro edge carrying that name.

```cpp
// in : g, cls (edge_class, name, dictionary)
// out: std::vector<W> name_weight   // per dictionary index; infinity() if no macro edge has it
for (EdgeId e : edges with cls.edge_class[e] == kMacro)
  name_weight[cls.name[e]] = min(name_weight[cls.name[e]], g.edge(e).weight);
```

- [Interpretation, Q-65] "physical weight" is the input weight `w_e` (Lemma 7.1 is stated for `w`), not the projected length `A_e`, and the minimum runs over macro edges only.
- Existing functions: none; a linear scan. Lengths are only compared.
- Invariants (testable): for every macro `e` with `a = name_weight[name[e]]`: `a ≤ w_e`; `w_e < 2P(2h+1)·a` (saturating product); `a` is the weight of some macro edge with that name. For every state `z` with finite θ and `a_z = name_weight[name[b0_edge[z]]]`: `a_z ≤ b0[z] ≤ theta[z]` and `theta[z] < F·a_z`.

### 2.2 Stage S1 — F

`f = outer_scale_factor_log2(p_log2, exact h)` (§1.3), once per graph. [Supplement] When there is no macro edge (`t = 1`), `f` is computed but never applied.

### 2.3 Stage S2 — block scales and attachment (Sec 8.1 l.1999; Sec 8.2 l.2083–2086)

```cpp
// in : states (A-1), od.radius, name_weight
// out: std::vector<OuterBlockScale<W>> blocks      // per block
struct OuterBlockScale {
  Kind kind;            // kSingleton, kAligned, kUnaligned, kUnattached
  std::size_t size;
  W s;                  // s_U = min over states x of U of od.radius[projection[x]]   (Q-35)
  Vertex s_vertex;      // a vertex attaining s_U (least id) -- needed for the separate grid
  StateId z;            // a state minimising theta over U (least id); kNoState if all infinite
  std::uint32_t name;   // dictionary index of the name of b0_edge[z]; a_U = name_weight[name]
  W t;                  // t_U, set in S5
};
```

- [Paper] A singleton block (every pivot copy, every special state, and any one-state ordinary block) keeps only its own finite θ join and skips the block-scale test, the grid, the MST classification and the atom substitution (l.2083–2085).
- [Paper] For a nonsingleton positive ordinary block with a finite attachment, pick a minimising state and its name `a_U` (l.2085–2086).
- [Interpretation, Q-35] `s_U` uses the original radius `r_v = b_T(v, p(v))`, not `J`.
- [Supplement] Ties (several minimisers of θ, several vertices of radius `s_U`) go to the least id. Any choice is valid: every minimiser `z` satisfies `a_z ≤ θ(z) < F·a_z`.
- Invariants: `0 < s_U ≤ r_x` for every state `x` of `U`; `s_U ≤ θ(z)` for every `z ∈ U` (Lemma 8.2 `r ≤ θ`); `θ(z_U) = min_U θ`; within a block one copy per physical vertex (the states of a block project injectively).

### 2.4 Stage S3 — aligned or unaligned (l.2086–2087; Lemma 8.4 proof (1), (2))

```cpp
// in : blocks (kind not yet split), name_weight, f
// out: blocks[U].kind in {kAligned, kUnaligned} for attached nonsingleton blocks
// group the attached nonsingleton blocks by blocks[U].name; in a group with name weight a:
//   s_star = pivots::median_prune(keys = {s_U}, feasible(s) = !(sat::mul_pow2(s, 2f) < a))
//   U is aligned  iff  s_star exists and !(s_U < *s_star)
```

- [Paper] Aligned iff `F²·s_U ≥ a_U`. The paper groups the tests by name and uses median pruning with `O(log F)` doublings per tested median (l.2130–2131).
- Existing functions: `pivots::median_prune` (`src/paper/pivots/name_routes.hpp`; the predicate is monotone in `s`), `sat::mul_pow2`. Grouping by name is a counting sort on the dictionary index.
- [Supplement] A per-block test gives the same classification with up to `2f` doublings per block; it changes only the work bound. If used, record it as a time-only deviation.
- Invariants:
  - Unaligned block: every block MST edge (indeed every within-window arc of `U`) has weight `< a_U/2`, tested as `sat::add(w, w) < a_U`. [Paper, l.2086–2087] Reason: `w ≤ 4K²·J_v ≤ 4K²h·s_U < a_U/2`.
  - Unaligned block: θ is the same for all states of `U` (Lemma 8.4 proof (1), "equal attachment threshold").
  - Aligned block: `s_U ≤ θ(z_U) < F·a_U` and `a_U ≤ F²·s_U`.
  - Differential: the grouped result equals the per-block test `!(sat::mul_pow2(s_U, 2f) < a_U)`.
  - Comparison count: at most `⌊log2 n_a⌋ + 1` predicate tests per name group of size `n_a`.

### 2.5 Stage S4 — supports and the global grid (l.2081–2083, l.2087–2091; Lemma 8.4 proof (2))

Keys: one per dictionary index that names some finite join, and one per corridor of `T_A` (`od.corridors`).

```cpp
// in : name_weight, od.corridors, aligned blocks, theta joins, f
// out: std::vector<W> grid;                       // increasing grid points of all components
//      std::vector<std::uint32_t> component_of;   // per grid point
//      std::vector<GridSupport<W>> support;       // per key: lower, upper, lower_index, upper_index
// 1. lower[key] = a (name weight, or the corridor weight); upper[key] = sat::mul_pow2(a, 3f)
// 2. for every aligned U: lower[name_U] = min(lower[name_U], s_U)
// 3. sort the supports by lower; merge two neighbours when they overlap or when
//    next.lower < 2 * current_component.upper (gap of ratio below two)
// 4. per merged component: grid points lo, 2lo, 4lo, ... until the first point strictly above the
//    component's upper end (one final doubling as sentinel); stop at infinity()
// 5. per key: lower_index = largest grid point <= lower; upper_index = first grid point > upper
```

- [Paper] Supports `[a, F³a]`; lowering the start to `s_U` for aligned blocks; merging overlaps and gaps below a factor two; indices from the grid floor of `l` through the first point strictly above `u`; `O(k log F)` octaves for `k = O(t log²(N+2))` names.
- [Interpretation, Q-65] A corridor join uses the support of its own weight, not of its macro name.
- [Interpretation, Q-68] "gap smaller than a factor two" is the strict test `next.lower < 2·upper`. After this merge no class can cross a component start, so every class is exactly `[g, 2g)` (unlike Sec 5.4, where a class may stop at the next component, l.1279).
- [Supplement, as in `assemble.cpp`] Grid points above the largest finite join weight need not be generated: nothing is classified there.
- Existing functions: none callable. The same steps exist in `src/paper/local/assemble.cpp` (anonymous namespace) with `N^16σ` in place of `F³a`; see M-6.
- Invariants:
  - Every finite `θ(z)` lies in `[lower, upper]` of its own key; every corridor weight lies in its own support.
  - For an aligned `U` with key `κ`: `lower[κ] ≤ s_U`, and every block MST edge `w ≥ t_U` of `U` satisfies `grid[lower_index[κ]] ≤ w < grid[upper_index[κ]]` (Lemma 8.4 proof (1): `w ≤ 4K²h·s_U < F²·a_U`).
  - For an aligned `U` and every `z ∈ U` with finite θ: the supports of `name(z)` and of `name_U` lie in the same merged component. [Implied] `s_U ≤ θ(z) ≤ max(θ(z_U), 4K²h·s_U) < F²·a_U`, so `θ(z)` lies in both supports.
  - Grid: strictly increasing; inside a component `grid[i+1] = 2·grid[i]`; between components `grid[first of next] ≥ 2·(upper end of previous)`.
  - Size: at most `3f + 2 + log2(lowering)` points per input support before merging, where the lowering is at most `F²` (aligned `s_U ≥ a_U/F²`), so at most `5f + 2`.

### 2.6 Stage S5 — the cutoff t_U (l.2091–2093, l.2100–2104)

```cpp
// aligned U (key k):   t_U = largest grid point <= s_U, searched in [lower_index[k], upper_index[k]]
// every other nonsingleton positive ordinary block (unaligned, or without finite attachment):
//   v  = blocks[U].s_vertex;  x = aux.aux_vertex[v];  C = cal.graph_piece[x]
//   i  = cal.radius_index[x]              // cal grid is 2^7-scaled: grid[i] <= 2^7 r_v < 2 grid[i]
//   t_U = sat::mul_pow2(cal.graph_pieces[C].origin, i - 7)     // origin * 2^(i-7), input units
```

- [Paper] Both choices give `s_U/2 < t_U ≤ s_U` (Q-36: strict on the left, unlike Sec 5.4).
- [Paper] For a separate grid: pick `v` attaining `s_U`; every `J` in the block is at most `K²·J_v`; the largest local edge is at most `4K²·R_v` and the least significant edge exceeds `r_v/2`; use the auxiliary grid pointer of `v`, with each piece's table extended once by a fixed number of octaves.
- [Interpretation, Q-67] "auxiliary grid pointer" is the Sec 7.3 graph-piece grid of `v` (`Calendars`). That grid is stored scaled by `2^7` with `grid[7] = 2^7·origin`, so the grid point below `r_v` in input units is `origin·2^(i−7)` with `i ≥ 7` (because `r_v ≥ origin`). It is reached by doublings; no halving is needed.
- [Interpretation, Q-66] "Every other nonsingleton positive ordinary block" includes blocks with no finite attachment (as in Sec 5.4 l.1213).
- [Supplement] Separate-grid classes of `U` are `[t_U·2^j, t_U·2^(j+1))`, `j ≥ 0`, up to the class of the largest block MST edge (`≤ 4K²·R_v`, so `j ≤ 10 + ⌈log2 h⌉`). These are the piece-grid points. The existing piece grid is only guaranteed to reach `4·max R`; sharing one table as the paper describes needs it extended by 9 more doublings at the top (8 for `K²`, one sentinel). Generating the points per block by doubling from `t_U` gives the same classes and costs `O(log(K²h))` doublings per block (time only; M-8).
- Invariants: `t_U ≤ s_U < 2·t_U` (tested as `!(s_U < t_U) && s_U < sat::add(t_U, t_U)`); aligned: `t_U` is a grid point of the component of `name_U` and `grid[lower_index] ≤ t_U`; separate: `t_U = origin·2^j` for the piece of `s_vertex`.

### 2.7 Stage S6 — atoms (l.2000–2001, l.2133–2142; Def 5.11)

```cpp
// in : block MST join records, blocks[U].t
// out: per state: part id; std::vector<ReplaceableAtom<W>> atoms
// DisjointSets over the states; unite the endpoints of every block MST edge with weight < t_U.
// A component with one state is a singleton part. A component Y with >= 2 states is an atom:
//   states        = members of Y
//   edges (H)     = every within-window arc pair with both endpoints in Y, once per input edge,
//                   with its physical weight (heavier chords included)
//   spanning_tree = the block MST edges < t_U inside Y
//   cutoff        = t_U
```

- [Paper] The components of edges lighter than `t_U` are replaced by local representations; each atom lies in one component of `T − A` and has fewer than `b` vertices; inside an atom the owner and the original radius are constant (Lemma 8.4 proof (3)).
- [Interpretation, Q-53 extended; Q-74] The internal graph holds every within-window edge with both ends in the atom. `substitute_atoms` rejects any arc of `D` with both endpoints in one atom, so these must all go to `H`.
- [Implied] No other arc kind can have both endpoints in one atom: transfers and switches leave the window or the phase, pivot and special arcs touch singleton blocks, and a macro edge joins two owner cells while an atom has one owner. Check this and throw `std::logic_error` on a violation.
- [Interpretation, Q-74] One-state components are singleton parts (Def 5.12, Lemma 5.15 "every nonsingleton part"); Alg 4 is not called for them.
- Existing functions: `DisjointSets` (`src/paper/primitives/dsu.hpp`); types `ReplaceableAtom<W>`, `kNoAtom` (`src/paper/local/atom_replacement.hpp`).
- Invariants: `|Y| < b`; one owner, one `r`, and `r ≥ s_U ≥ t_U` on `Y`; `Y` lies inside one block and one component of `T − A`; `Y` projects injectively; every spanning-tree edge is `< t_U` (so `< 2·cutoff`); `H` is connected; every finite `θ(z)`, `z ∈ Y`, is `≥ t_U`.

### 2.8 Stage S7 — classifying the events (l.2091, l.2094; Lemma 8.4 proof (2))

| Event | Class |
|---|---|
| block MST edge `w < t_U` | none (inside an atom) |
| block MST edge `w ≥ t_U`, `U` aligned | global class `i`: `grid[i] ≤ w < grid[i+1]`, binary search in the index range of `name_U` |
| block MST edge `w ≥ t_U`, `U` unaligned or unattached | separate class `j` of `U`: `t_U·2^j ≤ w < t_U·2^(j+1)` |
| state-to-owner join at finite θ(z) | global class in the index range of `name(z)` |
| corridor join | global class in the index range of its corridor key |

- [Paper] `O(N)` events, each `O(1 + log log F) = O(τ(N))` comparisons.
- Invariants: the class found brackets the weight; binary-search steps per event `≤ ⌈log2(upper_index − lower_index + 1)⌉`; for an unaligned `U`, with `p_top` the lower end of its top separate class and `g` the global class of `θ(z_U)`: `2·p_top ≤ g` [Implied: `2·p_top ≤ 2·w_max < 8K²h·a_U/F² < a_U/2 < g`], so the whole block is one node before its first global event (l.2093–2094).

### 2.9 Stage S8 — the partition hierarchy (l.2094–2099; Def 5.12; Lemma 5.10 post-processing)

```cpp
// in : events with classes, parts (singletons and atoms), the state graph
// out: PartitionHierarchy<W> ph   // graph = outside arcs only; nodes; leaf_atom; atoms; c0 = 2
// 1. leaves: one per singleton part and one per atom
// 2. join vertices: states and one auxiliary vertex per pivot; DisjointSets (Lemma 3.6)
// 3. per separate-grid block, classes in increasing j; then global classes in increasing index.
//    All events of one class are merged simultaneously; the node's height is the class lower
//    boundary. Auxiliary nodes without states are dropped and unary nodes compressed.
// 4. links: every outside arc is assigned to the pair of children it joins at the LCA of its
//    endpoints' parts; at a node of class [g, 2g) the assigned arcs of weight < 2g form the child
//    adjacency; a DFS takes one arc per discovery edge (endpoints are the actual states).
// 5. order: at each node, topological order of the children by the assigned arcs of weight < g;
//    roots in condensation order of the arcs between roots.
```

- [Paper] Processing unions by integer grid classes, suppressing auxiliary nodes and ordering arcs lighter than each node gives a partition hierarchy; connecting trees come from the child-quotient searches of Sec 5.4; for a class `[a, 2a)` the strict form of Lemma 8.3 gives connectivity by arcs strictly below `2a`.
- [Interpretation, Q-68] Height = the class lower boundary, links = arcs below the upper boundary, hence `c0 = 2`, exactly as Lemma 5.10 does for Sec 5.4. The order "separate grids first, then global classes" is ours; any order consistent with the heights along each root path is equivalent.
- Existing functions: `DisjointSets`, `ForestQueries` (LCA and level ancestor), `add_internal_node`, `collect_roots`; type `PartitionHierarchy<W>`. The merge builder and the link/order passes of `assemble.cpp` are not callable (M-6, M-7).
- Invariants:
  - `substitute_atoms` accepts `ph` (parts partition the states, no arc of `D` inside an atom, every atom's cutoff `≤` its parent's height).
  - For every internal node `X` of height `g` and class `[g, g')`: the states below `X` are one SCC of the state arcs of length `< g'`, and the children of `X` are the SCCs of the arcs of length `< g` among them (Lemma 8.3 strict form; checked against Tarjan on small inputs).
  - Heights do not decrease towards the root; every link has length `< 2·height`.
  - Every atom's parent height `≥ t_U`.

### 2.10 Stage S9 — atom substitution (l.2099; Lemma 5.15)

```cpp
AtomSubstitution<W> sub = substitute_atoms(ph, &diag, params.local);   // one Alg 4 call per atom
```

- Existing function: `substitute_atoms` (`src/paper/local/atom_replacement.hpp`). It calls `local_hierarchy` on every atom, which pushes `FallbackRecord{"local: k <= 2*L_k^4 (paper base case)", k, m}` when Alg 4 takes its base case.
- Invariants: `validate_ordered_hierarchy(sub.hierarchy, contracted graph)` reports no violation with `c0 ≤ 2` (the acceptance criterion in `docs/TASKS.md`); `sub.old_state` composed with `projection` gives the physical image of every new state; one `local_hierarchy` call per atom (`atom_replace_atoms`).

### 2.11 Complexity claimed

[Paper, Lemma 8.4] Classification of the significant local events and the finite joins, including all global grids, in `O(m + N·τ(N))` work. Lemma 5.15 adds `O(m + N log log(b+4))` for the substitutions.

---

## 3. Algorithm 8 l.8–17 (Sec 8.3)

### 3.1 Stages, inputs and outputs

`c` is the `outer::ZeroContraction<W>` from Alg 8 l.2 (positive weights, no loops, connected, parallel edges kept), `g = c.graph`, `N = g.num_vertices()`.

| Alg 8 | Step | Call (on `main` unless marked) | Input → output |
|---|---|---|---|
| l.8 | MST (Lemma 3.5), owner cells, deletion of strictly tree-dominated nontree edges | `owner_decomposition(g, b, &diag)` | `g`, `b` → `OwnerDecomposition<W> od` (`tree_edges`, `reduced`, `pivots`, `owner`, `radius`, `owner_distance`, `corridors`) |
| l.9 | NameRoutes on the projected interowner edges | `owner::classify_interowner_edges(g, od, rng, name_routes_params, &diag)` | → `OuterEdgeClassification<W> cls` (`name`, `dictionary`, `p_log2`, `detour` d_e) |
| l.10 | macro / radial; `G'`, `T'` | (classification is part of the call above) `owner::build_aux_graph(g, od, cls, &diag)` | → `AuxGraph<W> aux` |
| l.11 | OwnerDistances (Alg 7) | `owner::owner_distances(aux, &diag, {.local = local_params})` | → `J` (`distance`, `scaled_distance`), `calendars` |
| l.12 | delete radial `e = xy` with `w_e > J_x + J_y + d_e` | **P3-OUTER-STATES** | → surviving edge set |
| l.13a | outer states, θ, joins (Sec 8.1) | **P3-OUTER-STATES** | → A-1 |
| l.13b | names, F, scales, grids, partition hierarchy (Sec 8.2) | **P3-OUTER-SCALES** (§2, S0–S8) | → `PartitionHierarchy<W> ph` |
| l.13c | capped local hierarchies | `substitute_atoms(ph, &diag, local_params)` | → `AtomSubstitution<W> sub` |
| l.14 | seed | `redirect_offers(sub, {{native_inward[c.source], 0}})` | → one `StateOffer<W>` at the native inward copy of `s`, or at its native local subcopy if `s` lies in an atom |
| l.15 | Evaluate, minimum per vertex | `evaluate(SeededOrderedHierarchy<W>{sub.hierarchy, {{c.source, 0}}, offers})` | → `EvalResult<W>::value` per new state; `dist[v] = min` over the new states `x` with `projection[sub.old_state[x]] = v` |
| l.16–17 | restore, tight arcs, first-discovery tree | `outer::restore_predecessor_tree(g_original, s, c, dist)` (already called by `paper_sssp`) | → `SsspResult<W>` |

Notes:
- l.8: the paper's "retain a lightest edge per vertex pair" belongs to the base case (l.5) only. The outer branch keeps parallel edges; NameRoutes deduplicates its own instance (Alg 5 l.4).
- l.9 is the only randomised step (Lemma 6.8). The generator is a `model::WordRng&` passed in explicitly.
- l.15 [Interpretation, Q-73]: the minimum runs over every state that projects to `v`, including the special states and all local subcopies. Lemma 4.6 gives no underestimate and at least one exact copy. A vertex without a finite emission means a violated invariant: throw `std::logic_error`, never return infinity for a vertex of the source component.
- `outer_construction` returns the contracted distances only; the tree is not produced by Evaluate.

### 3.2 Lemma 8.5 (l.2145–2157)

[Paper] After contracting the zero-edge components, exact distances of the contracted source component determine all original distances and a shortest-path tree in `O(m + n)`: orient `uv` from `u` to `v` when `d(u) + w(uv) = d(v) < ∞`, traverse from `s`, keep the first-discovery arcs, and output the sums of original weights down the tree.

This is implemented and VERIFIED on the base branch (`restore_predecessor_tree`, PAPER_MAPPING row "Sec 8.3, Lemma 8.5"). It takes distances only, so Alg 8 l.16 needs no new code. What the outer branch must respect:
- The input labels must be exact. `restore_predecessor_tree` throws `std::invalid_argument` when a finite label is not attained by a tight original arc. With `double`, the outer distances are sums taken in another order than the tree sums, so non-dyadic inputs can fail here (Q-06, Q-55 option A, Q-63). In Mode B that is a **failure**, never a reason to substitute another search.
- Unreachable original vertices are outside `c.vertex_map` and get infinity in the restoration; the outer construction never sees them.

### 3.3 Where fallbacks are recorded

| Place | Record | Who pushes it |
|---|---|---|
| Alg 8 l.4–5 (`N ≤ L_N^64`) | `kOuterBaseCaseScanReason` / `kOuterBaseCaseD1Reason`, `(N, m')` | `outer_base_case*` (exists) |
| Alg 7 l.6: Alg 4 on every Sec 7.3 atom | `"local: k <= <c>*L_k^<e> (paper base case)"`, `(k, m_Y)` | `local_hierarchy` via `cap_local_hierarchy` (exists) |
| Alg 8 l.13 / Lemma 5.15: Alg 4 on every Sec 8.2 atom | same reason, `(k, m_Y)` | `local_hierarchy` via `substitute_atoms` (exists) |
| Alg 5 retry cap (D-10, Mode B only) | NameRoutes' record | `name_routes` (exists) |
| `outer_construction` itself | **none** | it never replaces a stage by a direct search; an unsupported input throws |

- [Measured] Alg 7 pushes a record for one-vertex atoms too (`k = 1`, `m = 0`): four such records on the 6-vertex example of §3.5.
- A vertex can be covered by several records: once in Alg 7 and once per Sec 8.2 atom containing one of its copies (up to four copies: `in_i`, `in_{i−1}`, `out_i`, `out_{i+1}`). On the 6-vertex example the pair `{x, y}` forms an atom in three blocks. The sum of `FallbackRecord::vertices` can therefore exceed `N`; the fallback share needs a definition (Q-69).
- [Supplement] Proposed counters (not fallbacks): `outer_constructions`, `outer.pivots`, `outer.macro_edges`, `outer.radial_edges_deleted`, `outer.states`, `outer.blocks_{singleton,aligned,unaligned,unattached}`, `outer.atoms`, `outer.atom_states`, `outer.grid_points`, `outer.grid_components`, `outer.global_classes`, `outer.separate_classes`, `outer.scale_factor_log2`.
- The degenerate case `t = 1` (no macro, every block unattached, `J` is the full distance from the pivot) is the paper's construction, not a fallback; report it through `outer.pivots = 1` (Q-72).

### 3.4 Mode A and Mode B

Mode A: `takes_outer_construction(N)` is false for every `N < 2^64` (`L_N ≥ 2`, so `L_N^64 ≥ 2^64 > N`). Alg 8 l.8–17 is unreachable through `paper_sssp`; the least `N` that reaches it is `589^64 + 1` (Q-04, `thresholds.md` §2). The code is therefore exercised only by (i) Mode B parameters and (ii) calling `outer_construction` directly.

[Implied] Calling `outer_construction` directly with the paper's constants on a small graph is well defined: `b` saturates above `N`, the decomposition has one pivot, there is no macro edge, every nonsingleton block is unattached, and Alg 7 computes `J` as the full distance from the pivot. This is a useful test of the Mode A instantiation (§6, PL-2), but it does not exercise Sec 8.2 beyond the unattached path.

Minimal Mode B parameter set that forces the outer branch with real work (`paper_experimental::PaperParams`):

| Parameter | Value | Paper | Deviation |
|---|---|---|---|
| `outer.threshold` | `SizeThreshold{0, 0}` (applies for every `N ≥ 1`) | `1·L^64` | D-7 |
| `outer.separator` | `SeparatorRule::constant(b)`, `b ∈ {1, 2, 3, …}` | `L^64` | D-7 |
| everything else | paper values: sample `16t`, residual `32t`, paper `P`, paper `H`, no retry cap, local threshold `2·L^4`, local `b = L^4` | — | none |

Derived: `h = b + 1`, `F` from §1 (e.g. `f = 48` for `b = 3`, `t = 2`). With these values every outer atom has fewer than `b` vertices and takes Alg 4's base case (recorded), and NameRoutes never samples (`|U| ≤ 32t` when `t ≤ 65`), so every interowner edge certifies its own weight, `d_e = P·A_e > w_e`, and **every interowner edge is macro**; radial edges are exactly the intraowner edges.

Optional knobs that widen coverage (each a listed deviation):

| Goal | Knob | Deviation |
|---|---|---|
| NameRoutes sampling rounds, hence radial interowner edges and Alg 8 l.12 deletions | `sample_factor = 1`, `residual_factor = 2`, small `p_log2`, `retry_cap` | D-17, D-10 |
| unaligned blocks on `uint64` | `p_log2 = 0` with `b ≤ 3` (`f ≤ 31`) | D-17 |
| Alg 4's radial branch inside an outer atom | outer `b ≥ 33` and a local threshold and local `b` with `k ≥ max(32, b_local + 1)` for some atom (Q-03) | D-7 |

### 3.5 Smallest N at which the outer construction runs

"Runs" has several meanings. `N` is the contracted size; the outer threshold is disabled; NameRoutes and local parameters are the paper's.

| Level | Meaning | Least N | b | Example |
|---|---|---|---|---|
| L0 | the branch is taken | 1 | any | one vertex; one pivot; nothing else |
| L1 | two pivots, one corridor, one macro edge, a NameRoutes instance with an edge | **2** | 1 | one edge; both ends are pivots |
| L2 | L1 and a nonpivot vertex (positive ordinary states, nontrivial `J`) | **4** | 2 | path 0–1–2–3 |
| L3 | L2 and a nonsingleton block with a finite attachment (the alignment test runs) | **4** | 2 | path 0–1–2–3 (weights 10, 3, 3) plus the edge 1–3 (weight 5) |
| L4 | L3 and a nonsingleton atom (Alg 4, capping and substitution do real work): **every stage of Sec 8.2 is nontrivial** | **6** | 3 | edges 0–1 (4), 1–2 (4), 0–3 (100), 3–4 (8), 4–5 (1) |

Lower bounds [Implied]: a nonsingleton atom lies in one component of `T − A`, which has fewer than `b` vertices, so `b ≥ 3`. Two pivots need two marked vertices (with one mark there is no branch vertex), and each mark consumes at least `b` vertices of its own (Sec 5.1 marking rule), so `N ≥ 2b`: `N ≥ 6` for L4. With `b = 1` every vertex is marked, so a nonpivot vertex needs `b ≥ 2`, hence `N ≥ 4` for L2 and L3.

[Measured] The components on `main` (`owner_decomposition`, `classify_interowner_edges`, `build_aux_graph`, `owner_distances`, and `build_radial_states` with coordinate `J`, `K = 16`, after the l.12 deletion) were run on the examples and on 60,000 random connected multigraphs per `(N, b)`, `N ≤ 7`, `b ≤ 4`:
- L4 example: pivots {0, 3}; one macro (0–3); `r = 8` at vertices 4 and 5, `J = 8, 9`; blocks `{4, 5}` in the inward window 0 and the outward windows 0 and 1, each with `s_U = 8` and the edge of weight 1 below `s_U/2`; blocks `{1, 2}` with `s_U = 4` and an edge of weight 4 (significant, no atom).
- L3 example: pivots {0, 2}; macro 0–1; `r = J = 3` at vertices 1 and 3; three blocks `{1, 3}` with `s_U = 3` and the edge of weight 5.
- Random search: no instance with two pivots and a nonsingleton block for `N ≤ 3`; none with two pivots and a within-window edge below `s_U` for `N ≤ 5`; both occur at `N = 4, b = 2` and `N = 6, b = 3`.
The Sec 8.1 special states, θ, and everything after them were not run (not on `main`).

Beyond L4 (not measured; bounds only):
- NameRoutes sampling with the paper's constants needs more than `32t` distinct pivot pairs, hence `t ≥ 66` (`C(66, 2) = 2145 > 2112`). With `b = 1` this is `N = 66` (complete graph; no nonpivot vertex). Together with L4 (`b = 3`): `102 ≤ N ≤ 198` (`t ≤ 2|A0| − 1` and `|A0| ≤ N/b` give the lower bound; 66 marks of 3 vertices give the upper bound). The exact minimum was not computed.
- With D-17 (`sample 1t`, `residual 2t`) sampling needs `t ≥ 6`: `12 ≤ N ≤ 18` together with L4.
- Alg 4's radial branch inside an outer atom needs an atom of at least 32 vertices, so outer `b ≥ 33` and `N ≥ 66` for two pivots.

---

## 4. What the Existing Components Lack

| ID | Missing | Needed by | Detail |
|---|---|---|---|
| M-1 | Output type of Sec 8.1 | scales, pipeline | P3-OUTER-STATES is in progress; assumption A-1 (§2.0) lists the fields. In particular `b0_edge` (or `b0_key` = dictionary index) per state and the block kind |
| M-2 | Exact `f` in `src/paper/` | scales | A multiword natural (the `Nat` of `experimental/paper_experimental/params.cpp` is file-local) and `outer_scale_factor_log2(p_log2, exact h)`; an exact `b = L_N^64` (today `outer_separator_size` saturates). Mode B's `OuterDerived::f_log2` should then call the same function (Q-04: one implementation) |
| M-3 | Parameters and generator for the outer branch | pipeline | `paper::outer_construction(const ZeroContraction<W>&, Diagnostics&)` has neither. Needed: a parameter struct in `src/paper/outer/` (separator rule, `pivots::NameRoutesParams`, `LocalHierarchyParams`) and a `model::WordRng&`; a seed in `PaperOptions`; the Mode B driver (P5-EXP-THRESHOLDS) that maps `PaperParams` to it |
| M-4 | New names | scales | `name_weight` per dictionary index (§2.1) |
| M-5 | Outer classification as `EdgeClassification` | outer states | `OuterEdgeClass` → `EdgeClass` with the l.12 deletions omitted and `macro_key = name` (probably inside P3-OUTER-STATES) |
| M-6 | Reusable grid / support / simultaneous-merge / link / child-order engine | scales | The Sec 5.4 code is in the anonymous namespace of `local/assemble.cpp` (`MergeBuilder`, `floor_index`, `bucket`) and is tied to `OwnerDecomposition` corridors and `N^c`, `N^16`. Mode B already has a second copy (`local_assembly.cpp`, Q-59, P5-LOCAL-UNIFY). Extract a shared engine or accept a third copy explicitly (lead decision) |
| M-7 | Partition-hierarchy output | scales | `assemble_local_hierarchy` returns an `OrderedHierarchy` with exact Kruskal trees for small components; Sec 8.2 needs `PartitionHierarchy<W>` with atom leaves (`leaf_atom`, `ReplaceableAtom`) |
| M-8 | Separate grids in input units | scales | `Calendars::graph_pieces[].grid` is scaled by `2^7` and reaches `4·max R`; Sec 8.2 needs points up to `4K²·R_v` in input units (§2.6). Also the map vertex → aux vertex (`aux.aux_vertex`) for the lookup |
| M-9 | Checker | tests | `check_outer_scales(...)` returning one message per violated invariant of §2, in the style of `check_aux_metric` / `check_calendars` |
| M-10 | Local variants for outer atoms | Mode B | `substitute_atoms` and `owner_distances` accept `LocalHierarchyParams` only (size test and separator); the D-8 and D-11 engines cannot be selected for atoms. Irrelevant while atoms take the base case |
| M-11 | Diagnostics plumbing | pipeline | `absorb` is file-local in `sssp.cpp`; `evaluate` returns its own `EvalResult::diag`, which must be merged |
| M-12 | Minimum over copies | pipeline | new states → `sub.old_state` → `projection` → vertex; a linear pass, plus the "no emission" check of Q-73 |

Existing pieces that fit as they are: `owner_decomposition`, `classify_interowner_edges`, `build_aux_graph`, `owner_distances`, `sat::add` / `sat::mul_pow2`, `pivots::median_prune`, `DisjointSets`, `ForestQueries`, `PartitionHierarchy` / `ReplaceableAtom` / `substitute_atoms` / `redirect_offers`, `validate_ordered_hierarchy`, `evaluate`, `contract_zero_edges`, `restore_predecessor_tree`.

---

## 5. Ambiguities: Proposed Open Questions (Q-65 onward)

These are proposals for `open-questions.md` (the lead edits that file). None is decided.

| ID | Location | Question | Candidates | Proposed |
|---|---|---|---|---|
| Q-65 | Sec 8.2 l.2077–2083 | (a) Is the "least physical weight of any macro with that name" the input weight `w_e` or the projected length `A_e`, and over which edges? (b) Does a corridor join use the support of its own weight or of its macro name? | (a1) `w_e` over macro edges; (a2) `A_e`. (b1) own weight, one key per corridor; (b2) the name | (a1): Lemma 7.1 bounds `w`. (b1): the text lists "every pivot corridor weight a" next to the names |
| Q-66 | Sec 8.2 l.2085–2093 | Which blocks does "every other nonsingleton positive ordinary block" cover? | (1) unaligned blocks and blocks without a finite attachment; (2) unaligned only, unattached blocks unspecified | (1), as Sec 5.4 l.1213 says for the local case |
| Q-67 | Sec 8.2 l.2100–2104 | What is the "auxiliary grid pointer of v", in which units, and how many octaves are added? | (1) the Sec 7.3 graph-piece grid, converted to input units by doubling from the origin, extended by 9 doublings at the top; (2) the same grid used in `2^7`-scaled units (every block weight scaled, D-20 range limit); (3) a per-block grid from `s_U` with `t_U = s_U` as in Sec 5.4 (not what Sec 8.2 says; a deviation) | (1). The number of octaves is ours; the paper says "a fixed number" |
| Q-68 | Sec 8.2 l.2089–2099; Lemma 8.4 proof (2) | (a) Is the gap test strict (`next.lower < 2·upper`)? (b) Node height and link bound for a class `[a, 2a)`. (c) Order of processing separate and global classes | (a) strict / non-strict. (b) height `a`, links `< 2a`, `c0 = 2` (as Lemma 5.10) / height `2a` with links `< 2a`. (c) separate grids first / interleaved by value | (a) strict ("smaller than a factor two"). (b) height `a`. (c) separate first; equivalent to any order consistent with heights |
| Q-69 | CLAUDE.md fallback policy; Alg 7 l.6; Lemma 5.15 | How is the fallback share defined when a vertex is covered by several local base cases (Alg 7 atom, and up to four Sec 8.2 atoms through its copies)? | (1) raw sums per reason, stated as "can exceed N"; (2) distinct physical vertices (extra bookkeeping); (3) states in base-case atoms divided by all states, reported per stage | (3) per stage plus (1) raw; needs the lead's agreement |
| Q-70 | Alg 8 l.9; `include/sssp/paper.hpp` | Where do the generator and the parameters of the outer branch enter the Mode A API? | (1) `PaperOptions::seed` with a fixed default and an internal `WordRng`; (2) an overload taking `WordRng&`; (3) no public way, tests call `outer::outer_construction` directly | (1) and (3) together; a public header change, so the lead decides |
| Q-71 | Sec 8.2 l.2080–2081; Sec 7.1 l.1644 | In Mode A `b = L_N^64` is held saturated. From which `h` is `F` derived? | (1) exact multiword `h` (M-2); (2) skip `F` when the saturated `b ≥ N` (then `t = 1` and `F` is unused) | (1); (2) adds a branch the paper does not have |
| Q-72 | Alg 8 l.8–13 with `t = 1` | Is the single-pivot run (no macro, `J` = full distances by Alg 7) a fallback? | (1) no, it is the paper's construction; counter only; (2) record a FallbackRecord | (1) |
| Q-73 | Alg 8 l.15 | Over which states is "the minimum emitted value for each vertex" taken, and what if a vertex has none? | (1) all states projecting to the vertex (ordinary, special, local subcopies); none → `std::logic_error`; (2) ordinary copies only | (1): special states project to physical walks, so they cannot underestimate |
| Q-74 | Sec 8.1 l.2000–2001; Def 5.11–5.12 | (a) Which edges form the internal graph of a Sec 8.2 atom? (b) Are one-state components atoms? | (a1) every within-window edge with both ends in the atom (heavier chords included, as Q-53); (a2) only the edges `< t_U`. (b1) singleton parts, no Alg 4 call; (b2) one-state atoms (as Alg 7 does today) | (a1), required by `substitute_atoms`; (b1), following Def 5.12 |

Places checked and found consistent (no question): the alignment bound of Lemma 8.4 proof (1) (`a_U/F² ≤ s_U ≤ F·a_U`) follows from Lemma 8.2 and Lemma 7.1 with a spare factor `8K²h`; `s_U/2 < t_U ≤ s_U` holds for both grids; the atom cutoff `t_U` satisfies Def 5.11 (spanning tree edges `< t_U < 2t_U`, parent height `≥ t_U`).

---

## 6. Proposed Acceptance Tests

Weights: `uint64` unless stated; `double` only with dyadic values (Q-55 option A). Seeds fixed and varied where NameRoutes samples. Oracles are independent: Bellman–Ford on the original graph, Tarjan SCC on the state arcs, Python big integers for `f`.

### 6.1 P3-OUTER-SCALES

Proposed files: `tests/unit/paper/outer/test_outer_scales.cpp`, `tests/property/paper/test_outer_scales_property.cpp`.

| Test | Input | Checks |
|---|---|---|
| SC-1 `scale_factor_log2 matches big-integer values` | `(p, b)` from §1.4, including `b = 589^64` (1808), `20^64` (871), `b = 3, p = 17` (48), `b = 160000, p = 1` (77) | `f` equal to the table; `2^f > X ≥ 2^(f−1)` recomputed in the test with a second multiword routine |
| SC-2 `scale factor is not derived from a saturated h` | `N = 10^6`, paper `b` | 871, not 233 |
| SC-3 `F^2 s >= a is decided exactly` | random `s`, `a`, `f` in `uint64` and dyadic `double`, around the saturation boundary | equals the exact answer (`unsigned __int128` shifts, or exponent arithmetic for powers of two) |
| SC-4 `new names are the least macro weight` | hand classification with two names and four macros | `name_weight`; `a ≤ w < 2P(2h+1)a`; `a ≤ B0 ≤ θ < F·a` on every finite join |
| SC-5 `hand example L3` | §3.5 L3, `b = 2` | three nonsingleton blocks `{1, 3}`, `s_U = 3`, aligned, `t_U ∈ (1.5, 3]`, no atom, the edge of weight 5 classified in a global class that brackets it |
| SC-6 `hand example L4` | §3.5 L4, `b = 3` | atoms `{4, 5}` in three blocks, cutoff `t_U ∈ (4, 8]`, `H` = one edge of weight 1; blocks `{1, 2}` without atom; `ph` accepted by `substitute_atoms`; `validate_ordered_hierarchy` clean with `c0 ≤ 2` |
| SC-7 `unaligned block on uint64` | L3 graph with corridor weight `2^60`, `p_log2 = 0`, `b = 2` (`f = 29`) | block `{1, 3}` unaligned; all its within-window edges `< a_U/2`; θ equal on the block; `t_U` on the piece grid; `2·p_top ≤ g` |
| SC-8 `unaligned block with an atom` | L4 graph with weights 0–3: `2^63 + 1`, 3–4: 2, 4–5: 1, `p_log2 = 0`, `b = 3` (`f = 31`) | the blocks `{4, 5}` are unaligned (`F²·2 = 2^63 < a_U`); each is one atom with cutoff 2; the blocks `{1, 2}` are aligned |
| SC-9 `several grid components (double)` | dyadic weights with two corridor weights more than `3f` octaves apart | two components; `grid[i+1] = 2·grid[i]` inside each; no class crosses a component start; supports with a gap below two are merged |
| SC-10 `unattached blocks` | one pivot (`b > N`) | every nonsingleton block unattached, separate grids only, several roots in condensation order |
| SC-11 property: Lemma 8.4 | random connected graphs, `N ≤ 200`, `b ∈ {1, 2, 3, 4, 8, ⌈log2 N⌉}`, generator families of the harness, several seeds | every invariant of §2 through `check_outer_scales`: unaligned MST edges `< a_U/2`; aligned significant weights inside their supports; `s_U` = least original radius; `s_U/2 < t_U ≤ s_U`; atoms `< b` vertices, one owner, one `r`; hierarchy validates Def 4.1 with `c0 ≤ 2` |
| SC-12 property: SCC equality | small `N` (≤ 40) | for every internal node: states below = one SCC of arcs `<` class upper bound; children = SCCs of arcs `<` height (Tarjan) |
| SC-13 `alignment by median pruning equals the per-block test` | random groups | same classification; predicate tests per group `≤ ⌊log2 n⌋ + 1` |
| SC-14 `distances on the hierarchy` | SC-11 inputs | `evaluate` on the substituted hierarchy with the seed of Alg 8 l.14: minimum per vertex equals Bellman–Ford |

### 6.2 P3-OUTER-PIPELINE

Proposed files: `tests/unit/paper/outer/test_shortest_paths.cpp`, `tests/differential/paper_experimental/test_outer_forced_vs_dijkstra.cpp` (the name already in `PAPER_MAPPING.md`).

| Test | Input | Checks |
|---|---|---|
| PL-1 `Mode A never takes the outer branch` | the differential families | `takes_outer_construction` false; exactly one base-case FallbackRecord (existing behaviour kept) |
| PL-2 `outer_construction with the paper's constants` | direct call on small graphs (`N ≤ 300`) | distances equal Bellman–Ford; `outer.pivots = 1`; no macro; no FallbackRecord with an `outer:` reason |
| PL-3 `forced outer, levels L0–L4` | the five examples of §3.5 | distances and tree equal the oracle; `verify_sssp`; the level's feature is present in the counters (pivots, macro edges, nonsingleton blocks, atoms) |
| PL-4 `forced outer vs oracle` | every generator family, `N ≤ 2000`, `b ∈ {1, 2, 3, 4, 8, 16, 64}`, seeds fixed and varied | distances equal Bellman–Ford and baseline Dijkstra; `verify_sssp` (root 0, every arc feasible, tree arcs tight, parents reach `s`, unreachable = infinity) |
| PL-5 `stage checkers pass on the intermediate states` | PL-4 inputs | `check_route_names`, `check_aux_metric`, `check_calendars`, `check_atoms`, the P3-OUTER-STATES checker, `check_outer_scales`, `validate_seeded_ordered_hierarchy` all empty |
| PL-6 `fallback records` | PL-4 inputs | no outer base-case record; every record is a local base-case reason (or the D-10 reason when a cap is set); one record per Alg 4 call that takes its base case; the parameter set is recorded (`PaperParams::record`) |
| PL-7 `zero edges, unreachable vertices, loops, parallel edges` | inputs with zero cycles and several components | restoration correct; contracted `N` is what the outer branch sees |
| PL-8 `seed independence` | 16 seeds on inputs where NameRoutes samples (D-17 parameters, `t ≥ 6`) | identical distances for every seed; at least one instance with `aux_metric.radial_interowner_edges > 0` and `outer.radial_edges_deleted > 0` |
| PL-9 `source position` | source at a pivot, at a nonpivot outside any atom, inside an atom | all equal the oracle (Alg 8 l.14 redirection) |
| PL-10 `double weights` | dyadic `double`; non-dyadic `double` | dyadic: exact equality with the exact reference; non-dyadic: either equal within the documented tolerance or an exception, never a silently wrong result (Mode B failure rule, Q-63) |
| PL-11 `unaligned end to end` | SC-7 and SC-8 inputs | distances equal the oracle |
| PL-12 `sanitizer` | PL-3, PL-4 at reduced size | clean under `make sanitizer` (P4-SANITIZER) |

Suggested additions to the acceptance criteria in `docs/TASKS.md` (for the lead):
- P3-OUTER-SCALES: "`f` equals the big-integer table of note 04 §1.4 and is never derived from a saturated `h`; the unaligned branch is covered by a `uint64` test with `p_log2 = 0` and by a dyadic `double` test; `check_outer_scales` is empty on every property input."
- P3-OUTER-PIPELINE: "covers levels L0–L4 of note 04 §3.5; a direct call with the paper's constants (single pivot) equals the oracle; no FallbackRecord other than local base cases; a vertex without emission throws."
