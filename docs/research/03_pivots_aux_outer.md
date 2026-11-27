# 03. Pivot Routes, Auxiliary Distances, and the Outer Hierarchy (Paper Sec 6–8)

English | [日本語](03_pivots_aux_outer.ja.md)

Target: arXiv 2611.34911v1 "Linear Single-Source Shortest Path, Almost" (Gholami, Hajiaghayi).
Written after reading l.1367–2217 (PDF p.29–47) of the text `.cache/paper/2611.34911v1.txt` in full. Sec 2–5 were also read as background (for the detailed analysis, see the documents by the other owners).
Notation: `[Supplement]` marks an addition for implementation that the paper does not state. `[Ambiguous]` marks a place where the paper's wording does not determine a unique reading. Line numbers refer to the text version.

The quantitative treatment of thresholds and constants (minimum N, which paths real inputs take, Mode B) is collected in `thresholds.md`.

---

## 0. Key Points (for Implementers)

1. **Randomness appears only in Sec 6 (NameRoutes and Estimate)**. It is Las Vegas; there is no Monte Carlo part. The correctness of a certificate (`d_S(x,y) ≤ Z ≤ Pλ`) holds independently of the random bits (the unconditional upper bound of Lemma 6.7(1)). Randomness affects only the number of trials until one is accepted (= the time).
2. The outer construction (Alg 8 line 8 onward) runs only when `N > L_N^64`. The smallest N satisfying this is `589^64 + 1 ≈ 1.94×10^177`, so every real input goes to the base case (scan-based Dijkstra) (`thresholds.md`).
3. The construction **does not recurse**. Its depth is fixed: Alg 8 (outer) → Alg 4 (local, per atom) inside Alg 7 and Sec 8.2 → the "exact Kruskal merge tree" inside Alg 4 (atom Y, |Y| < L_k^4). Alg 8 never calls itself.
4. Constants that matter for correctness (changing them can make distances wrong): Estimate's `T ≥ 4t`, the scale top `≥ tT²a`, the window `[a/T², T²a]`, the outer shell ratio `K = 16`, Sec 7.3's `32q ≤ r < 64q`, overlap test `[r/4, 4R]` and `16R`, Sec 8.2's `F > 2P(2h+1)(8K²h)²`, and the local `K ≥ 8h` and macro test `w ≤ N³σ`.
   Constants that affect only time: the outer and local thresholds, the size of `P`, the sample sizes `16t` and `32t`, `H`, and the block length of Sec 3.
5. In a floating-point implementation, underestimating `d_e = Pλ_e`, or rounding in the comparisons scaled up by `H`, **can make the radial-edge deletion test (Sec 8.1) wrong and delete a shortest path**. Integer weights with exact comparisons (saturating shift comparisons or multi-precision) are recommended (`thresholds.md` §6).

---

## 1. Sec 6 Short Pivot Routes from Random Samples (p.29–35)

### 1.1 Definition 6.1 (Pivot route certificate, p.29, l.1368)

- **Input (pivot-routing instance)**: a loopless undirected multigraph `G_P` whose vertices are t pivots. It has M edges `e = xy`. Each edge has a positive length `a_e`, which is the length of a **physical walk** joining x and y.
- **Route certificate**: for a common factor `P ≥ 1`, an input weight `λ_e ≤ a_e` such that a physical x–y walk of length `≤ Pλ_e` exists.
- **Dictionary**: a sorted list of the distinct weights used as certificates. Each edge holds one dictionary index.
- In the application (Sec 7.1), `t = O(1 + N/L_N^64)` (N is the number of physical vertices). Loops may be discarded beforehand.

### 1.2 Preprocessing in Sec 6.1 (l.1382–1392)

- Replace each class of parallel edges by one representative of minimum weight. Keep a pointer from each original edge to its representative. The representative's certificate is valid for every edge in the class (`a_e ≥ a_rep`).
- The representative graph is simple, with `≤ t²` edges. This serves only to bound the number of sampling rounds; the first and last scans are charged to the original M.
- **Median pruning**: find the least a* in a group with `P·a* ≥ Z`. Select the median m by BFPRT (Lemma 3.1), form `P·m` by `log2 P` doublings, and compare it with Z. If the comparison holds, keep m as a candidate and move to the keys **strictly smaller** than m. Otherwise move to the keys **strictly larger** than m (keys equal to m are discarded in both cases). The least remaining candidate is the answer.

### 1.3 Algorithm 5 `NameRoutes(G_P)` (p.30, l.1395–1418)

**Input** an instance of Definition 6.1. **Output** a sorted dictionary, a dictionary index per input edge, and the factor P.

```cpp
// Alg 5 l.2-3
u = ceil_log2(t + 2);  g = 4*u + 3;
P = pow2_ge(8 * 256 * g * u * u);          // smallest power of 2 at least 2048·g·u² (computed in integers)
if (t < 2) return {empty_dict, {}, P};     // l.1422: no edges exist
// l.4: choose representatives (stable counting sort by endpoint pair (min,max) → per pair, least a, ties by least id)
rep[o] = representative of original edge o; U = list of representatives;
vector<Rec> dict; vector<int> ptr(|U|);    // ptr: representative → dictionary record
while (U.size() > 32 * t) {                // l.5
  S = random_subset(U, 16 * t);            // l.6: partial Fisher-Yates in O(|U|) expected time
  auto [table, idx, disc] = Estimate(S, endpoint_pairs(U));   // l.7
  groups = bucket_by(idx, U where !disc);  // l.8: counting sort by index
  qualified = {};  astar_of_group = {};
  for (grp : groups) {                     // l.9
    Z = table[grp.index];
    optional<W> a = median_prune(grp.weights, Z, P);   // l.10
    if (a) for (e : grp) if (a_e >= *a) qualify(e, *a);  // l.11-12
  }
  if (2 * qualified.size() >= U.size()) {  // l.13: accept if at least half qualify
    for (grp with a*) { r = dict.push(a*); for (e qualified in grp) ptr[e] = r; }   // l.14
    U.remove(qualified);
  } // otherwise discard all proposals and retry with fresh random bits (l.1421)
}
for (e : U) ptr[e] = dict.push(a_e);       // l.15: the rest use their own weights as certificates
sort(dict by weight); coalesce ties → newidx;   // l.16
for (o : original edges) name[o] = newidx[ptr[rep[o]]];
return {dict, name, P};                    // l.17
```

- **Invariant**: after each accepted round, every edge removed from U holds a valid certificate (l.1421).
- **Randomization**: S is a uniform 16t-element subset of U. Inside Estimate, an independent uniform permutation and a uniform integer j are also drawn per scale. A failed round is discarded entirely and redrawn with fresh random bits. **Las Vegas**. The acceptance probability of one round is ≥ 1/2 (Lemma 6.8(2)).
- **Computational model**: weights use only comparison, addition, and doubling (`P·m` by doubling, BFPRT by comparison). Indices and grouping use integer counting sort. Random numbers use the rejection method of Sec 2 (draw `⌈log2 h⌉` bits and redraw if the value is at least h).
- **Degenerate cases**: if `t < 2`, empty. If `|U| ≤ 32t`, no sampling is done and every edge is treated as terminal. Median pruning works correctly even when all edges have the same weight (equal keys are discarded on both sides).
- [Supplement] The "at least half" in `2·|qualified| ≥ |U|` is the paper's "at least half" compared directly in integers.

### 1.4 Lemma 6.2 (Critical edges, p.30, l.1424)

- Claim: let `g = 4⌈log2(t+2)⌉ + 3`. An edge is **critical** if, after deleting it, its endpoints can no longer be joined by a path with "at most g edges, each of weight ≤ the weight of the deleted edge". There are **≤ 2t** critical edges.
- Proof outline: the critical subgraph is simple and has no cycle of length ≤ g+1. Repeatedly delete vertices of degree ≤ 2 and charge them. The remaining subgraph has minimum degree ≥ 3, so a BFS to depth `j = ⌊g/2⌋` produces `3·2^{j−1} > t` vertices, a contradiction.
- Check: `3·2^{2u} > t` always holds (verified exhaustively for t < 2×10^5).
- **Test**: on small random multigraphs, count critical edges by brute force and check `≤ 2t`.

### 1.5 Lemma 6.3 (Distances in a uniform sample, p.30, l.1439)

- Claim: let the residual U have `h > 32t` edges and let S be a uniform 16t-edge subset of it. The expected number of edges e with `d_S(x,y) > g·a_e` is ≤ h/8 (including disconnected pairs).
- Proof: backward analysis. If e is chosen uniformly from U\S, each edge of S∪{e} is equally likely to be e. There is no witness exactly when e is critical in S∪{e}. The expectation is ≤ `(h−16t)·2t/(16t+1) < h/8` (checked numerically).
- **Test**: over many random seeds, check empirically that the mean fraction with `d_S > g·a_e` is ≤ 1/8.

### 1.6 Sec 6.2 Algorithm 6 `Estimate(S, Q)` (p.31–33, l.1455–1508)

**Input** a positive-weight sample graph S (O(t) edges) and h requests with distinct endpoints. **Output** a table of upper bounds, and per request either a table index or a disconnected flag.

```cpp
T = pow2_ge(4 * t);                         // l.1473
sort S by (w, id);                          // merge sort of adjacent runs, O(t log t)
// band: the least unassigned weight a starts a band, which contains the weights a ≤ w < T·a (T·a by log2 T doublings)
// Kruskal sweep: each union creates a binary parent node (weight, band id) → merge forest (l.2)
for (q : Q) if (find(x) != find(y)) disc[q] = true;
            else band_of[q] = band(LCA_mergeforest(x, y));     // l.3 (Lemma 3.9)
for (a : anchors of bands with demand, in increasing order) {   // l.4
  // quotient Q_a (l.1481–1488): contract edges of weight < a/T² (tested as T²w < a).
  //   Insert the selected forest edges when they fall below the contraction cutoff.
  //   Labels are maintained by relabeling the smaller list. Two cursors maintain the window [a/T², T²a].
  //   Kept edges: those inside the window whose contracted endpoints differ. Vertices: only those incident to kept edges.
  build Q_a; map demand endpoints by component label;            // l.5
  ell = least ℓ with 2^ℓ ≥ t·T²;   // scales Δ_i = 2^i·a, i = 0..ell
  H = pow2_ge(pow(16*t, 5));
  for (i = ell; i >= 0; --i) {                                    // l.6-7
    perm = uniform_permutation(V(Q_a));  j = uniform_int(0, H-1);
    // radius R = Δ(1 + j/H). Compare as H·dist ≤ (H+j)·Δ (l.1493)
    best[] = +∞ (lazy initialization by timestamps); label[] = none;
    for (c : perm) radius_limited_dijkstra(c, R):                 // binary heap
        // When vertex z is settled at distance d, only if d < best[z] (strictly)
        // set best[z] = d and scan its edges. Otherwise discard that visit.
        // On the first improvement, label[z] = c (permanent).
  }
  // l.8-9: in each quotient component, build a partition tree from coarse to fine scales.
  //   A part at scale Δ = a pair (parent part, label). Group by stable counting sort.
  //   For each new part, write the value 4Δ + a (= a·2^{i+2} + a) to the global table.
  //   Attach one leaf per vertex below the finest parts. LCA preprocessing.
  for (q in band a) idx[q] = table_index(LCA(leaf(x'), leaf(y')));   // l.10
}
return {table, idx, disc};                                         // l.11
```

- **Invariant**: vertices in a part at scale Δ have the same center at that scale (l.1507). At the top scale `Δ ≥ tT²a`, the first center covers the whole quotient component (Lemma 6.7(1)).
- **Computational model**: weights use comparison, addition, and doubling (T², multiplication by H, `(H+j)Δ`). H and j are integers, multi-word integers of about `log2 H = 5·log2(16t)` bits (`thresholds.md` §5). Dijkstra uses a binary heap (Lemma 3.2, a log per operation).
- [Supplement] The definition of `ell`, "the least ℓ_a such that the last scale is `≥ tT²a`", is equivalent in integers to `2^ℓ ≥ t·T²` (since a > 0).

#### Lemma 6.4 (Lifting a quotient path, p.32, l.1509)

- Claim: if x, y are distinct connected sample vertices whose merge weight is β and whose band anchor is a, then their images in Q_a are distinct, exist, and are connected. For the quotient distance d, `a ≤ β ≤ d ≤ d_S(x,y) ≤ d + a`.
- A deterministic claim. The facts used are `(t−1)β < T²a` (`T > t−1`) and that lifting adds length `t²a/T² ≤ a` (`T ≥ t`).
- **Test**: on random S, check the four inequalities for every connected request.

#### Lemma 6.5 (Separation at one scale, p.33, l.1530)

The probability that two points at distance d receive different centers at scale Δ is `≤ H_t(d/Δ + 1/H)`, where `H_t = Σ_{i≤t} 1/i`. The proof is the same ordering argument as in FRT.

#### Lemma 6.6 (Work to assign one scale's labels, p.33, l.1542)

The pruned searches return each vertex's "first covering center" (correctness by induction, deterministic). The expected time is `O((v+e)·log²(t+2))`. The expected number of scans per vertex is ≤ H_t, with a log per heap operation. This also holds conditioned on a fixed radius.

#### Lemma 6.7 (Sample estimator, p.33–34, l.1560)

- Claim: t ≥ 2 pivots, O(t) edges. The table has `O(t log(t+2))` entries, and the expected time is `O(h + t log³(t+2))`. Disconnected pairs are identified. For connected requests,
  `d_S(x,y) ≤ Z` (**deterministic**), `E[Z | S] ≤ 256u²·d_S(x,y)`, `u = ⌈log2(t+2)⌉`.
- Proof essentials: (1) within a part, the distance to the common center is < 2Δ, so `d_S ≤ 4Δ + a` (by Lemma 6.4). (2) `Z ≤ 5a + 8ΣΔ·I_Δ`, `E[Z|S] ≤ 5a + 8H_t((ℓ_a+1)d + ΣΔ/H)`. Since `H_t ≤ 2u`, `T < 8t`, `ℓ_a+1 ≤ 7u` (checked numerically for t < 10^5), and `ΣΔ < 4tT²a`, `H ≥ (16t)^5` gives `ΣΔ/H ≤ a ≤ d`. (3) Each edge enters ≤ 5 windows, so the total size of the quotients is O(t).
- [Note] `ΣΔ/H ≤ a` needs only `H ≥ 4tT²`, so `(16t)^5` is a value with slack (the basis for lowering H in Mode B; H affects only time).
- **Test**: (a) Z ≥ d_S for every request (deterministic, so every time). (b) `mean(Z/d_S) ≤ 256u²` empirically. (c) The disc flags match the connected components of S. (d) Number of table entries ≤ c·t·log t.

### 1.7 Lemma 6.8 (Pivot route certificates, p.34–35, l.1605)

- Claim: on a loopless instance with `t = O(1 + N/L_N^64)`, Alg 5 gives every edge a certificate in the sense of Definition 6.1. The dictionary has `O(t log²(N+2))` entries. P is a power of 2 of order `log³(t+2)`. It **terminates almost surely**, in expected time `O(M + N)` (including the sorted dictionary).
- Proof essentials:
  1. For an accepted edge, `a* ≤ a_e` and `d_S ≤ Z ≤ Pa*`. Concatenating the physical walks represented by a shortest path in S gives the proof. A terminal edge is certified by its own weight.
  2. In one trial, for requests with `d_S ≤ g·a_e`, Markov's inequality and Lemma 6.7 give `Pr[Z > P·a_e | S] ≤ 1/8` (because `P ≥ 8·256·g·u²`). The expected number of failures is ≤ h/8 + h/8 = h/4. The probability of rejection (failures > h/2) is ≤ 1/2.
  3. The cost of one trial is `O(h + t log³ t)`. The probability of reaching the i-th trial is ≤ 2^{1−i}.
  4. Each accepted round halves the residual. There are `O(1 + log t)` accepted rounds, and ≤ 32t terminal edges. Sorting the dictionary takes `O(t log³ t)`. `O(t log⁴ t) = O(N)` uses the pivot-count assumption `t = O(N/L^64)`.
- **Randomization class**: Las Vegas. On failure (rejection), it retries on the same residual with fresh random bits. There is no bound on the number of retries (it terminates almost surely).
- **Test**: for every edge e, `λ_e ≤ a_e` and `dist_{G_P}(x,y) ≤ P·λ_e` (with Dijkstra on G_P as the reference). The dictionary is sorted and has no duplicates. Statistics of the number of trials and the acceptance rate (that it is ≥ 1/2). All invariants hold across random seeds.

---

## 2. Sec 7 Auxiliary Distances to the Owner Pivots (p.35–41)

### 2.1 Sec 7.1 Edge Classification and the Auxiliary Metric (l.1642–1660)

- The outer level uses the owner decomposition of Sec 5 with **separation size `b = L_N^64`**. `r_v = b_T(v, p(v))` (the maximum edge on the tree path), `R_v = D_T(v, p(v))` (the tree distance), `h = b + 1`. The macro rule of this section replaces the local one (`w ≤ N³σ`).
- For each edge `xy` crossing owners, the projected length `A_e = R_x + w + R_y ≤ (2h+1)w` (Lemma 5.2) is passed to Lemma 6.8 as the edge length of a pivot-routing instance. The vertices are the pivots, and `t = |A|`.
- From the name λ_e, form `d_e = P·λ_e` (once per name, by log2 P doublings).
- **Classification**: an edge crossing owners is **radial** if `2d_e ≤ w`, and **macro** otherwise. An edge within an owner is radial, with `d_e = 0`.

#### Lemma 7.1 (p.35, l.1653)

Every cut corridor edge is a macro. A macro of weight w has a name λ with `λ/(2h+1) ≤ w < 2Pλ`. For a cut edge f, `d_f ≥ w(f) > 0` follows from the cut property of the MST.
**Test**: the two inequalities for every macro edge. All corridor edges are macros.

### 2.2 Sec 7.2 The Owner Coordinate (l.1662–1690)

- Identify all pivots into one point `a`. Keep the radial edges with auxiliary length `w'_e = w_e + d_e`, and keep all spokes of the owner MSTs. These spokes form the MST `T'` of the new graph. Discard positive-length loops at a. `J_v = d_{w'}(v, a)`. At pivots, `J = 0`.

#### Lemma 7.2 (Physical interpretation of the owner coordinate, p.35, l.1668)

- Claim: `J_v` is an upper bound on the length of a physical route from v to its owner pivot. `r_v ≤ J_v ≤ R_v ≤ h·r_v`. On radial edges, `|J_x − J_y| ≤ w'_e ≤ 3w_e/2`. `T'` is the MST after the identification.
- **Test**: compute J by reference Dijkstra on G' and check the inequalities on every vertex and every radial edge. That `T'` is an MST (whether its weight sum matches Kruskal's).

### 2.3 Sec 7.3 Width Classes, Calendars, and Atoms (l.1692–1779)

Data structures and procedure (l.1693–1710):

```cpp
// r_v, R_v top-down on T': a child v of a has r=R=w(av). Otherwise r_v = max(r_u, w), R_v = R_u + w
// tree pieces P: keep the parent-child edges (u,v) of T'−a whose intervals [r/4, 4R] overlap
//   overlap ⇔ r_v ≤ 16R_u and r_u ≤ 16R_v (×16 is 4 doublings)
// graph pieces C: apply the same overlap test to the auxiliary edges of G'−a and group by union-find (each P lies in one C)
// root of P: the vertex closest to the common pivot
for (C : graph pieces) {
  origin = min_{v∈C, r_v>0} r_v;
  grid = { origin·2^i : i ≥ −7 (or below), ..., max R + 2 octaves };   // l.1697-1698
  // Lemma 7.3: root a spanning tree of the overlap edges of C at a vertex of least radius.
  //   Binary-search the child's r within O(1+log(h+2)) entries from the parent's grid index.
  idx_r[v] = grid index i (grid[i] ≤ r_v < grid[i+1]);
  q_v = grid[idx_r[v] − 5];                 // 32q ≤ r < 64q
  U_plus[v] = first grid value strictly greater than R_v;
}
for (each pair (C,q)) calendar: cells [iq,(i+1)q), i = 0..(largest U+)/q;   // Lemma 7.4
within C, merge the markers of all q from coarse to fine into one increasing chain;
// edge uv crossing graph pieces (oriented with r_v > 16R_u): marker arc U+_u → 16q_v
// consecutive-marker arcs of each chain + the arcs above → DAG. Process in topological order
// atoms: connected components Y (union-find) of the non-pivot edges with w' < 2·max(q_u,q_v)
//   the physical spanning tree of Y (defining edges) and its weight sum L_Y (≤ 2q|Y|)
// initialization: place the feasible upper bound R_v at each v, plus the offers from edges incident to s
```

#### Lemma 7.3 (Octave spans of overlapping ranges, p.36, l.1711)

Over all graph pieces, the number of dyadic grid points from the least r to the largest R (including a constant number at both ends) is O(N). Each vertex's width and upper marker index are found in `O(1 + log log(b+4))` time.
Proof essentials: at the root of a tree piece, `R_x ≤ 17r_x/16`. On a kept edge, `R_v ≤ 17R_u`. Merging by graph overlap adds only a constant number of octaves per joining edge.

#### Lemma 7.4 (Counting width cells, p.36, l.1736)

The total number of cells over the (C, q) calendars, and the time to merge the markers, is O(N). A tree edge of weight w contributes `Σ_{q>w/64} w/q = O(1)`.

#### Lemma 7.5 (Atom size at an auxiliary width, p.37, l.1751)

An atom belongs to one graph piece and one width, and has **< b vertices**. Its defining edges give a spanning tree with < 2q per edge. Reason: the endpoints of a defining edge have equal r, and the path in T' does not pass through a pivot, so it lies in one component with the pivots removed (< b vertices).

#### Lemma 7.6 (p.37, l.1762)

The grids and calendars have O(N) entries in total. Construction takes `O(m + N log log(b+4))`. Atoms have < b vertices. **Every marker arc increases the numerical time** (`U+_u ≤ 2R_u < r_v/8`, `16q_v > r_v/4`).

**Tests for Sec 7.3 (intermediate state)**: `32q_v ≤ r_v < 64q_v` and `R_v < U+_v ≤ 2R_v` at every vertex. `r_v > 16R_u` on edges crossing pieces. The marker DAG is acyclic and every arc points in the direction of increasing value. Total cells / N is bounded. For every atom, `|Y| < b`, q and C are uniform, `L_Y ≤ 2q|Y|`, and `T'[Y]` is connected.

### 2.4 Sec 7.4 Local Searches, Clipping, and External Arrivals (l.1780–1963)

**Families and activation lists** (l.1781–1800):
- An atom belongs to the family of its representative's (P, q). Within a family, order the representatives by preorder of T' (integer sort). The gap `g(i,j)` between consecutive representatives is computed as a tree distance (LCA and distances from the root).
- The family horizon is `max_{Y∈family}(R_rep(Y) + L_Y)`. For `v ∈ Y`, `R_v ≤ R_rep(Y) + L_Y`.
- Walk the family's (C, q) calendar from cell 0 until "start ≤ horizon", or to the end of the calendar. Push the family record onto the activation list of each cell's starting marker.

**Run management within a family** (l.1801–1815):
- A linked list of registered atoms, and contiguous runs of unregistered ones. In cell `[a0, b0)`, split the run at the atom i containing the minimum z that first becomes due. For each newly created subrun (at most 2), offer `u = max{b0, z + L_i + g(i,j)}` to the native inward copy of the representative v of the boundary atom j adjacent to i. Discard it if `u > R_v` or if it does not lower the pending key.
- Splitting and extracting minima (split–findmin) take `O(m + Nτ(N))` in total (Lemma 3.11).

**Local representation of an atom** (l.1816–1833):
- Each atom Y uses the whole induced graph `G'[Y]` with auxiliary lengths (including heavy chords). The exits leaving the atom are listed at their tails.
- Suppressing exits: do not emit from tail u if the emission time is `> R_u`. Discard an external offer to head v if it is `> R_v`. Of the edges crossing graph pieces, keep only the orientation from the side with the smaller range to the side with the larger range.
- The family's registration time z and L_Y determine the persistent deadline `z + L_Y` (inclusive).
- Y is a replaceable atom with cutoff q (Definition 5.11). By Lemmas 5.13 and 5.14, build a local hierarchy capped to a single root, and normalize it by Lemma 4.3. The internal units are ≤ q.

**Visiting active atoms** (l.1835–1841): for a single-vertex root, emit the due minimum once, relax the permitted exits, and retire the leaf. For an internal root, open it at the first due key D0 and set the cap `c = min(atom deadline, D0 + S(root))` (inclusive). In each cell, call `Scan(root, I, c)` on the nominal root slots whose intersection with the cell contains a time ≤ c. When resuming a clipped slot, start from the first live child rank, keeping registrations, offers, and deadlines intact.

#### Algorithm 7 `OwnerDistances(T', G')` (p.39, l.1843–1867)

```cpp
// l.2  r_v, R_v, tree pieces and graph pieces, calendars, marker arcs crossing pieces (Lemma 7.6)
// l.3  components Y of the non-root edges with w' < 2max(q_u,q_v). Discovery tree (weight sum L_Y) and representative
// l.4  group atoms by the representative's (tree piece, width)
// l.5  family gaps and activation lists (preorder of representatives, tree queries of Lemma 3.9)
// l.6  collect G'[Y] and the interatom exits in one edge scan. Run LocalHierarchy (Alg 4) on each G'[Y] and cap it by the width
// l.7  normalization of each hierarchy, child order, representatives, gaps (Lemma 4.3, 4.5)
// l.8  leaf layout: family order → within a family, T' preorder of atom representatives → local child order. All keys ∞, split at family boundaries
// l.9  place R_v at each native inward copy. Place all direct offers from the root a
for (marker : topo_order(marker_DAG))                 // l.10
  for (fam : activation[marker]) {                    // l.11
    register_due_atoms(fam);                          // l.12: run splitting (Lemma 3.11) + padded offers
    for (atom : active(fam)) visit_atom(atom, cell);  // l.13: emit for a single vertex, Scan the root slots for an internal one
  }
return the least emitted value at each vertex (0 at the root);   // l.14
```

#### Lemma 7.7 (Slots clipped by the calendar, p.39–40, l.1871)

- Assumptions: a seeded ordered hierarchy with one root and units ≤ q. External offers are feasible, arrive before the cell of their timestamp is visited, and do not enter below the left end of the current cell. The first cell starts at or below every finite initial timestamp. The visited cells cover all protected timestamps (including those equal to the deadline).
- Conclusion: every protected prefix is emitted in its own cell. The work is `O(M + Nτ(N) + Q + E + C)`. A nominal slot crosses an inner boundary of the external grid at most once, so each local scan is needed at most twice (Figure 2).

#### Lemma 7.8 (Time order of external offers, p.40–41, l.1919)

Keep an interatom offer `z + w'` only in a permitted orientation and only when `z ≤ R_u` and `z + w' ≤ R_v`. An effective decrease arrives before the starting marker of the destination cell. A tight interatom step passes all orientation and endpoint tests. The proof is by cases (across pieces, the arc `U+_u → 16q_v`; within the same piece, `w' ≥ q_u + q_v`).

#### Lemma 7.9 (p.41, l.1944)

Once the names are given, Alg 7 computes all J_v **deterministically** in `O(m + Nτ(N))` time. `τ(N) = 1 + log log log(N+16)`.
Correctness follows by induction over the marker order, with Lemmas 5.14, 7.7, and 7.8.

**Tests for Sec 7.4**: J matches reference Dijkstra on G'. Check "arrives before the starting marker of the destination cell" against a log of arrival times (Lemma 7.8). Family horizon `R_v ≤ R_rep + L_Y`. Emitted values never fall below the lower bound (the true J) (the "no underestimate" of Lemma 4.6). Degenerate inputs such as a single atom, a single family, or all vertices in one cell.

---

## 3. Sec 8 The Outer Hierarchy and the Single-Source Algorithm (p.41–47)

### 3.1 Sec 8.1 The Outer State Graph and Threshold Components (l.1966–2001)

```cpp
// Alg 8 l.12: delete radial edge e=xy if w_e > J_x + J_y + d_e (spokes + detour is strictly shorter)
//   owner tree edges remain (J_v ≥ r_v ≥ w_e). Surviving radial edges have w_e ≤ 2(J_x + J_y)
// two-shell construction (Sec 5.2) with coordinate J and K = 16:
//   near: radial edges whose endpoints both have positive J with ratio ≤ K² = 256. Per near component, a base-16 grid with the least positive J as its base
//   native shell i(v): a·16^i ≤ J_v < a·16^{i+1}
//   inward state j is shells {j, j+1}, outward state j is shells {j−1, j}. Near edges inside the window are bidirectional
//   lower the shell record when inward / raise it when outward → to a new native copy
//   far inward descends to native inward, far outward rises to native outward
//   zero switch: inward copy → native outward. A pivot (J=0) has two pivot copies
//   macro: leaves any ordinary copy and returns to the head's native inward
// special states (v_in, v_out for each non-pivot v):
//   ordinary inward copy(v) → v_in (0). v_in → parent(v)_in (length of the physical owner edge) …
//     → ordinary inward pivot copy
//   ordinary outward pivot copy → child_out → … → v_out → ordinary native outward(v) (0)
//   tails and heads of macro arcs are ordinary copies only. Special states are singleton blocks
// propagate I0/O0 (macro minimum ignoring weights) and I*/O* (bottleneck minimization) in block-DAG order → θ = max(I*,O*)
// joins: auxiliary pivot vertices, corridor joins (edge weights of T_A), MST edges of local blocks, state → owner joins (at θ)
// positive ordinary block U: s_U = min_{v∈U} r_v. The local cutoff t_U is determined in Sec 8.2
```

#### Lemma 8.1 (p.42, l.2002)

- Claim: surviving radial edges have `w ≤ 2(J_x + J_y)`. If a globally shortest radial run is split at its least coordinate, the inward side and the outward side each swing by a factor of at most 16. The ordinary graph realizes every geodesic. The augmented quotient without macro arcs is acyclic.
- Key: an interval reaching vertex M with both ends ≤ M/16 has length `≥ (2/3)(M − M/16) = 5M/8`. The detour costs `≤ M/8 + L/2 < L`.
- **K = 16 is required for correctness** (below 16 the realization argument fails).
- **Test**: for each radial run of a reference shortest path tree, the swing ratio is ≤ 16. The "minimum over copies" of Dijkstra on the state graph matches the physical distance. DFS on the block quotient without macros finds no cycle.

#### Lemma 8.2 (Bottleneck of access to a macro edge, p.43, l.2027)

`B0 = max(I0, O0)` and `θ = max(I*, O*)` are finite together, and when finite, `r ≤ θ` and `B0 ≤ θ ≤ 4K²h·B0` (1024h since K = 16).
**Test**: at every state, check the inequalities with θ and B0 computed from brute-force bottleneck-minimizing paths.

#### Lemma 8.3 (p.43–44, l.2053)

For any W > 0, the block MST edges, the state → owner joins at θ, and the corridor edges of weight ≤ W exactly describe the state components strongly connected by arcs of length ≤ W. The strict-inequality version (< W) also holds.
**Test**: for random W, Tarjan's SCC (arcs ≤ W) matches the union-find components of the joins restricted to states.

### 3.2 Sec 8.2 Numerical Scales and Atom Substitution (p.44–45, l.2076–2142)

```cpp
// replace each macro's name by "the least physical weight of any macro with that name" (ratio within one name ≤ 2P(2h+1))
F = pow2_gt(2*P*(2*h+1) * sq(8*K*K*h));      // once, with integer parameters (l.2080)
// supports: an interval [a, F³a] for each finite join name a and each corridor weight a
// singleton blocks (including pivots and special states) contribute only their own finite θ join. No internal merge
for (positive nonsingleton ordinary block U with a finite attachment) {
  z_U = state minimizing θ; a_U = its name;
  aligned = (F²·s_U ≥ a_U);                   // doubling and comparison (median pruning, O(log F) doublings)
  if (aligned) start of the name's interval = min(start, s_U);
}
// sort the interval endpoints and merge overlaps and gaps of ratio < 2 → a doubling grid per merged component
// aligned:   t_U = s_U rounded down on the global grid of names
// unaligned: t_U = rounded down on a separate grid (using pointers into the auxiliary grids of Sec 7.3, extended by a constant number of octaves at both ends)
//   in both cases s_U/2 < t_U ≤ s_U
// binary-search events (block MST edges ≥ t_U, θ joins, corridors) into integer grid classes (O(log log F) = O(τ(N)))
// simultaneous merging per class (Lemma 3.6) → partition hierarchy (Def 5.12)
//   suppress auxiliary nodes, and topologically sort the "arcs below the height" at each node
//   physical connection trees are built by the child quotient search of Sec 5.4 (connections of class [a,2a) are arcs < 2a, by the strict version of Lemma 8.3)
// components of edges < t_U within a block = atoms (< b vertices) → Lemma 5.15: run Alg 4 on each atom, cap it, and substitute
```

#### Lemma 8.4 (p.44–45, l.2105)

With `b = L_N^64` and names satisfying Lemma 6.8, the meaningful local events (block MST weights at or above the block cutoff) and the finite state → pivot joins can be classified in `O(m + Nτ(N))` (including the global grid). An atom below the threshold lies in one component with the pivots removed, and has < b vertices.
Proof essentials: (1) if unaligned, the edges in the window are `≤ 4K²h·s_U < a_U/2`. If aligned, `a_U/F² ≤ s_U ≤ F·a_U`. (2) The number of octaves is `O(k log F)`, `k = O(t log² N)`. (3) Edges with `w < t_U` change neither the owner nor the radius.

- **Test**: all MST edges of an unaligned block are `< a_U/2`. The meaningful weights of an aligned block fall in the support interval (plus the lower rounding boundary). `s_U/2 < t_U ≤ s_U`. Atoms have < b vertices, one owner, and one r. The resulting hierarchy satisfies Definition 4.1 (the quotient of light arcs between children is acyclic, connection length ≤ c0·h, c0 ≤ 2).

### 3.3 Sec 8.3 The Single-Source Algorithm and Main Proof (p.45–47)

#### Lemma 8.5 (Restoring a predecessor tree, p.45, l.2145)

From the distances obtained after contracting the components of zero edges, restore all original distances and a shortest path tree in `O(m + n)`. Follow the tight arcs `d(u) + w(uv) = d(v) < ∞` with orientation, and build the tree from the first discovered arcs. The output distances are recomputed by adding the original edge weights while descending the tree.
**Test**: every tree edge is tight, the tree spans all reachable vertices, and the sums from the root match the output distances. Inputs containing zero cycles.

#### Algorithm 8 `ShortestPaths(G, s)` (p.46, l.2167–2196) and the Overall Flow

```cpp
// l.2  contract the components of zero edges (BFS / union-find). Discard loops and keep only the component of s. Keep the original edge ids
// l.3  N = number of remaining vertices
if (N <= pow_sat(L(N), 64)) {      // l.4  ← always true for a 64-bit N (thresholds.md)
  // l.5 keep only the lightest edge per vertex pair, scan-based Dijkstra (Lemma 3.2, O(N² + m))
  // l.6 restore the original vertices. Unreachable is ∞. Follow finite tight arcs from s (Lemma 8.5)
  return;                          // l.7
}
// l.8  MST by the three-level growth of Lemma 3.5 → owner cells (b = L^64) → delete non-tree edges strictly dominated by the tree
// l.9  NameRoutes (Alg 5) on the projected edges crossing owners          ← the only randomized stage
// l.10 macro / radial classification. Identify the pivots into the root to get G' (radial, length w+d). Discard loops. T' = owner forest after identification
// l.11 J by OwnerDistances (Alg 7)      ← internally, LocalHierarchy (Alg 4) per atom
// l.12 delete radial edges with w_e > J_x + J_y + d_e
// l.13 outer states, threshold joins, grids, capped local hierarchies (Sec 8.1–8.2, Alg 4 per atom)
// l.14 seed 0 at the native inward outer copy of s (if substituted, at its native inward local copy)
// l.15 Evaluate (Alg 3). Minimum of the emitted values per vertex
// l.16 restore the original vertices. Unreachable is ∞. Follow tight arcs from s
// l.17 return the distances and the tree of first discovered arcs
```

**Places that enter a base case**

| Place | Condition | Base case contents | Complexity |
|---|---|---|---|
| Alg 8 l.4 | `N ≤ L_N^64` | lightest edge per vertex pair → scan-based Dijkstra → Lemma 8.5 | O(N² + m) (treated as N being bounded by an absolute constant) |
| Alg 4 l.2 (per atom in Alg 7 l.6 and Sec 8.2) | `k ≤ 2L_k^4` | Kruskal hierarchy with one state per vertex, arcs in both directions, and simultaneous merging of equal-weight MST edges (Sec 2) | O(m + k log k) (sorting; treated as k being bounded by an absolute constant) |
| Atoms inside Alg 4 (Sec 5.4) | components Y of edges < t_U within a block (< L_k^4 vertices) | exact local Kruskal merge tree (sorting < b MST edges) | O(N log(b+2)) in total |
| Alg 5 l.5 | `|U| ≤ 32t` | the rest use their own weights as certificates | linear |
| Alg 5 | `t < 2` | empty dictionary | O(1) |

**Breakdown of the expected complexity `O(m + n log log log n)`** (proof of the main theorem, l.2206–2216)

| Stage | Basis | Cost |
|---|---|---|
| l.2 contraction and reachability, l.16 restoration | graph traversal, Lemma 8.5 | O(m + n) |
| l.8 MST | Lemma 3.5 (Fibonacci heap, three levels, B1 = max(2,⌈log L⌉), B2 = L) | O(m + Nτ(N)) |
| l.8 owner cells and deletion of dominated edges | Sec 5.1, LCA | O(m + N + \|A\| log \|A\|) |
| l.9 NameRoutes | Lemma 6.8 (**expected**, the only randomness) | O(M + N), t log⁴ t = O(N) from t = O(N/L^64) |
| l.10 classification and d_e | log P doublings per name | O(m + t log² N · log log N) |
| l.11 J | Lemma 7.9 (deterministic). Local hierarchies O(m + N log log(b+4)), split–findmin O(m + Nτ(N)) | O(m + Nτ(N)) |
| l.12–13 outer representation | Sec 8.1 (the grid is linear since K = 16), Lemma 8.4, 5.15 | O(m + Nτ(N)) |
| l.15 evaluation | Lemma 4.7 | O(m + Nτ(N)) |

Here `log log(b+4)` with `b = L^64` is `6 + log log L + o(1) = O(log log log N) = O(τ(N))` (l.2211–2212). The triple logarithm comes from two sources: (i) `log log k` per vertex in the local construction (k ≤ polylog N), and (ii) binary search of the global events over a range of `O(log log N)` powers of 2.

---

## 4. List of Invariant Tests (Intermediate State)

| Target | Invariant | Source |
|---|---|---|
| pivot / owner | components with the pivots removed have < b vertices, ≤ 2 pivot neighbors. `R_v ≤ h·r_v` (h = b+1). `\|A\| = O(1+N/b)` | Def 5.1 |
| certificate | `λ_e ≤ a_e`, `dist_{G_P}(x,y) ≤ Pλ_e`. The dictionary is sorted with no duplicates | Def 6.1, Lemma 6.8 |
| estimator | `Z ≥ d_S` (every time). `a ≤ β ≤ d ≤ d_S ≤ d+a` | Lemma 6.4, 6.7 |
| classification | macro: `λ/(2h+1) ≤ w < 2Pλ`. corridor ⊂ macro | Lemma 7.1 |
| J | matches reference Dijkstra. `r ≤ J ≤ R ≤ hr`. `\|ΔJ\| ≤ 3w/2` | Lemma 7.2 |
| calendar | `32q ≤ r < 64q`, `R < U+ ≤ 2R`, the marker DAG is acyclic and increasing, number of cells = O(N) | Lemma 7.3–7.6 |
| atom | `\|Y\| < b`, (C, q) uniform, `L_Y ≤ 2q\|Y\|`, T'[Y] connected | Lemma 7.5 |
| deadline | the atom deadline `z + L_Y` covers every J in the atom. Family horizon `R_v ≤ R_rep + L_Y` | Sec 7.4, Lemma 7.9 |
| arrival order | kept interatom offers arrive before the marker of the destination cell | Lemma 7.8 |
| outer states | `w ≤ 2(J_x+J_y)` on surviving radial edges. Swing ≤ 16. Minimum over copies = physical distance | Lemma 8.1 |
| θ | `B0 ≤ θ ≤ 1024h·B0`, `r ≤ θ` | Lemma 8.2 |
| threshold SCC | union-find of the joins = SCC of arcs ≤ W | Lemma 8.3 |
| scale | if unaligned, MST edges < a_U/2. `s_U/2 < t_U ≤ s_U` | Lemma 8.4 |
| output | all distances = reference Dijkstra. The tree is tight and spans the reachable vertices | Lemma 8.5, Thm 1 |

In Mode A the outer stage never runs on a real input. Tests therefore need to build its intermediate states with Mode B arguments (force the outer stage and make b small) (`thresholds.md` §7).

---

## 5. Ambiguities and Caveats

1. [Ambiguous] Sec 8.1 sets `s_U = min_{v∈U} r_v` (l.1999), whereas the local version in Sec 5.4 sets `s_U` to "the lower end of the lowest shell of the window" (l.1184). Since Lemma 8.4(1) also uses "the v of least radius" at the outer level, reading the outer level as `min r` is consistent.
2. [Ambiguous] The grid of Sec 7.3 runs "from at least 7 octaves below the origin to 2 octaves above the largest R" (l.1698), but the proof of Lemma 7.3 says only "a constant number of octaves". Since q = grid index − 5, at least 5 octaves are needed below. We adopt 7 as is.
3. [Ambiguous] The role of the "2" in the local threshold `N > 2L^4` is not stated. What the inequalities of Sec 5.4 (N ≥ 32, h ≤ N, K < 16h) actually need is `N ≥ 28 ∧ h ≤ N` (⇔ `N > L^4`), which is weaker (`thresholds.md` §4).
4. [Ambiguous] The statement of Lemma 6.8 gives the dictionary size as `O(t log²(N+2))`, but the proof gives `O(t log²(t+2))` (no contradiction, since t ≤ N).
5. [Ambiguous] When building Alg 6's "scale tree per band" from the coarsest scale, the fact that one quotient component becomes exactly one part at the top scale follows from Lemma 6.7(1), but it is not stated in the procedure in the text. [Supplement] Top-level part = quotient component.
6. [Ambiguous] We interpret Alg 7 l.9 "direct offers from the root" as placing the `w'` of each radial edge incident to a in G' (originally an edge incident to a pivot) at the native inward copy of its other endpoint (J starts at (a, 0), l.1819). We read "offers from edges incident to the source" at the end of Sec 7.3 (l.1710) as referring to the root a in the context of J.
7. [Note] Lemma 6.7 uses `H ≥ (16t)^5`, but the proof needs only `H ≥ 4tT²`.
8. [Note] At the minimum N `589^64 + 1`, |A| = 1 is possible. In that case NameRoutes is empty, every edge is radial, and J becomes a full SSSP from the pivot, so the whole outer machinery is redundant work. The outer level is meaningful only when N/b is large enough.
9. [Note] The base case of Alg 8 l.5 specifies a "scan-based" Dijkstra (O(N²)). This choice is justified only at constant size, so when Mode A handles a real N, it is reasonable to use the binary-heap version (in the same Lemma 3.2). The output is the same. Record it as a deviation (`thresholds.md` §7).
