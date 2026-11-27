# Numerics: exact model vs. IEEE-754 doubles

The paper (arXiv 2611.34911v1, Sec 2) works in the **comparison-addition model**: edge lengths
are non-negative reals, and an algorithm may only compare and add them, exactly. A CPU `double`
rounds every addition. This document fixes how the project tells the two apart, what tolerance
a `double` result is checked with and why, and why integer weights are checked exactly.

Code: `include/sssp/verify.hpp` (`Tolerance`, `verify_sssp`, `write_result_json`),
`include/sssp/counted.hpp` (`Counted<W>`), `tests/python/reference/` (exact rational oracle).

## 1. Two weight types, two standards of correctness

| weight type     | arithmetic                      | correctness standard                          |
|-----------------|---------------------------------|-----------------------------------------------|
| `std::uint64_t` | exact (as in the model)         | **exact** equality, overflow is an error      |
| `double`        | IEEE-754 binary64, round-to-nearest | within the tolerance of Sec 3, or exact in the rounded algebra (Sec 2) |

Integers are the primary correctness target: on integer inputs every model operation is exact
as long as no sum reaches `WeightTraits<uint64_t>::infinity() = 2^64 - 1`, so an
implementation of the paper must reproduce the exact distances. `double` results are judged
with a rounding-error bound.

## 2. What relaxation algorithms compute in `double`

Write `fl(x)` for rounding to the nearest double and `a ⊕ w = fl(a + w)`. For `w ≥ 0`:

1. **monotone:** `a ≤ b ⇒ a ⊕ w ≤ b ⊕ w` (rounding is monotone);
2. **non-decreasing:** `a ⊕ w ≥ a`, but **not strictly** — absorption: `2^60 ⊕ 1 = 2^60`.

(1)+(2) make `(double, min, ⊕)` a monotone, non-decreasing path algebra, on which Dijkstra and
Bellman–Ford are correct. Both therefore compute exactly the same value

    d̂(v) = min over s–v paths P of the left-to-right rounded sum of P,

and the minimum is attained by a simple path (removing a cycle never increases a left fold, by
(1)+(2)). Consequences used by the tests:

* Any two relaxation-based algorithms (Bellman–Ford, every Dijkstra variant) must agree
  **bit for bit** on `dist`; only the parent tree may differ (ties).
* For such results the local checks of `verify_sssp` hold **exactly**:
  `dist[v] == fl(dist[parent] + w)` and `dist[v] ≤ fl(dist[u] + w)` on every arc. Hence
  `Tolerance::exact()` is the default even for `double`.
* With strict `<` relaxation, property (2) keeps the parent graph acyclic even under absorption.

An algorithm that adds in another association (e.g. precomputed shortcut lengths
`d(u) + (w1 + w2)`, Floyd–Warshall's `d[i][k] + d[k][j]`) computes a different rounded value,
so it can differ from the relaxation result by rounding. That is what the tolerance is for.

## 3. The tolerance

`Tolerance{rel, abs}` accepts two distance values `a`, `b` if

    |a − b| ≤ rel · max(|a|, |b|) + abs                  (one-sided: a ≤ b + same slack)

and `Tolerance::for_path_hops(k)` sets `rel = γ_k`, `abs = 0`, with

    u = 2^-53 (unit roundoff of binary64),     γ_k = k·u / (1 − k·u)     (requires k·u < 1/2).

`Tolerance::for_graph(g)` uses `k = n − 1`.

### 3.1 Why γ_k

Summing `m` non-negative doubles `x_1..x_m` with **any** parenthesisation (any binary tree of
`m − 1` additions) gives (Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed.,
Lemma 3.1 and eq. (4.4)):

    |ŝ − s| ≤ γ_{m−1} · Σ|x_i| = γ_{m−1} · s.

A shortest path is simple, so it has at most `n − 1` edges; with `dist[source] = 0` the first
addition `0 + w_1` is exact, so at most `n − 2` additions round. `γ_{n−1}` is therefore a safe
(one step conservative) bound for any evaluation order. Let `d*` be the exact distance and
`d̂ = ŝ(P)` the computed one, realised by some path `P` (the tree path):

* lower: `d̂ = ŝ(P) ≥ (1 − γ) L(P) ≥ (1 − γ) d*`;
* upper: `d̂ ≤ ŝ(P*) ≤ (1 + γ) L(P*) = (1 + γ) d*`  (the algorithm minimises over paths);

so **`|d̂ − d*| ≤ γ_{n−1} · d*`** — the bound `check_double.py` checks against the exact
rational `d*`.

### 3.2 Comparing two rounded results

If `a` and `b` are both rounded results (e.g. a paper implementation vs. double Floyd–Warshall),
each is within `γ_k d*`, so `|a − b| ≤ 2γ_k d* ≤ 2γ_k/(1 − γ_k) · max(a, b) = γ_{2k} · max(a, b)`
(the identity `2γ_k / (1 − γ_k) = γ_{2k}` is exact). Use `Tolerance::for_path_hops(2(n − 1))`.
Against the double nearest to `d*`, add half an ulp: `γ_k + u ≤ γ_{k+1}`.

### 3.3 Why no absolute term by default

IEEE-754 addition never loses accuracy to underflow: if the sum of two doubles is subnormal it
is exact. So for non-negative sums the error is purely relative and `abs = 0`. `abs` exists for
callers that knowingly leave the comparison-addition model (e.g. bucket keys computed by
subtraction, `d − base`) or want explicit slack; it must be justified where it is used.

### 3.4 What the local checks can and cannot conclude under a tolerance

Tree check `dist[v] ≈ fl(dist[p] + w)` and certificate `dist[v] ≲ fl(dist[u] + w)` use the same
slack per arc. With slack, the certificate argument (`dist[v] ≤ d*(v)`) accumulates one slack
per hop, so a tolerant local check alone proves only `dist ≤ d*(1 + γ)^{hops}`-type bounds. Global
accuracy of a non-relaxation algorithm is established by comparison with a reference: exact
rational (`check_double.py`), or a relaxation result with `for_path_hops(2(n − 1))`.

`check_double.py` also bounds the exact length `L(v)` of the reported tree path: since
`dist(v)` is a rounded sum of that path, `(1 − rel) L ≤ dist(v) ≤ (1 + rel) d* + abs`, hence
`d* ≤ L ≤ ((1 + rel) d* + abs) / (1 − rel)`.

## 4. Where the model and doubles really differ

Reasoning valid for reals that fails for doubles, and must not be relied on in `double` code:

* **strict growth:** `w > 0 ⇒ a + w > a` (false under absorption; see `extreme_ratio` in
  `tests/support/graphs.hpp`, weights `2^60` and `1`);
* **associativity:** `(a + b) + c = a + (b + c)` (false; makes results order-dependent);
* **cancellation:** `(a + w) − w = a` (no subtraction exists in the model; `Counted<W>` forbids it).

Algorithms whose termination or bucket arguments use strict growth must be tested on
`extreme_ratio`-style inputs with `double`, and are primarily validated on `uint64_t`.

## 5. Why fast-math is forbidden

`cmake/Warnings.cmake` passes `-fno-fast-math` to every target. `-ffast-math` would enable:

* `-ffinite-math-only`: the compiler may assume no infinities, so `dist[v] == inf` (our
  "unreachable" sentinel) may be folded to `false` — silently wrong results;
* `-fassociative-math` / reassociation: destroys the rounded path algebra of Sec 2, so
  relaxation algorithms no longer agree bit for bit and results depend on optimisation level;
* flush-to-zero / denormals-are-zero: subnormal sums are no longer exact, so Sec 3.3 fails;
* `-fno-signed-zeros`, reciprocal math: harmless here, but no benefit.

Floating-point contraction (FMA) is irrelevant because distances use only additions.
`FLT_EVAL_METHOD == 0` (SSE2 / AArch64, no x87 extended precision) is assumed.

## 6. Integer weights: exact, with overflow detection

* `infinity() = 2^64 − 1` is reserved; every finite distance is `< 2^64 − 1`.
* `WeightTraits<uint64_t>::checked_add` throws `std::overflow_error` if `a + b ≥ infinity()`.
* `bellman_ford` skips an overflowing candidate (it may lie on a non-shortest path) and throws
  `std::overflow_error` only if, after convergence, a vertex adjacent to a finite one is still
  infinite — i.e. iff some true distance is not representable.
* `verify_sssp` compares integers exactly and ignores `Tolerance`; an overflowing
  `dist[parent] + w` is reported as `arithmetic_overflow`; on the certificate an overflowing
  bound is vacuous (the true constraint lies beyond the range).

## 7. Exact verification of a C++ `double` result

1. Write the graph as the edge list `n m` / `u v w` with weights printed `%.17g` (so
   `Fraction(float(w))` recovers the exact double the C++ code used).
2. Dump the result with `sssp::write_result_json(os, source, result)`:
   `{"n":..,"source":..,"dist":[%.17g or null],"parent":[-1..],"parent_edge":[-1..]}`.
3. `uv run python tests/python/reference/check_double.py graph.txt result.json`
   (`--hops K`, `--rel R`, `--abs A`, `--weights exact` to use the decimal weights as reals).
   It recomputes `d*` with `fractions.Fraction`, reads each dumped number as the exact binary
   value of its double, and checks Sec 3.1 plus the tree (exit 0 / 1 / 2 = ok / rejected /
   unreadable).

The dump format is pinned on both sides: `tests/unit/verify/test_verify.cpp` asserts the exact
string, and `tests/python/test_exact_reference.py::test_cpp_dump_format_is_readable_and_exact`
parses the same string.

## 8. Counting model operations

`sssp::Counted<W>` (W = `double` or `uint64_t`) counts `+`/`+=` and the six comparisons in
thread-local counters (`counted_stats()`, `reset_counted_stats()`, `CountedScope`). It has no
`-`, `*`, `/` or conversion back to `W`, so code that compiles with `Counted<W>` provably uses
only comparison-addition operations on weights. Overflow guards in `checked_add` are not
counted (bookkeeping outside the model). Limitation: `validate_edges` only validates signed
or floating-point `W`, so a negative or NaN `Counted<double>` weight is not rejected at graph
construction; build the graph from validated raw weights first.
