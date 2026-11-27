---
name: correctness-review
description: >-
  Independently verify an SSSP implementation or paper component against reference oracles, shortest-path-tree checks, boundary cases and random graphs; minimise any counterexample and add it as a regression test. Use after GREEN on any algorithmic task, before marking it DONE, and when a result looks suspicious.
when_to_use: >-
  "review correctness", "verify the implementation", "find a counterexample", "is this right", "differential test", before closing any P2-/P3-/P5- algorithm task.
---
# Correctness review

Goal: find a counterexample, not confirm the implementation. Work from the spec / paper, not
from the implementer's explanation.

## 1. Oracles (independent of the code under test)
- `sssp::baseline::bellman_ford` (O(nm), small graphs), binary-heap Dijkstra, exact rational
  reference (Python, `fractions.Fraction`) for `double` inputs.
- Shortest-path-tree verifier for every result:
  `dist[s]==0`; for every arc `dist[v] <= dist[u] + w`; for every non-source reachable `v`,
  `parent_edge[v]` joins `parent[v]` and `v` and is tight (`dist[parent]+w == dist[v]`);
  parent pointers are acyclic and reach `s`; unreachable ⇔ `dist == infinity()` and
  `parent == kNoVertex`.
- For paper components, check the lemma's invariant directly (e.g. merge order, LCA answers
  against a naive walk, interval minima against a linear scan).

## 2. Inputs to cover
- Boundary: n=0 (if allowed), n=1, isolated source, disconnected graph, self-loops,
  multi-edges with different weights, zero-weight edges and zero-weight cycles, all-equal
  weights (ties), path/star/complete graphs, very large weights (uint64 overflow detection:
  expect `std::overflow_error`), subnormal and huge `double`s, invalid input (negative, NaN,
  inf, out-of-range endpoints, invalid source → `std::out_of_range`).
- Threshold edges: N exactly at `2*L_N^4` and `2*L_N^4 + 1`; L_N jumps (N+2 a power of two).
- Random: Erdős–Rényi, random geometric, grid, power-law, trees + few extra edges; integer
  weights in small ranges (many ties) and wide ranges; several seeds per family.
- Randomised algorithms: many seeds; the output must be exact for every terminating run.

## 3. Floating point
- `uint64_t`: exact equality, no tolerance.
- `double`: compare against the exact rational reference. Allowed error per vertex is
  `|d - d*| <= k * eps * d*` where `k` = hop count of the reference path (documented in the
  test); report the max observed error. Never widen a tolerance to make a test pass.
- Any build used here must not use fast-math.

## 4. Commands
```bash
make test-differential test-property test-regression
make sanitizer                      # ASan + UBSan
uv run pytest tests/python -k differential
```
Seeded fuzz loop (property tests use Catch2 generators / Hypothesis): raise the iteration count
locally, keep CI counts small.

## 5. When a counterexample is found
1. Save the failing input (edge list, source, seed) immediately.
2. Minimise: delete edges/vertices while it still fails (delta debugging); reduce weights to
   small integers if the failure persists.
3. Add it to `tests/regression/` as a named test case (`regression: <task-id> <short cause>`),
   run it to record RED via `/tdd`, then hand the fix to the implementer.
4. Report: input, expected vs actual, which oracle disagreed, suspected paper step.

## 6. Verdict
Report one of: PASS (with the list of checks and seeds run), FAIL (counterexample attached),
INCONCLUSIVE (what could not be checked and why). Never PASS without having run the checks.
