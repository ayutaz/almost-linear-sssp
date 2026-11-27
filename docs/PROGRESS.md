# Progress

Read this first when resuming work (then `docs/TASKS.md`, then `uv run python scripts/tdd.py status`).

## Current phase
Phase 3 (paper components, Mode A) in progress. Phase 2 foundations are on `main` except the
items listed under "In progress"; Phase 0 harness (hooks, `scripts/tdd.py`) is not yet committed.

## Implemented (on `main`, 9a58a44)
- Build skeleton: CMake presets (debug, release, release-diag, sanitizer, coverage, ci,
  ci-release), Catch2 v3.8.1 and Google Benchmark v1.9.4 pinned by SHA-256, Makefile with
  `make verify`. Shared interfaces `include/sssp/graph.hpp` and `include/sssp/result.hpp`.
- Scripts: `scripts/fetch_paper.sh` (pinned SHA-256), `scripts/bootstrap.sh` (exists, no TDD
  evidence yet), `scripts/check_public_safety.py`, `scripts/coverage.sh`.
- Baselines (`src/baseline/`): Dijkstra with binary, 4-ary, Fibonacci and radix heaps, and the
  four heaps; Bellman–Ford oracle.
- Verification: `verify_sssp` shortest-path-tree verifier (`src/verify/`), exact rational
  reference in Python (`tests/python/reference/`), `Counted<W>`, exhaustive small-graph
  differential test, differential test of every baseline on all generator families, `sssp_run`
  CLI with a cross-language check against the exact oracle.
- Graph library: `EdgeList`/`AdjacencyList` formats, DIMACS and edge-list I/O, reproducible PRNG,
  weight distributions, generators, small boundary-case suite, `sssp_gen` CLI.
- Benchmark infrastructure, Python side: record schema, environment capture, statistics, job
  plan, `scripts/run_benchmarks.sh`, `scripts/analyze_results.py`. The C++ runner is not on
  `main` yet (P2-BENCH-RUNNER).
- CI workflows in `.github/workflows/` (ci, sanitizer, benchmark, release).
- Paper, Sec 2 (`src/paper/model/`): thresholds, length type, integer log, randomness.
- Paper, Sec 3 (`src/paper/primitives/`): searches, stable sorts, BFPRT, radius Dijkstra
  (Lemma 3.1–3.2); Fibonacci heap, GrowthPass, three-pass MST (Lemma 3.3–3.5); DSU, LCA, level
  ancestors, forest queries (Lemma 3.6–3.9); split-findmin (Lemma 3.10–3.11).
- Paper, Sec 4 (`src/paper/hierarchy/`): ordered hierarchy types and H1–H6 validator
  (Def 4.1/4.2), threshold-SCC reference hierarchy (test support).
- Paper, Alg 4 base case (`src/paper/local/base_case`); Alg 8 l.2 zero-edge contraction and
  Lemma 8.5 predecessor-tree restoration (`src/paper/outer/contract_restore`).
- Research notes and mapping: `docs/PAPER_MAPPING.md` and all notes in `docs/research/`
  (Phase 1 DONE). ADRs 0001–0008.

## In progress
- Phase 0: P0-HARNESS-TDD-CLI and the six P0-HOOKS sub-tasks (not committed), P0-CI,
  P0-RULES-DOCS.
- Phase 2: P2-BENCH-RUNNER, P2-DIFF-HARNESS (open regressions P2-DIFF-DOUBLE-OVERFLOW,
  P2-DIFF-REGULAR-SWITCH), P2-DIJKSTRA-BINARY/4ARY/FIB (blocked on the double-overflow fix),
  P2-GRAPH-GEN.
- Phase 3: P3-NORMALIZE, P3-CHILD-ORDER, P3-EVAL-ORACLE (verification infrastructure for
  P3-EVALUATE), P3-OWNER-DECOMP, P3-RADIAL-STATES, P3-THRESHOLD-JOINS, P3-LOCAL-HIERARCHY
  (Alg 4 entry), P3-ESTIMATE, P3-NAME-ROUTES, P3-OUTER-BASE and P3-PAPER-SSSP (Mode A
  end-to-end path), and the remaining items of P3-INTERVAL-MIN and P3-HIERARCHY-TYPES (see
  TASKS.md).

## Test results
- `ctest --preset debug` on `main` at ced30ef: 100% passed, 556 of 556 (2026-11-27T11:00+0900).
  9a58a44 changes documentation only.
  Three of these are known-bug regression tests tagged `[!shouldfail]`; they count as passed
  because they fail as expected (see Known issues).
- pytest and `make verify` were not run for this update.

## Latest measurements
Not measured. No benchmark has been run yet.

## CI
Workflows exist but have never run: no git remote is configured yet.

## Known issues
- Mode A outer threshold `N > L_N^64` is false for every `N < 2^64`; the smallest N that
  satisfies it is 589^64 + 1 ≈ 1.94×10^177 (exact values in `docs/research/thresholds.md`).
  Mode A always takes the recorded outer base case on feasible inputs; Mode B
  (P5-EXP-THRESHOLDS) exists to study the outer construction.
- Known bug, fix in progress (P2-DIFF-DOUBLE-OVERFLOW): on double weights the binary, 4-ary
  and Fibonacci Dijkstra baselines return infinity when a path sum overflows, i.e. for a
  reachable vertex whose true distance exceeds DBL_MAX, instead of throwing
  `std::overflow_error` (Bellman–Ford already throws). The regression test
  `tests/regression/test_double_sum_overflow.cpp` is temporarily tagged `[!shouldfail]`.
- Known bug, fix in progress (P2-DIFF-REGULAR-SWITCH): `random_regular(n, d, kSwitchToSimple)`
  does not converge for small, dense (n, d), e.g. (3, 2) with seed 20 never leaves a
  three-self-loop state. The regression test
  `tests/regression/test_regular_switch_dead_end.cpp` is temporarily tagged `[!shouldfail]`.
- Evidence files P2-VERIFY-SPT, P2-EXACT-REFERENCE, P2-COUNTED-WEIGHT and P2-EXHAUSTIVE-DIFF
  record pre-rebase SHAs; TASKS.md gives the SHAs on `main`.
- Evidence of the heaps, Dijkstra baselines and P2-BENCH-* ends at REFACTOR; `tdd.py done` has
  not been recorded. Most other DONE evidence files have an empty `commit` field.
- P2-GRAPH-GEN: no random geometric generator yet (planned in the acceptance criteria).
- P3-EDGE-PAIRS and P3-MERGE-EVENTS are TODO; the MST and the Alg 4 base case contain their
  own lightest-edge-per-pair and equal-weight merge logic for now.

## Next tasks
1. Finish the fixes for P2-DIFF-DOUBLE-OVERFLOW and P2-DIFF-REGULAR-SWITCH (remove
   `[!shouldfail]` with the fix), then close P2-DIJKSTRA-BINARY/4ARY/FIB and P2-GRAPH-GEN.
2. Commit P0-HARNESS-TDD-CLI and P0-HOOKS, then set `.tdd/config.json` enforce=true.
3. Integrate P2-BENCH-RUNNER; decide whether P2-BENCH-CPP is still needed.
4. Phase 3: integrate P3-OUTER-BASE and P3-PAPER-SSSP (the Mode A path that runs on feasible
   inputs); P3-EDGE-PAIRS, P3-MERGE-EVENTS; P3-NORMALIZE → P3-CHILD-ORDER → P3-EVALUATE (with
   P3-EVAL-ORACLE).
5. Configure the remote and get a first green CI run (P0-CI).

Last updated: 2026-11-27T11:15+0900

<!-- AUTO:precompact -->
<!-- /AUTO -->
