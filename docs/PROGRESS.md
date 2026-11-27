# Progress

Read this first when resuming work (then `docs/TASKS.md`, then `uv run python scripts/tdd.py status`).

## Current phase
Phase 0 (setup/harness) and Phase 1 (paper research) running in parallel with the Phase 2
foundations (graph library, baselines, verification, benchmark infrastructure).

## Implemented
- Build skeleton: CMake presets (debug, release, release-diag, sanitizer, coverage, ci,
  ci-release), Catch2 v3.8.1 and Google Benchmark v1.9.4 pinned by SHA-256, Makefile with
  `make verify`.
- Shared interfaces: `include/sssp/graph.hpp` (CSR, weight traits, validation),
  `include/sssp/result.hpp` (`SsspResult`, `Diagnostics`, `FallbackRecord`),
  `include/sssp/baseline.hpp` (declarations only).
- Tasks DONE: P0-VERSION, P0-GRAPH-CSR (commit ce4c2ac).

## In progress
P0-TDD-TOOL, P0-HOOKS, P0-CI, P0-RULES-DOCS, P1-PAPER-SEC1-3, P1-PAPER-SEC4-5,
P1-PAPER-SEC6-8, P1-OSS-SURVEY, P2-GRAPH-*, P2-DIJKSTRA-*, P2-VERIFY-*, P2-BENCH-PY.

## Test results
- Last recorded: `ctest --preset debug` — 100% passed (6 tests) at P0-GRAPH-CSR GREEN
  (2026-11-27T10:15+0900). `make verify` not yet run on the full tree.

## Latest measurements
Not measured. No benchmark has been run yet.

## CI
Not run yet (workflow in progress, P0-CI).

## Known issues
- Mode A outer threshold `N > L_N^64` cannot hold for any `N < 2^64`; Mode A will always take
  the recorded outer base case on feasible inputs. Mode B (P5-EXP-THRESHOLDS) exists to study
  the outer construction.
- `docs/PAPER_MAPPING.md` and `docs/research/open-questions.md` not yet created.
- `scripts/fetch_paper.sh` and `scripts/bootstrap.sh` referenced but not yet implemented.

## Next tasks
1. Finish P0-TDD-TOOL and P0-HOOKS, then set `.tdd/config.json` enforce=true.
2. P1-PAPER-MAPPING once the three section analyses are done; refine Phase 3 rows.
3. P2-DIFF-HARNESS and P2-BENCH-CPP once generators, oracles and binary Dijkstra are DONE.

Last updated: 2026-11-27T10:20+0900

<!-- AUTO:precompact -->
<!-- /AUTO -->
