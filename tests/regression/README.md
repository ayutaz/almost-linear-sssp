# Regression cases

Every counterexample found by differential, property or random testing is **minimised** and
kept here forever, so the bug can never come back silently.

A case has two parts with the same id `<YYYYMMDD>-<short-name>` (e.g. `20261127-zero-cycle-parent`):

1. `tests/regression/cases/<id>.json` — the machine-readable record (format below).
2. `tests/regression/test_<short_name>.cpp` — a Catch2 `TEST_CASE("<id>: ...", "[regression]")`
   that inlines the graph, runs the algorithm and asserts `sssp::verify_sssp(...)` (with the
   exact reference distances) is `ok`. It must fail on the buggy version (RED) and pass after
   the fix; record both in the task's `docs/tdd/tasks/<ID>.json`.

## Record format (`cases/<id>.json`)

```json
{
  "id": "20261127-zero-cycle-parent",
  "summary": "parent pointers form a cycle on a zero-weight triangle",
  "algorithm": {
    "name": "sssp::paper::sssp",
    "version": "git 1a2b3c4 (or library version 0.1.0)",
    "params": {"variant": "...", "thresholds": {"L_N": 4}}
  },
  "weight_type": "uint64",
  "graph": {
    "n": 3,
    "edges": [[0, 1, "0"], [1, 2, "0"], [2, 0, "0"]]
  },
  "source": 0,
  "found_by": {
    "test": "differential:random multigraphs",
    "generator": "sssp::test::random_edges",
    "seed": 123,
    "params": {"n": 25, "m": 60, "lo": 0, "hi": 3}
  },
  "minimised": true,
  "failure": {"kind": "parent_cycle", "message": "vertex 1 lies on a parent cycle ..."},
  "expected_dist": ["0", "0", "0"],
  "tolerance": null,
  "fixed_by": "fix(paper): ... (commit sha)"
}
```

Field rules:

- `weight_type`: `"uint64"` or `"double"`.
- Edge weights and `expected_dist` are **strings**, so no JSON parser rounds them: integers
  for `uint64`; for `double` the `%.17g` text of the exact double used (round-trips).
  `expected_dist` entries are `"inf"` for unreachable vertices; compute them with
  `tests/python/reference/exact_sssp.py`, never from the implementation under test.
- Large graphs: instead of `graph`, give `"graph_file": "cases/<id>.txt"` in the edge-list
  format of `exact_sssp.py` (`n m` then `u v w` lines), after minimisation it should rarely
  be needed.
- `found_by.seed` and `found_by.params` must reproduce the original (unminimised) failure.
- `failure.kind` is a `sssp::to_string(FailureKind)` name (or `"exception"`, `"timeout"`).
- `tolerance`: `null` for exact checks, else `{"rel": "...", "abs": "..."}` with the
  justification in `summary` (see `docs/design/numerics.md`).

## Minimising

Delete edges, then vertices (renumbering), then lower weights, while the failure persists;
prefer the fewest vertices, then the fewest edges, then the smallest weights.
