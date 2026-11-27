# almost-linear-sssp — rules for Claude Code

## Purpose
Reproduce, verify, measure and optimise **"Linear Single-Source Shortest Path, Almost"**
(Gholami, Hajiaghayi), **arXiv:2611.34911 v1** — https://arxiv.org/abs/2611.34911 .
Scope: SSSP on undirected graphs with finite non-negative real weights; claimed expected time
O(m + n log log log n). We target **v1 only**; if a later version appears, record it in an ADR
before following it.

The paper is NOT redistributable (arXiv non-exclusive licence). It lives in `.cache/paper/`
(git-ignored) and is fetched with `scripts/fetch_paper.sh`. Never commit or quote long passages
of it; cite by section / lemma / algorithm line instead.

## Stack
C++20 (Apple Clang / LLVM Clang), CMake + Ninja via `CMakePresets.json`, Catch2 v3,
Google Benchmark, Python 3.12+ managed by **uv** (`uv add <pkg>`, `uv run <cmd>`; never pip),
pytest, Ruff, clang-format, clang-tidy, GitHub Actions.
Shared types: `include/sssp/graph.hpp` (CSR), `include/sssp/result.hpp` (`SsspResult`,
`Diagnostics`, `FallbackRecord`).

## Layout and implementation modes
| Path | Mode | Rule |
|---|---|---|
| `src/baseline/` | baselines | Dijkstra variants, Bellman–Ford oracle. Never "improved" to look like the paper. |
| `src/paper/` | **A: paper-faithful** | Only what the paper specifies. Thresholds exactly as written (below). |
| `src/experimental/paper_experimental/` | **B: paper-experimental** | Paper structure with tuned thresholds. |
| `src/experimental/practical/` | **C: practical-optimized** | Anything goes, as long as it is correct and labelled. |

Paper thresholds (Sec 2): `L_N = ceil(log2(N+2))`; outer construction requires `N > L_N^64`,
local construction requires `N > 2*L_N^4`. Evaluate them with overflow-safe arithmetic
(saturating power, early exit once the product exceeds N). Never use floating-point `log`/`pow`.

### Fallback policy
Whenever Mode A/B handles a region by direct search instead of the paper's construction (base
case, size threshold, unsupported input), push a `FallbackRecord{reason, vertices, edges}` into
`Diagnostics.fallbacks`. Silent fallbacks are bugs.
- Mode A: entering the paper's base case is legal but MUST be recorded.
- Mode B: an input whose correctness cannot be established is a **failure**, not a success.
  Substituting a reference algorithm and reporting success is forbidden.
- Reports must state the fallback share (vertices/edges handled by fallback) next to every timing.

## TDD is mandatory (see `/tdd` skill)
- Every feature or bug fix has a task ID from `docs/TASKS.md` and evidence in
  `docs/tdd/tasks/<ID>.json` (format: `docs/tdd/README.md`).
- Order: `tdd.py start` → write failing test → `tdd.py red` (must fail for the intended reason)
  → minimal implementation → `tdd.py green` → `tdd.py refactor` → commit → `tdd.py done`.
- With `.tdd/config.json` `enforce=true`, hooks block edits to `src/`, `include/`, `tools/`
  until a RED event exists for the active task. Do not bypass or disable the hooks.

## Correctness rules (see `/correctness-review` skill)
- Every SSSP implementation is checked against an **independent** oracle (Bellman–Ford, exact
  rational reference) and by shortest-path-tree verification (`dist[s]=0`, every arc satisfies
  `dist[v] <= dist[u]+w`, every tree arc is tight, parents reach `s`, unreachable = infinity).
- Integer weights: exact equality. `double`: compare to the exact reference with a documented
  tolerance; never loosen a tolerance to make a test pass.
- Every found counterexample is minimised and added to `tests/regression/` before the fix.

## Benchmark fairness (see `/benchmark` skill)
- Release preset (`release`, diagnostics OFF) only; same `CsrGraph`, same inputs, same machine,
  same warm-up and repetitions for every algorithm. Record environment metadata with results.
- Report median, mean, p95, stddev and n. Never mix measured and estimated numbers; anything
  not measured is written as "not measured".

## Git
- Conventional Commits (`feat:`, `fix:`, `test:`, `refactor:`, `perf:`, `docs:`, `build:`,
  `ci:`, `chore:`), one task per commit where practical. See `/git-workflow` skill.
- Never force-push, never rewrite pushed history, never add/change/remove remotes.
- Never hide failures: no `|| true`, skipped tests, `-DSSSP_WERROR=OFF` in CI, or deleted
  assertions to get green.
- Commit or push only when the user (or team lead acting on the user's request) asks.

## Quality gate
`make verify` (format-check, ruff, CI build with `-Werror`, ctest, clang-tidy, pytest) must pass
before every commit. `make sanitizer` before merging algorithmic changes.
`uv run python scripts/tdd.py check` must pass (every task file is consistent).

## Resuming work
1. Read `docs/PROGRESS.md` (current phase, known issues, next tasks).
2. Read `docs/TASKS.md` (status, dependencies, acceptance criteria).
3. Run `uv run python scripts/tdd.py status` (active task and phase).
4. `git status` / `git log --oneline -10`; then continue the active task's TDD cycle.

## Documentation rules
- Update `docs/TASKS.md` (Status, Test Evidence, Commit SHA) and `docs/PROGRESS.md` when a task
  changes state. Regenerate `docs/TDD_LOG.md` with `tdd.py log`; never edit it by hand.
- Paper-derived code: update `docs/PAPER_MAPPING.md` (see `/paper-reproduction` skill).
  Deviations from the paper go there and, if structural, into an ADR (`docs/adr/`).
- Decisions that constrain future work get an ADR. Superseding requires a new ADR.

## Parallel work
- Never let two agents edit the same file at the same time. Split work by file ownership at
  launch time; shared files (`docs/TASKS.md`, `docs/PROGRESS.md`, `CMakeLists.txt`) are
  edited only by the coordinating (lead) agent.
- Per-task evidence files and globbed CMake sources exist so that modules can be added in
  parallel without touching shared files.

## Forbidden
- Implementations that only make the current tests pass (hard-coded outputs, special-casing
  test inputs).
- Hiding an unimplemented paper component behind Dijkstra or another existing algorithm.
- Reporting performance that was not measured, or measured under different conditions.
- Claiming something is "the paper's algorithm" when it contains steps the paper does not have.
- Changing a test's expected value to match the implementation's output.
- Enabling `-ffast-math` (or `-Ofast`, `-funsafe-math-optimizations`) in any build used for
  correctness checks. `cmake/Warnings.cmake` forces `-fno-fast-math`; keep it.
- Committing the paper PDF/text, datasets under `datasets/raw/`, or raw benchmark output.

## Definition of Done (per task)
1. RED→GREEN(→REFACTOR) evidence recorded; `tdd.py check` passes.
2. Acceptance criteria in `docs/TASKS.md` met, including differential + tree verification
   for any SSSP implementation.
3. `make verify` passes locally; CI is green for the pushed commit (`gh run watch`).
4. `docs/PAPER_MAPPING.md`, `docs/TASKS.md`, `docs/PROGRESS.md` updated; deviations recorded.
5. Committed with a Conventional Commit message; `tdd.py done <ID> --commit <SHA>` recorded.
