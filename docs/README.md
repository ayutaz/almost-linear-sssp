# Documentation guide

**English** | [日本語](README.ja.md) | [简体中文](README.zh-CN.md)

This page explains what lives under `docs/`, who each part is for, and the order to read it in.
For the project overview, see the [top-level README](../README.md).

Some documents are currently available only in Japanese. Where a document exists in both
languages, the Japanese version is the `.ja.md` file next to the English one. File names, section
numbers, commands and identifiers are the same in every language.

## Where to start

| If you want to… | Read, in this order |
|---|---|
| understand what the project found so far | [../README.md](../README.md) → [research/thresholds.md](research/thresholds.md) §2–§3 → [ADR 0003](adr/0003-three-implementation-modes.md) → [ADR 0008](adr/0008-mode-a-outer-base-case-variants.md) |
| see what is done and what is next | [PROGRESS.md](PROGRESS.md) → [TASKS.md](TASKS.md) |
| implement a part of the paper | [research/README.md](research/README.md) → the note for your section → [PAPER_MAPPING.md](PAPER_MAPPING.md) → [research/open-questions.md](research/open-questions.md) → [tdd/README.md](tdd/README.md) |
| work on graphs, baselines or verification | [design/graph_library.md](design/graph_library.md), [design/baselines.md](design/baselines.md), [design/numerics.md](design/numerics.md) |
| run or review benchmarks | [design/benchmark_protocol.md](design/benchmark_protocol.md) → [design/benchmark_schema.md](design/benchmark_schema.md) → [design/benchmark_runner_cli.md](design/benchmark_runner_cli.md) → [research/prior_work_and_oss.md](research/prior_work_and_oss.md) |
| change CI or tooling | [design/ci.md](design/ci.md), then the ADRs it refers to |
| contribute code | [../CONTRIBUTING.md](../CONTRIBUTING.md) → [tdd/README.md](tdd/README.md) → [TASKS.md](TASKS.md) |

## What each part is

### `research/` — reading the paper

Notes that explain the paper (arXiv 2611.34911v1) section by section, in our own words: the
definitions, lemmas, algorithms and constants, what they depend on, and what an implementation
has to watch out for. The paper itself is not stored in the repository; run `make fetch-paper`
to read it locally. Line numbers such as `l.284` refer to the text extracted by that script.

- **Who reads it:** anyone implementing or reviewing a paper component, and anyone who wants to
  check our claims about the paper.
- **Index and reading order:** [research/README.md](research/README.md).
- **Most important note:** [research/thresholds.md](research/thresholds.md) — the exact sizes at
  which the paper's constructions start, and which constants affect correctness versus only
  running time.
- **Prior work:** [research/prior_work_and_oss.md](research/prior_work_and_oss.md) — related
  algorithms and existing open-source implementations, with licences.

### `PAPER_MAPPING.md` — paper ↔ code

One row per definition, lemma, algorithm and key constant of the paper, with the module that
implements it, its tests, its mode (A or B) and its status. It also lists every **deviation**
from the paper (`D-xx`). Read it before touching `src/paper/` or
`src/experimental/paper_experimental/`, and update it with every change there.

### `design/` — how the code is built

Design documents for the parts that are not taken from the paper.

| Document | Topic |
|---|---|
| [graph_library.md](design/graph_library.md) | Graph representations, text formats, reproducible generators, `sssp_gen` |
| [baselines.md](design/baselines.md) | Dijkstra variants and heaps, tie-breaking, operation counters |
| [numerics.md](design/numerics.md) | Comparison–addition model versus IEEE-754 `double`; the tolerance used for `double` results and why |
| [benchmark_protocol.md](design/benchmark_protocol.md) | How benchmarks are run: builds, repetitions, timing boundaries, verification, statistics, suites |
| [benchmark_schema.md](design/benchmark_schema.md) | The format of benchmark result records |
| [benchmark_runner_cli.md](design/benchmark_runner_cli.md) | Command-line interface of the C++ benchmark runner |
| [ci.md](design/ci.md) | GitHub Actions workflows and the setup scripts |

- **Who reads it:** contributors changing that area, and reviewers.

### `adr/` — decisions

Architecture Decision Records: one short file per decision that constrains future work, with its
context and consequences. A decision is never edited away; it is replaced by a new ADR. The index
is [adr/README.md](adr/README.md). Start with
[0003](adr/0003-three-implementation-modes.md) (the three modes) and
[0008](adr/0008-mode-a-outer-base-case-variants.md) (the paper-faithful base case).

- **Who reads it:** anyone who wonders "why is it done this way?" before proposing a change.

### `TASKS.md` and `PROGRESS.md` — tracking

- [TASKS.md](TASKS.md) is the single list of work items: ID, priority, status, dependencies,
  acceptance criteria, test evidence and closing commit.
- [PROGRESS.md](PROGRESS.md) is a short snapshot: current phase, what is implemented, latest
  test and benchmark state, known issues and next tasks. Read it first when resuming work.

Both are maintained by the coordinating maintainer; other contributors propose changes rather
than editing them in parallel.

### `tdd/` and `TDD_LOG.md` — test-first evidence

Every feature or bug fix has a task file `tdd/tasks/<ID>.json` that records the failing test run
(RED), the passing run (GREEN) and the refactoring run, with commands, exit codes and output.
The format is described in [tdd/README.md](tdd/README.md). `TDD_LOG.md` is a log generated from
these files by `scripts/tdd.py log`; it is not edited by hand (it does not exist until it is
first generated).

- **Who reads it:** reviewers checking that a change was really test-first.

### `reports/` — results

Benchmark and review reports generated from saved runs. It is empty for now: **nothing has been
measured yet**. Every number that appears here must come from a run saved under `results/`.

## Conventions

- The paper is cited by section, lemma, algorithm and line number, never quoted at length.
- Measured and estimated numbers are never mixed; anything not measured is written as
  "not measured".
- No absolute local paths or personal information: this is a public repository. Run
  `make public-safety` before publishing changes.
