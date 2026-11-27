# 0004. Pinned dependencies
- Status: Accepted
- Date: 2026-11-27

## Context
Reproducibility requires that tests and benchmarks run against the same library versions on
every machine and in CI.

## Decision
- C++: dependencies are fetched with CMake `FetchContent` by release tag **and** SHA-256
  (`cmake/Dependencies.cmake`): Catch2 v3.8.1, Google Benchmark v1.9.4.
- Python: managed by uv; dependencies are added with `uv add` (`uv add --dev` for tooling) and
  locked in `uv.lock`, which is committed. pip is not used.
- Bumping a pinned version requires a new ADR (or an amendment ADR) stating the reason.

## Consequences
- Builds are reproducible offline once the FetchContent cache exists.
- Version bumps are deliberate and reviewed.
