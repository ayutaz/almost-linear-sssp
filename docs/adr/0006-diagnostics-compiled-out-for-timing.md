# 0006. Diagnostics compiled out in timed builds
- Status: Accepted
- Date: 2026-11-27

## Context
Execution-path counters and operation counts (`Diagnostics` in `include/sssp/result.hpp`) are
needed to verify which paper constructions ran and how often fallbacks occur, but they perturb
timing.

## Decision
Counters are updated through `SSSP_COUNT(...)`, controlled by the `SSSP_DIAGNOSTICS` compile
definition. The `release` preset sets `SSSP_DIAGNOSTICS=OFF` and is the only preset used for
timing. Counts come from a separate `release-diag` run on the same inputs. `FallbackRecord`s
are always recorded (they are part of the result, not a counter).

## Consequences
- Timed runs carry no counter overhead.
- Reports combine timing (release) and counts (release-diag) from the same commit and inputs.
