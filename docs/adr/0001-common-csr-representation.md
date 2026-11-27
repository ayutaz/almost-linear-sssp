# 0001. Common CSR graph representation for all algorithms
- Status: Accepted
- Date: 2026-11-27

## Context
Benchmarks compare baselines, the paper-faithful algorithm and experimental variants. If each
used its own graph layout, timing differences would mix algorithmic effects with memory-layout
effects. The paper's input is an undirected multigraph with finite non-negative lengths (Sec 2).

## Decision
`sssp::CsrGraph<W>` (`include/sssp/graph.hpp`) is the single input type for every algorithm.
Each undirected edge is stored once in the edge list and as two arcs (one for a self-loop) in
the CSR; arcs carry the original edge id. Multi-edges, self-loops and isolated vertices are
kept. Weights are `std::uint64_t` (exact, overflow-checked) or `double` (IEEE-754). Invalid
input (out-of-range endpoints, negative/NaN/inf weights) is rejected at construction.
Algorithms needing other structures build them from the CSR inside the timed region and report
that cost.

## Consequences
- Fair comparisons on identical memory layouts and inputs.
- Results (`SsspResult`) can report parent edge ids for tree verification.
- Algorithms that want a different layout pay the conversion inside their own timing.
