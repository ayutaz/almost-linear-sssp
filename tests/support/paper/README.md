# Test support for Sec 4 (Evaluate)

Test-only headers for checking Algorithm 3 (`Evaluate`) and Algorithm 2 (`Scan`) of
arXiv 2611.34911 v1 against oracles that share no code with the implementation. Nothing in this
directory may be included from `src/`.

| Header | Contents |
|---|---|
| `eval_oracle.hpp` | `target_distances` (delta of Def 4.2, multi-source Dijkstra on the physical graph), `protected_states` (tight-arc BFS from tight initial offers), `state_graph_distances` (multi-source Dijkstra on the directed state graph), `compute_eval_reference` (all of them, rejects infeasible offers), `physical_minima`, `served_targets` |
| `eval_invariants.hpp` | `EvalObserver<W>` (hook contract), `EvalObserverFor` (same contract as a concept), `ForwardingObserver<Hooks, W>` (adapter from a src-side virtual hook type), `RecordingObserver`, `replay`, `state_dijkstra_trace`, `InvariantChecker<W>` (I1-I8 of research notes 02 sec 2.8 plus protocol checks) |
| `seeded_instances.hpp` | `reference_instance`, `zero_arc_instance`, `base_case_instance`, `reference_catalogue`, `base_case_catalogue` (fixed SplitMix64 seeds) |

## What Evaluate must produce (Lemma 4.6)

- Every protected state `x` emits exactly `delta(pi(x))` = `ref.state_target[x]`; for these
  states this also equals `ref.state_distance[x]` (dist_D from the initial offers).
- Every emission at any state is `>= delta(pi(x))`.
- Unprotected states may emit any value `>= delta` (it can be below dist_D because supplemental
  offers follow physical gap walks) or not emit at all (Q-29, research notes 02 A9). Do not
  compare them with dist_D.
- Per physical vertex `v` with `served_targets(...)[v]`, the minimum emission equals
  `ref.target[v]`.

## Hook contract for the Evaluate implementer

The full contract is in the header comment of `eval_invariants.hpp`. In short, `EvalOptions`
(research notes 02 sec 8: `{bool supplemental_offers; Checker* checks;}`) holds a pointer to a
src-side hook type with these member functions (plain argument types only):

```cpp
on_initial_offer(StateId x, W v)                       // Alg 3 l.3, before anything else
on_relax_offer(ArcId arc, W v)                         // every attempted decrease, arc of the INPUT D
on_supplemental_offer(StateId x, W v)                  // Alg 2 l.6
on_open(uint32 node, span<const StateId> states, W D0, W origin, W nu, W S)   // internal nodes
on_scan_begin(uint32 node, span<const StateId> states, W a', W b', W a, W b, W c)
on_scan_end(uint32 node, span<const StateId> states)
on_emit(StateId x, W v)                                // before relaxing the outgoing arcs
on_retire(uint32 node, span<const StateId> states)     // covers every state below the node
on_finish()
```

`states` is the node's leaf-layout interval; a normalised node is identified by its leaf set,
which equals that of the topmost input node it came from (Lemma 4.3 contraction), so no NodeIds
are shared. In a test, bind the src-side type to the checker with
`ForwardingObserver<SrcHooks, W> adapter(checker);` (it fails to compile if a hook is missing or
has a different signature), or call an `InvariantChecker` directly if `evaluate` takes its
observer as a template parameter (`EvalObserverFor`).

Interpretations built into the checker: the origin rule and the padding of supplemental offers
use the parent's NOMINAL slot end `b` (Q-28, A6/A7); D0 must be the least key below the node at
opening; emissions inside `Scan(X, I, c)` lie in `I` and are `<= c`; I3 is checked in its first
part only (new offers `>= a'` of the innermost portion), the rank part needs data the hooks do
not carry.

## Differential tests (`tests/differential/paper/test_evaluate_vs_dijkstra.cpp`)

Written with P3-EVALUATE (merged in 3ef315c); the table below is the original plan, and the file
also adds invariant sweeps and `double` checks beyond it. Each test runs `evaluate` with an `InvariantChecker`
attached, requires `report().ok()`, and then compares emissions with the oracle as below.
Integer weights (`std::uint64_t`) first; `double` only with integer-valued lengths.

| Test | Generator | Oracle / assertions |
|---|---|---|
| `BaseCaseSingleSource` | `base_case_catalogue(seed)`, the `*/single` instances, several seeds | `compute_eval_reference`; every native state emits `ref.target[v]`; `physical_minima == ref.target` at every vertex (all reachable vertices are served); additionally `verify_sssp`-style check against `sssp::baseline` Dijkstra on `inst.physical` |
| `MultiSource` | `base_case_catalogue(seed)` `*/multi` instances and `reference_catalogue(seed0, count)` instances with `seeded.sources.size() > 1` | protected states emit `state_target`; `physical_minima == ref.target` where `served_targets`; every emission `>= state_target`; compare `target_distances` with min_i c_i + Floyd-Warshall on small ones |
| `RandomDirectedStateGraph` | `reference_catalogue(seed0, count)` (threshold-SCC hierarchies, copies, shortcuts in every third instance, tie-heavy and spread lengths) | protected emissions `== ref.state_target`, all emissions `>= ref.state_target`, protected count `== ref.num_protected()` emitted; I1-I8 via the checker |
| `SupplementalOffersOffEquivalent` | same instances as `RandomDirectedStateGraph` and `BaseCaseSingleSource` | run twice, `supplemental_offers = true / false` (checker with the same flag, which also asserts that no supplemental offer is posted when disabled); the protected emissions are identical and equal `state_target` (D-14); unprotected emissions may differ |
| `ZeroLengthArcsAcyclic` | `zero_arc_instance(params, seed)` with `zero_arcs > 0`, several sources | as `RandomDirectedStateGraph`; zero-length arcs make many arcs tight, so most copies are protected; lengths are even, zero arcs never close a cycle (Q-21) |

Hand example (research notes 02 sec 2.11 item 5): the two fixtures `chain` and `align` of
`tests/unit/paper/hierarchy/test_eval_invariants.cpp` spell out a correct hook trace, including
slots, splits, supplemental offers and the second origin rule. Their child order is the one of
Sec 4.2: neither node has light arcs, so every child is a DFS root that finishes at once, and
the rank order is the reverse finishing order (Lemma 4.4): [c, P] at R with rep(P) = a and
g(c, P) = D_F(a, c) (5 in `chain`, 11 in `align`), [b, a] at P with g(b, a) = D_F(a, b)
(1 / 3); the leaf layout is c, b, a. Hence the supplemental offer of R's first slot in `align`
goes to rep(P) = a with max{6, 0 + 0 + 11} = 11. Do not write a child order down by hand:
`test_eval_invariants.cpp` checks the fixtures' order, representatives and gaps against
`order_children()`, recomputes every supplemental offer of the traces from that output, and
compares the traces with an `evaluate()` run recorded by `RecordingObserver`, event by event.
The checker alone cannot catch a wrong order: it does not see ranks, representatives or gaps.
