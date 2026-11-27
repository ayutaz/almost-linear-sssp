# Analysis Documents for the Paper (docs/research/)

English | [日本語](README.ja.md)

The subject is arXiv 2611.34911v1, "Linear Single-Source Shortest Path, Almost". The paper itself cannot be redistributed, so fetch it into `.cache/paper/` (not tracked by git) with `scripts/fetch_paper.sh`. The documents here paraphrase and explain the paper; they do not reproduce its text. Line numbers `l.NNN` refer to `.cache/paper/2611.34911v1.txt`.

## Index

| Document | Scope | Contents |
|---|---|---|
| [01_preliminaries_primitives.md](01_preliminaries_primitives.md) / [01_preliminaries_primitives.ja.md](01_preliminaries_primitives.ja.md) | Sec 1–3 | Main theorem, list of symbols, how the thresholds are used, base-case conventions, the computation model (length arithmetic, integers, randomness), places where double causes problems, the Sec 3 primitives P1–P13 (search, sorting, selection, Dijkstra, Fibonacci heap, growth paths and MST, union-find, LCA, level ancestor, split–findmin), correspondence with reference implementations, proposed modules |
| [02_ordered_hierarchy_local.md](02_ordered_hierarchy_local.md) / [02_ordered_hierarchy_local.ja.md](02_ordered_hierarchy_local.ja.md) | Sec 4–5 | Ordered hierarchies (Def 4.1–4.2), normalization, child order, Scan and Evaluate (Alg 2–3), checkable invariants I1–I8, owner decomposition, radial state representation, threshold components, local hierarchies (Alg 4), atom replacement, lemmas that break when thresholds are lowered, proposed modules and implementation order |
| [03_pivots_aux_outer.md](03_pivots_aux_outer.md) / [03_pivots_aux_outer.ja.md](03_pivots_aux_outer.ja.md) | Sec 6–8 | Pivot route certificates (Alg 5–6, Lemma 6.2–6.8), the auxiliary metric and owner coordinates, width classes, calendars and atoms, OwnerDistances (Alg 7), the outer state graph and numeric scales, ShortestPaths (Alg 8) and the breakdown of the main theorem's running time, list of intermediate-state invariants |
| [04_outer_scales_pipeline.md](04_outer_scales_pipeline.md) / [04_outer_scales_pipeline.ja.md](04_outer_scales_pipeline.ja.md) | Sec 8.2–8.3 | Implementation spec for the outer scales (Lemma 8.4) and the outer pipeline (Alg 8 l.8–17, Lemma 8.5): the factor F in exact (value, exponent) form, Lemma 8.4 stage by stage, the stages of Alg 8 and where fallbacks are recorded, Mode A / Mode B, what the existing components lack (M-1–M-12), proposed questions Q-65–Q-74, proposed acceptance tests SC-1–14 and PL-1–12 |
| [thresholds.md](thresholds.md) / [thresholds.ja.md](thresholds.ja.md) | Thresholds and constants | List and classification of all constants (whether they affect correctness or only time), the smallest N satisfying each threshold, the paths taken at realistic N, what breaks when values are changed in Mode B, recommendations for the implementation |
| [prior_work_and_oss.md](prior_work_and_oss.md) / [prior_work_and_oss.ja.md](prior_work_and_oss.ja.md) | Prior work and existing OSS | Licenses, verification methods and measurement conditions of existing implementations (including one of the target paper), prior work, recommendations on whether to incorporate or reimplement |
| [open-questions.md](open-questions.md) / [open-questions.ja.md](open-questions.ja.md) | Whole paper | The ambiguities from the documents above, collected in one table with IDs (Q-xx): provisional interpretation, affected modules, and resolution status (Decided by interpretation / Needs agreement / Open) |
| [../PAPER_MAPPING.md](../PAPER_MAPPING.md) / [../PAPER_MAPPING.ja.md](../PAPER_MAPPING.ja.md) | Whole paper | Mapping table from every paper item to modules, tests, Mode, status and deviations (D-xx). The record of unified module names is also here |

## Reading order

1. **§2–§3 of `thresholds.md`** — First take in that "with the paper's thresholds, the outer construction never runs on any realistic input, and everything falls into the base case of scan-based Dijkstra". The reasons for the Mode A / B / C split (ADR 0003) and for the two versions of the base case (ADR 0008, D-1) are here.
2. **§0 and §2 of `01_preliminaries_primitives.md`** — The computation model and notation. Everyone writing any module should read this.
3. **`../PAPER_MAPPING.md`** — Look at the rows for the modules you own, the naming unification (§0), and the list of deviations (§3).
4. The analysis document for your section:
   - Sec 3 primitives → §3 and §5 of `01`
   - Sec 4 evaluator → §2 of `02`
   - Sec 5 local hierarchies → §3–§4 of `02`
   - Sec 6–8 → `03` (for numeric magnitudes, §5–§7 of `thresholds.md`); Sec 8.2–8.3 implementation spec → `04`
5. **`open-questions.md`** — Filter by the column of the modules you own, and reach agreement on "Needs agreement" rows before implementing. Add new ambiguities here with their location (section, lemma, line number) and candidate interpretations.
6. `prior_work_and_oss.md` — Read when deciding comparison targets and measurement conditions (benchmark owners).

## Writing rules

- Keep separate what the paper states, what follows from it, and your own interpretation (the [Supplement], [Interpretation] and [Ambiguous] markers in each document).
- Do not reproduce long passages. Limit yourself to quoting short phrases and referring to sections, lemmas and algorithm line numbers.
- Do not write absolute paths or personal information (this is a public OSS project). Run `uv run python scripts/check_public_safety.py` before pushing.
