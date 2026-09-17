# RESS reproduction package — Information Propagation Algorithm (IPA)

Reproduction data, scripts, and the current manuscript draft for the RESS journal submission
"Information Propagation Algorithm for Exact Reliability Analysis in Directed Acyclic Process
Networks" (Ohiani, Patelli, Pyke). This folder is a staging copy assembled for a future Zenodo
deposit — **not yet minted, not yet a git repository.** The author will `git init` and push this
themselves.

- **Code**: the algorithm itself is released separately as the open-source Julia package
  `InformationPropagationAnalysis.jl`, registered in the Julia General package registry (current
  version 0.2.1): https://github.com/Temi-Tory/InformationPropagationAnalysis.jl
- **This deposit**: the benchmark networks, validation scripts, and result data behind every
  table and figure in the manuscript, so every number in the paper can be traced to how it was
  produced.

## Folder structure

- `manuscript/` — the current manuscript draft (`main.tex`/`main.pdf`, 44pp, compiles clean),
  its figures, the original diagram assets, a log of every correction applied during the
  pre-submission audit (`CORRECTIONS_APPLIED.md`), and a proposed (not yet applied) addition
  explaining p-box computational cost (`PBOX_COST_MECHANISM_DRAFT.md`).
- `networks/` — every network topology with real validated data supporting a manuscript claim
  (34 folders): the benchmark grid, power-network, KarlNetwork, counterexample-n15, the
  worked-example network (`test-decomp3s2t`), the applied drone case study (3 configurations + the
  K=6 comparison variant, + `drone-network-full` for the unrestricted-connectivity finding),
  real-infrastructure networks not yet individually named in the manuscript text
  (`mlgw-gas-network`, `metro_directed_dag_for_ipm`, `net3`), and the full 17-network bnlearn
  benchmark set (8–724 nodes) — see `MANIFEST_DECISIONS.md` for the include/exclude rule (every
  topology with real validated data is in; only networks with no data, invalidated inputs, a
  different algorithm's validation, or a redundant duplicate are out) and two open manuscript-text
  recommendations it surfaced (name the real-infrastructure networks explicitly; add a bnlearn
  breadth claim — a strong, already-validated answer to reviewer comment 1/R3.2's request for more
  benchmark diversity). The 129-graph synthetic corpus and the fanin-k/mesh-w adversarial families
  are regenerated on the fly by seeded scripts in `scripts/`, not stored as static files.
- `data/` — every CSV/log behind a manuscript table or figure:
  - `final_csvs/` — the "as of the applied revision" data pack (129-graph corpus, timing
    comparisons, p-box soundness/tightness sweeps, drone belief tables).
  - `grid_case_study/` — the benchmark grid's own dedicated accuracy/cost suite.
  - `fresh_this_session/` — artifacts from the final pre-submission confirmation pass (Net3
    feasibility check, the corrected p-box steps-scaling curve, the grid p-box soundness
    re-confirmation at both tested uncertainty widths, and the p-box-cost-vs-diamond-count
    experiment), each with a written summary of exactly what was run and how.
- `scripts/` — the Julia scripts that generate the above: corpus generation, the reference
  recursion used to validate the algorithm, independent oracles (BDD, Monte Carlo, path
  enumeration), and the case-study drivers (grid, drone K-sweep, ASCE reproduction, Net3,
  p-box timing).
- `literature/CITATION_STATUS.md` — final citation audit: what's confirmed and already in the
  manuscript, what's confirmed but not yet integrated (recommended additions), and the resolution
  of two previously-open citations.
- `reviewer_response/` — the point-by-point response tracker (all 3 reviewers), three flagged
  framing recommendations awaiting author sign-off, and a cover-letter draft.
- `MANIFEST_DECISIONS.md` — the include/exclude reasoning for this package (which of ~68 locally
  available network folders were bundled and why; a few genuinely excluded as superseded/
  exploratory and not referenced in the current manuscript).

## How to reproduce

1. Install Julia (version per `InformationPropagationAnalysis.jl`'s `Project.toml`, currently
   1.12) and the framework environment the scripts activate (`InfoPropFrmwrk/` in the source
   repository, or `InformationPropagationAnalysis.jl` directly for new work).
2. Each script in `scripts/` is self-contained and states its own network/parameters in its
   header comment. Run order for the core validation claims:
   - `graph_gen.jl` / `graph_families.jl` — regenerate the synthetic corpus (seeded, deterministic).
   - `full_regression_sifted.jl`, `families_validate.jl`, `large_graphs.jl` — exactness vs. the
     independent BDD oracle (§5.2 of the manuscript).
   - `complexity_validate.jl` — diamond/conditioning-set structural statistics per network.
   - `probability/asce_grid_reproduction.jl`, `probability/asce_power_reproduction.jl` — the
     ASCE benchmark reproductions (§5.1).
   - `probability/drone_ksweep_remeasure.jl`, `drone_bdd_comparison.jl` — the applied case study
     (§5.4).
   - `probability/task2_pbox_steps_v2_one.jl <steps>` and
     `probability/pbox_cost_vs_diamonds_v2_one.jl <network>` — p-box cost characterisation
     (§5.3.2). **Run each invocation in its own fresh process** — see the header comment for why
     (a documented same-process timing-pollution failure mode this project has hit twice).
3. Compile the manuscript: `pdflatex main.tex` (twice, for cross-references) from `manuscript/`.

## Headline validated claims (see `manuscript/main.pdf` for full context)

- Exactness: 129/129 synthetic corpus graphs, six topological families, and real infrastructure
  networks match an independent reduced-ordered BDD to floating-point precision (worst
  disagreement 1.1e-16).
- Interval propagation is exact by construction (monotonicity) and, for one-shot queries, 3.9x–
  99.9x faster than the equivalent BDD-corner-evaluation route across 8 corpus families.
- p-box (probability-box) propagation gives guaranteed, Fréchet-bound-backed distributional
  bounds; sound on every tested configuration (grid at two uncertainty widths, KarlNetwork,
  mlgw-gas-network, drone networks). Cost is controllable via the discretisation level (steps);
  the empirical growth exponent is approximately 2.6, and cost tracks the same realised-work
  quantity (diamond-cache hits, not raw conditioning-set-width) that governs the exact/interval
  cost story.
- Applied case study: a medical drone logistics network for Scotland, with every input traced to
  the cited source design study; a network-redundancy parameter locates the practical boundary of
  exact computation, beyond which the decision-diagram alternative does not complete.

## License

MIT, matching the package (`InformationPropagationAnalysis.jl`) and the companion thesis-data
Zenodo repository (`Temi-Tory/thesis-data`, DOI 10.5281/zenodo.22180227).

## Citation

Once published, cite the paper. For the software: Ohiani & Patelli,
`InformationPropagationAnalysis.jl` (Julia General registry). For this specific reproduction
package: this deposit's own DOI, once minted.

## Status note (for the author, not for the deposit's public README)

This package was assembled from the working repository as part of a pre-submission confirmation
pass; it has not yet had the Zenodo DOI minted, and a few open decisions from that pass are listed
in `MANIFEST_DECISIONS.md`'s closing section. See the working repository's
`InfoPropFrmwrk/Publications/My work/RESS_response/pre-write final/README.md` for the full session
history and remaining open items (judgment calls awaiting sign-off, the adversarial-data scope
question, the still-pending original submission/decision letter). Remove this section before
publishing the deposit publicly.
