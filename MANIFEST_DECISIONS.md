# Data package scope — final rule and decisions (2026-09-06)

**Rule, converged in discussion with the author**: include every network topology with real
validated data (belief/comparison/timing) that supports a claim this manuscript makes or is
adding. There is no meaningful "corpus vs. benchmark" distinction to exclude on — the 129-graph
synthetic corpus, the six named families, and KarlNetwork are *themselves* synthetic,
non-decision-relevant, topology-only evidence; bnlearn is the same kind of evidence at a larger
and more diverse scale, not a different category. Decision-relevant reliability context is the
job of the two dedicated case studies (the grid methodology demonstrator, the drone applied case
study) — everything else in the corpus exists to validate exactness/cost claims across topology
diversity, regardless of whether it is "real" or synthetic, named individually in prose or not.

Exclusion is reserved for networks with a **substantive** reason, not a categorical one:

## Included

- **Named/discussed case-study and reference networks**: `grid-graph` (benchmark, §5.1),
  `power-network` (§5.2), `KarlNetwork` (§5.2), `counterexample-n15` (§5.3.2),
  `test-decomp3s2t` (Figure~12 worked example, §4.2), the applied drone case study
  (`drone-network-fw-reliant-centralized`, `-vtol-dense-decentralized`, `-concentrated-minimal`
  + `-concentrated-minimal-k6-test`, §5.4), and `drone-network-full` (the unrestricted-
  connectivity network the manuscript states did not complete identification).
- **Real-infrastructure networks not yet individually named in the manuscript text but backed by
  real validated data**: `mlgw-gas-network` (Memphis gas network), `metro_directed_dag_for_ipm`
  (Berlin metro), `net3` (EPANET water benchmark, trimmed to the actual network files — the
  working repo's copy also contains a ~255MB unrelated server-response dump, excluded).
  **Open manuscript TODO**: §5.2 claims "several real infrastructure networks" but only
  `power-network` is currently named as such — recommend naming these three explicitly in the
  rewrite so the claim is concretely traceable.
- **The full bnlearn benchmark set (17 networks, `*-bnlearn` + `munin-dag`/`munin-sub1`/`water`)**:
  real, published Bayesian-network topologies (bnlearn.com/bnrepository), 8 to 724 nodes,
  synthetic (non-decision-relevant) reliability values assigned for structural-exactness testing.
  15/17 propagate cleanly; `diabetes-bnlearn` and `andes-bnlearn` are honestly documented
  boundary/intractable cases (memory/time limits respectively), not hidden. **Currently unclaimed
  in the manuscript** — recommended addition: a separately-framed breadth statement (e.g.
  "structural exactness additionally verified against 17 published Bayesian-network benchmark
  topologies up to 724 nodes"), distinct from the real-infrastructure and applied-case-study
  claims. This is a genuinely stronger answer to reviewer comment 1 / R3.2's request for
  "additional benchmark networks with varying sizes and topological complexities" than the
  current corpus alone, at zero extra validation cost (already done).
- **Synthetic/generated, no static files needed**: the 129-graph corpus and the adversarial
  families (`fanin-k`, `mesh-w`, pending the agreed manuscript addition) — seeded and regenerated
  by `scripts/graph_gen.jl` / `graph_families.jl` / `graph_adversarial.jl`.

## Excluded — each for a substantive reason, not a categorical one

- **Exploratory drone variants** (`drone-network-balanced-k3`, `-cost-optimal`, `-geographic-knn`,
  `-resilience-optimal-k5`, `-time-optimal-k2`, `drone-medical-delivery-network`,
  `single-mission-drone-network`, `regional_hub_drone_medical`) — **verified by grepping the
  entire working repo's `validation/` tree**: zero references in any script, log, or CSV. No
  validated data exists for these at all; nothing to include.
- **The superseded six-Pareto case study** (`pareto-point-1..6-*`) — has historical data, but the
  reliability *inputs* assigned to it were explicitly determined to be scientifically
  indefensible (an invented distance-decay probability curve, an unjustified node-prior
  heuristic) and the case study was rebuilt from scratch for exactly this reason
  (`RESS_edit_proposals.md`'s own scope decision). The exclusion reason is "the numbers were
  wrong," not "we didn't use this topology."
- **CPM/Flow-thesis-chapter networks** (`HB0_local_*`, `central_scotland_*`, `edinburgh_area`,
  `glasgow_area`, `glasgow_to_shetland_extreme`, `highland_to_lowland_full_network`,
  `military_multi_domain_network`, `continental_medical_network`, `hybrid_power_hierarchical`,
  `ergo-proxy-dag-network`, `power-network-scenarios`, `psplib-j301_1`, `water-highvdemo`,
  `grid-graph-5x5`) — these validate a **different algorithm** in the broader framework
  (critical-path/scheduling, max-flow/min-cut) — there is no reliability/belief data for them
  relative to anything this paper claims. Not applicable, not a withheld piece of evidence.
- **`net3-water`** (an older, less-realistic arbitrary-BFS conversion of the same real Net3
  system) — has data, but is a redundant, inferior duplicate of the same real system already
  represented, properly, by `net3`'s hydraulically-grounded conversion (307 diamonds/maxcond=12
  vs. the older conversion's 51/5, which the working repo's own notes flag as giving a "different,
  less-realistic diamond structure"). Including both would add noise, not evidence.

## Added since the above was written: ISCAS85 combinational circuits (adversarial family, not real-infrastructure)

Real, published, widely-cited digital-circuit benchmark suite (`pld.ttu.ee/~maksim/benchmarks/iscas85/`),
converted this session and tested under the **same non-degenerate uniform-0.9 prior convention**
used throughout the corpus — a deliberate worst-case stress test (prevents the zero-weight-skip
optimization from pruning any diamond). Classify with `fanin-k`/`mesh-w` (adversarial), not with
`mlgw`/`metro`/`net3` (real infrastructure): the priors are chosen to maximise structural
complexity, not to represent real circuit-reliability semantics. Result is a genuine, citable
finding: `c432` (196 nodes) alone reaches `maxcond=67`/`sum_2^C≈1.17e21`; `c1355`/`c1908`
(587/913 nodes) crash **identification itself** via memory exhaustion on a 16GB machine — a harder
wall than diabetes-bnlearn. `c17` (11 nodes) is the one fully-validated data point (Float64,
Interval, and p-box all confirmed; p-box cost matched this session's `measured_ops`-based
prediction almost exactly). Included: `iscas85-c17` (full validation), `iscas85-c432`/`c499`/`c880`
(identify-only, adversarial-boundary structural data), `iscas85-c1355`/`c1908` (identify-only,
crashed — the crash itself is the data point). Not included: `c2670`–`c6288` (not attempted, no
data — see `data/corpus_expansion/iscas85_identify_only_DNF.csv` for why).

## Pending

- Add bnlearn's breadth claim + the real-infrastructure network names (`mlgw`, `metro`, `net3`,
  now freshly re-confirmed this session with proper one-clean-process timing) to the manuscript
  text — a rewrite-session task, not resolved by this data-scoping pass.
- Decide where ISCAS85 goes in the manuscript (adversarial section placement is already an open
  item — this is additional content for that same decision, not a new one).
- Once the manuscript's adversarial addition exists, add the fanin-k/mesh-w fresh timed data
  (`validation/fresh_20260816/adversarial_timed.csv` + `adversarial_fit_summary.txt` in the
  working repo) here alongside the ISCAS85 data already present.
- Zenodo DOI not minted — author will handle repo creation, push, and minting directly (no Claude
  git authorship).
