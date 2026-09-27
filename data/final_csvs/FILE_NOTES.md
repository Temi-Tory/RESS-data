# File notes

- `structured_fresh.csv` — output of `make_structured.jl`, the source of the power-network and
  KarlNetwork rows of the structured-networks table. Its counterexample-n15 row uses the diagram's
  older 340-node count; the canonical figure used in the manuscript is 323 (see
  `counterexample_n15_canonical_timing.csv`).
- `paper_data.csv` — counterexample-n15 `sift_nodes` is 323, matching the manuscript.
- `robdd_named_networks_20260921.csv` — sifted CUDD diagram sizes and IPA-vs-ROBDD agreement for
  the gas, transit, and Net3 networks (uniform 0.9 inputs).
- `grid_ipa_vs_cudd_cost_20260921.csv` — IPA vs CUDD cost on the benchmark grid (time, Julia-side
  memory, and allocations).
- `bnlearn_vs_sifted_bdd_20260921.csv` — IPA (Float64 and interval) against sifted CUDD on 18
  bnlearn-family networks including munin-dag/munin-sub1/water; produced by `bnl_oracle.jl`.

Scripts in `scripts/probability/` that produced the above: `grid_bench3.jl`, `net_bdd.jl`,
`worked16.jl` (§4 worked example), `drone_C.jl` (drone conditioning sizes: centralised 8, dense 16,
minimal 16), `bnl_oracle.jl` (bnlearn comparison). Path constants in these scripts point at a local
scratch layout and need adjusting before rerunning.
