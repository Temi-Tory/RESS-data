# Grid case study — the methodology demonstrator (self-contained)

Purpose: on one featured graph (the paper's grid), show the method is correct, capable, and its
cost is understood. This is a method demonstration, so inputs are illustrative (the Float64 values
match the paper's `dPrPm`; the imprecise extensions are this project's own).

## Input construction

Float64 = the paper's grid values (to match `dPrPm`). Non-float extensions are built around those
values:

| component                          | Float64      | Interval        | p-box                                   |
|------------------------------------|--------------|-----------------|-----------------------------------------|
| uncertain (links 0.9, interior nd) | paper `v`    | `[v-w, v+w]`    | triangular(min=`v-w`, mode=`v`, max=`v+w`) |
| perfect (1.0 sources/nodes)        | 1.0          | 1.0 (exact)     | 1.0 (degenerate, mode=1)                |

- Half-widths **w in {0.05, 0.10}** (both are run), clamped to [0,1].
- Perfect nodes stay exactly 1.0 in every representation, so Float64 is the all-modes special case
  of the imprecise model.
- The p-box triangular distribution is mode-centred at the paper value, so the point estimate is
  the most-likely value.

## Oracles (each capability validated against an independent ground truth)

- Float64 exactness -> sifted ROBDD (CUDD). Confirmed to 1.1e-16.
- Interval exactness -> sifted ROBDD at the two corners (all-low / all-high; exact range by
  monotonicity). Confirmed exact to 1e-16.
- p-box soundness -> the cvxP/cvxF conditioning operator gives `worst_unsound=0.000e+00` at both
  w=0.05 and w=0.10 against a CUDD-based accuracy rerun of this suite — see `data/grid_accuracy.csv`.

## The "increasing capability" arc (paper section order)

1. `dPrPm` baseline: published grid numbers + accessibility caveat (not reproducible, motivating a
   reproducible exact method).
2. IPA vs sifted-CUDD: exact agreement (accuracy) + performance.
3. Interval: exact belief range vs naive over-widening (interval overhead ~1.2x).
4. p-box @ {50,200,800}: soundness vs MC + tightness/cost tradeoff.
5. Monte Carlo overlay: the ground truth once inputs are imprecise.

## Run order & scripts (two separate concerns)

DATA / CORRECTNESS run (can use threads) -> `run_grid_data.jl`: produces the case-study table —
per-node Float64/Interval/p-box beliefs + BDD / BDD-corner / MC / naive comparisons + exactness &
soundness columns. Output -> `data/`.

COST run (must run alone, single-thread, warmup first) -> `bench_grid.jl`: (1) one warmup call per
type (excludes JIT), then (2) the trio — `@benchmark` (median time + memory + allocs), a sampling
profiler (time % PBA vs IPA vs other), and an allocation profiler (bytes PBA vs IPA vs other).
Output -> `data/` + `notes/profile_breakdown.md`. The cost run measures wall-clock, so nothing else
(including the data run) should run concurrently with it.

## Complexity validation

IPA per-instance cost is exactly `Work = sum over diamonds d of 2^|C_d| * O(|E_d|)` (`|C_d|` =
conditioning set size, computed by `new_identify`). Worst case `max|C_d| <= treewidth`. The
validation table compares measured `ipa_ops` against this formula's prediction across the grid and
corpus, confirming the model is exact rather than a heuristic approximation (there is no closed
form in `n` alone: the underlying problem is #P-hard).
