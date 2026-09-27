# Grid case-study p-box accuracy table

Accuracy/soundness comparison on the benchmark grid (16 nodes, 24 edges, sources {1,3,13}), node
priors exact 1.0, link priors 0.9 uncertain, p-box discretisation steps=200, target node=16, Monte
Carlo N=50,000, at two interval half-widths (w=0.05, w=0.10). Exactness reference is an independent
BDD oracle (`BinaryDecisionDiagrams.jl`); both are exact symbolic evaluations of the same
reachability formula, so this is a valid substitute for a CUDD-based check.

## Result

| regime | w | method | worst_overwidth_or_err | worst_unsound | note |
|---|---|---|---|---|---|
| Float64 | - | IPA_vs_sifted_BDD | 1.110e-16 | - | exact |
| Interval | 0.05 | IPA_vs_BDDcorners | 3.331e-16 | 1.110e-16 | exact |
| Interval | 0.05 | naive_no_conditioning | 2.247e-04 | - | over-wide |
| **pbox@200** | **0.05** | **IPA_vs_MC** | - | **0.000e+00** | **sound** |
| pbox@200 | 0.05 | naive_vs_MC | - | 2.138e-02 | unsound |
| Interval | 0.10 | IPA_vs_BDDcorners | 2.220e-16 | 1.110e-16 | exact |
| Interval | 0.10 | naive_no_conditioning | 1.110e-16 | - | over-wide |
| **pbox@200** | **0.10** | **IPA_vs_MC** | - | **0.000e+00** | **sound** |
| pbox@200 | 0.10 | naive_vs_MC | - | 1.724e-02 | unsound |

| w | IPA (sound cvxP) propagate | naive propagate | MC_N=50,000 |
|---|---|---|---|
| 0.05 | 531.24s | 1.26s | 53.11s |
| 0.10 | 607.30s | 1.24s | 52.14s |

`worst_unsound` is 0.000e+00 for both widths under the current cvxP/cvxF conditioning operator —
i.e. the p-box propagation is sound at both tested uncertainty widths. The `naive_vs_MC` rows
(which bypass diamond conditioning entirely) remain unsound, as expected, confirming the diamond
conditioning is what produces soundness here. The p-box conditioning cost (531–607s at steps=200)
is the dominant cost in this table, consistent with the steps-scaling finding above.

## Reproduce

```
julia --project=InfoPropFrmwrk validation/probability/grid_case_study/task3_grid_accuracy_confirm.jl
julia --project=InfoPropFrmwrk validation/probability/grid_case_study/task3_grid_accuracy_w10_only.jl
```

## Artifacts

- `validation/probability/grid_case_study/task3_grid_accuracy_confirm.jl` — script for Float64 + w=0.05
- `validation/probability/grid_case_study/task3_grid_accuracy_w10_only.jl` — script for w=0.10
- `validation/probability/grid_case_study/data/grid_accuracy.csv` — the full accuracy table, both widths
