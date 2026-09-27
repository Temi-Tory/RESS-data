# Net3 (EPANET water benchmark) — identify-only feasibility check

Diamond identification (`new_identify`) on the flow-oriented Net3 conversion
(`net3.EDGES`, 97 nodes / 119 edges — the hydraulically-grounded conversion, not the older
arbitrary-BFS `net3-water` conversion, which gives a different, less-realistic diamond structure:
51 diamonds/maxcond 5 vs this one's 307/12). Priors: non-degenerate uniform 0.9 on all nodes and
edges. Identification completed well within budget, so Float64 exact propagation was run as well.

## Result

| Metric | Value |
|---|---|
| Nodes | 97 |
| Edges | 119 |
| Sources (indegree 0) | 3 — {1, 2, 4} = Lake, River, tank "1" |
| Sinks (outdegree 0) | 18 |
| Forks | 35 |
| Joins | 23 |
| Layers | 28 |
| **Unique diamonds (incl. nested)** | **307** |
| Maximal/root diamonds | 21 |
| **maxcond (max \|conditioning_nodes\|)** | **12** |
| sum(2^\|C\|) over unique diamonds | 5.4734e+04 |
| Identify time (`new_identify`) | 1.348 s |
| Float64 exact propagation time | 2.886 s (completed) |

p-box propagation was not attempted on Net3: `sum(2^|C|)=5.47e4` exceeds the cost at which the
comparably-costed drone K=8 configuration failed to complete.

## Reproduce

```
julia --project=InfoPropFrmwrk validation/net3/thesis_writing_pack/net3_task1_confirm_rerun.jl
```

## Artifacts

- `validation/net3/thesis_writing_pack/net3_task1_confirm_rerun.jl` — the script
- `validation/net3/thesis_writing_pack/net3_task1_identify_only.csv` — identify-only summary row
- `validation/net3/thesis_writing_pack/net3_task1_full_summary.csv` — full summary incl. propagation status
- `validation/net3/thesis_writing_pack/net3_task1_propagation_beliefs.csv` — per-node Float64 beliefs (97 rows)
