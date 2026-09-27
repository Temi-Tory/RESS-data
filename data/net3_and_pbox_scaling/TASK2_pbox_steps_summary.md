# p-box steps-scaling curve

p-box propagation timing on `counterexample-n15` (15 nodes, 23 edges) as a function of the
discretisation step count, using `PBA.uniform`-based p-box inputs rebuilt fresh at each steps
value (avoiding a stale-discretisation bug where p-box inputs built once and reused across
`PBA.setSteps()` calls silently keep the first call's resolution).

## Result

| steps | time |
|---|---|
| 25 | 1.36s |
| 50 | 6.97s |
| 100 | 45.75s |
| 200 | 300.75s |

Growth is markedly super-quadratic: steps 200/25 = 8x the steps, but ~221x the time, consistent
with an empirical exponent of ~2.8 rather than the quadratic (~64x) scaling assumed elsewhere.
steps=800 was not attempted (~4h extrapolated).

## Reproduce

```
julia --project=InfoPropFrmwrk validation/probability/task2_pbox_steps_confirm.jl
```

## Artifacts

- `validation/probability/task2_pbox_steps_confirm.jl` — the script
- `validation/probability/task2_pbox_steps_confirm_results.csv` — the timing table above
