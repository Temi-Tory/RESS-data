# RESS response Task 2: p-box steps-scaling curve, bug-fix reconfirmation, FRESH re-run this
# session, on the 15-node reference network (counterexample-n15).
#
# BACKGROUND (see validation/probability/notes/CORPUS_CAMPAIGN_HANDBACK.md sec 2.1 /
# validation/timing_imprecise.jl's own header comment): the ORIGINAL bug was that
# timing_imprecise.jl built p-box inputs ONCE outside the steps loop and only called
# PBA.setSteps() between legs -- but PBA.uniform() bakes the CURRENT global step count into the
# pbox's discretization at CONSTRUCTION time, so reusing inputs built before the loop meant every
# leg after the first silently kept measuring the FIRST leg's resolution, not the one the column
# name claimed. timing_imprecise.jl has SINCE been patched in-place (its own header carries a
# "BUGFIX (2026-08-17)" note and now rebuilds npb/lpb fresh after every PBA.setSteps() call) --
# but that script measures on random_n20/random_n25 grids, not the 15-node reference network, and
# its own steps-scaling curve section only sweeps one graph (random_n20) at steps
# (25,50,100,200,400,800), not the specific {25,50,100,200} x counterexample-n15 combination this
# task calls for.
#
# The actual script that FIRST produced the faithful 1.23/5.98/39.2/275.7s figures (steps
# 25/50/100/200 on counterexample-n15, the "15-node reference network") is
# validation/fresh_20260816/pbox_steps_probe.jl -- confirmed by grepping its own two log files
# (pbox_steps_probe.log: steps=50 -> 5.98s, steps=200 -> 275.69s; pbox_steps_curve_25_100.log:
# steps=25 -> 1.23s, steps=100 -> 39.16s). That script was ALREADY CORRECT (it builds fresh
# np/lp -- via a fresh call to run_steps(s), which calls tri_pbox(...) fresh -- inside each
# per-steps call, no reuse-before-loop bug), it was just run as two separate half-sweeps across
# two nights. This script adapts it verbatim (same tri_pbox construction, same
# PBOX_COND_BLEND[]=:positive setting, same network/seed) into ONE process covering all four
# steps values in sequence, with results persisted to CSV as they complete.
const REPO = raw"C:\Development\Info_Prop_Framework_Project\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(joinpath(REPO, "InfoPropFrmwrk"))
include(joinpath(REPO,"InfoPropFrmwrk","src","Algorithms","InfoPropFramework.jl")); using .InfoPropFramework
using Random, Printf, Dates, Distributions
include(joinpath(REPO,"validation","graph_gen.jl")); include(joinpath(REPO,"validation","oracles.jl"))

InfoPropFramework.InputProcessingModule.PBOX_COND_BLEND[] = :positive  # cvxP, tight -- matches original probe

const OUTDIR = joinpath(REPO, "validation", "probability")
const RESSDIR = raw"C:\Development\Info_Prop_Framework_Project\Info_Prop_Framework_Project\InfoPropFrmwrk\Publications\My work\RESS_response\pre-write final\data\task2_pbox_steps"
mkpath(RESSDIR)
log_path = joinpath(OUTDIR, "task2_pbox_steps_confirm.log")
io_log = open(log_path, "w")
function logmsg(s...)
    msg = string("[", now(), "] ", s...)
    println(msg); println(io_log, msg); flush(stdout); flush(io_log)
end

logmsg("=== Task 2: p-box steps-scaling curve, fresh reconfirm, counterexample-n15 ===")

g = load_edges("cex", joinpath(REPO,"dag_ntwrk_files","counterexample-n15","counterexample-n15.EDGES"))
logmsg("network: counterexample-n15, nodes=", g.n, " edges=", length(g.edges))
P = make_problem(g); D = draw_probs(P, MersenneTwister(1007))
fk, jn = identify_fork_and_join_nodes(P.outgoing, P.incoming)

tri_pbox(v, w, steps) = begin
    lo = max(0.0, v-w); hi = min(1.0, v+w)
    lo >= hi && return PBA.makepbox(PBA.interval(lo, lo))
    qs = quantile.(TriangularDist(lo, hi, clamp(v, lo, hi)), [(i-0.5)/steps for i in 1:steps])
    PBA.pbox(qs, qs)
end

function run_steps(steps)
    PBA.setSteps(steps)
    # fresh np/lp built AFTER setSteps, for EVERY call -- this is the fix; no np/lp is ever reused
    # across a different steps value.
    np = Dict{Int64,pbox}(n => tri_pbox(v, 0.05, steps) for (n,v) in D.node_priors)
    lp = Dict{Tuple{Int64,Int64},pbox}(e => tri_pbox(v, 0.05, steps) for (e,v) in D.link_probs)
    roots, uniq = new_identify(P.edgelist, np, lp, Set{Int64}(P.sources), fk, jn, P.anc, P.desc, P.itersets)
    t = @elapsed update_beliefs_iterative(P.edgelist, P.itersets, P.outgoing, P.incoming, P.sources, np, lp, P.desc, P.anc, roots, jn, fk, uniq, Dict{CacheKey,DiamondCacheEntry{pbox}}())
    t
end

logmsg("warmup steps=20 (JIT, discarded)")
run_steps(20)

const OLD_REF = Dict(25=>1.23, 50=>5.98, 100=>39.2, 200=>275.7)  # prior faithful figures, for comparison only
csv_path = joinpath(OUTDIR, "task2_pbox_steps_confirm_results.csv")
open(csv_path, "w") do f
    println(f, "steps,t1_s,t2_s,t_min_s,old_ref_s")
    for s in (25, 50, 100, 200)
        logmsg("steps=", s, " starting...")
        t1 = run_steps(s); t2 = run_steps(s)
        tmin = min(t1,t2)
        old = get(OLD_REF, s, NaN)
        @printf(f, "%d,%.4f,%.4f,%.4f,%.2f\n", s, t1, t2, tmin, old); flush(f)
        logmsg(@sprintf("steps=%d: %.2fs / %.2fs (min %.2fs; prior faithful ref: %.2fs)", s, t1, t2, tmin, old))
    end
end
logmsg("wrote ", csv_path)
cp(csv_path, joinpath(RESSDIR, "task2_pbox_steps_confirm_results.csv"); force=true)
logmsg("=== Task 2 done (steps=800 deliberately NOT attempted, per instructions: documented ~4h extrapolated / impractical) ===")
close(io_log)
cp(log_path, joinpath(RESSDIR, "task2_pbox_steps_confirm.log"); force=true)
