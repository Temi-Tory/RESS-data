# Task 2 v2 -- p-box steps-scaling curve, CLEAN methodology: ONE steps value per fresh process.
#
# Fixes a real measurement-pollution flaw in task2_pbox_steps_confirm.jl, which ran a warm-up
# call PLUS eight timed propagations (2 calls x 4 steps values, increasing size) all in ONE Julia
# process. This project's own prior campaign already documented this exact failure mode and named
# the rule (validation/probability/MASTER_FINDINGS.md): "float corner legs ran after the interval
# leg's multi-GB cache in the same process: 2,291s vs 81s fresh-process reference ... process-
# hygiene rule, now twice learned: one timing measurement per fresh process."
#
# Fix applied here: this script measures EXACTLY ONE (network, steps) timed propagation per
# process invocation. JIT is warmed up first on a TINY, DIFFERENT, throwaway network (never the
# real target, never re-run at a different size) so the timed call itself is not the first-ever
# call into the generic propagation code path, without building a large heap/cache on the actual
# target before measuring it. Run this script once per steps value (25, 50, 100, 200), each in
# its own fresh `julia` process invocation, SEQUENTIALLY (not concurrently -- concurrent Julia
# processes competing for CPU would itself distort wall-clock timing).

const REPO = raw"C:\Development\Info_Prop_Framework_Project\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(joinpath(REPO, "InfoPropFrmwrk"))
include(joinpath(REPO,"InfoPropFrmwrk","src","Algorithms","InfoPropFramework.jl")); using .InfoPropFramework
using Random, Printf, Dates, Distributions
include(joinpath(REPO,"validation","graph_gen.jl")); include(joinpath(REPO,"validation","oracles.jl"))

InfoPropFramework.InputProcessingModule.PBOX_COND_BLEND[] = :positive  # shipped default, unchanged

length(ARGS) >= 1 || error("usage: julia task2_pbox_steps_v2_one.jl <steps>")
const STEPS = parse(Int, ARGS[1])

const OUTDIR = joinpath(REPO, "validation", "probability")
const OUT_CSV = joinpath(OUTDIR, "data", "task2_pbox_steps_v2_results.csv")
const LOG_PATH = joinpath(OUTDIR, "task2_pbox_steps_v2_steps$(STEPS).log")
io_log = open(LOG_PATH, "w")
function logmsg(s...)
    msg = string("[", now(), "] ", s...)
    println(msg); println(io_log, msg); flush(stdout); flush(io_log)
end
logmsg("=== Task 2 v2: steps=$STEPS -- ONE clean process, warm-up on a TRIVIAL SEPARATE network ===")

tri_pbox(v, w, steps) = begin
    lo = max(0.0, v-w); hi = min(1.0, v+w)
    lo >= hi && return PBA.makepbox(PBA.interval(lo, lo))
    qs = quantile.(TriangularDist(lo, hi, clamp(v, lo, hi)), [(i-0.5)/steps for i in 1:steps])
    PBA.pbox(qs, qs)
end

# --- JIT warm-up: a tiny, throwaway 5-node random DAG, NEVER the real target ---
PBA.setSteps(STEPS)
gwarm = gen_random_dag(MersenneTwister(999); n=5, p=0.4)
Pwarm = make_problem(gwarm)
fkw, jnw = identify_fork_and_join_nodes(Pwarm.outgoing, Pwarm.incoming)
npw = Dict{Int64,pbox}(n => tri_pbox(0.9, 0.05, STEPS) for n in Pwarm.all_nodes)
lpw = Dict{Tuple{Int64,Int64},pbox}(e => tri_pbox(0.9, 0.05, STEPS) for e in keys(Pwarm.eid))
logmsg("warm-up (trivial 5-node network) starting...")
t_warm = @elapsed begin
    rootsw, uniqw = new_identify(Pwarm.edgelist, npw, lpw, Set{Int64}(Pwarm.sources), fkw, jnw, Pwarm.anc, Pwarm.desc, Pwarm.itersets)
    update_beliefs_iterative(Pwarm.edgelist, Pwarm.itersets, Pwarm.outgoing, Pwarm.incoming, Pwarm.sources, npw, lpw, Pwarm.desc, Pwarm.anc, rootsw, jnw, fkw, uniqw, Dict{CacheKey,DiamondCacheEntry{pbox}}())
end
logmsg(@sprintf("warm-up done: %.4fs (discarded, trivial network -- NOT the real target)", t_warm))

# --- the real target: counterexample-n15, the "15-node reference network", ONE timed call ---
g = load_edges("cex", joinpath(REPO,"dag_ntwrk_files","counterexample-n15","counterexample-n15.EDGES"))
P = make_problem(g); D = draw_probs(P, MersenneTwister(1007))
fk, jn = identify_fork_and_join_nodes(P.outgoing, P.incoming)
np = Dict{Int64,pbox}(n => tri_pbox(v, 0.05, STEPS) for (n,v) in D.node_priors)
lp = Dict{Tuple{Int64,Int64},pbox}(e => tri_pbox(v, 0.05, STEPS) for (e,v) in D.link_probs)

logmsg("target network: counterexample-n15, nodes=", g.n, " edges=", length(g.edges))
logmsg("SINGLE timed call starting (steps=$STEPS)...")
t = @elapsed begin
    roots, uniq = new_identify(P.edgelist, np, lp, Set{Int64}(P.sources), fk, jn, P.anc, P.desc, P.itersets)
    update_beliefs_iterative(P.edgelist, P.itersets, P.outgoing, P.incoming, P.sources, np, lp, P.desc, P.anc, roots, jn, fk, uniq, Dict{CacheKey,DiamondCacheEntry{pbox}}())
end
logmsg(@sprintf("timed call done: %.4fs", t))

header_needed = !isfile(OUT_CSV)
open(OUT_CSV, "a") do f
    header_needed && println(f, "steps,pbox_propagate_seconds,warmup_seconds_trivial_network")
    @printf(f, "%d,%.4f,%.4f\n", STEPS, t, t_warm)
end
logmsg("wrote row to $OUT_CSV")
close(io_log)
