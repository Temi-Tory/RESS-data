# RESS response follow-up v2: p-box propagation TIME at a FIXED discretisation (steps=50) vs.
# diamond-count/conditioning structure, cross-network -- CLEAN methodology fix.
#
# The v1 script (pbox_cost_vs_diamonds_one.jl) ran a warm-up call THEN a timed call on the SAME
# large target network, in the same process -- exactly the measurement-pollution anti-pattern
# already documented in this project (MASTER_FINDINGS.md: "one timing measurement per fresh
# process"; the mesh-timing incident showed a 2,291s vs 81s same-process-vs-fresh-process gap).
# Fix: JIT warm-up happens on a tiny, throwaway, DIFFERENT network every time -- never the real
# target -- so the target network gets exactly ONE propagation call in its process, ever.
#
# Run ONE network per process invocation (ARGS[1] = network name), exactly as v1 did (that part
# was already correct -- each network already got its own fresh process from the driver).

const REPO = raw"C:\Development\Info_Prop_Framework_Project\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(joinpath(REPO, "InfoPropFrmwrk"))
include(joinpath(REPO,"InfoPropFrmwrk","src","Algorithms","InfoPropFramework.jl")); using .InfoPropFramework
using Random, Printf, Dates, Distributions
include(joinpath(REPO,"validation","graph_gen.jl"))
include(joinpath(REPO,"validation","graph_families.jl"))
include(joinpath(REPO,"validation","oracles.jl"))

InfoPropFramework.InputProcessingModule.PBOX_COND_BLEND[] = :positive  # shipped default, unchanged

const OUTDIR = joinpath(REPO, "validation", "probability")
const OUT_CSV = joinpath(OUTDIR, "data", "pbox_cost_vs_diamonds_v2.csv")
const LOG_PATH = joinpath(OUTDIR, "pbox_cost_vs_diamonds_v2.log")

function logmsg(io_log, s...)
    msg = string("[", now(), "] ", s...)
    println(msg); println(io_log, msg); flush(stdout); flush(io_log)
end

# --- EXACT same 14-graph corpus construction as complexity_validate.jl / v1 ---
function build_graphs()
    graphs = Tuple{String,Graph}[
      ("grid_4x4", load_edges("grid",joinpath(REPO,"dag_ntwrk_files","grid-graph","grid-graph.EDGES"))),
      ("grid_5x5", gen_grid(5,5)),
      ("counterexample", load_edges("cex",joinpath(REPO,"dag_ntwrk_files","counterexample-n15","counterexample-n15.EDGES"))),
      ("bridge_5", gen_bridge(5)), ("complete_8", gen_complete(8)),
      ("layered_5x4", gen_layered(MersenneTwister(2); layers=5, width=4, p=0.5)),
    ]
    for n in (12,15,20,25), s in 1:2; push!(graphs, ("random_n$(n)_s$s", gen_random_dag(MersenneTwister(s); n=n, p=0.15))); end
    graphs
end

tri_pbox(v, w, steps) = begin
    lo = max(0.0, v-w); hi = min(1.0, v+w)
    lo >= hi && return PBA.makepbox(PBA.interval(lo, lo))
    qs = quantile.(TriangularDist(lo, hi, clamp(v, lo, hi)), [(i-0.5)/steps for i in 1:steps])
    PBA.pbox(qs, qs)
end

function main()
    length(ARGS) >= 1 || error("usage: julia pbox_cost_vs_diamonds_v2_one.jl <network_name> [steps]")
    target = ARGS[1]
    steps = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 50

    graphs = build_graphs()
    idx = findfirst(t -> t[1] == target, graphs)
    idx === nothing && error("unknown network name: $target (known: $(join(first.(graphs), ", ")))")
    nm, g = graphs[idx]

    io_log = open(LOG_PATH, "a")
    logmsg(io_log, "=== network=$nm steps=$steps (v2: clean, trivial-network warm-up) ===")

    P = make_problem(g)
    fk, jn = identify_fork_and_join_nodes(P.outgoing, P.incoming)

    np90 = Dict{Int64,Float64}(n=>0.9 for n in P.all_nodes)
    lp90 = Dict{Tuple{Int64,Int64},Float64}(e=>0.9 for e in keys(P.eid))
    r90, u90 = new_identify(P.edgelist, np90, lp90, Set{Int64}(P.sources), fk, jn, P.anc, P.desc, P.itersets)
    conds = [length(d.diamond.conditioning_nodes) for d in values(u90)]
    ndia = length(conds); maxc = isempty(conds) ? 0 : maximum(conds)
    sum2C = isempty(conds) ? 0.0 : sum(2.0^c for c in conds)
    logmsg(io_log, "structural: V=$(length(P.all_nodes)) E=$(length(g.edges)) n_diamonds=$ndia maxcond=$maxc sum_2C=$sum2C")

    PBA.setSteps(steps)

    # --- JIT warm-up on a tiny, DIFFERENT, throwaway network -- never the real target ---
    gwarm = gen_random_dag(MersenneTwister(999); n=5, p=0.4)
    Pwarm = make_problem(gwarm)
    fkw, jnw = identify_fork_and_join_nodes(Pwarm.outgoing, Pwarm.incoming)
    npw = Dict{Int64,pbox}(n => tri_pbox(0.9, 0.05, steps) for n in Pwarm.all_nodes)
    lpw = Dict{Tuple{Int64,Int64},pbox}(e => tri_pbox(0.9, 0.05, steps) for e in keys(Pwarm.eid))
    logmsg(io_log, "warm-up (trivial 5-node network) starting...")
    t_warm = @elapsed begin
        rootsw, uniqw = new_identify(Pwarm.edgelist, npw, lpw, Set{Int64}(Pwarm.sources), fkw, jnw, Pwarm.anc, Pwarm.desc, Pwarm.itersets)
        update_beliefs_iterative(Pwarm.edgelist, Pwarm.itersets, Pwarm.outgoing, Pwarm.incoming, Pwarm.sources, npw, lpw, Pwarm.desc, Pwarm.anc, rootsw, jnw, fkw, uniqw, Dict{CacheKey,DiamondCacheEntry{pbox}}())
    end
    logmsg(io_log, @sprintf("warm-up done: %.4fs (discarded, trivial network -- NOT %s)", t_warm, nm))

    # --- the real target: nm, ONE timed call, ever, in this process ---
    npb = Dict{Int64,pbox}(n => tri_pbox(0.9, 0.05, steps) for n in P.all_nodes)
    lpb = Dict{Tuple{Int64,Int64},pbox}(e => tri_pbox(0.9, 0.05, steps) for e in keys(P.eid))
    logmsg(io_log, "SINGLE timed call starting on $nm...")
    t = @elapsed begin
        roots, uniq = new_identify(P.edgelist, npb, lpb, Set{Int64}(P.sources), fk, jn, P.anc, P.desc, P.itersets)
        update_beliefs_iterative(P.edgelist, P.itersets, P.outgoing, P.incoming, P.sources, npb, lpb, P.desc, P.anc, roots, jn, fk, uniq, Dict{CacheKey,DiamondCacheEntry{pbox}}())
    end
    logmsg(io_log, @sprintf("timed call done: %.4fs", t))

    header_needed = !isfile(OUT_CSV)
    open(OUT_CSV, "a") do f
        header_needed && println(f, "name,V,E,n_diamonds,maxcond,sum_2^C,steps,pbox_propagate_seconds,warmup_seconds_trivial,status")
        @printf(f, "%s,%d,%d,%d,%d,%.0f,%d,%.4f,%.4f,%s\n", nm, length(P.all_nodes), length(g.edges), ndia, maxc, sum2C, steps, t, t_warm, "ok")
    end
    logmsg(io_log, "wrote row for $nm to $OUT_CSV")
    close(io_log)
end

main()
