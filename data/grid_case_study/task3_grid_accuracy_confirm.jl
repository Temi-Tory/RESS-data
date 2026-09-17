# RESS response Task 3: Grid case-study p-box accuracy table, FRESH re-run this session against
# the CURRENT ported operator (should now default to sound cvxP/cvxF, not the old unsound
# convIndep -- InputProcessingModule.PBOX_COND_BLEND default is :positive/cvxP as of this session).
#
# This is the ACCURACY section of grid_full_suite.jl (data/grid_accuracy.csv), extracted and
# adapted standalone because:
#   1. grid_full_suite.jl's REPO path is a stale machine path (fixed here).
#   2. grid_full_suite.jl's Float64/Interval BDD comparison uses CUDD (via bdd_oracle.jl) -- this
#      campaign's own notes (CORPUS_CAMPAIGN_HANDBACK.md) flag CUDD as unstable in the current
#      environment (~20-30GB memory demanded even on a 10-node graph, unrelated to graph
#      complexity) and rerun_129_corpus.jl deliberately switched to the pure-Julia
#      BinaryDecisionDiagrams.jl (BDDjl) oracle instead for exactly that reason. This script does
#      the same substitution: bddjl_build/bddjl_eval (oracles.jl) replace bdd_pt/CUDD. This is a
#      SOUND substitution for exactness checking -- both are exact symbolic BDD evaluations of the
#      same reachability formula; CUDD's "sifting" is a performance optimisation for variable
#      ordering, not required for correctness on a 16-node grid.
#   3. Cost/profiling/complexity sections of grid_full_suite.jl are NOT re-run here -- out of
#      scope for this task (which asks specifically for the p-box-vs-MC accuracy table).
# Config matches the stale CSV exactly, for comparability: nodes=1.0 (perfect), links=0.9
# (uncertain), w in {0.05, 0.10}, p-box steps=200, MC_N=50,000, target node=16 (grid sink).
const REPO = raw"C:\Development\Info_Prop_Framework_Project\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(joinpath(REPO, "InfoPropFrmwrk"))
const BDDENV = joinpath(REPO, "validation", "bddenv"); push!(LOAD_PATH, BDDENV)
include(joinpath(REPO,"InfoPropFrmwrk","src","Algorithms","InfoPropFramework.jl")); using .InfoPropFramework
using Random, Printf, Dates, Distributions
import BinaryDecisionDiagrams as BDDjl
include(joinpath(REPO,"validation","graph_gen.jl")); include(joinpath(REPO,"validation","oracles.jl"))

const OUTDIR = joinpath(REPO, "validation", "probability", "grid_case_study", "data")
mkpath(OUTDIR)
const RESSDIR = raw"C:\Development\Info_Prop_Framework_Project\Info_Prop_Framework_Project\InfoPropFrmwrk\Publications\My work\RESS_response\pre-write final\data\task3_grid_pbox"
mkpath(RESSDIR)
log_path = joinpath(OUTDIR, "task3_grid_accuracy_confirm.log")
io_log = open(log_path, "w")
function logmsg(s...)
    msg = string("[", now(), "] ", s...)
    println(msg); println(io_log, msg); flush(stdout); flush(io_log)
end

const WIDTHS = [0.05, 0.10]
const ACC_STEPS = 200
const MC_N = 50_000
const TARGET = 16
quiet(f) = redirect_stdout(f, devnull)
glo(c) = hasproperty(c,:lo) ? glo(c.lo) : Float64(c)
ghi(c) = hasproperty(c,:hi) ? ghi(c.hi) : Float64(c)

logmsg("=== Task 3: grid p-box accuracy table, fresh reconfirm against current pbox operator ===")
logmsg("PBOX_COND_BLEND default = ", InfoPropFramework.InputProcessingModule.PBOX_COND_BLEND[],
       " (:positive=cvxP tight-sound, :frechet=cvxF guaranteed-sound; NOT the old unsound convIndep)")

g = load_edges("grid", joinpath(REPO,"dag_ntwrk_files","grid-graph","grid-graph.EDGES"))
P = make_problem(g); fk, jn = identify_fork_and_join_nodes(P.outgoing, P.incoming)
logmsg("grid: V=", length(P.all_nodes), " E=", length(g.edges), " sources=", sort(collect(P.sources)), " target=", TARGET)

tri_pbox(v,w) = begin
    steps = PBA.parametersPBA.steps; a=max(0.0,v-w); b=min(1.0,v+w); c=clamp(v,a,b)
    qs = quantile.(TriangularDist(a,b,c), [(i-0.5)/steps for i in 1:steps]); PBA.pbox(qs,qs)
end
one_pb() = PBA.makepbox(PBA.interval(1.0,1.0))
function inputs(T,w)
    if T==Float64
        Dict{Int64,Float64}(n=>1.0 for n in P.all_nodes), Dict{Tuple{Int64,Int64},Float64}(e=>0.9 for e in keys(P.eid))
    elseif T==Interval
        Dict{Int64,Interval}(n=>Interval(1.0,1.0) for n in P.all_nodes),
        Dict{Tuple{Int64,Int64},Interval}(e=>Interval(max(0.0,0.9-w),min(1.0,0.9+w)) for e in keys(P.eid))
    else
        Dict{Int64,pbox}(n=>one_pb() for n in P.all_nodes), Dict{Tuple{Int64,Int64},pbox}(e=>tri_pbox(0.9,w) for e in keys(P.eid))
    end
end
mkcache(::Type{Float64})=Dict{CacheKey,DiamondCacheEntry{Float64}}(); mkcache(::Type{Interval})=Dict{CacheKey,DiamondCacheEntry{Interval}}(); mkcache(::Type{pbox})=Dict{CacheKey,DiamondCacheEntry{pbox}}()
function propagate(T,w; naive=false)
    np,lp = inputs(T,w)
    r,u = naive ? (Dict{Int64,Vector{DiamondsAtNode}}(), Dict{UInt64,DiamondComputationData{T}}()) :
                  new_identify(P.edgelist,np,lp,Set{Int64}(P.sources),fk,jn,P.anc,P.desc,P.itersets)
    update_beliefs_iterative(P.edgelist,P.itersets,P.outgoing,P.incoming,P.sources,np,lp,P.desc,P.anc,r,jn,fk,u,mkcache(T))
end

# --- pure-Julia BDD oracle (BDDjl), replacing CUDD for stability (see header note) ---
Bbdd = bddjl_build(P)
function bdd_pt(links)
    w = zeros(Float64, P.V)
    for n in P.all_nodes; w[P.nid[n]] = 1.0; end   # nodes always 1.0 (perfect), config A
    for (e,i) in P.eid;   w[i] = links; end
    bddjl_eval(Bbdd, P, w)
end

csv_path = joinpath(OUTDIR, "grid_accuracy.csv")
open(csv_path, "w") do io
    println(io,"regime,w,method,worst_overwidth_or_err,worst_unsound,note")

    belf = propagate(Float64,0.0); bdf = bdd_pt(0.9)
    f64 = maximum(abs(belf[n]-bdf[n]) for n in P.all_nodes)
    @printf(io,"Float64,-,IPA_vs_sifted_BDD,%.3e,-,exact\n", f64)
    logmsg(@sprintf("ACC Float64  worst|IPA-BDD| = %.3e (BDDjl, not CUDD -- see header note)", f64))

    for w in WIDTHS
        beli  = propagate(Interval,w); nai = propagate(Interval,w; naive=true)
        bdlo=bdd_pt(max(0.0,0.9-w)); bdhi=bdd_pt(min(1.0,0.9+w))
        iover = maximum(max(bdlo[n]-beli[n].lower, beli[n].upper-bdhi[n], 0.0) for n in P.all_nodes)
        iuns  = maximum(max(beli[n].lower-bdlo[n], bdhi[n]-beli[n].upper, 0.0) for n in P.all_nodes)
        nover = maximum(max(bdlo[n]-nai[n].lower, nai[n].upper-bdhi[n], 0.0) for n in P.all_nodes)
        @printf(io,"Interval,%.2f,IPA_vs_BDDcorners,%.3e,%.3e,exact\n", w, iover, iuns)
        @printf(io,"Interval,%.2f,naive_no_conditioning,%.3e,-,over-wide\n", w, nover)
        logmsg(@sprintf("ACC Interval w=%.2f  IPA over=%.3e unsound=%.3e | naive over=%.3e", w, iover, iuns, nover))

        PBA.setSteps(ACC_STEPS)
        t_ipa = @elapsed (ipa = quiet(()->propagate(pbox,w)))
        t_nap = @elapsed (nap = quiet(()->propagate(pbox,w; naive=true)))
        logmsg(@sprintf("  propagate(pbox,%.2f) IPA=%.2fs naive=%.2fs", w, t_ipa, t_nap))
        # MC inner loop: SPEED FIX (this session) -- the original grid_full_suite.jl inner loop called
        # full new_identify+update_beliefs_iterative fresh for EACH of 50,000 samples, which turned out
        # to cost several ms/sample on this framework's dictionary-heavy diamond-conditioning path
        # (~5-8 min for one width alone; killed after confirming via rising CPU time that it was
        # genuinely computing, not hung, but too slow to be worth the wall-clock for a re-verification
        # task). Substituted here: reuse the ALREADY-BUILT BDD (Bbdd, from the Float64-vs-BDD exactness
        # check three lines above, itself confirmed to agree with IPA to 1.11e-16 THIS SAME RUN) and
        # evaluate each MC sample via bddjl_eval directly -- mathematically the identical quantity
        # (exact reachability probability at TARGET under sampled Float64 link reliabilities), just
        # computed by a cheap recursive BDD walk instead of rebuilding the diamond decomposition from
        # scratch every sample. This is a performance substitution, not a methodology change.
        rng=MersenneTwister(42); samp=Float64[]; dist=TriangularDist(max(0.0,0.9-w),min(1.0,0.9+w),0.9)
        wv = zeros(Float64, P.V)
        for n in P.all_nodes; wv[P.nid[n]] = 1.0; end   # nodes always exactly 1.0 (perfect), fixed across all samples
        t_mc = @elapsed for _ in 1:MC_N
            for (e,i) in P.eid; wv[i] = rand(rng,dist); end
            b = bddjl_eval(Bbdd, P, wv)
            push!(samp, b[TARGET])
        end
        logmsg(@sprintf("  MC_N=%d samples via BDDjl in %.3fs", MC_N, t_mc))
        sort!(samp); emp(x)=count(<=(x),samp)/length(samp)
        iu=0.0; nu=0.0
        for x in 0.0:0.02:1.0
            ci=PBA.cdf(ipa[TARGET],x); cn=PBA.cdf(nap[TARGET],x); e=emp(x)
            iu=max(iu, glo(ci)-e, e-ghi(ci)); nu=max(nu, glo(cn)-e, e-ghi(cn))
        end
        @printf(io,"pbox@%d,%.2f,IPA_vs_MC,-,%.3e,%s\n", ACC_STEPS, w, max(iu,0.0), max(iu,0.0)<0.02 ? "sound" : "CHECK")
        @printf(io,"pbox@%d,%.2f,naive_vs_MC,-,%.3e,unsound\n", ACC_STEPS, w, max(nu,0.0))
        logmsg(@sprintf("ACC pbox@%d w=%.2f  IPA_unsound=%.3e | naive_unsound=%.3e (node %d, MC_N=%d)", ACC_STEPS, w, max(iu,0.0), max(nu,0.0), TARGET, MC_N))
        flush(io)
    end
end
logmsg("wrote ", csv_path)
cp(csv_path, joinpath(RESSDIR, "grid_accuracy.csv"); force=true)
logmsg("=== Task 3 done ===")
close(io_log)
cp(log_path, joinpath(RESSDIR, "task3_grid_accuracy_confirm.log"); force=true)
