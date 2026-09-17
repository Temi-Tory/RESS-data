# RESS response Task 3 -- FINISH JOB: just the missing w=0.10 p-box leg
# (pbox@200, w=0.10, IPA_vs_MC / naive_vs_MC).
#
# The full task3_grid_accuracy_confirm.jl run already produced and wrote to
# data/grid_accuracy.csv (and is NOT re-run here, to avoid re-paying the ~584s
# w=0.05 leg cost):
#   Float64,-,IPA_vs_sifted_BDD
#   Interval,0.05,IPA_vs_BDDcorners / naive_no_conditioning
#   pbox@200,0.05,IPA_vs_MC / naive_vs_MC
#   Interval,0.10,IPA_vs_BDDcorners / naive_no_conditioning
# Those 7 rows are correct/final and are NOT touched by this script.
#
# This script reproduces IDENTICAL setup/config/methodology to
# task3_grid_accuracy_confirm.jl (same graph, same input construction, same
# BDDjl oracle substitution for CUDD, same reused-BDD MC speed fix, same
# steps=200, MC_N=50000, seed=42, target=16) but runs ONLY the w=0.10 pbox
# leg, and APPENDS its two rows to the existing CSV rather than truncating it.
#
# No time budget / early cutoff is applied here -- this is allowed to run to
# genuine completion (expected ~8-10 minutes based on the w=0.05 leg's
# measured 531s propagate + 53s MC time).
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
log_path = joinpath(OUTDIR, "task3_grid_accuracy_w010_pbox_only.log")
io_log = open(log_path, "w")
function logmsg(s...)
    msg = string("[", now(), "] ", s...)
    println(msg); println(io_log, msg); flush(stdout); flush(io_log)
end

const ACC_STEPS = 200
const MC_N = 50_000
const TARGET = 16
const W = 0.10
quiet(f) = redirect_stdout(f, devnull)
glo(c) = hasproperty(c,:lo) ? glo(c.lo) : Float64(c)
ghi(c) = hasproperty(c,:hi) ? ghi(c.hi) : Float64(c)

logmsg("=== Task 3 FINISH: pbox@200 w=0.10 leg only (no time budget this run) ===")
logmsg("PBOX_COND_BLEND default = ", InfoPropFramework.InputProcessingModule.PBOX_COND_BLEND[])

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

# --- pure-Julia BDD oracle (BDDjl), replacing CUDD for stability (same as full script) ---
Bbdd = bddjl_build(P)

csv_path = joinpath(OUTDIR, "grid_accuracy.csv")
@assert isfile(csv_path) "expected existing grid_accuracy.csv with the 7 already-completed rows"

PBA.setSteps(ACC_STEPS)
logmsg("Starting propagate(pbox, w=$W) IPA (sound cvxP) ... this is the long pole, expect several minutes")
t_ipa = @elapsed (ipa = quiet(()->propagate(pbox,W)))
logmsg(@sprintf("  IPA propagate(pbox,%.2f) done in %.2fs", W, t_ipa))
t_nap = @elapsed (nap = quiet(()->propagate(pbox,W; naive=true)))
logmsg(@sprintf("  naive propagate(pbox,%.2f) done in %.2fs", W, t_nap))

rng=MersenneTwister(42); samp=Float64[]; dist=TriangularDist(max(0.0,0.9-W),min(1.0,0.9+W),0.9)
wv = zeros(Float64, P.V)
for n in P.all_nodes; wv[P.nid[n]] = 1.0; end
logmsg("Starting MC (N=$MC_N) via reused BDDjl oracle ...")
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
    global iu = max(iu, glo(ci)-e, e-ghi(ci))
    global nu = max(nu, glo(cn)-e, e-ghi(cn))
end
iu = max(iu,0.0); nu = max(nu,0.0)
note_ipa = iu < 0.02 ? "sound" : "CHECK"
logmsg(@sprintf("ACC pbox@%d w=%.2f  IPA_unsound=%.3e (%s) | naive_unsound=%.3e (node %d, MC_N=%d)", ACC_STEPS, W, iu, note_ipa, nu, TARGET, MC_N))

open(csv_path, "a") do io
    @printf(io,"pbox@%d,%.2f,IPA_vs_MC,-,%.3e,%s\n", ACC_STEPS, W, iu, note_ipa)
    @printf(io,"pbox@%d,%.2f,naive_vs_MC,-,%.3e,unsound\n", ACC_STEPS, W, nu)
end
logmsg("appended 2 rows to ", csv_path)
cp(csv_path, joinpath(RESSDIR, "grid_accuracy.csv"); force=true)
logmsg("=== Task 3 FINISH done ===")
close(io_log)
cp(log_path, joinpath(RESSDIR, "task3_grid_accuracy_w010_pbox_only.log"); force=true)
