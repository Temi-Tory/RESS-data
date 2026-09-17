# RESS response Task 3 (continuation): grid p-box accuracy, w=0.10 leg.
# w=0.05 was already confirmed this session (worst_unsound=0.000e+00, sound) and persisted --
# see task3_grid_accuracy_confirm.jl / .log / grid_accuracy.csv. This script computes the
# remaining w=0.10 rows and appends them to the existing grid_accuracy.csv.
# NOTE: the first run's w=0.10 Interval rows were computed and LOGGED (visible in
# task3_grid_accuracy_confirm.log) but never reached the CSV -- the original script only called
# flush(io) once at the END of each width's whole loop body (after the pbox section), so when the
# w=0.10 pbox propagate call was killed mid-flight, the already-buffered-but-unflushed w=0.10
# Interval CSV rows were lost with it, even though they had already been correctly computed and
# printed to the (separately, per-line, flushed) log. This script recomputes them (cheap, exact
# BDD-corner comparison, sub-second) so the final CSV is complete and consistent, then computes
# the w=0.10 pbox rows (the expensive part).
# Per explicit instruction: let this run to ACTUAL completion, however long it takes -- no
# external kill this time.
const REPO = raw"C:\Development\Info_Prop_Framework_Project\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(joinpath(REPO, "InfoPropFrmwrk"))
const BDDENV = joinpath(REPO, "validation", "bddenv"); push!(LOAD_PATH, BDDENV)
include(joinpath(REPO,"InfoPropFrmwrk","src","Algorithms","InfoPropFramework.jl")); using .InfoPropFramework
using Random, Printf, Dates, Distributions
import BinaryDecisionDiagrams as BDDjl
include(joinpath(REPO,"validation","graph_gen.jl")); include(joinpath(REPO,"validation","oracles.jl"))

const OUTDIR = joinpath(REPO, "validation", "probability", "grid_case_study", "data")
const RESSDIR = raw"C:\Development\Info_Prop_Framework_Project\Info_Prop_Framework_Project\InfoPropFrmwrk\Publications\My work\RESS_response\pre-write final\data\task3_grid_pbox"
log_path = joinpath(OUTDIR, "task3_grid_accuracy_w10_only.log")
io_log = open(log_path, "w")
function logmsg(s...)
    msg = string("[", now(), "] ", s...)
    println(msg); println(io_log, msg); flush(stdout); flush(io_log)
end

const W = 0.10
const ACC_STEPS = 200
const MC_N = 50_000
const TARGET = 16
quiet(f) = redirect_stdout(f, devnull)
glo(c) = hasproperty(c,:lo) ? glo(c.lo) : Float64(c)
ghi(c) = hasproperty(c,:hi) ? ghi(c.hi) : Float64(c)

logmsg("=== Task 3 continuation: grid p-box accuracy, w=0.10 leg ONLY (running to full completion, no kill) ===")
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
Bbdd = bddjl_build(P)
function bdd_pt(links)
    w = zeros(Float64, P.V)
    for n in P.all_nodes; w[P.nid[n]] = 1.0; end
    for (e,i) in P.eid;   w[i] = links; end
    bddjl_eval(Bbdd, P, w)
end

# --- recompute + append the w=0.10 Interval rows (lost from the first run, see header note) ---
csv_path_early = joinpath(OUTDIR, "grid_accuracy.csv")
beli  = propagate(Interval,W); nai = propagate(Interval,W; naive=true)
bdlo=bdd_pt(max(0.0,0.9-W)); bdhi=bdd_pt(min(1.0,0.9+W))
iover = maximum(max(bdlo[n]-beli[n].lower, beli[n].upper-bdhi[n], 0.0) for n in P.all_nodes)
iuns  = maximum(max(beli[n].lower-bdlo[n], bdhi[n]-beli[n].upper, 0.0) for n in P.all_nodes)
nover = maximum(max(bdlo[n]-nai[n].lower, nai[n].upper-bdhi[n], 0.0) for n in P.all_nodes)
open(csv_path_early, "a") do io
    @printf(io,"Interval,%.2f,IPA_vs_BDDcorners,%.3e,%.3e,exact\n", W, iover, iuns)
    @printf(io,"Interval,%.2f,naive_no_conditioning,%.3e,-,over-wide\n", W, nover)
end
logmsg(@sprintf("ACC Interval w=%.2f (recomputed+flushed)  IPA over=%.3e unsound=%.3e | naive over=%.3e", W, iover, iuns, nover))

PBA.setSteps(ACC_STEPS)
logmsg("starting propagate(pbox, w=", W, ") IPA (sound cvxP operator) -- this is the long pole, expect several minutes...")
t_ipa = @elapsed (ipa = quiet(()->propagate(pbox,W)))
logmsg(@sprintf("  IPA propagate DONE in %.2fs", t_ipa))
t_nap = @elapsed (nap = quiet(()->propagate(pbox,W; naive=true)))
logmsg(@sprintf("  naive propagate DONE in %.2fs", t_nap))

rng=MersenneTwister(42); samp=Float64[]; dist=TriangularDist(max(0.0,0.9-W),min(1.0,0.9+W),0.9)
wv = zeros(Float64, P.V)
for n in P.all_nodes; wv[P.nid[n]] = 1.0; end
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
    global iu = max(iu, glo(ci)-e, e-ghi(ci)); global nu = max(nu, glo(cn)-e, e-ghi(cn))
end
iu = max(iu,0.0); nu = max(nu,0.0)
note = iu < 0.02 ? "sound" : "CHECK"
logmsg(@sprintf("ACC pbox@%d w=%.2f  IPA_unsound=%.3e | naive_unsound=%.3e (node %d, MC_N=%d)", ACC_STEPS, W, iu, nu, TARGET, MC_N))

# --- merge into the existing grid_accuracy.csv (append the two new rows) ---
csv_path = joinpath(OUTDIR, "grid_accuracy.csv")
open(csv_path, "a") do io
    @printf(io,"pbox@%d,%.2f,IPA_vs_MC,-,%.3e,%s\n", ACC_STEPS, W, iu, note)
    @printf(io,"pbox@%d,%.2f,naive_vs_MC,-,%.3e,unsound\n", ACC_STEPS, W, nu)
end
logmsg("appended w=0.10 pbox rows to ", csv_path)
cp(csv_path, joinpath(RESSDIR, "grid_accuracy.csv"); force=true)
logmsg("=== Task 3 w=0.10 leg done (FULL COMPLETION, not killed) ===")
close(io_log)
cp(log_path, joinpath(RESSDIR, "task3_grid_accuracy_w10_only.log"); force=true)
