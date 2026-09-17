# Validate the additional DAG families (multi-source, grid, layered, bridge, series-parallel, complete)
# against the SIFTED ROBDD oracle, Float64 + Interval, using the FRAMEWORK's new_identify.
const REPO = raw"c:\Users\ohian\OneDrive - University of Strathclyde\Documents\Programmming Files\Julia Files\InformationPropagation\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(REPO); push!(LOAD_PATH, joinpath(REPO,"validation","bddenv"))
include(joinpath(REPO,"InfoPropFrmwrk","src","Algorithms","InfoPropFramework.jl")); using .InfoPropFramework
using CUDD, Random, Printf
include(joinpath(REPO,"validation","graph_gen.jl")); include(joinpath(REPO,"validation","oracles.jl"))
include(joinpath(REPO,"validation","graph_families.jl"))
include(joinpath(REPO,"validation","bdd_oracle.jl")); using .BDDReliabilityOracle

function check(nm, g)
    P = make_problem(g); edges=collect(P.edgelist); srcs=collect(P.sources)
    rng=MersenneTwister(1007); bnp=Dict(n=>0.3+0.69*rand(rng) for n in P.all_nodes); blp=Dict(e=>0.3+0.69*rand(rng) for e in keys(P.eid))
    fk,jn=identify_fork_and_join_nodes(P.outgoing,P.incoming)
    # Float64
    roots,uniq=new_identify(P.edgelist,bnp,blp,Set{Int64}(P.sources),fk,jn,P.anc,P.desc,P.itersets)
    belf=update_beliefs_iterative(P.edgelist,P.itersets,P.outgoing,P.incoming,P.sources,bnp,blp,P.desc,P.anc,roots,jn,fk,uniq)
    bd,_=bdd_reliability(edges, Dict(Int(n)=>bnp[n] for n in P.all_nodes), Dict{Tuple{Int,Int},Float64}(e=>blp[e] for e in keys(P.eid)), srcs; sift=true)
    fw=maximum(abs(belf[n]-bd[n]) for n in P.all_nodes)
    # Interval
    w=0.3; nplo=Dict(n=>max(0.0,bnp[n]-w/2) for n in P.all_nodes); nphi=Dict(n=>min(1.0,bnp[n]+w/2) for n in P.all_nodes)
    lplo=Dict(e=>max(0.0,blp[e]-w/2) for e in keys(P.eid)); lphi=Dict(e=>min(1.0,blp[e]+w/2) for e in keys(P.eid))
    npi=Dict{Int64,Interval}(n=>Interval(nplo[n],nphi[n]) for n in P.all_nodes); lpi=Dict{Tuple{Int64,Int64},Interval}(e=>Interval(lplo[e],lphi[e]) for e in keys(P.eid))
    ri,ui=new_identify(P.edgelist,npi,lpi,Set{Int64}(P.sources),fk,jn,P.anc,P.desc,P.itersets)
    beli=update_beliefs_iterative(P.edgelist,P.itersets,P.outgoing,P.incoming,P.sources,npi,lpi,P.desc,P.anc,ri,jn,fk,ui,Dict{CacheKey,DiamondCacheEntry{Interval}}())
    bdlo,_=bdd_reliability(edges, Dict(Int(n)=>nplo[n] for n in P.all_nodes), Dict{Tuple{Int,Int},Float64}(e=>lplo[e] for e in keys(P.eid)), srcs; sift=true)
    bdhi,_=bdd_reliability(edges, Dict(Int(n)=>nphi[n] for n in P.all_nodes), Dict{Tuple{Int,Int},Float64}(e=>lphi[e] for e in keys(P.eid)), srcs; sift=true)
    iu=maximum(max(beli[n].lower-bdlo[n], bdhi[n]-beli[n].upper, 0.0) for n in P.all_nodes)
    io=maximum(max(bdlo[n]-beli[n].lower, beli[n].upper-bdhi[n], 0.0) for n in P.all_nodes)
    @printf("%-18s V=%2d E=%2d src=%d  f64|d|=%.2e  iv:unsound=%.2e over=%.2e  %s\n",
            nm, length(P.all_nodes), length(g.edges), length(P.sources), fw, iu, io,
            (fw>1e-6 || iu>1e-6 || io>1e-6) ? "WRONG" : "ok")
end

# warm
check("warm", gen_grid(2,3))
println("=== additional DAG families vs sifted BDD (Float64 + Interval) ===")
for (n,nsrc,s) in [(15,2,1),(15,3,2),(20,2,3),(20,3,4)]; check("multisrc_n$(n)_k$(nsrc)_s$s", gen_multisource(MersenneTwister(s); n=n, p=0.15, nsrc=nsrc)); end
for (r,c) in [(3,4),(4,4),(4,5),(5,5)]; check("grid_$(r)x$(c)", gen_grid(r,c)); end
for (L,W,s) in [(4,4,1),(5,4,2),(4,5,3)]; check("layered_$(L)x$(W)_s$s", gen_layered(MersenneTwister(s); layers=L, width=W, p=0.5)); end
for k in [1,2,3,4]; check("bridge_$k", gen_bridge(k)); end
for (b,s) in [(3,1),(5,2)]; check("seriesparallel_$(b)_s$s", gen_series_parallel(MersenneTwister(s); blocks=b)); end
for n in [6,8,10]; check("complete_$n", gen_complete(n)); end
println("# done")
