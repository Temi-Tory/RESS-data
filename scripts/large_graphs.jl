# Larger graphs (n=30..50, single- and multi-source) to preempt "corpus too small". Validate IPA
# (framework new_identify) Float64 vs SIFTED BDD (exact, where tractable) AND vs Monte Carlo (always,
# statistical cross-check). Interval vs sifted corners where the BDD is tractable.
const REPO = raw"c:\Users\ohian\OneDrive - University of Strathclyde\Documents\Programmming Files\Julia Files\InformationPropagation\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(REPO); push!(LOAD_PATH, joinpath(REPO,"validation","bddenv"))
include(joinpath(REPO,"InfoPropFrmwrk","src","Algorithms","InfoPropFramework.jl")); using .InfoPropFramework
using Random, Printf
include(joinpath(REPO,"validation","graph_gen.jl")); include(joinpath(REPO,"validation","oracles.jl"))
include(joinpath(REPO,"validation","graph_families.jl"))
include(joinpath(REPO,"validation","bdd_oracle.jl")); using .BDDReliabilityOracle
include(joinpath(REPO,"validation","oracles_tiered.jl")); using .TieredOracles

function ipa_f64(P,np,lp)
    fk,jn=identify_fork_and_join_nodes(P.outgoing,P.incoming)
    roots,uniq=new_identify(P.edgelist,np,lp,Set{Int64}(P.sources),fk,jn,P.anc,P.desc,P.itersets)
    update_beliefs_iterative(P.edgelist,P.itersets,P.outgoing,P.incoming,P.sources,np,lp,P.desc,P.anc,roots,jn,fk,uniq)
end
function ipa_iv(P,nplo,nphi,lplo,lphi)
    npi=Dict{Int64,Interval}(n=>Interval(nplo[n],nphi[n]) for n in P.all_nodes); lpi=Dict{Tuple{Int64,Int64},Interval}(e=>Interval(lplo[e],lphi[e]) for e in keys(P.eid))
    fk,jn=identify_fork_and_join_nodes(P.outgoing,P.incoming)
    ri,ui=new_identify(P.edgelist,npi,lpi,Set{Int64}(P.sources),fk,jn,P.anc,P.desc,P.itersets)
    update_beliefs_iterative(P.edgelist,P.itersets,P.outgoing,P.incoming,P.sources,npi,lpi,P.desc,P.anc,ri,jn,fk,ui,Dict{CacheKey,DiamondCacheEntry{Interval}}())
end

graphs = Tuple{String,Graph}[]
for n in (30,40,50), p in (0.06,0.10), s in 1:2
    push!(graphs, ("rand_n$(n)_p$(replace(string(p),"."=>""))_s$s", gen_random_dag(MersenneTwister(s); n=n, p=p)))
    push!(graphs, ("multi_n$(n)_p$(replace(string(p),"."=>""))_s$s", gen_multisource(MersenneTwister(100+s); n=n, p=p, nsrc=3)))
end

println("name,V,E,src,ipa_t_ms,bdd_status,f64_vs_bdd,iv_vs_bdd_unsound,iv_vs_bdd_over,f64_vs_mc_stderrs")
worst_bdd=0.0; worst_iv=0.0; worst_mc=0.0; nbdd=0; nmc=0
for (nm,g) in graphs
    P=make_problem(g); edges=collect(P.edgelist); srcs=collect(P.sources)
    rng=MersenneTwister(7); bnp=Dict(n=>0.3+0.69*rand(rng) for n in P.all_nodes); blp=Dict(e=>0.3+0.69*rand(rng) for e in keys(P.eid))
    t=@elapsed belf=ipa_f64(P,bnp,blp)
    # sifted BDD (exact) guarded
    bstat="ok"; fb=-1.0; ivu=-1.0; ivo=-1.0
    try
        bd,_=bdd_reliability(edges, Dict(Int(n)=>bnp[n] for n in P.all_nodes), Dict{Tuple{Int,Int},Float64}(e=>blp[e] for e in keys(P.eid)), srcs; sift=true)
        fb=maximum(abs(belf[n]-bd[n]) for n in P.all_nodes); global nbdd+=1; global worst_bdd=max(worst_bdd,fb)
        w=0.3; nplo=Dict(n=>max(0.0,bnp[n]-w/2) for n in P.all_nodes); nphi=Dict(n=>min(1.0,bnp[n]+w/2) for n in P.all_nodes)
        lplo=Dict(e=>max(0.0,blp[e]-w/2) for e in keys(P.eid)); lphi=Dict(e=>min(1.0,blp[e]+w/2) for e in keys(P.eid))
        beli=ipa_iv(P,nplo,nphi,lplo,lphi)
        bdlo,_=bdd_reliability(edges, Dict(Int(n)=>nplo[n] for n in P.all_nodes), Dict{Tuple{Int,Int},Float64}(e=>lplo[e] for e in keys(P.eid)), srcs; sift=true)
        bdhi,_=bdd_reliability(edges, Dict(Int(n)=>nphi[n] for n in P.all_nodes), Dict{Tuple{Int,Int},Float64}(e=>lphi[e] for e in keys(P.eid)), srcs; sift=true)
        ivu=maximum(max(beli[n].lower-bdlo[n], bdhi[n]-beli[n].upper, 0.0) for n in P.all_nodes)
        ivo=maximum(max(bdlo[n]-beli[n].lower, beli[n].upper-bdhi[n], 0.0) for n in P.all_nodes)
        global worst_iv=max(worst_iv, ivu, ivo)
    catch; bstat="bdd_skip"; end
    # MC cross-check ONLY when the sifted BDD was skipped (exact oracle unavailable). N modest.
    mcstderr = -1.0
    if bstat == "bdd_skip"
        mcb,mce=monte_carlo_reliability(edges, Dict(Int(n)=>bnp[n] for n in P.all_nodes), Dict{Tuple{Int,Int},Float64}(e=>blp[e] for e in keys(P.eid)), srcs; N=200_000, seed=11)
        mcstderr=maximum(abs(belf[n]-mcb[n])/max(mce[n],1e-9) for n in P.all_nodes); global nmc+=1; global worst_mc=max(worst_mc,mcstderr)
    end
    @printf("%s,%d,%d,%d,%.1f,%s,%.2e,%.2e,%.2e,%.2f\n", nm, length(P.all_nodes), length(g.edges), length(P.sources), t*1e3, bstat, fb, ivu, ivo, mcstderr); flush(stdout)
end
@printf("# SUMMARY bdd_checked=%d worst_f64_vs_bdd=%.2e worst_iv=%.2e | mc_checked=%d worst_f64_vs_mc=%.2f stderrs\n",
        nbdd, worst_bdd, worst_iv, nmc, worst_mc)
