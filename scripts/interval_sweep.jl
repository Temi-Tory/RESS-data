# Full-corpus INTERVAL exactness sweep: for each graph, compare IPA interval reliability and NAIVE
# interval propagation (no diamond conditioning) against the EXACT range [CUDD(lower),CUDD(upper)]
# (exact by monotonicity). Reports per-graph: exact max range, IPA over-width (expect ~0 = EXACT),
# naive over-width (the dependency over-widening IPA removes). CSV for the imprecise-propagation artifact.
const REPO = raw"c:\Users\ohian\OneDrive - University of Strathclyde\Documents\Programmming Files\Julia Files\InformationPropagation\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(REPO); push!(LOAD_PATH, joinpath(REPO,"validation","bddenv"))
include(joinpath(REPO,"InfoPropFrmwrk","src","Algorithms","InfoPropFramework.jl")); using .InfoPropFramework
using CUDD, Random, Printf
include(joinpath(REPO,"validation","graph_gen.jl")); include(joinpath(REPO,"validation","oracles.jl")); include(joinpath(REPO,"validation","new_identify.jl"))

const WIDTH = 0.3
function run_iv(P,nplo,nphi,lplo,lphi; naive)
    np=Dict{Int64,Interval}(n=>Interval(nplo[n],nphi[n]) for n in P.all_nodes); lp=Dict{Tuple{Int64,Int64},Interval}(e=>Interval(lplo[e],lphi[e]) for e in keys(P.eid))
    fk,jn=identify_fork_and_join_nodes(P.outgoing,P.incoming)
    roots,uniq = naive ? (Dict{Int64,Vector{DiamondsAtNode}}(), Dict{UInt64,DiamondComputationData{Interval}}()) :
                          new_identify(P.edgelist,np,lp,Set{Int64}(P.sources),fk,jn,P.anc,P.desc,P.itersets)
    update_beliefs_iterative(P.edgelist,P.itersets,P.outgoing,P.incoming,P.sources,np,lp,P.desc,P.anc,roots,jn,fk,uniq,Dict{CacheKey,DiamondCacheEntry{Interval}}())
end
function measure(nm,g)
    P=make_problem(g); P.V>130 && return nothing
    rng=MersenneTwister(1007); bnp=Dict(n=>0.3+0.69*rand(rng) for n in P.all_nodes); blp=Dict(e=>0.3+0.69*rand(rng) for e in keys(P.eid))
    nplo=Dict(n=>max(0.0,bnp[n]-WIDTH/2) for n in P.all_nodes);nphi=Dict(n=>min(1.0,bnp[n]+WIDTH/2) for n in P.all_nodes)
    lplo=Dict(e=>max(0.0,blp[e]-WIDTH/2) for e in keys(P.eid));lphi=Dict(e=>min(1.0,blp[e]+WIDTH/2) for e in keys(P.eid))
    try
        Bc=cudd_build(P); Bc.total_nodes>2_000_000 && (cudd_free(Bc); return nothing)
        wlo=zeros(P.V);whi=zeros(P.V);for n in P.all_nodes;wlo[P.nid[n]]=nplo[n];whi[P.nid[n]]=nphi[n];end;for (e,i) in P.eid;wlo[i]=lplo[e];whi[i]=lphi[e];end
        exlo=cudd_eval(Bc,P,wlo);exhi=cudd_eval(Bc,P,whi);cudd_free(Bc)
        ipa=run_iv(P,nplo,nphi,lplo,lphi;naive=false); nai=run_iv(P,nplo,nphi,lplo,lphi;naive=true)
        exr=maximum(exhi[n]-exlo[n] for n in P.all_nodes)
        ipa_ow=maximum(max(exlo[n]-ipa[n].lower,ipa[n].upper-exhi[n],0.0) for n in P.all_nodes)
        ipa_un=maximum(max(ipa[n].lower-exlo[n],exhi[n]-ipa[n].upper,0.0) for n in P.all_nodes)
        nai_ow=maximum(max(exlo[n]-nai[n].lower,nai[n].upper-exhi[n],0.0) for n in P.all_nodes)
        return (nm,exr,ipa_ow,ipa_un,nai_ow)
    catch; return nothing; end
end

graphs = Tuple{String,Graph}[]
push!(graphs,("counterexample-n15", load_edges("cex",joinpath(REPO,"dag_ntwrk_files","counterexample-n15","counterexample-n15.EDGES"))))
let base=gen_random_dag(MersenneTwister(7);n=28,p=0.12); for (nm,gg) in scaling_mutants("rand28",base;seeds=1:8,adds=5,dels=2); push!(graphs,(nm,gg)); end; end
for n in (10,12,15,20,25), p in (0.10,0.15,0.20), s in 1:8; push!(graphs,("random_n$(n)_p$(replace(string(p),"."=>""))_s$s", gen_random_dag(MersenneTwister(s);n=n,p=p))); end
measure("warm", gen_random_dag(MersenneTwister(1);n=8,p=0.3))  # compile

println("name,exact_range,ipa_overwidth,ipa_unsound,naive_overwidth")
n_ex=0; worst_ipa_ow=0.0; worst_un=0.0; worst_nai=0.0; nchk=0
for (nm,g) in graphs
    r=measure(nm,g); r===nothing && continue
    global nchk+=1; global worst_ipa_ow=max(worst_ipa_ow,r[3]); global worst_un=max(worst_un,r[4]); global worst_nai=max(worst_nai,r[5])
    r[3]<=1e-9 && (global n_ex+=1)
    @printf("%s,%.4f,%.2e,%.2e,%.3e\n", r[1],r[2],r[3],r[4],r[5])
end
@printf("# SUMMARY checked=%d IPA_exact=%d worst_ipa_overwidth=%.2e worst_ipa_unsound=%.2e worst_naive_overwidth=%.3e\n",
        nchk,n_ex,worst_ipa_ow,worst_un,worst_nai)
