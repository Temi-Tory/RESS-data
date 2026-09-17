# FULL-corpus regression against the SIFTED ROBDD oracle (which builds all 131 graphs, unlike the naive
# order that blows up on 17 dense ones). Validates BOTH:
#   Float64  : IPA belief == sifted-BDD belief (exact)
#   Interval : IPA interval == [sifted-BDD(lower corner), sifted-BDD(upper corner)] (exact by monotonicity)
# Uses the FRAMEWORK's new_identify. No V/2M skip (sifted handles the whole corpus).
const REPO = raw"c:\Users\ohian\OneDrive - University of Strathclyde\Documents\Programmming Files\Julia Files\InformationPropagation\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(REPO); push!(LOAD_PATH, joinpath(REPO,"validation","bddenv"))
include(joinpath(REPO,"InfoPropFrmwrk","src","Algorithms","InfoPropFramework.jl")); using .InfoPropFramework
using CUDD, Random, Printf
include(joinpath(REPO,"validation","graph_gen.jl")); include(joinpath(REPO,"validation","oracles.jl"))
include(joinpath(REPO,"validation","bdd_oracle.jl")); using .BDDReliabilityOracle   # sifted oracle

const VCAP = 200   # ultra-safety only; sifted expected to handle all
sources_of(edges)=begin ns=Set{Int}();inc=Set{Int}();for (u,v) in edges;push!(ns,u);push!(ns,v);push!(inc,v);end;sort(collect(setdiff(ns,inc)));end

function ipa_f64(P, np, lp)
    fk,jn=identify_fork_and_join_nodes(P.outgoing,P.incoming)
    roots,uniq=new_identify(P.edgelist,np,lp,Set{Int64}(P.sources),fk,jn,P.anc,P.desc,P.itersets)
    update_beliefs_iterative(P.edgelist,P.itersets,P.outgoing,P.incoming,P.sources,np,lp,P.desc,P.anc,roots,jn,fk,uniq)
end
function ipa_iv(P, nplo,nphi,lplo,lphi)
    np=Dict{Int64,Interval}(n=>Interval(nplo[n],nphi[n]) for n in P.all_nodes); lp=Dict{Tuple{Int64,Int64},Interval}(e=>Interval(lplo[e],lphi[e]) for e in keys(P.eid))
    fk,jn=identify_fork_and_join_nodes(P.outgoing,P.incoming)
    roots,uniq=new_identify(P.edgelist,np,lp,Set{Int64}(P.sources),fk,jn,P.anc,P.desc,P.itersets)
    update_beliefs_iterative(P.edgelist,P.itersets,P.outgoing,P.incoming,P.sources,np,lp,P.desc,P.anc,roots,jn,fk,uniq,Dict{CacheKey,DiamondCacheEntry{Interval}}())
end

graphs = Tuple{String,Graph}[]
push!(graphs,("counterexample-n15", load_edges("cex",joinpath(REPO,"dag_ntwrk_files","counterexample-n15","counterexample-n15.EDGES"))))
let base=gen_random_dag(MersenneTwister(7);n=28,p=0.12); for (nm,g) in scaling_mutants("rand28",base;seeds=1:8,adds=5,dels=2); push!(graphs,(nm,g)); end; end
for n in (10,12,15,20,25), p in (0.10,0.15,0.20), s in 1:8; push!(graphs,("random_n$(n)_p$(replace(string(p),"."=>""))_s$s", gen_random_dag(MersenneTwister(s);n=n,p=p))); end

println("name,V,E,f64_worst,iv_unsound,iv_overwidth")
nchk=0; f64w=0.0; ivu=0.0; ivo=0.0; f64bad=String[]; ivbad=String[]
for (nm,g) in graphs
    P=make_problem(g); P.V>VCAP && continue
    edges=collect(P.edgelist); srcs=collect(P.sources)
    rng=MersenneTwister(1007); bnp=Dict(n=>0.3+0.69*rand(rng) for n in P.all_nodes); blp=Dict(e=>0.3+0.69*rand(rng) for e in keys(P.eid))
    # Float64
    belf=ipa_f64(P,bnp,blp)
    npI=Dict(Int(n)=>bnp[n] for n in P.all_nodes); lpI=Dict{Tuple{Int,Int},Float64}(e=>blp[e] for e in keys(P.eid))
    bd,_=bdd_reliability(edges, npI, lpI, srcs; sift=true)
    fw=maximum(abs(belf[n]-bd[n]) for n in P.all_nodes)
    # Interval
    w=0.3; nplo=Dict(n=>max(0.0,bnp[n]-w/2) for n in P.all_nodes); nphi=Dict(n=>min(1.0,bnp[n]+w/2) for n in P.all_nodes)
    lplo=Dict(e=>max(0.0,blp[e]-w/2) for e in keys(P.eid)); lphi=Dict(e=>min(1.0,blp[e]+w/2) for e in keys(P.eid))
    beli=ipa_iv(P,nplo,nphi,lplo,lphi)
    bdlo,_=bdd_reliability(edges, Dict(Int(n)=>nplo[n] for n in P.all_nodes), Dict{Tuple{Int,Int},Float64}(e=>lplo[e] for e in keys(P.eid)), srcs; sift=true)
    bdhi,_=bdd_reliability(edges, Dict(Int(n)=>nphi[n] for n in P.all_nodes), Dict{Tuple{Int,Int},Float64}(e=>lphi[e] for e in keys(P.eid)), srcs; sift=true)
    un=maximum(max(beli[n].lower-bdlo[n], bdhi[n]-beli[n].upper, 0.0) for n in P.all_nodes)
    ow=maximum(max(bdlo[n]-beli[n].lower, beli[n].upper-bdhi[n], 0.0) for n in P.all_nodes)
    global nchk+=1; global f64w=max(f64w,fw); global ivu=max(ivu,un); global ivo=max(ivo,ow)
    fw>1e-6 && push!(f64bad,nm); (un>1e-6||ow>1e-6) && push!(ivbad,nm)
    @printf("%s,%d,%d,%.2e,%.2e,%.2e\n", nm, length(P.all_nodes), length(g.edges), fw, un, ow)
end
@printf("# SUMMARY checked=%d (all, sifted oracle) | Float64 worst=%.2e wrong=%d | Interval unsound=%.2e overwidth=%.2e wrong=%d\n",
        nchk, f64w, length(f64bad), ivu, ivo, length(ivbad))
!isempty(f64bad) && println("# F64 WRONG: ", f64bad)
!isempty(ivbad) && println("# IV WRONG: ", ivbad)
