# Validate the convex-combination p-box operator (cvxI / cvxF) across a corpus vs Monte Carlo.
# Standalone (rc_core recursion + PBA). For each graph & sink target, report soundness + band for
# BOTH regimes: perfect nodes (1.0, links uncertain -> stresses belief near 1) and uncertain nodes.
# cvxF must stay SOUND everywhere to claim a rigorous analytic p-box contribution.
const REPO = raw"c:\Users\ohian\OneDrive - University of Strathclyde\Documents\Programmming Files\Julia Files\InformationPropagation\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(REPO)
using ProbabilityBoundsAnalysis, Distributions, Random, Printf
const PBA = ProbabilityBoundsAnalysis
include(joinpath(REPO,"validation","rc_core.jl")); using .RCCore
include(joinpath(REPO,"validation","graph_gen.jl")); include(joinpath(REPO,"validation","graph_families.jl"))
PBA.setSteps(50)
glo(c)=hasproperty(c,:lo) ? glo(c.lo) : Float64(c); ghi(c)=hasproperty(c,:hi) ? ghi(c.hi) : Float64(c)
ONE=PBA.makepbox(PBA.interval(1.0,1.0)); ZERO=PBA.makepbox(PBA.interval(0.0,0.0))
mul(a,b)=PBA.convIndep(a,b,op=*); comp(a)=PBA.convIndep(ONE,a,op=-)
function cvxcomb(W,A,B,blend)
    n=length(W.u); ps=fill(1.0/n,n)
    Mu=PBA.mixture([blend(W.u[i]*A,(1.0-W.u[i])*B,op=+) for i in 1:n], ps)
    Md=PBA.mixture([blend(W.d[i]*A,(1.0-W.d[i])*B,op=+) for i in 1:n], ps)
    PBA.env(Mu,Md)
end
function rc_pbox(edgelist,np,lp,sources,mode)
    out=Dict{Int,Vector{Int}}();inc=Dict{Int,Vector{Int}}();nodes=Set{Int}()
    for (u,v) in edgelist; push!(get!(out,u,Int[]),v);push!(get!(inc,v,Int[]),u);push!(nodes,u);push!(nodes,v);end
    nodes=sort(collect(nodes));srcset=Set(sources)
    indeg=Dict(v=>length(get(inc,v,Int[])) for v in nodes);q=sort([v for v in nodes if indeg[v]==0]);topo=Int[]
    while !isempty(q);u=popfirst!(q);push!(topo,u);for w in get(out,u,Int[]);indeg[w]-=1;indeg[w]==0&&push!(q,w);end;end
    topi=Dict(v=>i for (i,v) in enumerate(topo));forks=Set(v for v in nodes if length(get(out,v,Int[]))>=2)
    isone(p)=(glo(PBA.left(p))>=1-1e-9);iszero(p)=(ghi(PBA.right(p))<=1e-9)
    det=Set(v for v in nodes if (haskey(np,v)&&iszero(np[v]))||(haskey(np,v)&&isone(np[v])&&v in srcset))
    anc=Dict{Int,Set{Int}}();for v in topo;s=Set{Int}();for u in get(inc,v,Int[]);push!(s,u);union!(s,anc[u]);end;anc[v]=s;end
    function pick(v,cond)
        ps=get(inc,v,Int[]);length(ps)<2&&return nothing;cover=Dict{Int,Int}()
        for p in ps,a in anc[p];(a in forks)&&!(a in keys(cond))&&!(a in det)&&(cover[a]=get(cover,a,0)+1);end
        for p in ps,qq in ps;(p==qq)&&continue;(p in anc[qq])&&!(p in keys(cond))&&!(p in det)&&(cover[p]=get(cover,p,0)+1);end
        best=nothing;bk=(typemax(Int),);for (a,c) in cover;c<2&&continue;k=(topi[a],);(k<bk)&&(bk=k;best=a);end;best
    end
    memo=Dict{Tuple{Int,Vector{Tuple{Int,Bool}}},pbox}()
    function R(v,cond)
        haskey(cond,v)&&return cond[v] ? ONE : ZERO
        v in srcset&&return np[v]
        av=anc[v];key=(v,sort([(c,cond[c]) for c in keys(cond) if c in av]));haskey(memo,key)&&return memo[key]
        ps=get(inc,v,Int[]);f=pick(v,cond);local val
        if f===nothing
            surv=ONE;for u in ps;surv=mul(surv,comp(mul(R(u,cond),lp[(u,v)])));end;val=mul(np[v],comp(surv))
        else
            rf=R(f,cond);cu=copy(cond);cu[f]=true;cd=copy(cond);cd[f]=false;bup=R(v,cu);bdn=R(v,cd)
            val = mode==:cvxI ? cvxcomb(rf,bup,bdn,PBA.convIndep) : cvxcomb(rf,bup,bdn,PBA.convFrechet)
        end
        memo[key]=val;val
    end
    Dict(v=>(v in srcset ? np[v] : R(v,Dict{Int,Bool}())) for v in nodes)
end

edgesof(g)=[(Int(u),Int(v)) for (u,v) in g.edges]
srcof(E)=begin dst=Set(v for (_,v) in E); sort([u for u in Set(vcat([e[1] for e in E],[e[2] for e in E])) if !(u in dst)]) end
tri(v,w)=(a=max(0.0,v-w);b=min(1.0,v+w);qs=quantile.(TriangularDist(a,b,clamp(v,a,b)),[(i-0.5)/50 for i in 1:50]);PBA.pbox(qs,qs))

function test(nm, E, tgt, perfect; w=0.05, N=6000)
    alln=sort(collect(Set(vcat([e[1] for e in E],[e[2] for e in E])))); srcs=srcof(E)
    np = perfect ? Dict{Int,pbox}(n=>ONE for n in alln) : Dict{Int,pbox}(n=>tri(0.7,w) for n in alln)
    lp = Dict{Tuple{Int,Int},pbox}(e=>tri(0.9,w) for e in E)
    ri=rc_pbox(E,np,lp,srcs,:cvxI); rf=rc_pbox(E,np,lp,srcs,:cvxF)
    rng=MersenneTwister(7);samp=Float64[]
    dn = perfect ? nothing : TriangularDist(0.7-w,0.7+w,0.7); dl=TriangularDist(0.9-w,0.9+w,0.9)
    for _ in 1:N
        npf=Dict{Int,Float64}(n=> dn===nothing ? 1.0 : rand(rng,dn) for n in alln); lpf=Dict{Tuple{Int,Int},Float64}(e=>rand(rng,dl) for e in E)
        push!(samp,RCCore.rc_reliability(E,npf,lpf,srcs)[tgt])
    end
    sort!(samp);emp(x)=count(<=(x),samp)/length(samp)
    uns(pb)=begin m=0.0;for x in 0.0:0.02:1.0;c=PBA.cdf(pb[tgt],x);m=max(m,glo(c)-emp(x),emp(x)-ghi(c));end;max(m,0.0) end
    bnd(pb)=begin m=0.0;for x in 0.0:0.02:1.0;c=PBA.cdf(pb[tgt],x);m=max(m,ghi(c)-glo(c));end;m end
    @printf("%-22s %-9s tgt=%2d  cvxI: uns=%.3f bnd=%.2f | cvxF: uns=%.3f bnd=%.2f  %s\n",
            nm, perfect ? "perfect" : "uncertain", tgt, uns(ri),bnd(ri), uns(rf),bnd(rf),
            (uns(rf)<0.03) ? "cvxF SOUND" : "cvxF UNSOUND!")
end

graphs = [ ("bridge_3", edgesof(gen_bridge(3))), ("bridge_5", edgesof(gen_bridge(5))),
           ("grid_3x4", edgesof(gen_grid(3,4))), ("grid_4x4", edgesof(gen_grid(4,4))),
           ("seriespar_3", edgesof(gen_series_parallel(MersenneTwister(1);blocks=3))),
           ("layered_4x3", edgesof(gen_layered(MersenneTwister(1);layers=4,width=3,p=0.6))),
           ("random_n12", edgesof(gen_random_dag(MersenneTwister(1);n=12,p=0.2))),
           ("random_n15", edgesof(gen_random_dag(MersenneTwister(3);n=15,p=0.2))) ]
println("# convex-combination p-box operator vs MC (steps=50). cvxF must stay SOUND (uns<0.03).")
for (nm,E) in graphs
    tgt=maximum(vcat([e[1] for e in E],[e[2] for e in E]))
    for perf in (true,false); test(nm,E,tgt,perf); end
end
println("# done")
