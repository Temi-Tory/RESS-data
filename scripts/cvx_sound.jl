# BROAD SOUNDNESS sweep of cvxP (positive-dependence convex-combination operator). cvxP ONLY, steps=20
# (soundness ~step-independent), MC N=3000. Goal: is cvxP SOUND (uns<0.03) across ALL topologies+regimes?
# Flush per config. This is the "be sure" evidence for the analytic p-box soundness claim.
const REPO = raw"c:\Users\ohian\OneDrive - University of Strathclyde\Documents\Programmming Files\Julia Files\InformationPropagation\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(REPO)
using ProbabilityBoundsAnalysis, Distributions, Random, Printf
const PBA = ProbabilityBoundsAnalysis
include(joinpath(REPO,"validation","rc_core.jl")); using .RCCore
include(joinpath(REPO,"validation","graph_gen.jl")); include(joinpath(REPO,"validation","graph_families.jl"))
PBA.setSteps(20)
glo(c)=hasproperty(c,:lo) ? glo(c.lo) : Float64(c); ghi(c)=hasproperty(c,:hi) ? ghi(c.hi) : Float64(c)
ONE=PBA.makepbox(PBA.interval(1.0,1.0)); ZERO=PBA.makepbox(PBA.interval(0.0,0.0))
mul(a,b)=PBA.convIndep(a,b,op=*); comp(a)=PBA.convIndep(ONE,a,op=-)
bP(x,y;op=+)=PBA.env(PBA.convIndep(x,y,op=op), PBA.convPerfect(x,y,op=op))
function cvxcomb(W,A,B)
    n=length(W.u); ps=fill(1.0/n,n)
    Mu=PBA.mixture([bP(W.u[i]*A,(1.0-W.u[i])*B) for i in 1:n], ps)
    Md=PBA.mixture([bP(W.d[i]*A,(1.0-W.d[i])*B) for i in 1:n], ps)
    PBA.env(Mu,Md)
end
function rc_pbox(edgelist,np,lp,sources)
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
            rf=R(f,cond);cu=copy(cond);cu[f]=true;cd=copy(cond);cd[f]=false;val=cvxcomb(rf,R(v,cu),R(v,cd))
        end
        memo[key]=val;val
    end
    Dict(v=>(v in srcset ? np[v] : R(v,Dict{Int,Bool}())) for v in nodes)
end
edgesof(g)=[(Int(u),Int(v)) for (u,v) in g.edges]
srcof(E)=begin dst=Set(v for (_,v) in E); sort([u for u in Set(vcat([e[1] for e in E],[e[2] for e in E])) if !(u in dst)]) end
tri(v,w)=(a=max(0.0,v-w);b=min(1.0,v+w);qs=quantile.(TriangularDist(a,b,clamp(v,a,b)),[(i-0.5)/20 for i in 1:20]);PBA.pbox(qs,qs))
function test(nm,E,tgt,nodeval;linkval=0.9,w=0.05,N=3000)
    alln=sort(collect(Set(vcat([e[1] for e in E],[e[2] for e in E])))); srcs=srcof(E)
    np = nodeval==1.0 ? Dict{Int,pbox}(n=>ONE for n in alln) : Dict{Int,pbox}(n=>tri(nodeval,w) for n in alln)
    lp = Dict{Tuple{Int,Int},pbox}(e=>tri(linkval,w) for e in E)
    rp=rc_pbox(E,np,lp,srcs)
    rng=MersenneTwister(7);samp=Float64[]
    dn = nodeval==1.0 ? nothing : TriangularDist(max(0.,nodeval-w),min(1.,nodeval+w),nodeval); dl=TriangularDist(max(0.,linkval-w),min(1.,linkval+w),linkval)
    for _ in 1:N
        npf=Dict{Int,Float64}(n=> dn===nothing ? 1.0 : rand(rng,dn) for n in alln); lpf=Dict{Tuple{Int,Int},Float64}(e=>rand(rng,dl) for e in E)
        push!(samp,RCCore.rc_reliability(E,npf,lpf,srcs)[tgt])
    end
    sort!(samp);emp(x)=count(<=(x),samp)/length(samp)
    u=0.0;b=0.0;for x in 0.0:0.02:1.0;c=PBA.cdf(rp[tgt],x);u=max(u,glo(c)-emp(x),emp(x)-ghi(c));b=max(b,ghi(c)-glo(c));end
    @printf("%-24s uns=%.3f band=%.2f  %s\n", nm, max(u,0.0), b, max(u,0.0)<0.03 ? "SOUND" : "UNSOUND!!"); flush(stdout)
end
println("# cvxP broad soundness sweep (steps=20). All must be SOUND.")
graphs=[("bridge_3",gen_bridge(3)),("bridge_5",gen_bridge(5)),("grid_3x4",gen_grid(3,4)),("grid_4x4",gen_grid(4,4)),
        ("seriespar_3",gen_series_parallel(MersenneTwister(1);blocks=3)),("layered_4x3",gen_layered(MersenneTwister(1);layers=4,width=3,p=0.6)),
        ("multisrc_n12",gen_multisource(MersenneTwister(1);n=12,p=0.2,nsrc=3)),("random_n12",gen_random_dag(MersenneTwister(1);n=12,p=0.2)),
        ("random_n15",gen_random_dag(MersenneTwister(3);n=15,p=0.2)),("random_n12b",gen_random_dag(MersenneTwister(6);n=12,p=0.25))]
for (nm,g) in graphs
    E=edgesof(g); tgt=maximum(vcat([e[1] for e in E],[e[2] for e in E]))
    test("$nm perf",E,tgt,1.0); test("$nm unc0.7",E,tgt,0.7)
end
println("# done")
