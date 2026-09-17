# Full band-vs-threshold curve for BOTH p-box operators (cvxP default, cvxF proven), on the SAME grid
# network and scenario, to see WHERE they diverge across the whole [0,1] range rather than trusting a
# single x*=0.95 point (which turned out to be a favourable, non-representative spot -- see
# certified_bound_vignette.jl vs certified_bound_vignette_frechet.jl, 2026-07-28). One p-box propagation
# PER operator per scenario; the threshold sweep itself is free (reading the same p-box's CDF at many x).
const REPO = raw"c:\Users\ohian\OneDrive - University of Strathclyde\Documents\Programmming Files\Julia Files\InformationPropagation\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(REPO); push!(LOAD_PATH, joinpath(REPO,"validation","bddenv"))
include(joinpath(REPO,"InfoPropFrmwrk","src","Algorithms","InfoPropFramework.jl")); using .InfoPropFramework
using Random, Printf, Distributions
include(joinpath(REPO,"validation","graph_gen.jl")); include(joinpath(REPO,"validation","oracles.jl"))
quiet(f)=redirect_stdout(f, devnull)
glo(c)=hasproperty(c,:lo) ? glo(c.lo) : Float64(c); ghi(c)=hasproperty(c,:hi) ? ghi(c.hi) : Float64(c)
PBA.setSteps(50)
const STEPS = 50
const Z = 1.96
const XSTARS = (0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,
                 0.75, 0.8, 0.85, 0.9, 0.95, 0.99)

mkdist(kind,v,w) = begin
    a=max(0.0,v-w); b=min(1.0,v+w)
    kind==:uni  ? Uniform(a,b) : TriangularDist(a,b, clamp(v,a,b))
end
pbfromdist(d) = (qs=quantile.(d,[(i-0.5)/STEPS for i in 1:STEPS]); PBA.pbox(qs,qs))
ONEPB() = PBA.makepbox(PBA.interval(1.0,1.0))

g = load_edges("grid", joinpath(REPO,"dag_ntwrk_files","grid-graph","grid-graph.EDGES"))
P = make_problem(g); fk,jn = identify_fork_and_join_nodes(P.outgoing,P.incoming); TGT=16

function run(kind, w, perfect)
    dl = mkdist(kind, 0.9, w)
    np = perfect ? Dict{Int64,pbox}(n=>ONEPB() for n in P.all_nodes) :
                   Dict{Int64,pbox}(n=>pbfromdist(mkdist(kind,0.7,w)) for n in P.all_nodes)
    lp = Dict{Tuple{Int64,Int64},pbox}(e=>pbfromdist(dl) for e in keys(P.eid))
    quiet() do
        r,u=new_identify(P.edgelist,np,lp,Set{Int64}(P.sources),fk,jn,P.anc,P.desc,P.itersets)
        update_beliefs_iterative(P.edgelist,P.itersets,P.outgoing,P.incoming,P.sources,np,lp,P.desc,P.anc,r,jn,fk,u,Dict{CacheKey,DiamondCacheEntry{pbox}}())
    end[TGT]
end

function mc_samples(kind, w, perfect; N=6000)
    dl = mkdist(kind, 0.9, w); rng=MersenneTwister(42); samp=Float64[]
    dn = perfect ? nothing : mkdist(kind,0.7,w)
    for _ in 1:N
        npf=Dict{Int64,Float64}(n=> dn===nothing ? 1.0 : rand(rng,dn) for n in P.all_nodes)
        lpf=Dict{Tuple{Int64,Int64},Float64}(e=>rand(rng,dl) for e in keys(P.eid))
        r,u=new_identify(P.edgelist,npf,lpf,Set{Int64}(P.sources),fk,jn,P.anc,P.desc,P.itersets)
        b=update_beliefs_iterative(P.edgelist,P.itersets,P.outgoing,P.incoming,P.sources,npf,lpf,P.desc,P.anc,r,jn,fk,u); push!(samp,b[TGT])
    end
    samp
end

scenarios = [("perfect", true, 0.10), ("uncert0.7", false, 0.10)]

open(joinpath(REPO,"InfoPropFrmwrk","Publications","My work","RESS_response","data","certified_bound_threshold_sweep.csv"),"w") do io
    println(io,"regime,operator,xstar,ipa_lo,ipa_hi,band,mc_phat,N_worst")
    for (regname, perfect, w) in scenarios
        samp = mc_samples(:tri, w, perfect)
        for (opname, blend) in (("cvxP", :positive), ("cvxF", :frechet))
            InfoPropFramework.InputProcessingModule.PBOX_COND_BLEND[] = blend
            bel = run(:tri, w, perfect)
            for xstar in XSTARS
                c = PBA.cdf(bel, xstar); a,b = glo(c), ghi(c); band = max(b-a,0.0)
                phat = count(<=(xstar), samp)/length(samp)
                Nworst = band/2 < 1e-9 ? -1 : ceil(Int, Z^2*0.25/(band/2)^2)
                @printf(io, "%s,%s,%.2f,%.4f,%.4f,%.4f,%.4f,%d\n", regname, opname, xstar, a, b, band, phat, Nworst)
                @printf("%-9s %-4s x*=%.2f  [%.3f,%.3f] band=%.3f  phat=%.3f  N_worst=%d\n",
                        regname, opname, xstar, a, b, band, phat, Nworst); flush(stdout)
            end
        end
    end
end
println("# done -> data/certified_bound_threshold_sweep.csv")
