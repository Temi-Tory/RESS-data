# Table-1 data for the named structured networks: V, E, #sources, #root-diamonds, #unique-diamonds,
# max conditioning size, sifted ROBDD nodes, and worst |Δ| vs exact CUDD (IPA = new_identify -> prop).
const REPO = raw"c:\Users\ohian\OneDrive - University of Strathclyde\Documents\Programmming Files\Julia Files\InformationPropagation\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(REPO)
push!(LOAD_PATH, joinpath(REPO,"validation","bddenv"))
include(joinpath(REPO,"InfoPropFrmwrk","src","Algorithms","InfoPropFramework.jl")); using .InfoPropFramework
using CUDD, Random, Printf
include(joinpath(REPO,"validation","graph_gen.jl")); include(joinpath(REPO,"validation","oracles.jl"))
include(joinpath(REPO,"validation","new_identify.jl"))
include(joinpath(REPO,"validation","bdd_oracle.jl")); using .BDDReliabilityOracle

function ipa_run(P, D)
    fork, join = identify_fork_and_join_nodes(P.outgoing, P.incoming)
    roots, uniq = new_identify(P.edgelist, D.node_priors, D.link_probs, Set{Int64}(P.sources), fork, join, P.anc, P.desc, P.itersets)
    bel = update_beliefs_iterative(P.edgelist, P.itersets, P.outgoing, P.incoming, P.sources, D.node_priors, D.link_probs, P.desc, P.anc, roots, join, fork, uniq)
    maxc = isempty(uniq) ? 0 : maximum(length(cd.diamond.conditioning_nodes) for (_,cd) in uniq)
    (length(roots), length(uniq), maxc, bel)
end

names = ["power-network","grid-graph","KarlNetwork","counterexample-n15"]
println("% name & V & E & sources & rootD & uniqD & maxcond & sift_nodes & maxAbsDelta \\\\")
for nm in names
    path = joinpath(REPO,"dag_ntwrk_files",nm,"$nm.EDGES")
    isfile(path) || (println("% (skip $nm)"); continue)
    g = load_edges(nm, path); P = make_problem(g); D = draw_probs(P, MersenneTwister(1007))
    nr, nu, mc, bel = ipa_run(P, D)
    Bc = cudd_build(P); cud = cudd_eval(Bc, P, D.w); cudd_free(Bc)
    worst = maximum(abs(bel[n]-cud[n]) for n in P.all_nodes)
    np = Dict{Int,Float64}(n=>D.node_priors[n] for n in P.all_nodes)
    lp = Dict{Tuple{Int,Int},Float64}(e=>D.link_probs[e] for e in keys(P.eid))
    _, st = bdd_reliability(collect(P.edgelist), np, lp, collect(P.sources); sift=true)
    nsrc = length(P.sources)
    @printf("%-18s & %d & %d & %d & %d & %d & %d & %d & %.1e \\\\\n",
            nm, length(P.all_nodes), length(g.edges), nsrc, nr, nu, mc, st.bdd_total_nodes, worst)
end
