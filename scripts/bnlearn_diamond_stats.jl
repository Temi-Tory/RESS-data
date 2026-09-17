# Safe identify-only diagnostic (polynomial, no exponential enumeration) across the newly-converted
# bnlearn discrete-network corpus, before attempting any full BDD comparison. LINK in particular is a
# known high-treewidth benchmark in the PGM exact-inference literature despite its modest node count --
# do not assume any of these "just work" without checking maxcond first (lesson from the drone case study).
const REPO = raw"c:\Users\ohian\OneDrive - University of Strathclyde\Documents\Programmming Files\Julia Files\InformationPropagation\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(REPO)
include(joinpath(REPO,"InfoPropFrmwrk","src","Algorithms","InfoPropFramework.jl")); using .InfoPropFramework
using Printf

function diamond_stats(name)
    base = joinpath(REPO,"dag_ntwrk_files",name)
    edgelist, _, _, source_nodes_vec = read_graph_to_dict(joinpath(base,"$name.EDGES"))
    node_priors = read_node_priors_from_json(joinpath(base,"float","$name-nodepriors.json"))
    link_probabilities = read_edge_probabilities_from_json(joinpath(base,"float","$name-linkprobabilities.json"))
    source_nodes = Set(source_nodes_vec)
    all_nodes = sort!(collect(union(Set(first.(edgelist)), Set(last.(edgelist)))))
    outgoing = Dict{Int64,Set{Int64}}(); incoming = Dict{Int64,Set{Int64}}()
    for (u,v) in edgelist
        push!(get!(outgoing,u,Set{Int64}()), v); push!(get!(incoming,v,Set{Int64}()), u)
    end
    t_iter = @elapsed (itersets, anc, desc) = find_iteration_sets(edgelist, outgoing, incoming)
    fk, jn = identify_fork_and_join_nodes(outgoing, incoming)
    t_id = @elapsed (r, u) = new_identify(edgelist, node_priors, link_probabilities, source_nodes, fk, jn, anc, desc, itersets)
    nuniq = length(u)
    maxcond = isempty(u) ? 0 : maximum(length(cd.diamond.conditioning_nodes) for (_,cd) in u)
    sum2c = isempty(u) ? 0.0 : sum(2.0^length(cd.diamond.conditioning_nodes) for (_,cd) in u)
    @printf("%-14s V=%-4d E=%-5d src=%-3d forks=%-4d joins=%-4d  t_iter=%.2fs t_id=%.2fs  uniq=%-5d maxcond=%-3d sum2^C=%.2e\n",
            name, length(all_nodes), length(edgelist), length(source_nodes), length(fk), length(jn), t_iter, t_id, nuniq, maxcond, sum2c)
    flush(stdout)
end

for name in ("cancer-bnlearn","earthquake-bnlearn","sachs-bnlearn","survey-bnlearn","alarm-bnlearn",
             "child-bnlearn","insurance-bnlearn","barley-bnlearn","mildew-bnlearn","hailfinder-bnlearn",
             "hepar2-bnlearn","win95pts-bnlearn","andes-bnlearn","pathfinder-bnlearn","pigs-bnlearn",
             "diabetes-bnlearn","link-bnlearn")
    try
        diamond_stats(name)
    catch e
        println("$name: FAILED - $(sprint(showerror, e))"); flush(stdout)
    end
end
println("# bnlearn diamond stats done"); flush(stdout)
