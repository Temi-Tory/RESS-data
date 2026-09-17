# K-sweep: for the two K-bounded drone proxy networks, measure maxcond AND actual full-propagation wall
# time (not just new_identify) as K increases, entirely in-memory (no disk I/O per K). Answers: "if K=12
# runs in seconds, how far can we push K before it stops being seconds?"
const REPO = raw"c:\Users\ohian\OneDrive - University of Strathclyde\Documents\Programmming Files\Julia Files\InformationPropagation\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(REPO)
include(joinpath(REPO,"InfoPropFrmwrk","src","Algorithms","InfoPropFramework.jl")); using .InfoPropFramework
using Printf

# Reuse the network-construction logic (does NOT auto-run generation: guarded by abspath(PROGRAM_FILE)).
include(joinpath(REPO,"drone_network_to_dag_reliability.jl"))

function to_framework_types(edges, edge_probs_tuples, nodes_df, active_ids)
    id_to_row = Dict(nodes_df.numberID[i] => i for i in 1:nrow(nodes_df))
    node_priors = Dict{Int64,Interval}(nid => Interval(node_prior_interval(nodes_df[id_to_row[nid], :])...) for nid in active_ids)
    link_probabilities = Dict{Tuple{Int64,Int64},Interval}((s,d) => Interval(v...) for ((s,d), v) in
        Dict((parse.(Int, split(k[2:end-1], ",")) |> t -> (t[1], t[2])) => v for (k,v) in edge_probs_tuples))
    node_priors, link_probabilities
end

function run_at_K(builder, nodes_df, vtol_matrix, fw_matrix, order_lookup, K)
    edges, edge_probs = builder(nodes_df, vtol_matrix, fw_matrix, order_lookup; K=K)
    active_ids = active_ids_from_edges(edges)
    node_priors, link_probabilities = to_framework_types(edges, edge_probs, nodes_df, active_ids)
    source_nodes = Set(n for n in active_ids if !any(e -> e[2]==n, edges))
    outgoing = Dict{Int64,Set{Int64}}(); incoming = Dict{Int64,Set{Int64}}()
    for (u,v) in edges
        push!(get!(outgoing,u,Set{Int64}()), v); push!(get!(incoming,v,Set{Int64}()), u)
    end
    itersets, anc, desc = find_iteration_sets(edges, outgoing, incoming)
    fk, jn = identify_fork_and_join_nodes(outgoing, incoming)
    t_id = @elapsed (r, u) = new_identify(edges, node_priors, link_probabilities, source_nodes, fk, jn, anc, desc, itersets)
    maxcond = isempty(u) ? 0 : maximum(length(cd.diamond.conditioning_nodes) for (_,cd) in u)
    t_prop = @elapsed bel = update_beliefs_iterative(edges, itersets, outgoing, incoming, source_nodes,
                             node_priors, link_probabilities, desc, anc, r, jn, fk, u)
    n_sound = count(b -> 0.0 <= b.lower <= b.upper <= 1.0 + 1e-9, values(bel))
    (length(edges), maxcond, t_id, t_prop, n_sound, length(bel))
end

nodes_df, vtol_matrix, fw_matrix = load_drone_network_data()
order_lookup = hub_spoke_order(nodes_df)

for (label, builder) in (("vtol-dense-decentralized", build_vtol_dense_decentralized),
                          ("concentrated-minimal", build_concentrated_minimal))
    println("=== $label ==="); flush(stdout)
    # warmup at K=3 to exclude JIT from the timings that matter
    run_at_K(builder, nodes_df, vtol_matrix, fw_matrix, order_lookup, 3)
    for K in (12, 16, 20, 24)
        try
            ne, maxcond, t_id, t_prop, n_sound, n_bel = run_at_K(builder, nodes_df, vtol_matrix, fw_matrix, order_lookup, K)
            @printf("  K=%-3d edges=%-5d maxcond=%-3d identify=%.2fs propagate=%.2fs sound=%d/%d\n",
                    K, ne, maxcond, t_id, t_prop, n_sound, n_bel); flush(stdout)
        catch e
            println("  K=$K FAILED: $(sprint(showerror, e))"); flush(stdout)
            break
        end
    end
    println(); flush(stdout)
end
println("# K sweep done"); flush(stdout)
