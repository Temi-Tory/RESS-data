# Tiered exact/approximate oracles for reachability reliability, to validate IPA beyond where the ROBDD
# is tractable. Reliability model: a path works iff every NODE and every EDGE on it is up; belief(v) =
# P(v operational AND reachable from some source). Standalone (no framework).
#
#   path_enum_reliability   — EXACT, per node, via source->node path enumeration + inclusion-exclusion
#                             over path events (each intersection = product over the UNION of components).
#                             Cost ~ 2^(#paths to node); guarded by a path/term cap -> returns nothing.
#   monte_carlo_reliability — APPROXIMATE, per node, with per-node standard error. Always tractable.
#
# For IE we exploit: P(∪ path_i up) = Σ_{∅≠S⊆paths} (-1)^{|S|+1} ∏_{c ∈ ⋃_{i∈S} comp(path_i)} p(c),
# because a component is up-or-not once (shared components appear once in the union) — exact.

module TieredOracles
using Random
export path_enum_reliability, monte_carlo_reliability

# all simple source->target paths (as sets of node-ids and edge-tuples); capped.
function _paths_to(target, incoming, sources, cap)
    # returns Vector of (nodes::Set{Int}, edges::Vector{Tuple{Int,Int}}) or nothing if > cap
    result = Vector{Tuple{Set{Int},Vector{Tuple{Int,Int}}}}()
    stack = Tuple{Int,Set{Int},Vector{Tuple{Int,Int}}}[(target, Set([target]), Tuple{Int,Int}[])]
    while !isempty(stack)
        v, vis, eds = pop!(stack)
        if v in sources
            push!(result, (vis, eds)); length(result) > cap && return nothing; continue
        end
        ps = get(incoming, v, Int[])
        isempty(ps) && continue                 # dead end (non-source with no parents)
        for u in ps
            u in vis && continue
            push!(stack, (u, union(vis, Set([u])), vcat(eds, [(u, v)])))
        end
    end
    result
end

"""
    path_enum_reliability(edgelist, node_priors, link_probs, sources; targets, path_cap=18)
Exact per-node reliability by path enumeration + inclusion-exclusion. Returns Dict(node=>belief) for the
requested `targets` (default all nodes), or `nothing` for any target whose path count exceeds `path_cap`
(caller can fall back to Monte Carlo).
"""
function path_enum_reliability(edgelist::Vector{Tuple{Int,Int}}, node_priors::Dict{Int,Float64},
                               link_probs::Dict{Tuple{Int,Int},Float64}, sources;
                               targets=nothing, path_cap::Int=18)
    incoming = Dict{Int,Vector{Int}}(); nodes = Set{Int}()
    for (u,v) in edgelist; push!(get!(incoming,v,Int[]),u); push!(nodes,u); push!(nodes,v); end
    srcset = Set(sources)
    tgts = targets === nothing ? sort(collect(nodes)) : targets
    out = Dict{Int,Float64}()
    for t in tgts
        if t in srcset; out[t] = node_priors[t]; continue; end
        paths = _paths_to(t, incoming, srcset, path_cap)
        paths === nothing && return nothing           # too many paths -> signal fallback
        np = length(paths)
        np == 0 && (out[t] = 0.0; continue)
        # inclusion-exclusion over 2^np subsets
        np > 20 && return nothing
        rel = 0.0
        for mask in 1:(2^np - 1)
            comp_nodes = Set{Int}(); comp_edges = Set{Tuple{Int,Int}}(); bits = 0
            for i in 1:np
                if (mask >> (i-1)) & 1 == 1
                    bits += 1; union!(comp_nodes, paths[i][1]); union!(comp_edges, Set(paths[i][2]))
                end
            end
            pr = 1.0
            for n in comp_nodes; pr *= node_priors[n]; end
            for e in comp_edges; pr *= link_probs[e]; end
            rel += (iseven(bits) ? -pr : pr)
        end
        out[t] = rel
    end
    out
end

"""
    monte_carlo_reliability(edgelist, node_priors, link_probs, sources; N=2_000_000, seed=1)
Approximate per-node reliability by sampling N component-worlds. Returns (belief::Dict, stderr::Dict).
"""
function monte_carlo_reliability(edgelist::Vector{Tuple{Int,Int}}, node_priors::Dict{Int,Float64},
                                 link_probs::Dict{Tuple{Int,Int},Float64}, sources; N::Int=2_000_000, seed::Int=1)
    using_nodes = sort(collect(union(Set(u for (u,_) in edgelist), Set(v for (_,v) in edgelist))))
    srcset = Set(sources)
    outadj = Dict{Int,Vector{Tuple{Int,Tuple{Int,Int}}}}()   # u -> list of (child v, edge (u,v))
    for (u,v) in edgelist; push!(get!(outadj,u,Tuple{Int,Tuple{Int,Int}}[]), (v,(u,v))); end
    cnt = Dict(n => 0 for n in using_nodes)
    rng = MersenneTwister(seed)
    for _ in 1:N
        upn = Dict(n => (rand(rng) < node_priors[n]) for n in using_nodes)
        reached = Set{Int}()
        q = Int[]
        for s in srcset; if upn[s]; push!(q, s); push!(reached, s); end; end
        while !isempty(q)
            u = popfirst!(q)
            for (v,e) in get(outadj, u, Tuple{Int,Int}[])
                if v ∉ reached && upn[v] && rand(rng) < link_probs[e]
                    push!(reached, v); push!(q, v)
                end
            end
        end
        for v in reached; cnt[v] += 1; end
    end
    bel = Dict(n => cnt[n]/N for n in using_nodes)
    err = Dict(n => sqrt(max(bel[n]*(1-bel[n]),0.0)/N) for n in using_nodes)
    (bel, err)
end

end # module
