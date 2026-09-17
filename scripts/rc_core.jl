# ---------------------------------------------------------------------------------------------
# rc_core — reference correct-by-construction network-reachability reliability via RECURSIVE
# CONDITIONING (the "diamond" idea done right). Standalone, no framework dependency.
#
# belief(v) = P(v operational AND reachable from some source). For a join whose parents are
# correlated, we condition on a shared fork ancestor via EXACT total probability
#     P(reach v) = P(f reached)·P(reach v | f) + P(¬f)·P(reach v | ¬f)
# recursing (one fork at a time) until the parents share no un-conditioned fork ancestor — then
# they are conditionally independent and combine by inclusion-exclusion. Conditioning on one node
# at a time keeps the joint exact (no independence assumption on the conditioning set), which is
# what the old completeness-loop got wrong (Bug #2). Memoized on (v, cond ∩ anc(v)).
# ---------------------------------------------------------------------------------------------
module RCCore
export rc_reliability

function _build(edgelist)
    out = Dict{Int,Vector{Int}}(); inc = Dict{Int,Vector{Int}}()
    nodes = Set{Int}()
    for (u,v) in edgelist
        push!(get!(out,u,Int[]), v); push!(get!(inc,v,Int[]), u); push!(nodes,u); push!(nodes,v)
    end
    (out, inc, sort(collect(nodes)))
end
function _toposort(nodes, out, inc)
    indeg = Dict(v => length(get(inc,v,Int[])) for v in nodes)
    q = sort([v for v in nodes if indeg[v]==0]); order = Int[]
    while !isempty(q)
        u = popfirst!(q); push!(order,u)
        for w in get(out,u,Int[]); indeg[w]-=1; indeg[w]==0 && push!(q,w); end
    end
    order
end

"""
    rc_reliability(edgelist, node_priors, link_probs, sources) -> Dict{Int,Float64}

Exact per-node reachability reliability via recursive conditioning.
"""
function rc_reliability(edgelist::Vector{Tuple{Int,Int}}, node_priors::Dict{Int,Float64},
                        link_probs::Dict{Tuple{Int,Int},Float64}, sources)
    out, inc, nodes = _build(edgelist)
    topo = _toposort(nodes, out, inc)
    topi = Dict(v => i for (i,v) in enumerate(topo))
    srcset = Set(sources)
    forks = Set(v for v in nodes if length(get(out,v,Int[])) >= 2)
    # exclude from conditioning only nodes with DETERMINISTIC reachability: dead (prior 0) or certain
    # SOURCES (prior 1). A prior-1 non-source fork has uncertain reachability -> must stay conditionable.
    det   = Set(v for v in nodes if get(node_priors,v,0.5) == 0.0 || (get(node_priors,v,0.5) == 1.0 && v in srcset))

    # ancestors (all), in topo order
    anc = Dict{Int,Set{Int}}()
    for v in topo
        s = Set{Int}()
        for u in get(inc,v,Int[]); push!(s,u); union!(s, anc[u]); end
        anc[v] = s
    end

    # pick the topologically-highest un-conditioned shared fork/ancestor of J's parents (or nothing)
    function pick_cond(v, cond)
        ps = get(inc, v, Int[]); length(ps) < 2 && return nothing
        cover = Dict{Int,Int}()
        for p in ps
            # candidate ancestors of this parent that could induce correlation: its fork ancestors,
            # plus p itself if p is an ancestor of another parent (asymmetric case)
            for a in anc[p]
                (a in forks) && !(a in keys(cond)) && !(a in det) && (cover[a] = get(cover,a,0)+1)
            end
            if p in forks && !(p in keys(cond)) && !(p in det)
                # p as a shared ancestor of other parents
                cover[p] = get(cover,p,0)  # ensure key exists; counted below
            end
        end
        # also count parent-as-ancestor-of-parent
        for p in ps, q in ps
            p == q && continue
            if p in anc[q] && !(p in keys(cond)) && !(p in det)
                cover[p] = get(cover,p,0) + 1
            end
        end
        best = nothing; bestkey = (typemax(Int),)
        for (a,c) in cover
            c < 2 && continue
            k = (topi[a],)
            if k < bestkey; bestkey = k; best = a; end
        end
        best
    end

    memo = Dict{Tuple{Int,Vector{Tuple{Int,Bool}}}, Float64}()
    function R(v, cond)
        haskey(cond, v) && return cond[v] ? 1.0 : 0.0
        v in srcset && return node_priors[v]
        av = anc[v]
        key = (v, sort([(c, cond[c]) for c in keys(cond) if c in av]))
        haskey(memo, key) && return memo[key]
        ps = get(inc, v, Int[])
        f = pick_cond(v, cond)
        local val
        if f === nothing
            # parents independent given cond -> inclusion-exclusion (noisy-OR)
            surv = 1.0
            for u in ps; surv *= (1.0 - R(u, cond) * link_probs[(u,v)]); end
            val = node_priors[v] * (1.0 - surv)
        else
            rf = R(f, cond)
            cu = copy(cond); cu[f] = true
            cd = copy(cond); cd[f] = false
            val = rf * R(v, cu) + (1.0 - rf) * R(v, cd)
        end
        memo[key] = val
        val
    end

    Dict(v => (v in srcset ? node_priors[v] : R(v, Dict{Int,Bool}())) for v in nodes)
end

end # module
