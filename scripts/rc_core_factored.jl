# rc_core_factored — the recursive-conditioning reference (see rc_core.jl) EXTENDED with
# INDEPENDENT-DIAMOND FACTORIZATION. At a join, the parents' "reach" events are combined by
# inclusion-exclusion; when the parents partition into groups with DISJOINT un-conditioned ancestry
# (given the current conditioning), those groups are INDEPENDENT, so the OR factorizes and each group is
# conditioned separately. This takes fanin-k from 2^k to O(k) while staying exact. Independence is
# context-sensitive: fixing a shared ancestor can split a group (its influence is then blocked).
#
# Standalone; no framework. Returns (beliefs, op_count) where op_count = number of reach_via/R
# evaluations (a machine-independent proxy for work), to compare against the non-factored rc_core.
module RCCoreFactored
export rc_reliability_factored

function _build(edgelist)
    out = Dict{Int,Vector{Int}}(); inc = Dict{Int,Vector{Int}}(); nodes = Set{Int}()
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

function rc_reliability_factored(edgelist::Vector{Tuple{Int,Int}}, node_priors::Dict{Int,Float64},
                                 link_probs::Dict{Tuple{Int,Int},Float64}, sources)
    out, inc, nodes = _build(edgelist)
    topo = _toposort(nodes, out, inc); topi = Dict(v => i for (i,v) in enumerate(topo))
    srcset = Set(sources)
    forks = Set(v for v in nodes if length(get(out,v,Int[])) >= 2)
    # exclude from conditioning only DETERMINISTIC-reachability nodes: dead (prior 0) or certain sources.
    det   = Set(v for v in nodes if get(node_priors,v,0.5) == 0.0 || (get(node_priors,v,0.5) == 1.0 && v in srcset))
    anc = Dict{Int,Set{Int}}()
    for v in topo
        s = Set{Int}(); for u in get(inc,v,Int[]); push!(s,u); union!(s, anc[u]); end; anc[v] = s
    end
    ops = Ref(0)

    # un-conditioned "influence set" of a parent p: itself + its ancestors, minus conditioned nodes
    infl(p, cond) = begin
        s = Set{Int}(); (p in keys(cond)) || push!(s, p)
        for a in anc[p]; (a in keys(cond)) || push!(s, a); end
        s
    end
    # partition pset into independent groups by shared un-conditioned ancestry (union-find)
    function components(pset, cond)
        pv = collect(pset); n = length(pv); parent = collect(1:n)
        find(i) = (while parent[i] != i; parent[i]=parent[parent[i]]; i=parent[i]; end; i)
        infls = [infl(p, cond) for p in pv]
        for i in 1:n, j in i+1:n
            isempty(intersect(infls[i], infls[j])) || (parent[find(i)] = find(j))
        end
        groups = Dict{Int,Vector{Int}}()
        for i in 1:n; push!(get!(groups, find(i), Int[]), pv[i]); end
        collect(values(groups))
    end
    # highest un-conditioned shared fork/asymmetric-parent among pset (or nothing)
    function shared_fork(pset, cond)
        cover = Dict{Int,Int}()
        for p in pset
            for a in anc[p]; (a in forks) && !(a in keys(cond)) && !(a in det) && (cover[a]=get(cover,a,0)+1); end
        end
        for p in pset, q in pset
            p==q && continue
            (p in anc[q]) && !(p in keys(cond)) && !(p in det) && (cover[p]=get(cover,p,0)+1)
        end
        best=nothing; bestk=typemax(Int)
        for (a,c) in cover; c<2 && continue; topi[a]<bestk && (bestk=topi[a]; best=a); end
        best
    end

    memoR = Dict{Tuple{Int,Vector{Tuple{Int,Bool}}}, Float64}()
    function R(v, cond)
        haskey(cond, v) && return cond[v] ? 1.0 : 0.0
        v in srcset && return node_priors[v]
        key = (v, sort([(c,cond[c]) for c in keys(cond) if c in anc[v]]))
        haskey(memoR, key) && return memoR[key]
        ops[] += 1
        ps = get(inc, v, Int[])
        val = isempty(ps) ? node_priors[v] : node_priors[v] * reach_via(v, Set(ps), cond)
        memoR[key] = val; val
    end
    # P(v reached via >=1 parent in pset | cond), excluding v's own prior
    function reach_via(v, pset, cond)
        ops[] += 1
        length(pset) == 1 && (p = first(pset); return R(p, cond) * link_probs[(p, v)])
        groups = components(pset, cond)
        if length(groups) > 1                       # independent groups -> factorize the OR
            surv = 1.0
            for G in groups; surv *= (1.0 - reach_via(v, Set(G), cond)); end
            return 1.0 - surv
        end
        f = shared_fork(pset, cond)
        if f === nothing                            # independent parents -> inclusion-exclusion
            surv = 1.0; for p in pset; surv *= (1.0 - R(p,cond)*link_probs[(p,v)]); end; return 1.0 - surv
        end
        rf = R(f, cond)                             # condition on the shared fork
        cu = copy(cond); cu[f] = true; cd = copy(cond); cd[f] = false
        rf * reach_via(v, pset, cu) + (1.0 - rf) * reach_via(v, pset, cd)
    end

    beliefs = Dict(v => (v in srcset ? node_priors[v] : R(v, Dict{Int,Bool}())) for v in nodes)
    (beliefs, ops[])
end

end # module
