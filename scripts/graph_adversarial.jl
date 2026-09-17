# Adversarial DAG families that expose the complementary worst cases of IPA (recursive conditioning)
# vs ROBDD, for the complexity study. Each returns (edges::Vector{Tuple{Int,Int}}, sources::Vector{Int}).
#
#   gen_fanin_k(k)   — k INDEPENDENT fork gadgets reconverging at one sink. The sink's parents are
#                      correlated through k distinct forks, so IPA must condition on all k -> 2^k. The
#                      reachability function is an OR of k independent terms, so the ROBDD is LINEAR in k.
#                      => IPA's structural worst case, BDD's best case. (Multi-source: the k forks.)
#
#   gen_mesh(w,L)    — a directed w×L lattice: (r,c) -> (r,c+1) and (r,c) -> (r+1,c+1). Reconvergence at
#                      every interior node; the "active frontier" under any variable order grows with the
#                      width w, stressing the ROBDD. IPA sees a regular nesting. => probes BDD's weak case.
#                      (Multi-source: column 1.)

# k independent fork->join gadgets into a single sink.
function gen_fanin_k(k::Int)
    edges = Tuple{Int,Int}[]
    f(i) = i                       # forks 1..k (sources)
    a(i) = k + 2i - 1              # first child
    b(i) = k + 2i                  # second child
    t = 3k + 1                     # sink
    for i in 1:k
        push!(edges, (f(i), a(i))); push!(edges, (f(i), b(i)))
        push!(edges, (a(i), t));    push!(edges, (b(i), t))
    end
    (edges, collect(1:k))
end

# directed w×L lattice with straight + down-diagonal edges.
function gen_mesh(w::Int, L::Int)
    id(r,c) = (c-1)*w + r          # column-major numbering keeps a topological order (c increasing)
    edges = Tuple{Int,Int}[]
    for c in 1:L-1, r in 1:w
        push!(edges, (id(r,c), id(r,c+1)))                 # straight
        r < w && push!(edges, (id(r,c), id(r+1,c+1)))      # down-diagonal
    end
    sources = [id(r,1) for r in 1:w]
    (edges, sources)
end
