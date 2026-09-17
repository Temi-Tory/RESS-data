# BDD (sifted CUDD) comparison for the reliability-grounded drone case study -- answers R2.3/R3.2 directly:
# "IPA must also be benchmarked against state-of-the-art exact solvers... using the SAME grid AND DRONE
# networks" (Reviewer #2, comment 2). We only ever did this for the grid; this closes the gap for the
# drone case study. ONE network per process invocation (ARGS[1]) so a blowup on one doesn't cost the others.
#   (a) EXACTNESS: Float64 midpoint weights (mid of each Interval) -- IPA exact vs sifted-CUDD exact.
#   (b) TIMING + STRUCTURE: sifted-CUDD build time + node count vs IPA-Float64 propagation time.
#   (c) INTERVAL TIMING: IPA-interval (one pass, ACTUAL interval data) vs BDD-interval (build + 2 corner
#       evals at the actual lo/hi bounds), mirroring validation/interval_bdd_vs_ipa_timing.jl.
const NAME = ARGS[1]
const REPO = raw"c:\Users\ohian\OneDrive - University of Strathclyde\Documents\Programmming Files\Julia Files\InformationPropagation\Info_Prop_Framework_Project"
using Pkg; Pkg.activate(REPO); const BDDENV = joinpath(REPO,"validation","bddenv"); push!(LOAD_PATH, BDDENV)
include(joinpath(REPO,"InfoPropFrmwrk","src","Algorithms","InfoPropFramework.jl")); using .InfoPropFramework
using CUDD, Printf

const NODE_CAP = 5_000_000

function cudd_build_sifted(nid, eid, edgelist, all_nodes, incoming, sources)
    mgr = Cudd_Init(0, 0, 256, 262144, 0)
    Cudd_AutodynEnable(mgr, CUDD.CUDD_REORDER_SIFT)
    ith(id) = Cudd_bddIthVar(mgr, id-1)
    And(f,g) = (r = Cudd_bddAnd(mgr,f,g); Cudd_Ref(r); r)
    Or(f,g)  = (r = Cudd_bddOr(mgr,f,g);  Cudd_Ref(r); r)
    Zero = Cudd_ReadLogicZero(mgr)
    reach = Dict{Int,Ptr{Nothing}}()
    topo = all_nodes   # all_nodes is already topologically sorted by construction (DAG order)
    for n in topo
        base = ith(nid[n])
        if n in sources
            r = base
        else
            acc = Zero
            for u in get(incoming, n, Int[])
                acc = Or(acc, And(reach[u], ith(eid[(u,n)])))
            end
            r = And(base, acc)
        end
        Cudd_Ref(r); reach[n] = r
    end
    Cudd_ReduceHeap(mgr, CUDD.CUDD_REORDER_SIFT, 0)
    (; mgr, reach, total_nodes = Int(Cudd_ReadNodeCount(mgr)))
end
function cudd_eval(B, nid, all_nodes, w)
    one_reg = Ptr{Nothing}(Cudd_Regular(Cudd_ReadOne(B.mgr)))
    memo = Dict{Ptr{Nothing},Float64}()
    function prob(f)
        reg = Ptr{Nothing}(Cudd_Regular(f)); b = get(memo, reg, -1.0)
        if b < 0.0
            b = reg == one_reg ? 1.0 : (p = w[Int(Cudd_NodeReadIndex(reg))+1]; p*prob(Cudd_T(reg)) + (1-p)*prob(Cudd_E(reg)))
            memo[reg] = b
        end
        (Cudd_IsComplement(f) != 0) ? 1.0-b : b
    end
    Dict(n => prob(B.reach[n]) for n in all_nodes)
end

function main()
    base = joinpath(REPO,"dag_ntwrk_files",NAME)
    edgelist, _, _, source_nodes_vec = read_graph_to_dict(joinpath(base,"$NAME.EDGES"))
    node_priors = read_node_priors_from_json(joinpath(base,"interval","$NAME-nodepriors.json"))
    link_probabilities = read_edge_probabilities_from_json(joinpath(base,"interval","$NAME-linkprobabilities.json"))
    sources = Set(source_nodes_vec)
    all_nodes_unsorted = sort!(collect(union(Set(first.(edgelist)), Set(last.(edgelist)))))
    outgoing = Dict{Int64,Set{Int64}}(); incoming = Dict{Int64,Set{Int64}}()
    for (u,v) in edgelist
        push!(get!(outgoing,u,Set{Int64}()), v); push!(get!(incoming,v,Set{Int64}()), u)
    end
    itersets, anc, desc = find_iteration_sets(edgelist, outgoing, incoming)
    fk, jn = identify_fork_and_join_nodes(outgoing, incoming)
    topo = Int64[]; for s in itersets, n in sort(collect(s)); push!(topo, n); end

    println("$NAME: $(length(all_nodes_unsorted)) nodes, $(length(edgelist)) edges"); flush(stdout)

    # --- CUDD variable indexing ---
    nid = Dict(n => i for (i,n) in enumerate(topo))
    eid = Dict{Tuple{Int,Int},Int}(); for (k,e) in enumerate(edgelist); eid[e] = length(topo)+k; end
    nvar = length(topo) + length(edgelist)
    wlo = zeros(nvar); whi = zeros(nvar); wmid = zeros(nvar)
    node_mid = Dict{Int64,Float64}(); link_mid = Dict{Tuple{Int64,Int64},Float64}()
    for n in topo
        iv = node_priors[n]; wlo[nid[n]]=iv.lower; whi[nid[n]]=iv.upper; m=(iv.lower+iv.upper)/2
        wmid[nid[n]]=m; node_mid[n]=m
    end
    for (e,i) in eid
        iv = link_probabilities[e]; wlo[i]=iv.lower; whi[i]=iv.upper; m=(iv.lower+iv.upper)/2
        wmid[i]=m; link_mid[e]=m
    end

    # WARMUP-then-measure: a bare @elapsed on the FIRST call includes JIT/compilation, not runtime -- the
    # exact rigor gap Reviewer #2 already flagged once (R2.2: "runtime comparison is not fully convincing").
    ipa_f64_run() = begin
        r,u = new_identify(edgelist, node_mid, link_mid, sources, fk, jn, anc, desc, itersets)
        update_beliefs_iterative(edgelist, itersets, outgoing, incoming, sources, node_mid, link_mid, desc, anc, r, jn, fk, u)
    end
    ipa_iv_run() = begin
        r,u = new_identify(edgelist, node_priors, link_probabilities, sources, fk, jn, anc, desc, itersets)
        update_beliefs_iterative(edgelist, itersets, outgoing, incoming, sources, node_priors, link_probabilities, desc, anc, r, jn, fk, u)
    end
    bdd_build_once() = cudd_build_sifted(nid, eid, edgelist, topo, incoming, sources)

    # --- warmup pass (untimed, discarded) ---
    bel_ipa = ipa_f64_run()
    Bwarm = bdd_build_once(); Bwarm.total_nodes > NODE_CAP && (println("  SKIP: sifted bdd_nodes=$(Bwarm.total_nodes) exceeds cap $NODE_CAP"); flush(stdout); return)
    cudd_eval(Bwarm, nid, topo, wmid); Cudd_Quit(Bwarm.mgr)
    ipa_iv_run()

    # --- (a)+(b) exactness + Float64 timing (warm) ---
    t_ipa_f64 = @elapsed bel_ipa = ipa_f64_run()
    t_build = @elapsed Bc = bdd_build_once()
    t_eval_mid = @elapsed bel_cudd = cudd_eval(Bc, nid, topo, wmid)
    maxdiff = maximum(abs(bel_ipa[n] - bel_cudd[n]) for n in topo)
    @printf("  EXACT (warm): bdd_nodes=%d, t_build=%.3fs, t_eval=%.3fs, t_ipa_float64=%.3fs, max|IPA-BDD|=%.2e\n",
            Bc.total_nodes, t_build, t_eval_mid, t_ipa_f64, maxdiff); flush(stdout)
    Cudd_Quit(Bc.mgr)

    # --- (c) interval timing (warm): IPA-interval (one pass, real data) vs BDD-interval (build + 2 corner evals) ---
    t_ipa_iv = @elapsed ipa_iv_run()
    t_bdd_iv = @elapsed begin
        B2 = bdd_build_once()
        _ = cudd_eval(B2, nid, topo, wlo); _ = cudd_eval(B2, nid, topo, whi)
        Cudd_Quit(B2.mgr)
    end
    @printf("  INTERVAL TIMING (warm): t_ipa_interval=%.3fs, t_bdd_interval(build+2eval)=%.3fs, ratio(bdd/ipa)=%.2fx\n",
            t_ipa_iv, t_bdd_iv, t_bdd_iv/t_ipa_iv); flush(stdout)
end
main()
println("# done: $NAME"); flush(stdout)
