# Convert a bnlearn BIF file (plain text, already gunzipped) into our EDGES + Float64 float/ prior format.
# BIF grammar used: `variable NAME { ... }` declares nodes; `probability ( CHILD | P1, P2, ... )` or
# `probability ( CHILD )` declares the DAG structure (edges P_i -> CHILD). We reuse ONLY the DAG topology
# (a real, standard, widely-cited PGM benchmark structure) -- NOT the actual conditional probability
# tables, which are categorical CPTs and do not map onto our binary up/down reachability model. Synthetic
# Float64 reliabilities are assigned the same way as for power-network/KarlNetwork/the random corpus
# (arbitrary but valid draws) -- legitimate for a pure exactness/structural validation claim, not a
# decision-relevant one (see PBOX_HANDOFF-adjacent discussion, 2026-07-28).
using Random, JSON

function parse_bif(path)
    text = read(path, String)
    names = String[]
    for m in eachmatch(r"variable\s+(\S+)\s*\{"m, text)
        push!(names, m.captures[1])
    end
    id = Dict(n => i for (i, n) in enumerate(names))
    edges = Tuple{Int,Int}[]
    for m in eachmatch(r"probability\s*\(\s*([^|)]+?)\s*(?:\|\s*([^)]+))?\)\s*\{"m, text)
        child = strip(m.captures[1])
        haskey(id, child) || continue
        if m.captures[2] !== nothing
            for p in split(m.captures[2], ',')
                p = strip(p)
                haskey(id, p) || continue
                push!(edges, (id[p], id[child]))
            end
        end
    end
    unique!(edges)
    (names, id, edges)
end

function write_network(name, names, edges, out_dag_dir)
    outdir = joinpath(out_dag_dir, name)
    mkpath(outdir); mkpath(joinpath(outdir, "float"))
    open(joinpath(outdir, "$name.EDGES"), "w") do f
        println(f, "source,destination")
        for (s, d) in edges; println(f, "$s,$d"); end
    end
    rng = MersenneTwister(7)
    node_ids = 1:length(names)
    dest_set = Set(d for (_, d) in edges)
    priors = Dict(string(n) => (n in dest_set ? 0.3 + 0.69 * rand(rng) : 1.0) for n in node_ids)
    links = Dict("($s,$d)" => 0.3 + 0.69 * rand(rng) for (s, d) in edges)
    open(joinpath(outdir, "float", "$name-nodepriors.json"), "w") do f
        JSON.print(f, Dict("nodes" => priors, "data_type" => "Float64",
            "description" => "Synthetic reliabilities on the real bnlearn '$name' DAG topology -- structural exactness validation only, not a decision-relevant claim about the represented system."), 2)
    end
    open(joinpath(outdir, "float", "$name-linkprobabilities.json"), "w") do f
        JSON.print(f, Dict("links" => links, "data_type" => "Float64",
            "description" => "Synthetic reliabilities on the real bnlearn '$name' DAG topology."), 2)
    end
    println("  $name: $(length(names)) nodes, $(length(edges)) edges -> $outdir")
end

if abspath(PROGRAM_FILE) == @__FILE__
    bif_path, name, out_dir = ARGS[1], ARGS[2], ARGS[3]
    names, id, edges = parse_bif(bif_path)
    write_network(name, names, edges, out_dir)
end
