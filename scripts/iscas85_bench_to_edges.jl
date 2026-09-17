# ISCAS85 .bench -> project .EDGES + float priors/link-probs converter.
#
# .bench format (plain text, ttu.ee/filebox.ece.vt.edu mirrors, standard since 1985):
#   comment lines start with '#'
#   INPUT(id)              -- primary input signal, no incoming edges
#   OUTPUT(id)              -- primary output signal (informational only here: the id already
#                              exists as a gate output or an INPUT, so no separate node/edge needed)
#   id = GATETYPE(a1, a2, ...)   -- gate: each arg a_i is a directed edge a_i -> id (signal flows
#                              INTO the gate that consumes it), per the task's own instruction.
#   GATETYPE in {AND,NAND,OR,NOR,NOT,XOR,BUFF} across the 10 standard circuits (verified by grep
#   before writing this parser -- no XNOR/DFF/other types present in c17..c6288).
#
# Output .EDGES format matches InputProcessingModule.jl's read_graph_from_edgelist expectations
# (verified by reading that function directly): a header line containing "source"/"destination"
# (any text with those words is accepted, first line skipped), then one "source,destination" CSV
# edge per line, plain Int64 ids, no self-loops, must resolve to a DAG (checked here again
# defensively via a topological pass even though combinational .bench circuits are feedforward by
# construction).
#
# Priors: uniform 0.9 node priors + 0.9 link probabilities (illustrative-uniform convention already
# used for KarlNetwork/bnlearn/the synthetic corpus in this project -- see report for why no
# citable gate-failure-rate figure was used instead).

using Printf

function parse_bench(path::String)
    inputs = Int[]
    outputs = Int[]
    edges = Tuple{Int,Int}[]
    gate_types = Dict{Int,String}()
    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        startswith(line, "#") && continue
        m_in = match(r"^INPUT\((\d+)\)$", line)
        if m_in !== nothing
            push!(inputs, parse(Int, m_in.captures[1]))
            continue
        end
        m_out = match(r"^OUTPUT\((\d+)\)$", line)
        if m_out !== nothing
            push!(outputs, parse(Int, m_out.captures[1]))
            continue
        end
        m_gate = match(r"^(\d+)\s*=\s*([A-Za-z]+)\(([0-9,\s]+)\)$", line)
        if m_gate !== nothing
            gid = parse(Int, m_gate.captures[1])
            gtype = uppercase(m_gate.captures[2])
            args = [parse(Int, strip(a)) for a in split(m_gate.captures[3], ',')]
            gate_types[gid] = gtype
            for a in args
                push!(edges, (a, gid))
            end
            continue
        end
        error("unrecognized line in $path: '$line'")
    end
    (; inputs, outputs, edges, gate_types)
end

function has_cycle_simple(edges::Vector{Tuple{Int,Int}})
    outgoing = Dict{Int,Vector{Int}}()
    nodes = Set{Int}()
    for (u,v) in edges
        push!(nodes, u); push!(nodes, v)
        push!(get!(outgoing, u, Int[]), v)
    end
    WHITE, GRAY, BLACK = 0, 1, 2
    color = Dict(n => WHITE for n in nodes)
    function dfs(u)
        color[u] = GRAY
        for v in get(outgoing, u, Int[])
            if color[v] == GRAY
                return true
            elseif color[v] == WHITE
                dfs(v) && return true
            end
        end
        color[u] = BLACK
        return false
    end
    for n in nodes
        if color[n] == WHITE
            dfs(n) && return true
        end
    end
    return false
end

function convert_circuit(name::String, bench_path::String, out_root::String)
    parsed = parse_bench(bench_path)
    edges = unique(parsed.edges)
    nodes = sort(collect(union(Set(u for (u,_) in edges), Set(v for (_,v) in edges))))

    @assert !has_cycle_simple(edges) "circuit $name: cycle detected -- not a valid DAG"
    # sanity: every declared INPUT that participates in an edge must have no incoming edge (a true
    # source). A handful of ISCAS85 circuits declare an INPUT that is ALSO declared OUTPUT with no
    # gate in between (a direct input->output passthrough pin, e.g. c2670 node 143) -- such a node
    # never appears in any edge at all, so it is simply absent from `nodes`/the converted graph
    # (isolated, no structural contribution); that is expected and not an error.
    incoming_targets = Set(v for (_,v) in edges)
    isolated_passthrough = Int[]
    for i in parsed.inputs
        if i in nodes
            @assert !(i in incoming_targets) "declared INPUT($i) unexpectedly has an incoming edge in $name"
        else
            push!(isolated_passthrough, i)
        end
    end
    if !isempty(isolated_passthrough)
        @info "circuit $name: $(length(isolated_passthrough)) INPUT node(s) never used in any gate (direct input->output passthrough, excluded from graph): $isolated_passthrough"
    end
    # sanity: every declared OUTPUT must either be a node (gate output) or one of the isolated
    # input->output passthroughs identified above.
    for o in parsed.outputs
        @assert (o in nodes) || (o in isolated_passthrough) "declared OUTPUT($o) not found among circuit nodes or isolated passthroughs in $name"
    end

    netdir = joinpath(out_root, "iscas85-$name")
    mkpath(netdir)
    mkpath(joinpath(netdir, "float"))

    edges_path = joinpath(netdir, "iscas85-$name.EDGES")
    open(edges_path, "w") do f
        println(f, "source,destination")
        for (u,v) in sort(edges)
            println(f, "$u,$v")
        end
    end

    np_path = joinpath(netdir, "float", "iscas85-$name-nodepriors.json")
    open(np_path, "w") do f
        println(f, "{")
        println(f, "  \"nodes\": {")
        for (k,n) in enumerate(nodes)
            comma = k < length(nodes) ? "," : ""
            println(f, "    \"$n\": 0.9$comma")
        end
        println(f, "  },")
        println(f, "  \"data_type\": \"Float64\",")
        println(f, "  \"description\": \"Illustrative uniform 0.9 node priors for ISCAS85 $name (no citable per-gate failure-rate figure found in a time-boxed search that maps onto this framework's structural up/down reliability semantics; falls back to the same convention used for KarlNetwork/bnlearn/the synthetic corpus)\"")
        println(f, "}")
    end
    lp_path = joinpath(netdir, "float", "iscas85-$name-linkprobabilities.json")
    open(lp_path, "w") do f
        println(f, "{")
        println(f, "  \"links\": {")
        sedges = sort(edges)
        for (k,(u,v)) in enumerate(sedges)
            comma = k < length(sedges) ? "," : ""
            println(f, "    \"($u,$v)\": 0.9$comma")
        end
        println(f, "  },")
        println(f, "  \"data_type\": \"Float64\",")
        println(f, "  \"description\": \"Illustrative uniform 0.9 link probabilities for ISCAS85 $name (see nodepriors description for rationale)\"")
        println(f, "}")
    end

    n_inputs = length(parsed.inputs)
    n_outputs = length(parsed.outputs)
    n_gates = length(parsed.gate_types)
    gate_type_counts = Dict{String,Int}()
    for (_, t) in parsed.gate_types
        gate_type_counts[t] = get(gate_type_counts, t, 0) + 1
    end
    @printf("%-8s V=%-5d E=%-5d inputs=%-4d outputs=%-4d gates=%-5d types=%s\n",
            name, length(nodes), length(edges), n_inputs, n_outputs, n_gates, gate_type_counts)

    (; V=length(nodes), E=length(edges), n_inputs, n_outputs, n_gates, edges_path, np_path, lp_path)
end

function main()
    REPO = raw"C:\Development\Info_Prop_Framework_Project\Info_Prop_Framework_Project"
    raw_dir = joinpath(REPO, "validation", "probability", "iscas85", "raw")
    out_root = joinpath(REPO, "dag_ntwrk_files")
    circuits = ["c17","c432","c499","c880","c1355","c1908","c2670","c3540","c5315","c6288"]
    results = NamedTuple[]
    for c in circuits
        bench_path = joinpath(raw_dir, "$c.bench")
        isfile(bench_path) || error("missing bench file: $bench_path")
        r = convert_circuit(c, bench_path, out_root)
        push!(results, (name=c, r...))
    end
    csv_path = joinpath(REPO, "validation", "probability", "iscas85", "iscas85_conversion_summary.csv")
    open(csv_path, "w") do f
        println(f, "name,V,E,n_inputs,n_outputs,n_gates")
        for r in results
            @printf(f, "%s,%d,%d,%d,%d,%d\n", r.name, r.V, r.E, r.n_inputs, r.n_outputs, r.n_gates)
        end
    end
    println("wrote ", csv_path)
end

main()
