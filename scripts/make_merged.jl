# Merge the three measurement CSVs (IPA perf, CUDD naive+sifted complexity) into one paper dataset,
# tag each graph with a family + density, and print a by-family LaTeX summary. Pure CSV work.
const REPO = raw"c:\Users\ohian\OneDrive - University of Strathclyde\Documents\Programmming Files\Julia Files\InformationPropagation\Info_Prop_Framework_Project"
using Printf, Statistics

readcsv(p) = [split(strip(l), ',') for l in eachline(p) if !isempty(strip(l)) && !startswith(strip(l), "#")]
ipa = Dict{String,Vector{String}}(); for r in readcsv(joinpath(REPO,"validation","perf_ipa.csv"))[2:end]; ipa[r[1]] = r; end
cx  = Dict{String,Vector{String}}(); for r in readcsv(joinpath(REPO,"validation","cudd_complexity.csv"))[2:end]; cx[r[1]] = r; end
ops = Dict{String,Vector{String}}(); for r in readcsv(joinpath(REPO,"validation","ops.csv"))[2:end]; ops[r[1]] = r; end

famof(nm) = startswith(nm,"counterexample") ? "counterexample" :
            startswith(nm,"rand28")        ? "mutant-n28" :
            occursin(r"_n(\d+)_", nm)      ? "random-n"*match(r"_n(\d+)_", nm).captures[1] : "other"

names = sort(collect(keys(cx)))
open(joinpath(REPO,"validation","paper_data.csv"),"w") do io
    println(io, "name,family,V,E,vars,density,nroots,nuniq,maxcond,ipa_struct_ms,ipa_prop_ms,ipa_ops,worst_delta,naive_status,naive_nodes,naive_ms,sift_nodes,sift_ms")
    for nm in names
        c = cx[nm]; V=parse(Int,c[2]); E=parse(Int,c[3]); vars=parse(Int,c[4])
        dens = round(E/(V*(V-1)/2), digits=4)
        i = get(ipa, nm, nothing)
        nr = i===nothing ? "" : i[5]; nu = i===nothing ? "" : i[6]; mc = i===nothing ? "" : i[7]
        isms = i===nothing ? "" : i[9]; ipms = i===nothing ? "" : i[10]
        o = get(ops, nm, nothing)
        opc = o===nothing ? "" : o[2]; wd = o===nothing ? "" : o[3]
        @printf(io, "%s,%s,%d,%d,%d,%.4f,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
                nm, famof(nm), V, E, vars, dens, nr, nu, mc, isms, ipms, opc, wd, c[5], c[6], c[7], c[9], c[10])
    end
end

# by-family summary
fams = Dict{String,Vector{String}}(); for nm in names; push!(get!(fams,famof(nm),String[]), nm); end
println("\n% ==== Table 2: corpus summary by family ====")
println("% family & #graphs & V & E & density & med.sift.nodes & max.sift.nodes & naive.blowups \\\\")
for fam in sort(collect(keys(fams)))
    gs = fams[fam]
    Vs = [parse(Int,cx[n][2]) for n in gs]; Es = [parse(Int,cx[n][3]) for n in gs]
    dens = [parse(Int,cx[n][3])/(parse(Int,cx[n][2])*(parse(Int,cx[n][2])-1)/2) for n in gs]
    sn = [parse(Int,cx[n][9]) for n in gs if cx[n][8]=="ok"]
    nb = count(n->cx[n][5]=="blowup", gs)
    @printf("%-16s & %d & %d–%d & %d–%d & %.2f–%.2f & %d & %d & %d \\\\\n",
            fam, length(gs), minimum(Vs),maximum(Vs), minimum(Es),maximum(Es),
            minimum(dens),maximum(dens), round(Int,median(sn)), maximum(sn), nb)
end

# global headline numbers
allsn = [parse(Int,cx[n][9]) for n in names if cx[n][8]=="ok"]
naive_ok = [parse(Int,cx[n][6]) for n in names if cx[n][5]=="ok"]
ipap = [parse(Float64,ipa[n][10]) for n in names if haskey(ipa,n)]
@printf("\n%% headline: graphs=%d | sift nodes med=%d max=%d (all build) | naive blowups=%d (>2.5M) | IPA prop med=%.1fms max=%.1fms\n",
        length(names), round(Int,median(allsn)), maximum(allsn),
        count(n->cx[n][5]=="blowup", names), median(ipap), maximum(ipap))
println("wrote validation/paper_data.csv")
