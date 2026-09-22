# Three-proteome scatter: log2(A/B) vs log2 mean abundance, per strategy, coloured by species.
# Usage: julia --project=<Pioneer> hye_scatter.jl <tag: h50|h250> <out.png> [lambda]
include(joinpath(@__DIR__, "hye_quant.jl_body"))
using Plots, Printf
gr()

function scatter_panels(tag, outpng, λ)
    sp = species_map()
    files = Dict{String, Vector{String}}(); quant = Dict{String, Dict{UInt32, NamedTuple}}(); pioneer = Dict{String, Dict{UInt32, Float64}}()
    for C in ("A", "B")
        run = "$(tag)_cond$C"
        pl = DataFrame(Arrow.Table(joinpath(H, "out_$run", "precursors_long.arrow")))
        for R in 1:3
            f = "$(tag)_$(C)_REP$R"
            isfile(joinpath(H, "dump_$run", f * "_chrom_weights.arrow")) || continue
            push!(get!(files, C, String[]), f)
            quant[f] = integrate_file(joinpath(H, "dump_$run", f * "_chrom_weights.arrow"), joinpath(H, "out_$run", "temp_data", "passing_psms", f * ".arrow"), λ)
            sub = pl[pl.file_name .== f, :]
            pioneer[f] = Dict{UInt32, Float64}(r.precursor_idx => Float64(r.peak_area) for r in eachrow(sub) if !ismissing(r.peak_area) && r.peak_area > 0)
        end
    end
    common = intersect([Set(keys(quant[f])) for C in ("A", "B") for f in files[C]]...)
    common_p = intersect([Set(keys(pioneer[f])) for C in ("A", "B") for f in files[C]]...)
    cols = Dict("HUMAN" => :steelblue, "YEAST" => :seagreen, "ECOLI" => :indianred)
    strategies = [(:raw, "raw sum, whole window"), (:v8, "raw 5x5 core - flat baseline"), (:v10, "footprint >= 25% apex - flat base"),
                  (:v3, "2D v3 (seeded, descent, rim)"), (:imsum1d, "IM-summed 1D"), (:pioneer, "Pioneer peak_area")]
    panels = []
    for (key, name) in strategies
        pids = key == :pioneer ? common_p : common
        get_v(f, pid) = key == :pioneer ? pioneer[f][pid] : getfield(quant[f][pid], key)
        xs = Dict(s => Float64[] for s in keys(cols)); ys = Dict(s => Float64[] for s in keys(cols))
        for pid in pids
            s = get(sp, pid, ""); haskey(cols, s) || continue
            va = [get_v(f, pid) for f in files["A"]]; vb = [get_v(f, pid) for f in files["B"]]
            (all(>(0), va) && all(>(0), vb)) || continue
            push!(xs[s], log2((mean(va) + mean(vb)) / 2)); push!(ys[s], log2(mean(va) / mean(vb)))
        end
        hc = isempty(ys["HUMAN"]) ? 0.0 : median(ys["HUMAN"])
        p = plot(xlabel = "log2 mean abundance", ylabel = "log2 A/B", title = name, titlefontsize = 9, legend = :topright, legendfontsize = 6, ylims = (-4, 3))
        for s in ("HUMAN", "YEAST", "ECOLI")
            isempty(xs[s]) && continue
            m = median(ys[s]) - hc
            scatter!(p, xs[s], ys[s] .- hc, ms = 1.0, msw = 0, alpha = 0.25, color = cols[s],
                     label = @sprintf("%s n=%d med %.2f (exp %.2f)", s, length(xs[s]), m, EXPECTED[s]))
            hline!(p, [m], color = cols[s], ls = :solid, lw = 1.2, label = "")
        end
        for s in ("HUMAN", "YEAST", "ECOLI"); hline!(p, [EXPECTED[s]], color = cols[s], ls = :dash, lw = 1, label = ""); end
        push!(panels, p)
    end
    savefig(plot(panels..., layout = (3, 2), size = (1400, 1150), left_margin = 6Plots.mm, bottom_margin = 5Plots.mm), outpng)
    println("wrote ", outpng)
end
scatter_panels(ARGS[1], ARGS[2], length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1e-5)
