# One precursor per row, the 6 runs across the columns: does the mobility peak sit in the same place?
# Red dashed = the old +/-32-scan window around that run's best PSM; white dotted = that run's best-PSM IM scan.
# Usage: julia --project=<Pioneer> plot_replicates.jl <out.png> <pid> [pid ...]
using Arrow, DataFrames, Statistics, Printf, Plots
gr()
const STRIDE = 8
const H = expanduser("~/BrukerTims/pride_hye")
outpng = ARGS[1]; pids = parse.(UInt32, ARGS[2:end])
runs = [("A", r) for r in 1:3]; append!(runs, [("B", r) for r in 1:3])
data = Dict{String, Tuple{DataFrame, Dict{UInt32, Int}, Dict{Int, Int}}}()
for (C, R) in runs
    f = "h50_$(C)_REP$R"
    d = DataFrame(Arrow.Table(joinpath(H, "dump_h50_cond$(C)_wide", f * "_chrom_weights.arrow")))
    d = d[in.(d.precursor_idx, Ref(Set(pids))), :]
    p = DataFrame(Arrow.Table(joinpath(H, "out_h50_cond$(C)_wide", "temp_data", "passing_psms", f * ".arrow")))
    p = p[(p.target) .& (p.qval .<= 0.01), :]
    bs = Dict{UInt32, Int}(r.precursor_idx => Int(r.scan_idx) for r in eachrow(p) if r.precursor_idx in pids)
    full = DataFrame(Arrow.Table(joinpath(H, "dump_h50_cond$(C)_wide", f * "_chrom_weights.arrow")))
    sim = Dict{Int, Int}(); for r in eachrow(full); haskey(bs, r.precursor_idx) && (sim[Int(r.scan_idx)] = Int(r.im_scan)); end
    for r in eachrow(full); sim[Int(r.scan_idx)] = Int(r.im_scan); end
    data[f] = (d, bs, sim)
end
panels = []
for pid in pids
    for (C, R) in runs
        f = "h50_$(C)_REP$R"; (d, bs, sim) = data[f]
        sub = d[d.precursor_idx .== pid, :]
        if nrow(sub) == 0
            push!(panels, plot(title = "$pid $f: absent", titlefontsize = 7, framestyle = :box, grid = false, showaxis = false)); continue
        end
        cyc = Int.(sub.cycle_idx); ims = Int.(sub.im_scan)
        c0, c1 = extrema(cyc); i0, i1 = extrema(ims); nr = c1 - c0 + 1; nc = (i1 - i0) ÷ STRIDE + 1
        M = zeros(nr, nc)
        for k in eachindex(cyc); M[cyc[k]-c0+1, (ims[k]-i0) ÷ STRIDE + 1] += Float64(sub.weight[k]); end
        w = Float64.(sub.weight); cen = sum(w .* ims) / sum(w)
        h = heatmap(c0:c1, i0:STRIDE:i1, (M ./ maximum(M))', c = :viridis, clims = (0, 1), colorbar = false,
                    title = @sprintf("%s  %s  centroid %.0f  sum %.1g", pid, f[5:end], cen, sum(w)), titlefontsize = 7,
                    xtickfontsize = 6, ytickfontsize = 6)
        imc = get(sim, get(bs, pid, -1), round(Int, cen))
        hline!(h, [imc - 32, imc + 32], color = :red, ls = :dash, lw = 1, label = "")
        hline!(h, [imc], color = :white, ls = :dot, lw = 1, label = "")
        push!(panels, h)
    end
end
savefig(plot(panels..., layout = (length(pids), 6), size = (1800, 300 * length(pids))), outpng)
println("wrote ", outpng)
