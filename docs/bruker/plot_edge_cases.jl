# 2D chromatograms of apex-on-IM-edge cases vs well-behaved cases, from a WIDE (64-scan) dump.
# The old 32-scan window is drawn as a box around the precursor's best-PSM IM scan, so what the narrow
# window would have missed is visible.
# Usage: julia --project=<Pioneer> plot_edge_cases.jl <dump.arrow> <passing_psms.arrow> <out.png>
using Arrow, DataFrames, Statistics, Printf, Plots
gr()
const STRIDE = 8
dumpf, psmsf, outpng = ARGS[1], ARGS[2], ARGS[3]
d = DataFrame(Arrow.Table(dumpf))
p = DataFrame(Arrow.Table(psmsf)); p = p[(p.target) .& (p.qval .<= 0.01), :]
best_scan = Dict{UInt32, Int}(r.precursor_idx => Int(r.scan_idx) for r in eachrow(p))
scan_im = Dict{Int, Int}(); for r in eachrow(d); scan_im[Int(r.scan_idx)] = Int(r.im_scan); end
edge = Tuple{UInt32, Float64}[]; good = Tuple{UInt32, Float64}[]
for sub in groupby(d, :precursor_idx)
    w = Float64.(sub.weight); s = sum(w); s > 0 || continue
    ims = Int.(sub.im_scan); i0, i1 = extrema(ims); i1 > i0 || continue
    ai = ims[argmax(w)]
    ring = maximum(vcat(0.0, w[(ims .== i0) .| (ims .== i1)])) / maximum(w)
    if ai == i0 || ai == i1 || ring > 0.5
        push!(edge, (sub.precursor_idx[1], s))
    elseif ring < 0.1
        push!(good, (sub.precursor_idx[1], s))
    end
end
sort!(edge; by = x -> -x[2]); sort!(good; by = x -> -x[2])
picks = vcat([(e[1], "apex on edge / high rim") for e in edge[[1, max(1, length(edge) ÷ 20), max(1, length(edge) ÷ 4)]]],
             [(g[1], "well behaved") for g in good[[1, max(1, length(good) ÷ 20), max(1, length(good) ÷ 4)]]])
panels = []
for (pid, lab) in picks
    sub = d[d.precursor_idx .== pid, :]
    cyc = Int.(sub.cycle_idx); ims = Int.(sub.im_scan)
    c0, c1 = extrema(cyc); i0, i1 = extrema(ims); nr = c1 - c0 + 1; nc = (i1 - i0) ÷ STRIDE + 1
    M = zeros(nr, nc)
    for k in eachindex(cyc); M[cyc[k]-c0+1, (ims[k]-i0) ÷ STRIDE + 1] += Float64(sub.weight[k]); end
    xs = c0:c1; ys = i0:STRIDE:i1
    h = heatmap(xs, ys, (M ./ maximum(M))', c = :viridis, clims = (0, 1), colorbar = false,
                title = @sprintf("prec %d — %s (sum %.2g)", pid, lab, sum(sub.weight)), titlefontsize = 8,
                xlabel = "cycle", ylabel = "IM scan")
    # the old +/- 32-scan window around the best PSM's IM scan
    bs = get(best_scan, pid, 0); im_c = get(scan_im, bs, ims[argmax(sub.weight)])
    hline!(h, [im_c - 32, im_c + 32], color = :red, ls = :dash, lw = 1.5, label = "")
    hline!(h, [im_c], color = :white, ls = :dot, lw = 1, label = "")
    push!(panels, h)
end
savefig(plot(panels..., layout = (3, 2), size = (1300, 1100)), outpng)
println("wrote ", outpng)
