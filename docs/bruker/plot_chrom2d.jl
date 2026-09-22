# Side-by-side raw / smoothed / baseline-subtracted (+ region) 2D chromatograms for sampled precursors.
# Usage: julia --project=<Pioneer> plot_chrom2d.jl <dump.arrow> <out.png> [lambda_rt lambda_im]
include(joinpath(@__DIR__, "chrom2d.jl"))
using Plots; gr()
dumpf, outpng = ARGS[1], ARGS[2]
λr = length(ARGS) >= 4 ? parse(Float64, ARGS[3]) : 1e-6; λi = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1e-6
d = DataFrame(Arrow.Table(dumpf))
g = groupby(d, :precursor_idx)
tot = [(sub.precursor_idx[1], sum(sub.weight), length(unique(sub.cycle_idx)), length(unique(sub.im_scan))) for sub in g]
cands = sort([t for t in tot if t[3] >= 3 && t[4] >= 3]; by = t -> -t[2])
picks = [cands[max(1, round(Int, q * length(cands)))][1] for q in (0.02, 0.15, 0.4, 0.6, 0.8, 0.95)]
panels = []
for pid in picks
    sub = d[d.precursor_idx .== pid, :]
    cyc = Int.(sub.cycle_idx); ims = Int.(sub.im_scan) .÷ STRIDE
    c0, c1 = extrema(cyc); i0, i1 = extrema(ims); nr = c1 - c0 + 1; nc = i1 - i0 + 1
    Y = zeros(nr, nc); W = zeros(nr, nc)
    for k in eachindex(cyc); Y[cyc[k]-c0+1, ims[k]-i0+1] += Float64(sub.weight[k]); W[cyc[k]-c0+1, ims[k]-i0+1] = 1.0; end
    scale = maximum(Y); Z = wh2d(Y ./ scale, W, λr, λi); B = border_plane(Z); Zb = max.(Z .- B, 0.0)
    apex = argmax(Zb); reg = grow_region(Zb, apex, 0.03)
    xs = c0:c1; ys = (i0:i1) .* STRIDE
    cl = (0, 1)
    h1 = heatmap(xs, ys, (Y ./ scale)', c = :viridis, clims = cl, colorbar = false, title = @sprintf("prec %d raw (max %.2g)", pid, scale), titlefontsize = 8, ylabel = "IM scan")
    h2 = heatmap(xs, ys, Z', c = :viridis, clims = cl, colorbar = false, title = @sprintf("smoothed (λ %.2g / %.2g)", λr, λi), titlefontsize = 8)
    h3 = heatmap(xs, ys, Zb', c = :viridis, clims = cl, colorbar = false, title = @sprintf("baseline plane removed; region %d cells, area %.2g", count(reg), trapezoid2d(Zb, reg, 1.0, 1.0) * scale), titlefontsize = 8, xlabel = "cycle")
    # region outline: mark region cells
    rr = [(xs[r], ys[c]) for r in 1:nr, c in 1:nc if reg[r, c]]
    scatter!(h3, first.(rr), last.(rr), marker = (:square, 2.5, :white, stroke(0)), label = "", alpha = 0.35)
    push!(panels, h1, h2, h3)
end
savefig(plot(panels..., layout = (6, 3), size = (1500, 2100)), outpng)
println("wrote ", outpng)
