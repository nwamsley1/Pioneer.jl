histogram(prec_mzs[(prec_mzs.>398).&(prec_mzs.<410)], xlim=(400, 408), bins = LinRange(398, 410, 3000))


prec_mzs = copy(precursors[:mz])
prec_mzs = prec_mzs[precursors[:prec_charge].==2]
bind_idxs, n_bins = estimateKdeBins(prec_mzs)
x_kde = KernelDensity.kde(prec_mzs, npoints = 60000, bandwidth = 0.01)

kde_eval_points = LinRange(minimum(prec_mzs), maximum(prec_mzs), 60000)
density_atx = KernelDensity.pdf(x_kde, (kde_eval_points))
p = plot(kde_eval_points, density_atx)


local_maxima = findLocalMaxima(density_atx)

x_bin_boundaries = kde_eval_points[local_maxima]
x_bin_boundaries = x_bin_boundaries[density_atx[local_maxima].>0.001]
p = plot(kde_eval_points, density_atx, xlim = (600, 604))
vline!([x_bin_boundaries])

A = hcat(collect(range(1, length(x_bin_boundaries))), 
ones(length(x_bin_boundaries)))
a, b = A\x_bin_boundaries
fitted_bins = [a*x + b for x in range(1, length(x_bin_boundaries))]
plot(fitted_bins, x_bin_boundaries, seriestype=:scatter)


prec_mzs = copy(precursors[:mz])
prec_mzs = prec_mzs[precursors[:prec_charge].==2]
p = histogram(prec_mzs[(prec_mzs.>398).&(prec_mzs.<410)], xlim=(398, 402), bins = LinRange(398, 402, 3000))
vline!(p, x_bin_boundaries)
plot_bins = LinRange(398, 402, 100)
quad_func = getQuadTransmissionFunction(GeneralGaussModel(5.0f0, 0.0f0), 400.0f0, 2.0f0)
plot!(p, plot_bins, quad_func.(plot_bins), lw = 2, alpha = 0.5, show = true)

# Plot the histogram
prec_mzs = copy(precursors[:mz])
prec_mzs = prec_mzs[precursors[:prec_charge].==2]
p = stephist(prec_mzs[(prec_mzs.>398).&(prec_mzs.<410)], xlim=(398, 402), bins = LinRange(398, 402, 3000), label = "Histogram", ylabel = "Count")
prec_mzs = copy(precursors[:mz])
prec_mzs = prec_mzs[precursors[:prec_charge].==3]
stephist!(p, prec_mzs[(prec_mzs.>398).&(prec_mzs.<410)], xlim=(398, 402), bins = LinRange(398, 402, 3000), label = "Histogram", ylabel = "Count")
#vline!(p, x_bin_boundaries)
# Create a twin y-axis for the quadratic function
p2 = twinx(p)
# Plot the quadratic function on the twin y-axis
plot_bins = LinRange(398, 402, 100)
quad_func = getQuadTransmissionFunction(GeneralGaussModel(5.0f0, 0.0f0), 400.0f0, 2.0f0)
plot!(p2, plot_bins, quad_func.(plot_bins), lw = 2, alpha = 0.5, xlim = (398, 402), label = "Quadratic Function", ylabel = "Quadratic Function Value")
vline!(p2, x_bin_boundaries[1:4:end])

window_center_mzs = x_bin_boundaries[1:4:end]
window_center_mzs_a = window_center_mzs[1:2:end]
window_center_mzs_b = window_center_mzs[2:2:end]
CSV.write("/Users/n.t.wamsley/Desktop/AstralAlternatingMethod.csv", DataFrame((center_mz = vcat(window_center_mzs_a, window_center_mzs_b))))
CSV.write("/Users/n.t.wamsley/Desktop/AstralStandardMethod.csv", DataFrame((center_mz = window_center_mzs)))

plot(p, p2, legend = :topleft, show = true)
