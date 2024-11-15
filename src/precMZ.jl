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



gpsms = groupby(DataFrame(Arrow.Table(readdir(passing_psms_folder, join=true))),[:ms_file_idx,:precursor_idx])
area_ratios = Float32[]
minimum_weights = Float32[]
for (key, psms) in pairs(gpsms)
    if size(psms, 1) == 2
        row_order = sortperm([last(x) for x in psms[!,:isotopes_captured]])
        if abs(psms[1,:scan_idx] - psms[2,:scan_idx]) > 1
            continue
        end
        push!(area_ratios, psms[row_order[1],:weight]/psms[row_order[2],:weight])
        push!(minimum_weights, minimum(psms[!,:weight]))
    end
end
histogram(log2.(area_ratios), xlim = (-3, 3))
describe(log2.(area_ratios))
plot(log2.(minimum_weights), log2.(area_ratios), ylim = (-3, 3), seriestype=:scatter, alpha = 0.01,
xlabel = "Minimum Isotope Trace Intensity", ylabel = "Log2 Fold Change Between Trace Apex")


test2q = DataFrame(Tables.columntable(Arrow.Table(readdir(second_quant_folder, join=true))))
filter!(x->x.ms_file_idx==1, test2q)
g2q = groupby(test2q, [:ms_file_idx,:precursor_idx])
describe([size(x, 1) for x in g2q])
split_trace = [x for x in g2q if size(x, 1)==2]
area_ratios = zeros(Float32, length(split_trace))
ms_table = Arrow.Table(first(MS_TABLE_PATHS))
mz_offsets = zeros(Float32, length(area_ratios))
for (i, (key, psms)) in enumerate(pairs(split_trace))
    psms_perm = sortperm(psms[!,:isotopes_captured])
    area_ratios[i] = log2(psms[psms_perm[1],:peak_area]/psms[psms_perm[2],:peak_area])
    mz_offsets[i] = precursors[:mz][psms[1,:precursor_idx]] - ms_table[:centerMz][psms[1,:scan_idx]]
end
histogram(area_ratios)
vline!([median(area_ratios)])
plot(mz_offsets, area_ratios, alpha = 0.1, seriestype=:scatter)
if params_[:output_params]["delete_temp"]
    rm(passing_psms_folder,recursive=true)
end

unique_precursors = collect(unique(test2q[!,:precursor_idx]))
prec_mz = [precursors[:mz][pid] for pid in unique_precursors]
prec_charge = [precursors[:prec_charge][pid] for pid in unique_precursors ]
sulfur_count = [precursors[:sulfur_count][pid] for pid in unique_precursors ]
iso_ratios = zeros(Float32, length(unique_precursors))
for i in range(1, length(iso_ratios))
    iso_ratios[i] = iso_splines(min(Int64(sulfur_count[i]), 5), 0, prec_mz[i]*prec_charge[i])
end


median(prec_mz.*prec_charge)

