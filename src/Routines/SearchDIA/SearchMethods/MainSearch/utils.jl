"""
    recalibrate_rt!(search_context, ms_file_idx, best_psms, scores; kwargs...)

Per-file RT recalibration after initial LightGBM scoring.

Uses high-confidence PSMs (score > min_prob, target) to refit the RT→iRT spline
and update iRT error tolerance. Downstream steps naturally benefit via
`getRtIrtModel()` and `getIrtErrors()`.
"""
function recalibrate_rt!(
    search_context::SearchContext,
    ms_file_idx::Int64,
    best_psms::DataFrame,
    scores::Vector{Float32};
    min_prob::Float32 = 0.9f0,
    irt_tol_multiplier::Float32 = 4.0f0,
    min_calib_psms::Int = 30
)
    # 1. Filter to high-confidence target PSMs for calibration
    calib_mask = (scores .> min_prob) .& best_psms[!, :target]
    n_calib = count(calib_mask)

    if n_calib < min_calib_psms
        @user_warn "RT recalibration: only $n_calib high-confidence PSMs (need $min_calib_psms), skipping"
        return nothing
    end

    # 2. Build calibration DataFrame with columns expected by fit_irt_model
    calib_df = DataFrame(
        rt = best_psms[calib_mask, :rt],
        irt_predicted = best_psms[calib_mask, :irt_pred]
    )

    # 3. Fit new iRT model
    local model, rts, irts, mad
    try
        model, rts, irts, mad = fit_irt_model(calib_df)
    catch e
        @user_warn "RT recalibration: fit_irt_model failed ($e), skipping"
        return nothing
    end

    # 4. Store the improved model (overwrites coarse ParameterTuning model)
    setRtIrtMap!(search_context, model, ms_file_idx)

    # 5. Update iRT error tolerance
    new_irt_tol = mad * irt_tol_multiplier
    getIrtErrors(search_context)[ms_file_idx] = new_irt_tol

    @debug_l1 "RT spline fit: $n_calib high-prob target PSMs (prob > $min_prob), MAD=$(round(mad, digits=3)), new iRT tol=$(round(new_irt_tol, digits=2))"

    return nothing
end

"""
    fit_im_lines(scan, pred, charge, calib; min_calib=100)

Per-charge ion-mobility calibration lines. For the calibration rows (`calib`), fits
library-predicted 1/K0 = a + b * packet IM scan index by least squares, once pooled
over all charges (key 0) and once per charge with at least `min_calib` rows, with
sigma = 1.4826 * MAD of the residuals (floored at 1e-4). Returns a
`Dict{Int, NTuple{3,Float32}}` of (a, b, sigma); empty when the pooled set is too small.
"""
function fit_im_lines(
    scan::AbstractVector{Float32},
    pred::AbstractVector{Float32},
    charge::AbstractVector,
    calib::AbstractVector{Bool};
    min_calib::Int = 100
)
    models = Dict{Int, NTuple{3, Float32}}()
    function fit(idx)
        x = scan[idx]; y = pred[idx]
        xm = mean(x); ym = mean(y)
        vx = sum((x .- xm) .^ 2)
        b = vx > 0 ? sum((x .- xm) .* (y .- ym)) / vx : 0f0
        a = ym - b * xm
        r = y .- (a .+ b .* x)
        s = 1.4826f0 * median(abs.(r .- median(r)))
        return (Float32(a), Float32(b), max(Float32(s), 1f-4))
    end
    all_idx = findall(calib)
    length(all_idx) < min_calib && return models
    models[0] = fit(all_idx)
    for z in unique(charge[all_idx])
        idx = findall(i -> calib[i] && charge[i] == z, eachindex(calib))
        length(idx) >= min_calib && (models[Int(z)] = fit(idx))
    end
    return models
end

"""
    add_im_error!(best_psms, scores, spectra, precursors, ms_file_idx; min_prob=0.9, min_calib=100)

Per-file ion-mobility calibration after initial LightGBM scoring (ion-mobility packet
data). Fits per-charge lines of library-predicted 1/K0 against the packet's IM scan
index on high-confidence target PSMs (score > `min_prob`) via `fit_im_lines`, and writes
`im_error` = |predicted 1/K0 - line(scan)| / sigma for every row. Charges with fewer
than `min_calib` calibration PSMs use the pooled line. When the file has no IM scan
column or the library no mobility predictions, `im_error` is 0 everywhere — the column
always exists because ScoringSearch takes its feature list from the first file's schema.
Returns the fitted lines (`Dict{Int, NTuple{3,Float32}}`, empty when nothing was fit) so
the caller can store them on the SearchContext for downstream stages.
"""
function add_im_error!(
    best_psms::DataFrame,
    scores::Vector{Float32},
    spectra::MassSpecData,
    precursors,
    ms_file_idx::Int64;
    min_prob::Float32 = 0.9f0,
    min_calib::Int = 100
)
    n = nrow(best_psms)
    im_error = zeros(Float32, n)
    im_scans = getImScans(spectra)
    im_lib = getInvIonMobility(precursors)
    if im_scans === nothing || im_lib === nothing || n == 0
        best_psms[!, :im_error] = im_error
        return Dict{Int, NTuple{3, Float32}}()
    end
    scan = Float32[Float32(im_scans[si]) for si in best_psms[!, :scan_idx]]
    pred = Float32[Float32(im_lib[pid]) for pid in best_psms[!, :precursor_idx]]
    charge = best_psms[!, :charge]
    calib = (scores .> min_prob) .& best_psms[!, :target]
    models = fit_im_lines(scan, pred, charge, calib; min_calib = min_calib)
    if isempty(models)
        @debug_l1 "IM calibration (file $ms_file_idx): fewer than $min_calib high-confidence PSMs, im_error = 0"
        best_psms[!, :im_error] = im_error
        return models
    end
    pooled = models[0]
    @inbounds for i in 1:n
        a, b, s = get(models, Int(charge[i]), pooled)
        im_error[i] = abs(pred[i] - (a + b * scan[i])) / s
    end
    best_psms[!, :im_error] = im_error
    for (z, (a, b, s)) in sort(collect(models))
        n_z = z == 0 ? count(calib) : count(i -> calib[i] && Int(charge[i]) == z, eachindex(calib))
        @debug_l1 "  IM line " * (z == 0 ? "pooled" : "z=$z") * " (file $ms_file_idx): pred 1/K0 = " *
                  "$(round(a, digits=4)) + ($(round(b, digits=6))) * scan, sigma = $(round(s, digits=4)), n_calib = $n_z"
    end
    return models
end

"""
    plot_im_calibration(best_psms, scores, spectra, precursors, models, fname; min_prob=0.9)

QC plots for the per-file ion-mobility calibration (`add_im_error!`): page 1 scatters the
calibration PSMs (score > `min_prob`, target) as packet IM scan vs library 1/K0 per charge
with the fitted lines and +/- 3 sigma bands; page 2 overlays per-charge residual histograms
of the calibration targets and of all decoys. Returns a vector of two plots.
"""
function plot_im_calibration(
    best_psms::DataFrame,
    scores::AbstractVector{Float32},
    spectra::MassSpecData,
    precursors,
    models::Dict{Int, NTuple{3, Float32}},
    fname::AbstractString;
    min_prob::Float32 = 0.9f0
)
    im_scans = getImScans(spectra)
    im_lib = getInvIonMobility(precursors)
    scan = Float64[im_scans[si] for si in best_psms[!, :scan_idx]]
    pred = Float64[im_lib[pid] for pid in best_psms[!, :precursor_idx]]
    charge = Int.(best_psms[!, :charge])
    target = best_psms[!, :target]
    calib = (scores .> min_prob) .& target
    pooled = models[0]
    colors = Dict(1 => :gray, 2 => :steelblue, 3 => :darkorange, 4 => :seagreen, 5 => :purple)
    zs = sort(unique(charge[calib]))
    xs = range(minimum(scan), maximum(scan); length = 100)

    p1 = plot(xlabel = "packet IM scan", ylabel = "library 1/K0",
              title = "$fname\nIM calibration on targets with prob > $min_prob (dashed: +/- 3 sigma)",
              titlefontsize = 9, legend = :topright, legendfontsize = 7)
    for z in zs
        idx = findall(calib .& (charge .== z))
        a, b, s = get(models, z, pooled)
        col = get(colors, z, :black)
        own = haskey(models, z) ? "" : " (pooled line)"
        scatter!(p1, scan[idx], pred[idx], ms = 1.5, ma = 0.25, msw = 0, color = col,
                 label = "z=$z n=$(length(idx))$own")
        plot!(p1, xs, a .+ b .* xs, color = col, lw = 2,
              label = "z=$z: $(round(a, digits = 4)) + ($(round(b, digits = 6)))*scan, sigma $(round(s, digits = 4))")
        plot!(p1, xs, a .+ b .* xs .+ 3s, color = col, ls = :dash, lw = 1, label = "")
        plot!(p1, xs, a .+ b .* xs .- 3s, color = col, ls = :dash, lw = 1, label = "")
    end

    panels = Plots.Plot[]
    for z in zs
        a, b, s = get(models, z, pooled)
        resid(idx) = pred[idx] .- (a .+ b .* scan[idx])
        r_t = resid(findall(calib .& (charge .== z)))
        r_d = resid(findall(.!target .& (charge .== z)))
        lim = 6s
        p = plot(xlabel = "library 1/K0 - line (1/K0)", ylabel = "density", title = "z=$z",
                 titlefontsize = 9, legendfontsize = 7, xlims = (-lim, lim))
        histogram!(p, clamp.(r_t, -lim, lim), bins = 60, normalize = :pdf, alpha = 0.6, color = :steelblue,
                   lc = :match, label = "targets prob > $min_prob (n=$(length(r_t)))")
        isempty(r_d) || histogram!(p, clamp.(r_d, -lim, lim), bins = 60, normalize = :pdf, alpha = 0.5,
                                   color = :firebrick, lc = :match, label = "decoys, all (n=$(length(r_d)))")
        vline!(p, [-3s, 3s], color = :black, ls = :dash, label = "")
        push!(panels, p)
    end
    p2 = plot(panels..., layout = (1, length(panels)), size = (420 * length(panels), 380),
              plot_title = "$fname: IM residuals", plot_titlefontsize = 9,
              left_margin = 8 * Plots.mm, bottom_margin = 6 * Plots.mm)
    return Plots.Plot[p1, p2]
end

"""
    compute_rt_binned_tolerance!(search_context, rt_binned_tol, ms_data, n_files)

Store an RTBinnedTolerance for each non-failed file.
"""
function compute_rt_binned_tolerance!(
    search_context::SearchContext,
    rt_binned_tol::RTBinnedTolerance,
    ms_data,
    n_files::Int
)
    for ms_file_idx in 1:n_files
        getRtTolerances(search_context)[ms_file_idx] = rt_binned_tol
    end
end
