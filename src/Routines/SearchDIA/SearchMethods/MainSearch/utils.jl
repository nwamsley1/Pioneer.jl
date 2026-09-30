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
    fit_im_line(scan, pred, calib; min_calib=100)

Ion-mobility calibration line on the rows flagged by `calib`: library-predicted
1/K0 = a + b * packet IM scan index by least squares, with sigma = 1.4826 * MAD of the
residuals (floored at 1e-4). Returns `(a, b, sigma)`, or `nothing` with fewer than
`min_calib` rows. Callers pass the z2 rows: one line serves every charge, because the
library predictor's z3 / z4 values sit on the z2 line with wider scatter (see
IM_GATE_SIGMA_MULT), and z2 is the only charge with enough PSMs in every file.
"""
function fit_im_line(
    scan::AbstractVector{Float32},
    pred::AbstractVector{Float32},
    calib::AbstractVector{Bool};
    min_calib::Int = 100
)
    idx = findall(calib)
    length(idx) < min_calib && return nothing
    x = scan[idx]; y = pred[idx]
    xm = mean(x); ym = mean(y)
    vx = sum((x .- xm) .^ 2)
    b = vx > 0 ? sum((x .- xm) .* (y .- ym)) / vx : 0f0
    a = ym - b * xm
    r = y .- (a .+ b .* x)
    s = 1.4826f0 * median(abs.(r .- median(r)))
    return (Float32(a), Float32(b), max(Float32(s), 1f-4))
end

"""
    fill_missing_im_lines!(search_context, n_files)

After ParameterTuning: give every file without its own z2 line (fewer than
TUNING_IM_MIN_CALIB z2 PSMs, or tuning failed) the median of the other files' lines, so
the IM gate, `im_error` and `im_obs` exist for every file. When every file has an
instrument calibration (`.tdfs`), the median is taken in instrument 1/K0 units at two
anchor mobilities and mapped back to each receiving file's IM scan scale, so files with
different scan ranges combine correctly; otherwise it is taken in scan units. O(n_files),
no file I/O. Does nothing when no file has a line (no gate; `im_error` 0 everywhere).
"""
fill_missing_im_lines!(search_context::SearchContext, n_files::Integer) =
    fill_missing_im_lines!(search_context.im_models, search_context.im_cals, n_files)

function fill_missing_im_lines!(
    models::Dict{Int64, Dict{Int, NTuple{3, Float32}}},
    cals::Dict{Int64, NTuple{2, Float32}},
    n_files::Integer
)
    has_line(i) = haskey(get(models, Int64(i), Dict{Int, NTuple{3, Float32}}()), 2)
    donors = [i for i in 1:n_files if has_line(i)]
    receivers = [i for i in 1:n_files if !has_line(i)]
    (isempty(donors) || isempty(receivers)) && return nothing
    if all(i -> haskey(cals, Int64(i)), 1:n_files)
        k1, k2 = 0.8f0, 1.2f0                         # anchor 1/K0 values
        # line pred = a + b * scan, instrument k0 = c0 + m * scan  =>  pred at k0 = a + b * (k0 - c0) / m
        at(i, k) = ((a, b, _) = models[i][2]; (c0, m) = cals[i]; a + b * (k - c0) / m)
        y1 = median(at(i, k1) for i in donors)
        y2 = median(at(i, k2) for i in donors)
        B = (y2 - y1) / (k2 - k1); A = y1 - B * k1  # median line: pred = A + B * k0
        s = median(models[i][2][3] for i in donors)
        for i in receivers
            c0, m = cals[i]
            models[i] = Dict(2 => (Float32(A + B * c0), Float32(B * m), Float32(s)))
        end
    else
        med(j) = Float32(median(models[i][2][j] for i in donors))
        line = (med(1), med(2), med(3))
        for i in receivers
            models[i] = Dict(2 => line)
        end
    end
    @user_info "IM calibration: $(length(receivers)) file(s) without their own z2 line use the median of $(length(donors)) other file(s)"
    return nothing
end

"""
    add_im_error!(best_psms, scores, spectra, precursors, ms_file_idx; fallback=nothing, min_prob=0.9, min_calib=100)

Per-file ion-mobility calibration after initial LightGBM scoring (ion-mobility packet
data). Refits the z2 line of library-predicted 1/K0 against the packet's IM scan index on
high-confidence z2 target PSMs (score > `min_prob`, `fit_im_line`); with fewer than
`min_calib` it uses `fallback` (the file's ParameterTuning line, or the median line from
`fill_missing_im_lines!`). Writes, for every row and every charge,
`im_error` = (predicted 1/K0 - line(scan)) / sigma, SIGNED and in z2 sigma units (the
scoring model also sees :charge, so it learns each charge's offset and spread around the
line), and `im_obs` = line(scan), the observed mobility on the library's 1/K0 scale
(comparable between runs). Without an IM scan column or library mobility both are 0;
on mobility data with no line at all, `im_error` is 0 and `im_obs` NaN (MBR reads
non-finite as missing). Either way the value is the same for every file of the search;
the columns always exist because ScoringSearch takes its feature list from the first
file's schema. Returns `Dict(2 => line)`, empty when there is none, for the SearchContext.
"""
function add_im_error!(
    best_psms::DataFrame,
    scores::Vector{Float32},
    spectra::MassSpecData,
    precursors,
    ms_file_idx::Int64;
    fallback::Union{Nothing, NTuple{3, Float32}} = nothing,
    min_prob::Float32 = 0.9f0,
    min_calib::Int = 100
)
    n = nrow(best_psms)
    im_error = zeros(Float32, n)
    # Observed mobility on the library's 1/K0 scale, i.e. the calibration line evaluated at the PSM's
    # IM scan. Unlike im_error (a sigma-normalised residual) this is comparable BETWEEN runs, which is
    # what a donor/receiver mobility comparison needs. Zero on data without ion mobility (as before);
    # NaN on mobility data with no line, which MBR reads as missing.
    im_obs = zeros(Float32, n)
    im_scans = getImScans(spectra)
    im_lib = getInvIonMobility(precursors)
    if im_scans === nothing || im_lib === nothing || n == 0
        best_psms[!, :im_error] = im_error
        best_psms[!, :im_obs] = im_obs
        return Dict{Int, NTuple{3, Float32}}()
    end
    scan = Float32[Float32(im_scans[si]) for si in best_psms[!, :scan_idx]]
    pred = Float32[Float32(im_lib[pid]) for pid in best_psms[!, :precursor_idx]]
    charge = best_psms[!, :charge]
    calib = (scores .> min_prob) .& best_psms[!, :target] .& (charge .== 2)
    own = fit_im_line(scan, pred, calib; min_calib = min_calib)
    line = own !== nothing ? own : fallback
    if line === nothing
        @debug_l1 "IM calibration (file $ms_file_idx): no z2 line (fewer than $min_calib high-confidence z2 PSMs, no fallback), im_error = 0"
        best_psms[!, :im_error] = im_error
        best_psms[!, :im_obs] = fill(NaN32, n)
        return Dict{Int, NTuple{3, Float32}}()
    end
    a, b, s = line
    @inbounds for i in 1:n
        obs = a + b * scan[i]
        im_obs[i] = obs
        im_error[i] = (pred[i] - obs) / s
    end
    best_psms[!, :im_error] = im_error
    best_psms[!, :im_obs] = im_obs
    @debug_l1 "  IM line z=2 (file $ms_file_idx, " * (own !== nothing ? "own" : "fallback") * "): pred 1/K0 = " *
              "$(round(a, digits=4)) + ($(round(b, digits=6))) * scan, sigma = $(round(s, digits=4)), n_calib = $(count(calib))"
    return Dict(2 => line)
end

"""
    plot_im_calibration(best_psms, scores, spectra, precursors, models, fname; min_prob=0.9)

QC plots for the per-file ion-mobility calibration (`add_im_error!`): page 1 scatters the
calibration PSMs (score > `min_prob`, target) as packet IM scan vs library 1/K0 per charge
with the z2 line used for every charge and its +/- 3 sigma band; page 2 overlays per-charge
residual histograms (around that line) of the calibration targets and of all decoys.
Returns a vector of two plots.
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
    a, b, s = models[2]
    colors = Dict(1 => :gray, 2 => :steelblue, 3 => :darkorange, 4 => :seagreen, 5 => :purple)
    zs = sort(unique(charge[calib]))
    xs = range(minimum(scan), maximum(scan); length = 100)

    p1 = plot(xlabel = "packet IM scan", ylabel = "library 1/K0",
              title = "$fname\nIM calibration on targets with prob > $min_prob (dashed: +/- 3 sigma)",
              titlefontsize = 9, legend = :topright, legendfontsize = 7)
    for z in zs
        idx = findall(calib .& (charge .== z))
        scatter!(p1, scan[idx], pred[idx], ms = 1.5, ma = 0.25, msw = 0, color = get(colors, z, :black),
                 label = "z=$z n=$(length(idx))")
    end
    # the one z2 line used for every charge
    plot!(p1, xs, a .+ b .* xs, color = :black, lw = 2,
          label = "z2 line: $(round(a, digits = 4)) + ($(round(b, digits = 6)))*scan, sigma $(round(s, digits = 4))")
    plot!(p1, xs, a .+ b .* xs .+ 3s, color = :black, ls = :dash, lw = 1, label = "")
    plot!(p1, xs, a .+ b .* xs .- 3s, color = :black, ls = :dash, lw = 1, label = "")

    panels = Plots.Plot[]
    for z in zs
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
