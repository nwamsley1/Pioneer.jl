# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

#==========================================================
Per-file MS1 mass-error model for DDA MS1 quantification.

ParameterTuning fits one MS1 model per file (median offset +/- 5 MAD from ~1,200 tuning PSMs) that the DIA MS1
features and, by default, DDA MS1 quantification use. On timsTOF that window is ~+/-20 ppm: much of the spread is
an m/z- and RT-dependent bias and an intensity-dependent centroid precision (+/-25 ppm for weak peaks, ~+/-13 for
strong ones), and the wide window admits co-eluting ions into the precursor's isotope rows.

For DDA MS1 quantification only, the MS1 residuals of the calibration precursors (`ms1_calibration_rows`, the 3,000
most confident targets spread over RT, M0-M2 each) are fitted with the fragment model's machinery,
`fit_intensity_mass_error_model`: m/z + intensity + RT bias splines and an intensity-dependent spread, with a coverage
multiplier for MS1 quantification. The model is used for that file's MS1 extraction (`_ms1_match_peak`); the tuning
model stored in the SearchContext, and so the DIA path, is unchanged.
==========================================================#

"""Target coverage of the MS1-quantification window: the fitter sets k so the window covers this fraction of the
bias-corrected calibration residuals. The same rule as the fragment model (EMPK_QUANTILE, 95%). A narrower window
(85%) gave more accurate timsTOF DDA-PASEF ratios but dropped real signal inconsistently between runs (replicate CV:
a fixed-window sum of the raw signal degrades from 0.05 at +/-15 ppm to 0.19 at +/-5 ppm).
"""
const MS1_QUANT_MODEL_COVERAGE = EMPK_QUANTILE
"""Fewest MS1 residuals for which the intensity-aware model is fitted; below this the tuning model is kept."""
const MS1_QUANT_MODEL_MIN_SAMPLES = 500
"""Most residuals drawn in the QC plot (subsampled for plotting only)."""
const MS1_QUANT_MODEL_PLOT_MAX = 20_000

"""
    ms1_model_coverage(samples, model, fallback) -> NamedTuple

Per-sample bias-corrected residual (ppm), window half-width (ppm) and whether the window covers the theoretical m/z,
for the intensity-aware `model` (exactly the acceptance rule of `_ms1_match_peak`) and for the `fallback` simple model.
"""
function ms1_model_coverage(samples::AbstractVector{MassErrSample}, model::IntensityMassErrorModel,
                            fallback::SimpleMassErrorModel)
    n = length(samples)
    err_ppm = Vector{Float32}(undef, n); hw_ppm = Vector{Float32}(undef, n); covered = falses(n)
    raw_ppm = Vector{Float32}(undef, n)
    off = Float32(getMassOffset(fallback)); tol = Float32(getRightTol(fallback))
    @inbounds for i in 1:n
        s = samples[i]
        c, lo, hi = getCorrectedMzAndBounds(model, s.observed_mz, s.intensity, s.rt)
        err_ppm[i] = (c - s.theoretical_mz) / s.theoretical_mz * 1f6
        hw_ppm[i] = (hi - lo) / 2 / c * 1f6
        covered[i] = lo <= s.theoretical_mz <= hi &&
                     abs(s.observed_mz - (s.theoretical_mz + _mz_bias_da(model, s.theoretical_mz))) <= model.conservative_tol_da
        raw_ppm[i] = (s.observed_mz - s.theoretical_mz) / s.theoretical_mz * 1f6
    end
    fallback_cov = mean(abs.(raw_ppm .- off) .<= tol)
    return (err_ppm = err_ppm, hw_ppm = hw_ppm, covered = covered, raw_ppm = raw_ppm,
            coverage = mean(covered), median_hw = median(hw_ppm), fallback_coverage = fallback_cov,
            fallback_offset = off, fallback_tol = tol)
end

"""
    fit_ms1_quant_mass_model(spectra, calib_psms, search_context, ms_file_idx) -> AbstractMassErrorModel

The MS1 mass-error model for one DDA file's MS1 quantification: `fit_intensity_mass_error_model` on the M0-M2
residuals of `calib_psms` (`collect_ms1_samples`), at the coverage `MS1_QUANT_MODEL_COVERAGE`. Falls back to the tuning
model (`getMs1MassErrorModel`) with fewer than `MS1_QUANT_MODEL_MIN_SAMPLES` residuals, or when the fit fails. Logs the fit and writes a QC plot to
`qc_plots/ms1_quant_mass_error/<file>.pdf`.
"""
function fit_ms1_quant_mass_model(spectra, calib_psms::DataFrame, search_context, ms_file_idx::Integer)
    fallback = getMs1MassErrorModel(search_context, ms_file_idx)
    fallback isa SimpleMassErrorModel || return fallback
    t0 = time()
    samples = collect_ms1_samples(spectra, calib_psms, search_context, ms_file_idx)
    filter!(s -> s.intensity > 0f0 && isfinite(s.rt), samples)
    if length(samples) < MS1_QUANT_MODEL_MIN_SAMPLES
        @user_info "DDA MS1 mass model: $(length(samples)) residuals (< $MS1_QUANT_MODEL_MIN_SAMPLES), keeping the " *
                   "tuning model (offset $(round(getMassOffset(fallback), digits = 2)) ppm, +/-$(round(getRightTol(fallback), digits = 1)) ppm)"
        return fallback
    end
    target = MS1_QUANT_MODEL_COVERAGE
    model = fit_intensity_mass_error_model(samples, fallback; coverage = target)
    if !(model isa IntensityMassErrorModel)
        @user_info "DDA MS1 mass model: fit failed on $(length(samples)) residuals, keeping the tuning model"
        return fallback
    end
    cov = ms1_model_coverage(samples, model, fallback)
    @user_info "DDA MS1 mass model: $(length(samples)) residuals from $(nrow(calib_psms)) precursors, target coverage $(round(100target)) % (k = $(round(model.k, digits = 2))); " *
               "window median +/-$(round(cov.median_hw, digits = 1)) ppm (q10 $(round(quantile(cov.hw_ppm, 0.1), digits = 1)), " *
               "q90 $(round(quantile(cov.hw_ppm, 0.9), digits = 1))), coverage $(round(100cov.coverage, digits = 1))%; " *
               "tuning model +/-$(round(cov.fallback_tol, digits = 1)) ppm covers $(round(100cov.fallback_coverage, digits = 1))%; " *
               "corrected MAD $(round(mad(cov.err_ppm; normalize = true), digits = 2)) ppm (uncorrected " *
               "$(round(mad(cov.raw_ppm; normalize = true), digits = 2))); $(round(time() - t0, digits = 2)) s"
    try
        dir = joinpath(getDataOutDir(search_context), "qc_plots", "ms1_quant_mass_error")
        mkpath(dir)
        fname = getParsedFileName(search_context, ms_file_idx)
        plot_ms1_quant_mass_model(samples, model, cov, fname, joinpath(dir, fname * ".pdf"))
    catch e
        @user_warn "DDA MS1 mass model QC plot failed: $(sprint(showerror, e))"
    end
    return model
end

"""
    plot_ms1_quant_mass_model(samples, model, cov, fname, out_path)

QC for the DDA MS1-quantification mass model, four panels on the calibration residuals (uncorrected ppm):
residual vs m/z with the fitted bias (at the median intensity and RT) and the tuning offset; residual vs RT with
the fitted bias; residual vs log2 intensity with the fitted bias +/- window and the tuning window; and the histogram of
bias-corrected residuals with coverage of the intensity-aware and tuning windows.
"""
function plot_ms1_quant_mass_model(samples::AbstractVector{MassErrSample}, model::IntensityMassErrorModel, cov,
                                   fname::AbstractString, out_path::AbstractString)
    n = length(samples)
    idx = n <= MS1_QUANT_MODEL_PLOT_MAX ? collect(1:n) : sort!(randperm(n)[1:MS1_QUANT_MODEL_PLOT_MAX])
    mz = Float32[samples[i].theoretical_mz for i in idx]; rt = Float32[samples[i].rt for i in idx]
    l2 = Float32[log2(samples[i].intensity) for i in idx]; y = cov.raw_ppm[idx]
    mz0 = median(mz); rt0 = median(rt); l20 = median(l2)
    bias_ppm(m, l, r) = Float32(_total_bias_da(model, m, l, r)) / m * 1f6
    hw_ppm(m, l) = (bounds = getMzBoundsReverse(model, m, l); (bounds[2] - bounds[1]) / 2 / m * 1f6)
    ylim = (-30, 30)
    off, tol = cov.fallback_offset, cov.fallback_tol
    sc = (; markersize = 1.2, markerstrokewidth = 0, alpha = 0.25, color = :steelblue, label = nothing, ylim = ylim)

    gm = range(quantile(mz, 0.005), quantile(mz, 0.995), length = 100)
    p1 = Plots.scatter(mz, y; sc..., xlabel = "m/z", ylabel = "observed - theoretical (ppm)", title = "vs m/z")
    Plots.plot!(p1, gm, [bias_ppm(Float32(m), l20, rt0) for m in gm]; color = :red, lw = 2, label = "fitted bias")
    Plots.hline!(p1, [off]; color = :black, ls = :dash, label = "tuning offset")

    gr_ = range(quantile(rt, 0.005), quantile(rt, 0.995), length = 100)
    p2 = Plots.scatter(rt, y; sc..., xlabel = "RT (min)", title = "vs RT")
    Plots.plot!(p2, gr_, [bias_ppm(mz0, l20, Float32(r)) for r in gr_]; color = :red, lw = 2, label = "fitted bias")
    Plots.hline!(p2, [off]; color = :black, ls = :dash, label = "tuning offset")

    gl = range(quantile(l2, 0.005), quantile(l2, 0.995), length = 100)
    b = [bias_ppm(mz0, Float32(l), rt0) for l in gl]; h = [hw_ppm(mz0, Float32(l)) for l in gl]
    p3 = Plots.scatter(l2, y; sc..., xlabel = "log2 intensity", title = "vs intensity: window")
    Plots.plot!(p3, gl, b; color = :red, lw = 2, label = "fitted bias")
    Plots.plot!(p3, gl, b .+ h; color = :red, ls = :dash, label = "fitted window")
    Plots.plot!(p3, gl, b .- h; color = :red, ls = :dash, label = nothing)
    Plots.hline!(p3, [off - tol, off + tol]; color = :black, ls = :dash, label = "tuning window")

    p4 = Plots.histogram(clamp.(cov.err_ppm, -30f0, 30f0); bins = range(-30, 30, length = 121), color = :seagreen,
                         label = nothing, xlabel = "bias-corrected residual (ppm)", ylabel = "count",
                         title = "corrected: window covers $(round(100cov.coverage, digits = 1))% " *
                                 "(tuning +/-$(round(tol, digits = 1)): $(round(100cov.fallback_coverage, digits = 1))%)")
    Plots.vline!(p4, [-cov.median_hw, cov.median_hw]; color = :red, ls = :dash, label = "median window")

    p = Plots.plot(p1, p2, p3, p4; layout = (2, 2), size = (1300, 950), titlefontsize = 9, guidefontsize = 8,
                   legendfontsize = 7, plot_title = "$fname: DDA MS1 quant mass model (n = $n, k = $(round(model.k, digits = 2)))",
                   plot_titlefontsize = 10, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
    Plots.savefig(p, out_path)
    return p
end
