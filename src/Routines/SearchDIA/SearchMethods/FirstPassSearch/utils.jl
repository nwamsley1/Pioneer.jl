"""
    recalibrate_rt!(search_context, ms_file_idx, best_psms, scores; kwargs...)

Per-file RT recalibration after initial LightGBM scoring.

Uses high-confidence PSMs (score > min_prob, target) to refit the RT→iRT spline
and update iRT error tolerance. SecondPassSearch naturally benefits via
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
        @warn "RT recalibration: only $n_calib high-confidence PSMs (need $min_calib_psms), skipping"
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
        @warn "RT recalibration: fit_irt_model failed ($e), skipping"
        return nothing
    end

    # 4. Store the improved model (overwrites coarse ParameterTuning model)
    old_irt_tol = getIrtErrors(search_context)[ms_file_idx]
    setRtIrtMap!(search_context, model, ms_file_idx)

    # 5. Update iRT error tolerance
    new_irt_tol = mad * irt_tol_multiplier
    getIrtErrors(search_context)[ms_file_idx] = new_irt_tol

    @info "RT recalibration: $n_calib PSMs, MAD=$(round(mad, digits=3)), iRT tol $(round(old_irt_tol, digits=2)) → $(round(new_irt_tol, digits=2))"

    return nothing
end
