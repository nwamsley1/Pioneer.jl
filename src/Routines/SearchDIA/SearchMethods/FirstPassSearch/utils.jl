"""
    recalibrate_rt_and_rescore!(search_context, params, ms_file_idx, spectra, best_psms, scores; kwargs...)

Per-file RT recalibration after initial LightGBM scoring.

Uses high-confidence PSMs (score > min_prob, target) to refit the RT→iRT spline,
then re-deconvolves, recomputes features (with improved iRT), filters by tighter
iRT tolerance, and re-trains LightGBM.

Returns `(best_psms, scores, q_values, timings, irt_mad)` or `nothing` if
too few calibration PSMs are available (in which case caller keeps original results).
"""
function recalibrate_rt_and_rescore!(
    search_context::SearchContext,
    params::SecondPassSearchParameters,
    ms_file_idx::Int64,
    spectra::MassSpecData,
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

    @info "RT recalibration: fit from $n_calib PSMs, MAD=$(round(mad, digits=3))"

    # 4. Store the improved model (overwrites coarse ParameterTuning model)
    setRtIrtMap!(search_context, model, ms_file_idx)

    # 5. Reload fragment index matches
    frag_match_path = getFragmentIndexMatches(getMSData(search_context), ms_file_idx)
    scan_to_prec_idx, precursors_passed = load_fragment_index_matches(
        frag_match_path, length(spectra)
    )

    # 6. Re-deconvolve with the new RT model
    psms = perform_second_pass_search(
        spectra,
        scan_to_prec_idx,
        precursors_passed,
        search_context,
        params,
        ms_file_idx,
        MS2CHROM();
        n_frag_isotopes = params.prescore_n_frag_isotopes,
        max_frag_rank = params.prescore_max_frag_rank
    )

    n_psms_deconv = nrow(psms)
    if n_psms_deconv == 0
        @warn "RT recalibration: re-deconvolution produced 0 PSMs, skipping"
        return nothing
    end

    # 7. Recompute features (add_features! now uses the new RT spline for irt_obs/irt_error)
    prepare_psm_features!(psms, params, search_context, ms_file_idx, spectra, prescore_only=true)

    if nrow(psms) == 0
        @warn "RT recalibration: no PSMs after feature computation, skipping"
        return nothing
    end

    # 8. Filter by iRT tolerance
    irt_tol = irt_tol_multiplier * mad
    n_before_filter = nrow(psms)
    filter!(row -> row.irt_error <= irt_tol, psms)
    n_after_filter = nrow(psms)
    @info "RT recalibration: iRT filter (tol=$(round(irt_tol, digits=2))): $n_after_filter / $n_before_filter PSMs kept"

    if nrow(psms) == 0
        @warn "RT recalibration: no PSMs after iRT filter, skipping"
        return nothing
    end

    # 9. Re-train LightGBM and select best PSM per precursor
    best_psms2, scores2, q_values2, timings2 = train_lgbm_and_select_best(psms)

    return (best_psms2, scores2, q_values2, timings2, mad)
end
