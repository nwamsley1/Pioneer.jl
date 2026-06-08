# Fused scan loop for MainSearch. Replaces process_scans! when
# `MainSearchParameters.use_fused_scan == true`. Uses `run_fused!` for
# match+score+build in one pass, then runs the same post-design-matrix /
# distance-metric / scoring steps as classic.

@inline function _ms2_cycle_window_key(center_mz, isolation_width)
    center = ismissing(center_mz) ? 0.0 : Float64(center_mz)
    width = ismissing(isolation_width) ? 0.0 : Float64(isolation_width)
    return (round(Int32, center * 1000), round(Int32, width * 1000))
end

@inline _scan_ms2_cycle_window_key(spectra::MassSpecData, scan_idx::Integer) =
    _ms2_cycle_window_key(
        getCenterMz(spectra, scan_idx),
        getIsolationWidthMz(spectra, scan_idx),
    )

@inline function _advance_ms2_cycle_idx(
    cycle_idx::Integer,
    anchor_window_key::Tuple{Int32, Int32},
    has_anchor::Bool,
    ms_order::Integer,
    window_key::Tuple{Int32, Int32},
)
    ms_order == 2 || return Int64(cycle_idx), anchor_window_key, has_anchor
    if !has_anchor
        return Int64(1), window_key, true
    elseif window_key == anchor_window_key
        return Int64(cycle_idx) + 1, anchor_window_key, true
    else
        return Int64(cycle_idx), anchor_window_key, true
    end
end

function _build_scan_cycle_indices(spectra::MassSpecData)
    scan_to_cycle_idx = zeros(UInt32, length(spectra))
    anchor_window_key = (Int32(0), Int32(0))
    has_anchor = false
    cycle_idx = 0

    for scan_idx in 1:length(spectra)
        ms_order = getMsOrder(spectra, scan_idx)
        if ms_order == 2
            cycle_idx, anchor_window_key, has_anchor = _advance_ms2_cycle_idx(
                cycle_idx,
                anchor_window_key,
                has_anchor,
                ms_order,
                _scan_ms2_cycle_window_key(spectra, scan_idx),
            )
            scan_to_cycle_idx[scan_idx] = UInt32(cycle_idx)
        end
    end

    return scan_to_cycle_idx
end

"""
    process_scans_fused!(scan_range, spectra, prec_index,
                         scan_to_cycle_idx, ms_file_idx, search_data,
                         params, precursors, ion_list,
                         nce_model, qtm, mem, rt_to_irt_spline, irt_tol)
    -> DataFrame

MainSearch-only. Same I/O contract as `process_scans!` but uses
`run_fused!` to replace the `selectTransitions! + matchPeaks! + sort +
buildDesignMatrix! + sortSparse! + ScoreFragmentMatches!` chain with a
single per-precursor pass producing a `SparseArrayFused` in CSC order.

The caller must have:
1. Verified library m/z-sortedness (via `verify_mz_sorted`) — once per run.
2. Set `params.use_fused_scan = true` (dispatch decision happens upstream).
"""
function process_scans_fused!(
    scan_range::Vector{Int64},
    spectra::MassSpecData,
    prec_index::PI,
    scan_to_cycle_idx::Vector{UInt32},
    ms_file_idx::Int64,
    search_data::SearchDataStructures,
    params::P,
    precursors::LibraryPrecursors,
    ion_list::LibraryFragmentLookup,
    nce_model::NceModel{Float32},
    qtm::QuadTransmissionModel,
    mem::AbstractMassErrorModel,
    rt_to_irt_spline,
    irt_tol::AbstractFloat
    ;
    deconvolution_solver_override = nothing,
) where {P<:FragmentIndexSearchParameters, PI<:PrecursorIndex}

    Hs              = getHsFused(search_data)
    unscored_psms   = get_unscored_psms(search_data, params)
    id_to_col       = getIdToCol(search_data)
    fused_scratch   = getFusedScratch(search_data)
    corr_mz         = getScanCorrectedMz(search_data)
    obs_low         = getScanObsLow(search_data)
    obs_high        = getScanObsHigh(search_data)
    isotopes_buf    = getIsotopes(search_data)
    prec_trans_buf  = getPrecursorTransmission(search_data)

    prec_mzs     = getMz(precursors)
    prec_charges = getCharge(precursors)
    prec_sulfs   = getSulfurCount(precursors)
    prec_irts    = getIrt(precursors)
    prec_is_decoy = getIsDecoy(precursors)

    prec_estimation     = getPrecEstimation(params)
    n_frag_isotopes     = getNFragIsotopes(params)
    max_frag_rank       = getMaxFragRank(params)
    isotope_err_bounds  = getIsotopeErrBounds(params)
    irt_tol_f32         = Float32(irt_tol)

    # Trait that tags this as a FusedStandard search. Compiler specializes
    # run_fused! on K == FusedStandard{typeof(prec_estimation)}.
    kind = FusedStandard(prec_estimation, UInt8(max_frag_rank))

    last_val  = 0

    for scan_idx in scan_range
        (scan_idx < 1 || scan_idx > length(spectra)) && continue

        msn = getMsOrder(spectra, scan_idx)
        msn ∉ getSpecOrder(params) && continue
        ismissing(get_prec_range(prec_index, scan_idx)) && continue
        cycle_idx = Int64(scan_to_cycle_idx[scan_idx])

        scan_irt = Float32(rt_to_irt_spline(getRetentionTime(spectra, scan_idx)))

        # Pre-compute per-peak (corrected_mz, obs_low, obs_high) once per scan
        # using the 3-arg intensity-aware MEM API.
        scan_mz  = getMzArray(spectra, scan_idx)
        scan_int = getIntensityArray(spectra, scan_idx)
        peak_mz_len = prepare_scan_peaks!(corr_mz, obs_low, obs_high,
                                           mem, scan_mz, scan_int)

        quad_fn = getQuadTransmissionFunction(qtm,
            getCenterMz(spectra, scan_idx),
            getIsolationWidthMz(spectra, scan_idx))

        frag_mz_bounds = (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx))

        prec_range = get_prec_range(prec_index, scan_idx)
        precs_vec  = get_precursors(prec_index)

        # Run deconv pipeline for one prec subset (per-scan slice).
        # Returns updated last_val. The PSMs from this call accumulate in scored_psms.
        function _run_subset!(sub_indices)
            isempty(sub_indices) && return last_val
            nmatches, nmisses = run_fused!(
                kind,
                Hs, unscored_psms, id_to_col, fused_scratch,
                corr_mz, obs_low, obs_high, peak_mz_len,
                isotopes_buf, prec_trans_buf,
                ion_list, nce_model,
                precs_vec, sub_indices,
                prec_mzs, prec_charges, prec_sulfs, prec_irts,
                getIsoSplines(search_data), quad_fn, mem,
                scan_int, scan_irt, irt_tol_f32,
                frag_mz_bounds, n_frag_isotopes,
                isotope_err_bounds;
                m_rank = last(getMinTopNofM(params)),
                scan_idx = Int64(scan_idx)
            )
            if nmatches ≤ 2
                reset_scan_arrays!(id_to_col, Hs, unscored_psms)
                return last_val
            end
            resize_if_needed!(search_data, params)
            converged = if deconvolution_solver_override === nothing
                post_design_matrix!(search_data, Hs, params)
            else
                post_design_matrix!(
                    search_data,
                    Hs,
                    params;
                    deconvolution_solver_override = deconvolution_solver_override,
                )
            end
            if !converged
                reset_scan_arrays!(id_to_col, Hs, unscored_psms)
                return last_val
            end
            compute_distance_metrics!(Hs, search_data, params)
            new_last_val = score_psms!(search_data, params, Hs, scan_idx, nmatches, nmisses,
                                  spectra, last_val, ms_file_idx, cycle_idx; mem=mem)
            reset_scan_arrays!(id_to_col, Hs, unscored_psms)
            return new_last_val
        end

        last_val = _run_subset!(collect(Int64, prec_range))
    end

    return DataFrame(@view(get_scored_psms(search_data, params)[1:last_val]))
end
