# Fused scan loop for MainSearch. Replaces process_scans! when
# `MainSearchParameters.use_fused_scan == true`. Uses `run_fused!` for
# match+score+build in one pass, then runs the same post-design-matrix /
# distance-metric / scoring steps as classic.

"""
    process_scans_fused!(scan_range, spectra, prec_index,
                         search_data, params, precursors, ion_list,
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
# Set is_center[col] for each active design-matrix column: true iff the precursor
# at that column has its monoisotope m/z within the scan's center bin
# (|prec_mz − centerMz| ≤ isoWidth/2), matching filter_to_center_bin!. Columns not
# referenced by a candidate default to false (they carry no weight / aren't emitted).
@inline function fill_is_center!(is_center::Vector{Bool}, ncol::Int, id_to_col,
                                 precs, prec_range, prec_mzs::AbstractVector{Float32},
                                 cmz::Float32, hw::Float32)
    @inbounds for c in 1:ncol
        is_center[c] = false
    end
    @inbounds for ri in prec_range
        p = precs[ri]
        c = id_to_col[p]
        c == 0 && continue
        is_center[c] = abs(prec_mzs[p] - cmz) <= hw
    end
    return nothing
end

function process_scans_fused!(
    scan_range::Vector{Int64},
    spectra::MassSpecData,
    prec_index::PI,
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
    cycle_idxs = getCycleIdxs(spectra)

    # Sciex ZT lightweight neighbors: only the center bin of each (precursor,
    # cycle) meta-scan needs full spectral scoring; neighbor bins are kept only
    # for their weight. Gate scoring on a per-column center mask (reused buffer).
    _zt_lightweight = (params isa MainSearchParameters) &&
        something(tryparse(Int, get(ENV, "PIONEER_ZT_METASCAN_K", "")), 0) > 0
    is_center_buf = Bool[]

    for scan_idx in scan_range
        (scan_idx < 1 || scan_idx > length(spectra)) && continue

        msn = getMsOrder(spectra, scan_idx)
        msn ∉ getSpecOrder(params) && continue
        ismissing(get_prec_range(prec_index, scan_idx)) && continue
        cycle_idx = Int64(cycle_idxs[scan_idx])

        scan_rt = Float32(getRetentionTime(spectra, scan_idx))
        scan_irt = Float32(rt_to_irt_spline(scan_rt))

        # Pre-compute per-peak (corrected_mz, obs_low, obs_high) once per scan
        # using the intensity/RT-aware MEM API.
        scan_mz  = getMzArray(spectra, scan_idx)
        scan_int = getIntensityArray(spectra, scan_idx)
        peak_mz_len = prepare_scan_peaks!(corr_mz, obs_low, obs_high,
                                           mem, scan_mz, scan_int, scan_rt)

        quad_fn = getQuadTransmissionFunction(qtm,
            getCenterMz(spectra, scan_idx),
            getIsolationWidthMz(spectra, scan_idx))

        frag_mz_bounds = (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx))

        prec_range = get_prec_range(prec_index, scan_idx)
        precs_vec  = get_precursors(prec_index)

        # Run deconv pipeline for this scan's precursor subset. Inlined: this was
        # previously a closure (`_run_subset!`) defined inside the loop and called
        # once per scan, which allocated a capture object every iteration and boxed
        # the reassigned `last_val`. Inlining preserves the exact control flow.
        # Pass the precursor range directly. `prec_range` is a UnitRange{Int64}
        # (non-missing guaranteed by the ismissing guard above); run_fused! only
        # iterates it, so the previous `collect(Int64, ...)` per-scan Vector was
        # pure waste (~0.2 GB/file). Type-assert narrows the Union for dispatch.
        sub_indices = prec_range::UnitRange{Int64}
        if !isempty(sub_indices)
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
            else
                resize_if_needed!(search_data, params)
                converged = post_design_matrix!(search_data, Hs, params)
                if !converged
                    reset_scan_arrays!(id_to_col, Hs, unscored_psms)
                else
                    if _zt_lightweight
                        _cmz = getCenterMz(spectra, scan_idx)
                        _hw = getIsolationWidthMz(spectra, scan_idx)
                        if !ismissing(_cmz) && !ismissing(_hw)
                            ncol = Hs.n
                            ncol > length(is_center_buf) && resize!(is_center_buf, ncol)
                            fill_is_center!(is_center_buf, ncol, id_to_col, precs_vec,
                                            prec_range, prec_mzs, Float32(_cmz), Float32(_hw)/2)
                            compute_distance_metrics!(Hs, search_data, params;
                                is_center = view(is_center_buf, 1:ncol))
                        else
                            compute_distance_metrics!(Hs, search_data, params)
                        end
                    else
                        compute_distance_metrics!(Hs, search_data, params)
                    end
                    last_val = score_psms!(search_data, params, Hs, scan_idx, nmatches, nmisses,
                                          spectra, last_val, ms_file_idx, cycle_idx; mem=mem)
                    reset_scan_arrays!(id_to_col, Hs, unscored_psms)
                end
            end
        end
    end

    return DataFrame(@view(get_scored_psms(search_data, params)[1:last_val]))
end
