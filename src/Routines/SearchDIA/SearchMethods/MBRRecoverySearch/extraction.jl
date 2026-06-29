# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# Per-cell PMM extraction for MBRRecoverySearch.
#
# For one (donor_pid × receiver_file) cell, projects the donor's irt_obs
# through the receiver file's rt_to_irt spline to get a center RT, opens a
# calibrated rt_binned_tol window around that center, and runs the same
# per-scan pipeline IntegrateChromatogramsSearch uses
# (collect_rt_window_precursors! → run_fused! → solve_deconvolution!):
# the only difference is that the donor pid is manually prepended to the
# candidate list so it gets a column in the design matrix even though it's
# absent from `precursors_passing`.
#
# Returns the max-weight scan's features for use in inline MBR compute,
# or `nothing` if no scan in the window produced positive weight for the
# donor pid.

"""
    extract_max_weight_in_window!(scratch, ...) -> Union{Nothing, ReceiverExtraction}

Inputs:
  scratch              SearchDataStructures with Hs, weights, etc. (thread-local)
  spectra              MassSpecData for the receiver file
  donor_pid            UInt32 pid to extract
  precursors_passing   Set{UInt32} of already-IDed pids in this file (used as
                       allowlist for COMPETITORS; donor is prepended manually)
  rt_index             retentionTimeIndex built from already-IDed pids
  center_rt            Float32 center RT (donor irt_obs projected via rt_irt_model)
  rt_tol               Float32 half-width (from rt_binned_tol or floor)
  rt_irt_model         per-file rt_to_irt spline (for output irt_obs)
  mass_error_model, quad_model, nce_model    per-file calibrated models
  precursors           LibraryPrecursors (for mz/charge/sulfur)
  frag_lookup          FragmentLookupTable
  iso_splines          IsotopeSplineModel
  isotope_err_bounds   (UInt8, UInt8)
  min_fraction_transmitted   Float32
  irt_tol              Float32 fallback for per-precursor RT filter
  rt_binned_tol        per-file RT-binned tolerance struct (or nothing)
  max_iter_outer       Int   PMM outer iteration cap (e.g. 100)
  max_diff             Float32 PMM convergence threshold (e.g. 0.01)
"""
function extract_max_weight_in_window!(
    scratch::S,
    spectra,
    donor_pid::UInt32,
    precursors_passing::Set{UInt32},
    rt_index::retentionTimeIndex,
    center_rt::Float32,
    rt_tol::Float32,
    rt_irt_model,
    mass_error_model,
    quad_model,
    nce_model,
    precursors,
    frag_lookup,
    iso_splines,
    isotope_err_bounds::Tuple{UInt8, UInt8},
    min_fraction_transmitted::Float32,
    irt_tol::Float32,
    rt_binned_tol::Union{RTBinnedTolerance, Nothing},
    max_iter_outer::Int,
    max_diff::Float32,
    deconv_solver,
) where {S <: SearchDataStructures}
    # Pull scratch buffers from the per-thread SearchDataStructures
    Hs              = getHsFused(scratch)
    weights         = getTempWeights(scratch)
    colnorm2        = getColNorm2(scratch)
    precursor_weights = getPrecursorWeights(scratch)
    residuals       = getResiduals(scratch)
    fused_scratch   = getFusedScratch(scratch)
    corr_mz         = getScanCorrectedMz(scratch)
    obs_low         = getScanObsLow(scratch)
    obs_high        = getScanObsHigh(scratch)
    isotopes_buf    = getIsotopes(scratch)
    prec_trans_buf  = getPrecursorTransmission(scratch)
    id_to_col       = getIdToCol(scratch)
    unscored_psms   = getTuningUnscoredPsms(scratch)
    precs_temp      = getPrecIds(scratch)

    prec_mz_arr      = getMz(precursors)
    prec_charge_arr  = getCharge(precursors)
    prec_sulfur_arr  = getSulfurCount(precursors)
    prec_irt_arr     = getIrt(precursors)
    donor_mz         = prec_mz_arr[donor_pid]

    rt_lo = center_rt - rt_tol
    rt_hi = center_rt + rt_tol

    # Walk MS2 scans whose retention time is in window. Scans are not directly
    # indexed by RT here, so we sweep linearly within an expected range.
    # (Heuristic: get_ms2_scan_priority_order is RT-monotone for typical DIA;
    # for V0 we just loop all scans.)
    n_scans = length(spectra)

    best_scan_idx       = UInt32(0)
    best_weight         = 0.0f0
    best_log2_explained = 0.0f0
    best_log_by_ratio   = 0.0f0
    best_smoothed       = MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT
    best_rt             = 0.0f0

    kind = FusedRTIndexed(PartialCapture(), UInt8(255))

    for scan_idx in 1:n_scans
        getMsOrder(spectra, scan_idx) == 2 || continue
        rt = getRetentionTime(spectra, scan_idx)
        (rt < rt_lo || rt > rt_hi) && continue

        # Skip scans whose precursor isolation window doesn't cover donor_mz.
        quad_func = getQuadTransmissionFunction(
            quad_model,
            getCenterMz(spectra, scan_idx),
            getIsolationWidthMz(spectra, scan_idx),
        )
        min_qmz, max_qmz = getQuadrupoleBounds(quad_func)
        donor_mz < min_qmz && continue
        donor_mz > max_qmz && continue

        # 1. Enumerate already-IDed competitors at this scan from rt_index.
        rt_local_tol = rt_binned_tol !== nothing ? get_rt_tol(rt_binned_tol, Float32(rt)) : rt_tol
        rt_bin_start = max(searchsortedfirst(rt_index.rt_bins, rt - rt_local_tol,
                                             lt=(r,x)->r.lb<x) - 1, 1)
        rt_bin_stop  = min(searchsortedlast(rt_index.rt_bins, rt + rt_local_tol,
                                            lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))
        n_competitors = collect_rt_window_precursors!(
            precs_temp, rt_index, rt_bin_start, rt_bin_stop,
            precursors_passing, prec_mz_arr, prec_charge_arr, prec_sulfur_arr,
            iso_splines, quad_func, prec_trans_buf,
            isotope_err_bounds, min_fraction_transmitted,
            nothing, Float32(rt), rt_binned_tol, rt_local_tol,
        )

        # 2. Prepend donor pid (skip if it's already there — already-IDed pid case)
        donor_already_in = false
        @inbounds for j in 1:n_competitors
            if precs_temp[j] == donor_pid
                donor_already_in = true
                break
            end
        end
        if !donor_already_in
            if n_competitors + 1 > length(precs_temp)
                resize!(precs_temp, n_competitors + 1)
            end
            precs_temp[n_competitors + 1] = donor_pid
            n_competitors += 1
        end
        n_competitors == 0 && continue

        # 3. Pre-correct peak m/z.
        scan_mz  = getMzArray(spectra, scan_idx)
        scan_int = getIntensityArray(spectra, scan_idx)
        peak_mz_len = prepare_scan_peaks!(
            corr_mz, obs_low, obs_high,
            mass_error_model, scan_mz, scan_int,
            Float32(rt),
        )

        # 4. Fused fragment match + design matrix build.
        nmatches, nmisses = run_fused!(
            kind,
            Hs, unscored_psms, id_to_col, fused_scratch,
            corr_mz, obs_low, obs_high, peak_mz_len,
            isotopes_buf, prec_trans_buf,
            frag_lookup, nce_model,
            precs_temp, 1:n_competitors,
            prec_mz_arr, prec_charge_arr, prec_sulfur_arr, prec_irt_arr,
            iso_splines, quad_func, mass_error_model,
            scan_int, 0f0, Float32(Inf),
            (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx)),
            UInt8(1),                                    # n_frag_isotopes
            isotope_err_bounds,
        )

        if nmatches > 2
            n_active_cols = n_active(id_to_col)
            if n_active_cols > length(weights)
                new_entries = n_active_cols - length(weights) + 1000
                resize!(weights, length(weights) + new_entries)
                resize!(colnorm2, length(colnorm2) + new_entries)
            end
            initialize_weights!(id_to_col, weights, precursor_weights)
            solve_deconvolution!(
                deconv_solver,
                Hs, residuals, weights, colnorm2,
                getMu(scratch), getObserved(scratch),
                max_iter_outer, max_diff,
            )

            col = id_to_col[donor_pid]
            w = iszero(col) ? 0.0f0 : weights[col]
            if w > best_weight
                best_weight = w
                best_scan_idx = UInt32(scan_idx)
                best_rt = Float32(rt)
                # For the receiver-side payload features we approximate from
                # the per-scan deconv state. V0: log2_intensity_explained,
                # log_by_ratio, and smoothed_frag_sqrt are placeholder zeros;
                # they're computed properly downstream from frag intensities
                # at the chosen scan in V1.5. Track only the bare weight +
                # scan id for V0 to validate the loop end-to-end first.
                best_log2_explained = 0.0f0
                best_log_by_ratio   = 0.0f0
                best_smoothed       = MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT
            end
            update_precursor_weights!(id_to_col, weights, precursor_weights)
        end

        reset_scan_arrays!(id_to_col, Hs, unscored_psms)
    end

    best_weight > 0.0f0 || return nothing

    irt_obs = Float32(rt_irt_model(best_rt))
    return ReceiverExtraction(
        best_scan_idx,
        best_weight,
        best_log2_explained,
        irt_obs,
        best_rt,
        best_log_by_ratio,
        best_smoothed,
    )
end
