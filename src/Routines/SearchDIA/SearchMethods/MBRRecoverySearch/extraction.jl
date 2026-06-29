# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# Per-cell PMM extraction for MBRRecoverySearch.
#
# For one (donor_pid × receiver_file) cell, walks the receiver file's MS2
# scans in iRT space: each scan's RT is projected forward through the file's
# rt_to_irt spline and the scan is visited only when its iRT lands within
# `irt_tol` of the donor's `irt_obs`. This avoids needing an iRT→RT inverse
# (which Pioneer does not populate) and matches how the fragment-index search
# matches scans to precursors across files (_precompute_scan_properties).
#
# Within the window it runs the same per-scan pipeline
# IntegrateChromatogramsSearch uses (collect_rt_window_precursors! →
# run_fused! → solve_deconvolution!); the only difference is that the donor
# pid is manually prepended to the candidate list so it gets a column in the
# design matrix even though it's absent from `precursors_passing`.
#
# Returns the max-weight scan's features for use in inline MBR compute,
# or `nothing` if no scan in the window produced positive weight for the
# donor pid.

"""
    _recovery_smoothed_sqrt(u::MainUnscoredPSM) -> NTuple{8, Float32}

sqrt of the L1-normalized top-8 per-rank fragment trace intensities for one
scored row. Single-scan analog of the `frag{1..8}_smoothed_intensity` columns:
MainSearch's radius-1 chromatogram smoothing collapses to the raw per-rank
intensity when a scan has no flanking cycles, so the recovery seed uses the
raw `frag{r}_int` values directly. Mirrors `_read_smoothed_frag_sqrt`.
"""
@inline function _recovery_smoothed_sqrt(u::MainUnscoredPSM{Float32})
    f1 = u.frag1_int; f2 = u.frag2_int; f3 = u.frag3_int; f4 = u.frag4_int
    f5 = u.frag5_int; f6 = u.frag6_int; f7 = u.frag7_int; f8 = u.frag8_int
    s = f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8
    s > 0.0f0 || return MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT
    inv = 1.0f0 / s
    return (sqrt(f1 * inv), sqrt(f2 * inv), sqrt(f3 * inv), sqrt(f4 * inv),
            sqrt(f5 * inv), sqrt(f6 * inv), sqrt(f7 * inv), sqrt(f8 * inv))
end

"""
    extract_max_weight_in_window!(scratch, ...) -> Union{Nothing, ReceiverExtraction}

Inputs:
  scratch              SearchDataStructures with Hs, weights, etc. (thread-local)
  spectra              MassSpecData for the receiver file
  donor_pid            UInt32 pid to extract
  precursors_passing   Set{UInt32} of already-IDed pids in this file (used as
                       allowlist for COMPETITORS; donor is prepended manually)
  rt_index             retentionTimeIndex built from already-IDed pids
  target_irt           Float32 donor irt_obs — center of the iRT window
  rt_irt_model         per-file rt_to_irt spline (forward; used to project each
                       scan's RT into iRT space and for the output irt_obs)
  mass_error_model, quad_model, nce_model    per-file calibrated models
  precursors           LibraryPrecursors (for mz/charge/sulfur)
  frag_lookup          FragmentLookupTable
  iso_splines          IsotopeSplineModel
  isotope_err_bounds   (UInt8, UInt8)
  min_fraction_transmitted   Float32
  irt_tol              Float32 iRT window half-width (per-file getIrtErrors);
                       also the fallback for the per-scan competitor RT tol
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
    target_irt::Float32,
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
    # FusedStandard records full per-fragment scoring (b/y ion intensities,
    # per-rank trace intensities) into the MAIN unscored buffer via
    # apply_main_scoring! — needed for the secondary payload features. The
    # buffer auto-grows in record_match!'s ensure_unscored_capacity!.
    unscored_psms   = getMainUnscoredPsms(scratch)
    precs_temp      = getPrecIds(scratch)

    prec_mz_arr      = getMz(precursors)
    prec_charge_arr  = getCharge(precursors)
    prec_sulfur_arr  = getSulfurCount(precursors)
    prec_irt_arr     = getIrt(precursors)
    donor_mz         = prec_mz_arr[donor_pid]

    # Walk all MS2 scans; a scan is in-window when its forward-projected iRT
    # lands within irt_tol of the donor's target iRT. The forward rt_to_irt
    # model is always available (the iRT→RT inverse is not populated in
    # Pioneer), and iRT is the natural space for cross-file RT matching.
    n_scans = length(spectra)

    best_scan_idx       = UInt32(0)
    best_weight         = 0.0f0
    best_log2_explained = 0.0f0
    best_log_by_ratio   = 0.0f0
    best_smoothed       = MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT
    best_rt             = 0.0f0

    # FusedStandard (not FusedRTIndexed): same matched-fragment set and design
    # matrix (rank cap 255, candidates pre-filtered upstream so the inline
    # mz/iRT checks are no-ops with irt_tol=Inf below), but it ALSO populates
    # the per-precursor scoring buffer so we can read the donor's b/y ion and
    # per-rank fragment intensities at the chosen scan. Quant weights are
    # unchanged vs the FusedRTIndexed path.
    kind = FusedStandard(PartialPrecCapture(), UInt8(255))

    for scan_idx in 1:n_scans
        getMsOrder(spectra, scan_idx) == 2 || continue
        rt = getRetentionTime(spectra, scan_idx)
        scan_irt = Float32(rt_irt_model(rt))
        abs(scan_irt - target_irt) <= irt_tol || continue

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
        #    Competitor window stays in RT space (matches build_chromatograms);
        #    when no rt_binned_tol exists, convert irt_tol → RT via local slope.
        if rt_binned_tol !== nothing
            rt_local_tol = get_rt_tol(rt_binned_tol, Float32(rt))
        else
            h = 0.1f0
            local_slope = abs((Float32(rt_irt_model(rt + h)) - Float32(rt_irt_model(rt - h))) / (2f0 * h))
            rt_local_tol = irt_tol / max(local_slope, 0.01f0)
        end
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
            1,                                           # n_frag_isotopes (Int64)
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
                # Secondary payload features from the donor's scoring row at
                # this scan, mirroring MainSearch's Score! (ScoredPSMs.jl:199,
                # 242-243) exactly — col > 0 is guaranteed since w > 0.
                u = unscored_psms[col]
                spec_int_sum = Float32(sum(scan_int))
                best_log2_explained = log2(max(u.b_int + u.y_int, 1.0f-20) /
                                           max(spec_int_sum, 1.0f-20))
                best_log_by_ratio = log(u.b_int + 1.0f0) - log(u.y_int + 1.0f0)
                best_smoothed = _recovery_smoothed_sqrt(u)
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
