"""
Dispatch tag for the fragment correlation pass.
Used to dispatch a second getPSMS call that collects raw fragment match data
(selectTransitions! + matchPeaks!) without building a design matrix or scoring.
"""
struct FRAGCORR end

"""
    getPSMS(ms_file_idx, spectra, thread_task, precursors, ion_list, nce_model,
            scan_to_prec_idx, precursors_passed_scoring, search_data, params,
            qtm, mem, rt_to_irt_spline, irt_tol, ::FRAGCORR)

Second-pass getPSMS for fragment correlation features.
Runs selectTransitions! and matchPeaks! for each scan but skips buildDesignMatrix!
and scoring. Returns raw fragment match data grouped by precursor for downstream
correlation analysis.
"""
function getPSMS(
    ms_file_idx::UInt32,
    spectra::MassSpecData,
    thread_task::Vector{Int64},
    precursors::LibraryPrecursors,
    ion_list::LibraryFragmentLookup,
    nce_model::NceModel{Float32},
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed_scoring::Vector{UInt32},
    search_data::S,
    params::P,
    qtm::Q,
    mem::M,
    rt_to_irt_spline::Any,
    irt_tol::AbstractFloat,
    ::FRAGCORR
) where {M<:MassErrorModel, Q<:QuadTransmissionModel, S<:SearchDataStructures, P<:SearchParameters}

    isotopes = zeros(Float32, 5)
    precursor_transmission = zeros(Float32, 5)

    for scan_idx in thread_task
        (scan_idx == 0 || scan_idx > length(spectra)) && continue
        ismissing(scan_to_prec_idx[scan_idx]) && continue

        msn = getMsOrder(spectra, scan_idx)
        msn ∉ getSpecOrder(params) && continue

        # Ion Template Selection
        ion_idx, _ = selectTransitions!(
            getIonTemplates(search_data),
            StandardTransitionSelection(),
            getPrecEstimation(params),
            ion_list,
            nce_model,
            scan_to_prec_idx[scan_idx], precursors_passed_scoring,
            getMz(precursors),
            getCharge(precursors),
            getSulfurCount(precursors),
            getIrt(precursors),
            getIsoSplines(search_data),
            getQuadTransmissionFunction(qtm, getCenterMz(spectra, scan_idx), getIsolationWidthMz(spectra, scan_idx)),
            precursor_transmission, isotopes, getNFragIsotopes(params),
            getMaxFragRank(params),
            Float32(rt_to_irt_spline(getRetentionTime(spectra, scan_idx))),
            Float32(irt_tol),
            (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx));
            isotope_err_bounds = getIsotopeErrBounds(params)
        )

        ion_idx < 2 && continue

        # Match peaks
        nmatches, _ = matchPeaks!(
            getIonMatches(search_data),
            getIonMisses(search_data),
            getIonTemplates(search_data),
            ion_idx,
            getMzArray(spectra, scan_idx),
            getIntensityArray(spectra, scan_idx),
            mem,
            getHighMz(spectra, scan_idx),
            UInt32(scan_idx),
            ms_file_idx
        )

        nmatches < 2 && continue

        # Sort matches by (prec_id, peak_ind) for per-precursor grouping
        sort!(@view(getIonMatches(search_data)[1:nmatches]),
              by = x -> (x.prec_id, x.peak_ind), alg=QuickSort)

        # TODO: compute per-precursor fragment correlation features here
        # For now this is a placeholder — the match data is available in
        # getIonMatches(search_data)[1:nmatches] grouped by prec_id.
    end

    # TODO: return a DataFrame of per-precursor correlation features
    return DataFrame()
end
