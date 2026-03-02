"""
Dispatch tag for the fragment correlation pass.
Used to dispatch a second getPSMS call that collects raw fragment match data
(selectTransitions! + matchPeaks!) without building a design matrix or scoring.
"""
struct FRAGCORR end

"""
    scoreFragCorr!(results, last_val, matches, nmatches, scan_idx; block_size=50000)

Iterate sorted matches (by `(prec_id, peak_ind)`), detect precursor boundaries,
and for each precursor group write a `FragCorrScore` whose intensity slots are
indexed by library rank (1–5). Returns the updated write counter.
"""
function scoreFragCorr!(
    results::Vector{FragCorrScore},
    last_val::Int64,
    matches::AbstractVector{<:MatchIon{Float32}},
    nmatches::Int64,
    scan_idx::UInt32;
    block_size::Int64 = 50000
)
    nmatches < 1 && return last_val

    # Intensity slots for the current precursor group
    s1 = zero(Float32); s2 = zero(Float32); s3 = zero(Float32)
    s4 = zero(Float32); s5 = zero(Float32)

    current_prec = getPrecID(matches[1])

    # Assign first match
    rank = getRank(matches[1])
    if     rank == 0x01; s1 = getIntensity(matches[1])
    elseif rank == 0x02; s2 = getIntensity(matches[1])
    elseif rank == 0x03; s3 = getIntensity(matches[1])
    elseif rank == 0x04; s4 = getIntensity(matches[1])
    elseif rank == 0x05; s5 = getIntensity(matches[1])
    end

    for i in 2:nmatches
        prec_id = getPrecID(matches[i])

        if prec_id != current_prec
            # Flush previous precursor group
            last_val += 1
            if last_val > length(results)
                resize!(results, length(results) + block_size)
            end
            @inbounds results[last_val] = FragCorrScore(current_prec, scan_idx, s1, s2, s3, s4, s5)

            # Reset for next group
            current_prec = prec_id
            s1 = zero(Float32); s2 = zero(Float32); s3 = zero(Float32)
            s4 = zero(Float32); s5 = zero(Float32)
        end

        rank = getRank(matches[i])
        if     rank == 0x01; s1 = getIntensity(matches[i])
        elseif rank == 0x02; s2 = getIntensity(matches[i])
        elseif rank == 0x03; s3 = getIntensity(matches[i])
        elseif rank == 0x04; s4 = getIntensity(matches[i])
        elseif rank == 0x05; s5 = getIntensity(matches[i])
        end
    end

    # Flush last group
    last_val += 1
    if last_val > length(results)
        resize!(results, length(results) + block_size)
    end
    @inbounds results[last_val] = FragCorrScore(current_prec, scan_idx, s1, s2, s3, s4, s5)

    return last_val
end

"""
    FragCorrTransitionSelection <: TransitionSelectionStrategy

Lightweight transition selection strategy for the FRAGCORR pass.
Only populates prec_id, mz, and rank in DetailedFrag — skips intensity prediction,
isotope pattern calculation, and precursor transmission modeling. Emits monoisotopic
fragments only.
"""
struct FragCorrTransitionSelection <: TransitionSelectionStrategy end

function _select_transitions_impl!(
    transitions::Vector{DetailedFrag{Float32}},
    ::FragCorrTransitionSelection,
    ::PrecEstimation,
    transition_idx::Int64,
    lookup::LibraryFragmentLookup,
    scan_to_prec_idx::UnitRange{Int64},
    precursors_passed_scoring::Vector{UInt32},
    prec_mzs::AbstractArray{Float32},
    prec_charges::AbstractArray{UInt8},
    prec_irts::AbstractArray{Float32},
    quad_transmission_func::QuadTransmissionFunction,
    max_frag_rank::UInt8,
    iRT::Float32,
    iRT_tol::Float32,
    frag_mz_bounds::Tuple{Float32, Float32};
    isotope_err_bounds::Tuple{I,I} = (3, 1),
    block_size::Int64 = 10000
) where {I<:Integer}

    fragments = getFragments(lookup)

    for i in scan_to_prec_idx
        prec_idx = precursors_passed_scoring[i]
        prec_mz = prec_mzs[prec_idx]
        prec_charge = prec_charges[prec_idx]

        # NOTE: iRT filter removed — RT eligibility is handled by the caller
        # via ExpandedPrecursorSet. Only the quad window m/z filter applies here.

        # Quad window filter (same as Standard)
        mz_low = getPrecMinBound(quad_transmission_func) - NEUTRON * first(isotope_err_bounds) / prec_charge
        mz_high = getPrecMaxBound(quad_transmission_func) + NEUTRON * last(isotope_err_bounds) / prec_charge
        (prec_mz < mz_low) | (prec_mz > mz_high) && continue

        # Direct fragment loop — monoisotopic only, no intensity prediction
        for frag_idx in getPrecFragRange(lookup, prec_idx)
            frag = fragments[frag_idx]
            getRank(frag) > max_frag_rank && continue

            frag_mz = Float32(getMZ(frag))
            (frag_mz < first(frag_mz_bounds) || frag_mz > last(frag_mz_bounds)) && continue

            transition_idx += 1
            transitions[transition_idx] = DetailedFrag(
                frag.prec_id,
                frag_mz,
                Float16(0),
                UInt16(0),
                false,
                false,
                false,
                false,
                UInt8(0),
                UInt8(0),
                UInt8(0),
                getRank(frag),
                UInt8(0)
            )
            ensureTransitionCapacity!(transitions, transition_idx, block_size)
        end
    end

    return transition_idx, 0
end

"""
    getPSMS(ms_file_idx, spectra, thread_task, precursors, ion_list,
            expanded, search_data, params, qtm, mem,
            rt_to_irt_spline, irt_tol, ::FRAGCORR)

Second-pass getPSMS for fragment correlation features.
Processes ALL MS2 scans (not just those with fragment index hits), using the
`ExpandedPrecursorSet` to determine which precursors are eligible at each scan's RT.
Runs selectTransitions! (lightweight, no isotope/intensity prediction) and matchPeaks!
but skips buildDesignMatrix! and scoring. Returns raw fragment match data
grouped by precursor for downstream correlation analysis.
"""
function getPSMS(
    ms_file_idx::UInt32,
    spectra::MassSpecData,
    thread_task::Vector{Int64},
    precursors::LibraryPrecursors,
    ion_list::LibraryFragmentLookup,
    expanded::ExpandedPrecursorSet,
    search_data::S,
    params::P,
    qtm::Q,
    mem::M,
    rt_to_irt_spline::Any,
    irt_tol::AbstractFloat,
    ::FRAGCORR
) where {M<:MassErrorModel, Q<:QuadTransmissionModel, S<:SearchDataStructures, P<:SearchParameters}

    n_scans_processed = 0
    n_scans_with_matches = 0
    total_matches = 0
    last_val = 0
    frag_corr_scores = getFragCorrScores(search_data)

    # Pre-allocate reusable buffers (no per-scan allocations)
    eligible_buf = Vector{UInt32}(undef, max(length(expanded.prec_ids), 1))
    trans_prec_buf = Vector{UInt32}(undef, 10000)
    scored_prec_buf = Vector{UInt32}(undef, 10000)

    for scan_idx in thread_task
        (scan_idx == 0 || scan_idx > length(spectra)) && continue

        msn = getMsOrder(spectra, scan_idx)
        msn ∉ getSpecOrder(params) && continue
        n_scans_processed += 1

        # Find eligible precursors for this scan's RT via binary search
        scan_rt = Float32(getRetentionTime(spectra, scan_idx))
        n_eligible = 0
        hi = searchsortedlast(expanded.rt_lows, scan_rt)
        for j in 1:hi
            if expanded.rt_highs[j] >= scan_rt
                n_eligible += 1
                @inbounds eligible_buf[n_eligible] = expanded.prec_ids[j]
            end
        end
        n_eligible < 1 && continue

        # Ion Template Selection — lightweight, no isotope/intensity prediction
        ion_idx, _ = selectTransitions!(
            getIonTemplates(search_data),
            FragCorrTransitionSelection(),
            getPrecEstimation(params),
            ion_list,
            1:n_eligible, eligible_buf,
            getMz(precursors),
            getCharge(precursors),
            getIrt(precursors),
            getQuadTransmissionFunction(qtm, getCenterMz(spectra, scan_idx), getIsolationWidthMz(spectra, scan_idx)),
            getMaxFragRank(params),
            Float32(rt_to_irt_spline(getRetentionTime(spectra, scan_idx))),
            Float32(irt_tol),
            (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx));
            isotope_err_bounds = getIsotopeErrBounds(params)
        )

        ion_idx < 1 && continue  # truly zero transitions → nothing to do

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

        if nmatches > 0
            n_scans_with_matches += 1
        end
        total_matches += nmatches

        # Sort and score matched precursors
        prev_last_val = last_val
        if nmatches > 0
            # Sort matches by (prec_id, peak_ind) for per-precursor grouping
            sort!(@view(getIonMatches(search_data)[1:nmatches]),
                  by = x -> (x.prec_id, x.peak_ind), alg=QuickSort)

            # Extract rank-indexed intensities per precursor
            last_val = scoreFragCorr!(
                frag_corr_scores,
                last_val,
                @view(getIonMatches(search_data)[1:nmatches]),
                nmatches,
                UInt32(scan_idx)
            )
        end

        # --- Collect unique prec IDs from transitions ---
        # (transitions are sorted by m/z, NOT grouped by prec_id,
        #  so we must collect, sort, and deduplicate)
        n_trans = 0
        templates = getIonTemplates(search_data)
        for j in 1:ion_idx
            n_trans += 1
            if n_trans > length(trans_prec_buf)
                resize!(trans_prec_buf, length(trans_prec_buf) * 2)
            end
            @inbounds trans_prec_buf[n_trans] = getPrecID(templates[j])
        end
        sort!(@view(trans_prec_buf[1:n_trans]))
        # Deduplicate in-place
        n_unique = 0
        for j in 1:n_trans
            if j == 1 || @inbounds trans_prec_buf[j] != trans_prec_buf[j-1]
                n_unique += 1
                @inbounds trans_prec_buf[n_unique] = trans_prec_buf[j]
            end
        end

        # --- Collect scored prec IDs (already sorted from scoreFragCorr!) ---
        n_scored = 0
        for j in (prev_last_val + 1):last_val
            n_scored += 1
            if n_scored > length(scored_prec_buf)
                resize!(scored_prec_buf, length(scored_prec_buf) * 2)
            end
            @inbounds scored_prec_buf[n_scored] = frag_corr_scores[j].precursor_idx
        end

        # --- Two-pointer merge: emit zeros for unscored prec IDs ---
        i_t = 1; i_s = 1
        while i_t <= n_unique
            if i_s > n_scored || @inbounds trans_prec_buf[i_t] < scored_prec_buf[i_s]
                # Not scored — emit zero row
                last_val += 1
                if last_val > length(frag_corr_scores)
                    resize!(frag_corr_scores, length(frag_corr_scores) + 50000)
                end
                @inbounds frag_corr_scores[last_val] = FragCorrScore(
                    trans_prec_buf[i_t], UInt32(scan_idx),
                    0f0, 0f0, 0f0, 0f0, 0f0
                )
                i_t += 1
            elseif @inbounds trans_prec_buf[i_t] > scored_prec_buf[i_s]
                i_s += 1
            else
                # Equal — already scored, skip
                i_t += 1
                i_s += 1
            end
        end
    end

    @info "FRAGCORR thread done: $n_scans_processed scans processed, " *
          "$n_scans_with_matches with matches, $total_matches total fragment matches, " *
          "$last_val precursor-scan records"

    df = DataFrame(@view(frag_corr_scores[1:last_val]))
    if nrow(df) > 0
        df[!, :rt] = Float32[getRetentionTime(spectra, sid) for sid in df.scan_idx]
        is_decoy = getIsDecoy(precursors)
        df[!, :target] = Bool[!is_decoy[pid] for pid in df.precursor_idx]
    else
        df[!, :rt] = Float32[]
        df[!, :target] = Bool[]
    end
    return df
end
