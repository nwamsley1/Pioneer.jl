"""
    searchFragmentBinUnconditional!(counter, fragments, frag_id_range)

Score all fragments in the range unconditionally (no binary search on prec_mz).
This is the key performance change — partitioning guarantees all fragments in the
partition have prec_mz within ~partition_width of the query window.
"""
@inline function searchFragmentBinUnconditional!(
        prec_id_to_score::Pioneer.Counter{UInt32, UInt8},
        fragments::Vector{Pioneer.IndexFragment{T}},
        frag_id_range::UnitRange{UInt32}) where {T<:AbstractFloat}
    @inbounds for i in frag_id_range
        frag = fragments[i]
        Pioneer.inc!(prec_id_to_score, Pioneer.getPrecID(frag), Pioneer.getScore(frag))
    end
    return nothing
end

"""
    queryFragmentPartitioned!(counter, frag_bin_max_idx, lower_bound_guess, upper_bound_guess,
                              frag_bins, fragments, frag_mz_min, frag_mz_max)

Same frag-bin lookup as `queryFragment!` (exponential search + binary search to find
matching frag bins), but calls `searchFragmentBinUnconditional!` instead of
`searchFragmentBin!`. No prec_mz_min/prec_mz_max needed.
"""
function queryFragmentPartitioned!(
        prec_id_to_score::Pioneer.Counter{UInt32, UInt8},
        frag_bin_max_idx::UInt32,
        lower_bound_guess::UInt32,
        upper_bound_guess::UInt32,
        frag_bins::Vector{Pioneer.FragIndexBin{T}},
        fragments::Vector{Pioneer.IndexFragment{T}},
        frag_mz_min::Float32,
        frag_mz_max::Float32) where {T<:AbstractFloat}

    lower_bound_guess, upper_bound_guess = Pioneer.exponentialFragmentBinSearch(
        frag_bins,
        frag_bin_max_idx,
        lower_bound_guess,
        upper_bound_guess,
        frag_mz_min,
        frag_mz_max,
        UInt32(2048)
    )

    frag_bin_idx = Pioneer.findFirstFragmentBin(
        frag_bins,
        lower_bound_guess,
        upper_bound_guess,
        frag_mz_min
    )

    @inbounds @fastmath begin
        if iszero(frag_bin_idx)
            return lower_bound_guess, upper_bound_guess
        end

        while frag_bin_idx <= frag_bin_max_idx
            frag_bin = frag_bins[frag_bin_idx]
            if Pioneer.getLow(frag_bin) > frag_mz_max
                break
            else
                if frag_bin_max_idx === frag_bin_idx
                    if Pioneer.getHigh(frag_bin) < frag_mz_min
                        break
                    end
                end
                frag_id_range = Pioneer.getSubBinRange(frag_bin)
                searchFragmentBinUnconditional!(prec_id_to_score, fragments, frag_id_range)
                frag_bin_idx += 1
            end
        end
    end

    return lower_bound_guess, upper_bound_guess
end

"""
    searchScanPartitioned!(counter, pfi, rt_bin_idx, irt_high, masses, mass_err_model, quad_func, isotope_err_bounds)

For each scan, determine which partitions overlap the quad isolation window, then
search each relevant partition using the unconditional scoring loop.
"""
function searchScanPartitioned!(
        prec_id_to_score::Pioneer.Counter{UInt32, UInt8},
        pfi::PartitionedFragmentIndex{T},
        rt_bin_idx::Int64,
        irt_high::Float32,
        masses::AbstractArray{Union{Missing, U}},
        mass_err_model::Pioneer.MassErrorModel,
        quad_transmission_func::Pioneer.QuadTransmissionFunction,
        isotope_err_bounds::Tuple{UInt8, UInt8}
        ) where {T<:AbstractFloat, U<:AbstractFloat}

    prec_min = U(Pioneer.getPrecMinBound(quad_transmission_func) - Pioneer.NEUTRON * first(isotope_err_bounds) / 2)
    prec_max = U(Pioneer.getPrecMaxBound(quad_transmission_func) + Pioneer.NEUTRON * last(isotope_err_bounds) / 2)

    # Determine which partitions overlap the precursor window
    first_k = get_partition_idx(pfi, T(prec_min))
    last_k  = get_partition_idx(pfi, T(prec_max))

    for k in first_k:last_k
        partition = getPartition(pfi, k)
        p_rt_bins   = Pioneer.getRTBins(partition)
        p_frag_bins = Pioneer.getFragBins(partition)
        p_fragments = Pioneer.getFragments(partition)

        isempty(p_frag_bins) && continue

        # Walk RT bins for this partition (same logic as searchScan!)
        local_rt_bin_idx = rt_bin_idx
        if local_rt_bin_idx > length(p_rt_bins)
            local_rt_bin_idx = length(p_rt_bins)
        end

        @inbounds @fastmath while Pioneer.getLow(p_rt_bins[local_rt_bin_idx]) < irt_high
            sub_bin_range = Pioneer.getSubBinRange(p_rt_bins[local_rt_bin_idx])
            min_frag_bin = first(sub_bin_range)
            max_frag_bin = last(sub_bin_range)

            # Skip empty RT bins (min > max means empty range)
            if min_frag_bin <= max_frag_bin
                lower_bound_guess = min_frag_bin
                upper_bound_guess = min_frag_bin

                for mass in masses
                    corrected_mz = Pioneer.getCorrectedMz(mass_err_model, mass)
                    frag_min, frag_max = Pioneer.getMzBoundsReverse(mass_err_model, corrected_mz)

                    lower_bound_guess, upper_bound_guess = queryFragmentPartitioned!(
                        prec_id_to_score,
                        max_frag_bin,
                        lower_bound_guess,
                        upper_bound_guess,
                        p_frag_bins,
                        p_fragments,
                        frag_min,
                        frag_max
                    )
                end
            end

            local_rt_bin_idx += 1
            if local_rt_bin_idx > length(p_rt_bins)
                local_rt_bin_idx = length(p_rt_bins)
                break
            end
        end
    end

    return nothing
end

"""
    searchFragmentIndexPartitioned(scan_to_prec_idx, pfi, spectra, thread_task,
                                   search_data, params, qtm, mem, rt_to_irt_spline,
                                   irt_tol, precursor_mzs)

Per-thread entry point for partitioned fragment index search. Same scan loop as the
current `searchFragmentIndex`, but calls `searchScanPartitioned!` and applies a
post-filter to remove false-positive precursors outside the actual quad window.
"""
function searchFragmentIndexPartitioned(
        scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
        pfi::PartitionedFragmentIndex{Float32},
        spectra::Pioneer.MassSpecData,
        thread_task::Vector{Int64},
        search_data::S,
        params::P,
        qtm::Q,
        mem::M,
        rt_to_irt_spline::Any,
        irt_tol::AbstractFloat,
        precursor_mzs::AbstractVector{Float32}
        ) where {M<:Pioneer.MassErrorModel, Q<:Pioneer.QuadTransmissionModel,
                 S<:Pioneer.SearchDataStructures, P<:Pioneer.FragmentIndexSearchParameters}

    prec_id = 0
    precursors_passed_scoring = Vector{UInt32}(undef, 250000)
    rt_bin_idx = 1

    # Use the first partition's RT bins for RT-bin advancement (all partitions share
    # the same RT bin boundaries)
    reference_rt_bins = Pioneer.getRTBins(getPartition(pfi, 1))

    for scan_idx in thread_task
        (scan_idx <= 0 || scan_idx > length(spectra)) && continue
        Pioneer.getMsOrder(spectra, scan_idx) ∉ Pioneer.getSpecOrder(params) && continue

        irt_lo, irt_hi = Pioneer.getRTWindow(rt_to_irt_spline(Pioneer.getRetentionTime(spectra, scan_idx)), irt_tol)

        # Advance RT bin index (same logic as searchFragmentIndex)
        while rt_bin_idx < length(reference_rt_bins) && Pioneer.getHigh(reference_rt_bins[rt_bin_idx]) < irt_lo
            rt_bin_idx += 1
        end
        while rt_bin_idx > 1 && Pioneer.getLow(reference_rt_bins[rt_bin_idx]) > irt_lo
            rt_bin_idx -= 1
        end

        counter = Pioneer.getPrecursorScores(search_data)

        # Partitioned search (no binary search on prec_mz)
        searchScanPartitioned!(
            counter,
            pfi,
            rt_bin_idx,
            irt_hi,
            Pioneer.getMzArray(spectra, scan_idx),
            mem,
            Pioneer.getQuadTransmissionFunction(qtm, Pioneer.getCenterMz(spectra, scan_idx), Pioneer.getIsolationWidthMz(spectra, scan_idx)),
            Pioneer.getIsotopeErrBounds(params)
        )

        # Post-filter: zero out scores for precursors outside the actual quad window
        quad_func = Pioneer.getQuadTransmissionFunction(qtm, Pioneer.getCenterMz(spectra, scan_idx), Pioneer.getIsolationWidthMz(spectra, scan_idx))
        iso_bounds = Pioneer.getIsotopeErrBounds(params)
        actual_prec_min = Float32(Pioneer.getPrecMinBound(quad_func) - Pioneer.NEUTRON * first(iso_bounds) / 2)
        actual_prec_max = Float32(Pioneer.getPrecMaxBound(quad_func) + Pioneer.NEUTRON * last(iso_bounds) / 2)

        @inbounds for idx in 1:(Pioneer.getSize(counter) - 1)
            pid = Pioneer.getID(counter, idx)
            pid == 0 && continue
            pmz = precursor_mzs[pid]
            if pmz < actual_prec_min || pmz > actual_prec_max
                counter.counts[pid] = zero(UInt8)
            end
        end

        # Filter and collect results (same as searchFragmentIndex)
        _match_count, _prec_count = Pioneer.filterPrecursorMatches!(counter, Pioneer.getMinIndexSearchScore(params))

        if Pioneer.getID(counter, 1) > 0
            start_idx = prec_id + 1
            n = 1
            while n <= counter.matches
                prec_id += 1
                if prec_id > length(precursors_passed_scoring)
                    append!(precursors_passed_scoring,
                            Vector{eltype(precursors_passed_scoring)}(undef, length(precursors_passed_scoring)))
                end
                precursors_passed_scoring[prec_id] = Pioneer.getID(counter, n)
                n += 1
            end
            scan_to_prec_idx[scan_idx] = start_idx:prec_id
        else
            scan_to_prec_idx[scan_idx] = missing
        end

        Pioneer.reset!(counter)
    end

    return precursors_passed_scoring[1:prec_id]
end
