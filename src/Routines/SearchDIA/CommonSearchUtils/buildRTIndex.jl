# Construction of retentionTimeIndex from precursor RT/mz data.
#
# buildRtIndex groups precursors into RT bins of configurable width, with
# precursors sorted by m/z within each bin for efficient isolation window
# matching during chromatogram extraction.

"""
    buildRtIndex(rts, prec_mzs, prec_ids, bin_rt_size) -> retentionTimeIndex

Build an RT index from sorted RT values, precursor m/z values, and precursor IDs.

Precursors are grouped into bins where consecutive precursors span at most
`bin_rt_size` in RT. Within each bin, precursors are sorted by m/z.

# Arguments
- `rts`: Sorted retention time values (one per precursor)
- `prec_mzs`: Precursor m/z values
- `prec_ids`: Precursor IDs (UInt32)
- `bin_rt_size`: Maximum RT span per bin
"""
function buildRtIndex(rts::Vector{T}, prec_mzs::Vector{U}, prec_ids::Vector{I}, bin_rt_size::AbstractFloat) where {T,U<:AbstractFloat,I<:Integer}
    if isempty(rts)
        return retentionTimeIndex(T, U)
    end

    start_idx = 1
    start_rt = rts[start_idx]
    rt_index = retentionTimeIndex(T, U)
    i = 1
    while i < length(rts) + 1
        if ((rts[min(i + 1, length(rts))] - start_rt) > bin_rt_size) | (i == length(rts))
            push!(rt_index.rt_bins,
                rtIndexBin(rts[start_idx], rts[i],
                    [(zero(UInt32), zero(Float32)) for _ in 1:(i - start_idx + 1)]
                )
            )

            n = 1
            for idx in start_idx:(min(i, length(rts)))
                rt_index.rt_bins[end].prec[n] = (prec_ids[idx], prec_mzs[idx])
                n += 1
            end

            sort!(rt_index.rt_bins[end].prec, by = x -> last(x))
            i += 1
            start_idx = i
            start_rt = rts[min(start_idx, length(rts))]
            continue
        else
            i += 1
        end
    end

    # Final sort pass (redundant with in-loop sort but ensures correctness)
    for bin in rt_index.rt_bins
        sort!(bin.prec, by = x -> last(x))
    end

    return rt_index
end

"""DataFrame convenience: expects columns :irt, :prec_mz, :precursor_idx."""
buildRtIndex(PSMs::DataFrame; bin_rt_size::AbstractFloat = 0.1) =
    buildRtIndex(PSMs[:,:irt], PSMs[:,:prec_mz], PSMs[:,:precursor_idx], bin_rt_size)

buildRtIndex(PSMs::SubDataFrame; bin_rt_size::AbstractFloat = 0.1) =
    buildRtIndex(PSMs[:,:irt], PSMs[:,:prec_mz], PSMs[:,:precursor_idx], bin_rt_size)
