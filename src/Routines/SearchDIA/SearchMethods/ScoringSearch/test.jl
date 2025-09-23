function getIrtBins(irts::AbstractVector{R}) where {R<:Real}
    sort_idx = sortperm(irts)
    bin_idx, bin_count = zero(UInt32), zero(UInt32)
    bin_idxs = similar(irts, UInt32, length(irts))
    for idx in sort_idx
        bin_count += one(UInt32)
        bin_idxs[idx] = bin_idx
        if bin_count >= IRT_BIN_SIZE
            bin_idx += one(UInt32)
            bin_count = zero(UInt32)
        end
    end
    return bin_idxs 
end


function getIrtBins!(psms::AbstractDataFrame)
    psms[!, :irt_bin_idx] = getIrtBins(psms.irts)
    return psms
end

function assign_random_target_decoy_pairs!(psms::DataFrame)
    last_pair_id = zero(UInt32)
    psms[!,:pair_id] = zeros(UInt32, nrow(psms))  # Initialize pair_id column
    psms[!,:irt_bin_idx] = getIrtBins(psms.irts)  # Ensure irt_bin_idx column exists
    for (irt_bin_idx, sub_psms) in pairs(groupby(psms, :irt_bin_idx))
        last_pair_id = assignPairIds!(sub_psms, last_pair_id)
    end
end

function assignPairIds!(psms::AbstractDataFrame, last_pair_id::UInt32)
    psms[!,:pair_id], last_pair_id = assign_pair_ids(
        psms.target, psms.decoy, psms.precursor_idx, psms.irt_bin_idx, last_pair_id
    )
    return last_pair_id
end

function assign_pair_ids(
    target::AbstractVector{Bool}, decoy::AbstractVector{Bool},
    precursor_idx::AbstractVector{UInt32}, irt_bin_idx::AbstractVector{UInt32},
    last_pair_id::UInt32
)
    targets = unique(precursor_idx[target])
    decoys = unique(precursor_idx[decoy])
    target_perm = randperm(length(targets); rng=MersenneTwister(PAIRING_RANDOM_SEED))
    precursor_idx_to_pair_id = Dict{UInt32,UInt32}()  # Map from precursor_idx to pair_id
    pair_ids = similar(precursor_idx, UInt32)
    for i in range(1, min(length(targets), length(decoys)))
        last_pair_id += one(UInt32)
        precursor_idx_to_pair_id[targets[target_perm[i]]] = last_pair_id
        precursor_idx_to_pair_id[decoys[i]] = last_pair_id
    end
    if length(decoys) < length(targets)
        @debug_l2 "Fewer decoy precursors ($(length(decoys))) than target precursors ($(length(targets))) in iRT bin $(first(irt_bin_idx)). Some targets will remain unpaired."
        for i in range(length(decoys)+1, length(targets))
            last_pair_id += one(UInt32)
            precursor_idx_to_pair_id[targets[target_perm[i]]] = last_pair_id
        end
    elseif length(targets) < length(decoys)
        @debug_l2 "Fewer target precursors ($(length(targets))) than decoy precursors ($(length(decoys))) in iRT bin $(first(irt_bin_idx)). Some decoys will remain unpaired."
        for i in range(length(targets)+1, length(decoys))
            last_pair_id += one(UInt32)
            precursor_idx_to_pair_id[decoys[i]] = last_pair_id
        end
    end
    for (row_idx, precursor_idx) in enumerate(precursor_idx)
        pair_ids[row_idx] = precursor_idx_to_pair_id[precursor_idx]
    end
    return pair_ids, last_pair_id
end