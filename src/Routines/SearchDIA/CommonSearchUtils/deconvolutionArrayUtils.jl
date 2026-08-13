# Shared array management utilities for spectral deconvolution.
#
# Used by MainSearch, ParameterTuning, QuadTuning, and IntegrateChromatograms
# to manage working arrays (weights, precursor_weights, id_to_col) between
# scan iterations. All bulk ops on the precursor maps go through the
# AbstractPrecursorMap interface (active_keys / getindex / setindex!) so the
# same code works whether the maps are dense-backed or sparse-backed.

"""
    initialize_weights!(id_to_col, weights, precursor_weights)

Copy warm-start weights from `precursor_weights` (keyed by precursor id)
into `weights` (indexed by column number) for the current scan's active
precursors.
"""
function initialize_weights!(
    id_to_col::AbstractPrecursorMap{UInt16},
    weights::Vector{Float32},
    precursor_weights::AbstractPrecursorMap{Float32}
)
    @inbounds for (k, col) in active_keys(id_to_col)
        weights[col] = precursor_weights[k]
    end
end

"""
    update_precursor_weights!(id_to_col, weights, precursor_weights)

Copy solved weights back from `weights` (column-indexed) into
`precursor_weights` (precursor-id-keyed) so subsequent scans can warm-start.
"""
function update_precursor_weights!(
    id_to_col::AbstractPrecursorMap{UInt16},
    weights::Vector{Float32},
    precursor_weights::AbstractPrecursorMap{Float32}
)
    @inbounds for (k, col) in active_keys(id_to_col)
        precursor_weights[k] = weights[col]
    end
end

"""
    reset_scan_arrays!(id_to_col, Hs, unscored_psms)

Reset working arrays between scan iterations. Clears the sparse design
matrix, the id-to-column map, and zeroes out unscored PSM slots.
"""
function reset_scan_arrays!(
    id_to_col::AbstractPrecursorMap,
    Hs::AbstractSparseDesignMatrix,
    unscored_psms::AbstractVector
)
    @inbounds for i in 1:Hs.n
        unscored_psms[i] = eltype(unscored_psms)()
    end
    reset!(id_to_col)
    reset!(Hs)
end
