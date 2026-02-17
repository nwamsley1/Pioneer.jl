# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
    PrecursorBin{T<:AbstractFloat}

Container for precursor Items. 

### Fields

- precs::Vector{PrecursorBinItem{T}} -- Precursor items

### Examples

- PrecursorBin(T::DataType, N::Int) = PrecursorBin(Vector{PrecursorBinItem{T}}(undef, N)) -- Constructor

### GetterMethods

- getPrecursors(pb::PrecursorBin) = pb.precs

### Methods

- setPrecursor!(pb::PrecursorBin, index::Int, pbi::PrecursorBinItem)
"""
struct PrecursorBin{T<:AbstractFloat}
    precs::Vector{PrecursorBinFragment{T}}
end

getPrecursors(pb::PrecursorBin)  = pb.precs
getPrecursor(pb::PrecursorBin, i::Int64) = getPrecursors(pb)[i] 
getLength(pb::PrecursorBin)= length(pb.precs)

function setPrecursor!(pb::PrecursorBin, index::Int, pbi::PrecursorBinFragment)
    pb.precs[index] = pbi
end

PrecursorBin(T::DataType, N::Int) = PrecursorBin(Vector{PrecursorBinFragment{T}}(undef, N))

abstract type FragmentIndexBin{T<:AbstractFloat} end 

getLow(fb::FragmentIndexBin{T}) where {T<:AbstractFloat} = fb.lb
getHigh(fb::FragmentIndexBin{T}) where {T<:AbstractFloat} = fb.ub
getSubBinRange(fb::FragmentIndexBin{T}) where {T<:AbstractFloat} = fb.first_bin:fb.last_bin

struct FragIndexBin{T<:AbstractFloat} <: FragmentIndexBin{T}
    lb::T
    ub::T
    first_bin::UInt32
    last_bin::UInt32
end
ArrowTypes.arrowname(::Type{FragIndexBin{Float32}}) = :FragIndexBin
ArrowTypes.JuliaType(::Val{:FragIndexBin}) = FragIndexBin

struct IndexFragment{T<:AbstractFloat} <: LibraryFragmentIon{T}
    prec_id::UInt32
    prec_mz::T #Only need to tell if the peptide is in the quad isolation window
    score::UInt8 
    charge::UInt8
end
ArrowTypes.arrowname(::Type{IndexFragment{Float32}}) = :IndexFragment
ArrowTypes.JuliaType(::Val{:IndexFragment}) = IndexFragment

getScore(ind_frag::IndexFragment{T}) where {T<:AbstractFloat} = ind_frag.score

"""
    FragmentIndex{T<:AbstractFloat}

A fragment index for an MSFragger-style/fragment-centric search. Contains fragment bins 
that indicate the highest and lowest m/z of a fragment ion in the bin. Each `FragBin` links
to a `PrecursorBin`. A precursor bin has the precursor m/z's and peptide ID's for each fragment ion. 

### Fields

- fragment_bins::Vector{FragBin{T}} -- Lowest m/z of a fragment ion in the `FragBin`
- precursor_bins::Vector{PrecursorBin{T}} -- Highest m/z of a fragment ion in the `FragBin`

### Examples

- FragmentIndex(T::DataType, M::Int, N::Int) = FragmentIndex(fill(FragBin(), N), fill(PrecursorBin(T, M), N))

### GetterMethods

- getFragmentBin(fi::FragmentIndex, bin::Int) = fi.fragment_bins[bin]
- getPrecursorBin(fi::FragmentIndex, bin::Int64) = fi.precursor_bins[bin]

### Methods

- setFragmentBin!(fi::FragmentIndex, bin::Int64, frag_bin::FragBin)
- setPrecursorBinItem!(fi::FragmentIndex{T}, bin::Int64, index::Int64, prec_bin_item::PrecursorBinItem{T}) where {T<:AbstractFloat}
"""
struct FragmentIndex{T<:AbstractFloat}
    fragment_bins::Arrow.Struct{FragIndexBin, Tuple{Arrow.Primitive{T, Vector{T}}, Arrow.Primitive{T, Vector{T}}, Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{UInt32, Vector{UInt32}}}, (:lb, :ub, :first_bin, :last_bin)}
    rt_bins::Arrow.Struct{FragIndexBin, Tuple{Arrow.Primitive{T, Vector{T}}, Arrow.Primitive{T, Vector{T}}, Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{UInt32, Vector{UInt32}}}, (:lb, :ub, :first_bin, :last_bin)}
    fragments::Arrow.Struct{IndexFragment, Tuple{Arrow.Primitive{UInt32, Vector{UInt32}}, Arrow.Primitive{T, Vector{T}}, Arrow.Primitive{UInt8, Vector{UInt8}}, Arrow.Primitive{UInt8, Vector{UInt8}}}, (:prec_id, :prec_mz, :score, :charge)}
end
#struct FragmentIndex{T<:AbstractFloat}
#    fragment_bins::Vector{FragIndexBin{T}}
#    rt_bins::Vector{FragIndexBin{T}}
#    fragments::Vector{IndexFragment{T}}
#end

getFragBins(fi::FragmentIndex{T}) where {T<:AbstractFloat} = fi.fragment_bins
getRTBins(fi::FragmentIndex{T}) where {T<:AbstractFloat} = fi.rt_bins
getFragmentBin(fi::FragmentIndex{T}, frag_bin_idx::I) where {T<:AbstractFloat,I<:Integer} = getFragBins(fi)[frag_bin_idx]
getRTBin(fi::FragmentIndex{T}, rt_bin_idx::I) where {T<:AbstractFloat,I<:Integer} = getRTBins(fi)[rt_bin_idx]
getFragments(fi::FragmentIndex{T}) where {T<:AbstractFloat} = fi.fragments


abstract type SpectralLibrary end

struct FragmentIndexLibrary <: SpectralLibrary
    presearch_fragment_index::FragmentIndex{Float32}
    fragment_index::FragmentIndex{Float32}
    precursors::LibraryPrecursors
    proteins::LibraryProteins
    fragment_lookup_table::StandardFragmentLookup
end

struct SplineFragmentIndexLibrary <: SpectralLibrary
    presearch_fragment_index::FragmentIndex{Float32}
    fragment_index::FragmentIndex{Float32}
    precursors::LibraryPrecursors
    proteins::LibraryProteins
    fragment_lookup_table::SplineFragmentLookup
end


getPresearchFragmentIndex(sl::SpectralLibrary) = sl.presearch_fragment_index
getFragmentIndex(sl::SpectralLibrary) = sl.fragment_index
getPrecursors(sl::SpectralLibrary) = sl.precursors
getFragmentLookupTable(sl::SpectralLibrary) = sl.fragment_lookup_table
getProteins(sl::SpectralLibrary) = sl.proteins

#==========================================================
Batched Spectral Library Support
==========================================================#

"""
    BatchedLibraryConfig

Configuration for a batched spectral library.
Serialized to/from JSON for persistence.

# Indexing Strategy
Precursors are sorted by m/z and assigned sequential indices 1 to N.
Each batch contains exactly `batch_size` precursors (last batch may have fewer).
This allows O(1) computation of batch and local index from any precursor_idx:

    batch_idx = (prec_idx - 1) ÷ batch_size + 1
    local_idx = (prec_idx - 1) % batch_size + 1

No dictionary needed - just simple arithmetic.

# Fields
- `n_batches::Int` - Total number of batches
- `n_precursors::Int` - Total precursors across all batches
- `batch_size::Int` - Precursors per batch (last batch may have fewer)
- `batch_mz_ranges::Vector{Tuple{Float32, Float32}}` - (low, high) m/z for each batch
- `batch_min_charges::Vector{Int}` - Minimum precursor charge per batch (for scan skipping)
- `library_base_path::String` - Base path to the library directory
- `is_batched::Bool` - Whether this is a batched library
"""
struct BatchedLibraryConfig
    n_batches::Int
    n_precursors::Int
    batch_size::Int
    batch_mz_ranges::Vector{Tuple{Float32, Float32}}
    batch_min_charges::Vector{Int}
    library_base_path::String
    is_batched::Bool
end

# Constructor for backward-compatible single-batch mode
function BatchedLibraryConfig(library_path::String, n_precursors::Int)
    BatchedLibraryConfig(
        1,
        n_precursors,
        n_precursors,
        [(0.0f0, Inf32)],
        [1],  # Conservative default min charge
        library_path,
        false
    )
end

# Accessors for BatchedLibraryConfig
getBatchMzRange(config::BatchedLibraryConfig, batch_idx::Int) = config.batch_mz_ranges[batch_idx]
getBatchMinCharge(config::BatchedLibraryConfig, batch_idx::Int) = config.batch_min_charges[batch_idx]
getNPrecursors(config::BatchedLibraryConfig) = config.n_precursors
getBatchSize(config::BatchedLibraryConfig) = config.batch_size
getNBatches(config::BatchedLibraryConfig) = config.n_batches
getLibraryBasePath(config::BatchedLibraryConfig) = config.library_base_path
isBatched(config::BatchedLibraryConfig) = config.is_batched

"""
    get_batch_idx(config::BatchedLibraryConfig, prec_idx::UInt32) -> Int

Compute which batch contains the given precursor index. O(1) arithmetic.
"""
function get_batch_idx(config::BatchedLibraryConfig, prec_idx::UInt32)::Int
    return Int((prec_idx - 1) ÷ config.batch_size + 1)
end

"""
    get_local_idx(config::BatchedLibraryConfig, prec_idx::UInt32) -> Int

Compute the local index within the batch for a given precursor index. O(1) arithmetic.
"""
function get_local_idx(config::BatchedLibraryConfig, prec_idx::UInt32)::Int
    return Int((prec_idx - 1) % config.batch_size + 1)
end

"""
    get_prec_idx(config::BatchedLibraryConfig, batch_idx::Int, local_idx::Int) -> UInt32

Compute the global precursor index from batch and local indices. O(1) arithmetic.
"""
function get_prec_idx(config::BatchedLibraryConfig, batch_idx::Int, local_idx::Int)::UInt32
    return UInt32((batch_idx - 1) * config.batch_size + local_idx)
end

"""
    BatchedSpectralLibrary <: SpectralLibrary

A spectral library split into multiple batches for memory efficiency.
Only one batch is loaded at a time.

# Fields
- `config::BatchedLibraryConfig` - Batch configuration and file paths
- `current_batch_idx::Int` - Index of currently loaded batch (0 = none)
- `current_batch::Union{Nothing, SpectralLibrary}` - The loaded batch
- `proteins::LibraryProteins` - Shared protein data (always loaded)

# Indexing
No dictionary needed! Use config functions to compute batch/local indices:
- `get_batch_idx(config, prec_idx)` - which batch contains this precursor
- `get_local_idx(config, prec_idx)` - index within the batch
- `get_prec_idx(config, batch_idx, local_idx)` - global index from batch + local
"""
mutable struct BatchedSpectralLibrary <: SpectralLibrary
    config::BatchedLibraryConfig
    current_batch_idx::Int
    current_batch::Union{Nothing, SpectralLibrary}
    proteins::LibraryProteins
end

# Accessor methods for BatchedSpectralLibrary
getConfig(lib::BatchedSpectralLibrary) = lib.config
getCurrentBatchIdx(lib::BatchedSpectralLibrary) = lib.current_batch_idx
getCurrentBatch(lib::BatchedSpectralLibrary) = lib.current_batch
getProteins(lib::BatchedSpectralLibrary) = lib.proteins
getNBatches(lib::BatchedSpectralLibrary) = lib.config.n_batches
isBatchLoaded(lib::BatchedSpectralLibrary) = lib.current_batch !== nothing

# Indexing helpers that forward to config
get_batch_idx(lib::BatchedSpectralLibrary, prec_idx::UInt32) = get_batch_idx(lib.config, prec_idx)
get_local_idx(lib::BatchedSpectralLibrary, prec_idx::UInt32) = get_local_idx(lib.config, prec_idx)

# Forward methods to current batch when loaded
function getPresearchFragmentIndex(lib::BatchedSpectralLibrary)
    @assert isBatchLoaded(lib) "No batch currently loaded"
    return getPresearchFragmentIndex(lib.current_batch)
end

function getFragmentIndex(lib::BatchedSpectralLibrary)
    @assert isBatchLoaded(lib) "No batch currently loaded"
    return getFragmentIndex(lib.current_batch)
end

function getPrecursors(lib::BatchedSpectralLibrary)
    @assert isBatchLoaded(lib) "No batch currently loaded"
    return getPrecursors(lib.current_batch)
end

function getFragmentLookupTable(lib::BatchedSpectralLibrary)
    @assert isBatchLoaded(lib) "No batch currently loaded"
    return getFragmentLookupTable(lib.current_batch)
end

"""
    should_process_scan_for_batch(
        scan_isolation_center::Float32,
        scan_isolation_width::Float32,
        batch_mz_low::Float32,
        batch_mz_high::Float32,
        isotope_err_bounds::Tuple{UInt8, UInt8},
        min_charge::Int
    ) -> Bool

Check if a scan's isolation window could possibly contain precursors from this batch,
accounting for isotope error bounds. Returns false if the scan can be safely skipped.

# Arguments
- `scan_isolation_center`: Center m/z of the scan's isolation window
- `scan_isolation_width`: Width of the isolation window
- `batch_mz_low`: Lowest precursor m/z in this batch
- `batch_mz_high`: Highest precursor m/z in this batch
- `isotope_err_bounds`: (left, right) isotope error bounds, e.g., (1, 3)
- `min_charge`: Minimum precursor charge in batch (gives maximum m/z shift per isotope)

# Safety Margin
A 2 Da safety margin is added to account for mass calibration uncertainty and quad edge effects.
"""
function should_process_scan_for_batch(
    scan_isolation_center::Float32,
    scan_isolation_width::Float32,
    batch_mz_low::Float32,
    batch_mz_high::Float32,
    isotope_err_bounds::Tuple{UInt8, UInt8},
    min_charge::Int
)::Bool
    # Safety margin for calibration uncertainty and quad edge effects
    SCAN_SKIP_SAFETY_MARGIN = 2.0f0

    # Isotope spacing for the minimum charge (gives largest m/z shift)
    isotope_spacing = 1.00335f0 / min_charge

    # Extend batch m/z range to account for isotope errors
    # If isotope_err_bounds = (1, 3):
    #   - We look up to 3 isotopes to the RIGHT of precursor m/z
    #   - So a precursor could be UP TO 3*spacing BELOW the scan window and still match
    #   - We look up to 1 isotope to the LEFT of precursor m/z
    #   - So a precursor could be UP TO 1*spacing ABOVE the scan window and still match
    left_isotopes, right_isotopes = isotope_err_bounds
    extended_batch_low = batch_mz_low - right_isotopes * isotope_spacing - SCAN_SKIP_SAFETY_MARGIN
    extended_batch_high = batch_mz_high + left_isotopes * isotope_spacing + SCAN_SKIP_SAFETY_MARGIN

    # Calculate scan isolation window bounds
    scan_mz_low = scan_isolation_center - scan_isolation_width / 2
    scan_mz_high = scan_isolation_center + scan_isolation_width / 2

    # Check for overlap
    # No overlap if: scan ends before extended batch starts OR scan starts after extended batch ends
    return !(scan_mz_high < extended_batch_low || scan_mz_low > extended_batch_high)
end

#==========================================================
Reduced Spectral Library (for post-FirstPassSearch)
==========================================================#

"""
    ReducedSpectralLibrary <: SpectralLibrary

A reduced library containing only passing precursors from FirstPassSearch.
No fragment index needed - only used for stages after FirstPassSearch.

# Fields
- `precursors::LibraryPrecursors` - Precursor metadata (contiguous local indices)
- `proteins::LibraryProteins` - Shared protein data
- `fragment_lookup_table::StandardFragmentLookup` - Detailed fragments
- `original_prec_idx_to_local::Dictionary{UInt32, UInt32}` - Maps original precursor_idx to local index
- `local_to_original_prec_idx::Vector{UInt32}` - Maps local index to original precursor_idx

# Usage
After FirstPassSearch, the reduced library is created with only passing precursors.
Subsequent search stages (SecondPassSearch, ScoringSearch, etc.) use the reduced library,
which is much smaller and fits in memory.
"""
struct ReducedSpectralLibrary <: SpectralLibrary
    precursors::LibraryPrecursors
    proteins::LibraryProteins
    fragment_lookup_table::StandardFragmentLookup
    original_prec_idx_to_local::Dictionary{UInt32, UInt32}
    local_to_original_prec_idx::Vector{UInt32}
end

# Lookup by original precursor_idx
function get_local_idx(lib::ReducedSpectralLibrary, original_prec_idx::UInt32)::UInt32
    return lib.original_prec_idx_to_local[original_prec_idx]
end

function get_original_prec_idx(lib::ReducedSpectralLibrary, local_idx::Int)::UInt32
    return lib.local_to_original_prec_idx[local_idx]
end

# Check if a precursor is in the reduced library
function has_precursor(lib::ReducedSpectralLibrary, original_prec_idx::UInt32)::Bool
    return haskey(lib.original_prec_idx_to_local, original_prec_idx)
end

# Forward standard accessors
getPrecursors(lib::ReducedSpectralLibrary) = lib.precursors
getProteins(lib::ReducedSpectralLibrary) = lib.proteins
getFragmentLookupTable(lib::ReducedSpectralLibrary) = lib.fragment_lookup_table

# ReducedSpectralLibrary does not have fragment indices
# These will throw errors if called - they shouldn't be needed after FirstPassSearch
function getPresearchFragmentIndex(lib::ReducedSpectralLibrary)
    error("ReducedSpectralLibrary does not have a presearch fragment index. " *
          "This should only be used after FirstPassSearch.")
end

function getFragmentIndex(lib::ReducedSpectralLibrary)
    error("ReducedSpectralLibrary does not have a fragment index. " *
          "This should only be used after FirstPassSearch.")
end
