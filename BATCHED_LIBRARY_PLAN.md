# Batched Spectral Library Implementation Plan

## Overview

Split spectral libraries into N batches to handle very large libraries that exceed memory limits. During early search passes (ParameterTuning, FirstPass), loop through each batch. After FirstPassSearch, create a reduced library containing only passing precursors for subsequent searches.

## Partitioning Strategy

**Partition by precursor m/z ranges** (not RT):
- Each precursor belongs to exactly one batch (no duplicates)
- Fragment index queries filter by precursor m/z for quad isolation
- m/z distribution is relatively uniform across proteomes
- Preserves original precursor_idx values globally

### Key Optimization: Scan Skipping

Since batches are partitioned by precursor m/z, we can **skip entire MS/MS scans** where the isolation window doesn't overlap with the batch's m/z range. However, we must account for **isotope error bounds** - precursors can match scans via their isotope peaks, not just the monoisotopic peak.

**Isotope Error Consideration:**
- Typical isotope_err_bounds: `(-1, 3)` meaning look 1 isotope left, 3 isotopes right
- Isotope spacing ≈ 1.003 Da / charge
- For charge 2: spacing ≈ 0.5 Da
- A precursor at m/z 549 could match a scan at 550-560 via its +1 isotope peak

```julia
"""
    should_process_scan_for_batch(
        scan_isolation_center::Float32,
        scan_isolation_width::Float32,
        batch_mz_low::Float32,
        batch_mz_high::Float32,
        isotope_err_bounds::Tuple{Int, Int},
        min_charge::Int
    ) -> Bool

Check if a scan's isolation window could possibly contain precursors from this batch,
accounting for isotope error bounds. Returns false if the scan can be safely skipped.

# Arguments
- `scan_isolation_center`: Center m/z of the scan's isolation window
- `scan_isolation_width`: Width of the isolation window
- `batch_mz_low`: Lowest precursor m/z in this batch
- `batch_mz_high`: Highest precursor m/z in this batch
- `isotope_err_bounds`: (left, right) isotope error bounds, e.g., (-1, 3)
- `min_charge`: Minimum precursor charge in batch (gives maximum m/z shift per isotope)
"""
function should_process_scan_for_batch(
    scan_isolation_center::Float32,
    scan_isolation_width::Float32,
    batch_mz_low::Float32,
    batch_mz_high::Float32,
    isotope_err_bounds::Tuple{Int, Int},
    min_charge::Int
)::Bool
    # Isotope spacing for the minimum charge (gives largest m/z shift)
    isotope_spacing = 1.00335f0 / min_charge

    # Extend batch m/z range to account for isotope errors
    #
    # If isotope_err_bounds = (-1, 3):
    #   - We look up to 3 isotopes to the RIGHT of precursor m/z
    #   - So a precursor could be UP TO 3*spacing BELOW the scan window and still match
    #   - We look up to 1 isotope to the LEFT of precursor m/z
    #   - So a precursor could be UP TO 1*spacing ABOVE the scan window and still match
    #
    # Therefore:
    #   extended_low = batch_low - (right_isotopes * spacing)
    #   extended_high = batch_high + (left_isotopes * spacing)  [left is negative, so subtract]

    left_isotopes, right_isotopes = isotope_err_bounds
    extended_batch_low = batch_mz_low - right_isotopes * isotope_spacing
    extended_batch_high = batch_mz_high - left_isotopes * isotope_spacing  # left is negative

    # Calculate scan isolation window bounds
    scan_mz_low = scan_isolation_center - scan_isolation_width / 2
    scan_mz_high = scan_isolation_center + scan_isolation_width / 2

    # Check for overlap
    # No overlap if: scan ends before extended batch starts OR scan starts after extended batch ends
    return !(scan_mz_high < extended_batch_low || scan_mz_low > extended_batch_high)
end
```

**Example:**
- Batch covers precursors m/z 400-550
- isotope_err_bounds = (-1, 3), min_charge = 2
- isotope_spacing = 1.003/2 ≈ 0.5 Da
- Extended range: (400 - 3*0.5) to (550 + 1*0.5) = 398.5 to 550.5
- Scan with isolation window 395-405: **process** (overlaps extended range)
- Scan with isolation window 560-570: **skip** (no overlap)

**Additional Safety Margins:**
Consider adding a small buffer (e.g., 1-2 Da) to account for:
- Mass calibration uncertainty
- Quad transmission edge effects
- Rounding errors

```julia
const SCAN_SKIP_SAFETY_MARGIN = 2.0f0  # Da

extended_batch_low = batch_mz_low - right_isotopes * isotope_spacing - SCAN_SKIP_SAFETY_MARGIN
extended_batch_high = batch_mz_high - left_isotopes * isotope_spacing + SCAN_SKIP_SAFETY_MARGIN
```

**Implementation Location:**
This check should be added in the scan iteration loop in `LibrarySearch.jl`, before any expensive fragment matching operations. The batch m/z range and min_charge should be stored in `BatchedLibraryConfig` for quick access.

## File Structure

```
spectral_library/
├── library_config.json              # BatchedLibraryConfig
├── proteins_table.arrow             # Shared across batches
├── batch_1/
│   ├── f_index_fragments.arrow
│   ├── f_index_rt_bins.arrow
│   ├── f_index_fragment_bins.arrow
│   ├── presearch_f_index_*.arrow
│   ├── detailed_fragments.jld2
│   ├── precursor_to_fragment_indices.jld2
│   └── precursors_table.arrow       # Subset for this batch
├── batch_2/
│   └── ...
└── reduced_library/                 # Created after FirstPassSearch
    └── ... (standard single-batch structure)
```

---

## Phase 1: New Data Structures

### File: `src/structs/LibraryFragmentIndex.jl`

Add after existing SpectralLibrary types (around line 147):

```julia
#==========================================================
Batched Spectral Library Support
==========================================================#

"""
Configuration for a batched spectral library.
Serialized to/from JSON for persistence.

# Indexing Strategy
Precursors are sorted by m/z and assigned sequential indices 1 to N.
Each batch contains exactly `batch_size` precursors (last batch may have fewer).
This allows O(1) computation of batch and local index from any precursor_idx:

    batch_idx = (prec_idx - 1) ÷ batch_size + 1
    local_idx = (prec_idx - 1) % batch_size + 1

No dictionary needed - just simple arithmetic.
"""
struct BatchedLibraryConfig
    n_batches::Int
    n_precursors::Int                                  # Total precursors across all batches
    batch_size::Int                                    # Precursors per batch (last batch may have fewer)
    batch_mz_ranges::Vector{Tuple{Float32, Float32}}   # (low, high) m/z for each batch
    batch_min_charges::Vector{Int}                     # Minimum precursor charge per batch (for scan skipping)
    library_base_path::String
    is_batched::Bool
end

# Constructor for backward-compatible single-batch mode
BatchedLibraryConfig(library_path::String, n_precursors::Int) = BatchedLibraryConfig(
    1,
    n_precursors,
    n_precursors,
    [(0.0f0, Inf32)],
    [1],  # Conservative default min charge
    library_path,
    false
)

# Accessors
getBatchMzRange(config::BatchedLibraryConfig, batch_idx::Int) = config.batch_mz_ranges[batch_idx]
getBatchMinCharge(config::BatchedLibraryConfig, batch_idx::Int) = config.batch_min_charges[batch_idx]
getNPrecursors(config::BatchedLibraryConfig) = config.n_precursors
getBatchSize(config::BatchedLibraryConfig) = config.batch_size

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
```

---

## Phase 2: Library Loading

### File: `src/Routines/SearchDIA/ParseInputs/loadSpectralLibrary.jl`

Add these functions before the existing `loadSpectralLibrary`:

```julia
using JSON3

"""
    load_library_config(config_path::String) -> BatchedLibraryConfig

Load batched library configuration from JSON file.
"""
function load_library_config(config_path::String)::BatchedLibraryConfig
    config_data = JSON3.read(read(config_path, String))

    batch_ranges = [(Float32(r[1]), Float32(r[2])) for r in config_data.batch_mz_ranges]

    return BatchedLibraryConfig(
        Int(config_data.n_batches),
        batch_ranges,
        String(config_data.library_base_path),
        Bool(config_data.is_batched)
    )
end

"""
    save_library_config(config_path::String, config::BatchedLibraryConfig)

Save batched library configuration to JSON file.
"""
function save_library_config(config_path::String, config::BatchedLibraryConfig)
    config_data = Dict(
        "n_batches" => config.n_batches,
        "batch_mz_ranges" => config.batch_mz_ranges,
        "library_base_path" => config.library_base_path,
        "is_batched" => config.is_batched
    )
    open(config_path, "w") do f
        JSON3.write(f, config_data)
    end
end

"""
    get_batch_path(config::BatchedLibraryConfig, batch_idx::Int) -> String

Get the directory path for a specific batch.
"""
function get_batch_path(config::BatchedLibraryConfig, batch_idx::Int)::String
    return joinpath(config.library_base_path, "batch_$batch_idx")
end

"""
    loadBatchedSpectralLibrary(lib_dir::String, config::BatchedLibraryConfig, params::PioneerParameters)

Load a batched spectral library. Only loads shared data initially; batches are loaded on demand.

No dictionary needed for precursor->batch mapping! The config contains batch_size,
so we can compute batch_idx = (prec_idx - 1) ÷ batch_size + 1 in O(1).
"""
function loadBatchedSpectralLibrary(
    lib_dir::String,
    config::BatchedLibraryConfig,
    params::PioneerParameters
)::BatchedSpectralLibrary

    # Load shared protein data
    proteins = Arrow.Table(joinpath(lib_dir, "proteins_table.arrow"))

    return BatchedSpectralLibrary(
        config,
        0,  # No batch loaded initially
        nothing,
        SetProteins(proteins)
    )
end

"""
    load_batch!(lib::BatchedSpectralLibrary, batch_idx::Int, params::PioneerParameters)

Load a specific batch into memory. Unloads any previously loaded batch first.
"""
function load_batch!(lib::BatchedSpectralLibrary, batch_idx::Int, params::PioneerParameters)
    # Check if already loaded
    if lib.current_batch_idx == batch_idx && lib.current_batch !== nothing
        return nothing
    end

    # Unload current batch
    unload_current_batch!(lib)

    # Load new batch using existing loadSpectralLibrary for the batch directory
    batch_path = get_batch_path(lib.config, batch_idx)
    lib.current_batch = loadSpectralLibrary_internal(batch_path, params)
    lib.current_batch_idx = batch_idx

    @info "Loaded batch $batch_idx of $(lib.config.n_batches)"

    return nothing
end

"""
    unload_current_batch!(lib::BatchedSpectralLibrary)

Unload the current batch to free memory.
"""
function unload_current_batch!(lib::BatchedSpectralLibrary)
    if lib.current_batch !== nothing
        lib.current_batch = nothing
        lib.current_batch_idx = 0
        GC.gc()  # Trigger garbage collection to free memory
    end
    return nothing
end

# Rename existing function to internal version
function loadSpectralLibrary_internal(SPEC_LIB_DIR::String, params::PioneerParameters)
    # ... existing loadSpectralLibrary code ...
end
```

### Modify existing `loadSpectralLibrary`:

```julia
function loadSpectralLibrary(SPEC_LIB_DIR::String, params::PioneerParameters)
    # Check for batched library
    config_path = joinpath(SPEC_LIB_DIR, "library_config.json")

    if isfile(config_path)
        config = load_library_config(config_path)
        if config.is_batched
            @info "Detected batched spectral library with $(config.n_batches) batches"
            return loadBatchedSpectralLibrary(SPEC_LIB_DIR, config, params)
        end
    end

    # Original monolithic loading
    return loadSpectralLibrary_internal(SPEC_LIB_DIR, params)
end
```

---

## Phase 3: SearchContext Extensions

### File: `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl`

Add to SearchContext struct (around line 220):

```julia
mutable struct SearchContext{N,L<:SpectralLibrary,M<:MassSpecDataReference}
    # ... existing fields ...

    # Batched library support
    reduced_library::Base.Ref{Union{Nothing, SpectralLibrary}}

    # Constructor - add initialization
    function SearchContext(
        spec_lib::L,
        temp_structures::AbstractVector{<:SearchDataStructures},
        mass_spec_data_reference::M,
        n_threads::Int64,
        n_precursors::Int64,
        buffer_size::Int64
    ) where {L<:SpectralLibrary,M<:MassSpecDataReference}
        N = length(temp_structures)
        new{N,L,M}(
            # ... existing initializations ...
            Ref{Union{Nothing, SpectralLibrary}}(nothing),  # reduced_library
        )
    end
end
```

Add accessor methods (around line 420):

```julia
#==========================================================
Batched Library Support
==========================================================#

"""
    isBatchedLibrary(ctx::SearchContext) -> Bool

Check if the search context is using a batched spectral library.
"""
isBatchedLibrary(ctx::SearchContext) = isa(getSpecLib(ctx), BatchedSpectralLibrary)

"""
    getReducedLibrary(ctx::SearchContext) -> Union{Nothing, SpectralLibrary}

Get the reduced library created after FirstPassSearch, or nothing if not yet created.
"""
getReducedLibrary(ctx::SearchContext) = ctx.reduced_library[]

"""
    setReducedLibrary!(ctx::SearchContext, lib::SpectralLibrary)

Set the reduced library after FirstPassSearch.
"""
setReducedLibrary!(ctx::SearchContext, lib::SpectralLibrary) = (ctx.reduced_library[] = lib)

"""
    getActiveLibrary(ctx::SearchContext) -> SpectralLibrary

Get the library to use for searching. Returns reduced library if available,
otherwise returns the main spectral library.
"""
function getActiveLibrary(ctx::SearchContext)
    reduced = getReducedLibrary(ctx)
    return reduced !== nothing ? reduced : getSpecLib(ctx)
end
```

---

## Phase 4: Batched Search in ParameterTuningSearch

### File: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/ParameterTuningSearch.jl`

Add helper function:

```julia
"""
    search_all_batches(
        spectra::MassSpecData,
        search_context::SearchContext,
        params::ParameterTuningSearchParameters,
        ms_file_idx::Int64,
        search_func::Function
    ) -> DataFrame

Execute search across all batches of a batched library, merging results.
"""
function search_all_batches(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::ParameterTuningSearchParameters,
    ms_file_idx::Int64,
    search_func::Function
)
    lib = getSpecLib(search_context)

    if !isa(lib, BatchedSpectralLibrary)
        # Not batched, use normal search
        return search_func(spectra, search_context, params, ms_file_idx, lib)
    end

    all_psms = DataFrame[]
    pioneer_params = getParams(search_context)  # Get full params for batch loading

    for batch_idx in 1:getNBatches(lib)
        # Load batch
        load_batch!(lib, batch_idx, pioneer_params)

        # Search this batch
        batch_psms = search_func(spectra, search_context, params, ms_file_idx, lib.current_batch)

        if !isempty(batch_psms)
            push!(all_psms, batch_psms)
        end
    end

    # Unload final batch
    unload_current_batch!(lib)

    # Merge results
    if isempty(all_psms)
        return DataFrame()
    end

    return vcat(all_psms...)
end
```

Modify `collect_psms` or `library_search` to use batched search:

```julia
function library_search(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::ParameterTuningSearchParameters,
    ms_file_idx::Int64
)
    lib = getSpecLib(search_context)

    if isa(lib, BatchedSpectralLibrary)
        return search_all_batches(
            spectra, search_context, params, ms_file_idx,
            _library_search_single_batch
        )
    else
        return _library_search_single_batch(spectra, search_context, params, ms_file_idx, lib)
    end
end

# Rename existing implementation
function _library_search_single_batch(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::ParameterTuningSearchParameters,
    ms_file_idx::Int64,
    spec_lib::SpectralLibrary
)
    # ... existing library_search code, using spec_lib parameter ...
end
```

---

## Phase 5: Batched Search in FirstPassSearch

### File: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl`

Apply same pattern as ParameterTuningSearch:

```julia
function library_search(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::FirstPassSearchParameters,
    ms_file_idx::Int64
)
    lib = getSpecLib(search_context)

    if isa(lib, BatchedSpectralLibrary)
        return search_all_batches_first_pass(
            spectra, search_context, params, ms_file_idx
        )
    else
        return _library_search_single_batch(spectra, search_context, params, ms_file_idx, lib)
    end
end

function search_all_batches_first_pass(
    spectra::MassSpecData,
    search_context::SearchContext,
    params::FirstPassSearchParameters,
    ms_file_idx::Int64
)
    lib = getSpecLib(search_context)
    all_psms = DataFrame[]
    pioneer_params = getParams(search_context)

    for batch_idx in 1:getNBatches(lib)
        @info "FirstPassSearch: Processing batch $batch_idx of $(getNBatches(lib))"

        # Load batch
        load_batch!(lib, batch_idx, pioneer_params)

        # Search this batch
        batch_psms = _library_search_single_batch(
            spectra, search_context, params, ms_file_idx, lib.current_batch
        )

        if !isempty(batch_psms)
            push!(all_psms, batch_psms)
        end
    end

    # Unload final batch
    unload_current_batch!(lib)

    if isempty(all_psms)
        return DataFrame()
    end

    return vcat(all_psms...)
end
```

---

## Phase 6: Create Reduced Library

**Key Insight:** The fragment index is only used during ParameterTuning and FirstPassSearch. After FirstPassSearch, subsequent stages (SecondPassSearch, ScoringSearch, etc.) use RT indices and direct precursor lookups. Therefore, the reduced library **does not need a fragment index** - just precursor metadata and detailed fragments.

### Precursor Index Handling

The reduced library will have non-contiguous original precursor_idx values (e.g., [5, 100, 1050, 1900, ...]).
Since the reduced library is small (~5-10k precursors), we use a dictionary for lookup:

```julia
"""
    ReducedSpectralLibrary <: SpectralLibrary

A reduced library containing only passing precursors from FirstPassSearch.
No fragment index needed - only used for stages after FirstPassSearch.

# Fields
- `precursors::LibraryPrecursors` - Precursor metadata (contiguous local indices)
- `proteins::LibraryProteins` - Shared protein data
- `fragment_lookup_table::StandardFragmentLookup` - Detailed fragments
- `original_prec_idx_to_local::Dictionary{UInt32, Int}` - Maps original precursor_idx to local index
- `local_to_original_prec_idx::Vector{UInt32}` - Maps local index to original precursor_idx
"""
struct ReducedSpectralLibrary <: SpectralLibrary
    precursors::LibraryPrecursors
    proteins::LibraryProteins
    fragment_lookup_table::StandardFragmentLookup
    original_prec_idx_to_local::Dictionary{UInt32, Int}
    local_to_original_prec_idx::Vector{UInt32}
end

# Lookup by original precursor_idx
function get_local_idx(lib::ReducedSpectralLibrary, original_prec_idx::UInt32)::Int
    return lib.original_prec_idx_to_local[original_prec_idx]
end

function get_original_prec_idx(lib::ReducedSpectralLibrary, local_idx::Int)::UInt32
    return lib.local_to_original_prec_idx[local_idx]
end

# Forward standard accessors
getPrecursors(lib::ReducedSpectralLibrary) = lib.precursors
getProteins(lib::ReducedSpectralLibrary) = lib.proteins
getFragmentLookupTable(lib::ReducedSpectralLibrary) = lib.fragment_lookup_table
```

### New File: `src/Routines/SearchDIA/CommonSearchUtils/createReducedLibrary.jl`

```julia
"""
Functions for creating a reduced spectral library from passing precursors.
Used after FirstPassSearch to create a memory-efficient library for subsequent passes.

NOTE: The reduced library does NOT have a fragment index since it's only used
after FirstPassSearch, which is the last stage requiring the fragment index.
"""

export create_reduced_library

"""
    create_reduced_library(
        search_context::SearchContext,
        passing_precursor_ids::Vector{UInt32},
        params::PioneerParameters
    ) -> ReducedSpectralLibrary

Create a reduced spectral library containing only the specified precursors.
Extracts precursor data and detailed fragments from batched library.
No fragment index is built since subsequent stages don't need it.
"""
function create_reduced_library(
    search_context::SearchContext,
    passing_precursor_ids::Vector{UInt32},
    params::PioneerParameters
)::ReducedSpectralLibrary

    original_lib = getSpecLib(search_context)

    @info "Creating reduced library with $(length(passing_precursor_ids)) precursors (no fragment index needed)"

    if isa(original_lib, BatchedSpectralLibrary)
        return _create_reduced_from_batched(original_lib, passing_precursor_ids, params)
    else
        return _create_reduced_from_monolithic(original_lib, passing_precursor_ids, params)
    end
end

"""
    _create_reduced_from_batched(
        batched_lib::BatchedSpectralLibrary,
        passing_ids::Vector{UInt32},
        params::PioneerParameters
    ) -> ReducedSpectralLibrary

Extract passing precursors from batched library and build reduced library in memory.
"""
function _create_reduced_from_batched(
    batched_lib::BatchedSpectralLibrary,
    passing_ids::Vector{UInt32},
    params::PioneerParameters
)::ReducedSpectralLibrary

    config = getConfig(batched_lib)

    # Group passing precursors by batch using O(1) arithmetic
    batch_to_precs = Dict{Int, Vector{UInt32}}()
    for prec_id in passing_ids
        batch_idx = get_batch_idx(config, prec_id)
        if !haskey(batch_to_precs, batch_idx)
            batch_to_precs[batch_idx] = UInt32[]
        end
        push!(batch_to_precs[batch_idx], prec_id)
    end

    # Collect data from each batch
    all_precursor_rows = []
    all_detailed_frags = Vector{DetailedFrag{Float32}}()
    prec_frag_ranges = Vector{UnitRange{UInt64}}()

    # Build bidirectional index mapping
    original_prec_idx_to_local = Dictionary{UInt32, Int}()
    local_to_original_prec_idx = Vector{UInt32}()

    current_frag_offset = UInt64(1)
    current_local_idx = 1

    for batch_idx in sort(collect(keys(batch_to_precs)))
        batch_prec_ids = batch_to_precs[batch_idx]

        @info "Extracting $(length(batch_prec_ids)) precursors from batch $batch_idx"

        # Load batch
        load_batch!(batched_lib, batch_idx, params)
        batch = batched_lib.current_batch

        # Get batch data
        batch_precursors = getPrecursors(batch)
        batch_frag_lookup = getFragmentLookupTable(batch)
        batch_frags = batch_frag_lookup.frags
        batch_ranges = batch_frag_lookup.prec_frag_ranges

        # Extract data for each passing precursor in this batch
        for original_prec_id in batch_prec_ids
            # Get local index within batch using O(1) arithmetic
            batch_local_idx = get_local_idx(config, original_prec_id)

            # Build index mapping
            insert!(original_prec_idx_to_local, original_prec_id, current_local_idx)
            push!(local_to_original_prec_idx, original_prec_id)

            # Collect precursor row (copy all columns)
            # The precursor row keeps original_prec_id for reference
            push!(all_precursor_rows, batch_precursors.data[batch_local_idx])

            # Get fragment range for this precursor
            frag_range = batch_ranges[batch_local_idx]
            n_frags = length(frag_range)

            # Copy fragments
            for i in frag_range
                push!(all_detailed_frags, batch_frags[i])
            end

            # Record new range (using new local indexing)
            new_range = current_frag_offset:(current_frag_offset + n_frags - 1)
            push!(prec_frag_ranges, new_range)

            current_frag_offset += n_frags
            current_local_idx += 1
        end
    end

    # Unload final batch
    unload_current_batch!(batched_lib)

    # Create reduced library in memory (no disk I/O needed)
    fragment_lookup = StandardFragmentLookup(all_detailed_frags, prec_frag_ranges)

    # Create precursors table from collected rows
    prec_table = DataFrame(all_precursor_rows)
    precursors = SetPrecursors(Arrow.Table(prec_table))

    return ReducedSpectralLibrary(
        precursors,
        batched_lib.proteins,
        fragment_lookup,
        original_prec_idx_to_local,
        local_to_original_prec_idx
    )
end
```

### Integrate into FirstPassSearch

**File: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl`**

In `summarize_results!`, after `get_best_precursors_accross_runs`:

```julia
function summarize_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:FirstPassSearchParameters}

    # ... existing code for RT mapping, precursor selection ...

    # Create reduced library for batched searches
    if isBatchedLibrary(search_context)
        precursor_dict = getPrecursorDict(search_context)
        passing_ids = collect(keys(precursor_dict))

        @info "Creating reduced library from $(length(passing_ids)) passing precursors"

        pioneer_params = getParams(search_context)  # Get full params
        reduced_lib = create_reduced_library(search_context, passing_ids, pioneer_params)
        setReducedLibrary!(search_context, reduced_lib)

        @info "Reduced library created successfully"
    end

    # ... rest of existing code ...
end
```

---

## Phase 7: Use Reduced Library in Later Passes

### File: `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/SecondPassSearch.jl`

Replace `getSpecLib(search_context)` with `getActiveLibrary(search_context)`:

```julia
function process_file!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:SecondPassSearchParameters}

    # Use reduced library if available
    spec_lib = getActiveLibrary(search_context)

    # ... rest of existing code, using spec_lib ...
end
```

Apply same pattern to:
- `ScoringSearch/ScoringSearch.jl`
- `IntegrateChromatogramSearch/IntegrateChromatogramSearch.jl`
- `MaxLFQSearch/MaxLFQSearch.jl`

---

## Phase 8: Library Building (Future Enhancement)

### File: `src/Routines/BuildSpecLib/build/build_poin_lib.jl`

```julia
"""
    partition_precursors_by_mz(
        precursor_mzs::Vector{Float32},
        n_batches::Int
    ) -> Tuple{Vector{Tuple{Float32, Float32}}, Vector{Vector{UInt32}}}

Partition precursors into N batches based on m/z ranges.
Returns batch m/z ranges and precursor indices per batch.
"""
function partition_precursors_by_mz(precursor_mzs::Vector{Float32}, n_batches::Int)
    n_precs = length(precursor_mzs)
    sorted_indices = sortperm(precursor_mzs)
    batch_size = ceil(Int, n_precs / n_batches)

    batch_ranges = Vector{Tuple{Float32, Float32}}(undef, n_batches)
    batch_precursor_indices = Vector{Vector{UInt32}}(undef, n_batches)

    for batch_idx in 1:n_batches
        start_idx = (batch_idx - 1) * batch_size + 1
        end_idx = min(batch_idx * batch_size, n_precs)

        if start_idx > n_precs
            # Empty batch
            batch_ranges[batch_idx] = (Inf32, -Inf32)
            batch_precursor_indices[batch_idx] = UInt32[]
            continue
        end

        batch_indices = sorted_indices[start_idx:end_idx]
        mz_low = precursor_mzs[batch_indices[1]]
        mz_high = precursor_mzs[batch_indices[end]]

        batch_ranges[batch_idx] = (mz_low, mz_high)
        batch_precursor_indices[batch_idx] = UInt32.(batch_indices)
    end

    return batch_ranges, batch_precursor_indices
end

"""
    buildBatchedPionLib(
        spec_lib_path::String,
        n_batches::Int;
        kwargs...
    )

Build a batched Pioneer spectral library from preprocessed data.
"""
function buildBatchedPionLib(
    spec_lib_path::String,
    n_batches::Int;
    rt_bin_tol::Float64 = 0.1,
    frag_bin_tol_ppm::Float64 = 10.0,
    kwargs...
)
    @info "Building batched library with $n_batches batches"

    # Load precursor data
    precursors_table = Arrow.Table(joinpath(spec_lib_path, "precursors_table.arrow"))
    precursor_mzs = Vector{Float32}(precursors_table[:mz])

    # Partition precursors
    batch_ranges, batch_precursor_indices = partition_precursors_by_mz(precursor_mzs, n_batches)

    # Save library config
    config = BatchedLibraryConfig(n_batches, batch_ranges, spec_lib_path, true)
    save_library_config(joinpath(spec_lib_path, "library_config.json"), config)

    # Load fragment data
    fragments_table = Arrow.Table(joinpath(spec_lib_path, "fragments_table.arrow"))
    prec_to_frag = Arrow.Table(joinpath(spec_lib_path, "prec_to_frag.arrow"))

    # Build each batch
    for batch_idx in 1:n_batches
        batch_path = joinpath(spec_lib_path, "batch_$batch_idx")
        mkpath(batch_path)

        batch_prec_indices = batch_precursor_indices[batch_idx]

        if isempty(batch_prec_indices)
            @warn "Batch $batch_idx is empty, skipping"
            continue
        end

        @info "Building batch $batch_idx with $(length(batch_prec_indices)) precursors"

        # Filter and build batch
        _build_single_batch(
            batch_path,
            batch_prec_indices,
            precursors_table,
            fragments_table,
            prec_to_frag,
            rt_bin_tol,
            frag_bin_tol_ppm
        )
    end

    @info "Batched library build complete"
end
```

---

## Critical Files Summary

| File | Changes |
|------|---------|
| `src/structs/LibraryFragmentIndex.jl` | Add `BatchedLibraryConfig`, `BatchedSpectralLibrary` types and accessors |
| `src/Routines/SearchDIA/ParseInputs/loadSpectralLibrary.jl` | Add batch detection, `loadBatchedSpectralLibrary`, `load_batch!`, `unload_current_batch!` |
| `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl` | Add `reduced_library` field, `isBatchedLibrary`, `getActiveLibrary` |
| `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/ParameterTuningSearch.jl` | Add `search_all_batches` for batch iteration |
| `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl` | Add batch iteration + trigger reduced library creation |
| **NEW**: `src/Routines/SearchDIA/CommonSearchUtils/createReducedLibrary.jl` | Add `create_reduced_library` and helpers |
| `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/SecondPassSearch.jl` | Use `getActiveLibrary` instead of `getSpecLib` |
| `src/Routines/BuildSpecLib/build/build_poin_lib.jl` | Add `buildBatchedPionLib` (future) |

---

## Verification

1. **Unit test**: Create small batched library (2 batches), verify loading/unloading
2. **Integration test**: Run full SearchDIA with batched library, compare PSM counts to monolithic
3. **Memory test**: Profile memory usage during batch iteration vs. monolithic
4. **Reduced library test**: Verify reduced library contains only passing precursors
5. **Backward compatibility**: Run existing tests with monolithic libraries
