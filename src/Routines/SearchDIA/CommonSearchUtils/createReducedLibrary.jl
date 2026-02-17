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
Functions for creating a reduced spectral library from passing precursors.
Used after FirstPassSearch to create a memory-efficient library for subsequent passes.

NOTE: The reduced library does NOT have a fragment index since it's only used
after FirstPassSearch, which is the last stage requiring the fragment index.
"""

"""
    create_reduced_library(
        original_lib::SpectralLibrary,
        passing_precursor_ids::Vector{UInt32},
        params::PioneerParameters
    ) -> ReducedSpectralLibrary

Create a reduced spectral library containing only the specified precursors.
Extracts precursor data and detailed fragments from the original library.
No fragment index is built since subsequent stages don't need it.
"""
function create_reduced_library(
    original_lib::SpectralLibrary,
    passing_precursor_ids::Vector{UInt32},
    params::PioneerParameters
)::ReducedSpectralLibrary

    @info "Creating reduced library with $(length(passing_precursor_ids)) precursors (no fragment index needed)"

    if isa(original_lib, BatchedSpectralLibrary)
        return _create_reduced_from_batched(original_lib, passing_precursor_ids, params)
    else
        return _create_reduced_from_monolithic(original_lib, passing_precursor_ids)
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
    all_detailed_frags = Vector{DetailedFrag{Float32}}()
    prec_frag_ranges = Vector{UInt64}()

    # Build bidirectional index mapping
    original_prec_idx_to_local = Dictionary{UInt32, UInt32}()
    local_to_original_prec_idx = Vector{UInt32}()

    # Collect precursor data - we'll build the Arrow table at the end
    precursor_data_rows = NamedTuple[]

    current_frag_offset = UInt64(1)
    current_local_idx = UInt32(1)

    for batch_idx in sort(collect(keys(batch_to_precs)))
        batch_prec_ids = batch_to_precs[batch_idx]

        @info "Extracting $(length(batch_prec_ids)) precursors from batch $batch_idx"

        # Load batch
        load_batch!(batched_lib, batch_idx, params)
        batch = batched_lib.current_batch

        # Get batch data
        batch_precursors = getPrecursors(batch)
        batch_frag_lookup = getFragmentLookupTable(batch)
        batch_frags = getFragments(batch_frag_lookup)

        # Extract data for each passing precursor in this batch
        for original_prec_id in batch_prec_ids
            # Get local index within batch using O(1) arithmetic
            batch_local_idx = get_local_idx(config, original_prec_id)

            # Build index mapping
            insert!(original_prec_idx_to_local, original_prec_id, current_local_idx)
            push!(local_to_original_prec_idx, original_prec_id)

            # Collect precursor row data
            push!(precursor_data_rows, _extract_precursor_row(batch_precursors, batch_local_idx))

            # Get fragment range for this precursor
            frag_range = getPrecFragRange(batch_frag_lookup, batch_local_idx)
            n_frags = length(frag_range)

            # Record the start position for this precursor's fragments
            push!(prec_frag_ranges, current_frag_offset)

            # Copy fragments
            for i in frag_range
                push!(all_detailed_frags, batch_frags[i])
            end

            current_frag_offset += n_frags
            current_local_idx += one(UInt32)
        end
    end

    # Add final sentinel value for fragment ranges
    push!(prec_frag_ranges, current_frag_offset)

    # Unload final batch
    unload_current_batch!(batched_lib)

    # Create reduced library in memory
    fragment_lookup = StandardFragmentLookup(all_detailed_frags, prec_frag_ranges)

    # Create precursors table from collected rows
    precursors = _build_reduced_precursors(precursor_data_rows)

    return ReducedSpectralLibrary(
        precursors,
        batched_lib.proteins,
        fragment_lookup,
        original_prec_idx_to_local,
        local_to_original_prec_idx
    )
end

"""
    _create_reduced_from_monolithic(
        lib::SpectralLibrary,
        passing_ids::Vector{UInt32}
    ) -> ReducedSpectralLibrary

Extract passing precursors from a monolithic (non-batched) library.
"""
function _create_reduced_from_monolithic(
    lib::SpectralLibrary,
    passing_ids::Vector{UInt32}
)::ReducedSpectralLibrary

    precursors = getPrecursors(lib)
    frag_lookup = getFragmentLookupTable(lib)
    frags = getFragments(frag_lookup)

    # Collect data
    all_detailed_frags = Vector{DetailedFrag{Float32}}()
    prec_frag_ranges = Vector{UInt64}()

    original_prec_idx_to_local = Dictionary{UInt32, UInt32}()
    local_to_original_prec_idx = Vector{UInt32}()
    precursor_data_rows = NamedTuple[]

    current_frag_offset = UInt64(1)
    current_local_idx = UInt32(1)

    for original_prec_id in passing_ids
        # Build index mapping
        insert!(original_prec_idx_to_local, original_prec_id, current_local_idx)
        push!(local_to_original_prec_idx, original_prec_id)

        # Collect precursor row data
        push!(precursor_data_rows, _extract_precursor_row(precursors, Int(original_prec_id)))

        # Get fragment range for this precursor
        frag_range = getPrecFragRange(frag_lookup, Int(original_prec_id))
        n_frags = length(frag_range)

        # Record the start position for this precursor's fragments
        push!(prec_frag_ranges, current_frag_offset)

        # Copy fragments
        for i in frag_range
            push!(all_detailed_frags, frags[i])
        end

        current_frag_offset += n_frags
        current_local_idx += one(UInt32)
    end

    # Add final sentinel value for fragment ranges
    push!(prec_frag_ranges, current_frag_offset)

    # Create reduced library in memory
    fragment_lookup_new = StandardFragmentLookup(all_detailed_frags, prec_frag_ranges)

    # Create precursors table from collected rows
    reduced_precursors = _build_reduced_precursors(precursor_data_rows)

    return ReducedSpectralLibrary(
        reduced_precursors,
        getProteins(lib),
        fragment_lookup_new,
        original_prec_idx_to_local,
        local_to_original_prec_idx
    )
end

"""
    _extract_precursor_row(precursors::LibraryPrecursors, idx::Int) -> NamedTuple

Extract a single precursor's data as a NamedTuple for later table construction.
"""
function _extract_precursor_row(precursors::LibraryPrecursors, idx::Int)
    data = precursors.data

    return (
        sequence = data[:sequence][idx],
        structural_mods = data[:structural_mods][idx],
        isotopic_mods = data[:isotopic_mods][idx],
        charge = data[:charge][idx],
        mz = data[:mz][idx],
        irt = data[:irt][idx],
        prec_charge = data[:prec_charge][idx],
        sulfur_count = data[:sulfur_count][idx],
        missed_cleavages = data[:missed_cleavages][idx],
        is_decoy = data[:is_decoy][idx],
        entrapment_group_id = data[:entrapment_group_id][idx],
        accession_numbers = data[:accession_numbers][idx],
        gene_names = haskey(data, :gene_names) ? data[:gene_names][idx] : "",
        proteome_identifiers = haskey(data, :proteome_identifiers) ? data[:proteome_identifiers][idx] : ""
    )
end

"""
    _build_reduced_precursors(rows::Vector{NamedTuple}) -> LibraryPrecursors

Build a StandardLibraryPrecursors from collected row data.
"""
function _build_reduced_precursors(rows::Vector{NamedTuple})
    if isempty(rows)
        error("Cannot build precursors from empty row list")
    end

    # Build DataFrame from rows
    df = DataFrame(rows)

    # Convert to Arrow table and create StandardLibraryPrecursors
    arrow_table = Arrow.Table(df)
    return SetPrecursors(arrow_table)
end
