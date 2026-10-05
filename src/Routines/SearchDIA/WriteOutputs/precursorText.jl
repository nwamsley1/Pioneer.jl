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

# Intermediate PSM tables carry no library text: sequence, modifications, accessions, species and
# start positions are functions of :precursor_idx, and the file name of :ms_file_idx. The final
# precursor writers put those columns back here, at the positions they held when they were
# stored on the PSMs (around :peak_area_normalized), so the outputs are unchanged.

"""
    precursor_text_column(text, name, precursor_idx, ms_file_idx) -> Vector

One text column for the given rows.
"""
function precursor_text_column(text::PrecursorOutputText, name::Symbol,
                               precursor_idx::AbstractVector, ms_file_idx::AbstractVector)
    p = text.precursors
    if name === :accession_numbers
        accessions = getAccessionNumbers(p)
        return String[accessions[pid] for pid in precursor_idx]
    elseif name === :species
        ids = text.text_ids
        return String[ids.species_names[ids.species_id[pid]] for pid in precursor_idx]
    elseif name === :structural_mods
        mods = getStructuralMods(p)
        return Union{Missing, String}[mods[pid] for pid in precursor_idx]
    elseif name === :isotopic_mods
        mods = getIsotopicMods(p)
        return Union{Missing, String}[mods[pid] for pid in precursor_idx]
    elseif name === :sequence
        sequences = getSequence(p)
        return String[sequences[pid] for pid in precursor_idx]
    elseif name === :peptide_start_positions
        starts = getStartIdx(p)
        return String[_format_start_idx(starts[pid]) for pid in precursor_idx]
    elseif name === :file_name
        return String[text.file_names[idx] for idx in ms_file_idx]
    end
    throw(ArgumentError("Not a precursor text column: $name"))
end

"""
    with_precursor_text(cols::NamedTuple, text) -> NamedTuple

Insert the text columns into a precursor column table: `accession_numbers` and `species` just
before `:peak_area_normalized`, the rest just after it.
"""
function with_precursor_text(cols::NamedTuple, text::PrecursorOutputText)
    haskey(cols, :peak_area_normalized) ||
        throw(ArgumentError("Precursor table lacks :peak_area_normalized, the text column anchor"))
    precursor_idx, ms_file_idx = cols.precursor_idx, cols.ms_file_idx
    text_col(name) = precursor_text_column(text, name, precursor_idx, ms_file_idx)
    names = Symbol[]
    values = Any[]
    for (name, col) in pairs(cols)
        name in PRECURSOR_TEXT_COLUMNS && continue   # a stale copy; the library is authoritative
        if name === :peak_area_normalized
            for t in PRECURSOR_TEXT_COLUMNS_BEFORE_NORMALIZED
                push!(names, t); push!(values, text_col(t))
            end
            push!(names, name); push!(values, col)
            for t in PRECURSOR_TEXT_COLUMNS_AFTER_NORMALIZED
                push!(names, t); push!(values, text_col(t))
            end
        else
            push!(names, name); push!(values, col)
        end
    end
    return NamedTuple{Tuple(names)}(Tuple(values))
end

"""
    precursor_text_pools(chunk_refs, text, names) -> Dict{Symbol, Vector}

Sorted distinct values of the given text columns over every row of the chunks, read from the
chunks' `:precursor_idx` and `:ms_file_idx` alone.
"""
function precursor_text_pools(chunk_refs, text::PrecursorOutputText, names)
    pids = Set{UInt32}()
    files = Set{Int}()
    for chunk_ref in chunk_refs
        tbl = Arrow.Table(file_path(chunk_ref))
        union!(pids, tbl.precursor_idx)
        union!(files, tbl.ms_file_idx)
    end
    pid_vec, file_vec = collect(pids), collect(files)
    return Dict(name => sort!(collect(Set(precursor_text_column(
                    text, name,
                    name === :file_name ? UInt32[] : pid_vec,
                    name === :file_name ? file_vec : Int[]))))
                for name in names)
end
