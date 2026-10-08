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

# precursors.arrow -> the library's precursors_table.arrow, one record batch at a time: the search-side column names,
# partner_precursor_idx (add_pair_indices!) and, with entrapments, entrapment_target_idx (add_entrapment_indices!).

const FINAL_PRECURSOR_NAMES = Dict(:accession_number => :accession_numbers, :precursor_charge => :prec_charge,
                                   :decoy => :is_decoy, :mods => :structural_mods)

"""
    finalize_precursor_table(in_path, out_path; entrapment_targets) -> (n_precursors, n_decoys)

Write `out_path` from the sorted precursor table at `in_path` with the contents of the in-memory step (rename,
`add_pair_indices!`, and `add_entrapment_indices!` when `entrapment_targets`), streaming record batches. Rows keep
their order; the indices are row positions in it.
"""
function finalize_precursor_table(in_path::String, out_path::String; entrapment_targets::Bool)
    t = Arrow.Table(in_path)
    pair_ids = t.pair_id
    partner = pair_partners(pair_ids)
    etarget = entrapment_targets && hasproperty(t, :entrapment_pair_id) ?
        entrapment_target_rows(t.entrapment_pair_id, t.entrapment_group_id, t.decoy) : nothing
    n_decoys = count(t.decoy)
    n = length(pair_ids)
    t = nothing
    writer = open(Arrow.Writer, out_path)
    try
        row0 = 0
        for batch in Arrow.Stream(in_path)
            m = length(batch.pair_id)
            Arrow.write(writer, _final_precursor_batch(batch, partner, etarget, row0))
            row0 += m
        end
    finally
        close(writer)
    end
    return n, n_decoys
end

"One record batch in the final schema: renamed columns, materialized, then the index columns."
function _final_precursor_batch(batch, partner::Vector{Int64}, etarget::Union{Nothing, Vector{UInt32}}, row0::Int)
    names = Symbol[]; cols = AbstractVector[]
    for (name, col) in zip(Tables.columnnames(batch), Tables.columns(batch))
        push!(names, get(FINAL_PRECURSOR_NAMES, name, name))
        push!(cols, name === :start_idx ? [UInt32.(collect(s)) for s in col] : collect(col))
    end
    m = length(first(cols))
    push!(names, :partner_precursor_idx)
    push!(cols, Union{Missing, Int64}[partner[row0 + i] == 0 ? missing : partner[row0 + i] for i in 1:m])
    if etarget !== nothing
        push!(names, :entrapment_target_idx)
        push!(cols, Union{Missing, UInt32}[etarget[row0 + i] == 0 ? missing : etarget[row0 + i] for i in 1:m])
    end
    return NamedTuple{Tuple(names)}(Tuple(cols))
end

"""
    pair_partners(pair_ids) -> Vector{Int64}

add_pair_indices!: the row of the other member of each row's pair_id group when the group has exactly two rows,
else 0.
"""
function pair_partners(pair_ids::AbstractVector)
    n = length(pair_ids)
    any(ismissing, pair_ids) && error("precursor table has rows without a pair_id")
    max_id = n == 0 ? 0 : Int(maximum(pair_ids))
    first_row = zeros(UInt32, max_id); second_row = zeros(UInt32, max_id); n_rows = zeros(UInt8, max_id)
    for (i, p) in enumerate(pair_ids)
        c = n_rows[p]
        c == 0 ? (first_row[p] = UInt32(i)) : c == 1 && (second_row[p] = UInt32(i))
        n_rows[p] = min(c + 0x01, 0x03)
    end
    partner = zeros(Int64, n)
    for (i, p) in enumerate(pair_ids)
        n_rows[p] == 2 && (partner[i] = first_row[p] == i ? second_row[p] : first_row[p])
    end
    return partner
end

"""
    entrapment_target_rows(entrapment_pair_ids, entrapment_group_ids, decoys) -> Vector{UInt32}

add_entrapment_indices!: for each non-decoy row with an entrapment_pair_id, the row of the (last) group-0 non-decoy
row with that id, else 0.
"""
function entrapment_target_rows(epair::AbstractVector, group::AbstractVector, decoy::AbstractVector)
    n = length(epair)
    max_id = 0
    for e in epair; ismissing(e) || (max_id = max(max_id, Int(e))); end
    target = zeros(UInt32, max_id)
    for (i, (e, g, d)) in enumerate(zip(epair, group, decoy))
        (!ismissing(e) && g == 0 && !d) && (target[e] = UInt32(i))
    end
    rows = zeros(UInt32, n)
    for (i, (e, d)) in enumerate(zip(epair, decoy))
        (!ismissing(e) && !d) && (rows[i] = target[e])
    end
    return rows
end
