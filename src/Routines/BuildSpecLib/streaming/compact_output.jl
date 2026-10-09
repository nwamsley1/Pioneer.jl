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

# The streaming builder's schema-2 output (precursor_table_v2.jl): packed sequences and mod entries per row, written
# from chunk-level flat buffers (no per-row strings or arrays), and the side tables per base_pep_id.

"""
    RaggedColumn(data, offsets)

Row `i` is `view(data, offsets[i]:offsets[i+1]-1)`: a list column backed by one flat vector. Arrow writes it as a
list column without materializing a vector per row.
"""
struct RaggedColumn{T} <: AbstractVector{SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int}}, true}}
    data::Vector{T}
    offsets::Vector{Int}
end
Base.size(v::RaggedColumn) = (length(v.offsets) - 1,)
Base.@propagate_inbounds Base.getindex(v::RaggedColumn, i::Int) = view(v.data, v.offsets[i]:(v.offsets[i + 1] - 1))

"Schema-2 residue code (letter - 'A' + 1) of each SeqCode residue code (index into SEQ_ALPHABET)."
const _SEQCODE_TO_PACKED = UInt8[UInt8(c - 'A' + 1) for c in SEQ_ALPHABET]

"Write the packed form (pack_sequence) of `c` into `buf[off+1 : off+cld(5L,8)]` (zeroed)."
function pack_seqcode!(buf::Vector{UInt8}, off::Int, c::SeqCode)
    L = seq_length(c)
    for i in 1:L
        sc = i <= 25 ? (c.hi >> (5 * (25 - i) + 3)) & 0x1f : (c.lo >> (5 * (50 - i) + 3)) & 0x1f
        code = _SEQCODE_TO_PACKED[Int(sc)]
        bit = 5 * (i - 1)
        for b in 0:4
            ((code >> (4 - b)) & 0x01) == 0x01 || continue
            k = bit + b
            buf[off + k >> 3 + 1] |= 0x80 >> (k & 7)
        end
    end
    return nothing
end

"Number of residues of a packed SeqCode that are C or M (the schema-1 sulfur_count of a sequence)."
function sulfur_residues(c::SeqCode)
    n = 0
    for i in 1:seq_length(c)
        sc = i <= 25 ? (c.hi >> (5 * (25 - i) + 3)) & 0x1f : (c.lo >> (5 * (50 - i) + 3)) & 0x1f
        r = SEQ_ALPHABET[Int(sc)]
        n += (r == 'C') | (r == 'M')
    end
    return n
end

function _compact_columns_alloc(n::Int)
    return (sequence_packed = RaggedColumn(UInt8[], Vector{Int}(undef, n + 1)),
        mod_entries = RaggedColumn(UInt16[], Vector{Int}(undef, n + 1)),
        num_variable_modifications = Vector{UInt8}(undef, n), precursor_charge = Vector{UInt8}(undef, n),
        num_enzymatic_termini = Vector{UInt8}(undef, n), decoy = Vector{Bool}(undef, n),
        entrapment_group_id = Vector{UInt8}(undef, n), base_target_id = Vector{UInt32}(undef, n),
        base_pep_id = Vector{UInt32}(undef, n), pair_id = Vector{Union{Missing, UInt32}}(undef, n),
        mz = Vector{Float32}(undef, n), length = Vector{UInt8}(undef, n), missed_cleavages = Vector{UInt8}(undef, n),
        entrapment_pair_id = Vector{Union{Missing, UInt32}}(undef, n), irt = Vector{Float32}(undef, n),
        sulfur_count = Vector{UInt8}(undef, n))
end

"""
    chunk_columns_compact(src, chunk)

The schema-2 intermediate columns of rows `chunk` (the schema-1 columns without the text: protein columns, Koina
sequence, collision energy and isotope mods are gone, sequence and mods are packed), built in parallel: a pass for the
per-row sizes, a prefix sum, then a fill into the flat buffers.
"""
function chunk_columns_compact(src::TableSource, chunk::AbstractVector{UInt32})
    n = length(chunk)
    cols = _compact_columns_alloc(n)
    seq_off = cols.sequence_packed.offsets; mod_off = cols.mod_entries.offsets
    parts = collect(Iterators.partition(1:n, 4096))
    Threads.@threads :dynamic for js in parts
        _compact_sizes!(seq_off, mod_off, src, chunk, js)
    end
    seq_off[1] = 1; mod_off[1] = 1
    for j in 1:n
        seq_off[j + 1] += seq_off[j]; mod_off[j + 1] += mod_off[j]
    end
    resize!(cols.sequence_packed.data, seq_off[end] - 1); fill!(cols.sequence_packed.data, 0x00)
    resize!(cols.mod_entries.data, mod_off[end] - 1)
    Threads.@threads :dynamic for js in parts
        _fill_compact_rows!(cols, src, chunk, js, Tuple{UInt8, UInt8}[], UInt8[])
    end
    return cols
end

"Row sizes (stored at j + 1, prefix-summed by the caller): packed sequence bytes and mod count."
function _compact_sizes!(seq_off::Vector{Int}, mod_off::Vector{Int}, src::TableSource, chunk::AbstractVector{UInt32},
                         js::UnitRange{Int})
    units = src.units
    for j in js
        u = (Int(chunk[j]) - 1) ÷ src.nz + 1
        seq_off[j + 1] = cld(5 * seq_length(unit_code(units, u)), 8)
        f = unit_target(units, u)
        mod_off[j + 1] = Int(units.mod_off[f + 1] - units.mod_off[f])
    end
    return nothing
end

function _fill_compact_rows!(cols::NamedTuple, src::TableSource, chunk::AbstractVector{UInt32}, js::UnitRange{Int},
                             tbuf::Vector{Tuple{UInt8, UInt8}}, trev::Vector{UInt8})
    units = src.units; peps = units.peps
    seq_data = cols.sequence_packed.data; seq_off = cols.sequence_packed.offsets
    mod_data = cols.mod_entries.data; mod_off = cols.mod_entries.offsets
    for j in js
        r = chunk[j]
        u = (Int(r) - 1) ÷ src.nz + 1; zi = (Int(r) - 1) % src.nz + 1
        k = unit_peptide(units, u)
        code = unit_code(units, u)
        pack_seqcode!(seq_data, seq_off[j] - 1, code)
        mods, nvar = unit_mods!(tbuf, trev, units, u)
        o = mod_off[j] - 1
        for (i, (pos, id)) in enumerate(mods)
            mod_data[o + i] = mod_entry(pos, MOD_SITE_RESIDUE, id)
        end
        L = seq_length(code)
        cols.num_variable_modifications[j] = UInt8(nvar)
        cols.precursor_charge[j] = src.charges[zi]
        cols.num_enzymatic_termini[j] = peps.nte[k]
        decoy = is_decoy_unit(units, u)
        cols.decoy[j] = decoy
        cols.entrapment_group_id[j] = entrapment_group(units, u)
        cols.base_target_id[j] = base_target_id(units, u)
        cols.base_pep_id[j] = base_pep_id(units, u)
        cols.pair_id[j] = src.pair_id[r]
        cols.mz[j] = src.row_mz[r]
        cols.length[j] = UInt8(L)
        cols.missed_cleavages[j] = src.cleavage === nothing ? 0x00 : UInt8(count(src.cleavage, decode_seq(code)))
        cols.entrapment_pair_id[j] = decoy || src.epair_id[r] == 0 ? missing : src.epair_id[r]
        cols.irt[j] = src.unit_irt[u]
        cols.sulfur_count[j] = UInt8(sulfur_residues(code))
    end
    return nothing
end

"""
    write_precursor_side_tables(lib_dir, units, rows, nz)

The schema-2 side tables of the table's rows (precursor_table_v2.jl): one row per base_pep_id up to the largest one in
the table (ids without rows hold zeros), distinct accession sets and proteome strings in sort order, accessions in sort
order, and the mod names. A peptide's accession set and proteome string are its occurrences' accessions / proteomes
joined with ';' (most recent occurrence first), its starts the occurrences' start positions in the same order.
"""
function write_precursor_side_tables(lib_dir::AbstractString, units::StreamUnits, rows::Vector{UInt32}, nz::Int)
    peps = units.peps
    n_pep = 0
    used = falses(units.F)
    for r in rows
        b = base_pep_id(units, (Int(r) - 1) ÷ nz + 1)
        used[b] = true; n_pep = max(n_pep, Int(b))
    end
    # distinct protein lists -> accession / proteome strings (single-protein peptides need no joined string)
    set_of = Dict{String, Int}(); set_names = String[]; proteome_of_set = String[]
    pep_set_name = zeros(Int, n_pep)                                     # index into set_names, per base_pep_id
    for b in 1:n_pep
        used[b] || continue
        k = Int(units.unit_pep[b])
        occ = peps.occ_offsets[k]:(peps.occ_offsets[k + 1] - 1)
        name = join_occurrences(peps.protein_accession, peps.occ_protein, occ)
        pep_set_name[b] = get!(set_of, name) do
            push!(set_names, name)
            push!(proteome_of_set, join_occurrences(peps.protein_proteome, peps.occ_protein, occ))
            length(set_names)
        end
    end
    order = sortperm(set_names)
    rank = invperm(order)
    sorted_sets = set_names[order]
    acc_names = sort!(unique(String[a for s in sorted_sets for a in split(s, ';')]))
    acc_id = Dict(a => UInt32(k) for (k, a) in enumerate(acc_names))
    proteome_names = sort!(unique(proteome_of_set))
    length(proteome_names) <= typemax(UInt16) || error("more than $(typemax(UInt16)) distinct proteome strings")
    proteome_id = Dict(s => UInt16(k) for (k, s) in enumerate(proteome_names))
    pep_set = zeros(UInt32, n_pep); pep_prot = zeros(UInt16, n_pep)
    start_off = ones(Int, n_pep + 1); starts = UInt32[]
    for b in 1:n_pep
        if used[b]
            s = pep_set_name[b]
            pep_set[b] = UInt32(rank[s]); pep_prot[b] = proteome_id[proteome_of_set[s]]
            k = Int(units.unit_pep[b])
            for o in reverse(peps.occ_offsets[k]:(peps.occ_offsets[k + 1] - 1))
                push!(starts, peps.occ_start[o])
            end
        end
        start_off[b + 1] = length(starts) + 1
    end
    f = PRECURSOR_SIDE_FILES
    Arrow.write(joinpath(lib_dir, f.peptides), (accession_set_id = pep_set, proteome_id = pep_prot,
                                                start_idx = RaggedColumn(starts, start_off)))
    Arrow.write(joinpath(lib_dir, f.accession_sets), (accession_numbers = sorted_sets,
        members = [sort!(unique!(UInt32[acc_id[a] for a in split(s, ';')])) for s in sorted_sets]))
    Arrow.write(joinpath(lib_dir, f.accessions), (accession = acc_names,))
    Arrow.write(joinpath(lib_dir, f.proteomes), (proteome_identifiers = proteome_names,))
    Arrow.write(joinpath(lib_dir, f.mod_names), (name = units.mod_names,))
    return nothing
end
