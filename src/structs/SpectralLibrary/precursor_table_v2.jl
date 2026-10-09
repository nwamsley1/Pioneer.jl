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

# Compact precursor table, schema 2 (dev_docs/precursor_table_compaction/plan.pdf).
#
# precursors_table.arrow (metadata pioneer_precursor_schema = "2") replaces three kinds of per-row text:
#   sequence        -> sequence_packed: 5-bit residue codes (letter - 'A' + 1), left-aligned, zero padded, so byte order
#                      is String order and the length is the number of non-zero codes
#   structural_mods -> mod_entries: one UInt16 per mod, in the string's order: position << 10 | site << 8 | name id
#                      (site 0 = the residue at the position, 1 = 'n', 2 = 'c'; names in precursor_mod_names.arrow)
#   accession_numbers, proteome_identifiers, start_idx -> nothing per row: base_pep_id determines them, and
#                      precursor_peptides.arrow holds them per base_pep_id (accession-set and proteome IDs, start lists)
# Side tables: precursor_peptides.arrow, precursor_accession_sets.arrow (distinct sets in sort order, with member
# accession IDs), precursor_accessions.arrow, precursor_proteomes.arrow, precursor_mod_names.arrow.
# The getters decode back to the schema-1 values, so readers see the same columns.

const PRECURSOR_SCHEMA_KEY = "pioneer_precursor_schema"
const PRECURSOR_SIDE_FILES = (peptides = "precursor_peptides.arrow", accession_sets = "precursor_accession_sets.arrow",
                              accessions = "precursor_accessions.arrow", proteomes = "precursor_proteomes.arrow",
                              mod_names = "precursor_mod_names.arrow")

precursor_schema(tbl::Arrow.Table) = (m = Arrow.getmetadata(tbl); m === nothing ? 1 : parse(Int, get(m, PRECURSOR_SCHEMA_KEY, "1")))

# ── packed sequences ──────────────────────────────────────────────────────────────────────────────────────────────

"Residue code of an uppercase letter: 1-26 ('A' = 1), in alphabetical order."
@inline function residue_code(c::Char)
    'A' <= c <= 'Z' || error("packed sequences hold uppercase letters only, got '$c'")
    return UInt8(c - 'A' + 1)
end

"The 5-bit codes of `seq`, left-aligned and zero padded in ceil(5L/8) bytes."
function pack_sequence(seq::AbstractString)
    L = ncodeunits(seq)
    bytes = zeros(UInt8, cld(5 * L, 8))
    for (i, c) in enumerate(seq)
        code = residue_code(c)
        bit = 5 * (i - 1)                                   # first bit (from the most significant end)
        for b in 0:4
            ((code >> (4 - b)) & 0x01) == 0x01 || continue
            k = bit + b
            bytes[k >> 3 + 1] |= 0x80 >> (k & 7)
        end
    end
    return bytes
end

"Residue code at position `i` (1-based) of a packed sequence; 0 past its end."
@inline function packed_code(bytes::AbstractVector{UInt8}, i::Int)
    bit = 5 * (i - 1)
    bit + 5 > 8 * length(bytes) && return 0x00
    code = 0x00
    @inbounds for b in 0:4
        k = bit + b
        code = (code << 1) | ((bytes[k >> 3 + 1] >> (7 - (k & 7))) & 0x01)
    end
    return code
end

"Number of residues of a packed sequence."
function packed_length(bytes::AbstractVector{UInt8})
    n = 0
    while packed_code(bytes, n + 1) != 0x00
        n += 1
    end
    return n
end

"The sequence of a packed sequence."
function unpack_sequence(bytes::AbstractVector{UInt8})
    L = packed_length(bytes)
    v = Base.StringVector(L)
    for i in 1:L
        v[i] = UInt8('A') + packed_code(bytes, i) - 0x01
    end
    return String(v)
end

# ── mod entries ───────────────────────────────────────────────────────────────────────────────────────────────────

const MOD_SITE_RESIDUE = 0x00
const MOD_SITE_NTERM = 0x01
const MOD_SITE_CTERM = 0x02

@inline mod_entry(pos::Integer, site::UInt8, name_id::Integer) = UInt16(pos) << 10 | UInt16(site) << 8 | UInt16(name_id)
@inline mod_entry_position(e::UInt16) = Int(e >> 10)
@inline mod_entry_site(e::UInt16) = UInt8((e >> 8) & 0x03)
@inline mod_entry_name_id(e::UInt16) = Int(e & 0x00ff)

"""
    encode_mods(mods, seq, name_id) -> Vector{UInt16} or missing

The mod entries of a `structural_mods` string ("(pos,aa,name)..."), in the string's order. `name_id` maps mod names to
IDs and is extended with new names. Errors if a mod's residue is neither the sequence's residue at its position nor
'n' / 'c', or a position or name does not fit its field.
"""
encode_mods(::Missing, seq::AbstractString, name_id::Dict{String, UInt8}) = missing
function encode_mods(mods::AbstractString, seq::AbstractString, name_id::Dict{String, UInt8})
    entries = UInt16[]
    for m in eachmatch(r"\((\d+),(.),([^)]*)\)", mods)
        pos = parse(Int, m.captures[1]); aa = m.captures[2][1]; name = String(m.captures[3])
        1 <= pos <= 63 || error("mod position $pos does not fit the 6-bit field ($mods)")
        site = aa == 'n' ? MOD_SITE_NTERM : aa == 'c' ? MOD_SITE_CTERM :
               aa == seq[pos] ? MOD_SITE_RESIDUE : error("mod residue '$aa' at $pos is not $(seq[pos]) in $seq ($mods)")
        id = get!(name_id, name) do
            length(name_id) < 255 || error("more than 255 mod names")
            UInt8(length(name_id) + 1)
        end
        push!(entries, mod_entry(pos, site, id))
    end
    return entries
end

"The `structural_mods` string of mod entries on a packed sequence."
decode_mods(::Missing, ::AbstractVector{UInt8}, ::Vector{String}) = missing
function decode_mods(entries::AbstractVector{UInt16}, packed::AbstractVector{UInt8}, names::Vector{String})
    isempty(entries) && return ""
    io = IOBuffer()
    for e in entries
        pos = mod_entry_position(e); site = mod_entry_site(e)
        aa = site == MOD_SITE_NTERM ? 'n' : site == MOD_SITE_CTERM ? 'c' : Char(UInt8('A') + packed_code(packed, pos) - 0x01)
        print(io, '(', pos, ',', aa, ',', names[mod_entry_name_id(e)], ')')
    end
    return String(take!(io))
end

# ── list columns from flat buffers ────────────────────────────────────────────────────────────────────────────────

"""
    arrow_list_column(data, offsets) -> Arrow.List

A list column whose row `i` is `data[offsets[i]:offsets[i+1]-1]` (1-based `offsets`, length n + 1), built directly in
Arrow's layout: Arrow.write writes an Arrow.List as it is, while any other vector of vectors goes through Arrow's
ToList, which collects every row's element first (about 75 bytes per row).
"""
function arrow_list_column(data::Vector{T}, offsets::Vector{Int}) where {T}
    n = length(offsets) - 1
    total = offsets[end] - 1
    O = total <= typemax(Int32) ? Int32 : Int64
    zero_based = O[o - 1 for o in offsets]
    values = Arrow.Primitive(T, UInt8[], Arrow.ValidityBitmap(UInt8[], 1, total, 0), data, total, nothing)
    return Arrow.List{Vector{T}, O, typeof(values)}(UInt8[], Arrow.ValidityBitmap(UInt8[], 1, n, 0),
                                                   Arrow.Offsets(UInt8[], zero_based), values, n, nothing)
end

"The rows of a vector of vectors (e.g. a multi-batch Arrow list column) as one flat Arrow.List."
function flat_list_column(col::AbstractVector)
    T = eltype(eltype(col))
    offsets = Vector{Int}(undef, length(col) + 1); offsets[1] = 1
    for (i, x) in enumerate(col)
        offsets[i + 1] = offsets[i] + length(x)
    end
    data = Vector{T}(undef, offsets[end] - 1)
    for (i, x) in enumerate(col)
        copyto!(data, offsets[i], x, 1, length(x))
    end
    return arrow_list_column(data, offsets)
end

# ── lazy schema-1 columns over schema-2 data ──────────────────────────────────────────────────────────────────────

"Sequences decoded from a packed column."
struct PackedSequenceColumn{C} <: AbstractVector{String}
    packed::C
end
Base.size(v::PackedSequenceColumn) = (length(v.packed),)
Base.@propagate_inbounds Base.getindex(v::PackedSequenceColumn, i::Int) = unpack_sequence(v.packed[i])

"`structural_mods` strings decoded from mod entries and packed sequences."
struct ModStringColumn{E, C} <: AbstractVector{Union{Missing, String}}
    entries::E
    packed::C
    names::Vector{String}
end
Base.size(v::ModStringColumn) = (length(v.entries),)
Base.@propagate_inbounds Base.getindex(v::ModStringColumn, i::Int) = decode_mods(v.entries[i], v.packed[i], v.names)

"values[keys[i]]: a per-precursor view of a per-key table."
struct KeyedColumn{T, K, V} <: AbstractVector{T}
    keys::K
    values::V
end
KeyedColumn(keys::K, values::V) where {K, V} = KeyedColumn{eltype(values), K, V}(keys, values)
Base.size(v::KeyedColumn) = (length(v.keys),)
Base.@propagate_inbounds Base.getindex(v::KeyedColumn, i::Int) = v.values[v.keys[i]]

"""
    PrecursorSideTables

The schema-2 side tables of a library, memory-mapped from its directory.
"""
struct PrecursorSideTables{A, P, S, M, N}
    pep_accession_set::A        # per base_pep_id: accession-set ID (row of accession_sets)
    pep_proteome::P             # per base_pep_id: proteome ID (row of proteomes)
    pep_starts::S               # per base_pep_id: protein start positions
    accession_sets::N           # distinct accession_numbers strings, in sort order
    set_members::M              # per set: its accession IDs (rows of accessions), sorted
    accessions::Vector{String}
    proteomes::Vector{String}
    mod_names::Vector{String}
end

function load_precursor_side_tables(lib_dir::AbstractString)
    f = PRECURSOR_SIDE_FILES
    pep = Arrow.Table(joinpath(lib_dir, f.peptides))
    sets = Arrow.Table(joinpath(lib_dir, f.accession_sets))
    return PrecursorSideTables(pep.accession_set_id, pep.proteome_id, pep.start_idx, sets.accession_numbers, sets.members,
                               String.(Arrow.Table(joinpath(lib_dir, f.accessions)).accession),
                               String.(Arrow.Table(joinpath(lib_dir, f.proteomes)).proteome_identifiers),
                               String.(Arrow.Table(joinpath(lib_dir, f.mod_names)).name))
end

# ── schema 1 -> schema 2 ──────────────────────────────────────────────────────────────────────────────────────────

"""
    convert_precursor_table_v2(lib_dir)

Rewrite a library's schema-1 precursors_table.arrow as schema 2 plus its side tables (in place). Errors, leaving the
library unchanged, if the table breaks an assumption of schema 2 (one set of protein columns per base_pep_id, mod
residues matching the sequence).
"""
function convert_precursor_table_v2(lib_dir::AbstractString)
    path = joinpath(lib_dir, "precursors_table.arrow")
    t = Arrow.Table(path)
    precursor_schema(t) == 1 || error("$path is already schema $(precursor_schema(t))")
    n = length(t.mz)
    seqs = t.sequence; mods = t.structural_mods; bpid = t.base_pep_id
    accs = t.accession_numbers; proteomes = t.proteome_identifiers; starts = t.start_idx
    # per base_pep_id protein columns (must agree on every row of the id)
    n_pep = n == 0 ? 0 : Int(maximum(bpid))
    pep_row = zeros(Int, n_pep)
    for i in 1:n
        j = pep_row[bpid[i]]
        if j == 0
            pep_row[bpid[i]] = i
        elseif !(accs[j] == accs[i] && proteomes[j] == proteomes[i] && collect(starts[j]) == collect(starts[i]))
            error("rows $j and $i share base_pep_id $(bpid[i]) but not their protein columns")
        end
    end
    set_names = sort!(unique(String[accs[pep_row[p]] for p in 1:n_pep if pep_row[p] != 0]))
    set_id = Dict(s => UInt32(k) for (k, s) in enumerate(set_names))
    acc_names = sort!(unique(String[a for s in set_names for a in split(s, ';')]))
    acc_id = Dict(a => UInt32(k) for (k, a) in enumerate(acc_names))
    proteome_names = sort!(unique(String[proteomes[pep_row[p]] for p in 1:n_pep if pep_row[p] != 0]))
    length(proteome_names) <= typemax(UInt16) || error("more than $(typemax(UInt16)) distinct proteome strings")
    proteome_id = Dict(s => UInt16(k) for (k, s) in enumerate(proteome_names))
    pep_set = zeros(UInt32, n_pep); pep_prot = zeros(UInt16, n_pep); pep_starts = [UInt32[] for _ in 1:n_pep]
    for p in 1:n_pep
        j = pep_row[p]; j == 0 && continue
        pep_set[p] = set_id[accs[j]]; pep_prot[p] = proteome_id[proteomes[j]]; pep_starts[p] = collect(UInt32, starts[j])
    end
    # per-row compact columns, in schema-1 column order
    name_id = Dict{String, UInt8}()
    packed = [pack_sequence(seqs[i]) for i in 1:n]
    entries = Union{Missing, Vector{UInt16}}[encode_mods(mods[i], seqs[i], name_id) for i in 1:n]
    dropped = (:accession_numbers, :proteome_identifiers, :start_idx)
    names = Symbol[]; cols = AbstractVector[]
    for (name, col) in zip(Tables.columnnames(t), Tables.columns(t))
        name in dropped && continue
        push!(names, name === :sequence ? :sequence_packed : name === :structural_mods ? :mod_entries : name)
        push!(cols, name === :sequence ? packed : name === :structural_mods ? entries : col)
    end
    mod_names = Vector{String}(undef, length(name_id))
    for (name, id) in name_id; mod_names[id] = name; end
    tmp = path * ".v2tmp"
    Arrow.write(tmp, NamedTuple{Tuple(names)}(Tuple(cols)); metadata = [PRECURSOR_SCHEMA_KEY => "2"])
    f = PRECURSOR_SIDE_FILES
    Arrow.write(joinpath(lib_dir, f.peptides), (accession_set_id = pep_set, proteome_id = pep_prot, start_idx = pep_starts))
    Arrow.write(joinpath(lib_dir, f.accession_sets), (accession_numbers = set_names,
        members = [sort!(unique!(UInt32[acc_id[a] for a in split(s, ';')])) for s in set_names]))
    Arrow.write(joinpath(lib_dir, f.accessions), (accession = acc_names,))
    Arrow.write(joinpath(lib_dir, f.proteomes), (proteome_identifiers = proteome_names,))
    Arrow.write(joinpath(lib_dir, f.mod_names), (name = mod_names,))
    t = nothing; seqs = nothing; mods = nothing; accs = nothing; proteomes = nothing; starts = nothing; bpid = nothing
    GC.gc()
    mv(tmp, path; force = true)
    return path
end
