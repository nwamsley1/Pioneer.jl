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

#==========================================================
Streaming precursor table (BuildSpecLib stages A-D)

Builds the sorted precursor table `precursors.arrow` (what prepare_chronologer_input + predict_retention_times +
parse_chronologer_output produce) with memory bounded by compact per-peptide / per-precursor numeric arrays:
peptide sequences are packed integers (`SeqCode`), precursors are (unit, charge) pairs, and strings exist only for
one output chunk at a time. The content is identical to the in-memory path, including its row order and its
order-dependent id columns (base_pep_id, base_target_id, pair_id, entrapment_pair_id).

Units: a unit is one precursor without its charge.
  target unit  f in 1:F          = (peptide k, mod variant v), in add_mods order (base_pep_id = base_target_id = f)
  entrap unit  F+g, g in 1:G     = (target unit, entrapment sequence e), in add_entrapment_sequences_grouped order
  decoy unit   F+G+u             = the decoy of input unit u (when its group got a decoy)
Rows are (unit, charge); row id = (unit - 1) * n_charges + charge index.
==========================================================#

# ── packed sequences ──────────────────────────────────────────────────────────────────────────────────────────

const SEQ_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"             # codes 1..20 in alphabetical (= String isless) order
const SEQ_MAX_LENGTH = 50                               # 5 bits per residue in 2 x UInt128 (25 + 25 residues)
const _SEQ_CODE = let t = zeros(UInt8, 128); for (i, c) in enumerate(SEQ_ALPHABET); t[Int(c)] = UInt8(i); end; t end
const _CODE_I, _CODE_L = _SEQ_CODE[Int('I')], _SEQ_CODE[Int('L')]

"A peptide sequence packed 5 bits per residue, first residue in the high bits: compares like the String."
struct SeqCode
    hi::UInt128
    lo::UInt128
end
Base.isless(a::SeqCode, b::SeqCode) = a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo)

function encode_seq(s::AbstractString)
    hi = zero(UInt128); lo = zero(UInt128)
    n = 0
    for c in s
        code = UInt128(_SEQ_CODE[Int(c)])
        n += 1
        if n <= 25
            hi |= code << (5 * (25 - n) + 3)
        else
            lo |= code << (5 * (50 - n) + 3)
        end
    end
    return SeqCode(hi, lo)
end

function decode_seq(c::SeqCode)
    io = IOBuffer()
    for n in 1:SEQ_MAX_LENGTH
        code = n <= 25 ? (c.hi >> (5 * (25 - n) + 3)) & 0x1f : (c.lo >> (5 * (50 - n) + 3)) & 0x1f
        code == 0 && break
        write(io, SEQ_ALPHABET[Int(code)])
    end
    return String(take!(io))
end

"I/L-folded code (I -> L), the uniqueness key of combine_shared_peptides and PeptideSequenceSet."
function fold_seq(c::SeqCode)
    hi, lo = c.hi, c.lo
    for n in 1:SEQ_MAX_LENGTH
        if n <= 25
            sh = 5 * (25 - n) + 3
            ((hi >> sh) & 0x1f) == _CODE_I && (hi = (hi & ~(UInt128(0x1f) << sh)) | (UInt128(_CODE_L) << sh))
        else
            sh = 5 * (50 - n) + 3
            ((lo >> sh) & 0x1f) == _CODE_I && (lo = (lo & ~(UInt128(0x1f) << sh)) | (UInt128(_CODE_L) << sh))
        end
    end
    return SeqCode(hi, lo)
end

# ── stage A: digest + combine_shared_peptides ─────────────────────────────────────────────────────────────────

"Unique peptides in combine_shared_peptides output order (first occurrence), with their occurrences."
struct StreamPeptides
    code::Vector{SeqCode}          # first-seen spelling
    nte::Vector{UInt8}             # max num_enzymatic_termini over occurrences
    occ_offsets::Vector{Int}       # peptide k's occurrences: occ_*[occ_offsets[k]:occ_offsets[k+1]-1], encounter order
    occ_protein::Vector{UInt32}    # index into proteins
    occ_start::Vector{UInt32}
    protein_accession::Vector{String}
    protein_proteome::Vector{String}
end

function stream_digest(proteins::Vector{FastaEntry}, regex, max_length::Int, min_length::Int, missed_cleavages::Int,
                       specificity::AbstractString, nterm_met_excision::Bool)
    max_length <= SEQ_MAX_LENGTH || error("streaming build supports peptides up to $SEQ_MAX_LENGTH residues")
    codes = SeqCode[]; prot = UInt32[]; starts = UInt32[]; nte = UInt8[]
    for (p, e) in enumerate(proteins)
        peptides, st, termini = digest_sequence(get_sequence(e), regex, max_length, min_length, missed_cleavages,
                                                specificity; nterm_met_excision = nterm_met_excision)
        for (pep, s0, t) in zip(peptides, st, termini)
            all(aa -> aa ∈ VALID_AAS, pep) || continue
            push!(codes, encode_seq(pep)); push!(prot, UInt32(p)); push!(starts, UInt32(s0)); push!(nte, t)
        end
    end
    n = length(codes)
    folded = map(fold_seq, codes)
    perm = sortperm(folded; alg = Base.Sort.DEFAULT_STABLE)     # runs of one folded sequence, encounter order kept
    run_starts = Int[]
    for i in 1:n
        (i == 1 || folded[perm[i]] != folded[perm[i - 1]]) && push!(run_starts, i)
    end
    push!(run_starts, n + 1)
    folded = nothing
    n_runs = length(run_starts) - 1
    first_occ = [perm[run_starts[r]] for r in 1:n_runs]          # encounter index of each run's first occurrence
    order = sortperm(first_occ)                                   # combine_shared_peptides output order
    pep_code = Vector{SeqCode}(undef, n_runs); pep_nte = Vector{UInt8}(undef, n_runs)
    occ_offsets = Vector{Int}(undef, n_runs + 1); occ_offsets[1] = 1
    occ_protein = Vector{UInt32}(undef, n); occ_start = Vector{UInt32}(undef, n)
    o = 0
    for (k, r) in enumerate(order)
        rng = run_starts[r]:(run_starts[r + 1] - 1)
        pep_code[k] = codes[perm[first(rng)]]
        pep_nte[k] = maximum(i -> nte[perm[i]], rng)
        for i in rng
            o += 1; occ_protein[o] = prot[perm[i]]; occ_start[o] = starts[perm[i]]
        end
        occ_offsets[k + 1] = o + 1
    end
    return StreamPeptides(pep_code, pep_nte, occ_offsets, occ_protein, occ_start,
        String[get_id(e) for e in proteins], String[get_proteome(e) for e in proteins])
end

# ── modification variants (add_mods semantics) ────────────────────────────────────────────────────────────────

struct ModConfig
    fixed::Vector{@NamedTuple{p::Regex, r::String}}
    var::Vector{@NamedTuple{p::Regex, r::String}}
    max_var_mods::Int
    min_var_mods::Int
end

"The mod vectors of `seq`'s variants, in add_mods order (fixed mods then var mods, variants sorted)."
function mod_variants(seq::String, mc::ModConfig)
    fixed = PeptideMod[]
    for m in mc.fixed
        getFixedMods!(fixed, eachmatch(m.p, seq), m.r)
    end
    matches = matchVarMods(seq, mc.var)
    n = countVarModCombinations(matches, mc.max_var_mods, mc.min_var_mods)
    variants = Vector{Vector{PeptideMod}}(undef, n)
    fillVarModStrings!(variants, matches, fixed, mc.max_var_mods, mc.min_var_mods)
    return variants, length(fixed)
end

"""
Derived (entrapment or decoy) sequences: code, the shuffle's position map (`new_positions`, as
adjust_mod_positions takes it) and the group id (entrapment index within its target, 0 for decoys).
"""
struct DerivedSeqs
    code::Vector{SeqCode}
    pos_offsets::Vector{Int}
    positions::Vector{UInt8}
end
DerivedSeqs() = DerivedSeqs(SeqCode[], [1], UInt8[])
function push_derived!(d::DerivedSeqs, code::SeqCode, positions::AbstractVector{UInt8})
    push!(d.code, code); append!(d.positions, positions); push!(d.pos_offsets, length(d.positions) + 1)
    return length(d.code)
end
derived_positions(d::DerivedSeqs, i::Int) = view(d.positions, d.pos_offsets[i]:(d.pos_offsets[i + 1] - 1))

# ── the whole unit structure ──────────────────────────────────────────────────────────────────────────────────

struct StreamUnits
    peps::StreamPeptides
    mc::ModConfig
    var_offsets::Vector{Int}       # target units of peptide k: var_offsets[k]:var_offsets[k+1]-1 (f = unit id)
    unit_pep::Vector{UInt32}       # target unit f -> peptide k
    n_fixed::Vector{UInt8}         # per peptide: number of fixed mods
    entrap::DerivedSeqs            # entrapment sequences e
    entrap_pep::Vector{UInt32}     # e -> target peptide
    entrap_group::Vector{UInt8}    # e -> entrapment_group_id (index among the peptide's successful entrapments)
    eunit_target::Vector{UInt32}   # entrap unit g -> target unit f
    eunit_seq::Vector{UInt32}      # entrap unit g -> entrapment sequence e
    decoy::DerivedSeqs             # decoy sequences, one per decoyed group
    decoy_of_input::Vector{UInt32} # input unit u (1:F+G) -> decoy sequence index (0 = no decoy)
    F::Int
    G::Int
end

n_units(s::StreamUnits) = 2 * (s.F + s.G)
is_decoy_unit(s::StreamUnits, u::Int) = u > s.F + s.G
input_of(s::StreamUnits, u::Int) = is_decoy_unit(s, u) ? u - s.F - s.G : u
has_unit(s::StreamUnits, u::Int) = !is_decoy_unit(s, u) || s.decoy_of_input[u - s.F - s.G] != 0

function unit_code(s::StreamUnits, u::Int)
    if is_decoy_unit(s, u)
        return s.decoy.code[s.decoy_of_input[u - s.F - s.G]]
    elseif u > s.F
        return s.entrap.code[s.eunit_seq[u - s.F]]
    else
        return s.peps.code[s.unit_pep[u]]
    end
end
unit_target(s::StreamUnits, u::Int) = (i = input_of(s, u); i > s.F ? Int(s.eunit_target[i - s.F]) : i)
unit_peptide(s::StreamUnits, u::Int) = Int(s.unit_pep[unit_target(s, u)])
base_target_id(s::StreamUnits, u::Int) = UInt32(input_of(s, u))
base_pep_id(s::StreamUnits, u::Int) = UInt32(unit_target(s, u))
entrapment_group(s::StreamUnits, u::Int) = (i = input_of(s, u); i > s.F ? s.entrap_group[s.eunit_seq[i - s.F]] : 0x00)

"A unit's mods: target = add_mods order; entrapment / decoy = adjust_mod_positions (sorted)."
function unit_mods(s::StreamUnits, u::Int, cache::Dict{Int, Tuple{Vector{Vector{PeptideMod}}, Int}})
    f = unit_target(s, u); k = Int(s.unit_pep[f])
    variants, _ = get!(() -> mod_variants(decode_seq(s.peps.code[k]), s.mc), cache, k)
    mods = variants[f - s.var_offsets[k] + 1]
    i = input_of(s, u)
    L = UInt8(length(decode_seq(s.peps.code[k])))
    if i > s.F
        mods = adjust_mod_positions(mods, collect(derived_positions(s.entrap, Int(s.eunit_seq[i - s.F]))), L)
    end
    if is_decoy_unit(s, u)
        mods = adjust_mod_positions(mods, collect(derived_positions(s.decoy, Int(s.decoy_of_input[i]))), L)
    end
    return mods, length(variants[f - s.var_offsets[k] + 1]) - s.n_fixed[k]
end

function new_shuffler()
    ShuffleSeq("", Vector{Char}(undef, 255), Vector{UInt8}(undef, 255), Vector{UInt8}(undef, 255),
               zero(UInt8), zero(UInt8), Vector{Char}())
end

"First candidate with `method`, then up to `max_attempts` shuffles, avoiding `taken` (I/L folded). nothing = exhausted."
function _unique_shuffle!(ss, seq::String, method::String, rng, taken::Function, max_attempts::Int)
    cand = shuffle_sequence!(ss, seq; method = method, rng = rng)
    n = 0
    while taken(fold_seq(encode_seq(cand))) && n < max_attempts
        cand = shuffle_sequence!(ss, seq; method = "shuffle", rng = rng)
        n += 1
    end
    taken(fold_seq(encode_seq(cand))) && return nothing
    return cand
end

function stream_units(peps::StreamPeptides, mc::ModConfig, entrapment_r::Int, entrapment_method::String,
                      add_decoys::Bool, decoy_method::String, seed::Int)
    U = length(peps.code)
    # add_mods: variant counts per peptide (a peptide with no admissible variant emits no units)
    var_offsets = Vector{Int}(undef, U + 1); var_offsets[1] = 1
    n_fixed = Vector{UInt8}(undef, U)
    for k in 1:U
        variants, nf = mod_variants(decode_seq(peps.code[k]), mc)
        n_fixed[k] = UInt8(nf)
        var_offsets[k + 1] = var_offsets[k] + length(variants)
    end
    F = var_offsets[end] - 1
    unit_pep = Vector{UInt32}(undef, F)
    for k in 1:U, f in var_offsets[k]:(var_offsets[k + 1] - 1); unit_pep[f] = UInt32(k); end
    # only peptides with units take part in grouping (add_mods emitted nothing for the others)
    live = [k for k in 1:U if var_offsets[k + 1] > var_offsets[k]]
    target_folded = sort!([fold_seq(peps.code[k]) for k in live])
    in_targets(c) = (r = searchsorted(target_folded, c); !isempty(r))

    # entrapments: groups = exact target sequences, sorted; each slot shuffled from the peptide's own RNG stream
    entrap = DerivedSeqs(); entrap_pep = UInt32[]; entrap_group = UInt8[]
    reserved = Set{SeqCode}()
    sorted_live = sort(live, by = k -> peps.code[k])
    ss = new_shuffler()
    pep_entraps = [UInt32[] for _ in 1:(entrapment_r > 0 ? U : 0)]
    if entrapment_r > 0
        taken(c) = in_targets(c) || c in reserved
        for k in sorted_live
            seq = decode_seq(peps.code[k])
            rng = sequence_rng(seq, seed, 2)
            for _ in 1:entrapment_r
                cand = _unique_shuffle!(ss, seq, entrapment_method, rng, taken, 20)
                cand === nothing && continue
                e = push_derived!(entrap, encode_seq(cand), view(ss.new_positions, 1:length(seq)))
                push!(entrap_pep, UInt32(k)); push!(entrap_group, UInt8(length(pep_entraps[k]) + 1))
                push!(pep_entraps[k], UInt32(e))
                push!(reserved, fold_seq(entrap.code[e]))
            end
        end
    end
    # entrapment units: by sorted base sequence, then target unit, then entrapment index
    eunit_target = UInt32[]; eunit_seq = UInt32[]
    for k in sorted_live, f in var_offsets[k]:(var_offsets[k + 1] - 1), e in (entrapment_r > 0 ? pep_entraps[k] : UInt32[])
        push!(eunit_target, UInt32(f)); push!(eunit_seq, e)
    end
    G = length(eunit_target)

    # decoys: groups = exact sequences of targets and entrapments, in sorted order; one decoy sequence per group
    decoy = DerivedSeqs(); decoy_of_input = zeros(UInt32, F + G)
    if add_decoys && decoy_method != "diann_mutation"
        entrap_units_of = [UInt32[] for _ in 1:length(entrap.code)]
        for g in 1:G; push!(entrap_units_of[eunit_seq[g]], UInt32(F + g)); end
        groups = Tuple{SeqCode, Int, Int}[]        # (code, kind 0 = target peptide / 1 = entrapment, index)
        for k in live; push!(groups, (peps.code[k], 0, k)); end
        for e in eachindex(entrap.code); push!(groups, (entrap.code[e], 1, e)); end
        sort!(groups, by = first)
        entrap_folded = Set{SeqCode}(fold_seq(c) for c in entrap.code)
        taken_d(c) = in_targets(c) || c in entrap_folded || c in reserved
        empty!(reserved)
        for (code, kind, idx) in groups
            seq = decode_seq(code)
            rng = sequence_rng(seq, seed, 1)
            cand = _unique_shuffle!(ss, seq, decoy_method, rng, taken_d, 20)
            cand === nothing && continue
            d = push_derived!(decoy, encode_seq(cand), view(ss.new_positions, 1:length(seq)))
            push!(reserved, fold_seq(decoy.code[d]))
            members = kind == 0 ? (var_offsets[idx]:(var_offsets[idx + 1] - 1)) : entrap_units_of[idx]
            for u in members; decoy_of_input[u] = UInt32(d); end
        end
    end
    return StreamUnits(peps, mc, var_offsets, unit_pep, n_fixed, entrap, entrap_pep, entrap_group,
                       eunit_target, eunit_seq, decoy, decoy_of_input, F, G)
end

# ── stage B/C/D: rows, m/z filter, ids, retention times, final order, write ───────────────────────────────────

"CSR (offsets, items) of `owner[i]` -> the items i it owns, items in ascending order; owner 0 = none."
function owner_csr(owner::AbstractVector{<:Integer}, n_owners::Int)
    offsets = zeros(Int, n_owners + 2)
    for o in owner; o != 0 && (offsets[o + 2] += 1); end
    offsets[1] = 1; offsets[2] = 1
    for j in 3:(n_owners + 2); offsets[j] += offsets[j - 1]; end
    items = Vector{UInt32}(undef, offsets[end] - 1)
    for (i, o) in enumerate(owner)
        o == 0 && continue
        items[offsets[o + 1]] = UInt32(i); offsets[o + 1] += 1
    end
    return view(offsets, 1:(n_owners + 1)), items     # owner o: items[offsets[o]:offsets[o+1]-1]
end

"""
The units in the in-memory path's pre-sort row order (rows = these units x charges): with decoys,
add_decoy_sequences_grouped sorts all rows by sequence (stable: a group's rows keep their input order); without,
targets then entrapments in unit order.
"""
function presort_units(s::StreamUnits, has_decoy_step::Bool)
    has_decoy_step || return UInt32.(1:(s.F + s.G))
    e_off, e_items = owner_csr(s.eunit_seq, length(s.entrap.code))      # entrapment seq -> entrap units g
    d_off, d_items = owner_csr(s.decoy_of_input, length(s.decoy.code))  # decoy seq -> input units
    # groups as (kind << 40 | index): 0 = target peptide, 1 = entrapment sequence, 2 = decoy sequence
    keys = UInt64[]
    for k in eachindex(s.peps.code); s.var_offsets[k + 1] > s.var_offsets[k] && push!(keys, UInt64(k)); end
    for e in eachindex(s.entrap.code); e_off[e + 1] > e_off[e] && push!(keys, (UInt64(1) << 40) | UInt64(e)); end
    for d in eachindex(s.decoy.code); d_off[d + 1] > d_off[d] && push!(keys, (UInt64(2) << 40) | UInt64(d)); end
    group_code(x) = (kind = x >> 40; i = Int(x & 0xffffffffff);
        kind == 0 ? s.peps.code[i] : kind == 1 ? s.entrap.code[i] : s.decoy.code[i])
    sort!(keys, by = group_code)        # target / entrapment / decoy sequences are distinct: no ties
    out = Vector{UInt32}(undef, 0); sizehint!(out, 2 * (s.F + s.G))
    for x in keys
        kind = x >> 40; i = Int(x & 0xffffffffff)
        if kind == 0
            for f in s.var_offsets[i]:(s.var_offsets[i + 1] - 1); push!(out, UInt32(f)); end
        elseif kind == 1
            for j in e_off[i]:(e_off[i + 1] - 1); push!(out, UInt32(s.F + e_items[j])); end
        else
            for j in d_off[i]:(d_off[i + 1] - 1); push!(out, UInt32(s.F + s.G + d_items[j])); end
        end
    end
    return out
end

"Log a phase's time with the process peak RSS so far and the live heap after a full collection."
function _stream_phase(name::String, t0::Float64)
    GC.gc()
    @user_info @sprintf("Streaming build: %-34s %8.1f s  peak RSS %6.2f GB  live heap %6.2f GB",
                        name, time() - t0, peak_rss() / 1e9, Base.gc_live_bytes() / 1e9)
    return time()
end

"Precursor m/z exactly as getMZs computes it (Float64 residue sum, then mod sum in mods-string order)."
function unit_mz(seq::String, mods_string::String, charge::UInt8, mod_to_mass::Dict{String, Float64})
    return Float32(getMZ(seq, mods_string, charge, mod_to_mass))
end

"""
    build_precursor_table_streaming(params, prec_mz_min, prec_mz_max, out_path, proteins_out_path;
                                    chunk_rows = 2_000_000) -> out_path

Streaming equivalent of prepare_chronologer_input + predict_retention_times + parse_chronologer_output: writes the
sorted precursor table to `out_path` (record batches of `chunk_rows`) and the protein table to `proteins_out_path`.
Retention times come from the active Koina client, once per unit (they depend on the Koina sequence only).
"""
function build_precursor_table_streaming(params::Dict{String, Any}, prec_mz_min::Float32, prec_mz_max::Float32,
                                         out_path::String, proteins_out_path::String;
                                         chunk_rows::Int = 1_000_000, rt_bin_tol::Float32 = 3.0f0)
    t = time()
    dp = params["fasta_digest_params"]; lp = params["library_params"]
    isempty(String(get(lp, "im_model", ""))) || error("streaming build does not support im_model yet")
    min_len, max_len = clamp_digest_length_to_model(get(lp, "prediction_model", "altimeter"),
                                                    dp["min_length"], dp["max_length"])
    # modification config, as prepare_chronologer_input reads it (mass strings: Float32 text, parsed as Float64)
    var_mods = @NamedTuple{p::Regex, r::String}[]; fixed_mods = @NamedTuple{p::Regex, r::String}[]
    mod_to_mass_dict = Dict{String, String}()
    for i in eachindex(params["variable_mods"]["pattern"])
        push!(var_mods, (p = Regex(params["variable_mods"]["pattern"][i]), r = params["variable_mods"]["name"][i]))
        mod_to_mass_dict[params["variable_mods"]["name"][i]] = string(Float32(params["variable_mods"]["mass"][i]))
    end
    for i in eachindex(params["fixed_mods"]["pattern"])
        push!(fixed_mods, (p = Regex(params["fixed_mods"]["pattern"][i]), r = params["fixed_mods"]["name"][i]))
        mod_to_mass_dict[params["fixed_mods"]["name"][i]] = string(Float32(params["fixed_mods"]["mass"][i]))
    end
    mod_to_mass = Dict(k => parse(Float64, v) for (k, v) in mod_to_mass_dict)
    mc = ModConfig(fixed_mods, var_mods, dp["max_var_mods"], Int(get(dp, "min_var_mods", 0)))

    # proteins, in FASTA order
    rx(key) = haskey(params, key) ? [isempty(r) ? nothing : Regex(r) for r in params[key]] : fill(nothing, length(params["fasta_paths"]))
    proteins = FastaEntry[]
    for (name, fasta, a, g, p, o) in zip(params["fasta_names"], params["fasta_paths"], rx("fasta_header_regex_accessions"),
                                         rx("fasta_header_regex_genes"), rx("fasta_header_regex_proteins"), rx("fasta_header_regex_organisms"))
        append!(proteins, parse_fasta(fasta, name; accession_regex = a, gene_regex = g, protein_regex = p, organism_regex = o))
    end
    Arrow.write(proteins_out_path, build_protein_df(proteins))
    t = _stream_phase("proteins", t)

    regex = isnothing(dp["cleavage_regex"]) ? nothing : Regex(dp["cleavage_regex"])
    peps = stream_digest(proteins, regex, max_len, min_len, dp["missed_cleavages"], get(dp, "specificity", "full"),
                         get(dp, "nterm_met_excision", true))
    @user_info "Streaming build: $(length(proteins)) proteins, $(length(peps.code)) unique peptides"
    t = _stream_phase("digest", t)
    seed = Int(get(params, "seed", 1844))
    units = stream_units(peps, mc, Int(dp["entrapment_r"]), get(dp, "entrapment_method", "shuffle"),
                         dp["add_decoys"], get(dp, "decoy_method", "shuffle"), seed)
    proteins = nothing
    charges = UInt8(dp["min_charge"]):UInt8(dp["max_charge"])
    nz = length(charges)
    NU = n_units(units)
    @user_info "Streaming build: $(units.F) target, $(units.G) entrapment, $(count(!=(0), units.decoy_of_input)) decoy units"
    t = _stream_phase("mods, entrapments, decoys", t)

    # per unit: sequence string helpers, mods strings and the Koina sequence are recomputed on demand
    row_mz = zeros(Float32, NU * nz); present = falses(NU * nz)
    unit_irt = zeros(Float32, NU)
    cache = Dict{Int, Tuple{Vector{Vector{PeptideMod}}, Int}}()
    rt_model = String(lp["rt_model"])
    batch_units = Int[]; batch_seqs = String[]
    flush_rt!() = begin
        isempty(batch_units) && return
        rts = predict_rt_koina(DataFrame(koina_sequence = batch_seqs); rt_model = rt_model)
        for (u, rt) in zip(batch_units, rts); unit_irt[u] = rt; end
        empty!(batch_units); empty!(batch_seqs)
    end
    for u in 1:NU
        has_unit(units, u) || continue
        length(cache) > 100_000 && empty!(cache)
        seq = decode_seq(unit_code(units, u))
        mods, _ = unit_mods(units, u, cache)
        ms = getModString(mods)
        any_row = false
        for (zi, z) in enumerate(charges)
            mz = unit_mz(seq, ms, z, mod_to_mass)
            r = (u - 1) * nz + zi
            row_mz[r] = mz
            if mz >= prec_mz_min && mz <= prec_mz_max
                present[r] = true; any_row = true
            end
        end
        if any_row
            push!(batch_units, u); push!(batch_seqs, koina_sequence(seq, ms))
            length(batch_units) >= 1_000_000 && flush_rt!()
        end
    end
    flush_rt!()
    t = _stream_phase("m/z filter + retention times", t)

    # pre-sort order (sequence, then input order, then charge) -> pair_id and entrapment_pair_id
    presort = presort_units(units, dp["add_decoys"] && get(dp, "decoy_method", "shuffle") != "diann_mutation")
    pair_id = zeros(UInt32, NU * nz); epair_id = zeros(UInt32, NU * nz)
    next_pair = UInt32(1)
    row_of(u, zi) = (u - 1) * nz + zi
    for u in Iterators.map(Int, presort), zi in 1:nz
        r = row_of(u, zi)
        (present[r] && !is_decoy_unit(units, u)) || continue
        pair_id[r] = next_pair
        d = units.decoy_of_input[u]
        if d != 0
            dr = row_of(units.F + units.G + u, zi)
            present[dr] && (pair_id[dr] = next_pair)
        end
        next_pair += 1
    end
    for u in Iterators.map(Int, presort), zi in 1:nz
        r = row_of(u, zi)
        (present[r] && is_decoy_unit(units, u) && pair_id[r] == 0) || continue
        pair_id[r] = next_pair; next_pair += 1
    end
    next_e = UInt32(1)
    for u in Iterators.map(Int, presort), zi in 1:nz
        r = row_of(u, zi)
        (present[r] && !is_decoy_unit(units, u)) || continue
        tr = row_of(unit_target(units, u), zi)          # the entrapment's group-0 target row (itself for targets)
        present[tr] || continue
        if epair_id[tr] == 0
            epair_id[tr] = next_e; next_e += 1
        end
        epair_id[r] = epair_id[tr]
    end
    presort = nothing
    t = _stream_phase("pair ids", t)

    # final order: iRT then (sequence, mods, charge, decoy, entrapment group); m/z within 3-iRT blocks, same ties
    rows = UInt32[r for r in 1:(NU * nz) if present[r]]
    unit_of(r) = (Int(r) - 1) ÷ nz + 1
    zi_of(r) = (Int(r) - 1) % nz + 1
    mods_str_cache = Dict{Int, String}()
    mods_string(u) = get!(() -> getModString(first(unit_mods(units, u, cache))), mods_str_cache, u)
    function tie_less(a, b)       # (sequence, mods, charge, decoy, entrapment group)
        ua, ub = unit_of(a), unit_of(b)
        ca, cb = unit_code(units, ua), unit_code(units, ub)
        ca != cb && return isless(ca, cb)
        if ua != ub
            ma, mb = mods_string(ua), mods_string(ub)
            ma != mb && return isless(ma, mb)
        end
        zi_of(a) != zi_of(b) && return zi_of(a) < zi_of(b)
        da, db = is_decoy_unit(units, ua), is_decoy_unit(units, ub)
        da != db && return isless(da, db)
        return entrapment_group(units, ua) < entrapment_group(units, ub)
    end
    irt_lt(a, b) = (ia = unit_irt[unit_of(a)]; ib = unit_irt[unit_of(b)]; isequal(ia, ib) ? tie_less(a, b) : isless(ia, ib))
    sort!(rows; lt = irt_lt)
    function mz_lt(a, b)
        ma, mb = row_mz[a], row_mz[b]
        isequal(ma, mb) || return isless(ma, mb)
        return irt_lt(a, b)
    end
    start = 1; start_irt = unit_irt[unit_of(rows[1])]
    for i in 1:length(rows)
        irt_i = unit_irt[unit_of(rows[i])]
        if (irt_i - start_irt) > rt_bin_tol && i > start
            sort!(view(rows, start:(i - 1)); lt = mz_lt); start = i; start_irt = irt_i
        end
    end
    sort!(view(rows, start:length(rows)); lt = mz_lt)
    empty!(mods_str_cache)
    t = _stream_phase("final sort", t)

    # write in chunks
    nce = Float32(params["nce_params"]["nce"])
    cleavage = isnothing(dp["cleavage_regex"]) ? nothing : Regex(dp["cleavage_regex"])
    # One Arrow.write call per chunk: each call returns once its record batch is written, so one chunk's strings are
    # alive at a time (Arrow.write over a lazy partitioner processes later partitions in @async tasks, letting the
    # producer run ahead and every chunk accumulate).
    function chunk_columns(chunk)
        n = length(chunk)
        cols = (proteome_identifiers = Vector{String}(undef, n), accession_number = Vector{String}(undef, n),
            sequence = Vector{String}(undef, n), start_idx = Vector{Vector{UInt32}}(undef, n),
            mods = Vector{Union{Missing, String}}(undef, n), isotopic_mods = Vector{Union{Missing, String}}(missing, n),
            num_variable_modifications = Vector{UInt8}(undef, n), precursor_charge = Vector{UInt8}(undef, n),
            num_enzymatic_termini = Vector{UInt8}(undef, n), collision_energy = fill(nce, n), decoy = Vector{Bool}(undef, n),
            entrapment_group_id = Vector{UInt8}(undef, n), base_target_id = Vector{UInt32}(undef, n),
            base_pep_id = Vector{UInt32}(undef, n), pair_id = Vector{Union{Missing, UInt32}}(undef, n),
            koina_sequence = Vector{String}(undef, n), mz = Vector{Float32}(undef, n), length = Vector{UInt8}(undef, n),
            missed_cleavages = Vector{UInt8}(undef, n), entrapment_pair_id = Vector{Union{Missing, UInt32}}(undef, n),
            irt = Vector{Float32}(undef, n), sulfur_count = Vector{UInt8}(undef, n),
            isotope_mods = Vector{Union{Missing, String}}(missing, n))
        for (j, r) in enumerate(chunk)
            u = unit_of(r); zi = zi_of(r); k = unit_peptide(units, u)
            occ = peps.occ_offsets[k]:(peps.occ_offsets[k + 1] - 1)
            cols.proteome_identifiers[j] = join((peps.protein_proteome[peps.occ_protein[o]] for o in reverse(occ)), ';')
            cols.accession_number[j] = join((peps.protein_accession[peps.occ_protein[o]] for o in reverse(occ)), ';')
            seq = decode_seq(unit_code(units, u))
            mods, nvar = unit_mods(units, u, cache)
            ms = getModString(mods)
            cols.sequence[j] = seq
            cols.start_idx[j] = UInt32[peps.occ_start[o] for o in reverse(occ)]
            cols.mods[j] = ms
            cols.num_variable_modifications[j] = UInt8(nvar)
            cols.precursor_charge[j] = charges[zi]
            cols.num_enzymatic_termini[j] = peps.nte[k]
            cols.decoy[j] = is_decoy_unit(units, u)
            cols.entrapment_group_id[j] = entrapment_group(units, u)
            cols.base_target_id[j] = base_target_id(units, u)
            cols.base_pep_id[j] = base_pep_id(units, u)
            cols.pair_id[j] = pair_id[r]
            cols.koina_sequence[j] = koina_sequence(seq, ms)
            cols.mz[j] = row_mz[r]
            cols.length[j] = UInt8(length(seq))
            cols.missed_cleavages[j] = cleavage === nothing ? 0x00 : UInt8(count(cleavage, seq))
            cols.entrapment_pair_id[j] = cols.decoy[j] || epair_id[r] == 0 ? missing : epair_id[r]
            cols.irt[j] = unit_irt[u]
            cols.sulfur_count[j] = UInt8(count(c -> c == 'C' || c == 'M', seq))
        end
        length(cache) > 100_000 && empty!(cache)
        return cols
    end
    open(Arrow.Writer, out_path) do writer
        for chunk in Iterators.partition(rows, chunk_rows)
            Arrow.write(writer, chunk_columns(chunk))
        end
    end
    t = _stream_phase("write", t)
    return out_path
end
