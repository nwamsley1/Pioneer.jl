# Isomer (localization) decoy generation — move-one-random onto shadow sites.
#
# For each MODIFIED target peptidoform we build ONE decoy peptidoform by moving a
# single randomly-chosen modification to a random decoy "shadow" site of the same
# mod type (keeping the other mods on their real sites). This mirrors a real
# one-site-off localization error. The decoy is NOT predicted by Koina; its
# fragments are derived from the target's by permuting only the moved site's m/z.
# Same sequence/charge/mz (mass preserved) -> co-isolated sibling.
#
# Shadow sites come from shadow_decoy_sites(seq) (site_determining.jl) — the SAME
# deterministic per-sequence draw retention uses, so retention already covers the
# real-vs-decoy gaps. Unmodified peptidoforms get no decoy.
#
# See dev_docs/prosit_ptm_integration/isomer_decoy_shadow_sites.md.

"""
    fragment_contains(base_type, frag_index, L, pos) -> Bool

Does a fragment contain sequence position `pos`? b-ion: `pos <= frag_index`;
y-ion: `pos >= L-frag_index+1`; precursor/other spans all -> true.
"""
@inline function fragment_contains(base_type::Char, frag_index::Integer, L::Integer, pos::Integer)
    base_type == 'b' && return pos <= frag_index
    base_type == 'y' && return pos >= L - frag_index + 1
    return true
end

"""
    permute_fragment_mz(mz, base_type, frag_index, frag_charge, L, k_old, k_new, mod_mass) -> new_mz

m/z when a modification moves from `k_old` to `k_new`: shift only ions containing
exactly one of the two sites by +/- mod_mass/charge.
"""
@inline function permute_fragment_mz(mz::T, base_type::Char, frag_index::Integer,
                                     frag_charge::Integer, L::Integer,
                                     k_old::Integer, k_new::Integer, mod_mass::Real) where {T<:Real}
    has_old = fragment_contains(base_type, frag_index, L, k_old)
    has_new = fragment_contains(base_type, frag_index, L, k_new)
    has_old == has_new && return mz
    delta = T(mod_mass / frag_charge)
    return has_new ? mz + delta : mz - delta
end

# Parse ":mods" string "(pos,aa,name)(pos,aa,name)" <-> (pos, aa, name) tuples.
function _parse_mod_tuples(s::AbstractString)
    out = Tuple{Int, Char, String}[]
    for m in eachmatch(r"\((\d+),([A-Za-z]),([^)]+)\)", s)
        push!(out, (parse(Int, m.captures[1]), m.captures[2][1], String(m.captures[3])))
    end
    return out
end
_format_mod_tuples(mods) = join(["($(p),$(a),$(n))" for (p, a, n) in sort(mods, by = x -> x[1])])

"""
    generate_isomer_decoys!(precursors_df, fragments_df, var_mods, mod_masses; global_seed)
        -> (precursors_df, fragments_df)

Append one isomer decoy per MODIFIED target peptidoform (Prosit path). Groups
precursors by sequence; per peptide draws the deterministic shadow sites; per
modified peptidoform picks a random modified site and swaps it to a random shadow
site of the same mod type, keeping decoy peptidoforms unique within the peptide.
Emits a decoy precursor (mods rewritten, `is_loc_decoy=true`, `loc_base_prec_id`,
fresh `pair_id`, same mz) and its fragments (one-site m/z permutation).
`precursor_idx == row index`; decoys are appended after all reals. Adds
`is_loc_decoy::Bool` + `loc_base_prec_id::UInt32`. No-op without var mods.
Callers MUST use the return values (vcat, not in-place).
"""
function generate_isomer_decoys!(precursors_df::DataFrame, fragments_df::DataFrame,
                                 var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
                                 mod_masses::AbstractDict;
                                 global_seed::UInt64 = _SHADOW_GLOBAL_SEED)
    R = nrow(precursors_df)
    precursors_df[!, :is_loc_decoy]     = falses(R)
    precursors_df[!, :loc_base_prec_id] = zeros(UInt32, R)
    (isempty(var_mods) || R == 0) && return precursors_df, fragments_df

    fr_pid = fragments_df[!, :precursor_idx]
    frag_rows = Dict{UInt32, Vector{Int}}()
    for i in 1:length(fr_pid)
        push!(get!(frag_rows, UInt32(fr_pid[i]), Int[]), i)
    end

    seqs = precursors_df[!, :sequence]
    modcol = precursors_df[!, :mods]
    # charge column: intermediate table uses :precursor_charge (renamed to
    # :prec_charge later). Uniqueness is per (charge, decoy-mods) so the same
    # peptidoform at different charges each gets its own decoy.
    chargecol = hasproperty(precursors_df, :prec_charge) ? precursors_df[!, :prec_charge] :
                hasproperty(precursors_df, :precursor_charge) ? precursors_df[!, :precursor_charge] :
                nothing
    # group rows by sequence (a peptide's peptidoforms decoy together for uniqueness)
    by_seq = Dict{String, Vector{Int}}()
    for r in 1:R; push!(get!(by_seq, String(seqs[r]), Int[]), r); end

    base_rows     = Int[]          # target row each decoy copies
    decoy_modstrs = String[]
    decoy_move    = Tuple{Int,Int,Float64}[]   # (k_old, k_new, mod_mass) per decoy

    for (seq, rows) in by_seq
        real_by = real_sites_by_mod(seq, var_mods)
        isempty(real_by) && continue
        shadow = shadow_decoy_sites(seq, real_by; seed = seq_seed(seq; global_seed))
        used = Set{Tuple{Int,String}}()
        rng  = Random.MersenneTwister(seq_seed(seq; global_seed) ⊻ 0xD5C0FFEE)
        L = length(seq)
        for r in rows
            ms = modcol[r]
            (ismissing(ms) || isempty(ms)) && continue
            mods = _parse_mod_tuples(ms)
            isempty(mods) && continue
            z = chargecol === nothing ? 0 : Int(chargecol[r])
            # candidate decoys: move one mod (pos,name) to a shadow site of that name
            cand_str = String[]; cand_mods = Vector{Tuple{Int,Char,String}}[]
            cand_move = Tuple{Int,Int,String}[]
            for (pos, res, name) in mods
                haskey(shadow, name) || continue
                for d in shadow[name]
                    nm = [(p == pos && n == name) ? (d, seq[d], n) : (p, a, n) for (p, a, n) in mods]
                    str = _format_mod_tuples(nm)
                    (z, str) in used && continue
                    push!(cand_str, str); push!(cand_mods, nm); push!(cand_move, (pos, d, name))
                end
            end
            isempty(cand_str) && continue
            j = rand(rng, 1:length(cand_str))
            push!(used, (z, cand_str[j]))
            push!(base_rows, r)
            push!(decoy_modstrs, cand_str[j])
            k_old, k_new, name = cand_move[j]
            push!(decoy_move, (k_old, k_new, Float64(get(mod_masses, name, 0.0))))
        end
    end

    ndec = length(base_rows)
    ndec == 0 && return precursors_df, fragments_df

    # decoy precursor rows
    decoy_prec = precursors_df[base_rows, :]
    decoy_prec[!, :mods]             = decoy_modstrs
    decoy_prec[!, :is_loc_decoy]     = trues(ndec)
    decoy_prec[!, :loc_base_prec_id] = UInt32.(base_rows)
    if hasproperty(precursors_df, :pair_id)
        mx = UInt32(0)
        @inbounds for x in precursors_df[!, :pair_id]
            ismissing(x) || (mx = max(mx, UInt32(x)))
        end
        decoy_prec[!, :pair_id] = UInt32[mx + UInt32(i) for i in 1:ndec]
    end

    # decoy fragment rows: copy the target's fragments, permute the one moved site.
    src_rows = Int[]; new_mz = Float32[]; new_pid = UInt32[]
    ann_col = fragments_df[!, :annotation]
    mz_col  = fragments_df[!, :mz]
    for di in 1:ndec
        r = base_rows[di]; (k_old, k_new, mass) = decoy_move[di]
        L = length(String(seqs[r]))
        dpid = UInt32(R + di)
        for fi in get(frag_rows, UInt32(r), Int[])
            pa = parse_fragment_annotation(GenericFragAnnotation(String(ann_col[fi])))
            push!(src_rows, fi)
            push!(new_mz, permute_fragment_mz(Float32(mz_col[fi]), pa.base_type,
                          Int(pa.frag_index), Int(pa.charge), L, k_old, k_new, mass))
            push!(new_pid, dpid)
        end
    end
    decoy_frag = fragments_df[src_rows, :]
    decoy_frag[!, :mz]            = new_mz
    decoy_frag[!, :precursor_idx] = new_pid

    return vcat(precursors_df, decoy_prec), vcat(fragments_df, decoy_frag)
end
