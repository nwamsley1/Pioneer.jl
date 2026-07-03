# Isomer (localization) decoy generation.
#
# An isomer decoy is a real peptidoform with one modification moved to a wrong
# (non-acceptor) site. It is NOT predicted by Koina (the acceptor can't carry the
# mod); instead its fragments are DERIVED from the real isomer's predicted
# fragments by permuting only the site-determining ion m/z — site-determining ions
# shift by ±mod_mass/charge, everything else (intensities, RT, precursor m/z)
# unchanged. This makes the hardest possible decoy: identical spectrum except the
# masses of the ions that actually determine the site.
#
# See dev_docs/prosit_ptm_integration/localization_decoy_library_plan.md and
# site_determining_fragment_retention.md. Decoy sites here MUST match the
# retention candidate set (decoy_neighbor_positions), so the reals keep the ions
# that distinguish these decoys.

"""
    fragment_contains(base_type, frag_index, L, pos) -> Bool

Whether a fragment contains sequence position `pos`. b-ion (residues `1..i`):
`pos <= frag_index`. y-ion (residues `L-frag_index+1..L`): `pos >= L-frag_index+1`.
Precursor / other ions span the whole peptide -> always `true`.
"""
@inline function fragment_contains(base_type::Char, frag_index::Integer, L::Integer, pos::Integer)
    base_type == 'b' && return pos <= frag_index
    base_type == 'y' && return pos >= L - frag_index + 1
    return true
end

"""
    permute_fragment_mz(mz, base_type, frag_index, frag_charge, L, k_old, k_new, mod_mass) -> new_mz

m/z of a fragment when a modification is moved from site `k_old` to `k_new`
(keeping intensity, ion type, charge). Only ions containing **exactly one** of the
two sites change: subtract `mod_mass/charge` if the ion contained `k_old` but not
`k_new`, add it if it contains `k_new` but not `k_old`, otherwise unchanged.
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

"""
    choose_decoy_site(L, acc_positions, mod_site, spacing, prefer_plus) -> Int

Decoy position for a modification at `mod_site`: `mod_site ± spacing`, in range
`1..L` and **not** an acceptor (so the placement is wrong by construction). Tries
the preferred sign first (randomized per precursor for +/- balance), then the
other. Returns 0 if neither is valid (e.g. flanked by acceptors or at the
terminus) — the caller then skips the decoy for that peptidoform.
"""
function choose_decoy_site(L::Integer, acc_positions, mod_site::Integer,
                           spacing::Integer, prefer_plus::Bool)
    accset = Set{Int}(Int.(acc_positions))
    cands = prefer_plus ? (mod_site + spacing, mod_site - spacing) :
                          (mod_site - spacing, mod_site + spacing)
    for d in cands
        (1 <= d <= L && !(d in accset)) && return Int(d)
    end
    return 0
end

# Parse a ":mods" string "(pos,aa,name)(pos,aa,name)" into (pos, aa, name) tuples,
# and format the reverse (matches getModString, sorted by position).
function _parse_mod_tuples(s::AbstractString)
    out = Tuple{Int, Char, String}[]
    for m in eachmatch(r"\((\d+),([A-Za-z]),([^)]+)\)", s)
        push!(out, (parse(Int, m.captures[1]), m.captures[2][1], String(m.captures[3])))
    end
    return out
end
_format_mod_tuples(mods) = join(["($(p),$(a),$(n))" for (p, a, n) in sort(mods, by = x -> x[1])])

"""
    generate_isomer_decoys!(precursors_df, fragments_df, var_mods, mod_masses;
                            spacing=1, target_mod_name="Unimod:21") -> (precursors_df, fragments_df)

Append isomer (localization) decoys to the intermediate precursor + raw-fragment
tables (Prosit path). For each real precursor carrying `target_mod_name`, move one
occurrence to a non-acceptor site (`choose_decoy_site`, sign randomized per row for
+/- balance) and emit:

- a **decoy precursor** — a copy of the real row with `:mods` rewritten (mod moved
  to the decoy site), `is_loc_decoy=true`, `loc_base_prec_id` = the real row. Same
  sequence/charge/**mz** (mass preserved), so it stays a co-isolated sibling; and
- its **fragments** — the real precursor's fragments with only the site-determining
  ion m/z permuted (`permute_fragment_mz`), `precursor_idx` = the new row index.

`precursor_idx == row index`, so decoys are appended after all reals and their
fragments reference the new indices. Adds `is_loc_decoy::Bool` and
`loc_base_prec_id::UInt32` columns (false/0 for reals). No-op when `var_mods` is
empty or the target mass is unknown.
"""
function generate_isomer_decoys!(precursors_df::DataFrame, fragments_df::DataFrame,
                                 var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
                                 mod_masses::AbstractDict;
                                 spacing::Integer = 1,
                                 target_mod_name::AbstractString = "Unimod:21")
    R = nrow(precursors_df)
    add_flag_cols!(df, nrow_) = begin
        df[!, :is_loc_decoy]     = falses(nrow_)
        df[!, :loc_base_prec_id] = zeros(UInt32, nrow_)
    end
    if isempty(var_mods) || !haskey(mod_masses, target_mod_name)
        add_flag_cols!(precursors_df, R)
        return precursors_df, fragments_df
    end
    mod_mass = Float64(mod_masses[target_mod_name])

    fr_pid = fragments_df[!, :precursor_idx]::Vector{UInt32}
    frag_rows = Dict{UInt32, Vector{Int}}()
    @inbounds for i in 1:length(fr_pid)
        push!(get!(frag_rows, fr_pid[i], Int[]), i)
    end

    base_rows      = Int[]
    decoy_modstrs  = String[]
    decoy_frag_src = Vector{Int}[]
    decoy_frag_kd  = Tuple{Int, Int}[]
    acc_cache = Dict{String, Vector{UInt8}}()

    for r in 1:R
        ms = precursors_df[r, :mods]
        (ismissing(ms) || isempty(ms)) && continue
        mods = _parse_mod_tuples(ms)
        tsites = Int[p for (p, a, n) in mods if n == target_mod_name]
        isempty(tsites) && continue
        seq = String(precursors_df[r, :sequence])
        acc = get!(() -> acceptor_positions(matchVarMods(seq, var_mods)), acc_cache, seq)
        k = tsites[1]
        d = choose_decoy_site(length(seq), acc, k, spacing, isodd(r))
        d == 0 && continue
        newmods = [(p == k && n == target_mod_name) ? (d, seq[d], n) : (p, a, n)
                   for (p, a, n) in mods]
        push!(base_rows, r)
        push!(decoy_modstrs, _format_mod_tuples(newmods))
        push!(decoy_frag_src, get(frag_rows, UInt32(r), Int[]))
        push!(decoy_frag_kd, (k, d))
    end

    add_flag_cols!(precursors_df, R)
    ndec = length(base_rows)
    ndec == 0 && return precursors_df, fragments_df

    # decoy precursor rows: copy the base reals (incl the new flag cols), override.
    decoy_prec = precursors_df[base_rows, :]
    decoy_prec[!, :mods]             = decoy_modstrs
    decoy_prec[!, :is_loc_decoy]     = trues(ndec)
    decoy_prec[!, :loc_base_prec_id] = UInt32.(base_rows)
    # Fresh unique pair_id per loc-decoy: inheriting the base target's pair_id would
    # make 3 rows share it and break add_pair_indices! (which only partners pair_ids
    # with exactly 2 rows), destroying the real target<->decoy pairing. Unique ids
    # leave loc-decoys as unpaired singletons and the real pairs intact.
    if hasproperty(precursors_df, :pair_id)
        mx = UInt32(0)
        @inbounds for x in precursors_df[!, :pair_id]
            ismissing(x) || (mx = max(mx, UInt32(x)))
        end
        decoy_prec[!, :pair_id] = UInt32[mx + UInt32(i) for i in 1:ndec]
    end

    # decoy fragment rows: copy source frags, permute m/z, repoint precursor_idx.
    total = sum(length, decoy_frag_src)
    src_rows = Vector{Int}(undef, total)
    new_mz   = Vector{Float32}(undef, total)
    new_pid  = Vector{UInt32}(undef, total)
    ann_col  = fragments_df[!, :annotation]
    mz_col   = fragments_df[!, :mz]::Vector{Float32}
    pos = 0
    for di in 1:ndec
        r = base_rows[di]; (k, d) = decoy_frag_kd[di]
        L = length(String(precursors_df[r, :sequence]))
        dpid = UInt32(R + di)
        for fi in decoy_frag_src[di]
            pos += 1
            src_rows[pos] = fi
            pa = parse_fragment_annotation(GenericFragAnnotation(String(ann_col[fi])))
            new_mz[pos]  = permute_fragment_mz(mz_col[fi], pa.base_type, Int(pa.frag_index),
                                               Int(pa.charge), L, k, d, mod_mass)
            new_pid[pos] = dpid
        end
    end
    decoy_frag = fragments_df[src_rows, :]
    decoy_frag[!, :mz]            = new_mz
    decoy_frag[!, :precursor_idx] = new_pid

    append!(precursors_df, decoy_prec)
    append!(fragments_df, decoy_frag)
    return precursors_df, fragments_df
end
