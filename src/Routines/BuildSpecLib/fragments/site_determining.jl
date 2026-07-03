# Site-determining fragment retention + shadow decoy sites (phospho localization).
#
# Ensures the library keeps, per peptidoform, the fragments that distinguish
# competing placements: top-N by intensity PLUS the highest-intensity fragment
# crossing each gap between consecutive candidate sites. For the isomer-decoy
# scheme the candidate set is the UNION of real acceptor sites and the (random,
# count-matched) decoy "shadow" sites, so real-vs-real, real-vs-decoy and
# decoy-vs-decoy placements are all distinguishable.
#
# The shadow sites are a DETERMINISTIC function of the sequence (seeded hash), so
# retention (compute_gap_sites) and the decoy generator (generate_isomer_decoys!)
# independently recompute the SAME sites with no shared state.
#
# See dev_docs/prosit_ptm_integration/isomer_decoy_shadow_sites.md and
# site_determining_fragment_retention.md.

const _SHADOW_GLOBAL_SEED = 0x243f6a8885a308d3

"""
    acceptor_positions(var_mod_matches) -> Vector{UInt8}

Sorted, unique 1-based positions where a variable modification can be placed,
from `matchVarMods` regex-match offsets (the same source the generator uses).
"""
function acceptor_positions(var_mod_matches)
    isempty(var_mod_matches) && return UInt8[]
    pos = UInt8[UInt8(m.regex_match.offset) for m in var_mod_matches]
    return sort!(unique!(pos))
end

"""
    seq_seed(seq; global_seed) -> UInt64

Stable (Julia-version-independent) FNV-1a hash of the sequence, used to seed the
per-peptide shadow-site draw so it is reproducible across builds.
"""
function seq_seed(seq::AbstractString; global_seed::UInt64 = _SHADOW_GLOBAL_SEED)
    h = global_seed
    @inbounds for c in codeunits(String(seq))
        h = (h ⊻ UInt64(c)) * 0x100000001b3
    end
    return h
end

"""
    real_sites_by_mod(seq, var_mods) -> Dict{String, Vector{Int}}

Real acceptor positions per modification NAME, from the variable-mod regexes
(offsets). `var_mods` is the `(p=Regex, r=name)` list the generator enumerates.
"""
function real_sites_by_mod(seq::AbstractString, var_mods)
    d = Dict{String, Vector{Int}}()
    s = String(seq)
    for vm in var_mods
        for m in eachmatch(vm.p, s)
            push!(get!(d, String(vm.r), Int[]), m.offset)
        end
    end
    for k in collect(keys(d)); d[k] = sort!(unique!(d[k])); end
    return d
end

"""
    shadow_decoy_sites(seq, real_by; seed) -> Dict{String, Vector{Int}}

Per modification NAME, draw as many random decoy positions as there are real
acceptor sites of that mod, from positions that are NOT a real acceptor of ANY
mod and NOT already a decoy site (all decoy sites mutually distinct). Sorted,
deterministic given `seed`. If a mod has fewer free positions than real sites,
it gets as many as fit (logged as a shortfall by the caller if needed).
"""
function shadow_decoy_sites(seq::AbstractString, real_by::AbstractDict; seed::UInt64)
    L = length(seq)
    blocked = Set{Int}()
    for v in values(real_by); union!(blocked, v); end
    rng = Random.MersenneTwister(seed)
    out = Dict{String, Vector{Int}}()
    for name in sort(collect(keys(real_by)))          # deterministic mod order
        n = length(real_by[name])
        avail = [p for p in 1:L if !(p in blocked)]
        k = min(n, length(avail))
        chosen = k == 0 ? Int[] : sort!(Random.shuffle(rng, avail)[1:k])
        out[name] = chosen
        union!(blocked, chosen)
    end
    return out
end

"""
    compute_gap_sites(seq, var_mods, include_decoy; global_seed) -> Vector{UInt8}

Candidate-site set `S`. Without decoys: the real acceptor positions. With decoys:
the UNION of real acceptor positions and the shadow decoy positions (so retention
covers real-vs-decoy boundaries too). Empty when there are no acceptors (no-op).
"""
function compute_gap_sites(seq::AbstractString, var_mods, include_decoy::Bool;
                           global_seed::UInt64 = _SHADOW_GLOBAL_SEED)
    real_by = real_sites_by_mod(seq, var_mods)
    isempty(real_by) && return UInt8[]
    all_pos = Int[]
    for v in values(real_by); append!(all_pos, v); end
    if include_decoy
        shadow = shadow_decoy_sites(seq, real_by; seed = seq_seed(seq; global_seed))
        for v in values(shadow); append!(all_pos, v); end
    end
    return UInt8.(sort!(unique!(all_pos)))
end

"""
    gap_sites_for_library(sequences, var_mods, include_decoy; global_seed) -> Vector{Vector{UInt8}}

Per-precursor candidate-site set, carried into the filter ctx. Empty (per row and
overall) when there are no variable mods -> retention is a strict no-op.
"""
function gap_sites_for_library(sequences::AbstractVector,
                               var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
                               include_decoy::Bool;
                               global_seed::UInt64 = _SHADOW_GLOBAL_SEED)
    n = length(sequences)
    out = Vector{Vector{UInt8}}(undef, n)
    if isempty(var_mods)
        @inbounds for i in 1:n; out[i] = UInt8[]; end
        return out
    end
    @inbounds for i in 1:n
        out[i] = compute_gap_sites(String(sequences[i]), var_mods, include_decoy; global_seed)
    end
    return out
end

"""
    gap_cover_indices(gap_sites, L, base_type, frag_index, block, topn_cut) -> Vector{Int}

For each gap between consecutive `gap_sites`, the extra fragment-row index to
retain beyond the top-N: the highest-intensity fragment whose cleavage falls in
the gap (b-ion cleavage = frag_index, y-ion = L - frag_index), if not already in
the top-N. `block` = kept-fragment rows sorted by descending intensity; first
`topn_cut` already kept. Pure.
"""
function gap_cover_indices(gap_sites::AbstractVector{<:Integer},
                           L::Integer,
                           base_type::AbstractVector{Char},
                           frag_index::AbstractVector{<:Integer},
                           block::AbstractVector{<:Integer},
                           topn_cut::Integer)
    extra = Int[]
    length(gap_sites) < 2 && return extra
    n_top = min(Int(topn_cut), length(block))
    @inbounds for g in 1:(length(gap_sites) - 1)
        a = Int(gap_sites[g]); b = Int(gap_sites[g + 1])
        for pos in 1:length(block)
            k = Int(block[pos])
            bt = base_type[k]
            (bt == 'b' || bt == 'y') || continue
            c = bt == 'b' ? Int(frag_index[k]) : Int(L) - Int(frag_index[k])
            if a <= c < b
                pos > n_top && push!(extra, k)   # first (highest-intensity) crosser
                break
            end
        end
    end
    return extra
end
