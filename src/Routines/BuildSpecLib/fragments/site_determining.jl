# Site-determining fragment retention (phospho localization).
#
# Ensures the library keeps, per peptidoform, the fragments that distinguish
# positional isomers (and, later, isomer decoys): top-N by intensity PLUS the
# highest-intensity fragment crossing each gap between consecutive candidate
# modification sites. See dev_docs/prosit_ptm_integration/
# site_determining_fragment_retention.md (algorithm + worked VLSAGSPESIK example)
# and site_determining_retention_impl_plan.md (this implementation).
#
# The candidate-site set S is computed at peptidoform GENERATION from the same
# variable-mod matches the generator enumerates over (authoritative), optionally
# unioned with decoy-candidate neighbours, then carried to prediction. The
# retention itself (gap_cover_indices) is pure and consumes the carried S.

"""
    acceptor_positions(var_mod_matches) -> Vector{UInt8}

Sorted, unique 1-based positions where a variable modification can be placed,
taken from the regex-match offsets produced by `matchVarMods` (the same source
the peptidoform generator uses). Authoritative: honours the configured var-mod
motifs rather than re-scanning for residue identity.
"""
function acceptor_positions(var_mod_matches)
    isempty(var_mod_matches) && return UInt8[]
    pos = UInt8[UInt8(m.regex_match.offset) for m in var_mod_matches]
    return sort!(unique!(pos))
end

"""
    decoy_neighbor_positions(seq, acc_positions) -> Vector{UInt8}

Decoy-candidate positions: the immediate sequence neighbours (`p-1`, `p+1`) of
each acceptor that are **not themselves acceptors** (so a phospho placed there is
wrong by construction). This is the shared helper the isomer-decoy generator
must also use, so the reals retain the ions that distinguish the decoys.
"""
function decoy_neighbor_positions(seq::AbstractString, acc_positions::AbstractVector{<:Integer})
    L = length(seq)
    accset = Set{Int}(Int.(acc_positions))
    neigh = UInt8[]
    for p in acc_positions
        for j in (Int(p) - 1, Int(p) + 1)
            if 1 <= j <= L && !(j in accset)
                push!(neigh, UInt8(j))
            end
        end
    end
    return sort!(unique!(neigh))
end

"""
    compute_gap_sites(seq, var_mod_matches, include_decoy) -> Vector{UInt8}

The candidate-site set `S` for one base peptide: acceptor positions, optionally
unioned with decoy-candidate neighbours, sorted and unique. Consecutive pairs of
`S` are the "gaps" whose crossing ions must be retained. Empty when there are no
acceptors (retention is then a strict no-op).
"""
function compute_gap_sites(seq::AbstractString, var_mod_matches, include_decoy::Bool)
    acc = acceptor_positions(var_mod_matches)
    isempty(acc) && return UInt8[]
    include_decoy || return acc
    return sort!(unique!(vcat(acc, decoy_neighbor_positions(seq, acc))))
end

"""
    gap_cover_indices(gap_sites, L, base_type, frag_index, block, topn_cut) -> Vector{Int}

Given a precursor's carried `gap_sites` and its kept fragments (`base_type`,
`frag_index` indexed by `block`, where `block` is fragment-row indices sorted by
**descending intensity**, of which the first `topn_cut` are already retained as
top-N), return the **extra** fragment-row indices to retain: for each gap, the
single highest-intensity fragment whose cleavage falls in the gap, if it is not
already in the top-N.

Cleavage of a fragment: `c = frag_index` for a b-ion, `c = L - frag_index` for a
y-ion (a y-ion of ordinal `m` cleaves after residue `L-m`). A fragment crosses
gap `(a, b)` iff `a <= c < b`. Precursor/other ion types never cross. Because
each fragment's cleavage lies in exactly one gap, no fragment is added twice.
Pure: no ctx, no I/O.
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

"""
    gap_sites_for_library(sequences, var_mods, include_decoy) -> Vector{Vector{UInt8}}

Compute the carried `gap_sites` (candidate-site set `S`) once per precursor from
the authoritative variable-mod patterns — the same `matchVarMods` the generator
uses — so the retention filter never re-derives acceptor identity. Returns one
`Vector{UInt8}` per precursor (empty when there are no variable mods, making the
whole feature a strict no-op).
"""
function gap_sites_for_library(sequences::AbstractVector,
                               var_mods::Vector{NamedTuple{(:p, :r), Tuple{Regex, String}}},
                               include_decoy::Bool)
    n = length(sequences)
    out = Vector{Vector{UInt8}}(undef, n)
    if isempty(var_mods)
        @inbounds for i in 1:n; out[i] = UInt8[]; end
        return out
    end
    @inbounds for i in 1:n
        seq = String(sequences[i])
        out[i] = compute_gap_sites(seq, matchVarMods(seq, var_mods), include_decoy)
    end
    return out
end
