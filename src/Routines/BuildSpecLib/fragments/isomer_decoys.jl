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
