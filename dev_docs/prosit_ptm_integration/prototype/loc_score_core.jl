# Site-determining-ion localization scoring -- CORE algorithm (prototype).
#
# Design refs: dev_docs/prosit_ptm_integration/localization_flr_model_plan.pdf
#   - distinguishing ion = fragment whose CUMULATIVE mod content differs between two isomers (2b/general)
#   - s(A,B) = support(A)/(support(A)+support(B)) on presence/absence, over distinguishing ions
#   - S_vs_targets = min over competitors (target: other targets; decoy: all targets incl parent) (2b(i))
#   - per-site accumulation with joint-only flag (2b(ii))
#
# Correctness-first. Single mod type here (phospho); multi-type notes inline.
# The observed spectrum is abstracted as a lookup `obs[(ion_type, index, modcount)] -> intensity`;
# in real data (extract_sds_features.jl) the same logic runs by matching each isomer's library m/z
# to observed peaks (a fragment's modcount fixes its m/z), so this core is representation-agnostic.

# ---- fragment mod content -----------------------------------------------------------------------
# A b_i ion covers residues 1..i; a y_j ion covers (L-j+1)..L. Mod content = # mod positions inside.
frag_modcount(mods::AbstractVector{<:Integer}, T::Symbol, index::Int, L::Int) =
    T === :b ? count(p -> p <= index, mods) :
    T === :y ? count(p -> p >= L - index + 1, mods) :
    error("unknown ion type $T")

# ---- s(A,B): competitive ratio on distinguishing ions -------------------------------------------
# A, B are sorted Vectors of mod positions (same composition => same length).
# Returns (s, ndist, cap) where cap = supA+supB (0 => undefined / unlocalizable vs B).
function s_pair(A::AbstractVector{<:Integer}, B::AbstractVector{<:Integer}, L::Int, obs;
                ion_types = (:b, :y))
    supA = 0.0; supB = 0.0; ndist = 0
    @inbounds for T in ion_types, i in 1:(L - 1)
        mcA = frag_modcount(A, T, i, L)
        mcB = frag_modcount(B, T, i, L)
        mcA == mcB && continue                       # not site-determining for this pair
        ndist += 1
        supA += get(obs, (T, i, mcA), 0.0)
        supB += get(obs, (T, i, mcB), 0.0)
    end
    cap = supA + supB
    return (s = cap == 0 ? NaN : supA / cap, ndist = ndist, cap = cap)
end

# ---- S score = min over competitors (undefined pairs don't constrain) ----------------------------
function S_score(A, competitors, L, obs; ion_types = (:b, :y))
    smin = Inf; any_def = false; ndist_min = 0; cap_min = Inf
    for B in competitors
        A == B && continue
        r = s_pair(A, B, L, obs; ion_types = ion_types)
        isnan(r.s) && continue
        any_def = true
        if r.s < smin
            smin = r.s; ndist_min = r.ndist; cap_min = r.cap
        end
    end
    return (S = any_def ? smin : NaN, n_distinguishing_min = ndist_min, min_pair_capacity = cap_min,
            defined = any_def)
end

# competitor sets per the plan (2b(i)):
#   target A  -> the OTHER targets (drops itself)
#   decoy  A' -> ALL targets, INCLUDING its parent (kept)
S_vs_targets_target(A, targets, L, obs; kw...) = S_score(A, filter(!=(A), targets), L, obs; kw...)
S_vs_targets_decoy(Ap, targets, L, obs; kw...) = S_score(Ap, targets, L, obs; kw...)  # parent kept

# ---- per-site evidence (2b(ii)) ------------------------------------------------------------------
# For each assigned mod site s of A, accumulate distinguishing evidence that depends on s's occupancy
# (moving the mod off s changes the fragment's modcount). We aggregate against the competitor that is
# the eventual min (worst case); here we expose the per-site support vs a given competitor set.
#
# An ion (T,i) "depends on site s" iff s is inside the ion's residue set (its modcount includes s).
# vs a competitor B, that ion is *usable per-site* only if it distinguishes A from B AND the difference
# is attributable to s alone; if the ion also spans another differing site it is "joint-only".
function per_site_support(A::AbstractVector{<:Integer}, B::AbstractVector{<:Integer}, L::Int, obs;
                          ion_types = (:b, :y))
    diff_sites = symdiff(Set(A), Set(B))                 # sites where A and B differ
    site_supA = Dict{Int,Float64}(); site_supB = Dict{Int,Float64}()
    joint = Dict{Int,Bool}()
    for s in intersect(Set(A), diff_sites); site_supA[s] = 0.0; site_supB[s] = 0.0; joint[s] = true; end
    @inbounds for T in ion_types, i in 1:(L - 1)
        mcA = frag_modcount(A, T, i, L); mcB = frag_modcount(B, T, i, L)
        mcA == mcB && continue
        # which of A's differing sites lie inside this ion?
        inside = [s for s in keys(site_supA) if (T === :b ? s <= i : s >= L - i + 1)]
        length(inside) == 1 || continue                  # spans >1 differing site -> joint-only, skip per-site
        s = inside[1]
        joint[s] = false                                 # this site has an independent distinguishing ion
        site_supA[s] += get(obs, (T, i, mcA), 0.0)
        site_supB[s] += get(obs, (T, i, mcB), 0.0)
    end
    return (site_supA = site_supA, site_supB = site_supB, joint_only = joint)
end
