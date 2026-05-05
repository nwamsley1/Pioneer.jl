# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# 1:1 target↔decoy pair regeneration. The library's `pair_id` is many-to-one
# in places, which breaks the FTR controller's null distribution
# (target↔target "pairs" sneak in). This module regenerates a clean 1:1
# pair_id with stratification on (prec_mz × irt_pred). Within each stratum
# we pair `min(|T|, |D|)` precursors via a shuffled zip; the leftover side
# (when one outnumbers the other) gets standalone `pair_id`s with no
# partner. No cloning, no row duplication.
#
# Develop has one PSM per (precursor, file) after MainSearch's reduction, so
# the donor-tracking key is `pair_id` alone — no `isotopes_captured`
# composite key is needed (unlike v0.6.6).

using Random

"""
    regenerate_pair_ids!(psms::DataFrame, precursors::LibraryPrecursors;
                         rng_seed::Int=1776) -> NamedTuple

Regenerate `:pair_id` on `psms` in-place with a 1:1 target↔decoy invariant.
For each 10×10 (prec_mz × irt_pred) stratum, shuffle targets and decoys, zip
them to make `min(|T|, |D|)` paired pair_ids; remaining targets/decoys get
unique standalone pair_ids.

Adds one column to `psms`:
- `:pair_id::UInt32` — fresh unique pair identifier (overwrites library)

Returns a NamedTuple with diagnostic counts.
"""
function regenerate_pair_ids!(
    psms::DataFrame,
    precursors::LibraryPrecursors;
    rng_seed::Int = 1776,
)
    rng = MersenneTwister(rng_seed)
    t0 = time()

    # ── 1. Build precursor-level table from PSMs (one row per unique pid) ──
    prec_idx_col   = psms[!, :precursor_idx]::Vector{UInt32}
    target_col_psm = psms[!, :target]::Vector{Bool}
    seen = Set{UInt32}()
    plist_pids   = UInt32[]
    plist_target = Bool[]
    plist_mz     = Float32[]
    plist_irt    = Float32[]

    prec_mz_full  = getMz(precursors)
    prec_irt_full = getIrt(precursors)

    @inbounds for i in eachindex(prec_idx_col)
        pid = prec_idx_col[i]
        if pid in seen
            continue
        end
        push!(seen, pid)
        push!(plist_pids, pid)
        push!(plist_target, target_col_psm[i])
        push!(plist_mz, prec_mz_full[pid])
        push!(plist_irt, prec_irt_full[pid])
    end
    n_precs = length(plist_pids)

    # ── 2. 10×10 stratification on (prec_mz × irt_pred) ──
    function _bin_assign(values::Vector{Float32})
        finite = filter(isfinite, values)
        if isempty(finite)
            return fill(Int8(1), length(values))
        end
        edges = quantile(finite, 0:0.1:1)
        edges = unique(edges)
        if length(edges) < 11
            mn, mx = extrema(finite)
            edges = collect(LinRange(mn, mx, 11))
        end
        out = Vector{Int8}(undef, length(values))
        @inbounds for i in eachindex(values)
            v = values[i]
            if !isfinite(v)
                out[i] = Int8(5); continue
            end
            b = searchsortedlast(edges, v)
            out[i] = Int8(clamp(b, 1, 10))
        end
        return out
    end
    bin_mz  = _bin_assign(plist_mz)
    bin_irt = _bin_assign(plist_irt)

    stratum = Vector{Int}(undef, n_precs)
    @inbounds for i in 1:n_precs
        stratum[i] = 10 * Int(bin_mz[i]) + Int(bin_irt[i])
    end

    strata_t = Dict{Int, Vector{Int}}()
    strata_d = Dict{Int, Vector{Int}}()
    @inbounds for i in 1:n_precs
        s = stratum[i]
        d = plist_target[i] ? strata_t : strata_d
        push!(get!(d, s, Int[]), i)
    end

    # ── 3. Stratified 1:1 pairing — no cloning ──
    # pid_to_pair_id[pid] = pair_id (each precursor gets exactly one).
    # paired_pids = set of precursors that got a target↔decoy partner.
    next_pair_id = UInt32(1)
    pid_to_pair_id = Dict{UInt32, UInt32}()
    sizehint!(pid_to_pair_id, n_precs)
    n_pairs   = 0
    n_lone_t  = 0
    n_lone_d  = 0
    n_strata_total = 0

    for s in keys(merge(strata_t, strata_d))
        n_strata_total += 1
        ts = copy(get(strata_t, s, Int[]))
        ds = copy(get(strata_d, s, Int[]))
        shuffle!(rng, ts); shuffle!(rng, ds)
        n_pair_here = min(length(ts), length(ds))
        # Paired: each (target, decoy) gets the same fresh pair_id
        @inbounds for k in 1:n_pair_here
            pid_t = plist_pids[ts[k]]
            pid_d = plist_pids[ds[k]]
            pid_to_pair_id[pid_t] = next_pair_id
            pid_to_pair_id[pid_d] = next_pair_id
            next_pair_id += UInt32(1)
            n_pairs += 1
        end
        # Lone targets (excess on target side)
        @inbounds for k in (n_pair_here + 1):length(ts)
            pid_to_pair_id[plist_pids[ts[k]]] = next_pair_id
            next_pair_id += UInt32(1)
            n_lone_t += 1
        end
        # Lone decoys (excess on decoy side)
        @inbounds for k in (n_pair_here + 1):length(ds)
            pid_to_pair_id[plist_pids[ds[k]]] = next_pair_id
            next_pair_id += UInt32(1)
            n_lone_d += 1
        end
    end

    # ── 4. Apply :pair_id column to all PSM rows ──
    psms[!, :pair_id] = UInt32[get(pid_to_pair_id, p, UInt32(0)) for p in prec_idx_col]

    elapsed = time() - t0

    # ── 5. Invariant check: per (pair_id, ms_file_idx), at most 1 target + 1 decoy ──
    n_violations = 0
    if hasproperty(psms, :ms_file_idx) && hasproperty(psms, :target)
        gb = combine(groupby(psms, [:pair_id, :ms_file_idx]),
                     :target => sum => :n_target,
                     :target => (x -> sum(.!x)) => :n_decoy)
        n_violations = count(i -> gb.n_target[i] > 1 || gb.n_decoy[i] > 1, 1:nrow(gb))
    end

    @user_info "MBR Phase 1 — pair regeneration:"
    @user_info "  unique precursors: $n_precs"
    @user_info "  paired (T+D): $n_pairs   lone targets: $n_lone_t   lone decoys: $n_lone_d"
    @user_info "  strata visited: $n_strata_total"
    if n_violations > 0
        @user_warn "  pair_id × ms_file_idx invariant violations: $n_violations"
    else
        @user_info "  pair_id × ms_file_idx invariant: OK (no group has >1 target or >1 decoy)"
    end
    @user_info "  regenerate_pair_ids! elapsed: $(round(elapsed, digits=2))s"

    return (
        n_precs                = n_precs,
        n_pairs                = n_pairs,
        n_lone_targets         = n_lone_t,
        n_lone_decoys          = n_lone_d,
        n_strata               = n_strata_total,
        n_invariant_violations = n_violations,
        elapsed_s              = elapsed,
    )
end
