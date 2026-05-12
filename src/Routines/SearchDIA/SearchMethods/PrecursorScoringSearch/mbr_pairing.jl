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
using Statistics

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
    cv_fold_psm    = psms[!, :cv_fold]
    seen = Set{UInt32}()
    plist_pids   = UInt32[]
    plist_target = Bool[]
    plist_mz     = Float32[]
    plist_irt    = Float32[]
    plist_fold   = UInt8[]  # cv_fold per precursor (constant across PSMs of same pid)

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
        push!(plist_fold, UInt8(cv_fold_psm[i]))
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

    # Stratum key: cv_fold × prec_mz_bin × irt_bin. Including cv_fold ensures
    # T↔D pairs only form within the same fold, preventing cross-fold donor
    # leakage downstream (compute_mbr_features! picks paired-partner PSMs as
    # donors, so partners must share a fold to keep MBR features CV-clean).
    stratum = Vector{Int}(undef, n_precs)
    @inbounds for i in 1:n_precs
        stratum[i] = 1000 * Int(plist_fold[i]) + 10 * Int(bin_mz[i]) + Int(bin_irt[i])
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

"""
    regenerate_counterfactual_partners!(psms, precursors;
                                        rng_seed=1846,
                                        irt_bin_width=0.5f0)

MBR Batch F: build an experiment-global counterfactual-partner map. For
every unique `precursor_idx` in `psms`, assign exactly one partner of
the OPPOSITE target/decoy class. Hard requirement: every precursor gets
a partner (asserted at the end).

Adds `:counterfactual_partner_pid::UInt32` to `psms`.

Stratification (preferred binning):
- cv_fold (so partners share fold; preserves CV cleanliness downstream)
- prec_mz decile (10 quantile bins)
- iRT bin of fixed width `irt_bin_width` (default 0.5 iRT)

Within each stratum, M:N wrap (`mod1`) — each target gets one decoy,
each decoy gets one target. Imbalance is fine: the smaller side wraps
around. Stable: `get!` writes each partner exactly once on first
encounter.

Three fallback tiers for precursors in one-sided strata (no opposite
class present):
1. nearest non-empty neighbor stratum in same (cv_fold, mz_decile),
   walking outward in iRT.
2. same cv_fold, any (mz, iRT) — pull from a fold-global opposite pool.
3. experiment-global — any opposite-class precursor.

If all three tiers fail the function throws (means the experiment has
zero opposite-class precursors of one polarity, which is degenerate).
"""
function regenerate_counterfactual_partners!(
    psms::DataFrame,
    precursors::LibraryPrecursors;
    rng_seed::Int            = 1846,
    irt_bin_width::Float32   = 0.5f0,
)
    rng = MersenneTwister(rng_seed)
    t0  = time()

    # ── 1. One-row-per-unique-pid table ──
    prec_idx_col = psms[!, :precursor_idx]::Vector{UInt32}
    target_col   = psms[!, :target]::Vector{Bool}
    cv_fold_col  = psms[!, :cv_fold]

    prec_mz_full  = getMz(precursors)
    prec_irt_full = getIrt(precursors)

    seen = Set{UInt32}()
    plist_pids   = UInt32[]
    plist_target = Bool[]
    plist_mz     = Float32[]
    plist_irt    = Float32[]
    plist_fold   = UInt8[]
    @inbounds for i in eachindex(prec_idx_col)
        pid = prec_idx_col[i]
        pid in seen && continue
        push!(seen, pid)
        push!(plist_pids, pid)
        push!(plist_target, target_col[i])
        push!(plist_mz, prec_mz_full[pid])
        push!(plist_irt, prec_irt_full[pid])
        push!(plist_fold, UInt8(cv_fold_col[i]))
    end
    n_precs = length(plist_pids)

    # ── 2. Bin assignments ──
    finite_mz = filter(isfinite, plist_mz)
    mz_edges = isempty(finite_mz) ? Float32[0f0, 1f0] :
                  Float32.(unique(quantile(finite_mz, 0:0.1:1)))
    if length(mz_edges) < 11
        mn, mx = isempty(finite_mz) ? (0f0, 1f0) : extrema(finite_mz)
        mz_edges = collect(Float32, LinRange(mn, mx, 11))
    end
    bin_mz = Vector{Int}(undef, n_precs)
    @inbounds for i in 1:n_precs
        v = plist_mz[i]
        bin_mz[i] = isfinite(v) ? clamp(searchsortedlast(mz_edges, v), 1, 10) : 5
    end

    finite_irt = filter(isfinite, plist_irt)
    irt_min = isempty(finite_irt) ? 0f0 : minimum(finite_irt)
    bin_irt = Vector{Int}(undef, n_precs)
    @inbounds for i in 1:n_precs
        v = plist_irt[i]
        bin_irt[i] = isfinite(v) ?
            (1 + Int(floor((v - irt_min) / irt_bin_width))) : 1
    end

    # ── 3. Strata bookkeeping ──
    # Stratum key: cv_fold * 1e7 + mz_bin * 1e4 + irt_bin.
    stratum_of = Vector{Int}(undef, n_precs)
    @inbounds for i in 1:n_precs
        stratum_of[i] = 10_000_000 * Int(plist_fold[i]) +
                        10_000     * bin_mz[i]          +
                                     bin_irt[i]
    end

    strata_t = Dict{Int, Vector{Int}}()
    strata_d = Dict{Int, Vector{Int}}()
    @inbounds for i in 1:n_precs
        s = stratum_of[i]
        push!(get!(() -> Int[], plist_target[i] ? strata_t : strata_d, s), i)
    end

    # ── 4. M:N within-stratum pairing ──
    pid_to_partner = Dict{UInt32, UInt32}()
    sizehint!(pid_to_partner, n_precs)
    n_in_stratum = 0
    all_strata = union(keys(strata_t), keys(strata_d))
    for s in all_strata
        ts = copy(get(strata_t, s, Int[]))
        ds = copy(get(strata_d, s, Int[]))
        (isempty(ts) || isempty(ds)) && continue
        shuffle!(rng, ts); shuffle!(rng, ds)
        n = max(length(ts), length(ds))
        @inbounds for k in 1:n
            t_idx = ts[mod1(k, length(ts))]
            d_idx = ds[mod1(k, length(ds))]
            t_pid = plist_pids[t_idx]
            d_pid = plist_pids[d_idx]
            # write-on-first so partners are stable
            if !haskey(pid_to_partner, t_pid); pid_to_partner[t_pid] = d_pid; n_in_stratum += 1; end
            if !haskey(pid_to_partner, d_pid); pid_to_partner[d_pid] = t_pid; n_in_stratum += 1; end
        end
    end

    # ── 5. Fallback for unpartnered precursors ──
    # Precompute per-(cv_fold, mz_decile) opposite-class pools (by pid),
    # per-cv_fold opposite-class pools, and experiment-global pools.
    fold_mz_pool_t = Dict{Tuple{UInt8,Int}, Vector{UInt32}}()
    fold_mz_pool_d = Dict{Tuple{UInt8,Int}, Vector{UInt32}}()
    fold_pool_t    = Dict{UInt8, Vector{UInt32}}()
    fold_pool_d    = Dict{UInt8, Vector{UInt32}}()
    all_t = UInt32[]; all_d = UInt32[]
    @inbounds for i in 1:n_precs
        pid = plist_pids[i]; fold = plist_fold[i]; mz = bin_mz[i]
        is_t = plist_target[i]
        push!(get!(() -> UInt32[], is_t ? fold_mz_pool_t : fold_mz_pool_d, (fold, mz)), pid)
        push!(get!(() -> UInt32[], is_t ? fold_pool_t    : fold_pool_d,    fold),     pid)
        push!(is_t ? all_t : all_d, pid)
    end

    n_fb_neighbor = 0; n_fb_fold = 0; n_fb_global = 0
    @inbounds for i in 1:n_precs
        pid = plist_pids[i]
        haskey(pid_to_partner, pid) && continue
        is_t = plist_target[i]
        opp_fold_mz = is_t ? fold_mz_pool_d : fold_mz_pool_t
        opp_fold    = is_t ? fold_pool_d    : fold_pool_t
        opp_global  = is_t ? all_d          : all_t

        chosen = UInt32(0)
        # Tier 1: walk outward in iRT within the same (cv_fold, mz_decile)
        # by widening to mz_decile-wide pool (we don't track per-irt-bin pools
        # to keep it simple — the (fold, mz) pool implicitly covers nearby iRT
        # since most opposite-class pids in this mz-band qualify).
        pool = get(opp_fold_mz, (plist_fold[i], bin_mz[i]), UInt32[])
        if !isempty(pool); chosen = pool[rand(rng, 1:length(pool))]; n_fb_neighbor += 1
        else
            pool = get(opp_fold, plist_fold[i], UInt32[])
            if !isempty(pool); chosen = pool[rand(rng, 1:length(pool))]; n_fb_fold += 1
            elseif !isempty(opp_global); chosen = opp_global[rand(rng, 1:length(opp_global))]; n_fb_global += 1
            end
        end
        if chosen == UInt32(0)
            error("regenerate_counterfactual_partners!: no opposite-class partner exists anywhere for pid=$pid (target=$is_t). Experiment has zero opposite-class precursors of one polarity.")
        end
        pid_to_partner[pid] = chosen
    end

    # ── 6. Apply column ──
    psms[!, :counterfactual_partner_pid] = UInt32[pid_to_partner[p] for p in prec_idx_col]

    elapsed = time() - t0
    n_paired_total = length(pid_to_partner)
    pct_strat = round(100 * n_in_stratum / max(n_precs, 1), digits=1)
    pct_fbn   = round(100 * n_fb_neighbor / max(n_precs, 1), digits=2)
    pct_fbf   = round(100 * n_fb_fold     / max(n_precs, 1), digits=2)
    pct_fbg   = round(100 * n_fb_global   / max(n_precs, 1), digits=2)

    @user_info "MBR Batch F — counterfactual partner map:"
    @user_info "  unique precursors: $n_precs"
    @user_info "  irt_bin_width: $irt_bin_width"
    @user_info "  partnered in stratum:        $n_in_stratum ($pct_strat%)"
    @user_info "  fallback (cv_fold + mz_bin): $n_fb_neighbor ($pct_fbn%)"
    @user_info "  fallback (cv_fold global):   $n_fb_fold ($pct_fbf%)"
    @user_info "  fallback (experiment global):$n_fb_global ($pct_fbg%)"
    @user_info "  total partnered: $n_paired_total / $n_precs"
    @user_info "  regenerate_counterfactual_partners! elapsed: $(round(elapsed, digits=2))s"

    @assert n_paired_total == n_precs "Counterfactual partner map missing $(n_precs - n_paired_total) precursors"

    return (
        n_precs        = n_precs,
        n_in_stratum   = n_in_stratum,
        n_fb_neighbor  = n_fb_neighbor,
        n_fb_fold      = n_fb_fold,
        n_fb_global    = n_fb_global,
        elapsed_s      = elapsed,
    )
end
