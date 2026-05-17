# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

# MBR pairing.
#
# Provides the experiment-global counterfactual-partner map used by MBR Batch F.
# For every unique `precursor_idx`, assigns exactly one partner of the OPPOSITE
# target/decoy class via iRT-NN stratified within (cv_fold × mz decile). The
# partner is written into the PSM column `:counterfactual_partner_pid` and is
# the lookup the MBR streaming feature compute uses to find each row's _false
# donor (see mbr_streaming.jl).

using Random
using Statistics

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

"""
    regenerate_counterfactual_partners_irt_nn!(psms, precursors)

MBR Batch F variant: counterfactual partner = the opposite-class precursor
with the CLOSEST predicted iRT within the same (cv_fold, prec_mz_decile)
stratum. Replaces the random-within-fine-iRT-bin assignment.

Symmetric: each target's partner is its iRT-nearest decoy; each decoy's
partner is its iRT-nearest target. M:N tying allowed.

Fallback if a stratum is one-sided:
- Tier 1: same cv_fold, ALL mz_deciles — iRT-nearest opposite from the
  fold-global opposite pool.
- Tier 2: experiment-global iRT-nearest.

Hard requirement: every precursor gets a partner (asserted).
"""
function regenerate_counterfactual_partners_irt_nn!(
    psms::DataFrame,
    precursors::LibraryPrecursors,
)
    t0 = time()

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

    # ── 2. prec_mz deciles ──
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

    # ── 3. Sorted (pid, irt) lists per (cv_fold, mz_decile, target_class) ──
    # Use NamedTuple {pids, irts} where irts is sorted ascending.
    function _pool_key(fold::UInt8, mz::Int); (Int(fold), mz); end
    pools_t = Dict{Tuple{Int,Int}, NamedTuple{(:pids, :irts), Tuple{Vector{UInt32}, Vector{Float32}}}}()
    pools_d = Dict{Tuple{Int,Int}, NamedTuple{(:pids, :irts), Tuple{Vector{UInt32}, Vector{Float32}}}}()
    # First collect, then sort.
    tmp_t = Dict{Tuple{Int,Int}, Vector{Tuple{Float32, UInt32}}}()
    tmp_d = Dict{Tuple{Int,Int}, Vector{Tuple{Float32, UInt32}}}()
    @inbounds for i in 1:n_precs
        key = _pool_key(plist_fold[i], bin_mz[i])
        dst = plist_target[i] ? tmp_t : tmp_d
        push!(get!(() -> Tuple{Float32, UInt32}[], dst, key), (plist_irt[i], plist_pids[i]))
    end
    for (key, v) in tmp_t
        sort!(v; by = x -> (x[1], x[2]))   # by irt then pid for determinism
        pools_t[key] = (pids = UInt32[x[2] for x in v], irts = Float32[x[1] for x in v])
    end
    for (key, v) in tmp_d
        sort!(v; by = x -> (x[1], x[2]))
        pools_d[key] = (pids = UInt32[x[2] for x in v], irts = Float32[x[1] for x in v])
    end

    # Fold-global pools (any mz) — built lazily via concatenate-then-sort.
    fold_pool_t = Dict{Int, NamedTuple{(:pids, :irts), Tuple{Vector{UInt32}, Vector{Float32}}}}()
    fold_pool_d = Dict{Int, NamedTuple{(:pids, :irts), Tuple{Vector{UInt32}, Vector{Float32}}}}()
    function _fold_pool!(cache::Dict, src::Dict{Tuple{Int,Int}, Any}, fold::Int)
        get!(cache, fold) do
            v = Tuple{Float32, UInt32}[]
            for ((f, _), pool) in src
                f == fold || continue
                @inbounds for k in 1:length(pool.pids); push!(v, (pool.irts[k], pool.pids[k])); end
            end
            sort!(v; by = x -> (x[1], x[2]))
            (pids = UInt32[x[2] for x in v], irts = Float32[x[1] for x in v])
        end
    end
    # Pre-build fold pools for cv_fold ∈ {0, 1}
    function _build_fold_pool(src_pools::Dict, fold::Int)
        v = Tuple{Float32, UInt32}[]
        for ((f, _), pool) in src_pools
            f == fold || continue
            @inbounds for k in 1:length(pool.pids); push!(v, (pool.irts[k], pool.pids[k])); end
        end
        sort!(v; by = x -> (x[1], x[2]))
        return (pids = UInt32[x[2] for x in v], irts = Float32[x[1] for x in v])
    end
    fold_pool_t[0] = _build_fold_pool(pools_t, 0)
    fold_pool_t[1] = _build_fold_pool(pools_t, 1)
    fold_pool_d[0] = _build_fold_pool(pools_d, 0)
    fold_pool_d[1] = _build_fold_pool(pools_d, 1)

    # Experiment-global pools
    function _global_pool(fold_pool::Dict)
        v = Tuple{Float32, UInt32}[]
        for (_, pool) in fold_pool
            @inbounds for k in 1:length(pool.pids); push!(v, (pool.irts[k], pool.pids[k])); end
        end
        sort!(v; by = x -> (x[1], x[2]))
        return (pids = UInt32[x[2] for x in v], irts = Float32[x[1] for x in v])
    end
    global_pool_t = _global_pool(fold_pool_t)
    global_pool_d = _global_pool(fold_pool_d)

    # ── 4. Nearest-iRT lookup helper ──
    # Returns the pid whose irt is closest to `target_irt` in the sorted pool;
    # 0 if pool is empty.
    @inline function _nearest_pid(pool::NamedTuple{(:pids, :irts), Tuple{Vector{UInt32}, Vector{Float32}}},
                                  target_irt::Float32)
        n = length(pool.irts)
        n == 0 && return UInt32(0)
        # binary search for insertion point
        idx = searchsortedfirst(pool.irts, target_irt)
        if idx > n
            return pool.pids[n]
        elseif idx == 1
            return pool.pids[1]
        else
            d_left  = abs(target_irt - pool.irts[idx-1])
            d_right = abs(pool.irts[idx] - target_irt)
            return d_left <= d_right ? pool.pids[idx-1] : pool.pids[idx]
        end
    end

    # ── 5. Assign partner for every precursor ──
    pid_to_partner = Dict{UInt32, UInt32}()
    sizehint!(pid_to_partner, n_precs)
    n_stratum = 0; n_fb_fold = 0; n_fb_global = 0
    @inbounds for i in 1:n_precs
        my_pid = plist_pids[i]
        my_irt = plist_irt[i]
        my_fold = Int(plist_fold[i])
        my_mz = bin_mz[i]
        opp_pool_local = plist_target[i] ?
            get(pools_d, (my_fold, my_mz), (pids=UInt32[], irts=Float32[])) :
            get(pools_t, (my_fold, my_mz), (pids=UInt32[], irts=Float32[]))

        partner_pid = _nearest_pid(opp_pool_local, my_irt)
        if partner_pid != UInt32(0)
            pid_to_partner[my_pid] = partner_pid
            n_stratum += 1
            continue
        end

        # Tier 1: same cv_fold, any mz
        opp_fold = plist_target[i] ? fold_pool_d[my_fold] : fold_pool_t[my_fold]
        partner_pid = _nearest_pid(opp_fold, my_irt)
        if partner_pid != UInt32(0)
            pid_to_partner[my_pid] = partner_pid
            n_fb_fold += 1
            continue
        end

        # Tier 2: experiment-global
        opp_global = plist_target[i] ? global_pool_d : global_pool_t
        partner_pid = _nearest_pid(opp_global, my_irt)
        if partner_pid == UInt32(0)
            error("regenerate_counterfactual_partners_irt_nn!: no opposite-class partner exists anywhere for pid=$my_pid (target=$(plist_target[i])).")
        end
        pid_to_partner[my_pid] = partner_pid
        n_fb_global += 1
    end

    psms[!, :counterfactual_partner_pid] = UInt32[pid_to_partner[p] for p in prec_idx_col]

    elapsed = time() - t0
    pct_strat = round(100 * n_stratum / max(n_precs, 1), digits=2)
    pct_fbf   = round(100 * n_fb_fold / max(n_precs, 1), digits=3)
    pct_fbg   = round(100 * n_fb_global / max(n_precs, 1), digits=3)

    @user_info "MBR Batch F (iRT-NN) — counterfactual partner map:"
    @user_info "  unique precursors: $n_precs"
    @user_info "  partner = closest-pred-iRT opposite within (cv_fold × mz_decile)"
    @user_info "  in stratum:                     $n_stratum ($pct_strat%)"
    @user_info "  fallback (cv_fold, any mz):     $n_fb_fold ($pct_fbf%)"
    @user_info "  fallback (experiment-global):   $n_fb_global ($pct_fbg%)"
    @user_info "  total partnered: $(length(pid_to_partner)) / $n_precs"
    @user_info "  regenerate_counterfactual_partners_irt_nn! elapsed: $(round(elapsed, digits=2))s"

    @assert length(pid_to_partner) == n_precs "iRT-NN partner map missing $(n_precs - length(pid_to_partner)) precursors"

    return (
        n_precs       = n_precs,
        n_in_stratum  = n_stratum,
        n_fb_fold     = n_fb_fold,
        n_fb_global   = n_fb_global,
        elapsed_s     = elapsed,
    )
end
