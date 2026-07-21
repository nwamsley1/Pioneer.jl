# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.
#
# Sciex ZT scanning DIA: collapse the ±k-bin meta-scan of each precursor (one
# center bin — where the precursor m/z falls in the ±half-bin window — plus its
# 2k neighboring Q1 bins in the same cycle) into a single "meta-PSM" row before
# the per-file LightGBM. The meta-PSM keeps the CENTER bin's feature columns and
# gains (a) the raw 2k+1 meta-scan weights and (b) shape features describing the
# weight profile, which discriminate real precursors (sharp centered triangle)
# from decoys/interference (flat, spread). First (deliberately simple) version:
# store the whole profile so features are cheap to iterate on.

# Shape features derived from the meta-scan weight profile.
const ZT_PROFILE_FEATURES = Symbol[
    :zt_tri_cosine,     # cosine of the profile with an ideal triangle template
    :zt_center_frac,    # center-bin weight / total
    :zt_tail_frac,      # weight at |j|≥4 / total (edge heaviness)
    :zt_apex_offset,    # j of the max-weight bin (0 = peaked at center)
    :zt_centroid,       # weighted mean offset Σ j·w / Σw
    :zt_spread,         # weighted std of the offset (peak width)
    :zt_entropy,        # Shannon entropy of the normalized profile
    :zt_symmetry,       # 1 - |Σ_{j<0}w - Σ_{j>0}w| / Σw
    :zt_monotonicity,   # fraction of outward steps where weight decreases
]

# Within-metascan (shape) fragment features: across the 2k+1 bins of ONE cycle,
# how well does each fragment's intensity profile track the weight/transmission
# triangle? Mirrors the develop across-precursor frag-corr aggregation
# (_frag_pcor / _positive_corr_summary / _bitvec_pattern_rank in features.jl) but on
# the per-cycle profile. The elution (across-cycle) counterparts are added in a
# separate post-collapse pass. `_shape` suffix parallels the develop feature names.
const ZT_SHAPE_FEATURES = Symbol[
    :frag_corr_strength_shape,               # rank-weighted positive-corr strength
    :frag_corr_effective_n_shape,            # effective # of positively-correlated frags
    :frag_corr_best_shape,                   # best-consensus fragment's shape corr vs weight
    :frag_apex_dispersion_shape,             # std of fragment apex bin-offsets (0 = all centered)
    :n_correlated_fragments_shape,           # # fragments with shape corr > 0.7
    :n_correlated_fragments_bitvec_rank_shape,  # bitvec rank of the >0.7 pattern
    # weight-profile-vs-triangle match, mean-centered Pearson (not a fragment feature —
    # grouped here for model wiring). Complements the uncentered zt_tri_cosine; a top-3
    # feature in the 2nd-pass. 3-file: +2.9% prec / +3.0% PG over the shape-only baseline.
    :zt_tri_pcor,
]

# Experimental shape-feature candidates (env `PIONEER_ZT_SHAPE_EXP`, default off →
# columns absent → hasproperty-filtered out of both LGBMs → byte-identical). Tested
# one at a time in isolation; this vector holds only the candidate currently under
# test (see zt_shape_feature_experiments_plan.md). Confirmed winners graduate to the
# permanent lists above. Next up: (a) zt_tri_spectral_angle.
const ZT_SHAPE_EXP_FEATURES = Symbol[
]

# Minimal ZT feature set the MODELS consume (2026-07-21), from the 3-file Condition-A
# feature-gain analysis: only the ZT within-metascan features that carry signal. The collapse
# still COMPUTES + stores all ZT_PROFILE/ZT_SHAPE columns (lists above); this just prunes what
# the per-file main-search LGBM uses. Dropped as noise/dead: the 8 non-cosine profile
# descriptors (zt_center_frac/tail_frac/apex_offset/centroid/spread/entropy/symmetry/monotonicity)
# and the weak shape variants (frag_corr_effective_n_shape, frag_corr_best_shape,
# frag_apex_dispersion_shape, n_correlated_fragments_bitvec_rank_shape).
const ZT_MAINSEARCH_MODEL_FEATURES = Symbol[
    :zt_tri_cosine, :zt_tri_pcor, :frag_corr_strength_shape, :n_correlated_fragments_shape,
]

# The across-cycle ELUTION features (frag_corr_*, ms1_corr_*, n_scans) are NOT a
# separate feature set: for ZT they are exactly the develop chromatogram features,
# run on the collapsed one-point-per-cycle meta trace (add_chromatogram_features!
# called post-collapse in MainSearch.jl). Develop's math on the meta trace is
# identical to a bespoke "elution" implementation, so we reuse it directly for
# maximum compatibility (same names, same code). Only the within-metascan SHAPE
# features above are genuinely new.

# Compute the 9 shape features from a length-(2k+1) weight vector `w` (j=-k..k,
# center at index k+1), given the ideal triangle template `tri` (+ its norm).
@inline function _zt_profile_features(w::Vector{Float32}, k::Int, tri::Vector{Float32}, tnorm::Float32)
    L = 2k + 1
    W = zero(Float32); @inbounds for t in 1:L; W += w[t]; end
    if W <= 0
        return (0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0)
    end
    w0 = w[k+1]
    dotwt = zero(Float32); nw = zero(Float32)
    tail = zero(Float32); centroid = zero(Float32); sl = zero(Float32); sr = zero(Float32)
    amax = 1; wmax = w[1]
    @inbounds for t in 1:L
        j = t - (k+1)
        dotwt += w[t]*tri[t]; nw += w[t]*w[t]
        abs(j) >= 4 && (tail += w[t])
        centroid += Float32(j)*w[t]
        j < 0 && (sl += w[t]); j > 0 && (sr += w[t])
        w[t] > wmax && (wmax = w[t]; amax = t)
    end
    centroid /= W
    cos = nw > 0 ? dotwt/(sqrt(nw)*tnorm) : 0f0
    spread = zero(Float32); ent = zero(Float32)
    @inbounds for t in 1:L
        j = t - (k+1)
        spread += (Float32(j)-centroid)^2 * w[t]
        p = w[t]/W; p > 0 && (ent -= p*log(p))
    end
    spread = sqrt(spread/W)
    # monotonic decay outward from center (both sides)
    mok = 0; mtot = 0
    @inbounds for j in 1:k
        mtot += 2
        w[k+1+j] <= w[k+j]     && (mok += 1)   # right
        w[k+1-j] <= w[k+2-j]   && (mok += 1)   # left
    end
    mono = mtot > 0 ? Float32(mok)/Float32(mtot) : 0f0
    return (cos, w0/W, tail/W, Float32(amax-(k+1)), centroid, spread, ent, 1f0 - abs(sl-sr)/W, mono)
end

"""
    collapse_to_metascans(psms, spectra, precursors, k) -> DataFrame

Reduce the raw per-(precursor,scan) post-deconv table to one meta-PSM per
precursor meta-scan. A row is a meta-scan CENTER iff the precursor m/z lies in
that scan's ±(isolationWidth/2) window; the meta-scan is that center ±k bins
(consecutive scan indices within a cycle). Keeps center rows, adds `metascan_w_*`
(2k+1 raw weights) and the ZT_PROFILE_FEATURES shape columns.
"""
# `_zt_std`: population-free std used for the within-metascan apex dispersion
# (avoids a Statistics dependency in the hot collapse loop).
@inline function _zt_std(a::AbstractVector{Float32})
    m = length(a)
    m < 2 && return 0f0
    mu = 0f0; @inbounds for x in a; mu += x; end; mu /= m
    s = 0f0; @inbounds for x in a; s += (x - mu) * (x - mu); end
    return sqrt(s / (m - 1))
end

# Gather `col[perm]` into a fresh concrete Vector{Float32} in one pass (no intermediate
# Float32.(col) copy). Used to sort only the columns the collapse loop reads.
@inline function _permute_f32(col, perm::Vector{Int})
    out = Vector{Float32}(undef, length(perm))
    @inbounds for i in eachindex(perm); out[i] = Float32(col[perm[i]]); end
    return out
end

function collapse_to_metascans(psms::DataFrame, spectra::MassSpecData, precursors, k::Int;
                               bitvec_rank_table = nothing,
                               own_scans::Union{Nothing,Set{UInt32}} = nothing)
    n = nrow(psms)
    (n == 0 || k <= 0) && return psms
    # Concrete Float32 vectors (function barrier): the accessors return an
    # Arrow.Primitive (getMz) and a Vector{Union{Missing,Float32}} (getCenterMzs on
    # FilteredMassSpecData), whose element type boxes on every index — and these are
    # indexed 35M times in the per-row center guard below. Materialize once; MS1 scans
    # (missing center m/z) coalesce to NaN32 and are never indexed by the MS2 PSMs.
    prec_mz::Vector{Float32} = Vector{Float32}(getMz(precursors))
    cmzs::Vector{Float32} = Float32.(coalesce.(getCenterMzs(spectra), NaN32))
    hws::Vector{Float32}  = Float32.(coalesce.(getIsolationWidthMzs(spectra), NaN32))

    # sort by (precursor_idx, scan_idx): meta-scan bins become contiguous rows
    _diag = haskey(ENV, "PIONEER_DIAG_ALLOC")
    _b0 = Base.gc_bytes()
    pid0 = psms[!, :precursor_idx]::Vector{UInt32}; scn0 = psms[!, :scan_idx]::Vector{UInt32}
    perm = sortperm(1:n; by = i -> (pid0[i], scn0[i]), alg = QuickSort)
    # Sort only the columns the loop reads (pid/scn/weight/8 frags), NOT the full ~60-col
    # table; meta rows are pulled from the original `psms` via perm[center_rows] at the end.
    pid = pid0[perm]; scn = scn0[perm]
    wt  = _permute_f32(psms[!, :weight], perm)
    _diag && @user_print string("[ALLOC]   collapse.sort_cols: ", round((Base.gc_bytes()-_b0)/1e9, digits=3), " GB")
    _b0 = Base.gc_bytes()
    L = 2k + 1
    tri = Float32[max(0f0, 1f0 - abs(Float32(j))/Float32(k+1)) for j in -k:k]
    tnorm = sqrt(sum(x->x*x, tri))

    # Fragment intensity columns for the within-metascan shape features (guarded —
    # absent only in degenerate tables). Mirror the develop frag-corr aggregation.
    _have_frags = all(r -> hasproperty(psms, Symbol("frag$(r)_int")), 1:8)
    # Concrete NTuple{8,Vector{Float32}} (never Union{Nothing,…}): a Union here boxes
    # fcols[b] on every one of the ~680M inner-gather iterations. Gate use with the
    # `_have_frags` Bool instead of a `fcols !== nothing` check.
    fcols::NTuple{8,Vector{Float32}} = _have_frags ?
        ntuple(r -> _permute_f32(psms[!, Symbol("frag$(r)_int")], perm), 8) :
        ntuple(_ -> Float32[], 8)
    rank_weights = _fragment_rank_weights(8)

    center_rows = Int[]
    profiles = Vector{Vector{Float32}}()
    sh_str = Float32[]; sh_effn = Float32[]; sh_best = Float32[]
    sh_disp = Float32[]; sh_n70 = UInt8[]; sh_rank = UInt16[]
    Fbuf = [zeros(Float32, L) for _ in 1:8]
    cfw = zeros(Float32, 8)
    has_sig = falses(8)
    apx = Float32[]

    i = 1
    @inbounds while i <= n
        j0 = i
        while i <= n && pid[i] == pid[j0]; i += 1; end
        blk_lo, blk_hi = j0, i - 1
        pm = Float32(prec_mz[pid[j0]])
        for r in blk_lo:blk_hi
            (abs(pm - Float32(cmzs[scn[r]])) <= Float32(hws[scn[r]])/2) || continue
            # Batched collapse: emit a meta-PSM only if this batch OWNS the center scan.
            # Rows at scans outside own_scans are borrowed boundary neighbors — read (for
            # the ±k weight/frag gather below) but never emitted as centers. nothing = own all.
            own_scans === nothing || (scn[r] in own_scans) || continue
            c = Int(scn[r])
            w = zeros(Float32, L)
            if _have_frags
                for b in 1:8; fill!(Fbuf[b], 0f0); end
            end
            # neighbors: bins are ~contiguous rows around r (sorted by scan);
            # scan a small window and bin by scan-index offset.
            lo = max(blk_lo, r - L); hi = min(blk_hi, r + L)
            for rr in lo:hi
                d = Int(scn[rr]) - c
                if -k <= d <= k
                    w[d+k+1] += wt[rr]
                    if _have_frags
                        @inbounds for b in 1:8; Fbuf[b][d+k+1] += fcols[b][rr]; end
                    end
                end
            end
            push!(center_rows, r); push!(profiles, w)

            # ---- within-metascan shape features (fragment profile vs weight) ----
            str = 0f0; effn = 0f0; best = 0f0; disp = 0f0; n70 = UInt8(0); rnk = UInt16(0)
            if _have_frags
                mask = UInt16(0); empty!(apx)
                for b in 1:8
                    fp = Fbuf[b]
                    mx = 0f0; @inbounds for v in fp; v > mx && (mx = v); end
                    has_sig[b] = mx > 0f0
                    if has_sig[b]
                        cfw[b] = _frag_pcor(fp, w)
                        if cfw[b] > 0.7f0; n70 += UInt8(1); mask |= UInt16(1) << (b - 1); end
                        ai = 1; vm = fp[1]; @inbounds for t in 2:L; fp[t] > vm && (vm = fp[t]; ai = t); end
                        push!(apx, Float32(ai - (k + 1)))
                    else
                        cfw[b] = 0f0
                    end
                end
                str, effn = _positive_corr_summary(cfw, rank_weights)
                # best-consensus fragment (mean corr to the other signal fragments)
                best_r = 0; best_cons = typemin(Float32)
                for b in 1:8
                    has_sig[b] || continue
                    cons = 0f0; np = 0
                    for b2 in 1:8
                        (b2 == b || !has_sig[b2]) && continue
                        cons += _frag_pcor(Fbuf[b], Fbuf[b2]); np += 1
                    end
                    avg = np > 0 ? cons / np : typemin(Float32)
                    if avg > best_cons; best_cons = avg; best_r = b; end
                end
                best = best_r > 0 ? cfw[best_r] : 0f0
                disp = _zt_std(apx)
                rnk = bitvec_rank_table === nothing ? UInt16(0) :
                      _bitvec_pattern_rank(bitvec_rank_table, mask)
            end
            push!(sh_str, str); push!(sh_effn, effn); push!(sh_best, best)
            push!(sh_disp, disp); push!(sh_n70, n70); push!(sh_rank, rnk)
        end
    end
    _diag && @user_print string("[ALLOC]   collapse.main_loop: ", round((Base.gc_bytes()-_b0)/1e9, digits=3), " GB")
    _b0 = Base.gc_bytes()

    # center_rows are positions in the PERMED order; map back to original psms rows.
    meta = psms[perm[center_rows], :]
    nm = length(center_rows)
    _diag && @user_print string("[ALLOC]   collapse.center_subset: ", round((Base.gc_bytes()-_b0)/1e9, digits=3), " GB")
    _b0 = Base.gc_bytes()
    # raw meta-scan weight columns
    for (jj, j) in enumerate(-k:k)
        col = Symbol("metascan_w_", j < 0 ? "m$(abs(j))" : "$(j)")
        meta[!, col] = Float32[profiles[m][jj] for m in 1:nm]
    end
    # weight-profile shape features
    feats = [_zt_profile_features(profiles[m], k, tri, tnorm) for m in 1:nm]
    for (fi, fname) in enumerate(ZT_PROFILE_FEATURES)
        meta[!, fname] = Float32[feats[m][fi] for m in 1:nm]
    end
    # within-metascan fragment shape features (new; additive — existing frag_corr_*
    # columns kept unchanged, still computed on the pre-collapse trace for now)
    meta[!, :frag_corr_strength_shape]              = sh_str
    meta[!, :frag_corr_effective_n_shape]           = sh_effn
    meta[!, :frag_corr_best_shape]                  = sh_best
    meta[!, :frag_apex_dispersion_shape]            = sh_disp
    meta[!, :n_correlated_fragments_shape]          = sh_n70
    meta[!, :n_correlated_fragments_bitvec_rank_shape] = sh_rank
    # zt_tri_pcor (promoted): mean-centered Pearson of the weight profile vs the triangle.
    # Computed from the stored profiles so the hot per-meta-PSM loop is untouched.
    meta[!, :zt_tri_pcor] = Float32[_frag_pcor(profiles[m], tri) for m in 1:nm]
    # Experimental shape candidates under test (env-gated; additive — absent by default).
    if get(ENV, "PIONEER_ZT_SHAPE_EXP", "0") != "0"
        # (next candidate columns go here; keep in sync with ZT_SHAPE_EXP_FEATURES)
    end
    _diag && @user_print string("[ALLOC]   collapse.materialize_cols: ", round((Base.gc_bytes()-_b0)/1e9, digits=3), " GB")
    return meta
end

# ============================================================================
# 5b: Batched / streaming meta-scan collapse (deconvolve-once + borrow).
#
# Instead of vcat-ing every thread's raw per-(precursor,scan) PSMs into one 35M-row
# table and collapsing globally, deconvolve CONTIGUOUS scan batches (each scan exactly
# once, by its owning batch), then collapse each batch — BORROWING the ±k boundary rows
# from the adjacent batches' already-deconvolved tables (never recomputed). The vcat then
# joins ~3M meta-PSMs instead of ~35M raw rows. Byte-identical to the global collapse:
# per-scan deconv is independent, and the final sort by (precursor_idx, scan_idx) restores
# the exact global row order the LightGBM saw. Gated on PIONEER_ZT_METASCAN_K>0 (ZT main).
# ============================================================================

# Contiguous ZT deconv partition: sort the valid MS2 scan indices and split into
# `n_batches` non-overlapping contiguous chunks (one per thread). Non-overlapping ⇒ each
# scan is deconvolved by exactly one batch. Returns Vector of (batch_id, scan_indices).
function partition_scans_contiguous_zt(scan_idxs::AbstractVector{<:Integer}, n_batches::Int)
    s = sort!(Int.(collect(scan_idxs)))
    n = length(s)
    n == 0 && return Tuple{Int,Vector{Int}}[]
    n_batches = max(1, min(n_batches, n))
    per = cld(n, n_batches)
    tasks = Tuple{Int,Vector{Int}}[]
    b = 1; lo = 1
    while lo <= n
        hi = min(n, lo + per - 1)
        push!(tasks, (b, s[lo:hi]))
        b += 1; lo = hi + 1
    end
    return tasks
end

# Rows of `raw` whose :scan_idx is >= thr (dir=:ge, head slice) or <= thr (dir=:le, tail
# slice), as a fresh DataFrame with the same columns as `raw` (so it vcats with siblings).
# These are the ±k boundary neighbors BORROWED from an adjacent, already-deconvolved batch.
function _zt_boundary_rows(raw::DataFrame, dir::Symbol, thr::Int)
    scn = raw[!, :scan_idx]::Vector{UInt32}
    t = UInt32(max(0, thr))
    keep = dir === :ge ? findall(>=(t), scn) : findall(<=(t), scn)
    return raw[keep, :]
end

# batch_tbl = head ++ raw ++ tail (borrowed slices may be nothing at the file ends).
# The three scan ranges are disjoint (contiguous partition), so no neighbor is double-counted.
function _zt_concat_borrow(head, raw::DataFrame, tail)
    (head === nothing && tail === nothing) && return raw
    parts = DataFrame[]
    head === nothing || push!(parts, head)
    push!(parts, raw)
    tail === nothing || push!(parts, tail)
    return vcat(parts...)
end

# Unwrap a TaskFailedException so the real deconv/collapse error surfaces (mirrors
# the pattern in library_search's fetch).
function _zt_fetch(t)
    try
        return fetch(t)
    catch e
        while e isa TaskFailedException
            e = e.task.exception
        end
        rethrow(e)
    end
end

"""
    zt_batched_deconv_collapse(zt_tasks, nce_entries, spectra, prec_index, ms_file_idx,
                               search_context, params, precursors, ion_list,
                               qtm_deconv, mem, rt_to_irt, irt_tol, k) -> DataFrame

ZT main-search replacement for the `vcat(fetched...)` raw table + downstream global
collapse. For each NCE model:
  Phase 1 (parallel, per batch): deconvolve the batch's own scans ONCE, then run the
    per-scan-local `add_scan_competition_features!` and per-row `add_ms1_lookup_features!`
    (both batch-safe — a scan lives entirely in one batch).
  Phase 2 (parallel, per batch): collapse, borrowing the ±k boundary rows from the
    neighbors' Phase-1 tables; emit centers only for scans this batch owns.
Then vcat the per-batch meta tables and sort by (precursor_idx, scan_idx) to reproduce the
global order. Returns the collapsed meta-PSM DataFrame (one row per (precursor, cycle)).
"""
function zt_batched_deconv_collapse(zt_tasks, nce_entries, spectra, prec_index,
                                    ms_file_idx::Int64, search_context, params,
                                    precursors, ion_list, qtm_deconv, mem, rt_to_irt,
                                    irt_tol, k::Int)
    search_data       = getSearchData(search_context)
    lib_precs         = getPrecursors(getSpecLib(search_context))
    bitvec_rank_table = getBitVecExcessRanks(search_context, Int64(ms_file_idx))
    nb = length(zt_tasks)

    metas = map(nce_entries) do (nce_model, nce_tag)
        # --- Phase 1: deconvolve each contiguous batch ONCE + batch-safe feature passes ---
        raw_tasks = map(zt_tasks) do (bid, own_scan_idxs)
            Threads.@spawn begin
                raw = process_scans_fused!(own_scan_idxs, spectra, prec_index, ms_file_idx,
                        search_data[bid], params, precursors, ion_list,
                        nce_model, qtm_deconv, mem, rt_to_irt, irt_tol)
                add_scan_competition_features!(raw)                               # per-scan-local
                add_ms1_lookup_features!(raw, spectra, search_context, ms_file_idx) # per-row
                raw
            end
        end
        raws = map(_zt_fetch, raw_tasks)

        # --- Phase 2: collapse each batch, borrowing ±k boundary rows from neighbors ---
        meta_tasks = map(enumerate(zt_tasks)) do (t, (bid, own_scan_idxs))
            Threads.@spawn begin
                lo_s = first(own_scan_idxs); hi_s = last(own_scan_idxs)  # sorted within batch
                head = t > 1  ? _zt_boundary_rows(raws[t-1], :ge, lo_s - k) : nothing
                tail = t < nb ? _zt_boundary_rows(raws[t+1], :le, hi_s + k) : nothing
                batch_tbl = _zt_concat_borrow(head, raws[t], tail)
                own = Set{UInt32}(UInt32.(own_scan_idxs))
                collapse_to_metascans(batch_tbl, spectra, lib_precs, k;
                        bitvec_rank_table = bitvec_rank_table, own_scans = own)
            end
        end
        # Keep only batches that produced meta-PSMs. A batch with rows always yields the
        # full meta schema (even with 0 centers); the only 0-row/raw-schema case is a fully
        # empty batch (collapse early-returns its input), which we drop so vcat sees uniform
        # meta-schema tables. Real contiguous ranges are non-empty.
        parts_b = filter(m -> nrow(m) > 0, map(_zt_fetch, meta_tasks))
        m = isempty(parts_b) ? DataFrame() : vcat(parts_b...)

        # Decisive isolation diagnostic (PIONEER_ZT_DIAG_COMPARE=1): run the GLOBAL collapse
        # on the SAME concatenated batched-raw. If it matches the baseline center count, the
        # deconv rows are identical and any discrepancy is in the batched collapse/borrow
        # logic; if it also differs, the raw rows themselves differ (partition-dependent deconv).
        if get(ENV, "PIONEER_ZT_DIAG_COMPARE", "0") != "0"
            raw_all = vcat(raws...)
            mg = collapse_to_metascans(raw_all, spectra, lib_precs, k;
                                       bitvec_rank_table = bitvec_rank_table)
            @user_print string("[ZT DIAG] raw_rows=", nrow(raw_all),
                "  global_collapse_on_batched_raw=", nrow(mg),
                "  batched_collapse=", nrow(m))
        end

        nce_tag === nothing || isempty(m) || (m[!, :nce] .= nce_tag)
        m
    end

    result = length(metas) == 1 ? metas[1] : vcat(metas...)
    isempty(result) || sort!(result, [:precursor_idx, :scan_idx])   # restore global order
    return result
end

# ============================================================================
# 5b (chosen): DISK-SPILLED collapse. Keeps develop's round-robin deconv EXACTLY
# (byte-identical raw rows, cache-optimal), but instead of vcat-ing every thread's raw
# PSMs into one 34M-row table, scatters each thread table to on-disk contiguous scan-range
# bins (freeing the thread's scored-PSM scratch after), then collapses one bin at a time via
# a 3-bin sliding window with ±k halo borrow. The full raw never coexists in memory.
# Byte-identical: raw rows unchanged; per-thread scan_competition == global (scans are
# thread-disjoint); collapse-per-bin-with-borrow proven == global collapse (DIAG_COMPARE);
# final sort restores (pid,scn) order. Gated PIONEER_ZT_SPILL=1.
# ============================================================================

# Scatter one thread's raw table into per-(bin,thread) Arrow files. `binof(scan)->1..B`.
function _zt_scatter_to_bins!(spill_dir::String, thread_id::Int, tbl::DataFrame, binof, B::Int)
    scn = tbl[!, :scan_idx]::Vector{UInt32}
    bins = Vector{Int}(undef, length(scn))
    @inbounds for i in eachindex(scn); bins[i] = binof(Int(scn[i])); end
    for b in 1:B
        rows = findall(==(b), bins)
        isempty(rows) && continue
        writeArrow(joinpath(spill_dir, "bin$(b)_thr$(thread_id).arrow"), tbl[rows, :])
    end
    return nothing
end

"""
    zt_disk_spill_deconv_collapse(...) -> DataFrame

ZT main-search replacement for `vcat(fetched...)` + downstream global collapse that keeps peak
memory low. Deconvolves with develop's round-robin `thread_tasks` (UNCHANGED), runs per-thread
`scan_competition`/`ms1_lookup` (byte-identical), scatters each thread table to `n_bins`
on-disk contiguous scan-range bins and frees that thread's scored scratch, then collapses bins
one at a time (3-bin sliding window, ±k halo). Returns the (pid,scn)-sorted meta table.
"""
function zt_disk_spill_deconv_collapse(thread_tasks, nce_entries, spectra, prec_index,
                                       ms_file_idx::Int64, search_context, params, precursors,
                                       ion_list, qtm_deconv, mem, rt_to_irt, irt_tol, k::Int,
                                       all_scan_idxs::Vector{Int}, n_bins::Int)
    search_data       = getSearchData(search_context)
    lib_precs         = getPrecursors(getSpecLib(search_context))
    bitvec_rank_table = getBitVecExcessRanks(search_context, Int64(ms_file_idx))
    nthreads_tasks    = length(thread_tasks)

    spill_dir = joinpath(getDataOutDir(search_context), "temp_data", "zt_spill")
    isdir(spill_dir) && rm(spill_dir; recursive=true, force=true)
    mkpath(spill_dir)

    # B contiguous scan-range bins over the sorted MS2 scans.
    sorted_scans = sort(all_scan_idxs)
    n = length(sorted_scans)
    B = max(1, min(n_bins, n))
    per = cld(n, B)
    bin_scans = Vector{Vector{Int}}(undef, B)
    bin_lo = Vector{Int}(undef, B); bin_hi = Vector{Int}(undef, B)
    for b in 1:B
        lo_i = (b-1)*per + 1; hi_i = min(b*per, n)
        bin_scans[b] = sorted_scans[lo_i:hi_i]
        bin_lo[b] = sorted_scans[lo_i]; bin_hi[b] = sorted_scans[hi_i]
    end
    binof(s::Int) = clamp(searchsortedlast(bin_lo, s), 1, B)   # bins disjoint & contiguous

    metas = map(nce_entries) do (nce_model, nce_tag)
        # --- Deconv (round-robin, UNCHANGED) then per-thread features + scatter + free scratch ---
        tasks = map(thread_tasks) do thread_task
            Threads.@spawn process_scans_fused!(
                last(thread_task), spectra, prec_index, ms_file_idx,
                search_data[first(thread_task)], params, precursors, ion_list,
                nce_model, qtm_deconv, mem, rt_to_irt, irt_tol)
        end
        for (idx, t) in enumerate(tasks)
            tid = first(thread_tasks[idx])
            tbl = _zt_fetch(t)
            # process_scans_fused! now returns a StructArray-backed DataFrame with VIEW columns
            # (copycols=false). Materialize this thread's rows to concrete Vector columns so the
            # typed per-thread passes (scan_competition/ms1) and the scatter work, and so the
            # scratch can be freed immediately (the copy is independent of the buffer).
            nrow(tbl) > 0 && (tbl = DataFrame(tbl; copycols = true))
            if nrow(tbl) > 0
                add_scan_competition_features!(tbl)                                # per-thread
                add_ms1_lookup_features!(tbl, spectra, search_context, ms_file_idx) # per-row
                _zt_scatter_to_bins!(spill_dir, tid, tbl, binof, B)
            end
            tbl = nothing
            # Release this thread's ~0.6 GB raw scored buffer so it BECOMES reclaimable.
            # NOTE: resize!/empty! KEEP the Vector's capacity (no memory freed); reassigning
            # drops the reference so GC can collect it. We do NOT force GC — Julia reclaims it
            # when there's memory pressure (promptly on a constrained machine, lazily on a big
            # one), so the program's *minimum required* memory drops without paying a GC tax
            # here. Regrown next file by growScoredPSMs!.
            search_data[tid].main_search_scored_psms = StructArray(MainSearchScoredPSM{Float32,Float16}[])
        end

        # --- Collapse bins one at a time (3-bin sliding window; borrow ±k halo from neighbors) ---
        function read_bin(b)
            paths = String[]
            for ti in 1:nthreads_tasks
                p = joinpath(spill_dir, "bin$(b)_thr$(ti).arrow")
                isfile(p) && push!(paths, p)
            end
            isempty(paths) && return nothing
            DataFrame(Tables.columntable(Arrow.Table(paths)))   # materialize (detach from mmap)
        end
        out = DataFrame[]
        prev = nothing
        cur  = read_bin(1)
        for b in 1:B
            nxt = b < B ? read_bin(b+1) : nothing
            if cur !== nothing && nrow(cur) > 0
                head = (prev !== nothing && nrow(prev) > 0) ?
                    prev[prev[!, :scan_idx] .>= UInt32(max(0, bin_lo[b]-k)), :] : nothing
                tail = (nxt !== nothing && nrow(nxt) > 0) ?
                    nxt[nxt[!, :scan_idx] .<= UInt32(bin_hi[b]+k), :] : nothing
                batch = _zt_concat_borrow(head, cur, tail)
                own = Set{UInt32}(UInt32.(bin_scans[b]))
                mb = collapse_to_metascans(batch, spectra, lib_precs, k;
                        bitvec_rank_table = bitvec_rank_table, own_scans = own)
                nrow(mb) > 0 && push!(out, mb)
            end
            prev = cur; cur = nxt
        end
        m = isempty(out) ? DataFrame() : vcat(out...)
        nce_tag === nothing || isempty(m) || (m[!, :nce] .= nce_tag)
        m
    end

    rm(spill_dir; recursive=true, force=true)
    result = length(metas) == 1 ? metas[1] : vcat(filter(x -> nrow(x) > 0, metas)...)
    isempty(result) || sort!(result, [:precursor_idx, :scan_idx])
    return result
end
