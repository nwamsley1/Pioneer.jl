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
]

# Across-cycle ELUTION features, computed post-collapse on the one-point-per-cycle
# meta trace: the correctly-leveled counterparts of the develop MS1/fragment
# chromatogram features. The MS1 four are the FIX (develop's are broken — MS1 is one
# scan per cycle, so constant within a cycle → correlated vs the within-cycle weight
# triangle = noise). n_scans_metascan = number of cycles (real peak width vs the
# 13×cycles `n_scans` artifact).
const ZT_ELUTION_FEATURES = Symbol[
    :ms1_corr_weight_m0_elution, :ms1_corr_m0_m1_elution,
    :ms1_apex_offset_irt_elution, :ms1_weight_apex_to_m0_apex_irt_elution,
    :frag_corr_strength_elution, :frag_corr_effective_n_elution,
    :frag_corr_best_elution, :frag_apex_dispersion_elution,
    :n_correlated_fragments_elution, :n_correlated_fragments_bitvec_rank_elution,
    :n_scans_metascan,
]

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

function collapse_to_metascans(psms::DataFrame, spectra::MassSpecData, precursors, k::Int;
                               bitvec_rank_table = nothing)
    n = nrow(psms)
    (n == 0 || k <= 0) && return psms
    prec_mz = getMz(precursors)
    cmzs = getCenterMzs(spectra)
    hws  = getIsolationWidthMzs(spectra)

    # sort by (precursor_idx, scan_idx): meta-scan bins become contiguous rows
    pid0 = psms[!, :precursor_idx]; scn0 = psms[!, :scan_idx]
    perm = sortperm(1:n; by = i -> (pid0[i], scn0[i]), alg = QuickSort)
    psms = psms[perm, :]
    pid = psms[!, :precursor_idx]::Vector{UInt32}
    scn = psms[!, :scan_idx]::Vector{UInt32}
    wt  = Float32.(psms[!, :weight])
    L = 2k + 1
    tri = Float32[max(0f0, 1f0 - abs(Float32(j))/Float32(k+1)) for j in -k:k]
    tnorm = sqrt(sum(x->x*x, tri))

    # Fragment intensity columns for the within-metascan shape features (guarded —
    # absent only in degenerate tables). Mirror the develop frag-corr aggregation.
    _have_frags = all(r -> hasproperty(psms, Symbol("frag$(r)_int")), 1:8)
    fcols = _have_frags ? ntuple(r -> Float32.(psms[!, Symbol("frag$(r)_int")]), 8) : nothing
    rank_weights = _fragment_rank_weights(8)

    center_rows = Int[]
    profiles = Vector{Vector{Float32}}()
    sh_str = Float32[]; sh_effn = Float32[]; sh_best = Float32[]
    sh_disp = Float32[]; sh_n70 = UInt8[]; sh_rank = UInt16[]
    # triangle-weighted combined per-cycle intensities (weight + 8 fragments) for the
    # across-cycle elution trace (PIONEER_ZT_ELUTION_INTENSITY=triangle option).
    mw = Float32[]; mf = [Float32[] for _ in 1:8]
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
            c = Int(scn[r])
            w = zeros(Float32, L)
            if fcols !== nothing
                for b in 1:8; fill!(Fbuf[b], 0f0); end
            end
            # neighbors: bins are ~contiguous rows around r (sorted by scan);
            # scan a small window and bin by scan-index offset.
            lo = max(blk_lo, r - L); hi = min(blk_hi, r + L)
            for rr in lo:hi
                d = Int(scn[rr]) - c
                if -k <= d <= k
                    w[d+k+1] += wt[rr]
                    if fcols !== nothing
                        @inbounds for b in 1:8; Fbuf[b][d+k+1] += fcols[b][rr]; end
                    end
                end
            end
            push!(center_rows, r); push!(profiles, w)

            # ---- within-metascan shape features (fragment profile vs weight) ----
            str = 0f0; effn = 0f0; best = 0f0; disp = 0f0; n70 = UInt8(0); rnk = UInt16(0)
            if fcols !== nothing
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

            # triangle-weighted combined per-cycle intensities
            mwi = 0f0; @inbounds for t in 1:L; mwi += w[t] * tri[t]; end
            push!(mw, mwi)
            if fcols !== nothing
                for b in 1:8
                    mfi = 0f0; @inbounds for t in 1:L; mfi += Fbuf[b][t] * tri[t]; end
                    push!(mf[b], mfi)
                end
            else
                for b in 1:8; push!(mf[b], 0f0); end
            end
        end
    end

    meta = psms[center_rows, :]
    nm = length(center_rows)
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
    # triangle-weighted combined per-cycle intensities (consumed by the post-collapse
    # elution feature pass; scratch, dropped before writing the second-pass table)
    meta[!, :zt_meta_weight] = mw
    for b in 1:8
        meta[!, Symbol("zt_meta_frag", b)] = mf[b]
    end
    return meta
end

# Post-collapse across-cycle ELUTION features on the collapsed meta-PSMs (one row per
# (precursor, cycle), already sorted by (precursor_idx, scan_idx) ≈ cycle order).
# Mirrors the develop MS1/fragment chromatogram aggregation but on the clean
# one-point-per-cycle trace. `use_triangle` picks the per-cycle intensity: the
# triangle-weighted combined value (zt_meta_*) or the center-bin value
# (metascan_w_0 / frag*_int). Per-precursor features are replicated across the
# precursor's cycle rows; best-per-precursor selects one downstream.
function add_zt_elution_features!(meta::DataFrame; bitvec_rank_table = nothing,
                                  use_triangle::Bool = true)
    n = nrow(meta)
    n == 0 && return meta
    pid = meta[!, :precursor_idx]::Vector{UInt32}
    irt = hasproperty(meta, :irt_obs) ? Float32.(meta[!, :irt_obs]) : zeros(Float32, n)
    wcol = use_triangle ? Float32.(meta[!, :zt_meta_weight]) : Float32.(meta[!, :metascan_w_0])
    m0col = hasproperty(meta, :ms1_m0_intensity) ? Float32.(meta[!, :ms1_m0_intensity]) : zeros(Float32, n)
    m1col = hasproperty(meta, :ms1_m1_intensity) ? Float32.(meta[!, :ms1_m1_intensity]) : zeros(Float32, n)
    fcols = ntuple(r -> use_triangle ? Float32.(meta[!, Symbol("zt_meta_frag", r)]) :
                                       Float32.(meta[!, Symbol("frag", r, "_int")]), 8)
    rank_weights = _fragment_rank_weights(8)

    out_wm0 = zeros(Float32, n); out_m01 = zeros(Float32, n)
    out_aoff = zeros(Float32, n); out_w2m0 = zeros(Float32, n)
    out_fstr = zeros(Float32, n); out_feffn = zeros(Float32, n); out_fbest = zeros(Float32, n)
    out_fdisp = zeros(Float32, n); out_fn70 = zeros(UInt8, n); out_frank = zeros(UInt16, n)
    out_ncyc = ones(UInt32, n)

    Wg = Float32[]; M0g = Float32[]; M1g = Float32[]; IRTg = Float32[]
    Fg = [Float32[] for _ in 1:8]
    cfw = zeros(Float32, 8); has_sig = falses(8); apxg = Float32[]

    i = 1
    @inbounds while i <= n
        j0 = i
        while i <= n && pid[i] == pid[j0]; i += 1; end
        lo = j0; hi = i - 1; m = hi - lo + 1
        for c in lo:hi; out_ncyc[c] = UInt32(m); end
        m < 2 && continue
        resize!(Wg, m); resize!(M0g, m); resize!(M1g, m); resize!(IRTg, m)
        for b in 1:8; resize!(Fg[b], m); end
        for (kk, c) in enumerate(lo:hi)
            Wg[kk] = wcol[c]; M0g[kk] = m0col[c]; M1g[kk] = m1col[c]; IRTg[kk] = irt[c]
            for b in 1:8; Fg[b][kk] = fcols[b][c]; end
        end
        c_wm0 = _frag_pcor(Wg, M0g); c_m01 = _frag_pcor(M0g, M1g)
        ai_m0 = argmax(M0g); ai_w = argmax(Wg)
        w2m0 = abs(IRTg[ai_w] - IRTg[ai_m0])
        empty!(apxg); mask = UInt16(0)
        for b in 1:8
            mx = 0f0; for v in Fg[b]; v > mx && (mx = v); end
            has_sig[b] = mx > 0f0
            if has_sig[b]
                cfw[b] = _frag_pcor(Fg[b], Wg)
                cfw[b] > 0.7f0 && (mask |= UInt16(1) << (b - 1))
                push!(apxg, IRTg[argmax(Fg[b])])
            else
                cfw[b] = 0f0
            end
        end
        n70 = UInt8(0); for b in 1:8; (has_sig[b] && cfw[b] > 0.7f0) && (n70 += UInt8(1)); end
        strv, effnv = _positive_corr_summary(cfw, rank_weights)
        best_r = 0; best_cons = typemin(Float32)
        for b in 1:8
            has_sig[b] || continue
            cons = 0f0; np = 0
            for b2 in 1:8
                (b2 == b || !has_sig[b2]) && continue
                cons += _frag_pcor(Fg[b], Fg[b2]); np += 1
            end
            avg = np > 0 ? cons / np : typemin(Float32)
            if avg > best_cons; best_cons = avg; best_r = b; end
        end
        c_best = 0f0
        if best_r > 0
            mxm0 = 0f0; for v in M0g; v > mxm0 && (mxm0 = v); end
            mxm0 > 0f0 && (c_best = _frag_pcor(Fg[best_r], M0g))
        end
        disp = length(apxg) >= 2 ? _zt_std(apxg) : 0f0
        rnk = bitvec_rank_table === nothing ? UInt16(0) : _bitvec_pattern_rank(bitvec_rank_table, mask)
        for (kk, c) in enumerate(lo:hi)
            out_wm0[c] = c_wm0; out_m01[c] = c_m01
            out_aoff[c] = abs(IRTg[kk] - IRTg[ai_m0]); out_w2m0[c] = w2m0
            out_fstr[c] = strv; out_feffn[c] = effnv; out_fbest[c] = c_best
            out_fdisp[c] = disp; out_fn70[c] = n70; out_frank[c] = rnk
        end
    end

    meta[!, :ms1_corr_weight_m0_elution]              = out_wm0
    meta[!, :ms1_corr_m0_m1_elution]                  = out_m01
    meta[!, :ms1_apex_offset_irt_elution]             = out_aoff
    meta[!, :ms1_weight_apex_to_m0_apex_irt_elution]  = out_w2m0
    meta[!, :frag_corr_strength_elution]              = out_fstr
    meta[!, :frag_corr_effective_n_elution]           = out_feffn
    meta[!, :frag_corr_best_elution]                  = out_fbest
    meta[!, :frag_apex_dispersion_elution]            = out_fdisp
    meta[!, :n_correlated_fragments_elution]          = out_fn70
    meta[!, :n_correlated_fragments_bitvec_rank_elution] = out_frank
    meta[!, :n_scans_metascan]                        = out_ncyc
    # drop the per-cycle combined-intensity scratch (consumed above)
    select!(meta, Not(vcat([:zt_meta_weight], [Symbol("zt_meta_frag", b) for b in 1:8])))
    return meta
end
