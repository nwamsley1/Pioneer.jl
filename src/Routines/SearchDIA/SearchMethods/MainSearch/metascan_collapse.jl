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
function collapse_to_metascans(psms::DataFrame, spectra::MassSpecData, precursors, k::Int)
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

    center_rows = Int[]
    profiles = Vector{Vector{Float32}}()

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
            # neighbors: bins are ~contiguous rows around r (sorted by scan);
            # scan a small window and bin by scan-index offset.
            lo = max(blk_lo, r - L); hi = min(blk_hi, r + L)
            for rr in lo:hi
                d = Int(scn[rr]) - c
                (-k <= d <= k) && (w[d+k+1] += wt[rr])
            end
            push!(center_rows, r); push!(profiles, w)
        end
    end

    meta = psms[center_rows, :]
    nm = length(center_rows)
    # raw meta-scan weight columns
    for (jj, j) in enumerate(-k:k)
        col = Symbol("metascan_w_", j < 0 ? "m$(abs(j))" : "$(j)")
        meta[!, col] = Float32[profiles[m][jj] for m in 1:nm]
    end
    # shape features
    feats = [_zt_profile_features(profiles[m], k, tri, tnorm) for m in 1:nm]
    for (fi, fname) in enumerate(ZT_PROFILE_FEATURES)
        meta[!, fname] = Float32[feats[m][fi] for m in 1:nm]
    end
    return meta
end
