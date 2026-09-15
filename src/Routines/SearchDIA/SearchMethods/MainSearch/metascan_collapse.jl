# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# ============================================================================
# Scanning-quad (ZT) meta-scan collapse.
#
# The swept quadrupole spreads a precursor's ions across ~2k+1 adjacent Q1 bins, so after
# deconvolution the same precursor holds a separate weight in every bin of its meta-scan.
# This reduces that raw per-(precursor, scan) table to one meta-PSM per meta-scan, carrying
# features derived from the shape of the weight profile across the bins.
# ============================================================================

"""
Weight-profile features. Of the nine descriptors originally evaluated only these two are
consumed by any model: `zt_tri_cosine` (2nd-pass gain 56k) and `zt_entropy` (1.3k). The other
seven — center_frac, tail_frac, apex_offset, centroid, spread, symmetry, monotonicity — all
measured below 400 and are not computed.
"""
const ZT_PROFILE_FEATURES = Symbol[
    :zt_tri_cosine,   # cosine of the weight profile against an ideal triangle template
    :zt_entropy,      # Shannon entropy of the normalized profile
    :zt_fit_h_ratio,  # per-meta-scan triangle half-base (a/b) over the file's fitted h; 0 if unfit
    :zt_fit_r2,       # R² of that per-meta-scan linear fit; 0 if unfit
]

"""
    _zt_local_triangle_fit(w, k, bin_step, delta0, h_file; lim) -> (h_ratio, r2)

The SAME estimator QuadTuningSearch uses per meta-scan (`w = a - b*|Δ|`, `h = a/b`), applied
to one collapsed profile: bins with weight > 0 and |Δ| <= `lim` Da, where Δ is the bin's offset
from the precursor m/z (`delta0` = precursor m/z minus the anchor bin's centre). Returns the
local half-base relative to the file's fitted `h_file` (clamped to [0, 4]) and the fit R².
Zero when fewer than 4 usable bins or the slope has the wrong sign.
"""
@inline function _zt_local_triangle_fit(w::Vector{Float32}, k::Int, bin_step::Float32,
                                        delta0::Float32, h_file::Float32; lim::Float32 = 4f0)
    h_file > 0f0 || return (0f0, 0f0)
    n = 0; sx = 0f0; sy = 0f0; sxx = 0f0; sxy = 0f0
    @inbounds for t in 1:(2k + 1)
        y = w[t]; y > 0f0 || continue
        x = abs(Float32(t - k - 1) * bin_step - delta0)
        x <= lim || continue
        n += 1; sx += x; sy += y; sxx += x * x; sxy += x * y
    end
    n >= 4 || return (0f0, 0f0)
    den = n * sxx - sx * sx
    den > 0f0 || return (0f0, 0f0)
    b = -(n * sxy - sx * sy) / den
    a = (sy + b * sx) / n
    (a > 0f0 && b > 0f0) || return (0f0, 0f0)
    ym = sy / n; ssr = 0f0; sst = 0f0
    @inbounds for t in 1:(2k + 1)
        y = w[t]; y > 0f0 || continue
        x = abs(Float32(t - k - 1) * bin_step - delta0)
        x <= lim || continue
        pr = a - b * x; ssr += (y - pr)^2; sst += (y - ym)^2
    end
    r2 = sst > 0f0 ? max(0f0, 1f0 - ssr / sst) : 0f0
    return (clamp((a / b) / h_file, 0f0, 4f0), r2)
end

"""
Within-metascan (shape) fragment features: across the 2k+1 bins of ONE cycle, how well does
each fragment's intensity profile track the transmission triangle? Mirrors develop's
across-precursor frag-corr aggregation (`_frag_pcor` / `_positive_corr_summary` /
`_bitvec_pattern_rank` in features.jl) but on the per-cycle profile. The across-cycle elution
counterparts are develop's own chromatogram features, run post-collapse on the meta trace.
"""
const ZT_SHAPE_FEATURES = Symbol[
    :frag_corr_strength_shape,
    :frag_corr_effective_n_shape,
    :frag_corr_best_shape,
    :frag_apex_dispersion_shape,
    :n_correlated_fragments_shape,
    :n_correlated_fragments_bitvec_rank_shape,
    :zt_tri_pcor,     # mean-centered Pearson vs the template; complements the uncentered cosine
]

# zt_emp_cosine (cosine vs a template shifted to the precursor's in-bin m/z offset) was dropped
# 2026-09-10: measured r = 0.991 with zt_tri_cosine over 5.46M meta-PSMs, with no discrimination
# advantage (|AUC-0.5| 0.0213 vs 0.0211). The physical idea -- transmission peaks at the
# precursor's true m/z, not the bin centre -- is sound, but the shift is at most +/-0.5 bins on a
# profile ~6 bins wide, so the shifted and centred templates are collinear. It cost a fresh
# template allocation per meta-PSM.

"""
The subset the per-file main-search LGBM consumes. The 2nd-pass model takes all of
`ZT_PROFILE_FEATURES` + `ZT_SHAPE_FEATURES` via `ADVANCED_FEATURE_SET`.
"""
const ZT_MAINSEARCH_MODEL_FEATURES = Symbol[
    :zt_tri_cosine, :zt_tri_pcor, :frag_corr_strength_shape, :n_correlated_fragments_shape,
]

"""
    _zt_profile_features(w, k, tri, tnorm) -> (zt_tri_cosine, zt_entropy)

Weight-profile descriptors over the length-(2k+1) profile `w` (j = -k..k, center at k+1),
given the ideal triangle template `tri` and its norm.
"""
@inline function _zt_profile_features(w::Vector{Float32}, k::Int, tri::Vector{Float32}, tnorm::Float32)
    L = 2k + 1
    W = zero(Float32); dotwt = zero(Float32); nw = zero(Float32)
    @inbounds for t in 1:L
        W     += w[t]
        dotwt += w[t] * tri[t]
        nw    += w[t] * w[t]
    end
    W <= 0 && return (0f0, 0f0)
    cos = nw > 0 ? dotwt / (sqrt(nw) * tnorm) : 0f0
    ent = zero(Float32)
    @inbounds for t in 1:L
        p = w[t] / W
        p > 0 && (ent -= p * log(p))
    end
    return (cos, ent)
end

"""
    _zt_std(a) -> Float32

Sample standard deviation, used for the within-metascan fragment apex dispersion. Local so the
hot collapse loop does not reach into Statistics.
"""
@inline function _zt_std(a::AbstractVector{Float32})
    m = length(a)
    m < 2 && return 0f0
    mu = 0f0; @inbounds for x in a; mu += x; end; mu /= m
    s = 0f0;  @inbounds for x in a; s += (x - mu) * (x - mu); end
    return sqrt(s / (m - 1))
end

"""
    _permute_f32(col, perm) -> Vector{Float32}

Gather `col[perm]` into a fresh concrete `Vector{Float32}` in one pass, with no intermediate
`Float32.(col)` copy. Used to sort only the columns the collapse loop reads.
"""
@inline function _permute_f32(col, perm::Vector{Int})
    out = Vector{Float32}(undef, length(perm))
    @inbounds for i in eachindex(perm); out[i] = Float32(col[perm[i]]); end
    return out
end

"""
    _zt_dump_precollapse(psms, search_context, ms_file_idx)

TEMPORARY (Stage 2 bring-up). When `PIONEER_ZT_DUMP_PRECOLLAPSE=<dir>` is set, write the columns
`collapse_to_metascans` reads to `<dir>/precollapse_file<N>.arrow`, so the collapse can be
benchmarked and validated against real input outside a full search. Inert when unset; remove
once Stage 2 is verified.
"""
function _zt_dump_precollapse(psms::DataFrame, search_context, ms_file_idx; chunk::Int = 0)
    dir = get(ENV, "PIONEER_ZT_DUMP_PRECOLLAPSE", "")
    isempty(dir) && return nothing
    cols = Symbol[:precursor_idx, :scan_idx, :weight]
    for b in 1:8
        c = Symbol("frag$(b)_int")
        hasproperty(psms, c) && push!(cols, c)
    end
    mkpath(dir)
    path = joinpath(dir, chunk > 0 ? "precollapse_file$(ms_file_idx)_chunk$(chunk).arrow" :
                                     "precollapse_file$(ms_file_idx).arrow")
    writeArrow(path, psms[!, cols])
    @user_info "ZT: dumped $(nrow(psms)) pre-collapse rows, $(length(cols)) cols -> $path"
    return nothing
end

"""
    collapse_to_metascans(psms, spectra, precursors, k; bitvec_rank_table = nothing) -> DataFrame

Reduce the raw per-(precursor, scan) post-deconvolution table to one meta-PSM per precursor
meta-scan. A row is a meta-scan CENTER iff the precursor m/z lies inside that scan's
±(isolationWidth/2) window; the meta-scan is that center ±`k` bins (consecutive scan indices
within a cycle). Keeps the center rows and attaches `ZT_PROFILE_FEATURES` +
`ZT_SHAPE_FEATURES`.

The raw 2k+1 weight profile is NOT materialized as columns — it is consumed inline and a
single buffer is reused across centers.
"""
function collapse_to_metascans(psms::DataFrame, spectra::MassSpecData, precursors,
                               geom::ZTGeometry; bitvec_rank_table = nothing)
    k = Int(geom.metascan_k)
    n = nrow(psms)
    (n == 0 || k <= 0) && return psms

    # Function barrier: these accessors return an Arrow.Primitive and a
    # Vector{Union{Missing,Float32}} whose element types box on every index — and they are
    # indexed once per row in the center guard below. Materialize once. MS1 scans (missing
    # center m/z) coalesce to NaN32 and are never indexed by MS2 PSMs.
    prec_mz::Vector{Float32} = Vector{Float32}(getMz(precursors))
    cmzs::Vector{Float32} = Float32.(coalesce.(getCenterMzs(spectra), NaN32))
    hws::Vector{Float32}  = Float32.(coalesce.(getIsolationWidthMzs(spectra), NaN32))
    cycles::Vector{UInt32} = UInt32.(getCycleIdxs(spectra))

    # Sort by (precursor_idx, scan_idx) so a precursor's meta-scan bins become contiguous
    # rows. Both keys are UInt32, so pack them into one UInt64: identical lexicographic
    # order, a single integer compare rather than a tuple built per comparison.
    pid0 = psms[!, :precursor_idx]::Vector{UInt32}
    scn0 = psms[!, :scan_idx]::Vector{UInt32}
    sortkeys = Vector{UInt64}(undef, n)
    @inbounds for i in 1:n
        sortkeys[i] = (UInt64(pid0[i]) << 32) | UInt64(scn0[i])
    end
    perm = sortperm(sortkeys; alg = QuickSort)
    sortkeys = UInt64[]                       # release before the column gather below

    # Permute ONLY the columns the loop reads, not the ~60-column table; meta rows are pulled
    # from the original `psms` via perm[center_rows] at the end.
    pid = pid0[perm]
    scn = scn0[perm]
    wt  = _permute_f32(psms[!, :weight], perm)

    L = 2k + 1
    tri, tnorm = zt_transmission_template(geom, k)

    # Fragment intensity columns for the shape features. Concrete NTuple{8,Vector{Float32}}
    # (never a Union): a Union boxes fcols[b] on every inner-gather iteration. Gate use with
    # the `_have_frags` Bool instead of a `!== nothing` check.
    _have_frags = all(r -> hasproperty(psms, Symbol("frag$(r)_int")), 1:8)
    fcols::NTuple{8,Vector{Float32}} = _have_frags ?
        ntuple(r -> _permute_f32(psms[!, Symbol("frag$(r)_int")], perm), 8) :
        ntuple(_ -> Float32[], 8)
    rank_weights = _fragment_rank_weights(8)

    # Per-bin fitted (deconvolved) and shadow (raw observed) spectra. References only -- the 8
    # frag columns above are permuted for locality, but permuting 16 more would cost ~3.8 GB, so
    # these are indexed as fit_cols[b][perm[rr]] inside the <=13-row inner loop instead.
    _have_fs = all(r -> hasproperty(psms, Symbol("fitted_frag$(r)_int")) &&
                        hasproperty(psms, Symbol("shadow_frag$(r)_int")), 1:8)
    fit_cols::NTuple{8,Vector{Float32}} = _have_fs ?
        ntuple(r -> psms[!, Symbol("fitted_frag$(r)_int")]::Vector{Float32}, 8) :
        ntuple(_ -> Float32[], 8)
    shd_cols::NTuple{8,Vector{Float32}} = _have_fs ?
        ntuple(r -> psms[!, Symbol("shadow_frag$(r)_int")]::Vector{Float32}, 8) :
        ntuple(_ -> Float32[], 8)
    accf = zeros(Float32, 8)
    accs = zeros(Float32, 8)

    # Reusable per-center scratch — nothing here outlives one iteration.
    w       = zeros(Float32, L)
    Fbuf    = [zeros(Float32, L) for _ in 1:8]
    cfw     = zeros(Float32, 8)
    has_sig = falses(8)
    apx     = Float32[]

    # Per-center outputs (scalars only; the profile is never retained).
    hint = n ÷ L + 1
    center_rows = Int[];    sizehint!(center_rows, hint)
    f_tri_cos  = Float32[]; sizehint!(f_tri_cos, hint)
    f_entropy  = Float32[]; sizehint!(f_entropy, hint)
    f_tri_pcor = Float32[]; sizehint!(f_tri_pcor, hint)
    f_fit_hr   = Float32[]; sizehint!(f_fit_hr, hint)
    f_fit_r2   = Float32[]; sizehint!(f_fit_r2, hint)
    h_file = geom.template_h > 0f0 ? geom.template_h : geom.transmission_fwhm
    sh_str  = Float32[];    sh_effn = Float32[]; sh_best = Float32[]
    sh_disp = Float32[];    sh_n70  = UInt8[];   sh_rank = UInt16[]
    out_fit = ntuple(_ -> Float32[], 8)
    out_shd = ntuple(_ -> Float32[], 8)

    blk_centers = Int[]                       # center rows of one precursor's meta-scans

    i = 1
    @inbounds while i <= n
        j0 = i
        while i <= n && pid[i] == pid[j0]; i += 1; end
        blk_lo, blk_hi = j0, i - 1
        pm = prec_mz[pid[j0]]
        # ---- pick the center row of each meta-scan, one cycle at a time ----
        # Preferred (unchanged): a bin whose recorded isolation window contains the precursor
        # m/z. Added fallback: when NO bin of that cycle contains it, re-anchor to the bin
        # nearest the precursor m/z rather than discarding the whole meta-scan.
        #
        # The old behaviour dropped 41.1% of (precursor, cycle) groups and 26.2% of deconvolved
        # rows on Sciex ZT A_REP1, with essentially no target/decoy discrimination (49.9% vs
        # 48.5% decoy) -- the discard is a geometric accident, because NNLS routinely zeroes a
        # precursor at its own center bin (22% of expanded bins fall below the 1e-6 weight
        # floor). For the dropped groups the nearest surviving bin is a median of 1.46 Da off
        # center, ~86% transmission at the measured 6.3 Da FWHM, so the meta-scan is still
        # physically meaningful. Measured: +394 precursors, +46 protein groups on A_REP1.
        empty!(blk_centers)
        rr0 = blk_lo
        while rr0 <= blk_hi
            cy0 = cycles[Int(scn[rr0])]
            rr1 = rr0
            while rr1 < blk_hi && cycles[Int(scn[rr1 + 1])] == cy0
                rr1 += 1
            end
            n_contained = 0
            for r in rr0:rr1
                c = Int(scn[r])
                if abs(pm - cmzs[c]) <= hws[c] / 2
                    push!(blk_centers, r)
                    n_contained += 1
                end
            end
            if n_contained == 0
                rbest = 0
                dbest = Inf32
                for r in rr0:rr1
                    c = Int(scn[r])
                    d = abs(pm - cmzs[c])
                    if d < dbest            # NaN center m/z never wins: keeps MS1 rows out
                        dbest = d
                        rbest = r
                    end
                end
                rbest != 0 && push!(blk_centers, rbest)
            end
            rr0 = rr1 + 1
        end

        for r in blk_centers
            c = Int(scn[r])

            fill!(w, 0f0)
            if _have_frags
                for b in 1:8; fill!(Fbuf[b], 0f0); end
            end
            # Bins are ~contiguous rows around r (sorted by scan); scan a small window and bin
            # by scan-index offset. ±L rather than ±k, because a precursor need not have a row
            # in every bin of its meta-scan.
            if _have_fs
                fill!(accf, 0f0); fill!(accs, 0f0)
            end
            lo = max(blk_lo, r - L); hi = min(blk_hi, r + L)
            for rr in lo:hi
                d = Int(scn[rr]) - c
                if -k <= d <= k
                    w[d+k+1] += wt[rr]
                    if _have_frags
                        for b in 1:8; Fbuf[b][d+k+1] += fcols[b][rr]; end
                    end
                    if _have_fs
                        # Matched filter: weight each bin by its expected transmission.
                        tw = tri[d+k+1]
                        oi = perm[rr]
                        for b in 1:8
                            accf[b] += tw * fit_cols[b][oi]
                            accs[b] += tw * shd_cols[b][oi]
                        end
                    end
                end
            end
            if _have_fs
                for b in 1:8
                    push!(out_fit[b], accf[b]); push!(out_shd[b], accs[b])
                end
            end

            push!(center_rows, r)

            # ---- weight-profile features (inline; `w` is reused next iteration) ----
            cosv, entv = _zt_profile_features(w, k, tri, tnorm)
            push!(f_tri_cos, cosv)
            push!(f_entropy, entv)
            push!(f_tri_pcor, _frag_pcor(w, tri))
            hr, r2 = _zt_local_triangle_fit(w, k, geom.bin_step, pm - cmzs[c], h_file;
                                            lim = zt_fit_limit_da(h_file))
            push!(f_fit_hr, hr); push!(f_fit_r2, r2)

            # ---- within-metascan shape features (fragment profile vs weight profile) ----
            str = 0f0; effn = 0f0; best = 0f0; disp = 0f0; n70 = UInt8(0); rnk = UInt16(0)
            if _have_frags
                mask = UInt16(0); empty!(apx)
                for b in 1:8
                    fp = Fbuf[b]
                    mx = 0f0; for v in fp; v > mx && (mx = v); end
                    has_sig[b] = mx > 0f0
                    if has_sig[b]
                        cfw[b] = _frag_pcor(fp, w)
                        if cfw[b] > 0.7f0; n70 += UInt8(1); mask |= UInt16(1) << (b - 1); end
                        ai = 1; vm = fp[1]
                        for t in 2:L; fp[t] > vm && (vm = fp[t]; ai = t); end
                        push!(apx, Float32(ai - (k + 1)))
                    else
                        cfw[b] = 0f0
                    end
                end
                str, effn = _positive_corr_summary(cfw, rank_weights)
                # best-consensus fragment: highest mean correlation to the other signal fragments
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
            push!(sh_str, str);   push!(sh_effn, effn); push!(sh_best, best)
            push!(sh_disp, disp); push!(sh_n70, n70);   push!(sh_rank, rnk)
        end
    end

    # center_rows are positions in the PERMED order; map back to original psms rows.
    meta = psms[perm[center_rows], :]
    meta[!, :zt_tri_cosine] = f_tri_cos
    meta[!, :zt_entropy]    = f_entropy
    meta[!, :zt_tri_pcor]   = f_tri_pcor
    meta[!, :zt_fit_h_ratio] = f_fit_hr
    meta[!, :zt_fit_r2]      = f_fit_r2
    meta[!, :frag_corr_strength_shape]                 = sh_str
    meta[!, :frag_corr_effective_n_shape]              = sh_effn
    meta[!, :frag_corr_best_shape]                     = sh_best
    meta[!, :frag_apex_dispersion_shape]               = sh_disp
    meta[!, :n_correlated_fragments_shape]             = sh_n70
    meta[!, :n_correlated_fragments_bitvec_rank_shape] = sh_rank

    # Replace the anchor bin's fitted/shadow spectra with transmission-weighted sums over the
    # whole meta-scan. `smoothed_2d_shadow_hellinger` (scoring.jl:1152) then contrasts integrated
    # spectra rather than single-bin ones. Hellinger is invariant to a common scale on both
    # vectors, so only the conditioning changes. Motivation: DIA-NN-only misses carry 4.73 of 8
    # ranks at the anchor bin vs 6.90 across the whole meta-scan (controls 5.38 -> 7.20), and 82%
    # of a precursor's fragments are intermittent across the bins of one meta-scan.
    if _have_fs
        for b in 1:8
            meta[!, Symbol("fitted_frag$(b)_int")] = out_fit[b]
            meta[!, Symbol("shadow_frag$(b)_int")] = out_shd[b]
        end
    end
    return meta
end
