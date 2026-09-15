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
# Scanning-quad (ZT) quadrupole tuning.
#
# Fundamentally different from the standard fit and much simpler. The standard path probes the
# M0/M1 isotope ratio of a precursor seen in two abutting windows and fits Razo to it. A swept
# quad has no such pair: a precursor appears in ~2k+1 bins of ONE meta-scan, and because the
# deconvolution runs with a wide SQUARE box (installed in ensure_zt_geometry!, flat across the
# whole meta-scan) the fitted weight in each bin is already proportional to true transmission:
#
#     w(bin) ~ abundance * T(centerMz_bin - precursor_centre_of_mass)
#
# So this is an ordinary search, wide square window, read the weights.
#
# `h` is estimated WITHOUT normalising anything. Per meta-scan least squares of
# `w = a - b*|Δ|` gives `h = a/b`, scale-free because a and b share the unknown abundance.
# Every normalisation we tried pinned a bin and manufactured structure: dividing by the observed
# max forces the apex to 1.0, and dividing by the bin nearest the precursor pins
# |Δ| < bin_step/2 to exactly 1.0, creating a spurious "flat top" one Q1 step wide that then
# fits a trapezoid better than a Gaussian for purely circular reasons.
#
# Reference (Sciex ZT 5 Da/5 min, A_REP1, 68,646 meta-scans): h = 6.17 Da, h/bin_step = 6.04.
# The estimate is stable to ~1.3% with only 200 meta-scans.
# ============================================================================

"""
    zt_isotope_center_offset(prec_mz, charge) -> Float32

Offset (Da) from monoisotopic m/z to the intensity-weighted centre of the transmitted isotope
envelope. Averagine: neutral mass M carries ~M*0.000472 expected 13C, so with M0+M1 the centroid
sits `[λ/(1+λ)] * 1.00336 / z` above M0. This RISES with m/z because heavier peptides hold more
carbon — measured apex offset runs -0.02 Da at m/z 400 to +0.31 Da at m/z 900, which no
fixed-centre model reproduces.
"""
@inline function zt_isotope_center_offset(prec_mz::Real, charge::Integer)
    z = max(Int(charge), 1)
    M = Float32(prec_mz) * z - z * 1.007276f0
    λ = 0.000472f0 * M
    return Float32((λ / (1f0 + λ)) * 1.00336f0 / z)
end

"""
    zt_isotope_center_offset(iso_splines, prec_mz, charge, sulfur_count; n_iso = 4) -> Float32

Composition-aware version: the offset from monoisotopic m/z to the intensity-weighted centroid
of the first `n_iso` isotopes, using the library's Goldfarb isotope splines (indexed by neutral
mass and sulfur count) instead of the averagine expectation.

`n_iso = 2` by default, matching what the deconvolution measures: the search models two
precursor isotopes and one fragment isotope, so the fitted weight tracks M0+M1 transmission
only. Measured on A_REP1 (2,232 meta-scans, median R²): mono 0.916, averagine M0+M1 0.921,
spline n_iso=2 0.920, n_iso=3 0.901, n_iso=4 0.888 — the full-envelope centroid over-shifts.
"""
@inline function zt_isotope_center_offset(iso_splines::IsotopeSplineModel{Float32},
                                          prec_mz::Real, charge::Integer,
                                          sulfur_count::Integer; n_iso::Int = 2)
    z = max(Int(charge), 1)
    M = Float32(prec_mz) * z - z * 1.007276f0
    sidx = min(Int(sulfur_count), 5)
    num = 0f0; den = 0f0
    @inbounds for i in 0:(n_iso - 1)
        pi_ = max(iso_splines(sidx, i, M), 0f0)
        num += pi_ * i; den += pi_
    end
    den <= 0f0 && return zt_isotope_center_offset(prec_mz, charge)
    return Float32((num / den) * 1.00336f0 / z)
end

"""
    zt_center_offsets(precursors, iso_splines) -> Function

Per-precursor apex offset lookup `pid -> Float32`. Uses the splines when given, else averagine.
"""
function zt_center_offsets(precursors, iso_splines)
    prec_mz = Vector{Float32}(getMz(precursors))
    prec_z  = Vector{UInt8}(getCharge(precursors))
    if iso_splines === nothing
        return pid -> zt_isotope_center_offset(prec_mz[pid], prec_z[pid])
    end
    prec_s = Vector{UInt8}(getSulfurCount(precursors))
    return pid -> zt_isotope_center_offset(iso_splines, prec_mz[pid], prec_z[pid], prec_s[pid])
end

"""
    ZTQuadFitResult

Per-file triangle fit. `h` is apex-to-zero half-width in Da; `k_implied = round(h / bin_step)`
is the physically implied meta-scan expansion half-width.
"""
struct ZTQuadFitResult
    h::Float32
    h_iqr_lo::Float32
    h_iqr_hi::Float32
    median_r2::Float32
    n_metascans::Int
    k_implied::Int
end

"""
    fit_zt_triangle_from_psms(psms, spectra, precursors, geom; fit_limit_da, min_bins, min_pts)

Per-meta-scan least squares of `w = a - b*|Δ - μ_com|`, pooled by median of `h = a/b`.
Returns `nothing` when support is too thin.

`fit_limit_da` restricts each regression to the reliable core. The fit extrapolates to the zero
crossing, so it need not observe it: fitting only ±2 bins of a 6-bin profile still recovers h to
+8%, which is why a too-small `metascan_k` can still reveal that it should be larger.
"""
function fit_zt_triangle_from_psms(psms::DataFrame, spectra::MassSpecData, precursors,
                                   geom::ZTGeometry; iso_splines = nothing, kwargs...)
    # Two passes: the first fits the core of the geometry's own estimate of h (half the
    # expansion span); the second refits over +/-h/2 of that result, so the fitted fraction of
    # the profile is the same on every method (see ZT_FIT_CORE_FRACTION).
    fit1, hist = _fit_zt_triangle_core(psms, spectra, precursors, geom;
                                       iso_splines = iso_splines,
                                       fit_limit_da = zt_fit_limit_da(geom), kwargs...)
    fit1 === nothing && return nothing, hist
    fit2, hist2 = _fit_zt_triangle_core(psms, spectra, precursors, geom;
                                        iso_splines = iso_splines,
                                        fit_limit_da = zt_fit_limit_da(fit1.h), kwargs...)
    return fit2 === nothing ? (fit1, hist) : (fit2, hist2)
end

function _fit_zt_triangle_core(psms::DataFrame, spectra::MassSpecData, precursors,
                               geom::ZTGeometry;
                               iso_splines = nothing,
                               fit_limit_da::Float32,
                                   min_bins::Int = 5,
                                   min_pts::Int = 4,
                                   min_metascans::Int = 25)
    (nrow(psms) == 0 || geom.metascan_k <= 0) && return nothing, Int[]
    (hasproperty(psms, :precursor_idx) && hasproperty(psms, :scan_idx) &&
     hasproperty(psms, :weight)) || return nothing, Int[]

    prec_mz::Vector{Float32} = Vector{Float32}(getMz(precursors))
    offset_of = zt_center_offsets(precursors, iso_splines)
    cmzs::Vector{Float32}    = Float32.(coalesce.(getCenterMzs(spectra), NaN32))
    cycles::Vector{UInt32}   = UInt32.(getCycleIdxs(spectra))

    pid = Vector{UInt32}(psms[!, :precursor_idx])
    scn = Vector{UInt32}(psms[!, :scan_idx])
    wt  = Vector{Float32}(psms[!, :weight])

    ord = sortperm(collect(zip(pid, cycles[scn])))
    hs = Float32[]; r2s = Float32[]; hist = zeros(Int, 16)
    xs = Float32[]; ys = Float32[]
    i = 1; n = length(ord)
    @inbounds while i <= n
        j = i; p0 = pid[ord[i]]; c0 = cycles[scn[ord[i]]]
        while j <= n && pid[ord[j]] == p0 && cycles[scn[ord[j]]] == c0; j += 1; end
        hist[min(j - i, 16)] += 1
        if j - i >= min_bins
            μ = prec_mz[p0] + offset_of(p0)
            empty!(xs); empty!(ys)
            for t in i:(j - 1)
                d = cmzs[scn[ord[t]]] - μ
                if isfinite(d) && abs(d) <= fit_limit_da
                    push!(xs, abs(d)); push!(ys, wt[ord[t]])
                end
            end
            nx = length(xs)
            if nx >= min_pts
                sx = sum(xs); sy = sum(ys)
                sxx = zero(Float32); sxy = zero(Float32)
                for t in 1:nx; sxx += xs[t]*xs[t]; sxy += xs[t]*ys[t]; end
                den = nx*sxx - sx*sx
                if den > 0
                    b = -(nx*sxy - sx*sy)/den          # slope magnitude
                    a = (sy + b*sx)/nx
                    if a > 0 && b > 0
                        ymean = sy/nx; ssr = zero(Float32); sst = zero(Float32)
                        for t in 1:nx
                            pr = a - b*xs[t]
                            ssr += (ys[t]-pr)^2; sst += (ys[t]-ymean)^2
                        end
                        sst > 0 && (push!(hs, a/b); push!(r2s, 1f0 - ssr/sst))
                    end
                end
            end
        end
        i = j
    end
    length(hs) < min_metascans && return nothing, hist
    sort!(hs)
    q(p) = hs[clamp(round(Int, p*length(hs)), 1, length(hs))]
    h = q(0.5)
    return ZTQuadFitResult(h, q(0.25), q(0.75), median(r2s), length(hs),
                           max(1, round(Int, h / geom.bin_step))), hist
end

"""
    plot_zt_triangle(fit, psms, spectra, precursors, geom, fname) -> Vector{Plots.Plot}

QC for the ZT transmission triangle, returned as plot objects so they join the combined
`quad_transmission_plots.pdf` written by `summarize_results!`, exactly like the Razo plots.

1. pooled empirical profile with the fitted triangle. Points are per-Δm/z medians of `w / a`,
   where `a` is that meta-scan's OWN fitted intercept — scale-free, and it pins nothing.
   (Dividing by the observed max forces the apex to 1.0; dividing by the bin nearest the
   precursor pins |Δ| < bin_step/2 to exactly 1.0 and manufactures a flat top one Q1 step wide.)
2. distribution of per-meta-scan `h`, median and IQR marked.
3. `h` vs precursor m/z — flat means a uniform quadrupole.
"""
function plot_zt_triangle(fit::ZTQuadFitResult, psms::DataFrame, spectra::MassSpecData,
                          precursors, geom::ZTGeometry, fname::AbstractString;
                          iso_splines = nothing)
    out = Plots.Plot[]
    try
        prec_mz::Vector{Float32} = Vector{Float32}(getMz(precursors))
        offset_of = zt_center_offsets(precursors, iso_splines)
        cmzs::Vector{Float32}    = Float32.(coalesce.(getCenterMzs(spectra), NaN32))
        cycles::Vector{UInt32}   = UInt32.(getCycleIdxs(spectra))
        pid = Vector{UInt32}(psms[!, :precursor_idx]); scn = Vector{UInt32}(psms[!, :scan_idx])
        wt  = Vector{Float32}(psms[!, :weight])
        ord = sortperm(collect(zip(pid, cycles[scn])))
        dall = Float32[]; wall = Float32[]; hs = Float32[]; mzs = Float32[]
        xs = Float32[]; ys = Float32[]
        i = 1; n = length(ord)
        @inbounds while i <= n
            j = i; p0 = pid[ord[i]]; c0 = cycles[scn[ord[i]]]
            while j <= n && pid[ord[j]] == p0 && cycles[scn[ord[j]]] == c0; j += 1; end
            if j - i >= 5
                μ = prec_mz[p0] + offset_of(p0)
                empty!(xs); empty!(ys)
                for t in i:(j-1)
                    d = cmzs[scn[ord[t]]] - μ
                    isfinite(d) && abs(d) <= zt_fit_limit_da(fit.h) && (push!(xs, abs(d)); push!(ys, wt[ord[t]]))
                end
                nx = length(xs)
                if nx >= 4
                    sx=sum(xs); sy=sum(ys); sxx=zero(Float32); sxy=zero(Float32)
                    for t in 1:nx; sxx+=xs[t]*xs[t]; sxy+=xs[t]*ys[t]; end
                    den = nx*sxx - sx*sx
                    if den > 0
                        b = -(nx*sxy - sx*sy)/den; a = (sy + b*sx)/nx
                        if a > 0 && b > 0
                            push!(hs, a/b); push!(mzs, prec_mz[p0])
                            for t in i:(j-1)
                                d = cmzs[scn[ord[t]]] - μ
                                isfinite(d) && (push!(dall, d); push!(wall, wt[ord[t]]/a))
                            end
                        end
                    end
                end
            end
            i = j
        end
        isempty(dall) && return out
        lim = Float32(geom.metascan_k + 1) * geom.bin_step
        nb = 48; edges = collect(range(-lim, lim; length = nb+1))
        ctr = Float32[]; med = Float32[]
        for q in 1:nb
            sel = findall(t -> edges[q] <= dall[t] < edges[q+1], eachindex(dall))
            length(sel) >= 20 || continue
            push!(ctr, Float32((edges[q]+edges[q+1])/2)); push!(med, median(@view wall[sel]))
        end
        gx = collect(range(-lim, lim; length = 300))
        tri = [max(0f0, 1f0 - abs(x)/fit.h) for x in gx]
        p1 = Plots.scatter(ctr, med; ms=3, label="empirical median",
            xlabel="Δm/z from isotope centre of mass (Da)", ylabel="weight / intercept",
            title="$fname — ZT transmission  h=$(round(fit.h;digits=2)) Da  k_implied=$(fit.k_implied)")
        Plots.plot!(p1, gx, tri; lw=2, lc=:red, label="fitted triangle")
        Plots.vline!(p1, [-fit.h, fit.h]; lc=:black, ls=:dash, lw=1, label="±h")
        push!(out, p1)
        # Near-flat slopes give h in the hundreds of Da; clip so the bulk is readable.
        hclip = filter(x -> x <= 3f0 * fit.h, hs)
        p2 = Plots.histogram(hclip; bins=60, legend=false, xlabel="per-meta-scan h (Da)",
            ylabel="meta-scans", title="$fname — h distribution (n=$(fit.n_metascans), $(length(hs)-length(hclip)) >3h clipped, median R²=$(round(fit.median_r2;digits=3)))")
        Plots.vline!(p2, [fit.h]; lc=:red, lw=2)
        Plots.vline!(p2, [fit.h_iqr_lo, fit.h_iqr_hi]; lc=:red, ls=:dash, lw=1)
        push!(out, p2)
        o2 = sortperm(mzs); step = max(length(o2) ÷ 12, 30)
        bx = Float32[]; by = Float32[]
        for q in 1:step:(length(o2)-step+1)
            idx = o2[q:q+step-1]
            push!(bx, median(@view mzs[idx])); push!(by, median(@view hs[idx]))
        end
        if length(bx) >= 3
            p3 = Plots.plot(bx, by; marker=:circle, ms=3, legend=false, xlabel="precursor m/z",
                ylabel="h (Da)", title="$fname — h vs m/z (flat = uniform quadrupole)")
            Plots.hline!(p3, [fit.h]; lc=:red, ls=:dash)
            push!(out, p3)
        end
    catch err
        @user_warn "ZT transmission QC plot failed for $fname: $err"
    end
    return out
end

# ============================================================================
# PSM collection for the triangle fit.
#
# The Razo collection path (`progressive_quad_psm_collection!`) is the wrong tool here for two
# independent reasons, both fatal:
#   1. it reduces to ONE row per precursor (`combine(groupby(_, :precursor_idx), first)`) and
#      drops every charge but 2 — correct for one isotope-ratio observation per precursor, but
#      the triangle fit regresses on the ~2k+1 bin weights WITHIN each meta-scan;
#   2. it searches a TIC/m-z prioritised scan subset. `expand_to_metascans!` unions candidates
#      from neighbour scans that were themselves searched, and the fused scan loop only emits
#      rows for searched scans, so a scattered subset yields 1-2 bins per (precursor, cycle).
# So: search WHOLE cycles (contiguous MS2 ranges), keep the raw per-(precursor, scan) output of
# `library_search`, and decide confidence per (precursor, cycle) group — never per row, because
# outer-bin rows are legitimately weak and would be culled by a per-row FDR cut, leaving only
# the apex bin.
# ============================================================================

"""
    zt_count_fittable_metascans(psms, spectra, min_bins) -> Int

Number of (precursor, cycle) groups with at least `min_bins` rows — the support the triangle
fit actually consumes.
"""
function zt_count_fittable_metascans(psms::DataFrame, spectra::MassSpecData, min_bins::Int)
    nrow(psms) == 0 && return 0
    cycles = UInt32.(getCycleIdxs(spectra))
    pid = Vector{UInt32}(psms[!, :precursor_idx])
    scn = Vector{UInt32}(psms[!, :scan_idx])
    keys_ = collect(zip(pid, cycles[scn]))
    sort!(keys_)
    n = 0; i = 1; m = length(keys_)
    while i <= m
        j = i
        while j <= m && keys_[j] == keys_[i]; j += 1; end
        (j - i) >= min_bins && (n += 1)
        i = j
    end
    return n
end

"""
    collect_zt_quad_psms(spectra, search_context, params, ms_file_idx;
                         target_metascans, min_bins, fdr_threshold)
        -> (psms::DataFrame, n_fittable::Int, n_cycles_used::Int)

Collect raw per-(precursor, scan) PSMs across whole cycles for the ZT triangle fit.

Cycles are visited in a bit-reversed (van der Corput) order so any prefix is spread evenly over
the gradient; the count doubles until `target_metascans` fittable groups exist or the file is
exhausted. Confidence is decided per (precursor, cycle): the group's best-scoring row must be a
target at `fdr_threshold` (q-values computed over group representatives, one observation per
meta-scan). ALL rows of a passing group are kept.

Requires the ZT meta-scan view in `library_search` (`zt_qtune`): meta-scan expansion on and the
wide square deconvolution box, so weights track transmission. Both are keyed off
`params isa QuadTuningSearchParameters` there.
"""
function collect_zt_quad_psms(spectra::MassSpecData, search_context::SearchContext,
                              params::QuadTuningSearchParameters, ms_file_idx::Int64;
                              target_metascans::Int = ZT_QUAD_TARGET_METASCANS,
                              min_bins::Int = 5,
                              fdr_threshold::Float64 = 0.01,
                              min_score::Integer = first(TUNING_SCORE_TIERS),
                              n_required_top::Int = TUNING_N_REQUIRED_TOP)
    ranges = zt_cycle_scan_ranges(spectra)
    n_cyc = length(ranges)
    n_cyc == 0 && return DataFrame(), 0, 0
    precursors = getPrecursors(getSpecLib(search_context))
    fdr_scale  = getLibraryFdrScaleFactor(search_context)
    cycles     = UInt32.(getCycleIdxs(spectra))

    # Bit-reversed visiting order: prefixes are stratified across the run.
    nb = max(1, ceil(Int, log2(n_cyc)))
    visit = Int[]
    for q in 0:(2^nb - 1)
        r = 0
        for b in 0:(nb - 1); (q >> b) & 1 == 1 && (r |= 1 << (nb - 1 - b)); end
        r < n_cyc && push!(visit, r + 1)
    end

    setBitVecFilter!(search_context, ms_file_idx,
                     make_top_n_required_lut(n_required_top, Int(min_score)))
    raw = DataFrame(); kept = DataFrame()
    n_fit = 0; used = 0
    batch = max(8, cld(n_cyc, 64))
    try
        while used < n_cyc
            hi = min(used + batch, n_cyc)
            scan_idxs = Int[]
            for c in visit[(used + 1):hi]; append!(scan_idxs, ranges[c]); end
            chunk = library_search(spectra, search_context, params, ms_file_idx;
                                   scan_indices = scan_idxs)
            if !isempty(chunk)
                add_tuning_search_columns!(chunk, spectra,
                    getIsDecoy(precursors), getIrt(precursors),
                    getCharge(precursors), getRetentionTimes(spectra), getTICs(spectra))
                raw = isempty(raw) ? chunk : (append!(raw, chunk); raw)
            end
            used = hi
            batch *= 2

            if nrow(raw) >= 50
                scored = copy(raw)
                score_presearch!(scored)
                scored[!, :cycle] = cycles[Vector{UInt32}(scored[!, :scan_idx])]
                # One observation per meta-scan: its best row.
                sort!(scored, [:precursor_idx, :cycle, order(:prob, rev = true)])
                reps = combine(groupby(scored, [:precursor_idx, :cycle]),
                               :prob => first => :prob, :target => first => :target)
                reps[!, :q_value] = zeros(Float16, nrow(reps))
                get_qvalues!(reps[!, :prob], reps[!, :target], reps[!, :q_value];
                             fdr_scale_factor = fdr_scale)
                pass = Set(zip(reps.precursor_idx[(reps.q_value .<= Float16(fdr_threshold)) .& reps.target],
                               reps.cycle[(reps.q_value .<= Float16(fdr_threshold)) .& reps.target]))
                kept = scored[[(p, c) in pass for (p, c) in zip(scored.precursor_idx, scored.cycle)], :]
                n_fit = zt_count_fittable_metascans(kept, spectra, min_bins)
            end
            @debug_l1 "  ZT quad collect: $used/$n_cyc cycles, $(nrow(raw)) raw rows, " *
                      "$(nrow(kept)) confident rows, $n_fit fittable meta-scans"
            n_fit >= target_metascans && break
        end
    finally
        delete!(search_context.bitvec_filter, ms_file_idx)
    end
    return kept, n_fit, used
end
