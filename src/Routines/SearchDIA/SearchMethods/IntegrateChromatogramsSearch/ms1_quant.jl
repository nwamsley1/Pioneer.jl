# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

#==========================================================
MS1-level quantitation for DDA.

DDA gives no MS2 chromatogram: under dynamic exclusion a precursor is fragmented in ~1 MS2 scan
(LFQBench PXD028735: n_scans median 1), so the MS2 trace has no extent and integrates to zero. MS1 is
the only place a DDA precursor has a chromatogram. For each MS1 scan we fit every candidate precursor's
isotope envelope (M0..M2) to the observed peaks with a non-negative Huber/L2 solve, giving one MS1
chromatogram point per candidate per scan; those traces are integrated in place of the MS2 traces.

Ported from origin/feat/dda-ms1-quant (2026-08); design notes in its ms1_quant_dda_plan.tex.
==========================================================#

"""
Half-width (minutes) of the DDA MS1 extraction window around each precursor's identifying MS2 scan.

It replaces the DIA RT tolerance rather than flooring it: `RTBinnedTolerance` is derived from MS2
chromatographic FWHM, which is 0 in DDA, and its sparse-bin fallback can be minutes wide. The window must
also be generous because the identifying MS2 scan is triggered on a threshold crossing, not at the apex
(LFQBench: 53% of apexes > 5 s away). Measured there: MS1 FWHM 15 s median, integrated span ~2.3x FWHM,
only 3.5% of FWHMs above 30 s.
"""
const MS1_DDA_RT_TOL_DEFAULT = 0.5f0

"""
    ms1_smoothing_rt_scale(calibrated::Float32) -> Float32

RT divisor (minutes) for Whittaker-Henderson smoothing of DDA MS1 traces, so the smoothing is set
in time rather than relative to the extraction window (the order-3 penalty grows as window^6).
Normally the per-file value from `calibrate_ms1_window` (`ms1_smoothing_scale`); 0 or negative
means no calibration and falls back to MS1_SMOOTH_SCALE_DEFAULT.
"""
function ms1_smoothing_rt_scale(calibrated::Float32)::Float32
    return calibrated > 0f0 ? calibrated : MS1_SMOOTH_SCALE_DEFAULT
end

"""
    collect_ms1_window_precursors!(precs_temp, rt_index, rt_start_idx, rt_stop_idx,
                                   precursors_passing, precursor_rt_map, scan_rt, rt_tol) -> Int

Candidate precursors for one MS1 scan, written to `precs_temp[1:n]`. Unlike the MS2 version there is no
isolation window, so the RT window around each precursor's identifying scan is the only bound, plus, on ion-mobility
data (`precursor_im_map`), a mobility band of `im_window_scans` around the IM scan of its best PSM.
"""
function collect_ms1_window_precursors!(
    precs_temp::Vector{UInt32},
    rt_index::retentionTimeIndex,
    rt_start_idx::Int64, rt_stop_idx::Int64,
    precursors_passing::Set{UInt32},
    precursor_rt_map::Dict{UInt32, NTuple{2, Float32}},
    scan_rt::Float32,
    rt_tol::Float32,
    rt_binned::Union{Nothing, RTBinnedTolerance} = nothing,   # per-precursor window by its MS2 RT
    im_scan::Float32 = 0f0,
    precursor_im_map::Union{Nothing, Dict{UInt32, Float32}} = nothing,
    im_window_scans::Float32 = 0f0,
)
    has_im = precursor_im_map !== nothing && im_window_scans > 0f0
    n = 0
    for rt_bin_idx in rt_start_idx:rt_stop_idx
        precs = rt_index.rt_bins[rt_bin_idx].prec
        for i in eachindex(precs)
            prec_idx = first(precs[i])
            prec_idx in precursors_passing || continue
            entry = get(precursor_rt_map, prec_idx, nothing)
            entry === nothing && continue
            tol = rt_binned === nothing ? rt_tol : get_rt_tol(rt_binned, first(entry))
            d = scan_rt - first(entry)                  # > 0: scan after the precursor's RT
            ((d > tol * last(entry)) | (-d > tol)) && continue   # last(entry): right-tail multiplier (CHROM_RHS_MULT)
            if has_im
                pim = get(precursor_im_map, prec_idx, NaN32)
                !isnan(pim) && abs(im_scan - pim) > im_window_scans && continue
            end
            n += 1
            n > length(precs_temp) && append!(precs_temp, Vector{UInt32}(undef, length(precs_temp)))
            precs_temp[n] = prec_idx
        end
    end
    return n
end

#==========================================================
(m/z, charge) grouping
==========================================================#

"""
m/z rounding for the grouping key (1e-5 Th, ~0.02 ppm at m/z 500): far finer than the MS1 tolerance, so
only exact duplicates collapse — every target/decoy pair (decoys are sequence permutations with the same
composition) and true isobars. Near-collinear distinct precursors keep their own columns and are handled by
the ridge term.
"""
const MS1_MZ_GROUP_SCALE = 100_000.0f0

@inline function ms1_group_key(prec_mz::Float32, prec_charge::UInt8)::UInt64
    # charge is part of the key: equal m/z at different charge has a different isotope spacing
    return (UInt64(round(Int64, prec_mz * MS1_MZ_GROUP_SCALE)) << 8) | UInt64(prec_charge)
end

"""
Per-scan scratch mapping candidates onto design-matrix columns, one column per distinct (m/z, charge).
`col_rep[c]` is the precursor whose envelope column `c` models; `col_of_cand[j]` is candidate `j`'s column.
"""
mutable struct Ms1MzGroups
    key_to_col::Dict{UInt64, UInt32}
    col_of_cand::Vector{UInt32}
    col_rep::Vector{UInt32}
    n_cols::Int
end
Ms1MzGroups(capacity::Int = 4096) =
    Ms1MzGroups(Dict{UInt64, UInt32}(), Vector{UInt32}(undef, capacity), Vector{UInt32}(undef, capacity), 0)

"""
    keep_ms1_covered!(precs_temp, n_cands, prec_mzs, prec_charges, coverage, scan_idx) -> Int

Segmented or BoxCar MS1 (MS1 scans measuring different parts of the precursor m/z range): keeps the candidates whose
M0 or M1 lies in a range MS1 scan `scan_idx` measures, compacting `precs_temp[1:n_cands]` in place, and returns
their number.
"""
function keep_ms1_covered!(precs_temp::Vector{UInt32}, n_cands::Int, prec_mzs::AbstractVector{<:AbstractFloat},
                           prec_charges::AbstractVector{<:Integer}, coverage::Ms1Coverage, scan_idx::Integer)::Int
    k = 0
    @inbounds for j in 1:n_cands
        p = precs_temp[j]
        m0 = Float32(prec_mzs[p])
        m1 = m0 + C13_C12_MASS_DIFF_F32 / Float32(prec_charges[p])
        (covers(coverage, scan_idx, m0) || covers(coverage, scan_idx, m1)) || continue
        k += 1
        precs_temp[k] = p
    end
    return k
end

function build_ms1_mz_groups!(groups::Ms1MzGroups, precs_temp::Vector{UInt32}, n_cands::Int,
                              prec_mzs::AbstractArray{Float32}, prec_charges::AbstractArray{UInt8})
    empty!(groups.key_to_col)
    groups.n_cols = 0
    n_cands > length(groups.col_of_cand) && resize!(groups.col_of_cand, n_cands)
    n_cands > length(groups.col_rep) && resize!(groups.col_rep, n_cands)
    @inbounds for j in 1:n_cands
        prec_idx = precs_temp[j]
        key = ms1_group_key(prec_mzs[prec_idx], prec_charges[prec_idx])
        col = get(groups.key_to_col, key, UInt32(0))
        if iszero(col)
            groups.n_cols += 1
            col = UInt32(groups.n_cols)
            groups.key_to_col[key] = col
            groups.col_rep[groups.n_cols] = prec_idx
        end
        groups.col_of_cand[j] = col
    end
    return groups.n_cols
end

#==========================================================
Design matrix and solver
==========================================================#

"""Precursor isotopes modelled per column (M0, M1, M2)."""
const MS1_N_ISOTOPES = 3

"""
Ridge penalty for the MS1 solve, applied to raw (isotope-abundance) columns.

The solve minimises ½‖Hw − y‖² + λ‖w‖² (w ≥ 0), so weights scale linearly with intensity at any λ and no
per-instrument intensity scaling is needed. An isolated column is shrunk by `Σh² / (Σh² + 2λ)`, where `Σh²` is the
squared abundances of the isotopes matched in that spectrum (~0.2-0.4). At λ = 1 that factor (~0.09-0.16) depends
on how many isotopes matched, which varies between spectra and runs; at λ = 0.01 it is ~0.95-0.97 and the weights
are close to unbiased, at the cost of more redistribution between columns that share peaks. λ = 0 does not
converge.
"""
const MS1_LAMBDA_DEFAULT = 0.01f0

"""
Huber delta for the MS1 solve, high enough to make the loss quadratic (OLS + L2). The MS2 chromatogram value
(300) is far below MS1 intensities (1e3-1e9), which puts every residual in the linear regime where the
Newton step degenerates and all weights collapse to a plateau.
"""
const MS1_HUBER_DELTA_DEFAULT = 1.0f10

"""The MS1 solver: Huber (quadratic at MS1 intensities) with the ridge penalty above."""
ms1_solver() = HuberSolver(MS1_HUBER_DELTA_DEFAULT, MS1_LAMBDA_DEFAULT, Int64(50), Int64(100), 10.0f0, 10.0f0, L2Norm())

"""Per-file MS1 solve statistics (scans solved, columns, outer iterations, solves that hit the iteration cap)."""
const MS1_SOLVE_STATS = (scans = Threads.Atomic{Int}(0), cols = Threads.Atomic{Int}(0),
                         iters = Threads.Atomic{Int}(0), capped = Threads.Atomic{Int}(0))
reset_ms1_solve_stats!() = foreach(a -> (a[] = 0), MS1_SOLVE_STATS)

"""
    _ms1_match_peak(mass_err_model, target, scan_mz, scan_int, n_peaks, scan_rt) -> (peak, intensity)

The observed peak matched to an isotope's theoretical m/z `target`, or (0, 0) when none is in tolerance.
`SimpleMassErrorModel`: the target is shifted by the fitted ppm offset (getMzBounds itself is offset-free), as the DIA
MS1 features do, and the peak nearest that corrected target inside its fixed window wins. `IntensityMassErrorModel`:
each candidate peak is bias-corrected and given its own window from its m/z, intensity and the scan RT
(`getCorrectedMzAndBounds`, as fragment matching does); among peaks whose window covers the target, the one whose
corrected m/z is nearest wins.
"""
@inline function _ms1_match_peak(m::SimpleMassErrorModel, target::Float32, scan_mz::AbstractVector,
                                 scan_int::AbstractVector, n_peaks::Int, ::Float32)
    centre = target * (1f0 + Float32(getMassOffset(m)) * 1f-6)
    lo, hi = getMzBounds(m, centre)
    p = searchsortedfirst(view(scan_mz, 1:n_peaks), lo)
    (p > n_peaks || Float32(scan_mz[p]) > hi) && return 0, 0f0
    best_p = p; best_i = Float32(scan_int[p]); best_d = abs(Float32(scan_mz[p]) - centre)
    q = p + 1
    @inbounds while q <= n_peaks && Float32(scan_mz[q]) <= hi
        dq = abs(Float32(scan_mz[q]) - centre)
        if dq < best_d
            best_i = Float32(scan_int[q]); best_p = q; best_d = dq
        end
        q += 1
    end
    return best_p, best_i
end

@inline function _ms1_match_peak(m::IntensityMassErrorModel, target::Float32, scan_mz::AbstractVector,
                                 scan_int::AbstractVector, n_peaks::Int, scan_rt::Float32)
    w = Float32(m.conservative_tol_da)                     # the widest window any peak can have
    centre = target + _mz_bias_da(m, target)               # where the peak sits: the span is offset-free, the bias is not
    p = searchsortedfirst(view(scan_mz, 1:n_peaks), centre - w)
    best_p = 0; best_i = 0f0; best_d = Inf32
    @inbounds while p <= n_peaks && Float32(scan_mz[p]) <= centre + w
        ip = Float32(scan_int[p])
        if ip > 0f0
            c, lo, hi = getCorrectedMzAndBounds(m, Float32(scan_mz[p]), ip, scan_rt)
            if lo <= target <= hi
                d = abs(c - target)
                if d < best_d
                    best_i = ip; best_p = p; best_d = d
                end
            end
        end
        p += 1
    end
    return best_p, best_i
end

"""
    build_ms1_design_matrix!(Hs, groups, n_cols, scan_mz, scan_int, n_peaks, prec_mzs, prec_charges,
                             prec_sulfurs, iso_splines, mass_err_model, m0_obs) -> nnz

One CSC column per (m/z, charge) group with up to `n_iso_model` rows: isotope `j` at
`mz + j * C13_C12_MASS_DIFF / z`, value from the sulfur-aware isotope splines, matched by `_ms1_match_peak` (the peak
nearest the bias-corrected target within the MS1 mass tolerance). Rows are
observed peak indices. As in the MS2 design matrix (`finalize_column!`), an isotope with no matched peak whose target
lies within the spectrum's m/z range gets its own row `n_peaks + k` with observed intensity 0, so predicting signal
where none was observed is penalised; a column with no matched isotope stays empty. `Hs.m` counts those rows.
No quadrupole transmission term: MS1 observes the whole m/z range. `m0_obs[c]` records the matched M0 intensity,
`n_iso[c]` the number of matched isotopes. `scan_rt` (minutes) feeds the RT term of an `IntensityMassErrorModel`.
"""
function build_ms1_design_matrix!(
    Hs::SparseArrayFused{UInt32, Float32},
    groups::Ms1MzGroups,
    n_cols::Int,
    scan_mz::AbstractVector,
    scan_int::AbstractVector,
    n_peaks::Int,
    prec_mzs::AbstractArray{Float32},
    prec_charges::AbstractArray{UInt8},
    prec_sulfurs::AbstractArray{UInt8},
    iso_splines::IsotopeSplineModel,
    mass_err_model,
    m0_obs::Vector{Float32},
    col_norm::Vector{Float32},
    n_iso::Vector{UInt8} = UInt8[],
    scan_rt::Float32 = 0f0,
    n_iso_model::Int = MS1_N_ISOTOPES,
)
    needed = n_cols * n_iso_model
    if needed + 1 > length(Hs.rowval)
        grow = needed + 1 - length(Hs.rowval)
        append!(Hs.rowval, zeros(UInt32, grow)); append!(Hs.nzval, zeros(Float32, grow))
        append!(Hs.x, zeros(Float32, grow)); append!(Hs.meta, fill(UInt8(1), grow))
        append!(Hs.rank, zeros(UInt8, grow))
    end
    n_cols + 2 > length(Hs.colptr) && append!(Hs.colptr, zeros(UInt32, n_cols + 2 - length(Hs.colptr)))
    n_cols > length(m0_obs) && resize!(m0_obs, n_cols)
    n_cols > length(col_norm) && resize!(col_norm, n_cols)
    n_cols > length(n_iso) && resize!(n_iso, n_cols)

    nnz = 0
    n_miss = 0
    mz_lo = n_peaks > 0 ? Float32(scan_mz[1]) : 0f0          # spectrum's observed m/z range: an isotope outside it
    mz_hi = n_peaks > 0 ? Float32(scan_mz[n_peaks]) : 0f0    # cannot be observed, so it is not a miss
    @inbounds for col in 1:n_cols
        Hs.colptr[col] = UInt32(nnz + 1)
        m0_obs[col] = 0f0
        prec_idx = groups.col_rep[col]
        prec_mz = prec_mzs[prec_idx]
        charge = prec_charges[prec_idx]
        mono_mass = (prec_mz - Float32(PROTON)) * Float32(charge)
        inv_charge = 1f0 / Float32(charge)
        sulfur = Int64(min(prec_sulfurs[prec_idx], UInt8(5)))
        n_match = 0
        matched = UInt32(0)                                  # bit j: isotope j has a row
        for j in 0:(n_iso_model - 1)
            target = prec_mz + Float32(j) * C13_C12_MASS_DIFF_F32 * inv_charge
            best_p, best_i = _ms1_match_peak(mass_err_model, target, scan_mz, scan_int, n_peaks, scan_rt)
            best_p == 0 && continue
            alpha = iso_splines(sulfur, j, mono_mass)
            alpha <= 0f0 && continue
            n_match += 1; matched |= UInt32(1) << j
            j == 0 && (m0_obs[col] = best_i)
            nnz += 1
            Hs.rowval[nnz] = UInt32(best_p)
            Hs.nzval[nnz] = Float32(alpha)
            Hs.x[nnz] = best_i
            Hs.meta[nnz] = pack_meta(true, UInt8(j))
            Hs.rank[nnz] = UInt8(0)
        end
        n_iso[col] = UInt8(n_match)
        # misses after the matches (rows above every peak index keep the column in row order); skipped for a column
        # with no match, which stays empty and gets weight 0
        if n_match > 0
            for j in 0:(n_iso_model - 1)
                target = prec_mz + Float32(j) * C13_C12_MASS_DIFF_F32 * inv_charge
                ((matched >> j) & UInt32(1)) == 0 || continue
                (mz_lo <= target <= mz_hi) || continue
                alpha = iso_splines(sulfur, j, mono_mass)
                alpha <= 0f0 && continue
                n_miss += 1; nnz += 1
                Hs.rowval[nnz] = UInt32(n_peaks + n_miss)
                Hs.nzval[nnz] = Float32(alpha)
                Hs.x[nnz] = 0f0
                Hs.meta[nnz] = pack_meta(false, UInt8(j))
                Hs.rank[nnz] = UInt8(0)
            end
        end
        # column norm (reported as col_norm_raw)
        first_nz = Int(Hs.colptr[col])
        ss = 0f0
        for k in first_nz:nnz; ss += Hs.nzval[k]^2; end
        col_norm[col] = sqrt(ss)
    end
    Hs.colptr[n_cols + 1] = UInt32(nnz + 1)
    Hs.n_vals = nnz
    Hs.n = n_cols
    Hs.m = n_peaks + n_miss
    return nnz
end

#==========================================================
MS1 scan loop
==========================================================#

"""
    build_chromatograms(..., ::MS1CHROM; kwargs...) -> DataFrame

MS1 counterpart of the MS2 `build_chromatograms` (same call signature, so `extract_chromatograms` drives
both). For every MS1 scan: candidates are the passing precursors whose identifying scan lies within
`ms1_rt_tol` (per RT bin with `ms1_rt_binned`); exact (m/z, charge) duplicates share one column; the isotope-envelope design matrix is
solved with Huber + L2; each candidate gets its column's weight (duplicated, not split, across a group).
Returns `rt, intensity, m0_intensity, scan_idx, precursor_idx, precursor_fraction_transmitted` (the last is
1: there is no isolation window at MS1), the columns `integrate_precursors` reads in combined-trace mode.
"""
function build_chromatograms(
    spectra::MassSpecData,
    scan_range::Vector{Int64},
    precursors_passing::Set{UInt32},
    precursor_rt_map::Dict{UInt32, NTuple{2, Float32}},
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64,
    ::MS1CHROM;
    scan_tic = nothing,
    precursor_im_map = nothing,
    im_window_scans::Float32 = 0f0,
    ms1_rt_tol::Float32 = 0f0,
    ms1_rt_binned::Union{Nothing, RTBinnedTolerance} = nothing,
    ms1_mass_err_model::Union{Nothing, AbstractMassErrorModel} = nothing,   # per-file MS1-quant model; nothing = tuning's
    ms1_coverage::Union{Nothing, Ms1Coverage} = nothing,   # m/z ranges each MS1 scan measures; nothing = all
)
    precs_temp = getPrecIds(search_data)
    decode_buf = getDecodeBuffer(search_data)
    groups = Ms1MzGroups()
    Hs = SparseArrayFused(UInt32(5000))
    weights = zeros(Float32, 8192)
    residuals = zeros(Float32, 8192)
    colnorm2 = zeros(Float32, 8192)
    m0_obs = zeros(Float32, 8192)
    col_norm = zeros(Float32, 8192)
    n_iso = zeros(UInt8, 8192)
    solver = ms1_solver()
    iso_splines = getIsoSplines(search_data)
    ms1_mass_err = ms1_mass_err_model === nothing ? getMs1MassErrorModel(search_context, ms_file_idx) : ms1_mass_err_model
    rt_tol = ms1_rt_binned !== nothing ? maximum(ms1_rt_binned.rt_tols) :   # bounds the RT-bin search
             ms1_rt_tol > 0f0 ? ms1_rt_tol : MS1_DDA_RT_TOL_DEFAULT
    im_scans = getImScans(spectra)
    # right-tail extended precursors (rhs_extended) are candidates up to CHROM_RHS_MULT * tol after their RT, so the
    # RT-bin search must reach that far back; the per-precursor test applies the exact window
    lookback = rt_tol * (any(v -> last(v) > 1f0, values(precursor_rt_map)) ? CHROM_RHS_MULT : 1f0)

    precursors = getPrecursors(getSpecLib(search_context))
    prec_mzs = getMz(precursors)
    prec_charges = getCharge(precursors)
    prec_sulfurs = getSulfurCount(precursors)

    out_rt = Float32[]; out_int = Float32[]; out_m0 = Float32[]; out_scan = UInt32[]; out_prec = UInt32[]
    out_norm = Float32[]; out_niso = UInt8[]

    for scan_idx in scan_range
        (scan_idx < 1 || scan_idx > length(spectra)) && continue
        getMsOrder(spectra, scan_idx) == 1 || continue
        rt = Float32(getRetentionTime(spectra, scan_idx))
        bin_start = max(searchsortedfirst(rt_index.rt_bins, rt - lookback, lt = (r, x) -> r.lb < x) - 1, 1)
        bin_stop = min(searchsortedlast(rt_index.rt_bins, rt + rt_tol, lt = (x, r) -> r.ub > x) + 1,
                       length(rt_index.rt_bins))
        n_cands = collect_ms1_window_precursors!(precs_temp, rt_index, bin_start, bin_stop,
                                                 precursors_passing, precursor_rt_map, rt, rt_tol, ms1_rt_binned,
                                                 im_scans === nothing ? 0f0 : Float32(im_scans[scan_idx]),
                                                 precursor_im_map, im_window_scans)
        if ms1_coverage !== nothing && !ms1_coverage.uniform
            n_cands = keep_ms1_covered!(precs_temp, n_cands, prec_mzs, prec_charges, ms1_coverage, scan_idx)
        end
        n_cands == 0 && continue
        n_cols = build_ms1_mz_groups!(groups, precs_temp, n_cands, prec_mzs, prec_charges)

        scan_mz, scan_int = getPeaks!(decode_buf, spectra, scan_idx)
        n_peaks = length(scan_mz)
        nnz = build_ms1_design_matrix!(Hs, groups, n_cols, scan_mz, scan_int, n_peaks,
                                       prec_mzs, prec_charges, prec_sulfurs, iso_splines, ms1_mass_err, m0_obs,
                                       col_norm, n_iso, rt, MS1_N_ISOTOPES)
        nnz == 0 && continue
        if n_cols > length(weights)
            resize!(weights, n_cols + 1000); resize!(colnorm2, n_cols + 1000)
        end
        Hs.m > length(residuals) && resize!(residuals, Hs.m + 1000)
        @inbounds for c in 1:n_cols; weights[c] = 0f0; end
        sres = solve_deconvolution!(solver, Hs, residuals, weights, colnorm2, getMu(search_data),
                                    getObserved(search_data), params.max_iter_outer, params.max_diff)
        Threads.atomic_add!(MS1_SOLVE_STATS.scans, 1); Threads.atomic_add!(MS1_SOLVE_STATS.cols, n_cols)
        if sres isa Tuple
            Threads.atomic_add!(MS1_SOLVE_STATS.iters, Int(sres[2]))
            sres[1] || Threads.atomic_add!(MS1_SOLVE_STATS.capped, 1)
        end
        scan32 = UInt32(scan_idx)
        @inbounds for j in 1:n_cands
            col = groups.col_of_cand[j]
            push!(out_rt, rt)
            push!(out_int, max(weights[col], 0f0))
            push!(out_norm, col_norm[col]); push!(out_niso, n_iso[col])
            push!(out_m0, m0_obs[col])
            push!(out_scan, scan32)
            push!(out_prec, precs_temp[j])
        end
    end
    return DataFrame(rt = out_rt, intensity = out_int, m0_intensity = out_m0, scan_idx = out_scan,
                     precursor_idx = out_prec, precursor_fraction_transmitted = ones(Float32, length(out_rt)),
                     col_norm_raw = out_norm, n_iso_matched = out_niso)
end

"""
    collapse_ms1_im_slices(chroms, spectra, prec_im, band) -> DataFrame

Ion-mobility data: the RT marginal of each precursor's MS1 grid over its mobility band. An MS1 frame arrives as
several mobility slices; the slices within `band` IM scans of the precursor's IM scan (`prec_im`, from its best PSM)
are summed per MS1 frame (one cycle) into one point at the frame's first in-band slice (`scan_idx`, `rt`). Used by
the per-file calibration (MS1 spacing, FWHM, apex offset), the counterpart of the RT marginal the 2D integration
takes over the same band. Without mobility data the input is returned.
"""
function collapse_ms1_im_slices(chroms::DataFrame, spectra::MassSpecData, prec_im::Dict{UInt32, Int}, band::Int)
    ims = getImScans(spectra)
    (nrow(chroms) == 0 || ims === nothing) && return chroms
    cyc = getCycleIdxs(spectra)
    keep = [abs(Int(ims[chroms.scan_idx[i]]) - get(prec_im, chroms.precursor_idx[i], typemin(Int) ÷ 2)) <= band
            for i in 1:nrow(chroms)]
    chroms = chroms[keep, :]
    chroms[!, :_cycle] = UInt32[UInt32(cyc[s]) for s in chroms.scan_idx]
    out = combine(groupby(chroms, [:precursor_idx, :_cycle]),
                  :rt => minimum => :rt, :intensity => sum => :intensity, :m0_intensity => sum => :m0_intensity,
                  :scan_idx => minimum => :scan_idx, :precursor_fraction_transmitted => first => :precursor_fraction_transmitted,
                  :col_norm_raw => maximum => :col_norm_raw, :n_iso_matched => maximum => :n_iso_matched)
    return select!(out, [:rt, :intensity, :m0_intensity, :scan_idx, :precursor_idx, :precursor_fraction_transmitted,
                         :col_norm_raw, :n_iso_matched])
end

#==========================================================
Per-run extraction window
==========================================================#

"""Calibration: number of most confident target precursors."""
const MS1_WINDOW_CALIB_N = 3000
"""Calibration extraction half-width h0 (minutes): the fitted iRT tolerance converted to RT, clamped to
[MS1_WINDOW_CALIB_MIN_CYCLES MS1 spacings, MS1_WINDOW_CALIB_MAX_HALF_WIDTH]."""
const MS1_WINDOW_CALIB_MIN_CYCLES = 20
const MS1_WINDOW_CALIB_MAX_HALF_WIDTH = 5.0f0
"""Calibration retry: h0 is doubled (up to MS1_WINDOW_CALIB_MAX_HALF_WIDTH, at most MS1_WINDOW_CALIB_MAX_RETRIES
times) when the window the peaks call for exceeds h0/2 or more than MS1_WINDOW_CALIB_MAX_CUT_FRAC of the calibration
peaks run into the trace ends."""
const MS1_WINDOW_CALIB_MAX_RETRIES = 3
const MS1_WINDOW_CALIB_MAX_CUT_FRAC = 0.01
"""Half-window = q99(|apex - MS2 scan|) + MS1_WINDOW_FWHM_MULT * q90(FWHM), clamped to [MS1_WINDOW_MIN_CYCLES MS1
cycles, h0]."""
const MS1_WINDOW_FWHM_MULT = 2.0f0
const MS1_WINDOW_MIN_CYCLES = 5
"""Fewest measurable calibration peaks for the window / smoothing statistics. At 500 peaks the window
(q99/q90 statistics) varies by about +/-6% and S (medians) by +/-2% (Astral DDA subsampling); below it the file
falls back to the iRT-based first guess."""
const MS1_WINDOW_CALIB_MIN_PEAKS = 500
"""Fallback smoothing scale as a fraction of the iRT-based first-guess half-width h0 (Astral 15 min: h0 ~37 s,
calibrated S ~12.4 s)."""
const MS1_SMOOTH_FALLBACK_H0_FRAC = 0.33f0
"""Smoothing scale S = max(MS1_SMOOTH_FWHM_MULT * median FWHM, MS1_SMOOTH_MIN_CYCLES * median MS1 spacing), the
same form and constants as DIA's chromatogram window (2 x max(2 FWHM, 10 cycles)). End to end, the equal-species
ratio error is lowest at S ~3.5-4 x FWHM on Astral DDA 15 min, Astral DDA 5 min and QE HF-X 150 min; more smoothing
improves precision but compresses ratios (7 x FWHM, tuned on Astral 15 min alone, compressed E. coli on the other
two). The factor is tied to the FWHM estimator in `_ms1_local_peak` (3-point smoothed, no baseline subtraction);
re-check it end to end if that estimator changes. Fallback (no calibration) MS1_SMOOTH_SCALE_DEFAULT = 20 s, in
minutes."""
const MS1_SMOOTH_FWHM_MULT = 4.0f0
const MS1_SMOOTH_MIN_CYCLES = 20
const MS1_SMOOTH_SCALE_DEFAULT = 20f0 / 60f0

"""Extraction half-width the calibration peaks call for: q99(|apex - MS2|) + MS1_WINDOW_FWHM_MULT * q90(FWHM)
(minutes); Inf when there are too few peaks to measure."""
function ms1_window_need(offsets::AbstractVector{<:Real}, fwhms::AbstractVector{<:Real})::Float32
    length(offsets) < MS1_WINDOW_CALIB_MIN_PEAKS && return Inf32
    return Float32(quantile(offsets, 0.99) + MS1_WINDOW_FWHM_MULT * quantile(fwhms, 0.9))
end

"""Whether the calibration window h0 was too narrow: the peaks call for more than h0/2, or more than
MS1_WINDOW_CALIB_MAX_CUT_FRAC of them run into the trace ends."""
ms1_calib_too_narrow(need::Real, h0::Real, cut_frac::Real) = need > h0 / 2 || cut_frac > MS1_WINDOW_CALIB_MAX_CUT_FRAC

"""
    ms1_calib_action(n_peaks, n_candidates, need, h0, cut_frac, can_widen) -> Symbol

Decision after one calibration extraction: `:accept` the measured peaks, `:widen` the window and re-extract, or
`:fallback` to the iRT-based first guess. Too few measurable peaks means too narrow a window only when there were
enough candidate precursors to begin with; a file with few confident IDs has nothing to calibrate.
"""
function ms1_calib_action(n_peaks::Integer, n_candidates::Integer, need::Real, h0::Real, cut_frac::Real,
                          can_widen::Bool)::Symbol
    if n_peaks < MS1_WINDOW_CALIB_MIN_PEAKS
        return (n_candidates >= MS1_WINDOW_CALIB_MIN_PEAKS && can_widen) ? :widen : :fallback
    end
    return (ms1_calib_too_narrow(need, h0, cut_frac) && can_widen) ? :widen : :accept
end

"""Per-file WH smoothing scale (minutes) from the calibration peaks: max(MS1_SMOOTH_FWHM_MULT * median FWHM,
MS1_SMOOTH_MIN_CYCLES * median MS1 spacing). Medians over all calibration peaks (~3,000), so files with the same
chromatography get nearly the same S (a 5% change in S already moves ~2% of areas by > 0.5 log2)."""
function ms1_smoothing_scale(fwhms::AbstractVector{<:Real}, spacings::AbstractVector{<:Real})::Float32
    return Float32(max(MS1_SMOOTH_FWHM_MULT * median(fwhms), MS1_SMOOTH_MIN_CYCLES * median(spacings)))
end

"""
    ms1_calib_initial_half_width(spectra, search_context, ms_file_idx, ms2_rts) -> Float32 (minutes)

First-guess MS1 extraction half-width for the calibration pass: the run's fitted iRT tolerance (3-4 MAD of the
RT-vs-library-iRT residual) divided by the median local slope d(iRT)/d(RT) at the calibration precursors' MS2
times. It reflects how far RTs scatter around the alignment, which grows with gradient length like peak width does,
and on Astral 15 min DDA is ~3.7x the window the calibration finally chooses. The RT-binned tolerance is not used:
it is built from MS2 peak widths, which DDA's dynamic exclusion makes ~0.
"""
function ms1_calib_initial_half_width(spectra::MassSpecData, search_context::SearchContext, ms_file_idx::Int64,
                                      ms2_rts::AbstractVector{<:Real})
    irt_tol = getIrtErrors(search_context)[ms_file_idx]
    rt_irt_model = getRtIrtModel(search_context, ms_file_idx)
    h = 0.1f0
    slopes = Float32[abs((rt_irt_model(Float32(rt) + h) - rt_irt_model(Float32(rt) - h)) / (2f0 * h)) for rt in ms2_rts]
    slope = isempty(slopes) ? 1f0 : max(median(slopes), 0.01f0)
    ms1_rts = Float32[getRetentionTime(spectra, i) for i in 1:length(spectra) if getMsOrder(spectra, i) == 1]
    spacing = length(ms1_rts) < 2 ? 0f0 : Float32(median(diff(ms1_rts)))
    h0 = Float32(irt_tol / slope)
    return clamp(isfinite(h0) ? h0 : MS1_WINDOW_CALIB_MAX_HALF_WIDTH,
                 MS1_WINDOW_CALIB_MIN_CYCLES * spacing, MS1_WINDOW_CALIB_MAX_HALF_WIDTH)
end

"""
    ms1_calibration_rows(passing_psms) -> Vector{Int}

Rows of the most confident target PSMs used to calibrate a DDA file's MS1 quantification (window, smoothing scale,
MS1 mass-error model): `MS1_WINDOW_CALIB_N` of them, spread evenly over RT.
Tens of thousands of targets tie at q = 0, so "the top N by q-value" would be whichever tied rows come first (an
early-RT subset in practice); every target at or below the N-th q-value is taken and N of them spread over RT, so the
calibration sees the whole gradient (peaks broaden late).
"""
function ms1_calibration_rows(passing_psms::DataFrame)
    tg = findall(passing_psms.target)
    isempty(tg) && return Int[]
    n_calib = MS1_WINDOW_CALIB_N
    order = sortperm(passing_psms.qval[tg])
    sel = tg[order[1:min(n_calib, length(order))]]
    q_cut = passing_psms.qval[last(sel)]
    pool = filter(i -> passing_psms.qval[i] <= q_cut, tg)
    if length(pool) > n_calib
        sort!(pool, by = i -> passing_psms.rt[i])
        sel = pool[round.(Int, range(1, length(pool), length = n_calib))]
    end
    return sel
end

"""
    calibrate_ms1_window(spectra, passing_psms, rt_index, search_context, params, ms_file_idx; ms1_mass_err_model)
        -> (half_width::Float32, smooth_scale::Float32, fwhm_q90::Float32)
           (minutes; fwhm_q90 = 0 when the file falls back to the iRT-based window)

Per-run MS1 extraction half-width for DDA. MS1 peak width and the offset between a precursor's identifying MS2
scan and its MS1 apex depend on gradient and instrument (Astral 15 min: FWHM ~2 s; QE HF-X 150 min: ~15 s), and
an oversized window admits co-eluting peaks and enlarges the MS1 deconvolution.
MS1 traces are built with a wide window h0 (`ms1_calib_initial_half_width`) for the most confident target
precursors (`ms1_calibration_rows`), matched with `ms1_mass_err_model` (nothing = tuning's MS1 model); for each,
the trace is climbed
from the point nearest its MS2 scan to its own local apex (not the trace's global maximum, which in a wide window
is often a different peak), and the apex offset and FWHM (interpolated half-height crossings) are recorded.
The same peaks give the file's smoothing scale (`ms1_smoothing_scale`). With fewer than MS1_WINDOW_CALIB_MIN_PEAKS
measurable peaks (`ms1_calib_action`) the file falls back to the iRT-based h0 and S = MS1_SMOOTH_FALLBACK_H0_FRAC * h0.
"""
function calibrate_ms1_window(spectra::MassSpecData, passing_psms::DataFrame, rt_index::retentionTimeIndex,
                              search_context::SearchContext, params::IntegrateChromatogramSearchParameters,
                              ms_file_idx::Int64;
                              ms1_mass_err_model::Union{Nothing, AbstractMassErrorModel} = nothing)
    t_start = time()
    sel = ms1_calibration_rows(passing_psms)
    isempty(sel) && return (half_width = MS1_DDA_RT_TOL_DEFAULT, smooth_scale = 0f0, fwhm_q90 = 0f0)
    calib = passing_psms[sel, :]
    h0 = ms1_calib_initial_half_width(spectra, search_context, ms_file_idx, calib.rt)
    h0_init = h0
    calib_im = Dict{UInt32, Int}()
    band = 0
    let ims = getImScans(spectra)
        if ims !== nothing
            band = im_half_width_scans(CHROM_IM_BAND_K0, spectra, search_context, ms_file_idx)
            for i in 1:nrow(calib); calib_im[calib.precursor_idx[i]] = Int(ims[calib.scan_idx[i]]); end
        end
    end
    ms2_rt = Dict(calib.precursor_idx[i] => Float32(calib.rt[i]) for i in 1:nrow(calib))
    offsets = Float32[]; fwhms = Float32[]; spacings = Float32[]
    attempts = String[]
    action = :fallback
    for attempt in 0:MS1_WINDOW_CALIB_MAX_RETRIES
        if nrow(calib) < MS1_WINDOW_CALIB_MIN_PEAKS   # few confident IDs: nothing to calibrate, skip the extraction
            push!(attempts, "only $(nrow(calib)) confident targets")
            break
        end
        chroms, _ = extract_chromatograms(spectra, calib, rt_index, search_context, params, ms_file_idx, MS1CHROM();
                                          ms1_rt_tol = h0, ms1_mass_err_model = ms1_mass_err_model)
        chroms = collapse_ms1_im_slices(chroms, spectra, calib_im, band)   # mobility data: RT marginal over the band
        sort!(chroms, [:precursor_idx, :rt])
        empty!(offsets); empty!(fwhms); empty!(spacings)
        n_trunc = 0
        for sub in groupby(chroms, :precursor_idx)
            r = get(ms2_rt, first(sub.precursor_idx), nothing); r === nothing && continue
            o, w, dt, cut = _ms1_local_peak(sub.rt, sub.intensity, r)
            isfinite(o) || continue
            push!(offsets, abs(o)); push!(fwhms, w); push!(spacings, dt); n_trunc += cut
        end
        need = ms1_window_need(offsets, fwhms)
        cut_frac = isempty(offsets) ? 1.0 : n_trunc / length(offsets)
        action = ms1_calib_action(length(offsets), nrow(calib), need, h0, cut_frac,
                                  h0 < MS1_WINDOW_CALIB_MAX_HALF_WIDTH && attempt < MS1_WINDOW_CALIB_MAX_RETRIES)
        push!(attempts, "h0 +/-$(round(60h0, digits = 1)) s: n=$(length(offsets)), cut $(round(100cut_frac, digits = 1))%, " *
                        "need +/-$(isfinite(need) ? round(60need, digits = 1) : "?") s -> $action")
        action === :widen || break
        h0 = min(2f0 * h0, MS1_WINDOW_CALIB_MAX_HALF_WIDTH)
    end
    if action === :fallback
        S0 = MS1_SMOOTH_FALLBACK_H0_FRAC * h0_init
        msg = "DDA MS1 window: +/-$(round(60h0_init, digits = 1)) s, smoothing scale $(round(60S0, digits = 1)) s " *
              "from the iRT tolerance (too few calibration peaks: $(join(attempts, "; ")))"
        nrow(calib) < MS1_WINDOW_CALIB_MIN_PEAKS ? (@user_info msg) : (@user_warn msg)
        return (half_width = h0_init, smooth_scale = S0, fwhm_q90 = 0f0)
    end
    lo = MS1_WINDOW_MIN_CYCLES * median(spacings)
    w = clamp(ms1_window_need(offsets, fwhms), Float32(lo), h0)
    S = ms1_smoothing_scale(fwhms, spacings)
    fq90 = Float32(quantile(fwhms, 0.9))
    @user_info "DDA MS1 window: +/-$(round(60w, digits = 1)) s, smoothing scale $(round(60S, digits = 1)) s " *
               "(apex offset q99 $(round(60quantile(offsets, 0.99), digits = 1)) s, " *
               "FWHM median $(round(60median(fwhms), digits = 2)) / q90 $(round(60quantile(fwhms, 0.9), digits = 1)) s, " *
               "MS1 spacing $(round(60median(spacings), digits = 2)) s; " *
               "calibration from iRT tol: $(join(attempts, "; ")); $(round(time() - t_start, digits = 2)) s)"
    return (half_width = w, smooth_scale = S, fwhm_q90 = fq90)
end

"""Offset (min) of the local apex reached by climbing from the point nearest `ms2_rt`, its FWHM (min), the median
point spacing (min), and whether the peak is cut off (a half-height crossing would lie beyond the trace end), on a
3-point smoothed trace. NaN offset when there is no peak."""
function _ms1_local_peak(rt::AbstractVector{Float32}, y::AbstractVector{Float32}, ms2_rt::Float32)
    n = length(y)
    n < 5 && return (NaN32, NaN32, NaN32, false)
    ys = Vector{Float32}(undef, n)
    @inbounds for i in 1:n
        lo = max(1, i - 1); hi = min(n, i + 1)
        s = 0f0; for k in lo:hi; s += y[k]; end
        ys[i] = s / (hi - lo + 1)
    end
    a = argmin(abs.(rt .- ms2_rt))
    while true
        b = a
        a > 1 && ys[a - 1] > ys[b] && (b = a - 1)
        a < n && ys[a + 1] > ys[b] && (b = a + 1)
        b == a && break
        a = b
    end
    ys[a] > 0f0 || return (NaN32, NaN32, NaN32, false)
    h = ys[a] / 2
    l = a; while l > 1 && ys[l - 1] > h; l -= 1; end
    r = a; while r < n && ys[r + 1] > h; r += 1; end
    # interpolate the half-height crossings
    tl = l > 1 ? rt[l - 1] + (rt[l] - rt[l - 1]) * (h - ys[l - 1]) / max(ys[l] - ys[l - 1], eps(Float32)) : rt[l]
    tr = r < n ? rt[r] + (rt[r + 1] - rt[r]) * (ys[r] - h) / max(ys[r] - ys[r + 1], eps(Float32)) : rt[r]
    return (rt[a] - ms2_rt, tr - tl, median(diff(rt)), (l == 1) | (r == n))
end

#==========================================================
Second pass: candidacy re-centred on each precursor's MS1 apex
==========================================================#

"""
Whether DDA MS1 quantification runs a second extraction re-centred on the first pass's MS1 apices (default on;
`PIONEER_MS1_TWO_PASS=0` turns it off). Against one pass on eight DDA sets (Astral / HF-X / timsTOF, 50 ng and
250 pg): ratios as accurate or closer to truth everywhere (timsTOF 50 ng E. coli -1.46 -> -1.58, -1.14 -> -1.24),
CV within +/-0.005 except timsTOF 50 ng 5 min (0.189 -> 0.176); quantified precursors within 1%, except timsTOF 50 ng
(-2.8% / -3.8%).

In DIA the best-scoring MS2 scan sits near the elution apex, so windows centred on the PSM are centred on the peak.
In DDA the identifying MS2 scan is wherever the instrument triggered: up to ~1 FWHM from the MS1 apex in RT and
up to ~30 IM scans in mobility. The first-pass candidacy window is therefore wide (calibrated RT half-width =
q99 apex offset + 2 x q90 FWHM; +/-CHROM_IM_WINDOW_K0 in mobility), and a precursor holds a column in many slices
where its own peak is absent, competing for shared peaks.
"""
ms1_two_pass() = get(ENV, "PIONEER_MS1_TWO_PASS", "1") != "0"

"""Second-pass RT half-width as a multiple of the calibrated q90 FWHM (the first-pass formula without the apex-offset
term, which only covered the PSM-to-apex distance)."""
const MS1_PASS2_RT_FWHM_MULT = 2.0f0

"""
    ms1_second_pass(ms1_chromatograms, spectra, passing_psms, search_context, params, ms_file_idx; kwargs...)
        -> (chromatograms, seed_scan::Vector{UInt32}, apex_im)

Integrates the first-pass MS1 chromatograms, then re-extracts with each precursor centred on its first-pass apex
cell: RT window apex +/- MS1_PASS2_RT_FWHM_MULT x q90 FWHM (at most the first-pass half-width; CHROM_RHS_MULT times
that on the right for `rhs_extended` precursors, as in DIA), mobility window = the integration band
(+/-CHROM_IM_BAND_K0) around the apex slice, so a precursor is deconvolved where it is integrated. Precursors without
a first-pass apex keep their PSM-centred RT and seed. Returns the second-pass chromatograms (with cycle_idx / im_scan
on mobility data), the integration seed per PSM (the apex scan) and the seed's IM scan (nothing without mobility).
"""
function ms1_second_pass(ms1_chromatograms::DataFrame, spectra::MassSpecData, passing_psms::DataFrame,
                         search_context::SearchContext, params::IntegrateChromatogramSearchParameters, ms_file_idx::Int64;
                         ms1_mass_err_model, fwhm_q90::Float32, pass1_rt_tol::Float32, rhs_extended::Set{UInt32},
                         im_half_scans::Int, rt_scale::Float32, apex_im)
    t0 = time()
    n = nrow(passing_psms)
    area = zeros(Float32, n); best = zeros(UInt32, n); npts = zeros(UInt32, n)
    s0 = zeros(UInt32, n); s1 = zeros(UInt32, n); withheld = falses(n)
    sort_chromatograms_for_integration!(ms1_chromatograms, CombineTraces(params.min_fraction_transmitted))
    integrate_precursors(ms1_chromatograms, CombineTraces(params.min_fraction_transmitted), params.min_fraction_transmitted,
                         passing_psms[!, :precursor_idx], passing_psms[!, :scan_idx], area, best, npts, s0, s1, withheld;
                         λ = params.wh_smoothing_strength, im_half_scans = im_half_scans, rt_scale = rt_scale, apex_im = apex_im)
    seed = Vector{UInt32}(passing_psms[!, :scan_idx]); rt = Vector{Float32}(passing_psms[!, :rt])
    n_apex = 0
    @inbounds for i in 1:n
        if area[i] > 0f0 && best[i] != 0
            seed[i] = best[i]; rt[i] = Float32(getRetentionTime(spectra, Int(best[i]))); n_apex += 1
        end
    end
    pids = Vector{UInt32}(passing_psms[!, :precursor_idx])
    order = sortperm(rt)
    prec_mzs = getMz(getPrecursors(getSpecLib(search_context)))
    rt_index2 = buildRtIndex(rt[order], Float32[prec_mzs[p] for p in pids[order]], pids[order], 0.1)
    h2 = fwhm_q90 > 0f0 ? min(MS1_PASS2_RT_FWHM_MULT * fwhm_q90, pass1_rt_tol) : pass1_rt_tol
    psms2 = DataFrame(precursor_idx = pids, rt = rt, scan_idx = seed)
    chroms, _ = extract_chromatograms(spectra, psms2, rt_index2, search_context, params, ms_file_idx, MS1CHROM();
                                      ms1_rt_tol = h2, ms1_mass_err_model = ms1_mass_err_model,
                                      ms1_im_window_k0 = CHROM_IM_BAND_K0, rhs_extended = rhs_extended)
    ims = getImScans(spectra)
    apex_im2 = nothing
    if ims !== nothing && nrow(chroms) > 0
        cyc = getCycleIdxs(spectra)
        chroms[!, :cycle_idx] = UInt32[UInt32(cyc[s]) for s in chroms.scan_idx]
        chroms[!, :im_scan] = UInt16[UInt16(ims[s]) for s in chroms.scan_idx]
        apex_im2 = UInt16[UInt16(ims[s]) for s in seed]
    end
    @user_info "DDA MS1 second pass: $(n_apex) of $(n) precursors re-centred on their MS1 apex; RT +/-$(round(60h2, digits = 1)) s " *
               "(first pass +/-$(round(60pass1_rt_tol, digits = 1)) s, $(length(rhs_extended)) right-tail extended)" *
               (ims === nothing ? "" : ", mobility +/-$(CHROM_IM_BAND_K0) 1/K0 (first pass +/-$(CHROM_IM_WINDOW_K0))") *
               "; $(nrow(chroms)) points (first pass $(nrow(ms1_chromatograms))); $(round(time() - t0, digits = 2)) s"
    return chroms, seed, apex_im2
end

#==========================================================
DIA: MS1 quantification next to MS2 (maxLFQ.ms1_quantification)
==========================================================#

"""Fallback MS1 candidacy half-width for DIA, in MS2 integration widths (q90), when the file has no per-RT-bin
chromatogram tolerance."""
const DIA_MS1_WINDOW_WIDTHS = 1.0f0

"""
    dia_ms1_quantification(spectra, passing_psms, rt_index, search_context, params, ms_file_idx)
        -> (area::Vector{Float32}, withheld::Vector{Bool})

MS1 peak area of every passing precursor of a DIA file, next to its MS2 `peak_area` (which, with `new_best_scan` and
the integration bounds, is left untouched: MBR reads them), and whether the 5x baseline rule withheld it. Reuses the
DDA MS1 machinery -- per-file MS1 mass model, M0-M2 isotope-envelope deconvolution per MS1 scan, chromatogram
integration (RT x IM on mobility data) -- but treats each MS1 trace like the precursor's MS2 trace: candidacy is
centred on the MS2 apex (`new_best_scan`, else the PSM scan) with the MS2 chromatogram window (per-RT-bin tolerance
from the MS2 peak widths), a mobility band of +/-CHROM_IM_BAND_K0, and window-relative smoothing (`rt_scale = 0`), so
no MS1-specific smoothing constant. Chosen on 250 pg three-proteome data (timsTOF diaPASEF and Astral, 5 min) over
MS1-calibrated windows, the DDA two-pass rule, MS2-bounds summation and nearest-maximum seeding: as accurate as the
best alternative on Astral with the lowest CV, and the least compressed MS1 ratios on timsTOF.
"""
function dia_ms1_quantification(spectra::MassSpecData, passing_psms::DataFrame, rt_index, search_context,
                                params, ms_file_idx::Integer)
    n = nrow(passing_psms)
    area = zeros(Float32, n)
    # true where QUANT_MIN_AREA_SURVIVING_RATIO withheld the MS1 area (1D integration only; the 2D timsTOF
    # integrator has no such rule)
    withheld = zeros(Bool, n)
    n == 0 && return area, withheld
    t0 = time()
    # MS2 apex (the PSM scan when there was no MS2 peak) and integration widths (minutes)
    apex = UInt32[s > 0 ? s : UInt32(passing_psms.scan_idx[i]) for (i, s) in enumerate(passing_psms.new_best_scan)]
    apex_rt = Float32[Float32(getRetentionTime(spectra, Int(s))) for s in apex]
    widths = Float32[]
    for i in 1:n
        a, b = passing_psms.integration_start_scan[i], passing_psms.integration_stop_scan[i]
        (a > 0 && b > a) || continue
        push!(widths, Float32(getRetentionTime(spectra, Int(b)) - getRetentionTime(spectra, Int(a))))
    end
    # one time per MS1 frame: on mobility data each IM slice of a frame carries its own, slightly later, RT
    ms1_idx = [i for i in 1:length(spectra) if getMsOrder(spectra, i) == 1]
    if getImScans(spectra) !== nothing && !isempty(ms1_idx)
        cyc = getCycleIdxs(spectra)
        ms1_idx = ms1_idx[[k == 1 || cyc[ms1_idx[k]] != cyc[ms1_idx[k - 1]] for k in eachindex(ms1_idx)]]
    end
    ms1_rts = Float32[Float32(getRetentionTime(spectra, i)) for i in ms1_idx]
    spacing = length(ms1_rts) > 1 ? median(diff(ms1_rts)) : 0f0
    w90 = isempty(widths) ? MS1_DDA_RT_TOL_DEFAULT : quantile(widths, 0.9)
    rt_tol = max(DIA_MS1_WINDOW_WIDTHS * w90, 3 * spacing)

    ms1_mem = fit_ms1_quant_mass_model(spectra, passing_psms[ms1_calibration_rows(passing_psms), :],
                                       search_context, ms_file_idx)
    reset_ms1_solve_stats!()
    # candidacy centred on the MS2 apex (extract_chromatograms reads :precursor_idx, :rt and, for mobility, :scan_idx)
    centred = DataFrame(precursor_idx = passing_psms.precursor_idx, rt = apex_rt, scan_idx = apex)
    chroms, _ = extract_chromatograms(spectra, centred, rt_index, search_context, params, ms_file_idx, MS1CHROM();
                                      ms1_rt_tol = Float32(rt_tol), ms1_mass_err_model = ms1_mem,
                                      ms1_rt_binned = haskey(getRtTolerances(search_context), ms_file_idx) ?
                                                      getRtTolerance(search_context, Int64(ms_file_idx)) : nothing,
                                      ms1_im_window_k0 = CHROM_IM_BAND_K0)
    if nrow(chroms) == 0
        @user_info "DIA MS1 quantification: no MS1 signal extracted ($(length(ms1_rts)) MS1 scans); peak_area_ms1 left missing"
        return area, withheld
    end
    im_half = 0; apex_im = nothing
    let ims = getImScans(spectra)
        if ims !== nothing
            im_half = im_half_width_scans(CHROM_IM_BAND_K0, spectra, search_context, ms_file_idx)
            apex_im = UInt16[UInt16(ims[s]) for s in apex]
            cyc = getCycleIdxs(spectra)
            chroms[!, :cycle_idx] = UInt32[UInt32(cyc[s]) for s in chroms.scan_idx]
            chroms[!, :im_scan] = UInt16[UInt16(ims[s]) for s in chroms.scan_idx]
        end
    end
    sort_chromatograms_for_integration!(chroms, CombineTraces(params.min_fraction_transmitted))
    # MS1-owned outputs: the MS2 apex, points and bounds in passing_psms stay as integrated from MS2
    best = zeros(UInt32, n); pts = zeros(UInt32, n); lo = zeros(UInt32, n); hi = zeros(UInt32, n)
    integrate_precursors(chroms, CombineTraces(params.min_fraction_transmitted), params.min_fraction_transmitted,
                         passing_psms.precursor_idx, apex, area, best, pts, lo, hi, withheld;
                         λ = params.wh_smoothing_strength, im_half_scans = im_half, rt_scale = 0f0,
                         apex_im = apex_im)
    st = MS1_SOLVE_STATS
    @user_info "DIA MS1 quantification: $(count(>(0f0), area)) of $n precursors quantified; window +/-" *
               "$(round(60rt_tol, digits = 1)) s (q90 MS2 width $(round(60w90, digits = 1)) s); $(st.scans[]) MS1 " *
               "scans, $(round(st.cols[] / max(st.scans[], 1), digits = 1)) cols/scan; $(nrow(chroms)) points; " *
               "$(round(time() - t0, digits = 2)) s"
    return area, withheld
end
