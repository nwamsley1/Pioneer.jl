# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
# Licensed under AGPL v3+; see LICENSE.

#==========================================================
MS1-level quantitation support.

DDA gives no MS2 chromatogram: under dynamic exclusion a precursor is
fragmented in ~1 MS2 scan (measured on LFQBench PXD028735: n_scans median 1),
so the MS2 trace has no extent and integrates to exactly zero. MS1 is the only
place a DDA precursor has a real chromatogram.

Two pieces live here:

  1. `collect_ms1_window_precursors!` — the MS1 analogue of
     `collect_rt_window_precursors!`. The MS2 version prunes candidates by the
     quadrupole isolation window; at MS1 there is no isolation window, so the
     RT window is the ONLY thing bounding the candidate set.

  2. `Ms1MzGroups` — collapses exactly-collinear candidates before the solve.
     At MS1 there is no fragmentation, so two precursors with the same m/z AND
     the same charge produce identical design-matrix columns and their
     individual weights are unidentifiable. Pioneer's decoys are sequence
     permutations, so every decoy shares its target's elemental composition,
     m/z, charge and isotope envelope — i.e. every target/decoy pair is an
     exact duplicate column. Those must be collapsed structurally; no penalty
     can separate identical columns. Near-collinear columns (distinct peptides
     within the MS1 ppm tolerance) are left to the ridge term instead.

See ms1_quant_dda_plan.pdf for the full design and test plan.
==========================================================#

"""
    collect_ms1_window_precursors!(precs_temp, rt_index, rt_start_idx, rt_stop_idx,
                                   precursors_passing, precursor_rt_map, scan_rt,
                                   rt_binned_tol, rt_tol_fallback) -> Int

Collect candidate precursors for one MS1 scan into `precs_temp[1:n]`, returning `n`.

The MS2 counterpart (`collect_rt_window_precursors!`) additionally filters on the
quadrupole transmission window, the transmitted-isotope set, and
`min_fraction_transmitted`. None of those exist at MS1 — the whole precursor m/z
range is observed — so this walks the full m/z extent of each RT bin.
"""
function collect_ms1_window_precursors!(
    precs_temp::Vector{UInt32},
    rt_index::retentionTimeIndex,
    rt_start_idx::Int64, rt_stop_idx::Int64,
    precursors_passing::Union{Set{UInt32}, Nothing},
    precursor_rt_map::Union{Dict{UInt32, Float32}, Nothing},
    scan_rt::Float32,
    rt_binned_tol::Union{RTBinnedTolerance, Nothing},
    rt_tol_fallback::Float32,
    rt_tol_floor::Float32 = 0f0)

    size = 0
    has_rt_filter = precursor_rt_map !== nothing

    for rt_bin_idx in rt_start_idx:rt_stop_idx
        precs = rt_index.rt_bins[rt_bin_idx].prec

        # No quadrupole bounds at MS1: take the whole bin rather than the
        # getPrecursorWindowRange sub-range the MS2 path can use.
        for i in eachindex(precs)
            prec_idx = first(precs[i])
            (!isnothing(precursors_passing) && prec_idx ∉ precursors_passing) && continue

            if has_rt_filter
                prec_rt_val = get(precursor_rt_map, prec_idx, NaN32)
                if !isnan(prec_rt_val)
                    prec_tol = rt_binned_tol !== nothing ?
                        get_rt_tol(rt_binned_tol, prec_rt_val) : rt_tol_fallback
                    # Floor: the DIA-derived tolerance is 0 on DDA data.
                    prec_tol = max(prec_tol, rt_tol_floor)
                    abs(scan_rt - prec_rt_val) > prec_tol && continue
                end
            end

            size += 1
            if size > length(precs_temp)
                append!(precs_temp, Vector{UInt32}(undef, length(precs_temp)))
            end
            precs_temp[size] = prec_idx
        end
    end
    return size
end

#==========================================================
m/z + charge grouping
==========================================================#

"""
Rounding precision for the m/z grouping key: m/z is rounded to `1/MS1_MZ_GROUP_SCALE`.

1e-5 Th is ~0.02 ppm at m/z 500 — far finer than the fitted MS1 tolerance
(~±9 ppm), so this collapses ONLY exact duplicates (target/decoy pairs and true
isobars). Near-collinear-but-distinct precursors stay in separate columns and
are handled by the ridge penalty. That division of labour is deliberate: see the
module comment above.
"""
const MS1_MZ_GROUP_SCALE = 100_000.0f0

#==========================================================
DDA retention-time tolerance
==========================================================#

"""
Default half-width (minutes) of the DDA extraction window: 15 s.

The DIA-derived RT tolerance cannot be reused. `RTBinnedTolerance` is built in
`prescore_aggregation.jl` from the observed MS2 chromatographic FWHM:

    fwhm_half  = 2.0 * median(selected_fwhm)
    cycle_half = min_n_scans * rt_cycle_time
    rt_tols[bin] = max(fwhm_half, cycle_half)

Both terms assume many MS2 scans per precursor. In DDA there is ~1, so the
measured FWHM is 0 and the tolerance collapses to 0.0 across the entire
gradient — verified on LFQBench PXD028735, where mid-gradient scans with 138
precursors in the RT bin returned 0 candidates purely because tol was 0.

Two DDA-specific reasons the window must also be generous rather than merely
non-zero:

  1. The MS2 scan used for identification is NOT acquired at the chromatographic
     apex — DDA triggers on threshold crossing, not on peak maximum — so the
     anchor RT is offset from the true apex by an unknown amount.
  2. That offset is one-sided and variable, so a window centred on the anchor
     must be wide enough to contain the apex on either side.

A fixed floor is a placeholder. Determining this automatically for DDA (e.g.
from the MS1 XIC width of confident IDs, once MS1 traces exist) is future work.
Override with PIONEER_DDA_RT_TOL_MIN (minutes).
"""
const MS1_DDA_RT_TOL_DEFAULT = 0.25f0

function ms1_dda_rt_tol_floor()::Float32
    s = get(ENV, "PIONEER_DDA_RT_TOL_MIN", "")
    isempty(s) && return MS1_DDA_RT_TOL_DEFAULT
    v = tryparse(Float32, s)
    return (v === nothing || v < 0f0) ? MS1_DDA_RT_TOL_DEFAULT : v
end

@inline function ms1_group_key(prec_mz::Float32, prec_charge::UInt8)::UInt64
    # Pack (rounded m/z, charge) into one integer key. m/z <= ~2e3 => the scaled
    # value needs <= 28 bits, so an 8-bit charge shift is safe in UInt64.
    mz_key = UInt64(round(Int64, prec_mz * MS1_MZ_GROUP_SCALE))
    return (mz_key << 8) | UInt64(prec_charge)
end

"""
    Ms1MzGroups

Reusable per-scan scratch mapping candidate precursors onto shared design-matrix
columns. One column per distinct `(rounded m/z, charge)` key.

Fields (all sized to the candidate count, reused across scans):
- `key_to_col`  : packed key -> column index
- `col_of_cand` : candidate slot j -> its column index
- `col_rep`     : column -> representative precursor_idx (the column actually solved)
- `n_cols`      : number of distinct columns

Charge is part of the key. The historical `MzGroupingMap` keyed on m/z alone,
which is wrong: two precursors with equal m/z but different charge have isotope
spacing 1.00336/z, so their M1/M2 rows land at different m/z and their columns
are NOT identical. Merging them would introduce a real error.
"""
mutable struct Ms1MzGroups
    key_to_col::Dict{UInt64, UInt32}
    col_of_cand::Vector{UInt32}
    col_rep::Vector{UInt32}
    n_cols::Int
end

function Ms1MzGroups(initial_capacity::Int = 4096)
    return Ms1MzGroups(
        Dict{UInt64, UInt32}(),
        Vector{UInt32}(undef, initial_capacity),
        Vector{UInt32}(undef, initial_capacity),
        0,
    )
end

"""
    build_ms1_mz_groups!(groups, precs_temp, n_cands, prec_mzs, prec_charges) -> Int

Assign each candidate in `precs_temp[1:n_cands]` to a column keyed on
`(rounded m/z, charge)`, returning the number of distinct columns.

After this, `groups.col_rep[1:n_cols]` holds the representative precursor for each
column — those are the only precursors that need design-matrix rows — and
`groups.col_of_cand[j]` maps candidate `j` back to its column.
"""
function build_ms1_mz_groups!(
    groups::Ms1MzGroups,
    precs_temp::Vector{UInt32},
    n_cands::Int,
    prec_mzs::AbstractArray{Float32},
    prec_charges::AbstractArray{UInt8},
)
    empty!(groups.key_to_col)
    groups.n_cols = 0

    if n_cands > length(groups.col_of_cand)
        resize!(groups.col_of_cand, n_cands)
    end
    if n_cands > length(groups.col_rep)
        resize!(groups.col_rep, n_cands)
    end

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
Dimension probe (risk 1 in the plan)

At MS2 the isolation window prunes candidates hard. At MS1 nothing does, so the
RT window alone bounds the design matrix. Before building the matrix we measure
what it would actually be, because if it is 100x the MS2 case the whole approach
needs rethinking. Enabled with PIONEER_MS1_DIMS=1.
==========================================================#

const MS1_DIM_STATS = Dict{Symbol, Int}(
    :n_scans => 0, :n_cand_total => 0, :n_cand_max => 0,
    :n_col_total => 0, :n_col_max => 0,
    :n_visited => 0, :n_empty_bins => 0, :n_zero_cand => 0,
)
# build_chromatograms runs one task per thread; the probe accumulates per-thread
# and merges here.
const MS1_DIM_LOCK = ReentrantLock()

function reset_ms1_dim_stats!()
    for k in keys(MS1_DIM_STATS); MS1_DIM_STATS[k] = 0; end
    return nothing
end

function report_ms1_dim_stats()
    n = max(1, MS1_DIM_STATS[:n_scans])
    @user_info string(
        "[MS1-DIMS] scans=", MS1_DIM_STATS[:n_scans],
        "  candidates/scan mean=", round(MS1_DIM_STATS[:n_cand_total]/n, digits=1),
        " max=", MS1_DIM_STATS[:n_cand_max],
        "  columns/scan mean=", round(MS1_DIM_STATS[:n_col_total]/n, digits=1),
        " max=", MS1_DIM_STATS[:n_col_max],
        "  visited=", MS1_DIM_STATS[:n_visited],
        " empty_bin_range=", MS1_DIM_STATS[:n_empty_bins],
        " zero_cand=", MS1_DIM_STATS[:n_zero_cand],
        "  collapse=", MS1_DIM_STATS[:n_cand_total] == 0 ? "n/a" :
            string(round(100*(1 - MS1_DIM_STATS[:n_col_total]/MS1_DIM_STATS[:n_cand_total]), digits=1), "%"),
    )
    return nothing
end

"""
    distribute_ms1_group_weights!(precursor_weights, group_weights, groups,
                                  precs_temp, n_cands)

Write each column's fitted weight back to EVERY candidate in that column.

The weight is duplicated, not split. Splitting a group weight n ways would
systematically underestimate every real peptide in proportion to how many
isobars happen to share its m/z — an abundance bias driven by library
composition rather than by chemistry. Duplication gives each candidate the
correct conditional quantity ("the abundance IF this is the peptide"); the MS2
identification decides which member that is, and only the identified member
survives FDR.

Consequence: a target and its decoy receive identical MS1 weights by
construction. That is safe for FDR — the value cannot discriminate, so it adds
no artificial separation — but it does mean MS1 quant can never contribute a
discriminating feature to scoring. It is a quantitation mechanism only.
"""
function distribute_ms1_group_weights!(
    precursor_weights::Vector{Float32},
    group_weights::AbstractVector{Float32},
    groups::Ms1MzGroups,
    precs_temp::Vector{UInt32},
    n_cands::Int,
)
    @inbounds for j in 1:n_cands
        col = groups.col_of_cand[j]
        iszero(col) && continue
        prec_idx = precs_temp[j]
        if prec_idx <= length(precursor_weights) && col <= length(group_weights)
            precursor_weights[prec_idx] = group_weights[col]
        end
    end
    return nothing
end

#==========================================================
Design matrix
==========================================================#

"""Number of precursor isotopes modelled per column (M0, M1, M2)."""
const MS1_N_ISOTOPES = 3

"""
Ridge penalty for the MS1 solve.

Grouping removes exactly-collinear columns; this handles the NEAR-collinear ones
(distinct peptides inside the MS1 ppm tolerance, whose columns share support and
differ only slightly in their isotope ratios). Override with PIONEER_MS1_LAMBDA.
"""
const MS1_LAMBDA_DEFAULT = 1.0f0

function ms1_lambda()::Float32
    s = get(ENV, "PIONEER_MS1_LAMBDA", "")
    isempty(s) && return MS1_LAMBDA_DEFAULT
    v = tryparse(Float32, s)
    return (v === nothing || v < 0f0) ? MS1_LAMBDA_DEFAULT : v
end

"""
Huber delta for the MS1 solve — set high enough to be effectively OLS.

The MS2 chromatogram default is 300, calibrated against Pioneer's MS2
deconvolution scale. MS1 peak intensities on this data run 2.7e3 - 1.7e9
(median 6e5), so delta=300 puts EVERY residual in the Huber *linear* regime,
where the second derivative vanishes and the per-column Newton/bisection solve
stalls. Measured consequence: fitted weights collapsed into 16.8-148.6 while the
M0 intensity driving them spanned six orders of magnitude.

delta=1e10 exceeds the largest observed residual, so the loss is purely
quadratic everywhere and this behaves as OLS with an L2 penalty. That is the
intended starting point; a properly scaled delta (or a real MS1 Huber tuning
step) is a later optimisation. Override with PIONEER_MS1_HUBER_DELTA.
"""
const MS1_HUBER_DELTA_DEFAULT = 1.0f10

function ms1_huber_delta()::Float32
    s = get(ENV, "PIONEER_MS1_HUBER_DELTA", "")
    isempty(s) && return MS1_HUBER_DELTA_DEFAULT
    v = tryparse(Float32, s)
    return (v === nothing || v <= 0f0) ? MS1_HUBER_DELTA_DEFAULT : v
end

"""
    ms1_huber_solver() -> HuberSolver

The MS2 chromatogram solver runs unregularised (`lambda=0`, `NoNorm()`); MS1
needs its own instance with an L2 term. OLSSolver/PoissonMMSolver carry no
penalty fields at all, so Huber is currently the only option here (plan T7).
"""
function ms1_huber_solver()
    return HuberSolver(
        ms1_huber_delta(),  # delta — effectively OLS, see above
        ms1_lambda(),   # lambda
        Int64(50),      # max_iter_newton
        Int64(100),     # max_iter_bisection
        10.0f0,         # accuracy_newton
        10.0f0,         # accuracy_bisection
        L2Norm(),
    )
end

"""
    build_ms1_design_matrix!(Hs, groups, n_cols, scan_mz, scan_int, n_peaks,
                             prec_mzs, prec_charges, prec_sulfurs, iso_splines,
                             mass_err_model) -> Int

Fill `Hs` (CSC) with one column per m/z group and up to `MS1_N_ISOTOPES` rows per
column, returning the number of non-zeros.

Column i, isotope j contributes at m/z = `prec_mz + j*C13_C12_MASS_DIFF/charge`
with value `iso_splines(sulfur, j, mono_mass)` — the sulfur-aware Goldfarb
splines, i.e. NOT plain averagine. Row index is the observed peak index, so
unmatched peaks are simply empty rows (residual 0, no loss contribution).

Unlike the MS2 path there is no quadrupole transmission factor: at MS1 the whole
precursor m/z range is observed.
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
    m0_obs::Vector{Float32} = Float32[],
)
    needed = n_cols * MS1_N_ISOTOPES
    if needed + 1 > length(Hs.rowval)
        grow = needed + 1 - length(Hs.rowval)
        append!(Hs.rowval, zeros(UInt32, grow))
        append!(Hs.nzval,  zeros(Float32, grow))
        append!(Hs.x,      zeros(Float32, grow))
        append!(Hs.meta,   fill(UInt8(1), grow))
        append!(Hs.rank,   zeros(UInt8, grow))
    end
    if n_cols + 2 > length(Hs.colptr)
        append!(Hs.colptr, zeros(UInt32, n_cols + 2 - length(Hs.colptr)))
    end

    want_m0 = !isempty(m0_obs)
    if want_m0 && n_cols > length(m0_obs)
        resize!(m0_obs, n_cols)
    end

    nnz = 0
    @inbounds for col in 1:n_cols
        Hs.colptr[col] = UInt32(nnz + 1)
        want_m0 && (m0_obs[col] = 0f0)

        prec_idx = groups.col_rep[col]
        prec_mz = prec_mzs[prec_idx]
        charge  = prec_charges[prec_idx]
        sulfur  = prec_sulfurs[prec_idx]
        mono_mass = (prec_mz - Float32(PROTON)) * Float32(charge)
        inv_charge = 1f0 / Float32(charge)

        for j in 0:(MS1_N_ISOTOPES - 1)
            target = prec_mz + Float32(j) * C13_C12_MASS_DIFF_F32 * inv_charge
            lo, hi = getMzBounds(mass_err_model, target)

            # Peaks are m/z-sorted: first peak >= lo.
            p = searchsortedfirst(view(scan_mz, 1:n_peaks), lo)
            (p > n_peaks || Float32(scan_mz[p]) > hi) && continue

            # Within tolerance take the most intense peak.
            best_p = p; best_i = Float32(scan_int[p])
            q = p + 1
            while q <= n_peaks && Float32(scan_mz[q]) <= hi
                iq = Float32(scan_int[q])
                if iq > best_i
                    best_i = iq; best_p = q
                end
                q += 1
            end

            alpha = iso_splines(Int64(min(sulfur, UInt8(5))), j, mono_mass)
            alpha <= 0f0 && continue

            want_m0 && j == 0 && (m0_obs[col] = best_i)

            nnz += 1
            Hs.rowval[nnz] = UInt32(best_p)
            Hs.nzval[nnz]  = Float32(alpha)
            Hs.x[nnz]      = best_i
            Hs.meta[nnz]   = pack_meta(true, UInt8(j))
            Hs.rank[nnz]   = UInt8(0)
        end
    end
    Hs.colptr[n_cols + 1] = UInt32(nnz + 1)

    Hs.n_vals = nnz
    Hs.n = n_cols
    Hs.m = n_peaks
    return nnz
end

#==========================================================
MS1 scan loop
==========================================================#

"""
    build_chromatograms(..., ::MS1CHROM) -> DataFrame

MS1 counterpart of the MS2 `build_chromatograms`. Mirrors its RT-window logic
exactly (same `rt_tol_local` derivation, same RT-bin walk) but drops everything
tied to the isolation window, which does not exist at MS1.

Currently implements candidate collection + (m/z, charge) grouping and, under
PIONEER_MS1_DIMS=1, records the design-matrix dimensions that the solve would
face. Returns an empty chromatogram frame until the matrix build lands.
"""
function build_chromatograms(
    spectra::MassSpecData,
    scan_range::Vector{Int64},
    precursors_passing::Set{UInt32},
    precursor_rt_map::Dict{UInt32, Float32},
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64,
    ::MS1CHROM;
    scan_tic::Union{Nothing, Vector{Float32}} = nothing,
)
    collect_dims = get(ENV, "PIONEER_MS1_DIMS", "0") == "1"

    precs_temp = getPrecIds(search_data)
    groups = Ms1MzGroups()

    Hs = SparseArrayFused(UInt32(5000))
    weights = zeros(Float32, 8192)
    residuals = zeros(Float32, 8192)
    colnorm2 = zeros(Float32, 8192)
    solver = ms1_huber_solver()
    iso_splines = getIsoSplines(search_data)
    ms1_mass_err = getMs1MassErrorModel(search_context, ms_file_idx)

    m0_obs = zeros(Float32, 8192)
    chrom_rt = Float32[]; chrom_int = Float32[]; chrom_m0 = Float32[]
    chrom_scan = UInt32[]; chrom_prec = UInt32[]

    irt_tol = getIrtErrors(search_context)[ms_file_idx]
    has_rt_tol = haskey(getRtTolerances(search_context), ms_file_idx)
    rt_binned_tol = has_rt_tol ? getRtTolerance(search_context, ms_file_idx) : nothing
    rt_irt_model = getRtIrtModel(search_context, ms_file_idx)
    rt_tol_floor = ms1_dda_rt_tol_floor()

    precursors = getPrecursors(getSpecLib(search_context))
    prec_mz_arr = getMz(precursors)
    prec_charge_arr = getCharge(precursors)
    prec_sulfur_arr = getSulfurCount(precursors)

    n_scans = 0; n_cand_total = 0; n_cand_max = 0; n_col_total = 0; n_col_max = 0
    n_visited = 0; n_empty_bins = 0; n_zero_cand = 0

    for scan_idx in scan_range
        ((scan_idx < 1) | (scan_idx > length(spectra))) && continue
        getMsOrder(spectra, scan_idx) == 1 || continue
        n_visited += 1

        rt = getRetentionTime(spectra, scan_idx)

        # Same local RT tolerance derivation as the MS2 path.
        if rt_binned_tol !== nothing
            rt_tol_local = get_rt_tol(rt_binned_tol, Float32(rt))
        else
            h = 0.1f0
            local_slope = abs((rt_irt_model(rt + h) - rt_irt_model(rt - h)) / (2f0 * h))
            rt_tol_local = irt_tol / max(local_slope, 0.01f0)
        end
        rt_tol_local = max(rt_tol_local, rt_tol_floor)

        rt_bin_start = max(searchsortedfirst(rt_index.rt_bins, rt - rt_tol_local,
                                             lt=(r,x)->r.lb<x) - 1, 1)
        rt_bin_stop  = min(searchsortedlast(rt_index.rt_bins, rt + rt_tol_local,
                                            lt=(x,r)->r.ub>x) + 1, length(rt_index.rt_bins))

        n_cands = collect_ms1_window_precursors!(
            precs_temp, rt_index, rt_bin_start, rt_bin_stop,
            precursors_passing, precursor_rt_map, Float32(rt),
            rt_binned_tol, rt_tol_local, rt_tol_floor)

        if rt_bin_stop < rt_bin_start
            n_empty_bins += 1
        end
        if n_cands == 0
            n_zero_cand += 1
            continue
        end

        n_cols = build_ms1_mz_groups!(groups, precs_temp, n_cands,
                                      prec_mz_arr, prec_charge_arr)

        scan_mz  = getMzArray(spectra, scan_idx)
        scan_int = getIntensityArray(spectra, scan_idx)
        n_peaks  = length(scan_mz)

        nnz = build_ms1_design_matrix!(
            Hs, groups, n_cols, scan_mz, scan_int, n_peaks,
            prec_mz_arr, prec_charge_arr, prec_sulfur_arr,
            iso_splines, ms1_mass_err, m0_obs)

        if nnz > 0
            if n_cols > length(weights)
                resize!(weights, n_cols + 1000)
                resize!(colnorm2, n_cols + 1000)
            end
            @inbounds for c in 1:n_cols
                weights[c] = 0f0
            end
            solve_deconvolution!(solver, Hs, residuals, weights, colnorm2,
                                 getMu(search_data), getObserved(search_data),
                                 params.max_iter_outer, params.max_diff)

            rt32 = Float32(rt); scan32 = UInt32(scan_idx)
            @inbounds for j in 1:n_cands
                col = groups.col_of_cand[j]
                push!(chrom_rt, rt32)
                push!(chrom_int, iszero(col) ? 0f0 : max(weights[col], 0f0))
                push!(chrom_m0, iszero(col) ? 0f0 : m0_obs[col])
                push!(chrom_scan, scan32)
                push!(chrom_prec, precs_temp[j])
            end
        end

        if collect_dims
            n_scans += 1
            n_cand_total += n_cands; n_cand_max = max(n_cand_max, n_cands)
            n_col_total  += n_cols;  n_col_max  = max(n_col_max, n_cols)
        end
    end

    if collect_dims
        lock(MS1_DIM_LOCK) do
            MS1_DIM_STATS[:n_scans]      += n_scans
            MS1_DIM_STATS[:n_cand_total] += n_cand_total
            MS1_DIM_STATS[:n_cand_max]    = max(MS1_DIM_STATS[:n_cand_max], n_cand_max)
            MS1_DIM_STATS[:n_col_total]  += n_col_total
            MS1_DIM_STATS[:n_col_max]     = max(MS1_DIM_STATS[:n_col_max], n_col_max)
            MS1_DIM_STATS[:n_visited]    += n_visited
            MS1_DIM_STATS[:n_empty_bins] += n_empty_bins
            MS1_DIM_STATS[:n_zero_cand]  += n_zero_cand
        end
    end

    return DataFrame(
        rt = chrom_rt,
        intensity = chrom_int,
        m0_intensity = chrom_m0,
        scan_idx = chrom_scan,
        precursor_idx = chrom_prec,
    )
end
