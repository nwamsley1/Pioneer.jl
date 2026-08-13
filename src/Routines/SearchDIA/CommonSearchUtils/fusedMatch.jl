# Fused per-precursor match+score+build.
#
# Loaded AFTER selectTransitions/ (which defines PrecEstimation and
# getFragIsotopes!). Uses the low-level primitives (FusedScratch,
# bsearch_hybrid, finalize_column!, etc.) defined in fusedScan.jl.

#==========================================================
FusedSearchKind — trait that factors the per-strategy variation points
out of the shared `run_fused!` body via multiple dispatch.

Variation points (all small, @inline'd — the compiler emits one
specialization of run_fused! per concrete kind):
  - n_iso_passes(kind) : outer iso-pass count (1 for normal searches, 3
    for Quad which fits per-isotope weights by one-hot transmission).
  - max_frag_rank_for(kind) : rank filter threshold.
  - prec_estimation_for(kind) : PartialPrecCapture or FullPrecCapture.
  - match_column_id(kind, prec_idx, iso_pass) : how (prec_idx, iso_pass)
    maps to a design-matrix column. Quad encodes iso in the column key so
    deconvolution fits a separate weight per isotope.

Concrete kinds land in separate commits as search methods migrate:
  - FusedStandard   : MainSearch, NceTuning (this step).
  - FusedRTIndexed  : IntegrateChromatogramsSearch (step 2).
  - FusedQuadEst    : QuadTuning (step 3).
==========================================================#

abstract type FusedSearchKind end

"""
    FusedStandard{P<:PrecEstimation}(prec_estimation, max_frag_rank)

Default fused search kind — walks `precursors_passed[prec_range]`, uses
the full quadrupole transmission per precursor, one pass per fragment
× isotope. Used by MainSearch (and NceTuning when/if migrated).
"""
struct FusedStandard{P<:PrecEstimation} <: FusedSearchKind
    prec_estimation::P
    max_frag_rank::UInt8
end

"""
    FusedRTIndexed{P<:PrecEstimation}(prec_estimation, max_frag_rank)

Fused kind for IntegrateChromatogramsSearch. The RT-bin + quad-window +
per-precursor-RT + min_fraction_transmitted filters are all applied
**upstream** (`collect_rt_window_precursors!`), producing a pre-filtered
`precs_temp::Vector{UInt32}`. `run_fused!` walks that list without
repeating the iRT / isotope_err_bounds checks (`check_prec_filters = false`).

Inline PSM scoring is a no-op — chromatogram points are recorded
post-deconv from `weights + id_to_col` (not from `MainUnscoredPSM`).
"""
struct FusedRTIndexed{P<:PrecEstimation} <: FusedSearchKind
    prec_estimation::P
    max_frag_rank::UInt8
end

"""
    FusedQuadEst()

Fused kind for QuadTuning. Runs three isotope-passes per precursor, each
with a one-hot `precursor_transmission` (1.0 for that iso, 0 elsewhere)
and a separate design-matrix column so the deconvolution fits an
independent weight per (precursor, iso). Fragment intensities are
normalised by the precursor's expected iso abundance (`iso_fac`) so the
fitted weights represent relative transmission.

Hardcoded: PartialPrecCapture, rank ≤ 7, fragment iso range 0:3.
Column key: `(prec_idx - 1) * 3 + iso_pass`.
No inline PSM scoring — results are read from `Hs + weights` post-deconv.
"""
struct FusedQuadEst <: FusedSearchKind end

#==========================================================
Dispatch points — all @inline; compiler folds them into each K
specialization of run_fused!. Defaults live on FusedSearchKind; each
concrete subtype overrides what it needs.
==========================================================#

# --- Outer iso-pass count (Quad = 3, others = 1) ----------------------
@inline n_iso_passes(::FusedSearchKind) = 1
@inline n_iso_passes(::FusedQuadEst)    = 3

# --- Rank filter threshold -------------------------------------------
@inline max_frag_rank_for(k::FusedStandard)  = k.max_frag_rank
@inline max_frag_rank_for(k::FusedRTIndexed) = k.max_frag_rank
@inline max_frag_rank_for(::FusedQuadEst)    = UInt8(7)

# --- Precursor estimation strategy -----------------------------------
@inline prec_estimation_for(k::FusedStandard)  = k.prec_estimation
@inline prec_estimation_for(k::FusedRTIndexed) = k.prec_estimation
@inline prec_estimation_for(::FusedQuadEst)    = PartialPrecCapture()

# --- Column ID mapping -----------------------------------------------
@inline match_column_id(::FusedSearchKind, prec_idx::UInt32, ::Int) = prec_idx
@inline match_column_id(::FusedQuadEst, prec_idx::UInt32, iso_pass::Int) =
    UInt32((prec_idx - UInt32(1)) * UInt32(3) + UInt32(iso_pass))

"""
    setup_transmission!(kind, prec_trans_buf, prec_mz, prec_charge, qfunc, iso_pass)
      -> prec_isotope_set::Tuple{Int64, Int64}

Populate `prec_trans_buf` with per-isotope quadrupole transmission for
this (precursor, iso_pass). Returns the allowable iso-set used by
`getFragIsotopes!`.

Standard: full qfunc across all isos (iso_pass ignored).
QuadEst: one-hot for `prec_iso_idx = iso_pass - 1`, set = `(iso, iso)`.
"""
@inline function setup_transmission!(::FusedSearchKind,
        prec_trans_buf::Vector{Float32}, prec_mz::Float32,
        prec_charge::UInt8, qfunc::QuadTransmissionFunction, ::Int)
    getPrecursorIsotopeTransmission!(prec_trans_buf, prec_mz, prec_charge, qfunc)
    return getPrecursorIsotopeSet(prec_mz, prec_charge, qfunc)
end

@inline function setup_transmission!(::FusedQuadEst,
        prec_trans_buf::Vector{Float32}, ::Float32,
        ::UInt8, ::QuadTransmissionFunction, iso_pass::Int)
    prec_iso_idx = iso_pass - 1
    fill!(prec_trans_buf, 0f0)
    prec_trans_buf[prec_iso_idx + 1] = 1f0
    return (prec_iso_idx, prec_iso_idx)
end

"""
    max_frag_iso_idx(kind, n_frag_isotopes, prec_isotope_set) -> Int

Largest fragment-isotope index to iterate. Standard respects
`n_frag_isotopes` and the quad-determined precursor iso upper bound;
Quad hardcodes 3 (iterates 0:3 per fragment).
"""
@inline max_frag_iso_idx(::FusedSearchKind, n_frag_isotopes::Int,
        prec_isotope_set::Tuple{<:Integer, <:Integer}) =
    min(n_frag_isotopes - 1, Int(last(prec_isotope_set)))
@inline max_frag_iso_idx(::FusedQuadEst, ::Int, ::Tuple{<:Integer, <:Integer}) = 3

"""
    frag_isotope_scale(kind, prec_sulfur, prec_mz, prec_charge, iso_pass, iso_splines)
      -> Float32

Scalar multiplier applied to each fragment isotope abundance. Standard
returns 1.0 (no scaling). Quad returns
`1 / iso_splines(min(sulfur,5), prec_iso_idx, prec_mz * prec_charge)` so
fitted weights represent relative precursor-iso transmission.
"""
@inline frag_isotope_scale(::FusedSearchKind, ::UInt8, ::Float32,
        ::UInt8, ::Int, ::IsotopeSplineModel) = 1f0

@inline function frag_isotope_scale(::FusedQuadEst, prec_sulfur::UInt8,
        prec_mz::Float32, prec_charge::UInt8, iso_pass::Int,
        iso_splines::IsotopeSplineModel)
    prec_iso_idx = iso_pass - 1
    return 1f0 / iso_splines(min(Int64(prec_sulfur), 5), prec_iso_idx,
                              prec_mz * prec_charge)
end

"""
    check_prec_filters(kind) -> Bool

Whether `run_fused!` applies the per-precursor iRT and
`isotope_err_bounds` filters. Standard: yes. Quad: no — the
scan→precursor mapping is already pre-filtered upstream, and the
classic Quad path skips these.
"""
@inline check_prec_filters(::FusedSearchKind) = true
@inline check_prec_filters(::FusedRTIndexed)  = false  # already done upstream
@inline check_prec_filters(::FusedQuadEst)    = false

const FRAGMENT_TRACE_ISOTOPE_REL_THRESHOLD = 0.25f0

@inline captures_fragment_trace_intensity(
        ::FusedSearchKind,
        ::AbstractVector{<:UnscoredPSM{Float32}}) = false
@inline captures_fragment_trace_intensity(
        ::FusedStandard,
        ::AbstractVector{MainUnscoredPSM{Float32}}) = true

@inline function fragment_trace_max_pred_intensity(
        isotopes_buf::Vector{Float32},
        frag_iso_idx_range::UnitRange{Int64},
        iso_scale::Float32)
    max_pred = 0f0
    @inbounds for iso_idx in frag_iso_idx_range
        pred_int = isotopes_buf[iso_idx + 1] * iso_scale
        if pred_int > max_pred
            max_pred = pred_int
        end
    end
    return max_pred
end

@inline function fragment_trace_isotope_contributes(
        pred_int::Float32,
        max_pred::Float32)
    return max_pred > 0f0 &&
           pred_int >= FRAGMENT_TRACE_ISOTOPE_REL_THRESHOLD * max_pred
end

"""
    record_match!(kind, unscored_psms, col, frag, iso_idx, intensity,
                   ppm_err, m_rank, prec_idx, pred_int, trace_intensity_match)

Per-match inline scoring hook, dispatched on `(kind, eltype(unscored_psms))`.

- `FusedStandard` + `MainUnscoredPSM`  → `apply_main_scoring!` (MainSearch)
- `FusedStandard` + `TuningUnscoredPSM` → `apply_tuning_scoring!` (Tuning paths)
- `FusedQuadEst` / `FusedRTIndexed`     → no-op (tuning reads `Hs + weights`
  after deconv; no PSM accumulation)
"""
@inline function ensure_unscored_capacity!(
        unscored::Vector{T}, col::Integer) where {T<:UnscoredPSM}
    needed = Int(col)
    needed <= length(unscored) && return nothing

    old_len = length(unscored)
    new_len = max(needed, old_len + 1000)
    resize!(unscored, new_len)
    default = T()
    @inbounds for i in (old_len + 1):new_len
        unscored[i] = default
    end
    return nothing
end

@inline function ensure_unscored_capacity!(
        unscored::AbstractVector{<:UnscoredPSM}, col::Integer)
    needed = Int(col)
    needed <= length(unscored) || throw(BoundsError(unscored, needed))
    return nothing
end

@inline function record_match!(::FusedStandard,
        unscored_psms::Vector{MainUnscoredPSM{Float32}}, col::Int, frag,
        iso_idx::UInt8, intensity::Float32, ppm_err::Float32,
        m_rank::Int64, prec_idx::UInt32, pred_int::Float32,
        trace_intensity_match::Bool)
    ensure_unscored_capacity!(unscored_psms, col)
    apply_main_scoring!(unscored_psms, col, frag, iso_idx,
                            intensity, ppm_err, m_rank, prec_idx, pred_int,
                            trace_intensity_match)
    return nothing
end

@inline function record_match!(::FusedStandard,
        unscored_psms::Vector{TuningUnscoredPSM{Float32}}, col::Int, frag,
        iso_idx::UInt8, intensity::Float32, ppm_err::Float32,
        m_rank::Int64, prec_idx::UInt32, ::Float32, ::Bool)
    ensure_unscored_capacity!(unscored_psms, col)
    apply_tuning_scoring!(unscored_psms, col, frag, iso_idx,
                            intensity, ppm_err, m_rank, prec_idx)
    return nothing
end

@inline record_match!(::FusedQuadEst, ::AbstractVector{<:UnscoredPSM},
        ::Int, _, ::UInt8, ::Float32, ::Float32, ::Int64, ::UInt32,
        ::Float32, ::Bool) = nothing

@inline record_match!(::FusedRTIndexed, ::AbstractVector{<:UnscoredPSM},
        ::Int, _, ::UInt8, ::Float32, ::Float32, ::Int64, ::UInt32,
        ::Float32, ::Bool) = nothing

#==========================================================
Inline scoring for the fused path — mirrors ModifyFeatures!(MainUnscoredPSM, …)
byte-for-byte, but takes a CompactFrag + iso_idx directly instead of a
FragmentMatch object (since the fused path never constructs matches).
==========================================================#

"""
    apply_main_scoring!(unscored, col, frag, iso_idx, intensity, ppm_err,
                        m_rank, prec_idx, pred_int, trace_intensity_match)

Update `unscored[col]` with one (fragment, isotope) match. Equivalent to
`ModifyFeatures!(score::MainUnscoredPSM, …)` called from classic's
`ScoreFragmentMatches!` — same branching on isotope vs mono, same rank /
topn / longest-y / longest-b / count accumulation, same `abs(ppm_err)`
error accumulation.
"""
@inline function apply_main_scoring!(unscored::Vector{MainUnscoredPSM{Float32}},
                                         col::Int,
                                         frag,
                                         iso_idx::UInt8,
                                         intensity::Float32,
                                         ppm_err::Float32,
                                         m_rank::Int64,
                                         prec_idx::UInt32,
                                         pred_int::Float32,
                                         trace_intensity_match::Bool)
    @inbounds score = unscored[col]

    best_rank            = score.best_rank
    longest_y            = score.longest_y
    longest_b            = score.longest_b
    longest_y_iso        = score.longest_y_iso
    isotope_count        = score.isotope_count
    b_count              = score.b_count
    b_int                = score.b_int
    y_count              = score.y_count
    y_int                = score.y_int
    y_count_iso          = score.y_count_iso
    p_count              = score.p_count
    non_cannonical_count = score.non_cannonical_count
    error                = score.error
    ms_file_idx          = score.ms_file_idx

    rank    = getRank(frag)
    ion_pos = getIonPosition(frag)
    frag1_int = score.frag1_int
    frag2_int = score.frag2_int
    frag3_int = score.frag3_int
    frag4_int = score.frag4_int
    frag5_int = score.frag5_int
    frag6_int = score.frag6_int
    frag7_int = score.frag7_int
    frag8_int = score.frag8_int
    top3_abs_ppm_err_sum = score.top3_abs_ppm_err_sum
    top3_ppm_err_count   = score.top3_ppm_err_count
    pred_int_sum_m0 = score.pred_int_sum_m0

    if iso_idx != UInt8(0)
        isotope_count += UInt8(1)
        if isY(frag)
            y_count_iso += UInt8(1)
            if ion_pos > longest_y_iso
                longest_y_iso = ion_pos
            end
        end
    else
        # E7: top-3 M0 fragment ppm-error capture
        if rank <= UInt8(3)
            top3_abs_ppm_err_sum += abs(ppm_err)
            top3_ppm_err_count   += UInt8(1)
        end
        # E3: matched-fragment predicted-intensity sum (M0 only, all ranks)
        pred_int_sum_m0 += pred_int
        if rank < best_rank
            best_rank = rank
        end
        if isB(frag)
            b_count += UInt8(1)
            b_int += intensity
            if ion_pos > longest_b
                longest_b = ion_pos
            end
        elseif isY(frag)
            y_count += UInt8(1)
            y_int += intensity
            if ion_pos > longest_y
                longest_y = ion_pos
            end
        elseif isP(frag)
            p_count += UInt8(1)
        else
            non_cannonical_count += UInt8(1)
        end
    end

    if trace_intensity_match
        if rank == UInt8(1);     frag1_int += intensity
        elseif rank == UInt8(2); frag2_int += intensity
        elseif rank == UInt8(3); frag3_int += intensity
        elseif rank == UInt8(4); frag4_int += intensity
        elseif rank == UInt8(5); frag5_int += intensity
        elseif rank == UInt8(6); frag6_int += intensity
        elseif rank == UInt8(7); frag7_int += intensity
        elseif rank == UInt8(8); frag8_int += intensity
        end
    end

    error += abs(ppm_err)

    @inbounds unscored[col] = MainUnscoredPSM{Float32}(
        best_rank,
        longest_y, longest_b, longest_y_iso,
        min(isotope_count, UInt8(255)),
        min(b_count, UInt8(255)), b_int,
        min(y_count, UInt8(255)), y_int,
        min(y_count_iso, UInt8(255)),
        p_count, non_cannonical_count,
        error,
        frag1_int, frag2_int, frag3_int, frag4_int, frag5_int, frag6_int,
        frag7_int, frag8_int,
        top3_abs_ppm_err_sum, top3_ppm_err_count,
        pred_int_sum_m0,
        prec_idx, ms_file_idx)
    return nothing
end

"""
    apply_tuning_scoring!(unscored::Vector{TuningUnscoredPSM}, col, frag, iso_idx,
                           intensity, ppm_err, m_rank, prec_idx)

Slim variant of `apply_main_scoring!` for tuning paths. Identical
accumulation except it skips the MainSearch-only `matched_rank_mask`
update and `frag1..8_int` captures (those fields don't exist on
`TuningUnscoredPSM`).
"""
@inline function apply_tuning_scoring!(unscored::Vector{TuningUnscoredPSM{Float32}},
                                         col::Int,
                                         frag,
                                         iso_idx::UInt8,
                                         intensity::Float32,
                                         ppm_err::Float32,
                                         m_rank::Int64,
                                         prec_idx::UInt32)
    @inbounds score = unscored[col]

    best_rank            = score.best_rank
    longest_y            = score.longest_y
    longest_b            = score.longest_b
    longest_y_iso        = score.longest_y_iso
    isotope_count        = score.isotope_count
    b_count              = score.b_count
    b_int                = score.b_int
    y_count              = score.y_count
    y_int                = score.y_int
    y_count_iso          = score.y_count_iso
    p_count              = score.p_count
    non_cannonical_count = score.non_cannonical_count
    error                = score.error
    ms_file_idx          = score.ms_file_idx

    rank    = getRank(frag)
    ion_pos = getIonPosition(frag)

    if iso_idx != UInt8(0)
        isotope_count += UInt8(1)
        if isY(frag)
            y_count_iso += UInt8(1)
            if ion_pos > longest_y_iso
                longest_y_iso = ion_pos
            end
        end
    else
        if rank < best_rank
            best_rank = rank
        end
        if isB(frag)
            b_count += UInt8(1)
            b_int += intensity
            if ion_pos > longest_b
                longest_b = ion_pos
            end
        elseif isY(frag)
            y_count += UInt8(1)
            y_int += intensity
            if ion_pos > longest_y
                longest_y = ion_pos
            end
        elseif isP(frag)
            p_count += UInt8(1)
        else
            non_cannonical_count += UInt8(1)
        end
    end

    error += abs(ppm_err)

    @inbounds unscored[col] = TuningUnscoredPSM{Float32}(
        best_rank,
        longest_y, longest_b, longest_y_iso,
        min(isotope_count, UInt8(255)),
        min(b_count, UInt8(255)), b_int,
        min(y_count, UInt8(255)), y_int,
        min(y_count_iso, UInt8(255)),
        p_count, non_cannonical_count,
        error, prec_idx, ms_file_idx)
    return nothing
end

"""
    prepare_scan_peaks!(corrected, obs_low, obs_high, mem, scan_mz, scan_int, scan_rt) -> Int

Per-scan precompute of the SoA peak table used by `run_fused!`. For each
peak calls the 4-arg `getCorrectedMzAndBounds(mem, mz, intensity, rt)` so
the intensity- and RT-dependent bias and tolerance (IntensityMassErrorModel)
are fully applied. Missing m/z or intensity → `Inf` / impossible window.

Three parallel `Vector{Float32}`s (not a tuple-packed AoS) so SIMD
bsearch primitives can load 8 adjacent values per operation.

Returns the peak count so callers can size views without rechecking.
"""
function prepare_scan_peaks!(corrected::Vector{Float32},
                              obs_low::Vector{Float32},
                              obs_high::Vector{Float32},
                              mem::AbstractMassErrorModel,
                              scan_mz::AbstractArray{Union{Missing, Float32}},
                              scan_int::AbstractArray{Union{Missing, Float32}},
                              scan_rt::Float32)
    n = length(scan_mz)
    if length(corrected) < n
        resize!(corrected, n)
        resize!(obs_low, n)
        resize!(obs_high, n)
    end
    @inbounds for i in 1:n
        m = scan_mz[i]
        it = scan_int[i]
        if ismissing(m) || ismissing(it)
            corrected[i] = Float32(Inf)
            obs_low[i]   = Float32(Inf)
            obs_high[i]  = -Float32(Inf)  # never matches any theoretical
        else
            c, lo, hi = getCorrectedMzAndBounds(mem, Float32(m), Float32(it), scan_rt)
            corrected[i] = c
            obs_low[i]   = lo
            obs_high[i]  = hi
        end
    end
    return n
end

"""
    scan_for_nearest_in_window(corrected_mz, obs_low, obs_high,
                                start_idx, n_peaks, iso_mz, conservative_high)
    -> (best_peak, best_abs, best_mz)

SIMD 8-wide forward scan over the per-peak window arrays, starting at
`start_idx`. Returns the peak that minimises `|corrected_mz[j] - iso_mz|`
among those whose tight window `[obs_low[j], obs_high[j]]` covers `iso_mz`.
Stops once `corrected_mz[j] > conservative_high` (monotonic sort order).

Three broadcast registers + three 32-byte SIMD loads per 8-peak block:
- `cmz`   = scan_corrected_mz[j:j+7]
- `lo`    = scan_obs_low[j:j+7]
- `hi`    = scan_obs_high[j:j+7]

Match mask =
  `(cmz ≤ conservative_high) & (lo ≤ iso_mz) & (hi ≥ iso_mz)`.
Set bits are iterated scalar-side to find the argmin. Scalar tail handles
partial blocks at the end of the array.

Returns `(0, typemax(Float32), 0f0)` on no match.
"""
@inline function scan_for_nearest_in_window(
    corrected_mz::Vector{Float32},
    obs_low::Vector{Float32},
    obs_high::Vector{Float32},
    start_idx::Int, n_peaks::Int,
    iso_mz::Float32, conservative_high::Float32)

    best_peak = 0
    best_abs  = typemax(Float32)
    best_mz   = Float32(0)

    iso_vec  = _fused_vbroadcast8(iso_mz)
    stop_vec = _fused_vbroadcast8(conservative_high)

    j = start_idx
    # SIMD main loop — process 8 peaks per iteration.
    @inbounds while j + 7 <= n_peaks
        cmz_vec   = _fused_vload8(corrected_mz, j)
        stop_mask = _fused_vcmple_mask(cmz_vec, stop_vec)
        stop_mask == 0x00 && return best_peak, best_abs, best_mz

        lo_vec = _fused_vload8(obs_low, j)
        hi_vec = _fused_vload8(obs_high, j)
        lo_mask = _fused_vcmple_mask(lo_vec, iso_vec)
        hi_mask = _fused_vcmpge_mask(hi_vec, iso_vec)
        match_mask = stop_mask & lo_mask & hi_mask

        # Iterate set bits; scalar comparison of |cmz - iso_mz| per match.
        m = match_mask
        while m != 0x00
            bit = trailing_zeros(m)
            peak = j + bit
            cmz_p = corrected_mz[peak]
            diff = abs(cmz_p - iso_mz)
            if diff < best_abs
                best_abs = diff
                best_peak = peak
                best_mz = cmz_p
            end
            m &= (m - UInt8(1))  # clear lowest set bit
        end

        # Any lane that failed the stop mask means we've crossed
        # conservative_high — done (corrected_mz is monotonic).
        stop_mask != 0xFF && return best_peak, best_abs, best_mz
        j += 8
    end

    # Scalar tail for the final <8 peaks.
    @inbounds while j <= n_peaks && corrected_mz[j] <= conservative_high
        if obs_low[j] <= iso_mz && iso_mz <= obs_high[j]
            cmz_j = corrected_mz[j]
            diff = abs(cmz_j - iso_mz)
            if diff < best_abs
                best_abs = diff
                best_peak = j
                best_mz = cmz_j
            end
        end
        j += 1
    end
    return best_peak, best_abs, best_mz
end

#==========================================================
Per-precursor pre-match filters.

Extracted as small, single-purpose helpers so each can be unit-tested in
isolation. All marked `@inline` because they live in the hot path and the
compiler should fold them into `run_fused!` with no call overhead.
==========================================================#

"""
    passes_irt_filter(prec_irt, scan_irt, irt_tol) -> Bool

True iff `|prec_irt - scan_irt| <= irt_tol`. Boundary inclusive.
"""
@inline passes_irt_filter(prec_irt::Float32, scan_irt::Float32, irt_tol::Float32) =
    abs(prec_irt - scan_irt) <= irt_tol

"""
    quad_window_with_iso_bounds(qfunc, prec_charge, iso_err_bounds) -> (low, high)

Returns the precursor-m/z acceptance window: the quadrupole bounds widened
by ±C13_C12_MASS_DIFF/charge × isotope-error bounds (so a precursor whose monoisotopic
m/z sits outside the literal quad window can still be admitted via an
allowed isotope offset).
"""
@inline function quad_window_with_iso_bounds(qfunc::QuadTransmissionFunction,
        prec_charge::Integer, iso_err_bounds::Tuple{<:Integer, <:Integer})
    low  = getPrecMinBound(qfunc) -
           Float32(C13_C12_MASS_DIFF * first(iso_err_bounds) / prec_charge)
    high = getPrecMaxBound(qfunc) +
           Float32(C13_C12_MASS_DIFF * last(iso_err_bounds)  / prec_charge)
    return low, high
end

"""
    passes_prec_mz_filter(prec_mz, low, high) -> Bool

True iff `low <= prec_mz <= high`. Boundary inclusive.
"""
@inline passes_prec_mz_filter(prec_mz::Float32, low::Float32, high::Float32) =
    low <= prec_mz <= high

"""
    iso_mz_for(frag_mz, iso_idx, frag_charge_inv) -> Float32

m/z of the `iso_idx`-th isotope of a fragment with monoisotopic m/z `frag_mz`
and `frag_charge_inv = 1 / frag_charge`. iso_idx=0 returns the monoisotopic
m/z. Caller hoists the reciprocal so this stays an FMA per call.
"""
@inline iso_mz_for(frag_mz::Float32, iso_idx::Integer, frag_charge_inv::Float32) =
    Float32(frag_mz + Float32(iso_idx) * C13_C12_MASS_DIFF_F32 * frag_charge_inv)

"""
    in_frag_mz_window(iso_mz, lo, hi) -> Bool

True iff `lo <= iso_mz <= hi`. Used to drop fragment isotopes whose m/z
falls outside the scan's empirical m/z range.
"""
@inline in_frag_mz_window(iso_mz::Float32, lo::Float32, hi::Float32) =
    lo <= iso_mz <= hi

"""
    conservative_half_width(mem, frag_mz) -> Float32

Half-width of the mass-error tolerance at the fragment's monoisotopic m/z.
For SimpleMEM (ppm-scaled) this is `ppm * frag_mz * 1e-6`; for IntensityMEM
it is m/z-independent. Hoisted out of the inner iso loop because iso m/z
values lie within ~2 Da of `frag_mz`, so the widest-tol over that range is
well-approximated by the value at `frag_mz` itself (<0.5% variation for
SimpleMEM at typical ppms).
"""
@inline function conservative_half_width(mem::AbstractMassErrorModel, frag_mz::Float32)
    _, hi = getMzBounds(mem, frag_mz)
    return hi - frag_mz
end

"""
    match_window(iso_mz, half_width) -> (low, high)

Symmetric m/z window `(iso_mz - half_width, iso_mz + half_width)` used for
peak-search bounds. Returned as a tuple so callers can destructure.
"""
@inline match_window(iso_mz::Float32, half_width::Float32) =
    (iso_mz - half_width, iso_mz + half_width)

"""
    compute_ppm_err(theoretical_mz, observed_mz) -> Float32

Signed ppm error between a theoretical iso m/z and the observed peak m/z.
Positive values mean the observed peak is below the theoretical mass.
Matches the inline `(iso_mz - best_mz) / (iso_mz * 1f-6)` form so the value
written into the scored PSM stays bit-identical.
"""
@inline compute_ppm_err(theoretical_mz::Float32, observed_mz::Float32) =
    (theoretical_mz - observed_mz) / (theoretical_mz * 1f-6)

"""
    push_match!(scratch, peak_row, pred_int, int_obs, iso_idx, rank)

Append a matched-peak entry into `FusedScratch`, growing the backing arrays
if necessary. Mirrors the inline form: row=UInt32(peak), nzval=pred_int,
x=int_obs, isotope=UInt8(iso_idx), rank=UInt8(rank).
"""
@inline function push_match!(scratch::FusedScratch, peak_row::Integer,
        pred_int::Float32, int_obs::Float32, iso_idx::Integer, rank::Integer)
    if scratch.n + 1 > length(scratch.row)
        grow_fused_scratch!(scratch, scratch.n + 1)
    end
    scratch.n += 1
    @inbounds begin
        scratch.row[scratch.n]     = UInt32(peak_row)
        scratch.nzval[scratch.n]   = pred_int
        scratch.x[scratch.n]       = int_obs
        scratch.isotope[scratch.n] = UInt8(iso_idx)
        scratch.rank[scratch.n]    = UInt8(rank)
    end
    return nothing
end

"""
    push_miss!(scratch, pred_int, iso_idx, rank)

Append a miss entry (no observed peak) into `FusedScratch`, growing the
backing arrays if necessary. Stores only `nzval = pred_int` and
`isotope = UInt8(iso_idx)` plus rank — misses don't carry a row index.
"""
@inline function push_miss!(scratch::FusedScratch, pred_int::Float32, iso_idx::Integer,
        rank::Integer)
    if scratch.miss_n + 1 > length(scratch.miss_nzval)
        grow_fused_scratch!(scratch, scratch.miss_n + 1)
    end
    scratch.miss_n += 1
    @inbounds begin
        scratch.miss_nzval[scratch.miss_n]   = pred_int
        scratch.miss_isotope[scratch.miss_n] = UInt8(iso_idx)
        scratch.miss_rank[scratch.miss_n]    = UInt8(rank)
    end
    return nothing
end

"""
    run_fused!(Hs, unscored_psms, id_to_col, scratch,
               corrected_peak_mz, peak_mz_len, isotopes_buf, prec_trans_buf,
               fragments, spline_data_fn, precursors_passed, prec_range,
               prec_mzs, prec_charges, prec_sulfur_counts, prec_irts,
               fragment_ranges_fn, iso_splines, quad_transmission_func, mem,
               scan_int, scan_irt, irt_tol,
               frag_mz_bounds, n_frag_isotopes, max_frag_rank,
               isotope_err_bounds, prec_estimation;
               m_rank) -> (nmatches, nmisses)

Fused per-precursor match+score+build. Replaces classic's
`selectTransitions!` + `matchPeaks!` + `buildDesignMatrix!` +
`sortSparse!` + `ScoreFragmentMatches!` in one pass, writing `Hs` in
CSC order.

**Pre-conditions**:
- Library fragments are m/z-sorted within each precursor (see
  `verify_mz_sorted`).
- `corrected_peak_mz` holds the scan's bias-corrected peak m/z values in
  positions `1:peak_mz_len` (call `corrected_peak_mz!` once before the
  first precursor on this scan).
- `Hs`, `id_to_col`, and the entries of `unscored_psms` we'll touch have
  been reset for this scan.
- `scratch` is reset between precursors inside the loop.

**Outputs**:
- `Hs` — CSC design matrix with `rowval/nzval/x/matched/isotope/colptr`
  set and `Hs.n/n_vals/m` finalized.
- `unscored_psms[col]` for each column populated via `apply_main_scoring!`.
- `id_to_col[prec_idx] = col` for precursors with ≥1 match.
- `(nmatches, nmisses)` — match ratio used by `Score!`.
"""
function run_fused!(
    kind::K,
    Hs::SparseArrayFused{UInt32, Float32},
    unscored_psms::AbstractVector{<:UnscoredPSM{Float32}},
    id_to_col::AbstractPrecursorMap{UInt16},
    scratch::FusedScratch,
    scan_corrected_mz::Vector{Float32},
    scan_obs_low::Vector{Float32},
    scan_obs_high::Vector{Float32},
    peak_mz_len::Int,
    isotopes_buf::Vector{Float32},
    prec_trans_buf::Vector{Float32},
    ion_list::LibraryFragmentLookup,
    nce_model::NceModel{Float32},
    precursors_passed::Vector{UInt32},
    prec_range::AbstractVector{Int64},
    prec_mzs::AbstractArray{Float32},
    prec_charges::AbstractArray{UInt8},
    prec_sulfur_counts::AbstractArray{UInt8},
    prec_irts::AbstractArray{Float32},
    iso_splines::IsotopeSplineModel,
    quad_transmission_func::QuadTransmissionFunction,
    mem::AbstractMassErrorModel,
    scan_int::AbstractArray{Union{Missing, Float32}},
    scan_irt::Float32,
    irt_tol::Float32,
    frag_mz_bounds::Tuple{Float32, Float32},
    n_frag_isotopes::Int64,
    isotope_err_bounds::Tuple{I, I};
    m_rank::Int64 = 3,
    scan_idx::Int64 = 0   # for blacklist lookup; 0 disables
) where {K<:FusedSearchKind, I<:Integer}

    reset!(Hs)

    fragments = getFragments(ion_list)

    n_peaks = peak_mz_len
    entry = 0
    col = UInt16(0)
    miss_row = UInt32(n_peaks)
    nmatches = 0
    nmisses = 0

    frag_mz_lo = first(frag_mz_bounds)
    frag_mz_hi = last(frag_mz_bounds)

    # Dispatched-once constants — compiler folds for each concrete K.
    n_passes      = n_iso_passes(kind)
    rank_thr      = max_frag_rank_for(kind)
    prec_est      = prec_estimation_for(kind)
    do_prec_check = check_prec_filters(kind)
    capture_trace_intensity = captures_fragment_trace_intensity(kind, unscored_psms)

    @inbounds for i in prec_range
        prec_idx = precursors_passed[i]
        prec_charge = prec_charges[prec_idx]
        prec_mz = prec_mzs[prec_idx]

        if do_prec_check
            passes_irt_filter(prec_irts[prec_idx], scan_irt, irt_tol) || continue
            mz_low, mz_high = quad_window_with_iso_bounds(
                quad_transmission_func, prec_charge, isotope_err_bounds)
            passes_prec_mz_filter(prec_mz, mz_low, mz_high) || continue
        end

        prec_sulfur = prec_sulfur_counts[prec_idx]
        spline_data = getSplineData(ion_list, nce_model, prec_charge, prec_mz)

        # Outer iso-pass loop: 1 iteration for FusedStandard (compiler
        # removes the wrapper), 3 for FusedQuadEst (one pass per isotope
        # column, fitted as a separate design-matrix column).
        for iso_pass in 1:n_passes
            # Populate quad transmission + get the allowable iso set
            # (Standard: full qfunc; Quad: one-hot for iso_pass-1).
            prec_isotope_set = setup_transmission!(kind, prec_trans_buf,
                prec_mz, prec_charge, quad_transmission_func, iso_pass)
            max_iso_idx = max_frag_iso_idx(kind, n_frag_isotopes, prec_isotope_set)
            max_iso_idx < 0 && continue
            frag_iso_idx_range = 0:max_iso_idx

            # Per-iso-pass fragment intensity scaling (Standard: 1.0,
            # Quad: 1 / iso_splines(...) so weights represent relative
            # precursor iso transmission).
            iso_scale = frag_isotope_scale(kind, prec_sulfur, prec_mz,
                prec_charge, iso_pass, iso_splines)

            reset_fused_scratch!(scratch)
            col_started = false
            this_col = UInt16(0)
            lower = 1

            frag_range = getPrecFragRange(ion_list, prec_idx)

            @inbounds for frag_idx in frag_range
                frag = fragments[frag_idx]
                rank = getRank(frag)
                rank > rank_thr && continue

                getFragIsotopes!(
                    prec_est, isotopes_buf, prec_trans_buf,
                    prec_isotope_set, frag_iso_idx_range, iso_splines,
                    prec_mz, prec_charge, prec_sulfur, frag, spline_data)
                trace_max_pred = capture_trace_intensity ?
                    fragment_trace_max_pred_intensity(isotopes_buf, frag_iso_idx_range, iso_scale) :
                    0f0

                mono_landing = lower
                frag_mz = getMz(frag)
                frag_charge_inv = 1f0 / Float32(getFragCharge(frag))

                half_width = conservative_half_width(mem, frag_mz)

                # Tighten per-iso bsearch `lo`: iso k's target is strictly greater
                # than iso k-1's, so k's first-matching-peak ≥ k-1's anchor.
                iso_anchor = lower

                @inbounds for iso_idx in frag_iso_idx_range
                    iso_mz = iso_mz_for(frag_mz, iso_idx, frag_charge_inv)
                    in_frag_mz_window(iso_mz, frag_mz_lo, frag_mz_hi) || continue

                    pred_int = isotopes_buf[iso_idx + 1] * iso_scale

                    conservative_low, conservative_high = match_window(iso_mz, half_width)

                    # Bsearch from iso_anchor (≥ lower, monotonic across isos).
                    start_idx = bsearch_hybrid(scan_corrected_mz, conservative_low, iso_anchor, n_peaks)
                    iso_anchor = start_idx  # tighten lo for the next iso

                    best_peak, best_abs, best_mz = scan_for_nearest_in_window(
                        scan_corrected_mz, scan_obs_low, scan_obs_high,
                        start_idx, n_peaks, iso_mz, conservative_high)

                    if best_peak > 0
                        nmatches += 1
                        raw_int = scan_int[best_peak]
                        int_obs = ismissing(raw_int) ? 0f0 : Float32(raw_int)

                        if !col_started
                            col += UInt16(1)
                            this_col = col
                            id_to_col[match_column_id(kind, prec_idx, iso_pass)] = this_col
                            col_started = true
                        end

                        push_match!(scratch, best_peak, pred_int, int_obs, iso_idx, rank)

                        ppm_err = compute_ppm_err(iso_mz, best_mz)
                        trace_intensity_match = capture_trace_intensity &&
                            fragment_trace_isotope_contributes(pred_int, trace_max_pred)
                        record_match!(kind, unscored_psms, Int(this_col),
                            frag, UInt8(iso_idx), int_obs, ppm_err, m_rank, prec_idx,
                            pred_int, trace_intensity_match)

                        if iso_idx == 0
                            mono_landing = best_peak
                        end
                    else
                        nmisses += 1
                        push_miss!(scratch, pred_int, iso_idx, rank)
                    end
                end  # iso_idx

                lower = mono_landing
            end  # frag_idx

            if col_started
                entry, miss_row = finalize_column!(Hs, this_col, scratch, entry, miss_row)
            end
        end  # iso_pass
    end

    grow_colptr_if_needed!(Hs, Int(col))
    @inbounds Hs.colptr[col + 1] = UInt32(entry + 1)
    Hs.n = Int(col)
    Hs.n_vals = entry
    Hs.m = Int(miss_row)

    return nmatches, nmisses
end

#==========================================================
Fused mass-error fragment collection (ParameterTuning).

Replaces the classic chain
  selectTransitions!(MassErrEstimationStrategy)
  → matchPeaks!
  → collectFragErrs(raw_masses=…)

with a single per-(precursor, fragment) walk that writes FragmentMatch
entries directly to the cumulative output buffer. No isotopes, no
scoring, no design matrix — just `(theoretical_mz, raw observed_mz,
intensity)` triples consumed downstream by `fit_mass_err_model`.

Filter (matches classic MassErrEstimationStrategy):
  - y-ions only
  - getFragCharge(frag) == 1
  - getRank(frag) ≤ max_rank   (default 6)
  - getIonPosition(frag) ≥ min_ion_position (default 5)
==========================================================#

"""
    run_fused_masserr!(fmatches, fmatches_idx,
                        scan_corrected_mz, scan_obs_low, scan_obs_high, peak_mz_len,
                        ion_list,
                        precursors_passed, prec_range,
                        mem, scan_mz, scan_int;
                        max_rank, min_ion_position) -> new_samples_idx

Walks the filtered fragment set for each precursor in `prec_range`,
matches the mono peak in `mem`'s window, writes a `MassErrSample`
to `samples[samples_idx + 1 : ...]` per match.

`observed_mz` is set to the **raw** peak m/z (not the bias-corrected
value) — the consumer (`fit_mass_err_model`) computes the bias from
`(theoretical - raw_observed)` ppm errors.

Returns the updated write head.
"""
function run_fused_masserr!(
    samples::Vector{MassErrSample},
    samples_idx::Int,
    scan_corrected_mz::Vector{Float32},
    scan_obs_low::Vector{Float32},
    scan_obs_high::Vector{Float32},
    peak_mz_len::Int,
    ion_list::LibraryFragmentLookup,
    precursors_passed::AbstractVector{UInt32},
    prec_range::UnitRange{Int64},
    mem::AbstractMassErrorModel,
    scan_mz::AbstractArray{Union{Missing, Float32}},
    scan_int::AbstractArray{Union{Missing, Float32}},
    scan_rt::Float32;
    max_rank::Int = 6,
    min_ion_position::Int = 4)

    fragments = getFragments(ion_list)
    n_peaks = peak_mz_len

    @inbounds for i in prec_range
        prec_idx = precursors_passed[i]
        frag_range = getPrecFragRange(ion_list, prec_idx)
        lower = 1  # running monotonic peak-index lower bound per precursor

        for frag_idx in frag_range
            frag = fragments[frag_idx]

            # Classic MassErr filter.
            if !(isY(frag) && getFragCharge(frag) == UInt8(1) &&
                 getRank(frag) <= UInt8(max_rank) &&
                 getIonPosition(frag) >= UInt8(min_ion_position))
                continue
            end

            frag_mz = getMz(frag)

            _, frag_conservative_upper = getMzBounds(mem, frag_mz)
            half_width = frag_conservative_upper - frag_mz
            conservative_low  = frag_mz - half_width
            conservative_high = frag_mz + half_width

            start_idx = bsearch_hybrid(scan_corrected_mz, conservative_low, lower, n_peaks)

            best_peak, _, _ = scan_for_nearest_in_window(
                scan_corrected_mz, scan_obs_low, scan_obs_high,
                start_idx, n_peaks, frag_mz, conservative_high)

            if best_peak > 0
                # Grow output buffer if needed (doubling).
                if samples_idx + 1 > length(samples)
                    extra = max(length(samples), 1024)
                    append!(samples, [MassErrSample() for _ in 1:extra])
                end
                samples_idx += 1

                raw_mz = Float32(scan_mz[best_peak])
                raw_int = scan_int[best_peak]
                int_obs = ismissing(raw_int) ? 0f0 : Float32(raw_int)

                samples[samples_idx] = MassErrSample(frag_mz, raw_mz, int_obs, scan_rt)

                # Monotonic advance: next fragment's mono mz ≥ current,
                # so its match peak index can't precede this one.
                lower = best_peak
            end
        end
    end

    return samples_idx
end
