# PSM feature computation for MainSearch.
#
# `add_psm_features!` is a single parallel pass that adds every per-PSM
# feature column consumed by the per-file LightGBM. Merged 2026-05-19 from
# two separate passes (`add_search_columns!` + `add_features!`) — same total
# work, but one parallel_foreach! task spawn / join instead of two and one
# cache-warm read of `precursor_idx` / `scan_idx` per row instead of two.

"""
    add_psm_features!(psms, search_context, spectra, ms_file_idx)

Single parallel pass that allocates and fills every per-PSM feature column
needed by the per-file LightGBM. Reads the deconv-emit columns
(:precursor_idx, :scan_idx, :error, :total_ions) and produces:

- `:decoy`, `:target`, `:cv_fold`, `:charge`, `:charge2`
- `:rt`, `:err_norm`
- `:irt_obs`, `:irt_pred`, `:irt_error`
- `:missed_cleavage`, `:num_enzymatic_termini`, `:num_variable_modifications`
- `:Mox`, `:sequence_length`, `:prec_mz`
- `:spectrum_peak_count`

Mox is computed lazily once per precursor via a length-`n_library` cache;
the count race is benign (same value written from any thread).
"""
function add_psm_features!(psms::DataFrame,
                           search_context::SearchContext,
                           spectra::MassSpecData,
                           ms_file_idx::Integer)
    # Function barrier: rt_to_irt is an abstract RtConversionModel (per-file
    # polymorphic), so calling it inside the per-PSM closure dynamically dispatches
    # and boxes the result every row. Resolve it once here, then run the body in a
    # method specialized on the concrete model type.
    rt_to_irt = getRtIrtModel(search_context, ms_file_idx)
    return _add_psm_features!(psms, search_context, spectra, ms_file_idx, rt_to_irt)
end

function _add_psm_features!(psms::DataFrame,
                            search_context::SearchContext,
                            spectra::MassSpecData,
                            ms_file_idx::Integer,
                            rt_to_irt::RtConversionModel)
    precursors_lib   = getPrecursors(getSpecLib(search_context))
    prec_is_decoy    = getIsDecoy(precursors_lib)
    prec_charge      = getCharge(precursors_lib)
    prec_mz          = getMz(precursors_lib)
    prec_irt         = getIrt(precursors_lib)
    prec_length      = getLength(precursors_lib)
    structural_mods  = getStructuralMods(precursors_lib)
    prec_missed_clv  = getMissedCleavages(precursors_lib)
    prec_enzymatic_termini = getNumEnzymaticTermini(precursors_lib)
    prec_num_var_mods = getNumVariableModifications(precursors_lib)
    scan_rts         = getRetentionTimes(spectra)
    masses           = getMzArrays(spectra)

    N = nrow(psms)

    # Inputs — type-asserted so the closure captures concrete Vector{T}.
    prec_idxs::Vector{UInt32}  = psms[!, :precursor_idx]
    scan_idxs::Vector{UInt32}  = psms[!, :scan_idx]
    error::Vector{Float32}     = psms[!, :error]
    total_ions::Vector{UInt16} = psms[!, :total_ions]

    # Allocate every output column up-front (so the closure writes typed buffers).
    decoys              = zeros(Bool,    N)
    targets             = zeros(Bool,    N)
    rt                  = zeros(Float32, N)
    err_norm            = zeros(Float16, N)
    prec_charges        = zeros(UInt8,   N)
    cv_fold             = zeros(UInt8,   N)
    irt_obs             = zeros(Float32, N)
    irt_pred            = zeros(Float32, N)
    irt_error           = zeros(Float32, N)
    missed_cleavage     = zeros(UInt8,   N)
    num_enzymatic_termini = zeros(UInt8, N)
    num_variable_modifications = zeros(UInt8, N)
    Mox                 = zeros(UInt8,   N)
    spectrum_peak_count = zeros(Float16, N)
    sequence_length     = zeros(UInt8,   N)
    prec_mzs            = zeros(Float32, N)

    # Per-precursor Mox cache (lazy; benign race).
    n_lib = length(prec_irt)
    _mox_vals     = Vector{UInt8}(undef, n_lib)
    _mox_computed = falses(n_lib)

    parallel_foreach!(N) do chunk
        @inbounds for i in chunk
            prec_idx = prec_idxs[i]
            scan_idx = scan_idxs[i]

            # Search-column work.
            is_decoy = prec_is_decoy[prec_idx]
            decoys[i]       = is_decoy
            targets[i]      = !is_decoy
            scan_rt_i       = Float32(scan_rts[scan_idx])
            rt[i]           = scan_rt_i
            prec_charges[i] = prec_charge[prec_idx]
            err_norm[i]     = Float16(min((2^min(error[i], 15f0)) / max(total_ions[i], one(UInt16)), Float32(6e4)))
            cv_fold[i]      = getCvFold(precursors_lib, prec_idx)

            # Feature work (uses scan_rt_i instead of re-reading rt[i]).
            irt_obs_i        = rt_to_irt(scan_rt_i)
            irt_pred_i       = prec_irt[prec_idx]
            irt_obs[i]       = irt_obs_i
            irt_pred[i]      = irt_pred_i
            irt_error[i]     = abs(irt_obs_i - irt_pred_i)
            missed_cleavage[i]     = prec_missed_clv[prec_idx]
            num_enzymatic_termini[i] = prec_enzymatic_termini[prec_idx]
            sequence_length[i]     = prec_length[prec_idx]
            prec_mzs[i]            = prec_mz[prec_idx]
            spectrum_peak_count[i] = length(masses[scan_idx])

            # Lazy per-precursor Mox.
            if !_mox_computed[prec_idx]
                _mox_vals[prec_idx]     = _count_mox(structural_mods[prec_idx])
                _mox_computed[prec_idx] = true
            end
            Mox[i] = _mox_vals[prec_idx]
            num_variable_modifications[i] = _num_variable_modifications_at(
                prec_num_var_mods,
                structural_mods,
                prec_idx
            )
        end
    end

    psms[!, :decoy]               = decoys
    psms[!, :target]              = targets
    psms[!, :rt]                  = rt
    psms[!, :err_norm]            = err_norm
    psms[!, :charge]              = prec_charges
    psms[!, :cv_fold]             = cv_fold
    psms[!, :charge2]             = Vector{UInt8}(prec_charges .== 2)
    psms[!, :irt_obs]             = irt_obs
    psms[!, :irt_pred]            = irt_pred
    psms[!, :irt_error]           = irt_error
    psms[!, :missed_cleavage]     = missed_cleavage
    psms[!, :num_enzymatic_termini] = num_enzymatic_termini
    psms[!, :num_variable_modifications] = num_variable_modifications
    psms[!, :Mox]                 = Mox
    psms[!, :spectrum_peak_count] = spectrum_peak_count
    psms[!, :sequence_length]     = sequence_length
    psms[!, :prec_mz]             = prec_mzs
    return nothing
end


#==========================================================
LightGBM Feature Set
==========================================================#

# Lean feature set for prescore LightGBM (fast per-file ranking)
const PRESCORE_FEATURES = [
    # ====================================================================
    # Smart-composite reduced feature set (2026-05-13).
    # Built from family-aware ablation campaign on 2-file Olsen Exploris +
    # 2-file MTAC 3P with cross-dataset reconciliation. Dropped 40 features
    # whose drop_all or keep-only-winner was positive on BOTH datasets.
    # Result on full pipeline (with MBR):
    #   Olsen: 91,193 → 91,576 final IDs (+383), 12,458 → 12,555 PGs (+97)
    #   MTAC:  180,428 → 180,672 IDs (+244),  20,107 → 20,335 PGs (+228)
    # ====================================================================

    # Core PSM / sequence metrics
    :fitted_manhattan_distance, :irt_error, :poisson, :err_norm,
    :total_ions, :missed_cleavage, :num_enzymatic_termini,
    :y_count, :weight, :gof,
    :Mox, :spectrum_peak_count, :sequence_length,
    :fitted_hellinger,

    # Cross-precursor at-scan
    :weight_ratio_at_scan, :weight_rank_at_scan,

    # MS1 features (chromatogram-correlations dropped, only m0_corr kept)
    :ms1_m0_mass_err_ppm,
    :ms1_weight_apex_to_m0_apex_irt,
    :ms1_m0_intensity, :ms1_m1_intensity,
    :ms1_m1_to_m0_ratio, :ms1_m1_to_m0_pred,
    :ms1_isotope_dotp_m0_m1_m2,
    :ms1_m0_m1_m2_window_fraction, :ms1_ms2_explained_delta,
    :ms1_m0_m1_m2_window_fraction_pc, :ms1_ms2_explained_delta_pc,

    # Per-rank top-8 fragment trace intensities; each rank sums matched isotope
    # peaks predicted at >=25% of the fragment's most abundant isotope.
    :frag1_int, :frag2_int, :frag3_int, :frag4_int,
    :frag5_int, :frag6_int, :frag7_int, :frag8_int,

    # Fragment-chromatogram correlations (frag_w_corr family dropped; pairs dropped)
    # frag_corr_mean_pairwise (Spearman) dropped 2026-05-13 — cross-dataset test
    # showed Olsen +471 IDs / MTAC +649 IDs. Saves 15 rank-sorts per precursor.
    :frag_apex_dispersion_irt,
    :n_correlated_fragments,
    :n_correlated_fragments_bitvec_rank,
    :frag_corr_strength,
    :frag_corr_effective_n,
    :frag_corr_best_m0,

    # Batch E features (E7, E14, E6 M0 kept; E1/E2 pred_obs dropped via composite)
    :top3_ms2_mass_error_mean,
    :delta_frame_peak_center,
    :log_by_ratio_m0,

    # n_scans (cross-experiment PSM count)
    :n_scans,

    # Backported from ADVANCED (verified +1k/+3k final IDs on Olsen/MTAC).
    # Tier-2 drop-all-5 (2026-05-13) removed: rt_fwhm, num_scans, irt_pred,
    # best_rank_iso, total_ions_iso — individually each lost 176–1,143 IDs
    # but dropped together they're net-neutral (Olsen +138, MTAC −405).
    :prec_mz,
    :irt_fwhm,
    :smoothness,
    :log2_intensity_explained, :longest_y,
    # 2026-05-21: min_* features removed (added 2026-05-20, ~+1% IDs not
    # worth the compute cost).

    # Neighborhood (windowed PSM-quality) family DROPPED 2026-05-18:
    # best_max_residual_3scan, best_gof_5scan, best_manhattan_5scan,
    # best_max_residual_5scan, irt_dist_best_gof_5scan,
    # irt_dist_best_manhattan_5scan, irt_dist_best_max_residual_5scan,
    # worst_max_residual_11scan, worst_manhattan_11scan. Overnight sweep
    # on Astral/Exploris/MTAC (6 files each) showed all 9 features net-
    # harmful: Astral +0.86%, Exploris +1.26%, MTAC +2.47% precursor IDs
    # when removed. add_neighborhood_features! deleted in the same commit.
]

"""
    add_ms1_features!(psms, spectra, search_context, ms_file_idx;
                       ms1_ppm_tol=10.0f0)

Per-PSM MS1 point-lookup features. For each PSM at (precursor_idx, MS2 scan_idx):

1. Find the nearest MS1 scan in time (binary search on RT).
2. Search the MS1 spectrum for M0, M+1, and M+2 of the precursor (charge-aware
   m/z = prec_mz + iso * C13_C12_MASS_DIFF / charge), within +/-ms1_ppm_tol ppm.
3. Populate MS1 point-lookup features per PSM:
   - `ms1_m0_mass_err_ppm` (Float32): |observed − theoretical| in ppm; 0 if no M0
   - `ms1_m0_intensity` (Float32): observed M0 intensity; 0 if not matched
   - `ms1_m1_intensity` (Float32): observed M+1 intensity; 0 if not matched
   - `ms1_isotope_dotp_m0_m1_m2` (Float32): cosine dot product between
     observed M0/M1/M2 intensities and isotope-spline expected M0/M1/M2
   - `ms1_m0_m1_m2_window_fraction` (Float32): M0/M1/M2 intensity fraction
     of the paired MS2 isolation window in the nearest MS1 scan
   - `ms1_ms2_explained_delta` (Float32): MS1 window fraction minus the MS2
     explained-intensity fraction
   - `ms1_m0_m1_m2_window_fraction_pc` and `ms1_ms2_explained_delta_pc`:
     pseudocounted variants using the nearest-MS1 scan noise floor

The intensities are also consumed by `_add_ms1_chromatogram_features!` to build
per-precursor M0 / M+1 / weight chromatograms and compute correlation features.

Earlier versions of this function also computed `ms1_iso_count`, `ms1_m{0,2,3}_matched`,
`ms1_log_iso_obs_pred[_m{2,3}]` (sulfur-aware predicted/observed isotope ratios), and
M+2/M+3 lookups. Feature-importance analysis (2026-05-10) showed the LGBM model recovers
all signal from raw `m0/m1_intensity` (presence ⇔ intensity > 0; ratio implicitly via
tree splits); the explicit features were redundant. Removing them costs ~174 IDs at
q≤.01 within run-to-run noise and slightly improves q≤.001 (+296).

If no MS1 scan is available (file has only MS2), all features are zeroed.
"""
# Top-level helpers for the MS1 per-PSM peak search. Operate on clean
# Vector{Float32} arrays (missings stripped during cache refresh — see
# _ms1_refresh_cache!). Use the same `bsearch_hybrid` primitive as the
# MS2 fused-match path so we get a branchless binary search to the window
# lower bound, followed by a tight scalar walk for the in-window peaks. M0
# lookups can reuse scan-run min/max bounds; M1/M2 reuse the previous isotope
# lower-bound cursor instead of paying another binary search.
@inline function _ms1_find_peak_at_lower_bound(mz::Vector{Float32}, intens::Vector{Float32},
                                               mz_target::Float32, hi::Float32,
                                               i0::Int, last_idx::Int)
    best_diff = Inf32
    best_int  = 0f0
    best_mz   = 0f0
    best_idx  = 0
    found = false
    last = min(last_idx, length(mz))
    i0 = max(i0, 1)
    i0 > last && return (false, 0f0, 0f0, 0, i0)
    @inbounds @fastmath for j in i0:last
        m = mz[j]
        m > hi && break
        diff = abs(m - mz_target)
        if diff < best_diff
            best_diff = diff
            best_int  = intens[j]
            best_mz   = m
            best_idx  = j
            found = true
        end
    end
    return found, best_int, best_mz, best_idx, i0
end

@inline function _ms1_find_peak_bsearch(mz::Vector{Float32}, intens::Vector{Float32},
                                        mz_target::Float32, ms1_ppm_tol::Float32)
    half_tol = mz_target * ms1_ppm_tol * 1f-6
    lo = mz_target - half_tol
    hi = mz_target + half_tol
    n_peaks = length(mz)
    n_peaks == 0 && return (false, 0f0, 0f0, 0, 1)
    i0 = bsearch_hybrid(mz, lo, 1, n_peaks)
    return _ms1_find_peak_at_lower_bound(mz, intens, mz_target, hi, i0, n_peaks)
end

@inline function _ms1_find_peak_from_cursor(mz::Vector{Float32}, intens::Vector{Float32},
                                            mz_target::Float32, ms1_ppm_tol::Float32,
                                            start_cursor::Int)
    return _ms1_find_peak_from_cursor(mz, intens, mz_target, ms1_ppm_tol,
                                      start_cursor, length(mz))
end

@inline function _ms1_find_peak_from_cursor(mz::Vector{Float32}, intens::Vector{Float32},
                                            mz_target::Float32, ms1_ppm_tol::Float32,
                                            start_cursor::Int, last_idx::Int)
    half_tol = mz_target * ms1_ppm_tol * 1f-6
    lo = mz_target - half_tol
    hi = mz_target + half_tol
    n_peaks = length(mz)
    n_peaks == 0 && return (false, 0f0, 0f0, 0, 1)
    last = min(last_idx, n_peaks)
    last < 1 && return (false, 0f0, 0f0, 0, 1)

    i0 = if start_cursor < 1
        1
    elseif start_cursor > last
        last + 1
    else
        start_cursor
    end

    # In the hot path, isotope targets are monotone increasing, so this is a
    # short forward walk. The fallback keeps the helper exact for accidental
    # non-monotone use.
    if i0 > 1 && i0 <= n_peaks
        @inbounds if mz[i0 - 1] >= lo
            i0 = bsearch_hybrid(mz, lo, 1, last)
        end
    elseif i0 == last + 1
        @inbounds if mz[last] >= lo
            i0 = bsearch_hybrid(mz, lo, 1, last)
        end
    end

    @inbounds while i0 <= last && mz[i0] < lo
        i0 += 1
    end
    return _ms1_find_peak_at_lower_bound(mz, intens, mz_target, hi, i0, last)
end

@inline function _ms1_find_peak_from_bounds(mz::Vector{Float32}, intens::Vector{Float32},
                                            mz_target::Float32, ms1_ppm_tol::Float32,
                                            floor_cursor::Int, ceiling_cursor::Int)
    half_tol = mz_target * ms1_ppm_tol * 1f-6
    lo = mz_target - half_tol
    hi = mz_target + half_tol
    n_peaks = length(mz)
    n_peaks == 0 && return (false, 0f0, 0f0, 0, 1)
    first = clamp(floor_cursor, 1, n_peaks)
    last = min(ceiling_cursor, n_peaks)
    first > last && return (false, 0f0, 0f0, 0, first)
    if first > 1
        @inbounds if mz[first - 1] >= lo
            first = bsearch_hybrid(mz, lo, 1, last)
        end
    end
    i0 = first <= last ? bsearch_hybrid(mz, lo, first, last) : first
    return _ms1_find_peak_at_lower_bound(mz, intens, mz_target, hi, i0, last)
end

@inline function _ms1_find_peak(mz::Vector{Float32}, intens::Vector{Float32},
                                mz_target::Float32, ms1_ppm_tol::Float32)
    found, best_int, best_mz, best_idx, _ =
        _ms1_find_peak_bsearch(mz, intens, mz_target, ms1_ppm_tol)
    return found, best_int, best_mz, best_idx
end

@inline function _ms1_sum_window_intensity(mz::Vector{Float32}, intens::Vector{Float32},
                                           low_mz::Float32, high_mz::Float32)
    (isfinite(low_mz) && isfinite(high_mz) && high_mz >= low_mz) || return 0f0
    n_peaks = length(mz)
    n_peaks == 0 && return 0f0
    i0 = bsearch_hybrid(mz, low_mz, 1, n_peaks)
    total = 0f0
    @inbounds @fastmath for j in i0:n_peaks
        m = mz[j]
        m > high_mz && break
        total += intens[j]
    end
    return total
end

@inline function _ms1_m0_m1_m2_window_fraction(m0_int, m1_int, m2_int,
                                               window_intensity)::Float32
    denom = Float32(window_intensity)
    denom > 0f0 || return 0f0
    numer = max(Float32(m0_int), 0f0) + max(Float32(m1_int), 0f0) +
            max(Float32(m2_int), 0f0)
    return clamp(numer / denom, 0f0, 1f0)
end

@inline function _ms1_noise_floor(intens::Vector{Float32})::Float32
    floor = Inf32
    @inbounds @fastmath for it in intens
        if it > 0f0 && it < floor
            floor = it
        end
    end
    return isfinite(floor) ? floor : 0f0
end

@inline function _ms1_m0_m1_m2_window_fraction_pseudocount(m0_int, m1_int, m2_int,
                                                           m0_hit::Bool, m1_hit::Bool, m2_hit::Bool,
                                                           window_intensity,
                                                           noise_floor)::Float32
    pc = max(Float32(noise_floor), 0f0)
    m0_eff = m0_hit ? max(Float32(m0_int), 0f0) : pc
    m1_eff = m1_hit ? max(Float32(m1_int), 0f0) : pc
    m2_eff = m2_hit ? max(Float32(m2_int), 0f0) : pc
    denom = max(Float32(window_intensity), 0f0) + 3f0 * pc
    denom > 0f0 || return 0f0
    return clamp((m0_eff + m1_eff + m2_eff) / denom, 0f0, 1f0)
end

@inline function _ms1_isotope_dotp_m0_m1_m2(m0_int, m1_int, m2_int,
                                            pred0, pred1, pred2)::Float32
    obs0 = max(Float32(m0_int), 0f0)
    obs1 = max(Float32(m1_int), 0f0)
    obs2 = max(Float32(m2_int), 0f0)
    exp0 = max(Float32(pred0), 0f0)
    exp1 = max(Float32(pred1), 0f0)
    exp2 = max(Float32(pred2), 0f0)
    obs_norm2 = obs0 * obs0 + obs1 * obs1 + obs2 * obs2
    exp_norm2 = exp0 * exp0 + exp1 * exp1 + exp2 * exp2
    denom = sqrt(obs_norm2 * exp_norm2)
    denom > 0f0 || return 0f0
    dotp = (obs0 * exp0 + obs1 * exp1 + obs2 * exp2) / denom
    return clamp(Float32(dotp), 0f0, 1f0)
end

@inline function _ms2_explained_fraction(log2_intensity_explained)::Float32
    x = Float32(log2_intensity_explained)
    isfinite(x) || return 0f0
    return clamp(exp2(x), 0f0, 1f0)
end

@inline function _ms1_ms2_explained_delta(ms1_fraction, log2_intensity_explained)::Float32
    return Float32(ms1_fraction) - _ms2_explained_fraction(log2_intensity_explained)
end

function _initialize_ms1_isolation_window_features!(psms)
    n = nrow(psms)
    psms[!, :ms1_m0_m1_m2_window_fraction] = zeros(Float32, n)
    psms[!, :ms1_ms2_explained_delta] = zeros(Float32, n)
    psms[!, :ms1_m0_m1_m2_window_fraction_pc] = zeros(Float32, n)
    psms[!, :ms1_ms2_explained_delta_pc] = zeros(Float32, n)
    psms[!, :ms1_isotope_dotp_m0_m1_m2] = zeros(Float32, n)
    return nothing
end

struct _M0PeakCompetitionScratch
    totals::Dict{UInt64,Float32}
    peak_counts::Dict{UInt64,UInt16}
    seen_peak_precursors::Set{UInt128}
    same_mz_counts::Dict{UInt32,UInt16}
    seen_mz_precursors::Set{UInt64}
end

function _M0PeakCompetitionScratch()
    return _M0PeakCompetitionScratch(
        Dict{UInt64,Float32}(),
        Dict{UInt64,UInt16}(),
        Set{UInt128}(),
        Dict{UInt32,UInt16}(),
        Set{UInt64}(),
    )
end

@inline _inc_u16_saturating(x::UInt16) = x == typemax(UInt16) ? x : x + UInt16(1)

@inline function _fragment_sum6(frag_cols, row::Int)::Float32
    return Float32(frag_cols[1][row]) + Float32(frag_cols[2][row]) +
           Float32(frag_cols[3][row]) + Float32(frag_cols[4][row]) +
           Float32(frag_cols[5][row]) + Float32(frag_cols[6][row])
end

function _add_m0_peak_fragment_competition_run!(
        frac_out::Vector{Float32},
        peak_count_out::Vector{UInt16},
        same_mz_count_out::Vector{UInt16},
        precursor_idx::Vector{UInt32},
        prec_mzs,
        m0_peak_keys::AbstractVector{UInt64},
        frag_cols,
        rows::UnitRange{Int},
        scratch::_M0PeakCompetitionScratch)
    isempty(rows) && return nothing

    totals = scratch.totals
    peak_counts = scratch.peak_counts
    seen_peak_precursors = scratch.seen_peak_precursors
    same_mz_counts = scratch.same_mz_counts
    seen_mz_precursors = scratch.seen_mz_precursors
    empty!(totals)
    empty!(peak_counts)
    empty!(seen_peak_precursors)
    empty!(same_mz_counts)
    empty!(seen_mz_precursors)

    first_row = first(rows)
    have_prec_mzs = prec_mzs !== nothing
    @inbounds for row in rows
        pid = precursor_idx[row]
        key = m0_peak_keys[row - first_row + 1]
        if key != 0
            seen_key = (UInt128(key) << 32) | UInt128(pid)
            if !(seen_key in seen_peak_precursors)
                push!(seen_peak_precursors, seen_key)
                peak_counts[key] = _inc_u16_saturating(get(peak_counts, key, UInt16(0)))
            end
            frag_sum = _fragment_sum6(frag_cols, row)
            if frag_sum > 0f0
                totals[key] = get(totals, key, 0f0) + frag_sum
            end
        end
        if have_prec_mzs && 1 <= Int(pid) <= length(prec_mzs)
            mz_key = reinterpret(UInt32, Float32(prec_mzs[pid]))
            seen_key = (UInt64(mz_key) << 32) | UInt64(pid)
            if !(seen_key in seen_mz_precursors)
                push!(seen_mz_precursors, seen_key)
                same_mz_counts[mz_key] = _inc_u16_saturating(get(same_mz_counts, mz_key, UInt16(0)))
            end
        end
    end

    @inbounds for row in rows
        pid = precursor_idx[row]
        key = m0_peak_keys[row - first_row + 1]
        if key != 0
            peak_count_out[row] = get(peak_counts, key, UInt16(0))
            denom = get(totals, key, 0f0)
            if denom > 0f0
                frag_sum = _fragment_sum6(frag_cols, row)
                frac_out[row] = frag_sum > 0f0 ? frag_sum / denom : 0f0
            end
        end
        if have_prec_mzs && 1 <= Int(pid) <= length(prec_mzs)
            mz_key = reinterpret(UInt32, Float32(prec_mzs[pid]))
            same_mz_count_out[row] = get(same_mz_counts, mz_key, UInt16(0))
        end
    end
    return nothing
end

# Filter the (possibly Union{Missing,Float32}) MS1 spectrum into per-task
# scratch Vector{Float32} buffers, dropping any peaks where either m/z or
# intensity is missing. Result: dense, type-stable arrays that bsearch_hybrid
# and the SIMD-friendly walk can consume directly.
@inline function _ms1_refresh_cache!(scratch_mz::Vector{Float32}, scratch_int::Vector{Float32},
                                      raw_mz, raw_int)
    n_raw = length(raw_mz)
    n_clean = 0
    @inbounds for j in 1:n_raw
        if !ismissing(raw_mz[j]) && !ismissing(raw_int[j])
            n_clean += 1
        end
    end
    resize!(scratch_mz, n_clean)
    resize!(scratch_int, n_clean)
    k = 0
    @inbounds for j in 1:n_raw
        m = raw_mz[j]; ismissing(m) && continue
        it = raw_int[j]; ismissing(it) && continue
        k += 1
        scratch_mz[k]  = Float32(m)
        scratch_int[k] = Float32(it)
    end
    return scratch_mz, scratch_int
end

function _scan_run_starts(scan::Vector{UInt32})
    n = length(scan)
    starts = Int[1]
    n == 0 && return starts
    sizehint!(starts, max(1, n ÷ 4))
    @inbounds for i in 2:n
        if scan[i] != scan[i-1]
            push!(starts, i)
        end
    end
    push!(starts, n + 1)
    return starts
end

function _ms1_m0_peak_competition_inputs(psms, prec_mzs)
    required = (:precursor_idx, :frag1_int, :frag2_int, :frag3_int,
                :frag4_int, :frag5_int, :frag6_int, :frag7_int, :frag8_int)
    all(c -> hasproperty(psms, c), required) || return nothing
    return (
        psms[!, :ms1_m0_peak_frag_intensity_fraction]::Vector{Float32},
        psms[!, :ms1_m0_peak_n_precursors]::Vector{UInt16},
        psms[!, :scan_prec_mz_n_precursors]::Vector{UInt16},
        psms[!, :precursor_idx]::Vector{UInt32},
        prec_mzs,
        (psms[!, :frag1_int], psms[!, :frag2_int], psms[!, :frag3_int],
         psms[!, :frag4_int], psms[!, :frag5_int], psms[!, :frag6_int],
         psms[!, :frag7_int], psms[!, :frag8_int]),
    )
end

# Per-scan-run worker. The input PSM table is contiguous by :scan_idx at this
# point in MainSearch, so each run shares one nearest MS1 scan and one MS2
# isolation window. That lets us compute window intensity/noise once per scan
# run instead of once per candidate PSM.
function _ms1_lookup_scan_runs!(psms, spectra,
                                prec_mzs, prec_charges, prec_sulfurs, iso_splines,
                                scan_to_ms1::Vector{Int32},
                                scan_window_low::Vector{Float32},
                                scan_window_high::Vector{Float32},
                                ms1_ppm_tol::Float32, ms1_ppm_offset::Float32,
                                isotope_spacing::Float32,
                                starts::Vector{Int}, run_chunk,
                                competition_inputs)
    scan_idxs::Vector{UInt32} = psms[!, :scan_idx]
    precursor_idxs::Vector{UInt32} = psms[!, :precursor_idx]
    log2_intensity_explained::Vector{Float16} = psms[!, :log2_intensity_explained]

    ms1_m0_mass_err_ppm::Vector{Float32} = psms[!, :ms1_m0_mass_err_ppm]
    ms1_m0_intensity::Vector{Float32} = psms[!, :ms1_m0_intensity]
    ms1_m1_intensity::Vector{Float32} = psms[!, :ms1_m1_intensity]
    ms1_m1_to_m0_ratio::Vector{Float32} = psms[!, :ms1_m1_to_m0_ratio]
    ms1_m1_to_m0_pred::Vector{Float32} = psms[!, :ms1_m1_to_m0_pred]
    ms1_isotope_dotp::Vector{Float32} = psms[!, :ms1_isotope_dotp_m0_m1_m2]
    ms1_window_fraction::Vector{Float32} = psms[!, :ms1_m0_m1_m2_window_fraction]
    ms1_delta::Vector{Float32} = psms[!, :ms1_ms2_explained_delta]
    ms1_window_fraction_pc::Vector{Float32} = psms[!, :ms1_m0_m1_m2_window_fraction_pc]
    ms1_delta_pc::Vector{Float32} = psms[!, :ms1_ms2_explained_delta_pc]

    cached_mz  = Vector{Float32}()
    cached_int = Vector{Float32}()
    cached_ms1_idx::Int = -1   # force first miss
    cached_noise_floor = 0f0
    m0_peak_keys = UInt64[]
    competition_scratch = _M0PeakCompetitionScratch()

    @inbounds for r in run_chunk
        run_start = starts[r]
        run_end = starts[r+1] - 1
        scan_id = Int(scan_idxs[run_start])
        ms1_idx = Int(scan_to_ms1[scan_id])

        if ms1_idx != cached_ms1_idx
            cached_ms1_idx = ms1_idx
            _ms1_refresh_cache!(cached_mz, cached_int,
                                getMzArray(spectra, ms1_idx),
                                getIntensityArray(spectra, ms1_idx))
            cached_noise_floor = _ms1_noise_floor(cached_int)
        end

        window_low = scan_window_low[scan_id]
        window_high = scan_window_high[scan_id]
        window_intensity = _ms1_sum_window_intensity(cached_mz, cached_int, window_low, window_high)
        resize!(m0_peak_keys, run_end - run_start + 1)
        fill!(m0_peak_keys, UInt64(0))

        m0_min_lo = Inf32
        m0_max_hi = -Inf32
        ms1_bias_factor = 1f0 + ms1_ppm_offset * 1f-6
        @inbounds for i in run_start:run_end
            pid = UInt32(precursor_idxs[i])
            target_m0 = Float32(prec_mzs[pid]) * ms1_bias_factor
            half_tol = target_m0 * ms1_ppm_tol * 1f-6
            m0_min_lo = min(m0_min_lo, target_m0 - half_tol)
            m0_max_hi = max(m0_max_hi, target_m0 + half_tol)
        end
        n_cached_peaks = length(cached_mz)
        m0_floor_cursor = 1
        m0_ceiling_cursor = 0
        if n_cached_peaks > 0 && isfinite(m0_min_lo) && isfinite(m0_max_hi)
            m0_floor_cursor = bsearch_hybrid(cached_mz, m0_min_lo, 1, n_cached_peaks)
            if m0_floor_cursor <= n_cached_peaks
                first_after_m0 = bsearch_hybrid(
                    cached_mz,
                    nextfloat(m0_max_hi),
                    m0_floor_cursor,
                    n_cached_peaks,
                )
                m0_ceiling_cursor = min(first_after_m0 - 1, n_cached_peaks)
            end
        end
        prev_m0_lo = -Inf32
        prev_m0_cursor = m0_floor_cursor

        for i in run_start:run_end
            pid = UInt32(precursor_idxs[i])
            prec_mz = Float32(prec_mzs[pid])
            prec_chg = Int(prec_charges[pid])
            prec_chg == 0 && (prec_chg = 1)
            inv_charge = 1f0 / Float32(prec_chg)

            # Shift theoretical target by the fitted per-file ppm bias so peaks
            # are searched where they actually land. Residual feature
            # ms1_m0_mass_err_ppm is computed against the shifted target.
            target_m0 = prec_mz * ms1_bias_factor
            target_m1 = target_m0 + isotope_spacing * inv_charge
            target_m2 = target_m0 + 2f0 * isotope_spacing * inv_charge

            target_m0_half_tol = target_m0 * ms1_ppm_tol * 1f-6
            target_m0_lo = target_m0 - target_m0_half_tol
            m0_hit, m0_int, m0_mz, m0_peak_idx, m0_cursor = if target_m0_lo >= prev_m0_lo
                _ms1_find_peak_from_cursor(
                    cached_mz, cached_int, target_m0, ms1_ppm_tol,
                    prev_m0_cursor, m0_ceiling_cursor,
                )
            else
                _ms1_find_peak_from_bounds(
                    cached_mz, cached_int, target_m0, ms1_ppm_tol,
                    m0_floor_cursor, m0_ceiling_cursor,
                )
            end
            prev_m0_lo = target_m0_lo
            prev_m0_cursor = m0_cursor
            m1_hit, m1_int, _, _, m1_cursor =
                _ms1_find_peak_from_cursor(cached_mz, cached_int, target_m1, ms1_ppm_tol, m0_cursor)
            m2_hit, m2_int, _, _, _ =
                _ms1_find_peak_from_cursor(cached_mz, cached_int, target_m2, ms1_ppm_tol, m1_cursor)

            sulf = clamp(Int(prec_sulfurs[pid]), 0, 5)
            neutral_mass = (prec_mz - Float32(PROTON)) * Float32(prec_chg)
            p0 = max(iso_splines(Int64(sulf), Int64(0), neutral_mass), 0f0)
            p1 = max(iso_splines(Int64(sulf), Int64(1), neutral_mass), 0f0)
            p2 = max(iso_splines(Int64(sulf), Int64(2), neutral_mass), 0f0)

            obs_m0 = m0_hit ? m0_int : 0f0
            obs_m1 = m1_hit ? m1_int : 0f0
            obs_m2 = m2_hit ? m2_int : 0f0
            ms1_m0_intensity[i] = obs_m0
            ms1_m1_intensity[i] = obs_m1
            ms1_isotope_dotp[i] = _ms1_isotope_dotp_m0_m1_m2(obs_m0, obs_m1, obs_m2, p0, p1, p2)

            m0_in_window = m0_hit && target_m0 >= window_low && target_m0 <= window_high
            m1_in_window = m1_hit && target_m1 >= window_low && target_m1 <= window_high
            m2_in_window = m2_hit && target_m2 >= window_low && target_m2 <= window_high
            ms1_fraction = _ms1_m0_m1_m2_window_fraction(
                m0_in_window ? m0_int : 0f0,
                m1_in_window ? m1_int : 0f0,
                m2_in_window ? m2_int : 0f0,
                window_intensity,
            )
            ms1_window_fraction[i] = ms1_fraction
            ms1_delta[i] = _ms1_ms2_explained_delta(ms1_fraction, log2_intensity_explained[i])

            ms1_fraction_pc = _ms1_m0_m1_m2_window_fraction_pseudocount(
                m0_int, m1_int, m2_int,
                m0_in_window, m1_in_window, m2_in_window,
                window_intensity,
                cached_noise_floor,
            )
            ms1_window_fraction_pc[i] = ms1_fraction_pc
            ms1_delta_pc[i] = _ms1_ms2_explained_delta(ms1_fraction_pc, log2_intensity_explained[i])

            if m0_hit
                m0_peak_keys[i - run_start + 1] =
                    (UInt64(UInt32(ms1_idx)) << 32) | UInt64(UInt32(m0_peak_idx))
                ms1_m0_mass_err_ppm[i] = abs(m0_mz - target_m0) / target_m0 * 1f6
            end
            if m0_hit && m0_int > 0f0 && m1_hit
                obs_ratio = m1_int / m0_int
                prec_mass_approx = prec_mz * Float32(prec_chg)
                ratio_p0 = max(iso_splines(Int64(sulf), Int64(0), prec_mass_approx), 1f-12)
                ratio_p1 = max(iso_splines(Int64(sulf), Int64(1), prec_mass_approx), 1f-12)
                ms1_m1_to_m0_ratio[i] = Float32(obs_ratio)
                ms1_m1_to_m0_pred[i] = Float32(ratio_p1 / ratio_p0)
            end
        end

        if competition_inputs !== nothing
            frac_out, peak_count_out, same_mz_count_out,
                comp_precursor_idx, comp_prec_mzs, frag_cols = competition_inputs
            _add_m0_peak_fragment_competition_run!(
                frac_out,
                peak_count_out,
                same_mz_count_out,
                comp_precursor_idx,
                comp_prec_mzs,
                m0_peak_keys,
                frag_cols,
                run_start:run_end,
                competition_scratch,
            )
        end
    end
    return nothing
end

"""
    add_ms1_lookup_features!(psms, spectra, search_context, ms_file_idx)

Per-PSM MS1 spectrum lookup. Populates:
- `:ms1_m0_intensity`
- `:ms1_m1_intensity`
- `:ms1_m0_mass_err_ppm`
- `:ms1_m1_to_m0_ratio`
- `:ms1_m1_to_m0_pred`
- `:ms1_isotope_dotp_m0_m1_m2`
- `:ms1_m0_m1_m2_window_fraction`
- `:ms1_ms2_explained_delta`
- `:ms1_m0_m1_m2_window_fraction_pc`
- `:ms1_ms2_explained_delta_pc`
- `:ms1_m0_peak_frag_intensity_fraction`
- `:ms1_m0_peak_n_precursors`
- `:scan_prec_mz_n_precursors`

For each PSM: finds the nearest MS1 scan by RT, extracts M0/M1 intensities
within a ppm tolerance window around the predicted precursor m/z, and
computes the observed/predicted M1:M0 ratio plus M0/M1/M2 isolation-window
features. Scan-local M0 peak competition is computed per contiguous MS2 scan
run using thread-local scratch.

**PRECONDITION (perf-only)**: the per-task MS1-spectrum cache hits ~once
per unique `scan_idx` in each run chunk, so this is dramatically faster when
the input PSMs are contiguous-by-:scan_idx — the natural ordering of the
deconv output. Call BEFORE `permute_psms_by_precursor_idx!`.
"""
function add_ms1_lookup_features!(psms::DataFrame,
                            spectra,
                            search_context,
                            ms_file_idx::Integer)
    # Per-file MS1 tolerance + bias fitted in ParameterTuningSearch
    # (collect_ms1_residuals → setMs1MassErrorModel!). Falls back to the
    # ±30 ppm default in getMs1MassErrorModel if no fit was installed.
    mem = getMs1MassErrorModel(search_context, ms_file_idx)
    ms1_ppm_tol    = max(getLeftTol(mem), getRightTol(mem))
    ms1_ppm_offset = getMassOffset(mem)
    n = nrow(psms)
    psms[!, :ms1_m0_mass_err_ppm]   = zeros(Float32, n)
    psms[!, :ms1_m0_intensity]      = zeros(Float32, n)
    psms[!, :ms1_m1_intensity]      = zeros(Float32, n)
    psms[!, :ms1_m1_to_m0_ratio]     = zeros(Float32, n)
    psms[!, :ms1_m1_to_m0_pred]      = zeros(Float32, n)
    _initialize_ms1_isolation_window_features!(psms)
    psms[!, :ms1_m0_peak_frag_intensity_fraction] = zeros(Float32, n)
    psms[!, :ms1_m0_peak_n_precursors] = zeros(UInt16, n)
    psms[!, :scan_prec_mz_n_precursors] = zeros(UInt16, n)
    n == 0 && return

    # 1. Build MS1 scan index (sorted by RT) for fast nearest-MS1 lookup
    n_scans = length(spectra)
    ms1_scan_idxs = Int[]
    ms1_scan_rts  = Float32[]
    for s in 1:n_scans
        if getMsOrder(spectra, s) == 1
            push!(ms1_scan_idxs, s)
            push!(ms1_scan_rts, Float32(getRetentionTime(spectra, s)))
        end
    end
    if isempty(ms1_scan_idxs)
        @debug_l1 "add_ms1_lookup_features!: no MS1 scans found, features all zero"
        return
    end

    # 1b. Precompute `scan_to_ms1[scan_id]` = nearest MS1 scan id, once per
    # file. Replaces a per-PSM `getRetentionTime` + `searchsortedfirst` (was
    # ~1-1.5 s/file on Astral) with a single array indexing per PSM.
    n_ms1 = length(ms1_scan_rts)
    scan_to_ms1 = Vector{Int32}(undef, n_scans)
    @inbounds for s in 1:n_scans
        scan_rt = Float32(getRetentionTime(spectra, s))
        pos = searchsortedfirst(ms1_scan_rts, scan_rt)
        scan_to_ms1[s] = if pos == 1
            Int32(ms1_scan_idxs[1])
        elseif pos > n_ms1
            Int32(ms1_scan_idxs[end])
        else
            d_after  = abs(ms1_scan_rts[pos]   - scan_rt)
            d_before = abs(ms1_scan_rts[pos-1] - scan_rt)
            d_before <= d_after ? Int32(ms1_scan_idxs[pos-1]) : Int32(ms1_scan_idxs[pos])
        end
    end

    scan_window_low = fill(Inf32, n_scans)
    scan_window_high = fill(-Inf32, n_scans)
    @inbounds for s in 1:n_scans
        center_mz = getCenterMz(spectra, s)
        isolation_width = getIsolationWidthMz(spectra, s)
        if !ismissing(center_mz) && !ismissing(isolation_width)
            width = Float32(isolation_width)
            if width > 0f0
                center = Float32(center_mz)
                half_width = width / 2f0
                scan_window_low[s] = center - half_width
                scan_window_high[s] = center + half_width
            end
        end
    end

    # 2. Per-precursor info
    precursors = getPrecursors(getSpecLib(search_context))
    prec_mzs       = getMz(precursors)
    prec_charges   = getCharge(precursors)
    prec_sulfurs   = getSulfurCount(precursors)
    iso_splines    = getIsoSplines(getSearchData(search_context)[1])

    starts = _scan_run_starts(psms[!, :scan_idx]::Vector{UInt32})
    n_runs = length(starts) - 1
    competition_inputs = _ms1_m0_peak_competition_inputs(psms, prec_mzs)
    if competition_inputs === nothing
        @debug_l1 "add_ms1_lookup_features!: missing frag1_int..frag8_int, skipping M0 peak competition"
    end

    parallel_foreach!(n_runs) do chunk
        _ms1_lookup_scan_runs!(
            psms, spectra,
            prec_mzs, prec_charges, prec_sulfurs, iso_splines,
            scan_to_ms1,
            scan_window_low, scan_window_high,
            Float32(ms1_ppm_tol), Float32(ms1_ppm_offset),
            C13_C12_MASS_DIFF_F32,
            starts, chunk,
            competition_inputs,
        )
    end
    return
end

"""
    add_chromatogram_features!(psms)

Per-precursor chromatogram-feature passes (MS1 + MS2-fragment). Calls
`_add_ms1_chromatogram_features!` (uses :ms1_m0_intensity, :ms1_m1_intensity)
and `_add_fragment_chromatogram_features!` (uses :frag1..8_int) over the
shared group structure built once. When `spectra` is supplied, groups are split
by both precursor and isolation window.

**PRECONDITION**: requires psms to be sorted by :precursor_idx (the
`_build_precursor_groups` fast path), AND requires
`add_ms1_lookup_features!` to have already populated the :ms1_m0_intensity
/ :ms1_m1_intensity columns.
"""
function add_chromatogram_features!(
    psms::DataFrame,
    spectra::Union{Nothing, MassSpecData} = nothing;
    bitvec_rank_table = nothing,
)
    n = nrow(psms)
    n == 0 && return
    # Build grouping (perm + starts/ends) once; reuse across both passes.
    t_groups = @elapsed groups = hasproperty(psms, :precursor_idx) ?
        (spectra === nothing ?
            _build_precursor_groups(psms.precursor_idx) :
            _build_precursor_window_groups(
                psms.precursor_idx,
                psms.scan_idx,
                getCenterMzs(spectra),
                getIsolationWidthMzs(spectra),
            )) :
        nothing
    t_ms1_chrom = @elapsed _add_ms1_chromatogram_features!(psms; groups=groups)
    t_frag_chrom = @elapsed _add_fragment_chromatogram_features!(
        psms; groups=groups, bitvec_rank_table=bitvec_rank_table)
    @debug_l1 "  chrom-feature passes (n_psms=$n): " *
               "groups=$(round(t_groups, digits=2))s  " *
               "ms1_chrom=$(round(t_ms1_chrom, digits=2))s  " *
               "frag_chrom=$(round(t_frag_chrom, digits=2))s"
    return
end

"""
    _add_ms1_chromatogram_features!(psms)

Phase 2 (B1): per-precursor MS1 chromatogram correlation features. Operates on
the per-PSM `ms1_m0_intensity`, `ms1_m1_intensity`, `weight`, and `irt` columns
populated by the deconv pipeline + Phase 1 MS1 lookup. Adds:

- `ms1_corr_weight_m0` — Pearson correlation between the per-precursor weight
  chromatogram and the MS1 M0 intensity chromatogram (both sampled at the MS2
  scans where the precursor has a PSM). Real precursors: high (~1).
- `ms1_corr_m0_m1` — Pearson(M0 chrom, M+1 chrom). Real isotope envelopes
  co-elute; chimeric MS1 noise does not.
- `ms1_apex_offset_irt` — |this_psm_irt − arg-irt-max(M0 chrom)|; how far this
  scan is from the M0 elution apex.
- `ms1_weight_apex_to_m0_apex_irt` — |arg-irt-max(weight) − arg-irt-max(M0)|;
  agreement between MS2 and MS1 apex location. Real: 0; chimeric: large.
"""
function _add_ms1_chromatogram_features!(psms::DataFrame;
        groups::Union{Nothing,Tuple{Vector{Int},Vector{UInt32},Vector{UInt32}}} = nothing)
    n = nrow(psms)
    psms[!, :ms1_corr_weight_m0]            = zeros(Float32, n)
    psms[!, :ms1_corr_m0_m1]                = zeros(Float32, n)
    psms[!, :ms1_apex_offset_irt]           = zeros(Float32, n)
    psms[!, :ms1_weight_apex_to_m0_apex_irt]= zeros(Float32, n)
    n == 0 && return

    # Required columns
    if !all(c -> hasproperty(psms, c), (:precursor_idx, :ms1_m0_intensity, :ms1_m1_intensity, :weight, :irt_obs))
        @debug_l1 "_add_ms1_chromatogram_features!: missing required columns, skipping"
        return
    end

    prec = psms.precursor_idx
    m0   = psms.ms1_m0_intensity
    m1   = psms.ms1_m1_intensity
    w    = psms.weight
    irt  = psms.irt_obs

    # Reuse the shared precursor grouping (perm + starts/ends) if the caller
    # has already computed it; otherwise build it locally.
    perm, starts, ends = groups === nothing ?
        _build_precursor_groups(prec) :
        groups
    n_prec = length(starts)

    # Per-thread scratch reused across precursors (matches the fragment path):
    # avoids 4 fresh Vector allocations per precursor (millions of small allocs).
    # maxthreadid() (not nthreads()) because threadid() can exceed the default
    # pool size; safe under :static where threadid() is stable within the body.
    nthr = Threads.maxthreadid()
    vm0_scratch  = [Float32[] for _ in 1:nthr]
    vm1_scratch  = [Float32[] for _ in 1:nthr]
    vw_scratch   = [Float32[] for _ in 1:nthr]
    virt_scratch = [Float32[] for _ in 1:nthr]

    Threads.@threads :static for p in 1:n_prec
        @inbounds begin
            i_start = Int(starts[p])
            i_end   = Int(ends[p])
            npts    = i_end - i_start + 1
            npts < 2 && continue

            # Extract per-precursor chromatograms (m0, m1, weight, iRT) into
            # reused per-thread scratch.
            tid = Threads.threadid()
            v_m0  = vm0_scratch[tid];  resize!(v_m0, npts)
            v_m1  = vm1_scratch[tid];  resize!(v_m1, npts)
            v_w   = vw_scratch[tid];   resize!(v_w, npts)
            v_irt = virt_scratch[tid]; resize!(v_irt, npts)
            for k in 1:npts
                i_orig = perm[i_start + k - 1]
                v_m0[k]  = m0[i_orig]
                v_m1[k]  = m1[i_orig]
                v_w[k]   = Float32(w[i_orig])
                v_irt[k] = Float32(irt[i_orig])
            end

            # Two Pearson correlations (uses the same top-level _frag_pcor
            # helper defined for the fragment-chromatogram path). The third
            # variant (c_wm1 = weight-vs-M1) was dropped 2026-05-19 — never
            # in any feature list since 143d6b87.
            c_wm0 = _frag_pcor(v_w,  v_m0)
            c_m01 = _frag_pcor(v_m0, v_m1)

            # Apex of M0 and weight (arg-max). First-occurrence on ties.
            ai_m0 = 1; vmax_m0 = v_m0[1]
            ai_w  = 1; vmax_w  = v_w[1]
            for k in 2:npts
                if v_m0[k] > vmax_m0; vmax_m0 = v_m0[k]; ai_m0 = k; end
                if v_w[k]  > vmax_w;  vmax_w  = v_w[k];  ai_w  = k; end
            end
            irt_apex_m0 = v_irt[ai_m0]
            irt_apex_w  = v_irt[ai_w]
            weight_apex_to_m0 = abs(irt_apex_w - irt_apex_m0)

            # Scatter outputs back to original row indices. Bind the columns to
            # typed locals first — `psms.col` is type-unstable, so writing through
            # it boxes a value per row.
            out_wm0  = psms.ms1_corr_weight_m0::Vector{Float32}
            out_m01  = psms.ms1_corr_m0_m1::Vector{Float32}
            out_aoff = psms.ms1_apex_offset_irt::Vector{Float32}
            out_w2m0 = psms.ms1_weight_apex_to_m0_apex_irt::Vector{Float32}
            for k in 1:npts
                i_orig = perm[i_start + k - 1]
                out_wm0[i_orig]  = c_wm0
                out_m01[i_orig]  = c_m01
                out_aoff[i_orig] = abs(v_irt[k] - irt_apex_m0)
                out_w2m0[i_orig] = weight_apex_to_m0
            end
        end
    end
    return
end

"""
    _add_fragment_chromatogram_features!(psms)

Per-precursor MS2 fragment chromatogram features. For each precursor, builds
8 fragment chromatograms (`frag1_int .. frag8_int` indexed by MS2 scan, captured
in `Score!` from `MainUnscoredPSM`) plus the deconv weight chromatogram, then
computes:

- `frag_corr_top1_top2`        Pearson(rank-1 frag chrom, rank-2 frag chrom)
- `frag_corr_top1_top3`        Pearson(rank-1, rank-3)
- `frag_corr_top1_weight`      Pearson(rank-1, weight chrom)
- `frag_corr_mean_pairwise`    Mean **Spearman** over all 15 pairs of (frag_i, frag_j) chromatograms (won the Pearson A/B; intensity-scale is not informative for frag-vs-frag)
- `frag_corr_min_pairwise`     Min Pearson over the same pairs (catches single contaminated fragment)
- `frag_corr_top3_weight`      Mean Pearson(rank-1..3 chrom, weight chrom)
- `frag_apex_dispersion_irt`   Std-dev of arg-max iRT across the 8 fragments (real: tight; chimeric: wide)
- `n_correlated_fragments`     Count of fragments with Pearson(frag, weight) > 0.7

Validated 2026-05-10 to add ~+2,088 IDs at q≤.01 vs MS1-only baseline (Olsen
Exploris one-file, entrap1, paired EFDR ~0.0107). Mechanism is the same as
MS1Corr — co-elution of multiple ions of the same precursor is the strongest
discriminative signal in DIA, only here we exploit it among MS2 fragments
rather than MS1 isotopes.
"""
# Group rows by :precursor_idx via a single sortperm + boundary enumeration.
# Returns (perm, starts, ends) all UInt32-keyed so downstream code can iterate
# `for p in 1:length(starts)` over contiguous per-precursor row ranges.
# Both _add_ms1_chromatogram_features! and _add_fragment_chromatogram_features!
# need this; computing once in prepare_psm_features! and passing the result
# saves one sortperm pass per file (HIGH-confidence cleanup, IDs unchanged).
function _build_precursor_groups(prec_col::AbstractVector)
    n = length(prec_col)
    if n == 0
        return Vector{Int}(), Vector{UInt32}(), Vector{UInt32}()
    end

    # Fast path: input already sorted by precursor_idx (the post-2026-05-19
    # MainSearch invariant, set by permute_psms_by_precursor_idx! in
    # process_file!). Skip the sortperm; perm is the identity permutation
    # and we enumerate boundaries in a single sequential scan.
    if issorted(prec_col)
        perm = collect(1:n)
        n_prec = 1
        @inbounds for i in 2:n
            if prec_col[i] != prec_col[i-1]
                n_prec += 1
            end
        end
        starts = Vector{UInt32}(undef, n_prec)
        ends   = Vector{UInt32}(undef, n_prec)
        @inbounds let p = 1
            starts[1] = UInt32(1)
            for i in 2:n
                if prec_col[i] != prec_col[i-1]
                    ends[p] = UInt32(i - 1)
                    p += 1
                    starts[p] = UInt32(i)
                end
            end
            ends[n_prec] = UInt32(n)
        end
        return perm, starts, ends
    end

    # Slow path (input not sorted): full sortperm + boundary enumeration.
    pkeys = Vector{UInt32}(undef, n)
    @inbounds for i in 1:n
        pkeys[i] = UInt32(prec_col[i])
    end
    perm = sortperm(pkeys)

    n_prec = 0
    @inbounds let i = 1
        while i <= n
            n_prec += 1
            cur = pkeys[perm[i]]
            while i < n && pkeys[perm[i+1]] == cur
                i += 1
            end
            i += 1
        end
    end
    starts = Vector{UInt32}(undef, n_prec)
    ends   = Vector{UInt32}(undef, n_prec)
    @inbounds let i = 1, p = 0
        while i <= n
            p += 1
            starts[p] = UInt32(i)
            cur = pkeys[perm[i]]
            while i < n && pkeys[perm[i+1]] == cur
                i += 1
            end
            ends[p] = UInt32(i)
            i += 1
        end
    end
    return perm, starts, ends
end

@inline function _scan_window_key(
    scan_idx::Integer,
    center_mzs::AbstractVector,
    isolation_widths::AbstractVector,
)
    return (
        Float32(coalesce(center_mzs[scan_idx], 0f0)),
        Float32(coalesce(isolation_widths[scan_idx], 0f0)),
    )
end

function _build_precursor_window_groups(
    prec_col::AbstractVector,
    scan_col::AbstractVector{UInt32},
    center_mzs::AbstractVector,
    isolation_widths::AbstractVector,
)
    n = length(prec_col)
    if n == 0
        return Vector{Int}(), Vector{UInt32}(), Vector{UInt32}()
    end

    base_perm, base_starts, base_ends = _build_precursor_groups(prec_col)
    perm = Vector{Int}(undef, n)
    starts = UInt32[]
    ends = UInt32[]
    sizehint!(starts, length(base_starts))
    sizehint!(ends, length(base_starts))

    out_pos = 1
    @inbounds for group_idx in eachindex(base_starts)
        group_start = Int(base_starts[group_idx])
        group_stop = Int(base_ends[group_idx])
        first_row = base_perm[group_start]
        first_key = _scan_window_key(scan_col[first_row], center_mzs, isolation_widths)
        single_window = true
        for i in (group_start + 1):group_stop
            row = base_perm[i]
            if _scan_window_key(scan_col[row], center_mzs, isolation_widths) != first_key
                single_window = false
                break
            end
        end

        if single_window
            start_out = out_pos
            for i in group_start:group_stop
                perm[out_pos] = base_perm[i]
                out_pos += 1
            end
            push!(starts, UInt32(start_out))
            push!(ends, UInt32(out_pos - 1))
        else
            window_to_group = Dict{Tuple{Float32, Float32}, Int}()
            window_rows = Vector{Vector{Int}}()
            for i in group_start:group_stop
                row = base_perm[i]
                key = _scan_window_key(scan_col[row], center_mzs, isolation_widths)
                local_group = get(window_to_group, key, 0)
                if local_group == 0
                    push!(window_rows, Int[])
                    local_group = length(window_rows)
                    window_to_group[key] = local_group
                end
                push!(window_rows[local_group], row)
            end

            for rows in window_rows
                start_out = out_pos
                for row in rows
                    perm[out_pos] = row
                    out_pos += 1
                end
                push!(starts, UInt32(start_out))
                push!(ends, UInt32(out_pos - 1))
            end
        end
    end

    return perm, starts, ends
end

# Pearson correlation helper (Float32-ized, returns 0 on zero variance).
# Top-level (not closure) so the compiler specializes without boxing.
@inline function _frag_pcor(x::Vector{Float32}, y::Vector{Float32})
    m = length(x)
    m < 2 && return 0f0
    mx = mean(x); my = mean(y)
    sx = 0f0; sy = 0f0; sxy = 0f0
    @inbounds for j in 1:m
        dx = x[j]-mx; dy = y[j]-my
        sx += dx*dx; sy += dy*dy; sxy += dx*dy
    end
    d = sqrt(sx*sy)
    d > 0 ? Float32(sxy/d) : 0f0
end

function _fragment_rank_weights(n_frags::Integer)
    weights = Vector{Float32}(undef, n_frags)
    sum_w = 0f0
    @inbounds for r in 1:n_frags
        w = 1f0 / sqrt(Float32(r))
        weights[r] = w
        sum_w += w
    end
    scale = sum_w > 0f0 ? Float32(n_frags) / sum_w : 0f0
    @inbounds for r in 1:n_frags
        weights[r] *= scale
    end
    return weights
end

@inline function _positive_corr_summary(corrs::Vector{Float32}, rank_weights::Vector{Float32})
    strength = 0f0
    sumsq = 0f0
    @inbounds for i in eachindex(corrs)
        pos = min(max(corrs[i], 0f0), 1f0)
        weighted = rank_weights[i] * pos
        strength += weighted
        sumsq += weighted * weighted
    end
    effective_n = sumsq > 0f0 ? Float32((strength * strength) / sumsq) : 0f0
    return strength, effective_n
end

@inline function _bitvec_pattern_rank(rank_table, mask::Integer)
    rank_table === nothing && return UInt16(0)
    length(rank_table) == 256 || throw(ArgumentError("bitvec rank table must have length 256"))
    idx = Int(UInt16(mask) & 0x00ff) + 1
    @inbounds rank = UInt16(rank_table[idx])
    return rank
end

function _add_fragment_chromatogram_features!(psms::DataFrame;
        groups::Union{Nothing,Tuple{Vector{Int},Vector{UInt32},Vector{UInt32}}} = nothing,
        bitvec_rank_table = nothing)
    n = nrow(psms)
    # Only the features actually consumed by PRESCORE_FEATURES / ADVANCED_FEATURE_SET
    # are computed and written. Earlier prototypes emitted many more outputs
    # (frag_corr_mean_pairwise — 15 Spearman per precursor; frag_corr_min_pairwise
    # — 15 Pearson; frag_corr_top1_top2/_top3/_weight, frag_corr_top3_weight,
    # frag_corr_top5_weight, n_correlated_fragments / _50, frag_corr_best_weight,
    # frag_corr_to_median_mean — 1 sort + 6 Pearson per precursor; pred_obs_max/area
    # spectral_contrast/scribe). Dropping them removes ~50 correlation calls + 1
    # sort + ~10 vector allocations per precursor.
    psms[!, :frag_apex_dispersion_irt]    = zeros(Float32, n)
    psms[!, :n_correlated_fragments]      = zeros(UInt8,  n)  # threshold 0.7
    psms[!, :n_correlated_fragments_bitvec_rank] = zeros(UInt16, n)
    psms[!, :frag_corr_strength]          = zeros(Float32, n)
    psms[!, :frag_corr_effective_n]       = zeros(Float32, n)
    psms[!, :frag_corr_best_m0]           = zeros(Float32, n)
    psms[!, :delta_frame_peak_center]     = zeros(Float32, n)
    # :n_scans (per-precursor PSM count) — also a PRESCORE_FEATURES feature.
    # Populated inside the same threaded per-precursor loop below (saves a
    # second pass and the Dict-build that train_lgbm_and_select_best would
    # otherwise do over ~14M rows).
    psms[!, :n_scans]                     = ones(UInt32, n)   # default 1 for single-PSM precs
    n == 0 && return

    if !all(c -> hasproperty(psms, c), (:precursor_idx, :frag1_int, :frag2_int, :frag3_int,
                                        :frag4_int, :frag5_int, :frag6_int, :frag7_int,
                                        :frag8_int, :weight, :irt_obs))
        @debug_l1 "_add_fragment_chromatogram_features!: missing required columns, skipping"
        return
    end

    prec   = psms.precursor_idx
    weight = psms.weight
    irt    = psms.irt_obs
    f      = (psms.frag1_int, psms.frag2_int, psms.frag3_int,
              psms.frag4_int, psms.frag5_int, psms.frag6_int,
              psms.frag7_int, psms.frag8_int)
    has_m0 = hasproperty(psms, :ms1_m0_intensity)
    m0_int = has_m0 ? psms.ms1_m0_intensity : nothing
    n_scans_col = psms.n_scans::Vector{UInt32}

    # Reuse the shared precursor grouping (perm + starts/ends) if the caller
    # has already computed it; otherwise build it locally. Per-precursor row
    # order within the run does not affect any of the output features
    # (correlations + apex stats are order-independent).
    perm, starts, ends = groups === nothing ?
        _build_precursor_groups(prec) :
        groups
    n_prec = length(starts)
    rank_weights = _fragment_rank_weights(8)

    # Per-thread reusable scratch (resized per precursor) so the per-precursor
    # chromatogram extraction below does not allocate ~10 small vectors on every
    # iteration. Safe under :static scheduling: threadid() is stable within an
    # iteration and no two concurrent iterations share a thread. Size by
    # maxthreadid() (not nthreads()) because threadid() can exceed the default
    # pool size when an interactive threadpool is present.
    nthr = Threads.maxthreadid()
    F_scratch    = [[Float32[] for _ in 1:8] for _ in 1:nthr]
    W_scratch    = [Float32[] for _ in 1:nthr]
    IRT_scratch  = [Float32[] for _ in 1:nthr]
    cfw_scratch  = [Vector{Float32}(undef, 8) for _ in 1:nthr]
    vm0_scratch  = [Float32[] for _ in 1:nthr]
    apex_scratch = [Float32[] for _ in 1:nthr]

    # Parallel per-precursor walk. Each precursor writes to disjoint row indices
    # in the output columns; all input arrays are read-only.
    Threads.@threads :static for p in 1:n_prec
        @inbounds begin
            tid = Threads.threadid()
            i_start = Int(starts[p])
            i_end   = Int(ends[p])
            npts    = i_end - i_start + 1
            # Write :n_scans for every row in this precursor's group — done
            # BEFORE the npts<2 early exit so all groups get a value (default 1
            # from initialization is the right answer for the npts==1 case,
            # but we still write to be explicit and correct under any ordering).
            len_u32 = UInt32(npts)
            for k in 0:(npts-1)
                n_scans_col[perm[i_start + k]] = len_u32
            end
            npts < 2 && continue

            # Extract chromatograms for the 8 fragments + weight + iRT into
            # per-thread scratch (resized to this precursor's npts, fully
            # overwritten below — no stale data is read).
            F = F_scratch[tid]
            for r in 1:8
                v = F[r]
                resize!(v, npts)
                for k in 1:npts
                    v[k] = Float32(f[r][perm[i_start + k - 1]])
                end
            end
            W   = W_scratch[tid];   resize!(W, npts)
            IRT = IRT_scratch[tid]; resize!(IRT, npts)
            for k in 1:npts
                i_orig = perm[i_start + k - 1]
                W[k]   = Float32(weight[i_orig])
                IRT[k] = Float32(irt[i_orig])
            end

            has_signal = ntuple(r -> maximum(F[r]) > 0, 8)

            # Apex dispersion across fragments with signal (also feeds delta_frame).
            apex_irts = apex_scratch[tid]; empty!(apex_irts)
            for r in 1:8
                has_signal[r] || continue
                ai = 1; vmax = F[r][1]
                for k in 2:npts; if F[r][k] > vmax; vmax = F[r][k]; ai = k; end; end
                push!(apex_irts, IRT[ai])
            end
            apex_disp = length(apex_irts) >= 2 ? Float32(std(apex_irts)) : 0f0

            # E14: median fragment apex iRT − midpoint of precursor scan-window
            delta_frame = 0f0
            if !isempty(apex_irts)
                sorted_apex = sort(apex_irts)
                m = length(sorted_apex)
                median_apex = isodd(m) ?
                    sorted_apex[(m + 1) ÷ 2] :
                    (sorted_apex[m ÷ 2] + sorted_apex[m ÷ 2 + 1]) / 2f0
                irt_lo = IRT[1]; irt_hi = IRT[1]
                for k in 2:npts
                    vk = IRT[k]
                    if vk < irt_lo; irt_lo = vk; end
                    if vk > irt_hi; irt_hi = vk; end
                end
                delta_frame = Float32(median_apex - (irt_lo + irt_hi) / 2f0)
            end

            # Per-fragment weight correlation — feeds n_correlated_fragments (>0.7).
            # 2026-05-15: replaced n_correlated_fragments_90 (>0.9) — 0.7 threshold
            # was ~11× more informative in ScoringSearch Pass-1 LGBM gain on
            # 23-file Olsen.
            c_fw = cfw_scratch[tid]   # length 8, fully overwritten below
            for r in 1:8
                c_fw[r] = has_signal[r] ? _frag_pcor(F[r], W) : 0f0
            end
            n_corr_70 = UInt8(0)
            corr_mask = UInt16(0)
            for r in 1:8
                has_signal[r] || continue
                if c_fw[r] > 0.7f0
                    n_corr_70 += UInt8(1)
                    corr_mask |= UInt16(1) << (r - 1)
                end
            end
            corr_rank = _bitvec_pattern_rank(bitvec_rank_table, corr_mask)
            corr_strength, corr_effective_n = _positive_corr_summary(c_fw, rank_weights)

            # DIA-NN-style best fragment: rank r with the highest mean correlation
            # to the other top-8 fragments. 56 Pearson calls per precursor.
            # Anchors frag_corr_best_m0 = Pearson(best_frag, MS1 m0 chrom).
            best_r = 0
            best_consensus = typemin(Float32)
            for r in 1:8
                has_signal[r] || continue
                consensus = 0f0; npairs = 0
                for r2 in 1:8
                    (r2 == r || !has_signal[r2]) && continue
                    consensus += _frag_pcor(F[r], F[r2]); npairs += 1
                end
                avg = npairs > 0 ? consensus / npairs : typemin(Float32)
                if avg > best_consensus
                    best_consensus = avg
                    best_r = r
                end
            end
            c_best_m0 = 0f0
            if best_r > 0 && has_m0
                v_m0 = vm0_scratch[tid]; resize!(v_m0, npts)
                for k in 1:npts
                    v_m0[k] = Float32(m0_int[perm[i_start + k - 1]])
                end
                if maximum(v_m0) > 0
                    c_best_m0 = _frag_pcor(F[best_r], v_m0)
                end
            end

            # Scatter outputs to original row indices. Bind the columns to typed
            # locals first — `psms.col` is type-unstable, so writing through it
            # boxes a value per row.
            out_disp   = psms.frag_apex_dispersion_irt::Vector{Float32}
            out_ncorr  = psms.n_correlated_fragments::Vector{UInt8}
            out_rank   = psms.n_correlated_fragments_bitvec_rank::Vector{UInt16}
            out_str    = psms.frag_corr_strength::Vector{Float32}
            out_effn   = psms.frag_corr_effective_n::Vector{Float32}
            out_bestm0 = psms.frag_corr_best_m0::Vector{Float32}
            out_dframe = psms.delta_frame_peak_center::Vector{Float32}
            for k in 1:npts
                i_orig = perm[i_start + k - 1]
                out_disp[i_orig]   = apex_disp
                out_ncorr[i_orig]  = n_corr_70
                out_rank[i_orig]   = corr_rank
                out_str[i_orig]    = corr_strength
                out_effn[i_orig]   = corr_effective_n
                out_bestm0[i_orig] = c_best_m0
                out_dframe[i_orig] = delta_frame
            end
        end
    end
    return
end


"""
    add_scan_competition_features!(psms::DataFrame)

Add per-scan competition features:
- `weight_ratio_at_scan`: weight / max(weight) over all PSMs at the same scan
- `weight_rank_at_scan`: rank (1 = highest) of this PSM's weight at its scan

Scans with a single PSM get `weight_ratio_at_scan=1.0` and `weight_rank_at_scan=1`.

**PRECONDITION**: assumes `psms` is **contiguous-by-:scan_idx** — all rows
sharing a scan_idx form one contiguous block. This is the natural ordering
of the deconv output (see process_scans_fused.jl per-scan emit loop). Call
this BEFORE any sort that reorders rows (e.g. `permute_psms_by_precursor_idx!`).

This invariant lets us replace the Dict-of-vectors grouping (~7s/8 files)
with a single linear sweep for run boundaries, and lets the per-run loop
parallelize trivially over disjoint contiguous row ranges.
"""
function add_scan_competition_features!(psms::DataFrame)
    n = nrow(psms)
    weight_ratio = Vector{Float32}(undef, n)
    weight_rank  = Vector{UInt16}(undef, n)
    if n == 0
        psms[!, :weight_ratio_at_scan] = weight_ratio
        psms[!, :weight_rank_at_scan]  = weight_rank
        return
    end

    scan::Vector{UInt32}  = psms[!, :scan_idx]
    w::Vector{Float32}    = psms[!, :weight]

    # 1. Sweep for contiguous-run boundaries (O(n), single pass).
    t_runs = @elapsed begin
        starts = Int[1]
        sizehint!(starts, n ÷ 4)
        @inbounds for i in 2:n
            if scan[i] != scan[i-1]
                push!(starts, i)
            end
        end
        n_runs = length(starts)
        # Append sentinel so run k covers starts[k]:(starts[k+1]-1).
        push!(starts, n + 1)
    end

    # 2. Per-run rank/ratio. Each run writes disjoint output rows, so
    #    threading is contention-free. Use a thread-local sortperm buffer.
    t_loop = @elapsed parallel_foreach!(n_runs) do chunk
        @inbounds begin
            buf_w   = Float32[]
            buf_ord = Int[]
            for r in chunk
                s = starts[r]
                e = starts[r+1] - 1
                len = e - s + 1
                if len == 1
                    weight_ratio[s] = 1f0
                    weight_rank[s]  = UInt16(1)
                    continue
                end
                resize!(buf_w,   len)
                resize!(buf_ord, len)
                max_w = 0f0
                for k in 1:len
                    wk = w[s + k - 1]
                    buf_w[k]   = wk
                    buf_ord[k] = k
                    if wk > max_w
                        max_w = wk
                    end
                end
                # sortperm descending by buf_w (closure-key is fine for small len).
                sort!(view(buf_ord, 1:len), by = k -> -buf_w[k])
                for rank in 1:len
                    k = buf_ord[rank]
                    row = s + k - 1
                    weight_ratio[row] = max_w > 0 ? buf_w[k] / max_w : 1f0
                    weight_rank[row]  = UInt16(min(rank, typemax(UInt16)))
                end
            end
        end
    end

    psms[!, :weight_ratio_at_scan] = weight_ratio
    psms[!, :weight_rank_at_scan]  = weight_rank
    @debug_l1 "  scan_competition (n_psms=$n n_scans=$n_runs): " *
               "runs=$(round(t_runs*1000, digits=0))ms  per_run_loop=$(round(t_loop*1000, digits=0))ms"
    return
end

# LGBM_RECOVERY_FEATURES const removed 2026-05-18 — it was a stale legacy
# list (23 entries) that claimed to be "the full Phase-2 feature set" but
# was never referenced anywhere. The active feature lists are
# PRESCORE_FEATURES (above) and ADVANCED_FEATURE_SET in
# PrecursorScoringSearch/model_config.jl.
#
# Computed-but-currently-unwired columns retained for potential
# re-introduction (have tested individual lift in past A/B):
#   :tic           — Float16, log2(spectrum_tic). +4.2k IDs on
#                    MTAC_3P_Standard solo (commit 4287cee7), dropped
#                    in smart-composite reduction 143d6b87.
#   :matched_ratio — Float16, log2((b+y+1)/(pred_int_sum_m0+1)). Batch E
#                    feature (c734e62a); dropped during 143d6b87.

"""
    prepare_psm_features!(psms, params, search_context, ms_file_idx, spectra)

Per-file PSM feature computation pipeline. Adds every prescore feature
column consumed by the per-file LightGBM in MainSearch via a single
parallel pass.

The historical two-pass form (`add_search_columns!` + `add_features!`)
and the optional `prescore_only=false` isotope-capture-and-filter path
were collapsed 2026-05-19. The isotope-capture work now runs only later
in MainSearch's `phase2`, on the much smaller best-per-precursor frame
(~400k rows). Measurement showed pre-filter saved <1% of rows and cost
more than it saved due to an extra sort.
"""
function prepare_psm_features!(
    psms::DataFrame,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData,
) where {P<:MainSearchParameters}
    t0 = time()
    add_psm_features!(psms, search_context, spectra, ms_file_idx)
    r = s -> round(s, digits=3)
    @debug_l2 "  prepare_psm_features! ($(nrow(psms)) PSMs): $(r(time()-t0))s"
    return psms
end
