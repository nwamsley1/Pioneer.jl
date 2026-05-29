# PSM feature computation for MainSearch.
#
# `add_psm_features!` is a single parallel pass that adds every per-PSM
# feature column consumed by the per-file LightGBM. Merged 2026-05-19 from
# two separate passes (`add_search_columns!` + `add_features!`) — same total
# work, but one parallel_foreach! task spawn / join instead of two and one
# cache-warm read of `precursor_idx` / `scan_idx` per row instead of two.

@inline _count_mox(seq::AbstractString) = UInt8(count("Unimod:35", seq))
@inline _count_mox(::Missing)            = zero(UInt8)

const FRAGMENT_PEAK_INDEX_COLUMNS = (
    :frag1_peak_idx, :frag2_peak_idx, :frag3_peak_idx,
    :frag4_peak_idx, :frag5_peak_idx, :frag6_peak_idx,
)

"""
    add_psm_features!(psms, search_context, spectra, ms_file_idx)

Single parallel pass that allocates and fills every per-PSM feature column
needed by the per-file LightGBM. Reads the deconv-emit columns
(:precursor_idx, :scan_idx, :error, :total_ions) and produces:

- `:decoy`, `:target`, `:cv_fold`, `:charge`, `:charge2`
- `:rt`, `:err_norm`
- `:irt_obs`, `:irt_pred`, `:irt_error`
- `:missed_cleavage`, `:Mox`, `:sequence_length`, `:prec_mz`
- `:spectrum_peak_count`

Mox is computed lazily once per precursor via a length-`n_library` cache;
the count race is benign (same value written from any thread).
"""
function add_psm_features!(psms::DataFrame,
                           search_context::SearchContext,
                           spectra::MassSpecData,
                           ms_file_idx::Integer)
    precursors_lib   = getPrecursors(getSpecLib(search_context))
    prec_is_decoy    = getIsDecoy(precursors_lib)
    prec_charge      = getCharge(precursors_lib)
    prec_mz          = getMz(precursors_lib)
    prec_irt         = getIrt(precursors_lib)
    prec_length      = getLength(precursors_lib)
    structural_mods  = getStructuralMods(precursors_lib)
    prec_missed_clv  = getMissedCleavages(precursors_lib)
    scan_rts         = getRetentionTimes(spectra)
    masses           = getMzArrays(spectra)
    rt_to_irt        = getRtIrtModel(search_context, ms_file_idx)

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
            sequence_length[i]     = prec_length[prec_idx]
            prec_mzs[i]            = prec_mz[prec_idx]
            spectrum_peak_count[i] = length(masses[scan_idx])

            # Lazy per-precursor Mox.
            if !_mox_computed[prec_idx]
                _mox_vals[prec_idx]     = _count_mox(structural_mods[prec_idx])
                _mox_computed[prec_idx] = true
            end
            Mox[i] = _mox_vals[prec_idx]
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
    :total_ions, :missed_cleavage, :y_count, :weight, :gof,
    :charge, :Mox, :spectrum_peak_count, :sequence_length,
    :fitted_hellinger, :spectral_contrast,

    # Cross-precursor at-scan
    :weight_ratio_at_scan, :weight_rank_at_scan,

    # MS1 features (chromatogram-correlations dropped, only m0_corr kept)
    :ms1_m0_mass_err_ppm,
    :ms1_weight_apex_to_m0_apex_irt,
    :ms1_m0_intensity, :ms1_m1_intensity,
    :ms1_m1_to_m0_ratio, :ms1_m1_to_m0_pred,

    # Per-rank top-6 fragment trace intensities (kept; used by chromatogram features).
    # Each rank sums matched isotope peaks predicted at >=25% of the fragment's
    # most abundant isotope.
    # frag5_int, frag6_int dropped 2026-05-14 (Tier-5 drop_5to6) — 8-file Olsen
    # +794 IDs / +29 PGs. Still captured per-scan for chromatogram-correlation features.
    :frag1_int, :frag2_int, :frag3_int, :frag4_int,

    # Fragment-chromatogram correlations (frag_w_corr family dropped; pairs dropped)
    # frag_corr_mean_pairwise (Spearman) dropped 2026-05-13 — cross-dataset test
    # showed Olsen +471 IDs / MTAC +649 IDs. Saves 15 rank-sorts per precursor.
    :frag_apex_dispersion_irt,
    :n_correlated_fragments,
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
2. Search the MS1 spectrum for M0 and M+1 of the precursor (charge-aware
   m/z = prec_mz + iso × C13_C12_MASS_DIFF / charge), within ±ms1_ppm_tol ppm.
3. Populate MS1 point-lookup features per PSM:
   - `ms1_m0_mass_err_ppm` (Float32): |observed − theoretical| in ppm; 0 if no M0
   - `ms1_m0_intensity` (Float32): observed M0 intensity; 0 if not matched
   - `ms1_m1_intensity` (Float32): observed M+1 intensity; 0 if not matched

The intensities are also consumed by `_add_ms1_chromatogram_features!` to build
per-precursor M0 / M+1 / weight chromatograms and compute correlation features.

Earlier versions of this function also computed `ms1_iso_count`, `ms1_m{0,2,3}_matched`,
`ms1_log_iso_obs_pred[_m{2,3}]` (sulfur-aware predicted/observed isotope ratios), and
M+2/M+3 lookups. Feature-importance analysis (2026-05-10) showed the LGBM model recovers
all signal from raw `m0/m1_intensity` (presence ⇔ intensity > 0; ratio implicitly via
tree splits); the explicit features were redundant. Removing them costs ~174 IDs at
q≤.01 within run-to-run noise and slightly improves q≤.001 (+296).

If no MS1 scan is available (file has only MS2), intensity features remain zeroed.
"""
# Top-level helper for the MS1 per-PSM peak search. Operates on clean
# Vector{Float32} arrays (missings stripped during cache refresh — see
# _ms1_refresh_cache!). Uses the same `bsearch_hybrid` primitive as the
# MS2 fused-match path so we get a branchless binary search to the window
# lower bound, followed by a tight scalar walk for the in-window peaks.
@inline function _ms1_find_peak(mz::Vector{Float32}, intens::Vector{Float32},
                                 mz_target::Float32, ms1_ppm_tol::Float32)
    half_tol = mz_target * ms1_ppm_tol * 1f-6
    lo = mz_target - half_tol
    hi = mz_target + half_tol
    n_peaks = length(mz)
    n_peaks == 0 && return (false, 0f0, 0f0, 0)
    i0 = bsearch_hybrid(mz, lo, 1, n_peaks)
    best_diff = Inf32
    best_int  = 0f0
    best_mz   = 0f0
    best_idx  = 0
    found = false
    @inbounds @fastmath for j in i0:n_peaks
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
    return found, best_int, best_mz, best_idx
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

# Per-chunk worker. Each task processes its row range with a thread-local
# MS1-spectrum cache so the (getMzArray, getIntensityArray) lookups amortize
# across consecutive PSMs that share the same nearest MS1 scan.
#
# `scan_to_ms1` is a precomputed mapping (scan_id → nearest MS1 scan id),
# built once in the caller. This eliminates the per-PSM
# `getRetentionTime` + `searchsortedfirst` work.
function _ms1_lookup_chunk!(psms, spectra,
                             prec_mzs, prec_charges, prec_sulfurs, iso_splines,
                             scan_to_ms1::Vector{Int32},
                             ms1_ppm_tol::Float32, ms1_ppm_offset::Float32, isotope_spacing::Float32,
                             m0_peak_keys::Vector{UInt64},
                             chunk_start::Int, chunk_end::Int)
    # Per-task scratch buffers; reused (resize-only) across cache misses
    # within this chunk. Type-stable Vector{Float32}.
    cached_mz  = Vector{Float32}()
    cached_int = Vector{Float32}()
    cached_ms1_idx::Int = -1   # force first miss

    @inbounds for i in chunk_start:chunk_end
        ms1_idx = Int(scan_to_ms1[Int(psms.scan_idx[i])])

        if ms1_idx != cached_ms1_idx
            cached_ms1_idx = ms1_idx
            _ms1_refresh_cache!(cached_mz, cached_int,
                                getMzArray(spectra, ms1_idx),
                                getIntensityArray(spectra, ms1_idx))
        end

        pid = UInt32(psms.precursor_idx[i])
        prec_mz   = Float32(prec_mzs[pid])
        prec_chg  = Int(prec_charges[pid])
        prec_chg == 0 && (prec_chg = 1)

        # Shift theoretical target by the fitted per-file ppm bias so peaks
        # are searched where they actually land. Residual feature
        # ms1_m0_mass_err_ppm is computed against the shifted (bias-corrected)
        # target so its zero point matches the calibration.
        target_m0 = prec_mz * (1f0 + ms1_ppm_offset * 1f-6)
        target_m1 = target_m0 + isotope_spacing / Float32(prec_chg)

        m0_hit, m0_int, m0_mz, m0_peak_idx = _ms1_find_peak(cached_mz, cached_int, target_m0, ms1_ppm_tol)
        m1_hit, m1_int, _, _               = _ms1_find_peak(cached_mz, cached_int, target_m1, ms1_ppm_tol)

        psms.ms1_m0_intensity[i] = m0_hit ? m0_int : 0f0
        psms.ms1_m1_intensity[i] = m1_hit ? m1_int : 0f0
        if m0_hit
            m0_peak_keys[i] = (UInt64(UInt32(ms1_idx)) << 32) | UInt64(UInt32(m0_peak_idx))
            psms.ms1_m0_mass_err_ppm[i] = abs(m0_mz - target_m0) / target_m0 * 1f6
        end
        if m0_hit && m0_int > 0f0 && m1_hit
            obs_ratio = m1_int / m0_int
            sulf_raw = Int(prec_sulfurs[pid])
            sulf = clamp(sulf_raw, 0, 5)
            prec_mass_approx = prec_mz * Float32(prec_chg)
            p0 = max(iso_splines(Int64(sulf), Int64(0), prec_mass_approx), 1f-12)
            p1 = max(iso_splines(Int64(sulf), Int64(1), prec_mass_approx), 1f-12)
            pred_ratio = Float32(p1 / p0)
            psms.ms1_m1_to_m0_ratio[i]    = Float32(obs_ratio)
            psms.ms1_m1_to_m0_pred[i]     = pred_ratio
            # ms1_envelope_dev_log2 = log2(obs_ratio / pred_ratio) was dropped
            # in 143d6b87 (smart-composite reduction); compute deleted 2026-05-19.
        end
    end
    return
end

"""
    add_ms1_lookup_features!(psms, spectra, search_context, ms_file_idx)

Per-PSM MS1 spectrum lookup. Populates:
- `:ms1_m0_intensity`
- `:ms1_m1_intensity`
- `:ms1_m0_mass_err_ppm`
- `:ms1_m1_to_m0_ratio`
- `:ms1_m1_to_m0_pred`
- `:ms1_m0_peak_frag_intensity_fraction`
- `:ms1_m0_peak_n_precursors`
- `:scan_prec_mz_n_precursors`

For each PSM: finds the nearest MS1 scan by RT, extracts M0/M1 intensities
within a ppm tolerance window around the predicted precursor m/z, and
computes isotope ratios/sentinel values. The peak-fragment fraction groups
PSMs within the same MS2 scan by the exact matched MS1 M0 peak and reports
each PSM's top-six matched M0 fragment-intensity share within that group.
The companion counts report how many unique candidate precursors share that
same matched M0 peak, and how many unique candidate precursors in the same
MS2 spectrum have the exact same library precursor m/z.

**PRECONDITION (perf-only)**: the per-task MS1-spectrum cache hits ~once
per unique `scan_idx` in each chunk, so this is dramatically faster when
the input PSMs are contiguous-by-:scan_idx — the natural ordering of the
deconv output. Call BEFORE `permute_psms_by_precursor_idx!`. (Correctness
is unchanged either way; the cache just thrashes when precursor-sorted.)
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

    # 2. Per-precursor info
    precursors = getPrecursors(getSpecLib(search_context))
    prec_mzs       = getMz(precursors)
    prec_charges   = getCharge(precursors)
    prec_sulfurs   = getSulfurCount(precursors)
    iso_splines    = getIsoSplines(getSearchData(search_context)[1])

    isotope_spacing = C13_C12_MASS_DIFF_F32
    m0_peak_keys = zeros(UInt64, n)

    # Parallel per-PSM lookup. Each PSM writes to disjoint indices in the output
    # columns (ms1_m0_intensity[i], etc.) so no synchronization needed. We use
    # Threads.@spawn per chunk so each task keeps its own ms1-spectrum cache.
    # When input is contiguous-by-scan, the cache hits ~once per unique scan
    # in each chunk; when input is precursor-sorted the cache thrashes.
    nt = Threads.nthreads()
    chunk_size = max(1, cld(n, nt))
    @sync for t in 1:nt
        chunk_start = (t - 1) * chunk_size + 1
        chunk_start > n && break
        chunk_end = min(t * chunk_size, n)
        Threads.@spawn _ms1_lookup_chunk!(
            psms, spectra,
            prec_mzs, prec_charges, prec_sulfurs, iso_splines,
            scan_to_ms1,
            Float32(ms1_ppm_tol), Float32(ms1_ppm_offset), isotope_spacing,
            m0_peak_keys,
            chunk_start, chunk_end
        )
    end
    _add_m0_peak_fragment_competition_feature!(psms, m0_peak_keys, prec_mzs)
    return
end

function _add_m0_peak_fragment_competition_feature!(psms::DataFrame,
                                                    m0_peak_keys::AbstractVector{UInt64},
                                                    prec_mzs = nothing)
    n = nrow(psms)
    length(m0_peak_keys) == n || throw(ArgumentError("m0_peak_keys length must match psms rows"))
    frac_out = zeros(Float32, n)
    peak_count_out = zeros(UInt16, n)
    same_mz_count_out = zeros(UInt16, n)
    psms[!, :ms1_m0_peak_frag_intensity_fraction] = frac_out
    psms[!, :ms1_m0_peak_n_precursors] = peak_count_out
    psms[!, :scan_prec_mz_n_precursors] = same_mz_count_out
    n == 0 && return

    required = (:precursor_idx, :scan_idx, :frag1_int, :frag2_int, :frag3_int, :frag4_int, :frag5_int, :frag6_int)
    if !all(c -> hasproperty(psms, c), required)
        @debug_l1 "_add_m0_peak_fragment_competition_feature!: missing required columns, skipping"
        return
    end

    precursor_idx::Vector{UInt32} = psms[!, :precursor_idx]
    scan::Vector{UInt32} = psms[!, :scan_idx]
    f = (psms.frag1_int, psms.frag2_int, psms.frag3_int,
         psms.frag4_int, psms.frag5_int, psms.frag6_int)
    have_prec_mzs = prec_mzs !== nothing

    starts = Int[1]
    sizehint!(starts, n ÷ 4)
    @inbounds for i in 2:n
        if scan[i] != scan[i-1]
            push!(starts, i)
        end
    end
    n_runs = length(starts)
    push!(starts, n + 1)

    parallel_foreach!(n_runs) do chunk
        totals = Dict{UInt64, Float32}()
        peak_counts = Dict{UInt64, UInt16}()
        seen_peak_precursors = Set{UInt128}()
        same_mz_counts = Dict{UInt32, UInt16}()
        seen_mz_precursors = Set{UInt64}()
        @inbounds for r in chunk
            empty!(totals)
            empty!(peak_counts)
            empty!(seen_peak_precursors)
            empty!(same_mz_counts)
            empty!(seen_mz_precursors)
            s = starts[r]
            e = starts[r+1] - 1
            for row in s:e
                pid = precursor_idx[row]
                key = m0_peak_keys[row]
                if key != 0
                    seen_key = (UInt128(key) << 32) | UInt128(pid)
                    if !(seen_key in seen_peak_precursors)
                        push!(seen_peak_precursors, seen_key)
                        prev = get(peak_counts, key, UInt16(0))
                        peak_counts[key] = prev == typemax(UInt16) ? prev : prev + UInt16(1)
                    end
                end
                if have_prec_mzs && 1 <= Int(pid) <= length(prec_mzs)
                    mz_key = reinterpret(UInt32, Float32(prec_mzs[pid]))
                    seen_key = (UInt64(mz_key) << 32) | UInt64(pid)
                    if !(seen_key in seen_mz_precursors)
                        push!(seen_mz_precursors, seen_key)
                        prev = get(same_mz_counts, mz_key, UInt16(0))
                        same_mz_counts[mz_key] = prev == typemax(UInt16) ? prev : prev + UInt16(1)
                    end
                end
                key == 0 && continue
                frag_sum = f[1][row] + f[2][row] + f[3][row] + f[4][row] + f[5][row] + f[6][row]
                frag_sum > 0f0 || continue
                totals[key] = get(totals, key, 0f0) + frag_sum
            end
            for row in s:e
                pid = precursor_idx[row]
                key = m0_peak_keys[row]
                if key != 0
                    peak_count_out[row] = get(peak_counts, key, UInt16(0))
                    denom = get(totals, key, 0f0)
                    if denom > 0f0
                        frag_sum = f[1][row] + f[2][row] + f[3][row] + f[4][row] + f[5][row] + f[6][row]
                        frac_out[row] = frag_sum > 0f0 ? frag_sum / denom : 0f0
                    end
                end
                if have_prec_mzs && 1 <= Int(pid) <= length(prec_mzs)
                    mz_key = reinterpret(UInt32, Float32(prec_mzs[pid]))
                    same_mz_count_out[row] = get(same_mz_counts, mz_key, UInt16(0))
                end
            end
        end
    end
    return
end

"""
    _add_fragment_peak_competition_features!(psms)

Add scan-local MS2 fragment peak competition features for the top-6 matched
fragment trace anchors:

- `:frag_competition_num_unique_fragments` — number of distinct observed MS2
  peaks matched by this PSM's top-6 fragment trace anchors.
- `:frag_competition_mean_candidates` — mean number of unique precursor
  candidates in the same MS2 scan that also matched those observed peaks.

The raw `frag*_peak_idx` columns are transient: `drop_fragment_peak_index_columns!`
removes them after these features are computed.
"""
function _add_fragment_peak_competition_features!(psms::DataFrame)
    n = nrow(psms)
    num_unique_out = zeros(UInt8, n)
    mean_candidates_out = zeros(Float32, n)
    psms[!, :frag_competition_num_unique_fragments] = num_unique_out
    psms[!, :frag_competition_mean_candidates] = mean_candidates_out
    n == 0 && return

    required = (:precursor_idx, :scan_idx, FRAGMENT_PEAK_INDEX_COLUMNS...)
    if !all(c -> hasproperty(psms, c), required)
        @debug_l1 "_add_fragment_peak_competition_features!: missing required columns, skipping"
        return
    end

    precursor_idx::Vector{UInt32} = psms[!, :precursor_idx]
    scan::Vector{UInt32} = psms[!, :scan_idx]
    peak_cols = (
        psms[!, :frag1_peak_idx]::Vector{UInt32},
        psms[!, :frag2_peak_idx]::Vector{UInt32},
        psms[!, :frag3_peak_idx]::Vector{UInt32},
        psms[!, :frag4_peak_idx]::Vector{UInt32},
        psms[!, :frag5_peak_idx]::Vector{UInt32},
        psms[!, :frag6_peak_idx]::Vector{UInt32},
    )

    scan_order::Vector{Int} = sortperm(scan)
    starts = Int[1]
    sizehint!(starts, n ÷ 4)
    @inbounds for i in 2:n
        if scan[scan_order[i]] != scan[scan_order[i-1]]
            push!(starts, i)
        end
    end
    n_runs = length(starts)
    push!(starts, n + 1)

    parallel_foreach!(n_runs) do chunk
        peak_counts = Dict{UInt32, UInt16}()
        seen_peak_precursors = Set{UInt64}()
        row_peaks = UInt32[]
        sizehint!(row_peaks, 6)
        @inbounds for r in chunk
            empty!(peak_counts)
            empty!(seen_peak_precursors)
            s = starts[r]
            e = starts[r+1] - 1

            for pos in s:e
                row = scan_order[pos]
                pid = precursor_idx[row]
                for peak_col in peak_cols
                    peak = peak_col[row]
                    peak == UInt32(0) && continue
                    seen_key = (UInt64(peak) << 32) | UInt64(pid)
                    if !(seen_key in seen_peak_precursors)
                        push!(seen_peak_precursors, seen_key)
                        prev = get(peak_counts, peak, UInt16(0))
                        peak_counts[peak] = prev == typemax(UInt16) ? prev : prev + UInt16(1)
                    end
                end
            end

            for pos in s:e
                row = scan_order[pos]
                empty!(row_peaks)
                for peak_col in peak_cols
                    peak = peak_col[row]
                    peak == UInt32(0) && continue
                    seen = false
                    for existing in row_peaks
                        if existing == peak
                            seen = true
                            break
                        end
                    end
                    seen || push!(row_peaks, peak)
                end

                n_unique = length(row_peaks)
                n_unique == 0 && continue

                sum_candidates = UInt32(0)
                for peak in row_peaks
                    count_i = get(peak_counts, peak, UInt16(0))
                    sum_candidates += UInt32(count_i)
                end

                num_unique_out[row] = UInt8(min(n_unique, typemax(UInt8)))
                mean_candidates_out[row] = Float32(sum_candidates) / Float32(n_unique)
            end
        end
    end
    return
end

function drop_fragment_peak_index_columns!(psms::DataFrame)
    keep_cols = [col for col in propertynames(psms)
                 if !(col in FRAGMENT_PEAK_INDEX_COLUMNS)]
    length(keep_cols) == length(propertynames(psms)) && return nothing
    select!(psms, keep_cols)
    return nothing
end

"""
    add_chromatogram_features!(psms)

Per-precursor chromatogram-feature passes (MS1 + MS2-fragment). Calls
`_add_ms1_chromatogram_features!` (uses :ms1_m0_intensity, :ms1_m1_intensity)
and `_add_fragment_chromatogram_features!` (uses :frag1..6_int) over the
shared per-precursor group structure built once via `_build_precursor_groups`.

**PRECONDITION**: requires psms to be sorted by :precursor_idx (the
`_build_precursor_groups` fast path), AND requires
`add_ms1_lookup_features!` to have already populated the :ms1_m0_intensity
/ :ms1_m1_intensity columns.
"""
function add_chromatogram_features!(psms::DataFrame)
    n = nrow(psms)
    n == 0 && return
    # Build precursor grouping (perm + starts/ends) once; reuse across both passes.
    t_groups = @elapsed groups = hasproperty(psms, :precursor_idx) ?
        _build_precursor_groups(psms.precursor_idx) :
        nothing
    t_ms1_chrom = @elapsed _add_ms1_chromatogram_features!(psms; groups=groups)
    t_frag_chrom = @elapsed _add_fragment_chromatogram_features!(psms; groups=groups)
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

    Threads.@threads :static for p in 1:n_prec
        @inbounds begin
            i_start = Int(starts[p])
            i_end   = Int(ends[p])
            npts    = i_end - i_start + 1
            npts < 2 && continue

            # Extract per-precursor chromatograms (m0, m1, weight, iRT)
            v_m0  = Vector{Float32}(undef, npts)
            v_m1  = Vector{Float32}(undef, npts)
            v_w   = Vector{Float32}(undef, npts)
            v_irt = Vector{Float32}(undef, npts)
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

            # Scatter outputs back to original row indices
            for k in 1:npts
                i_orig = perm[i_start + k - 1]
                psms.ms1_corr_weight_m0[i_orig]             = c_wm0
                psms.ms1_corr_m0_m1[i_orig]                 = c_m01
                psms.ms1_apex_offset_irt[i_orig]            = abs(v_irt[k] - irt_apex_m0)
                psms.ms1_weight_apex_to_m0_apex_irt[i_orig] = weight_apex_to_m0
            end
        end
    end
    return
end

"""
    _add_fragment_chromatogram_features!(psms)

Per-precursor MS2 fragment chromatogram features. For each precursor, builds
6 fragment chromatograms (`frag1_int .. frag6_int` indexed by MS2 scan, captured
in `Score!` from `MainUnscoredPSM`) plus the deconv weight chromatogram, then
computes:

- `frag_corr_top1_top2`        Pearson(rank-1 frag chrom, rank-2 frag chrom)
- `frag_corr_top1_top3`        Pearson(rank-1, rank-3)
- `frag_corr_top1_weight`      Pearson(rank-1, weight chrom)
- `frag_corr_mean_pairwise`    Mean **Spearman** over all 15 pairs of (frag_i, frag_j) chromatograms (won the Pearson A/B; intensity-scale is not informative for frag-vs-frag)
- `frag_corr_min_pairwise`     Min Pearson over the same pairs (catches single contaminated fragment)
- `frag_corr_top3_weight`      Mean Pearson(rank-1..3 chrom, weight chrom)
- `frag_apex_dispersion_irt`   Std-dev of arg-max iRT across the 6 fragments (real: tight; chimeric: wide)
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

function _add_fragment_chromatogram_features!(psms::DataFrame;
        groups::Union{Nothing,Tuple{Vector{Int},Vector{UInt32},Vector{UInt32}}} = nothing)
    n = nrow(psms)
    # Only the 4 features actually consumed by PRESCORE_FEATURES / ADVANCED_FEATURE_SET
    # are computed and written. Earlier versions emitted 15 more outputs
    # (frag_corr_mean_pairwise — 15 Spearman per precursor; frag_corr_min_pairwise
    # — 15 Pearson; frag_corr_top1_top2/_top3/_weight, frag_corr_top3_weight,
    # frag_corr_top5_weight, n_correlated_fragments / _50, frag_corr_best_weight,
    # frag_corr_to_median_mean — 1 sort + 6 Pearson per precursor; pred_obs_max/area
    # spectral_contrast/scribe). Dropping them removes ~50 correlation calls + 1
    # sort + ~10 vector allocations per precursor.
    psms[!, :frag_apex_dispersion_irt]    = zeros(Float32, n)
    psms[!, :n_correlated_fragments]      = zeros(UInt8,  n)  # threshold 0.7
    psms[!, :frag_corr_best_m0]           = zeros(Float32, n)
    psms[!, :delta_frame_peak_center]     = zeros(Float32, n)
    # :n_scans (per-precursor PSM count) — also a PRESCORE_FEATURES feature.
    # Populated inside the same threaded per-precursor loop below (saves a
    # second pass and the Dict-build that train_lgbm_and_select_best would
    # otherwise do over ~14M rows).
    psms[!, :n_scans]                     = ones(UInt32, n)   # default 1 for single-PSM precs
    n == 0 && return

    if !all(c -> hasproperty(psms, c), (:precursor_idx, :frag1_int, :frag2_int, :frag3_int,
                                        :frag4_int, :frag5_int, :frag6_int, :weight, :irt_obs))
        @debug_l1 "_add_fragment_chromatogram_features!: missing required columns, skipping"
        return
    end

    prec   = psms.precursor_idx
    weight = psms.weight
    irt    = psms.irt_obs
    f      = (psms.frag1_int, psms.frag2_int, psms.frag3_int,
              psms.frag4_int, psms.frag5_int, psms.frag6_int)
    has_m0 = hasproperty(psms, :ms1_m0_intensity)
    m0_int = has_m0 ? psms.ms1_m0_intensity : nothing
    n_scans_col = psms.n_scans::Vector{UInt32}

    # Reuse the shared precursor grouping (perm + starts/ends) if the caller
    # has already computed it; otherwise build it locally. Per-precursor row
    # order within the run does not affect any of the 4 output features
    # (correlations + apex stats are order-independent).
    perm, starts, ends = groups === nothing ?
        _build_precursor_groups(prec) :
        groups
    n_prec = length(starts)

    # Parallel per-precursor walk. Each precursor writes to disjoint row indices
    # in the 5 output columns; all input arrays are read-only.
    Threads.@threads :static for p in 1:n_prec
        @inbounds begin
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

            # Extract chromatograms for the 6 fragments + weight + iRT
            F = Vector{Vector{Float32}}(undef, 6)
            for r in 1:6
                v = Vector{Float32}(undef, npts)
                for k in 1:npts
                    v[k] = Float32(f[r][perm[i_start + k - 1]])
                end
                F[r] = v
            end
            W   = Vector{Float32}(undef, npts)
            IRT = Vector{Float32}(undef, npts)
            for k in 1:npts
                i_orig = perm[i_start + k - 1]
                W[k]   = Float32(weight[i_orig])
                IRT[k] = Float32(irt[i_orig])
            end

            has_signal = ntuple(r -> maximum(F[r]) > 0, 6)

            # Apex dispersion across fragments with signal (also feeds delta_frame).
            apex_irts = Float32[]
            for r in 1:6
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
            c_fw = Vector{Float32}(undef, 6)
            for r in 1:6
                c_fw[r] = has_signal[r] ? _frag_pcor(F[r], W) : 0f0
            end
            n_corr_70 = UInt8(0)
            for r in 1:6
                has_signal[r] || continue
                if c_fw[r] > 0.7f0; n_corr_70 += UInt8(1); end
            end

            # DIA-NN-style best fragment: rank r with the highest mean correlation
            # to the other top-6 fragments. 30 Pearson calls per precursor.
            # Anchors frag_corr_best_m0 = Pearson(best_frag, MS1 m0 chrom).
            best_r = 0
            best_consensus = typemin(Float32)
            for r in 1:6
                has_signal[r] || continue
                consensus = 0f0; npairs = 0
                for r2 in 1:6
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
                v_m0 = Vector{Float32}(undef, npts)
                for k in 1:npts
                    v_m0[k] = Float32(m0_int[perm[i_start + k - 1]])
                end
                if maximum(v_m0) > 0
                    c_best_m0 = _frag_pcor(F[best_r], v_m0)
                end
            end

            # Scatter outputs to original row indices
            for k in 1:npts
                i_orig = perm[i_start + k - 1]
                psms.frag_apex_dispersion_irt[i_orig] = apex_disp
                psms.n_correlated_fragments[i_orig]    = n_corr_70
                psms.frag_corr_best_m0[i_orig]         = c_best_m0
                psms.delta_frame_peak_center[i_orig]   = delta_frame
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
