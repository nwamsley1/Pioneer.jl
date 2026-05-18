# PSM feature computation for MainSearch.
#
# Functions that add feature columns to PSM DataFrames for LightGBM scoring.

"""
    add_search_columns!(psms::DataFrame, scan_retention_time, prec_charge, prec_is_decoy, precursors; prescore_only=false)

Add essential columns to PSM DataFrame for search analysis.

# Added Columns
- Retention time
- Charge state
- Target/decoy status
- Ion counts
- Error metrics
- CV fold assignments
"""
function add_search_columns!(psms::DataFrame,
                        scan_retention_time::AbstractVector{Float32},
                        prec_charge::AbstractVector{UInt8},
                        prec_is_decoy::AbstractVector{Bool},
                        precursors::LibraryPrecursors;
                        prescore_only::Bool=false
)
    # Allocate new columns

    N = size(psms, 1);
    decoys = zeros(Bool, N);
    rt = zeros(Float32, N);
    total_ions::Vector{UInt16} = psms[!,:total_ions]
    err_norm = zeros(Float16, N);
    targets = zeros(Bool, N);
    prec_charges = zeros(UInt8, N);
    cv_fold = zeros(UInt8, N);
    scan_idxs::Vector{UInt32} = psms[!,:scan_idx]
    prec_idxs::Vector{UInt32} = psms[!,:precursor_idx]
    error::Vector{Float32} = psms[!,:error]

    parallel_foreach!(size(psms, 1)) do chunk
        for i in chunk
            prec_idx = prec_idxs[i]
            scan_idx = scan_idxs[i]

            decoys[i] = prec_is_decoy[prec_idx];
            targets[i] = decoys[i] == false
            rt[i] = Float32(scan_retention_time[scan_idx]);
            prec_charges[i] = prec_charge[prec_idx];
            err_norm[i] = Float16(min((2^min(error[i], 15f0))/max(total_ions[i], one(UInt16)), Float32(6e4)))
            cv_fold[i] = getCvFold(precursors, prec_idx)
        end
    end
    psms[!,:decoy] = decoys
    psms[!,:rt] = rt
    psms[!,:err_norm] = err_norm
    psms[!,:target] = targets
    psms[!,:charge] = prec_charges
    psms[!,:cv_fold] = cv_fold
    psms[!,:charge2] = Vector{UInt8}(psms[!, :charge] .== 2)
    if !prescore_only
        sort!(psms,:rt); #Sorting before grouping is critical.
    end
    return nothing
end

"""
    add_features!(psms::DataFrame, search_context::SearchContext, tic, masses, ms_file_idx, rt_to_irt_interp; prescore_only=false)

Add feature columns to PSMs for scoring and analysis.

# Added Features
- RT and iRT metrics
- Sequence properties
- Intensity metrics
- Spectrum characteristics
"""
function add_features!(psms::DataFrame,
                        search_context::SearchContext,
                                    tic::AbstractVector{Float32},
                                    masses::AbstractArray,
                                    ms_file_idx::Integer,
                                    rt_to_irt_interp::RtConversionModel;
                                    prescore_only::Bool=false
                                    )

    precursors_lib = getPrecursors(getSpecLib(search_context))
    structural_mods = getStructuralMods(precursors_lib)
    prec_sequence = getSequence(precursors_lib)
    prec_mz = getMz(precursors_lib)
    prec_irt = getIrt(precursors_lib)
    prec_charge = getCharge(precursors_lib)
    entrap_group_ids = getEntrapmentGroupId(precursors_lib)
    precursor_missed_cleavage = getMissedCleavages(precursors_lib)
    prec_length = getLength(precursors_lib)
    ###########################
    #Allocate new columns
    N = size(psms, 1)
    irt_obs = zeros(Float32, N)
    irt_pred = zeros(Float32, N)
    irt_error = zeros(Float32, N)
    missed_cleavage = zeros(UInt8, N);
    spectrum_peak_count = zeros(Float16, N);
    Mox = zeros(UInt8, N);
    sequence_length = zeros(UInt8, N);
    # Sequence composition counts. Histidine + proline known to influence
    # ionization/fragmentation behavior. Lazy-populated like Mox.
    his_count = zeros(UInt8, N)
    pro_count = zeros(UInt8, N)
    lys_count = zeros(UInt8, N)
    arg_count = zeros(UInt8, N)

    # prec_mz is in PRESCORE_FEATURES (per-file LGBM) and the ScoringSearch
    # feature set, so it's allocated unconditionally.
    prec_mzs = zeros(Float32, N)

    # Columns only needed for Phase 2 (full feature set)
    if !prescore_only
        entrap_group_id = zeros(UInt8, N)
        adjusted_intensity_explained = zeros(Float16, N);
        prec_charges = zeros(UInt8, N)
        TIC = zeros(Float16, N);
    end

    precursor_idx::Vector{UInt32} = psms[!,:precursor_idx]
    scan_idx::Vector{UInt32} = psms[!,:scan_idx]
    rt::Vector{Float32} = psms[!,:rt]
    log2_intensity_explained = psms[!,:log2_intensity_explained]::Vector{Float16}

    function countMOX(seq::String)
        return UInt8(count("Unimod:35", seq))
    end

    function countMOX(seq::Missing)
        return zero(UInt8)
    end

    function _count_aa(seq::AbstractString, aa::Char)
        c = 0
        for ch in seq
            if ch == aa; c += 1; end
        end
        return UInt8(min(c, 255))
    end

    # Lazy-populate per-precursor caches: Vector-backed for O(1) lookup.
    n_lib = length(prec_irt)
    _mox_vals = Vector{UInt8}(undef, n_lib)
    _mox_computed = falses(n_lib)
    _his_vals = Vector{UInt8}(undef, n_lib)
    _pro_vals = Vector{UInt8}(undef, n_lib)
    _lys_vals = Vector{UInt8}(undef, n_lib)
    _arg_vals = Vector{UInt8}(undef, n_lib)
    _aa_computed = falses(n_lib)

    parallel_foreach!(size(psms, 1)) do chunk
        for i in chunk
            prec_idx = precursor_idx[i]
            irt_obs[i] = rt_to_irt_interp(rt[i])
            irt_pred[i] = prec_irt[prec_idx]
            irt_error[i] = abs(irt_obs[i] - irt_pred[i])
            missed_cleavage[i] = precursor_missed_cleavage[prec_idx]
            # Lazy Mox: compute once per precursor (benign race — same value)
            if !_mox_computed[prec_idx]
                _mox_vals[prec_idx] = countMOX(structural_mods[prec_idx])
                _mox_computed[prec_idx] = true
            end
            Mox[i] = _mox_vals[prec_idx]
            if !_aa_computed[prec_idx]
                seq = prec_sequence[prec_idx]
                _his_vals[prec_idx] = _count_aa(seq, 'H')
                _pro_vals[prec_idx] = _count_aa(seq, 'P')
                _lys_vals[prec_idx] = _count_aa(seq, 'K')
                _arg_vals[prec_idx] = _count_aa(seq, 'R')
                _aa_computed[prec_idx] = true
            end
            his_count[i] = _his_vals[prec_idx]
            pro_count[i] = _pro_vals[prec_idx]
            lys_count[i] = _lys_vals[prec_idx]
            arg_count[i] = _arg_vals[prec_idx]
            spectrum_peak_count[i] = length(masses[scan_idx[i]])
            sequence_length[i] = prec_length[prec_idx]

            prec_mzs[i] = prec_mz[prec_idx]

            if !prescore_only
                entrap_group_id[i] = entrap_group_ids[prec_idx]
                TIC[i] = Float16(log2(tic[scan_idx[i]]))
                adjusted_intensity_explained[i] = Float16(log2(TIC[i]) + log2_intensity_explained[i]);
                prec_charges[i] = prec_charge[prec_idx]
            end
        end
    end

    psms[!,:irt_obs] = irt_obs
    psms[!,:irt_pred] = irt_pred
    psms[!,:irt_error] = irt_error
    psms[!,:missed_cleavage] = missed_cleavage
    psms[!,:Mox] = Mox
    psms[!,:his_count] = his_count
    psms[!,:pro_count] = pro_count
    psms[!,:lys_count] = lys_count
    psms[!,:arg_count] = arg_count
    psms[!,:spectrum_peak_count] = spectrum_peak_count
    psms[!,:sequence_length] = sequence_length
    psms[!,:prec_mz] = prec_mzs

    if !prescore_only
        psms[!,:tic] = TIC
        psms[!,:adjusted_intensity_explained] = adjusted_intensity_explained
        psms[!,:charge] = prec_charges
        psms[!,:entrapment_group_id] = entrap_group_id
        psms[!,:ms_file_idx] .= ms_file_idx
    end
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
    :Mox, :spectrum_peak_count, :sequence_length,
    :fitted_hellinger,

    # Cross-precursor at-scan
    :weight_ratio_at_scan, :weight_rank_at_scan,

    # MS1 features (chromatogram-correlations dropped, only m0_corr kept)
    :ms1_m0_mass_err_ppm,
    :ms1_weight_apex_to_m0_apex_irt,
    :ms1_m0_intensity, :ms1_m1_intensity,
    :ms1_m1_to_m0_ratio, :ms1_m1_to_m0_pred,

    # Per-rank M0 fragment intensities (kept; used by chromatogram features)
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
   m/z = prec_mz + iso × NEUTRON_MASS / charge), within ±ms1_ppm_tol ppm.
3. Populate three features per PSM:
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

If no MS1 scan is available (file has only MS2), all features are zeroed.
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
    n_peaks == 0 && return (false, 0f0, 0f0)
    i0 = bsearch_hybrid(mz, lo, 1, n_peaks)
    best_diff = Inf32
    best_int  = 0f0
    best_mz   = 0f0
    found = false
    @inbounds @fastmath for j in i0:n_peaks
        m = mz[j]
        m > hi && break
        diff = abs(m - mz_target)
        if diff < best_diff
            best_diff = diff
            best_int  = intens[j]
            best_mz   = m
            found = true
        end
    end
    return found, best_int, best_mz
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
function _ms1_lookup_chunk!(psms, spectra,
                             prec_mzs, prec_charges, prec_sulfurs, iso_splines,
                             ms1_scan_idxs, ms1_scan_rts,
                             ms1_ppm_tol::Float32, NEUTRON::Float32,
                             chunk_start::Int, chunk_end::Int)
    # Per-task scratch buffers; reused (resize-only) across cache misses
    # within this chunk. Type-stable Vector{Float32}.
    cached_mz  = Vector{Float32}()
    cached_int = Vector{Float32}()
    cached_ms1_idx::Int = -1   # force first miss

    @inbounds for i in chunk_start:chunk_end
        scan_idx = Int(psms.scan_idx[i])
        scan_rt  = Float32(getRetentionTime(spectra, scan_idx))

        pos = searchsortedfirst(ms1_scan_rts, scan_rt)
        ms1_idx = if pos == 1
            ms1_scan_idxs[1]
        elseif pos > length(ms1_scan_rts)
            ms1_scan_idxs[end]
        else
            d_after  = abs(ms1_scan_rts[pos]   - scan_rt)
            d_before = abs(ms1_scan_rts[pos-1] - scan_rt)
            d_before <= d_after ? ms1_scan_idxs[pos-1] : ms1_scan_idxs[pos]
        end

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

        target_m0 = prec_mz
        target_m1 = prec_mz + NEUTRON / Float32(prec_chg)

        m0_hit, m0_int, m0_mz = _ms1_find_peak(cached_mz, cached_int, target_m0, ms1_ppm_tol)
        m1_hit, m1_int, _    = _ms1_find_peak(cached_mz, cached_int, target_m1, ms1_ppm_tol)

        psms.ms1_m0_intensity[i] = m0_hit ? m0_int : 0f0
        psms.ms1_m1_intensity[i] = m1_hit ? m1_int : 0f0
        if m0_hit
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
            psms.ms1_envelope_dev_log2[i] = log2(max(obs_ratio, 1f-6) / max(pred_ratio, 1f-6))
        end
    end
    return
end

function add_ms1_features!(psms::DataFrame,
                            spectra,
                            search_context,
                            ms_file_idx::Integer;
                            ms1_ppm_tol::Float32 = Float32(parse(Float64, get(ENV, "PIONEER_MS1_PPM_TOL", "100.0"))))
    n = nrow(psms)
    psms[!, :ms1_m0_mass_err_ppm]   = zeros(Float32, n)
    psms[!, :ms1_m0_intensity]      = zeros(Float32, n)
    psms[!, :ms1_m1_intensity]      = zeros(Float32, n)
    psms[!, :ms1_m1_to_m0_ratio]     = zeros(Float32, n)
    psms[!, :ms1_m1_to_m0_pred]      = zeros(Float32, n)
    psms[!, :ms1_envelope_dev_log2]  = zeros(Float32, n)
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
        @debug_l1 "add_ms1_features!: no MS1 scans found, features all zero"
        return
    end

    # 2. Per-precursor info
    precursors = getPrecursors(getSpecLib(search_context))
    prec_mzs       = getMz(precursors)
    prec_charges   = getCharge(precursors)
    prec_sulfurs   = getSulfurCount(precursors)
    iso_splines    = getIsoSplines(getSearchData(search_context)[1])

    NEUTRON = Float32(1.00335)

    # Parallel per-PSM lookup. Each PSM writes to disjoint indices in the output
    # columns (ms1_m0_intensity[i], etc.) so no synchronization needed. We use
    # Threads.@spawn per chunk so each task keeps its own ms1-spectrum cache
    # (consecutive PSMs in a chunk often map to the same nearest MS1 scan).
    nt = Threads.nthreads()
    chunk_size = max(1, cld(n, nt))
    @sync for t in 1:nt
        chunk_start = (t - 1) * chunk_size + 1
        chunk_start > n && break
        chunk_end = min(t * chunk_size, n)
        Threads.@spawn _ms1_lookup_chunk!(
            psms, spectra,
            prec_mzs, prec_charges, prec_sulfurs, iso_splines,
            ms1_scan_idxs, ms1_scan_rts,
            ms1_ppm_tol, NEUTRON,
            chunk_start, chunk_end
        )
    end

    # B1: per-precursor MS1 chromatogram features
    t_ms1_chrom = @elapsed _add_ms1_chromatogram_features!(psms)
    # B2: per-precursor MS2-fragment chromatogram features (uses frag1..6_int captured by Score!)
    t_frag_chrom = @elapsed _add_fragment_chromatogram_features!(psms)
    @user_info "  chrom-feature passes (file_idx=$ms_file_idx): " *
               "ms1_chrom=$(round(t_ms1_chrom, digits=2))s  " *
               "frag_chrom=$(round(t_frag_chrom, digits=2))s  " *
               "(n_psms=$(nrow(psms)))"
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
function _add_ms1_chromatogram_features!(psms::DataFrame)
    n = nrow(psms)
    psms[!, :ms1_corr_weight_m0]            = zeros(Float32, n)
    psms[!, :ms1_corr_m0_m1]                = zeros(Float32, n)
    psms[!, :ms1_corr_weight_m1]            = zeros(Float32, n)
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

    # Single sortperm on precursor_idx + UInt32 boundary enumeration —
    # mirrors the neighborhood / frag_chrom parallelization pattern.
    pkeys = Vector{UInt32}(undef, n)
    @inbounds for i in 1:n
        pkeys[i] = UInt32(prec[i])
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

            # Three Pearson correlations (uses the same top-level _frag_pcor
            # helper defined for the fragment-chromatogram path).
            c_wm0 = _frag_pcor(v_w,  v_m0)
            c_m01 = _frag_pcor(v_m0, v_m1)
            c_wm1 = _frag_pcor(v_w,  v_m1)

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
                psms.ms1_corr_weight_m1[i_orig]             = c_wm1
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

# n_correlated_fragments threshold mode.
#   raw (default): count fragments with Pearson r > FRAG_R_CUT (legacy 0.7).
#   tstat:         count fragments with t = r·√((n−2)/(1−r²)) > FRAG_T_CUT.
# Threshold mode + values are env-gated for A/B sweeps:
#   PIONEER_FRAG_NCORR_MODE  ∈ {"raw","tstat"}  default "raw"
#   PIONEER_FRAG_RCUT        Float           default 0.7  (raw mode)
#   PIONEER_FRAG_TCUT        Float           default 2.5  (tstat mode)
# Evaluated at call time (per file) so the env can change between Julia
# sessions without busting the precompile cache.
_frag_ncorr_mode() = get(ENV, "PIONEER_FRAG_NCORR_MODE", "raw")
_frag_rcut()       = parse(Float32, get(ENV, "PIONEER_FRAG_RCUT", "0.7"))
_frag_tcut()       = parse(Float32, get(ENV, "PIONEER_FRAG_TCUT", "2.5"))

# Returns true if the fragment-vs-weight Pearson `r` clears the active
# threshold for the active mode. `npts` is the chromatogram length.
@inline function _frag_corr_passes(r::Float32, npts::Int, mode::String,
                                    rcut::Float32, tcut::Float32)
    if mode == "tstat"
        npts < 3 && return false
        # r ≈ ±1 → t diverges; treat as pass for r > 0.
        if r >= 0.999f0; return true; end
        if r <= -0.999f0; return false; end
        denom = 1f0 - r*r
        denom <= 0f0 && return false
        t = r * sqrt(Float32(npts - 2) / denom)
        return t > tcut
    else
        return r > rcut
    end
end

function _add_fragment_chromatogram_features!(psms::DataFrame)
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

    # Single sortperm on precursor_idx (UInt32 key) — replaces the Dict + per-
    # precursor Vector{Int} bucketing. Per-precursor row order within the run
    # does not affect any of the 4 output features (correlations + apex stats
    # are order-independent).
    pkeys = Vector{UInt32}(undef, n)
    @inbounds for i in 1:n
        pkeys[i] = UInt32(prec[i])
    end
    perm = sortperm(pkeys)

    # Enumerate precursor-run boundaries: count, then exact-size UInt32 fill.
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

    # Snapshot env-gated n_correlated_fragments mode + thresholds once per call
    # so the inner parallel loop sees plain locals (not env-dict reads per row).
    ncorr_mode = _frag_ncorr_mode()
    rcut       = _frag_rcut()
    tcut       = _frag_tcut()

    # Parallel per-precursor walk. Each precursor writes to disjoint row indices
    # in the 4 output columns; all input arrays are read-only.
    Threads.@threads :static for p in 1:n_prec
        @inbounds begin
            i_start = Int(starts[p])
            i_end   = Int(ends[p])
            npts    = i_end - i_start + 1
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

            # Per-fragment weight correlation — feeds n_correlated_fragments.
            # Two threshold modes (env-gated, evaluated per-call):
            #   raw   (default, legacy): count fragments with r > 0.7.
            #   tstat: count fragments with t = r·√((n−2)/(1−r²)) > T_CUT.
            # The t-statistic mode normalizes by degrees of freedom so low-scan
            # precursors aren't free-counted (Pearson at n=2,3 carries no info).
            c_fw = Vector{Float32}(undef, 6)
            for r in 1:6
                c_fw[r] = has_signal[r] ? _frag_pcor(F[r], W) : 0f0
            end
            n_corr_70 = UInt8(0)
            for r in 1:6
                has_signal[r] || continue
                if _frag_corr_passes(c_fw[r], npts, ncorr_mode, rcut, tcut)
                    n_corr_70 += UInt8(1)
                end
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
    add_apex_distance_feature!(psms::DataFrame)

For each PSM, compute the iRT distance from this scan to the scan where the
SAME precursor had its highest deconv weight (its "apex"). PSMs at the apex
get value 0; PSMs far from the apex get a larger value. Real precursors
should have all their PSMs clustered around the apex; noise hits may be
scattered.

Adds column `:irt_dist_to_weight_apex`.
"""
function add_apex_distance_feature!(psms::DataFrame)
    n = nrow(psms)
    psms[!, :irt_dist_to_weight_apex] = zeros(Float32, n)
    # Batch A — edge-distance features (kept the 2 that scored; dropped 4 zero-gain).
    psms[!, :relative_position]      = zeros(Float32, n)
    psms[!, :dist_to_relative_center]= zeros(Float32, n)
    # Batch B — weight-chromatogram peak-shape stats (per-precursor, broadcast to all PSMs).
    psms[!, :weight_apex_relative_pos]    = zeros(Float32, n)
    psms[!, :weight_chrom_skewness]       = zeros(Float32, n)
    psms[!, :weight_chrom_kurtosis]       = zeros(Float32, n)
    psms[!, :apex_to_edge_weight_ratio]   = zeros(Float32, n)
    psms[!, :n_above_hm_left_of_apex]     = zeros(UInt8,   n)
    psms[!, :n_above_hm_right_of_apex]    = zeros(UInt8,   n)
    psms[!, :hm_asymmetry]                = zeros(Float32, n)
    # Batch C — log-quadratic Gaussian fit to weight chromatogram.
    psms[!, :weight_chrom_gaussian_r2]       = zeros(Float32, n)
    psms[!, :weight_chrom_gaussian_sigma]    = zeros(Float32, n)
    psms[!, :weight_chrom_gaussian_apex_irt] = zeros(Float32, n)
    # Batch D — DIA-NN ports. shape_{m2..p2} = 5-bin normalized weight profile
    # centered on each PSM's scan position k within its precursor's scan list.
    psms[!, :shape_m2] = zeros(Float32, n)
    psms[!, :shape_m1] = zeros(Float32, n)
    psms[!, :shape_0]  = zeros(Float32, n)
    psms[!, :shape_p1] = zeros(Float32, n)
    psms[!, :shape_p2] = zeros(Float32, n)
    # pBestCorrDelta analog: gof - max(gof over this precursor's scans). ≤ 0.
    psms[!, :gof_minus_max_gof_precursor] = zeros(Float32, n)
    # pTotCorrSum analog: log(gof / (sum(gof over precursor) + 1)). Fraction of total signal.
    psms[!, :gof_fraction_of_total_precursor] = zeros(Float32, n)
    n == 0 && return
    @assert hasproperty(psms, :irt_obs) ":irt_obs required"
    by_prec = Dict{eltype(psms.precursor_idx), Vector{Int}}()
    @inbounds for i in 1:n
        push!(get!(() -> Int[], by_prec, psms.precursor_idx[i]), i)
    end
    @inbounds for (_, idxs) in by_prec
        # Sort by irt for position features
        if length(idxs) > 1
            sort!(idxs, by = j -> Float32(psms.irt_obs[j]))
        end
        N = length(idxs)
        # Find apex (max weight) and apex position
        apex_i = idxs[1]; apex_w = Float32(psms.weight[apex_i])
        apex_k = 0
        for (k0, j) in enumerate(idxs)
            w = Float32(psms.weight[j])
            if w > apex_w; apex_w = w; apex_i = j; apex_k = k0 - 1; end
        end
        apex_irt = Float32(psms.irt_obs[apex_i])
        apex_rel = N > 1 ? Float32(apex_k) / Float32(N - 1) : 0.5f0

        # Batch B sub-pass: moments of weight chromatogram + half-max asymmetry.
        # Skewness = E[((w-μ)/σ)^3]; excess kurtosis = E[((w-μ)/σ)^4] - 3.
        # Zero variance → 0. Single-PSM precursors → 0 for all.
        wsum = 0f0; w2sum = 0f0
        for j in idxs; wj = Float32(psms.weight[j]); wsum += wj; w2sum += wj*wj; end
        mu  = wsum / Float32(N)
        var = max(w2sum / Float32(N) - mu*mu, 0f0)
        sigma = sqrt(var)
        skew = 0f0; kurt = 0f0
        if sigma > 1f-12 && N >= 3
            s3 = 0f0; s4 = 0f0
            for j in idxs
                z = (Float32(psms.weight[j]) - mu) / sigma
                s3 += z^3; s4 += z^4
            end
            skew = s3 / Float32(N)
            kurt = s4 / Float32(N) - 3f0
        end
        # Apex-to-edge ratio: max(w) / max(w_first, w_last). Single-scan → 1.
        edge_w = N >= 2 ? max(Float32(psms.weight[idxs[1]]), Float32(psms.weight[idxs[end]])) : apex_w
        apex_edge_ratio = edge_w > 1f-12 ? apex_w / edge_w : Float32(N == 1 ? 1f0 : 0f0)
        # Half-max asymmetry around apex
        half_max = 0.5f0 * apex_w
        n_left_hm = UInt8(0); n_right_hm = UInt8(0)
        for (k0, j) in enumerate(idxs)
            (Float32(psms.weight[j]) >= half_max) || continue
            k = k0 - 1
            if k < apex_k; n_left_hm  += UInt8(1)
            elseif k > apex_k; n_right_hm += UInt8(1)
            end
        end
        hm_total = Int(n_left_hm) + Int(n_right_hm)
        asym = hm_total > 0 ? Float32(abs(Int(n_left_hm) - Int(n_right_hm))) / Float32(hm_total) : 0f0

        # Batch C: log-quadratic Gaussian fit. y = log(w + eps) = a + b·x + c·x²
        # where x = irt. If c < 0, this is a valid Gaussian: σ = sqrt(-1/(2c)),
        # apex_irt = -b/(2c). Skip if N < 3 (system underdetermined).
        gauss_r2 = 0f0; gauss_sigma = 0f0; gauss_apex_irt = 0f0
        if N >= 3 && apex_w > 1f-12
            # Normal equations for [a,b,c]:
            # Σ1·a + Σx·b + Σx²·c = Σy
            # Σx·a + Σx²·b + Σx³·c = Σxy
            # Σx²·a + Σx³·b + Σx⁴·c = Σx²y
            S0 = Float64(N); S1 = 0.0; S2 = 0.0; S3 = 0.0; S4 = 0.0
            Sy = 0.0; Sxy = 0.0; Sx2y = 0.0
            Y = Vector{Float64}(undef, N)
            X = Vector{Float64}(undef, N)
            for (k0, j) in enumerate(idxs)
                xi = Float64(psms.irt_obs[j])
                yi = log(max(Float64(psms.weight[j]), 1e-6))
                X[k0] = xi; Y[k0] = yi
                S1 += xi; xi2 = xi*xi; S2 += xi2
                S3 += xi2*xi; S4 += xi2*xi2
                Sy += yi; Sxy += xi*yi; Sx2y += xi2*yi
            end
            # Solve 3x3 via Cramer's rule (cheap; N is tiny)
            M11=S0; M12=S1; M13=S2
            M21=S1; M22=S2; M23=S3
            M31=S2; M32=S3; M33=S4
            det_M = M11*(M22*M33 - M23*M32) - M12*(M21*M33 - M23*M31) + M13*(M21*M32 - M22*M31)
            if abs(det_M) > 1e-20
                det_a = Sy*(M22*M33 - M23*M32) - M12*(Sxy*M33 - M23*Sx2y) + M13*(Sxy*M32 - M22*Sx2y)
                det_b = M11*(Sxy*M33 - M23*Sx2y) - Sy*(M21*M33 - M23*M31) + M13*(M21*Sx2y - Sxy*M31)
                det_c = M11*(M22*Sx2y - Sxy*M32) - M12*(M21*Sx2y - Sxy*M31) + Sy*(M21*M32 - M22*M31)
                a = det_a / det_M; b = det_b / det_M; c = det_c / det_M
                # R² of the log-space fit
                ymean = Sy / S0
                ss_tot = 0.0; ss_res = 0.0
                for k0 in 1:N
                    yhat = a + b*X[k0] + c*X[k0]*X[k0]
                    ss_res += (Y[k0] - yhat)^2
                    ss_tot += (Y[k0] - ymean)^2
                end
                gauss_r2 = Float32(ss_tot > 1e-20 ? 1.0 - ss_res/ss_tot : 0.0)
                if c < 0
                    gauss_sigma    = Float32(sqrt(-1.0/(2.0*c)))
                    gauss_apex_irt = Float32(-b/(2.0*c))
                end
            end
        end

        # Batch D pre-computation: precursor-wide max_gof and sum_gof
        has_gof = hasproperty(psms, :gof)
        max_gof_prec = 0f0; sum_gof_prec = 0f0
        if has_gof
            for j in idxs
                g = Float32(psms.gof[j])
                if g > max_gof_prec; max_gof_prec = g; end
                sum_gof_prec += g
            end
        end

        for (k0, i) in enumerate(idxs)
            k = k0 - 1
            rel_pos = N > 1 ? Float32(k) / Float32(N - 1) : 0.5f0
            # Batch D shape_*: take 5 weights at positions k-2..k+2; pad with 0 if outside;
            # normalize so the 5 sum to 1.
            sh_total = 0f0
            sh_v1 = 0f0; sh_v2 = 0f0; sh_v3 = 0f0; sh_v4 = 0f0; sh_v5 = 0f0
            for off in -2:2
                kk = k + off
                if 0 <= kk <= N - 1
                    w = Float32(psms.weight[idxs[kk + 1]])
                    if     off == -2; sh_v1 = w
                    elseif off == -1; sh_v2 = w
                    elseif off ==  0; sh_v3 = w
                    elseif off ==  1; sh_v4 = w
                    else;             sh_v5 = w
                    end
                    sh_total += w
                end
            end
            if sh_total > 1f-12
                sh_v1 /= sh_total; sh_v2 /= sh_total; sh_v3 /= sh_total
                sh_v4 /= sh_total; sh_v5 /= sh_total
            end
            # Batch D pBestCorrDelta / pTotCorrSum analogs (using gof)
            this_gof = has_gof ? Float32(psms.gof[i]) : 0f0
            gof_minus_max = has_gof ? this_gof - max_gof_prec : 0f0
            gof_frac = (has_gof && sum_gof_prec > 0) ?
                Float32(log(max(Float64(this_gof) / (Float64(sum_gof_prec) + 1.0), 1e-12))) : 0f0

            psms.irt_dist_to_weight_apex[i]       = abs(Float32(psms.irt_obs[i]) - apex_irt)
            psms.relative_position[i]             = rel_pos
            psms.dist_to_relative_center[i]       = abs(0.5f0 - rel_pos)
            psms.shape_m2[i] = sh_v1
            psms.shape_m1[i] = sh_v2
            psms.shape_0[i]  = sh_v3
            psms.shape_p1[i] = sh_v4
            psms.shape_p2[i] = sh_v5
            psms.gof_minus_max_gof_precursor[i]     = gof_minus_max
            psms.gof_fraction_of_total_precursor[i] = gof_frac
            psms.weight_apex_relative_pos[i]      = apex_rel
            psms.weight_chrom_skewness[i]         = skew
            psms.weight_chrom_kurtosis[i]         = kurt
            psms.apex_to_edge_weight_ratio[i]     = apex_edge_ratio
            psms.n_above_hm_left_of_apex[i]       = n_left_hm
            psms.n_above_hm_right_of_apex[i]      = n_right_hm
            psms.hm_asymmetry[i]                  = asym
            psms.weight_chrom_gaussian_r2[i]       = gauss_r2
            psms.weight_chrom_gaussian_sigma[i]    = gauss_sigma
            psms.weight_chrom_gaussian_apex_irt[i] = gauss_apex_irt
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
"""
function add_scan_competition_features!(psms::DataFrame)
    n = nrow(psms)
    weight_ratio       = Vector{Float32}(undef, n)
    weight_rank        = Vector{UInt16}(undef, n)
    # New: ratio to the 2nd-best and 3rd-best competitor (∈ (0, 1] when present,
    # else 1.0). For the apex PSM at a scan, ratio_to_2nd = max_w / w_2nd ≥ 1
    # so we flip and report w_other / w_apex for symmetric interpretation
    # (1.0 = strongly dominant; <1 = competing precursors exist).
    weight_ratio_2nd   = Vector{Float32}(undef, n)
    weight_ratio_3rd   = Vector{Float32}(undef, n)
    n_competitors_50pct = Vector{UInt16}(undef, n)
    if n == 0
        psms[!, :weight_ratio_at_scan]   = weight_ratio
        psms[!, :weight_rank_at_scan]    = weight_rank
        psms[!, :weight_ratio_to_2nd_best] = weight_ratio_2nd
        psms[!, :weight_ratio_to_3rd_best] = weight_ratio_3rd
        psms[!, :n_competitors_50pct]    = n_competitors_50pct
        return
    end
    scan = psms[!, :scan_idx]
    w    = psms[!, :weight]
    by_scan = Dict{eltype(scan), Vector{Int}}()
    @inbounds for i in 1:n
        push!(get!(() -> Int[], by_scan, scan[i]), i)
    end
    @inbounds for (_, idxs) in by_scan
        weights = Float32[Float32(w[i]) for i in idxs]
        max_w = maximum(weights)
        order = sortperm(weights, rev=true)
        sorted_w = weights[order]   # descending
        w2 = length(sorted_w) >= 2 ? sorted_w[2] : 0f0
        w3 = length(sorted_w) >= 3 ? sorted_w[3] : 0f0
        # Count competitors at this scan with weight ≥ 50% of apex
        n_50 = UInt16(0)
        if max_w > 0
            half = 0.5f0 * max_w
            for wt in weights
                if wt >= half; n_50 += UInt16(1); end
            end
        end
        for (rank, j) in enumerate(order)
            i = idxs[j]
            weight_ratio[i] = max_w > 0 ? weights[j] / max_w : Float32(1)
            weight_rank[i]  = UInt16(min(rank, typemax(UInt16)))
            # The "ratio to 2nd/3rd best" is computed from this PSM's perspective:
            # for apex (rank=1) it's w_apex / w_2nd (or 1 if no 2nd); for others
            # it's their weight / w_2nd (the runner-up's perspective).
            weight_ratio_2nd[i] = w2 > 0 ? weights[j] / w2 : Float32(1)
            weight_ratio_3rd[i] = w3 > 0 ? weights[j] / w3 : Float32(1)
            n_competitors_50pct[i] = n_50
        end
    end
    psms[!, :weight_ratio_at_scan]    = weight_ratio
    psms[!, :weight_rank_at_scan]     = weight_rank
    psms[!, :weight_ratio_to_2nd_best] = weight_ratio_2nd
    psms[!, :weight_ratio_to_3rd_best] = weight_ratio_3rd
    psms[!, :n_competitors_50pct]     = n_competitors_50pct
    return
end

# Full feature set used in Phase 2 (ScoringSearch gets these via fold Arrow files)
const LGBM_RECOVERY_FEATURES = [
    # Core spectral quality
    :fitted_manhattan_distance, :max_matched_residual, :gof,
    :max_unmatched_residual, :poisson, :spectral_contrast, :err_norm,
    # Scores / weights
    :scribe, :weight, :log2_intensity_explained,
    # Ion counts
    :y_count, :b_count, :isotope_count, :total_ions, :total_ions_iso,
    # RT
    :irt_error,
    # Peptide properties
    :charge, :sequence_length, :missed_cleavage, :Mox, :prec_mz,
    # Other
    :tic, :best_rank, :matched_ratio, :spectrum_peak_count,
]

"""
    prepare_psm_features!(psms, params, search_context, ms_file_idx, spectra; prescore_only=false)

Reusable feature computation pipeline for PSMs. Called in both
prescore (per-file) and full feature computation modes.

Steps:
1. add_search_columns! — RT, charge, target, cv_fold, err_norm, total_ions
2. get_isotopes_captured! — precursor_fraction_transmitted (skipped for prescore)
3. Filter by fraction_transmitted (skipped for prescore)
4. add_features! — irt_error, tic, prec_mz, sequence_length, etc.
"""
function prepare_psm_features!(
    psms::DataFrame,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData;
    prescore_only::Bool=false
) where {P<:MainSearchParameters}
    t0 = time()

    # 1. Add basic search columns (RT, charge, target/decoy status, cv_fold)
    add_search_columns!(psms,
        getRetentionTimes(spectra),
        getCharge(getPrecursors(getSpecLib(search_context))),
        getIsDecoy(getPrecursors(getSpecLib(search_context))),
        getPrecursors(getSpecLib(search_context));
        prescore_only=prescore_only
    )
    t1 = time()

    if prescore_only
        # Skip isotope computation and fraction_transmitted filter for prescore path
        t2 = t1
        t3 = t1
    else
        # 2. Determine which precursor isotopes are captured in each scan's isolation window
        get_isotopes_captured!(
            psms,
            getQuadTransmissionModel(search_context, ms_file_idx),
            getSearchData(search_context),
            psms[!, :scan_idx],
            getCharge(getPrecursors(getSpecLib(search_context))),
            getMz(getPrecursors(getSpecLib(search_context))),
            getSulfurCount(getPrecursors(getSpecLib(search_context))),
            getCenterMzs(spectra),
            getIsolationWidthMzs(spectra)
        )
        t2 = time()

        # 3. Filter by fraction_transmitted
        to_remove = findall(psms[!, :precursor_fraction_transmitted] .< params.min_fraction_transmitted)
        deleteat!(psms, to_remove)
        t3 = time()
    end

    # 4. Add ML features (irt_error, tic, prec_mz, sequence_length, etc.)
    add_features!(
        psms,
        search_context,
        getTICs(spectra),
        getMzArrays(spectra),
        ms_file_idx,
        getRtIrtModel(search_context, ms_file_idx);
        prescore_only=prescore_only
    )
    t4 = time()

    r = s -> round(s, digits=3)
    @debug_l2 "  prepare_psm_features! ($(nrow(psms)) PSMs): $(r(t4-t0))s"

    return psms
end
