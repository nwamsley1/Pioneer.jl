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

    # prec_mz is in PRESCORE_FEATURES (per-file LGBM) and the ScoringSearch
    # feature set, so it's allocated unconditionally.
    prec_mzs = zeros(Float32, N)

    # Columns only needed for Phase 2 (full feature set)
    if !prescore_only
        entrap_group_id = zeros(UInt8, N)
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
            spectrum_peak_count[i] = length(masses[scan_idx[i]])
            sequence_length[i] = prec_length[prec_idx]

            prec_mzs[i] = prec_mz[prec_idx]

            if !prescore_only
                entrap_group_id[i] = entrap_group_ids[prec_idx]
                TIC[i] = Float16(log2(tic[scan_idx[i]]))
                prec_charges[i] = prec_charge[prec_idx]
            end
        end
    end

    psms[!,:irt_obs] = irt_obs
    psms[!,:irt_pred] = irt_pred
    psms[!,:irt_error] = irt_error
    psms[!,:missed_cleavage] = missed_cleavage
    psms[!,:Mox] = Mox
    psms[!,:spectrum_peak_count] = spectrum_peak_count
    psms[!,:sequence_length] = sequence_length
    psms[!,:prec_mz] = prec_mzs

    if !prescore_only
        # :tic is intentionally kept — historically tested with +4.2k IDs on
        # MTAC_3P_Standard (commit 4287cee7) before being dropped from feature
        # lists in the smart-composite reduction (143d6b87). Compute preserved
        # for re-introduction without code surgery.
        psms[!,:tic] = TIC
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

    # Build precursor grouping (perm + starts/ends) once; reuse across B1+B2.
    # Saves one sortperm + boundary-enum pass per file.
    t_groups = @elapsed groups = hasproperty(psms, :precursor_idx) ?
        _build_precursor_groups(psms.precursor_idx) :
        nothing
    # B1: per-precursor MS1 chromatogram features
    t_ms1_chrom = @elapsed _add_ms1_chromatogram_features!(psms; groups=groups)
    # B2: per-precursor MS2-fragment chromatogram features (uses frag1..6_int captured by Score!)
    t_frag_chrom = @elapsed _add_fragment_chromatogram_features!(psms; groups=groups)
    @user_info "  chrom-feature passes (file_idx=$ms_file_idx): " *
               "groups=$(round(t_groups, digits=2))s  " *
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
function _add_ms1_chromatogram_features!(psms::DataFrame;
        groups::Union{Nothing,Tuple{Vector{Int},Vector{UInt32},Vector{UInt32}}} = nothing)
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

    # Reuse the shared precursor grouping (perm + starts/ends) if the caller
    # has already computed it; otherwise build it locally. Per-precursor row
    # order within the run does not affect any of the 4 output features
    # (correlations + apex stats are order-independent).
    perm, starts, ends = groups === nothing ?
        _build_precursor_groups(prec) :
        groups
    n_prec = length(starts)

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
        for (rank, j) in enumerate(order)
            i = idxs[j]
            weight_ratio[i] = max_w > 0 ? weights[j] / max_w : Float32(1)
            weight_rank[i]  = UInt16(min(rank, typemax(UInt16)))
        end
    end
    psms[!, :weight_ratio_at_scan] = weight_ratio
    psms[!, :weight_rank_at_scan]  = weight_rank
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
