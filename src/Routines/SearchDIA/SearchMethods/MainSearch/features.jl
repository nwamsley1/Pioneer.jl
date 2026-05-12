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
    precursor_pair_idxs = getPairIdx(precursors_lib)
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

    # Columns only needed for Phase 2 (full feature set)
    if !prescore_only
        irt_diff = zeros(Float32, N)
        pair_idxs = zeros(UInt32, N)
        entrap_group_id = zeros(UInt8, N)
        adjusted_intensity_explained = zeros(Float16, N);
        prec_charges = zeros(UInt8, N)
        prec_mzs = zeros(Float32, N);
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

            if !prescore_only
                entrap_group_id[i] = entrap_group_ids[prec_idx]
                irt_diff[i] = abs(irt_obs[i] - prec_irt[prec_idx])
                TIC[i] = Float16(log2(tic[scan_idx[i]]))
                adjusted_intensity_explained[i] = Float16(log2(TIC[i]) + log2_intensity_explained[i]);
                prec_charges[i] = prec_charge[prec_idx]
                pair_idxs[i] = extract_pair_idx(precursor_pair_idxs, prec_idx)
                prec_mzs[i] = prec_mz[prec_idx];
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

    if !prescore_only
        psms[!,:irt_diff] = irt_diff
        psms[!,:tic] = TIC
        psms[!,:adjusted_intensity_explained] = adjusted_intensity_explained
        psms[!,:charge] = prec_charges
        psms[!,:pair_id] = pair_idxs
        psms[!,:prec_mz] = prec_mzs
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
    :fitted_manhattan_distance, :irt_error, :poisson, :err_norm,
    :total_ions, :total_ions_iso, :missed_cleavage, :y_count, :weight, :gof,
    :max_unmatched_residual, :max_matched_residual, :Mox, :spectrum_peak_count,
    :sequence_length,
    :fitted_hellinger,
    :weight_ratio_at_scan, :weight_rank_at_scan,
    :best_gof_3scan, :best_manhattan_3scan, :best_max_residual_3scan,
    :irt_dist_best_gof_3scan, :irt_dist_best_manhattan_3scan, :irt_dist_best_max_residual_3scan,
    :worst_max_residual_3scan, :worst_manhattan_3scan,
    :best_gof_5scan, :best_manhattan_5scan, :best_max_residual_5scan,
    :irt_dist_best_gof_5scan, :irt_dist_best_manhattan_5scan, :irt_dist_best_max_residual_5scan,
    :worst_max_residual_11scan, :worst_manhattan_11scan,
    :irt_dist_to_weight_apex,
    :ms1_m0_mass_err_ppm,
    :ms1_corr_weight_m0, :ms1_corr_m0_m1, :ms1_corr_weight_m1,
    :ms1_apex_offset_irt, :ms1_weight_apex_to_m0_apex_irt,
    :ms1_m0_intensity, :ms1_m1_intensity,
    # Isotope envelope features (2026-05-11)
    :ms1_m1_to_m0_ratio, :ms1_m1_to_m0_pred, :ms1_envelope_dev_log2,
    # Sequence composition
    :his_count, :pro_count, :lys_count, :arg_count,
    # Cross-precursor competition
    :weight_ratio_to_2nd_best, :weight_ratio_to_3rd_best, :n_competitors_50pct,
    # Per-rank M0 fragment intensities (from MainUnscoredPSM)
    :frag1_int, :frag2_int, :frag3_int, :frag4_int, :frag5_int, :frag6_int,
    # Per-precursor fragment chromatogram features
    :frag_corr_top1_top2, :frag_corr_top1_top3, :frag_corr_top1_weight,
    :frag_corr_mean_pairwise, :frag_corr_min_pairwise, :frag_corr_top3_weight,
    :frag_corr_top5_weight,
    :frag_apex_dispersion_irt, :n_correlated_fragments,
    :n_correlated_fragments_50, :n_correlated_fragments_90,
    :frag_corr_best_weight, :frag_corr_best_m0,
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
function add_ms1_features!(psms::DataFrame,
                            spectra,
                            search_context,
                            ms_file_idx::Integer;
                            ms1_ppm_tol::Float32 = Float32(parse(Float64, get(ENV, "PIONEER_MS1_PPM_TOL", "10.0"))))
    n = nrow(psms)
    psms[!, :ms1_m0_mass_err_ppm]   = zeros(Float32, n)
    # Per-PSM intensities (also consumed by the per-precursor chromatogram pass)
    psms[!, :ms1_m0_intensity]      = zeros(Float32, n)
    psms[!, :ms1_m1_intensity]      = zeros(Float32, n)
    # Isotope envelope features (re-added 2026-05-11). The predicted M+1/M0
    # ratio comes from the sulfur-aware iso_splines (averagine-like). Real
    # precursors match prediction; noise/chimeric peaks deviate.
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

    # Cache MS1 spectra m/z + intensity per scan we visit (avoids repeated lookups)
    cached_ms1_idx::Int = 0
    cached_mz       = nothing
    cached_int      = nothing

    @inbounds for i in 1:n
        scan_idx = Int(psms.scan_idx[i])
        scan_rt  = Float32(getRetentionTime(spectra, scan_idx))

        # Binary-search nearest MS1 scan by RT
        pos = searchsortedfirst(ms1_scan_rts, scan_rt)
        if pos == 1
            ms1_idx = ms1_scan_idxs[1]
        elseif pos > length(ms1_scan_rts)
            ms1_idx = ms1_scan_idxs[end]
        else
            d_after  = abs(ms1_scan_rts[pos]   - scan_rt)
            d_before = abs(ms1_scan_rts[pos-1] - scan_rt)
            ms1_idx = d_before <= d_after ? ms1_scan_idxs[pos-1] : ms1_scan_idxs[pos]
        end

        if ms1_idx != cached_ms1_idx
            cached_ms1_idx = ms1_idx
            cached_mz  = getMzArray(spectra, ms1_idx)
            cached_int = getIntensityArray(spectra, ms1_idx)
        end

        pid = UInt32(psms.precursor_idx[i])
        prec_mz   = Float32(prec_mzs[pid])
        prec_chg  = Int(prec_charges[pid])
        prec_chg == 0 && (prec_chg = 1)

        # Per-isotope: target m/z, search ± ms1_ppm_tol
        function _find_peak(mz_target::Float32)
            half_tol = mz_target * ms1_ppm_tol * 1f-6
            lo = mz_target - half_tol
            hi = mz_target + half_tol
            best_j = 0
            best_diff = Inf32
            best_int  = 0f0
            best_mz   = 0f0
            n_peaks = length(cached_mz)
            i0 = searchsortedfirst(cached_mz, lo)
            j = i0
            while j <= n_peaks && cached_mz[j] <= hi
                m = cached_mz[j]; ismissing(m) && (j+=1; continue)
                it = cached_int[j]; ismissing(it) && (j+=1; continue)
                diff = abs(Float32(m) - mz_target)
                if diff < best_diff
                    best_diff = diff
                    best_j = j
                    best_int = Float32(it)
                    best_mz  = Float32(m)
                end
                j += 1
            end
            return best_j > 0, best_int, best_mz
        end

        target_m0 = prec_mz
        target_m1 = prec_mz + NEUTRON / Float32(prec_chg)

        m0_hit, m0_int, m0_mz = _find_peak(target_m0)
        m1_hit, m1_int, _    = _find_peak(target_m1)

        psms.ms1_m0_intensity[i] = m0_hit ? m0_int : 0f0
        psms.ms1_m1_intensity[i] = m1_hit ? m1_int : 0f0
        if m0_hit
            psms.ms1_m0_mass_err_ppm[i] = abs(m0_mz - target_m0) / target_m0 * 1f6
        end
        # Isotope envelope ratio features
        if m0_hit && m0_int > 0f0 && m1_hit
            obs_ratio = m1_int / m0_int
            sulf = prec_sulfurs[pid]
            # iso_splines(sulfur, isotope, mass) gives P(isotope | sulfurs, mass).
            # Use prec_mz × charge as an approximate neutral mass for the spline.
            prec_mass_approx = prec_mz * Float32(prec_chg)
            p0 = max(iso_splines(sulf, UInt8(0), prec_mass_approx), 1f-12)
            p1 = max(iso_splines(sulf, UInt8(1), prec_mass_approx), 1f-12)
            pred_ratio = Float32(p1 / p0)
            psms.ms1_m1_to_m0_ratio[i]    = Float32(obs_ratio)
            psms.ms1_m1_to_m0_pred[i]     = pred_ratio
            psms.ms1_envelope_dev_log2[i] = log2(max(obs_ratio, 1f-6) / max(pred_ratio, 1f-6))
        end
    end

    # B1: per-precursor MS1 chromatogram features
    _add_ms1_chromatogram_features!(psms)
    # B2: per-precursor MS2-fragment chromatogram features (uses frag1..6_int captured by Score!)
    _add_fragment_chromatogram_features!(psms)
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

    # Index PSMs by precursor_idx (no DataFrame groupby — just a per-row pass)
    prec = psms.precursor_idx
    m0   = psms.ms1_m0_intensity
    m1   = psms.ms1_m1_intensity
    w    = psms.weight
    irt  = psms.irt_obs

    # Build per-precursor index buckets
    buckets = Dict{UInt32, Vector{Int}}()
    @inbounds for i in 1:n
        push!(get!(() -> Int[], buckets, UInt32(prec[i])), i)
    end

    # For each precursor, compute chrom features once, broadcast to all PSMs
    @inbounds for (_pid, idxs) in buckets
        length(idxs) < 2 && continue   # correlation undefined with <2 points
        # Extract per-precursor chromatograms
        npts = length(idxs)
        v_m0  = Vector{Float32}(undef, npts)
        v_m1  = Vector{Float32}(undef, npts)
        v_w   = Vector{Float32}(undef, npts)
        v_irt = Vector{Float32}(undef, npts)
        for (k, i) in enumerate(idxs)
            v_m0[k]  = m0[i]
            v_m1[k]  = m1[i]
            v_w[k]   = Float32(w[i])
            v_irt[k] = Float32(irt[i])
        end

        # Pearson correlations (need variance > 0 in each vector)
        function _pcor(x::Vector{Float32}, y::Vector{Float32})
            mx = mean(x); my = mean(y)
            sx = 0f0; sy = 0f0; sxy = 0f0
            @inbounds for j in 1:length(x)
                dx = x[j]-mx; dy = y[j]-my
                sx  += dx*dx; sy += dy*dy; sxy += dx*dy
            end
            d = sqrt(sx*sy)
            d > 0 ? sxy/d : 0f0
        end
        c_wm0 = _pcor(v_w,  v_m0)
        c_m01 = _pcor(v_m0, v_m1)
        c_wm1 = _pcor(v_w,  v_m1)

        # Apex of M0 and weight (arg-max). Use first-occurrence on ties.
        ai_m0 = 1; vmax_m0 = v_m0[1]
        ai_w  = 1; vmax_w  = v_w[1]
        @inbounds for k in 2:npts
            if v_m0[k] > vmax_m0; vmax_m0 = v_m0[k]; ai_m0 = k; end
            if v_w[k]  > vmax_w;  vmax_w  = v_w[k];  ai_w  = k; end
        end
        irt_apex_m0 = v_irt[ai_m0]
        irt_apex_w  = v_irt[ai_w]
        weight_apex_to_m0 = abs(irt_apex_w - irt_apex_m0)

        # Broadcast to all PSMs of this precursor
        for (k, i) in enumerate(idxs)
            psms.ms1_corr_weight_m0[i]             = c_wm0
            psms.ms1_corr_m0_m1[i]                 = c_m01
            psms.ms1_corr_weight_m1[i]             = c_wm1
            psms.ms1_apex_offset_irt[i]            = abs(v_irt[k] - irt_apex_m0)
            psms.ms1_weight_apex_to_m0_apex_irt[i] = weight_apex_to_m0
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
- `frag_corr_mean_pairwise`    Mean Pearson over all 15 pairs of (frag_i, frag_j) chromatograms
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
function _add_fragment_chromatogram_features!(psms::DataFrame)
    n = nrow(psms)
    psms[!, :frag_corr_top1_top2]      = zeros(Float32, n)
    psms[!, :frag_corr_top1_top3]      = zeros(Float32, n)
    psms[!, :frag_corr_top1_weight]    = zeros(Float32, n)
    psms[!, :frag_corr_mean_pairwise]  = zeros(Float32, n)
    psms[!, :frag_corr_min_pairwise]   = zeros(Float32, n)
    psms[!, :frag_corr_top3_weight]    = zeros(Float32, n)
    psms[!, :frag_apex_dispersion_irt] = zeros(Float32, n)
    psms[!, :n_correlated_fragments]      = zeros(UInt8,  n)  # threshold 0.7 (original)
    psms[!, :n_correlated_fragments_50]   = zeros(UInt8,  n)  # threshold 0.5
    psms[!, :n_correlated_fragments_90]   = zeros(UInt8,  n)  # threshold 0.9
    # DIA-NN-style "best fragment" reference: pick the fragment with the highest
    # sum of correlations to the other top-6 fragments (consensus). Then compute
    # weight↔best-frag and m0↔best-frag correlations (if MS1 m0 column is
    # available). Real precursors: weight and m0 chromatograms both track the
    # consensus best fragment closely; chimeric: best-frag is decoupled.
    psms[!, :frag_corr_best_weight]    = zeros(Float32, n)
    psms[!, :frag_corr_best_m0]        = zeros(Float32, n)
    psms[!, :frag_corr_top5_weight]    = zeros(Float32, n)  # mean of top-5 vs weight (top-K sweep)
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

    buckets = Dict{UInt32, Vector{Int}}()
    @inbounds for i in 1:n
        push!(get!(() -> Int[], buckets, UInt32(prec[i])), i)
    end

    # Pearson correlation helper (Float32-ized, returns 0 on zero variance)
    function _pcor(x::Vector{Float32}, y::Vector{Float32})
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

    @inbounds for (_pid, idxs) in buckets
        npts = length(idxs)
        npts < 2 && continue

        # Extract chromatograms for the 6 fragments + weight + iRT
        F = Vector{Vector{Float32}}(undef, 6)
        for r in 1:6
            v = Vector{Float32}(undef, npts)
            for (k, i) in enumerate(idxs); v[k] = Float32(f[r][i]); end
            F[r] = v
        end
        W   = Vector{Float32}(undef, npts)
        IRT = Vector{Float32}(undef, npts)
        for (k, i) in enumerate(idxs)
            W[k]   = Float32(weight[i])
            IRT[k] = Float32(irt[i])
        end

        has_signal = ntuple(r -> maximum(F[r]) > 0, 6)

        # Pairwise frag-vs-frag correlations (skip pairs where either side has no signal)
        pair_count = 0
        pair_sum   = 0f0
        pair_min   = 1f0
        any_pair = false
        for r1 in 1:6, r2 in (r1+1):6
            (has_signal[r1] && has_signal[r2]) || continue
            c = _pcor(F[r1], F[r2])
            pair_sum += c
            pair_count += 1
            if c < pair_min; pair_min = c; end
            any_pair = true
        end
        mean_pairwise = pair_count > 0 ? pair_sum / pair_count : 0f0
        min_pairwise  = any_pair ? pair_min : 0f0

        c_top1_top2  = (has_signal[1] && has_signal[2]) ? _pcor(F[1], F[2]) : 0f0
        c_top1_top3  = (has_signal[1] && has_signal[3]) ? _pcor(F[1], F[3]) : 0f0
        c_top1_weight = has_signal[1] ? _pcor(F[1], W) : 0f0

        # Mean correlation of top-3 fragments to weight chromatogram
        n_top3 = 0; s_top3 = 0f0
        for r in 1:3
            has_signal[r] || continue
            s_top3 += _pcor(F[r], W); n_top3 += 1
        end
        c_top3_w = n_top3 > 0 ? s_top3 / n_top3 : 0f0

        # Apex dispersion across fragments with signal
        apex_irts = Float32[]
        for r in 1:6
            has_signal[r] || continue
            ai = 1; vmax = F[r][1]
            for k in 2:npts; if F[r][k] > vmax; vmax = F[r][k]; ai = k; end; end
            push!(apex_irts, IRT[ai])
        end
        apex_disp = length(apex_irts) >= 2 ? Float32(std(apex_irts)) : 0f0

        # Per-fragment weight correlation (reused below for n_corr_* and best frag).
        c_fw = Vector{Float32}(undef, 6)
        for r in 1:6
            c_fw[r] = has_signal[r] ? _pcor(F[r], W) : 0f0
        end
        n_corr_50 = UInt8(0); n_corr_70 = UInt8(0); n_corr_90 = UInt8(0)
        for r in 1:6
            has_signal[r] || continue
            if c_fw[r] > 0.5f0; n_corr_50 += UInt8(1); end
            if c_fw[r] > 0.7f0; n_corr_70 += UInt8(1); end
            if c_fw[r] > 0.9f0; n_corr_90 += UInt8(1); end
        end

        # Mean correlation of top-5 fragments (rank 1..5) to weight chromatogram
        n_top5 = 0; s_top5 = 0f0
        for r in 1:5
            has_signal[r] || continue
            s_top5 += c_fw[r]; n_top5 += 1
        end
        c_top5_w = n_top5 > 0 ? s_top5 / n_top5 : 0f0

        # DIA-NN-style best fragment: pick the rank r with the highest sum of
        # pairwise correlations to the other top-6 fragments (consensus).
        best_r = 0
        best_consensus = typemin(Float32)
        for r in 1:6
            has_signal[r] || continue
            consensus = 0f0; npairs = 0
            for r2 in 1:6
                (r2 == r || !has_signal[r2]) && continue
                consensus += _pcor(F[r], F[r2]); npairs += 1
            end
            avg = npairs > 0 ? consensus / npairs : typemin(Float32)
            if avg > best_consensus
                best_consensus = avg
                best_r = r
            end
        end
        c_best_w = 0f0
        c_best_m0 = 0f0
        if best_r > 0
            c_best_w = c_fw[best_r]
            if has_m0
                v_m0 = Vector{Float32}(undef, npts)
                for (k, i) in enumerate(idxs); v_m0[k] = Float32(m0_int[i]); end
                if maximum(v_m0) > 0
                    c_best_m0 = _pcor(F[best_r], v_m0)
                end
            end
        end

        for i in idxs
            psms.frag_corr_top1_top2[i]      = c_top1_top2
            psms.frag_corr_top1_top3[i]      = c_top1_top3
            psms.frag_corr_top1_weight[i]    = c_top1_weight
            psms.frag_corr_mean_pairwise[i]  = mean_pairwise
            psms.frag_corr_min_pairwise[i]   = min_pairwise
            psms.frag_corr_top3_weight[i]    = c_top3_w
            psms.frag_corr_top5_weight[i]    = c_top5_w
            psms.frag_apex_dispersion_irt[i] = apex_disp
            psms.n_correlated_fragments[i]     = n_corr_70
            psms.n_correlated_fragments_50[i]  = n_corr_50
            psms.n_correlated_fragments_90[i]  = n_corr_90
            psms.frag_corr_best_weight[i]      = c_best_w
            psms.frag_corr_best_m0[i]          = c_best_m0
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
    n == 0 && return
    @assert hasproperty(psms, :irt_obs) ":irt_obs required"
    by_prec = Dict{eltype(psms.precursor_idx), Vector{Int}}()
    @inbounds for i in 1:n
        push!(get!(() -> Int[], by_prec, psms.precursor_idx[i]), i)
    end
    @inbounds for (_, idxs) in by_prec
        # Find apex (max weight)
        apex_i = idxs[1]; apex_w = Float32(psms.weight[apex_i])
        for j in idxs
            w = Float32(psms.weight[j])
            if w > apex_w; apex_w = w; apex_i = j; end
        end
        apex_irt = Float32(psms.irt_obs[apex_i])
        for i in idxs
            psms.irt_dist_to_weight_apex[i] = abs(Float32(psms.irt_obs[i]) - apex_irt)
        end
    end
    return
end

"""
    add_neighborhood_features!(psms::DataFrame)

For each PSM at (precursor_idx, scan_idx), look at the PSMs of the SAME
precursor within ±k scans (sorted by scan_idx). For each of:
- gof (higher = better, take max)
- fitted_manhattan_distance (lower = better, take min)
- max_matched_residual (lower = better, take min)

compute the BEST value in the window and the |iRT| distance from this PSM to
the scan that contributed the best value (a chromatographic-quality signal:
real precursors should have their best neighbor close in iRT).

Currently emits TWO window sizes for A/B comparison via LightGBM gain:
- 3-scan (`half_window=1`): self + 1 before + 1 after — the original window.
- 5-scan (`half_window=2`): self + 2 before + 2 after — wider chromatographic
  context, includes flanking shoulders of an elution peak.

Adds columns (12 total):
- `:best_gof_3scan`, `:best_manhattan_3scan`, `:best_max_residual_3scan`,
  `:irt_dist_best_gof_3scan`, `:irt_dist_best_manhattan_3scan`,
  `:irt_dist_best_max_residual_3scan`,
- `:best_gof_5scan`, `:best_manhattan_5scan`, `:best_max_residual_5scan`,
  `:irt_dist_best_gof_5scan`, `:irt_dist_best_manhattan_5scan`,
  `:irt_dist_best_max_residual_5scan`.

PSMs that are the only PSM for a precursor get values from themselves; iRT
distance = 0.
"""
function add_neighborhood_features!(psms::DataFrame)
    # `best_*` + `irt_dist_best_*` emitted only for the narrow 3-scan + 5-scan
    # windows (the strict elution-shoulder context). Window sweep showed wider
    # best_gof_*scan variants dominated individual feature importance, but
    # adding 7/9/11/15 didn't translate into a corresponding per-file q≤.01 lift
    # — within LGBM run-to-run noise. Trim to 3+5 to keep the feature count
    # lean.
    #
    # `worst_*` features stay at 3-scan + 11-scan (chosen because the 11-scan
    # worst_max_residual was the dominant feature when first tried). worst is
    # restricted to max_residual + manhattan; worst_gof was uninformative.
    _add_neighborhood_features_window!(psms, 1, "_3scan";  emit_worst=(:max_residual, :manhattan))
    _add_neighborhood_features_window!(psms, 2, "_5scan")
    _add_neighborhood_features_window!(psms, 5, "_11scan"; emit_best=false,
                                        emit_worst=(:max_residual, :manhattan))
    return
end

"""
    _add_neighborhood_features_window!(psms, half_window, suffix)

Worker for `add_neighborhood_features!`. Writes columns
`best_gof<suffix>`, `best_manhattan<suffix>`, `best_max_residual<suffix>`,
and their `irt_dist_*<suffix>` companions, computed over the window
`-half_window .. +half_window` around each PSM's scan within its precursor.
"""
function _add_neighborhood_features_window!(psms::DataFrame, half_window::Int, suffix::String;
                                             emit_best::Bool = true,
                                             emit_worst::NTuple = ())
    n = nrow(psms)
    if emit_best
        gof_col      = Symbol("best_gof", suffix)
        man_col      = Symbol("best_manhattan", suffix)
        mr_col       = Symbol("best_max_residual", suffix)
        dgof_col     = Symbol("irt_dist_best_gof", suffix)
        dman_col     = Symbol("irt_dist_best_manhattan", suffix)
        dmr_col      = Symbol("irt_dist_best_max_residual", suffix)
        psms[!, gof_col]  = Vector{Float32}(undef, n)
        psms[!, man_col]  = Vector{Float32}(undef, n)
        psms[!, mr_col]   = Vector{Float32}(undef, n)
        psms[!, dgof_col] = zeros(Float32, n)
        psms[!, dman_col] = zeros(Float32, n)
        psms[!, dmr_col]  = zeros(Float32, n)
    end
    # Worst-in-window features are computed only for the metrics specified by
    # `emit_worst` (subset of `:gof`, `:manhattan`, `:max_residual`).
    want_worst_gof = :gof in emit_worst
    want_worst_man = :manhattan in emit_worst
    want_worst_mr  = :max_residual in emit_worst
    worst_gof_col = want_worst_gof ? Symbol("worst_gof",          suffix) : nothing
    worst_man_col = want_worst_man ? Symbol("worst_manhattan",    suffix) : nothing
    worst_mr_col  = want_worst_mr  ? Symbol("worst_max_residual", suffix) : nothing
    if want_worst_gof; psms[!, worst_gof_col] = Vector{Float32}(undef, n); end
    if want_worst_man; psms[!, worst_man_col] = Vector{Float32}(undef, n); end
    if want_worst_mr;  psms[!, worst_mr_col]  = Vector{Float32}(undef, n); end
    n == 0 && return
    @assert hasproperty(psms, :irt_obs) "neighborhood features require :irt_obs column"

    # Group rows by precursor_idx; within each group sort indices by scan_idx
    by_prec = Dict{eltype(psms.precursor_idx), Vector{Int}}()
    @inbounds for i in 1:n
        push!(get!(() -> Int[], by_prec, psms.precursor_idx[i]), i)
    end

    @inbounds for (_, idxs) in by_prec
        if length(idxs) > 1
            sort!(idxs, by = j -> psms.scan_idx[j])
        end
        for (k, i) in enumerate(idxs)
            cur_irt = Float32(psms.irt_obs[i])
            cur_gof = Float32(psms.gof[i])
            cur_man = Float32(psms.fitted_manhattan_distance[i])
            cur_mr  = Float32(psms.max_matched_residual[i])
            best_gof = cur_gof; gof_idx = i
            best_man = cur_man; man_idx = i
            best_mr  = cur_mr;  mr_idx  = i
            worst_gof = cur_gof
            worst_man = cur_man
            worst_mr  = cur_mr
            for off in -half_window:half_window
                off == 0 && continue
                kk = k + off
                if 1 <= kk <= length(idxs)
                    j = idxs[kk]
                    g = Float32(psms.gof[j])
                    if g > best_gof; best_gof = g; gof_idx = j; end
                    if want_worst_gof && g < worst_gof; worst_gof = g; end
                    m = Float32(psms.fitted_manhattan_distance[j])
                    if m < best_man; best_man = m; man_idx = j; end
                    if want_worst_man && m > worst_man; worst_man = m; end
                    r = Float32(psms.max_matched_residual[j])
                    if r < best_mr;  best_mr  = r;  mr_idx  = j; end
                    if want_worst_mr && r > worst_mr; worst_mr = r; end
                end
            end
            if emit_best
                psms[!, gof_col][i]  = best_gof
                psms[!, man_col][i]  = best_man
                psms[!, mr_col][i]   = best_mr
                psms[!, dgof_col][i] = abs(Float32(psms.irt_obs[gof_idx]) - cur_irt)
                psms[!, dman_col][i] = abs(Float32(psms.irt_obs[man_idx]) - cur_irt)
                psms[!, dmr_col][i]  = abs(Float32(psms.irt_obs[mr_idx])  - cur_irt)
            end
            if want_worst_gof; psms[!, worst_gof_col][i] = worst_gof; end
            if want_worst_man; psms[!, worst_man_col][i] = worst_man; end
            if want_worst_mr;  psms[!, worst_mr_col][i]  = worst_mr;  end
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
    :irt_error, :irt_diff,
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
4. add_features! — irt_error, irt_diff, tic, prec_mz, sequence_length, etc.
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

    # 4. Add ML features (irt_error, irt_diff, tic, prec_mz, sequence_length, etc.)
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
