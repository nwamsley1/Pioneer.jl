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

    # Lazy-populate Mox cache: Vector-backed for O(1) lookup, computed on first access per precursor
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
    :n_scans,
]

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
