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

#==========================================================
QuadTuning fused-path dispatches.

Defined here (not in process_scans.jl) because QuadTuningSearchParameters
is loaded after process_scans.jl by the SearchMethods dir walker. These
mirror the ParameterTuningSearchParameters dispatches — Complex PSM
output, OLS-or-PMM deconv per params, MainSearch-style distance metrics.
==========================================================#

get_scored_psms(sd::SearchDataStructures, ::QuadTuningSearchParameters) = getTuningScoredPsms(sd)
get_unscored_psms(sd::SearchDataStructures, ::QuadTuningSearchParameters) = getTuningUnscoredPsms(sd)

function resize_if_needed!(search_data::SearchDataStructures, params::QuadTuningSearchParameters)
    weights = getTempWeights(search_data)
    if n_active(getIdToCol(search_data)) > length(weights)
        resize_arrays!(search_data, weights)
    end
end

function post_design_matrix!(search_data::SearchDataStructures,
                              Hs::AbstractSparseDesignMatrix,
                              params::QuadTuningSearchParameters)
    weights = getTempWeights(search_data)
    initialize_weights!(getIdToCol(search_data), weights, getPrecursorWeights(search_data))
    converged = first(solve_deconvolution!(
        params.deconvolution_solver,
        Hs, getResiduals(search_data), weights, getColNorm2(search_data),
        getMu(search_data), getObserved(search_data),
        params.max_iter_outer, params.max_diff))
    if converged
        update_precursor_weights!(getIdToCol(search_data), weights, getPrecursorWeights(search_data))
        zero_negligible_weights!(weights, Hs.n)
    end
    return converged
end

function compute_distance_metrics!(Hs::AbstractSparseDesignMatrix,
                                    search_data::SearchDataStructures,
                                    params::QuadTuningSearchParameters;
                                    lod_intensity_threshold::Float32 = 0.0f0)
    getDistanceMetrics(getTempWeights(search_data), getResiduals(search_data),
        Hs, getMainSearchSpectralScores(search_data);
        lod_intensity_threshold=lod_intensity_threshold)
end

function score_psms!(
    search_data::SearchDataStructures,
    params::QuadTuningSearchParameters,
    Hs::AbstractSparseDesignMatrix,
    scan_idx::Int64,
    nmatches::Int64,
    nmisses::Int64,
    spectra::MassSpecData,
    last_val::Int64,
    ms_file_idx::Int64,
    cycle_idx::Int64;
    mem::AbstractMassErrorModel = SimpleMassErrorModel(0f0, (0f0, 0f0))
)
    score_result = Score!(
        getTuningScoredPsms(search_data),
        getTuningUnscoredPsms(search_data),
        getMainSearchSpectralScores(search_data),
        getTempWeights(search_data),
        getIdToCol(search_data),
        ms_file_idx,
        cycle_idx,
        nmatches / (nmatches + nmisses),
        last_val,
        Hs.n,
        Float32(sum(getIntensityArray(spectra, scan_idx))),
        scan_idx;
        block_size = 500000,
        default_top3_ll = get_default_top3_ll(mem))
    return score_result.last_val
end

"""
    check_window_widths(spectra::MassSpecData) -> Set{String}

Examine MS2 isolation window widths across all spectra.

# Arguments
- `spectra`: Table containing MS data with isolationWidthMz and msOrder columns

# Returns
Set of unique isolation window widths (as strings) for MS2 spectra.
Used to verify consistent window settings across runs.
"""
function check_window_widths(spectra::MassSpecData)
    first_width = 0.0
    all_same = true
    for i in 1:length(spectra)
        if getMsOrder(spectra, i) == 2 && !ismissing(getIsolationWidthMz(spectra, i))
            w = round(Float64(getIsolationWidthMz(spectra, i)), digits=1)
            if first_width == 0.0
                first_width = w
            elseif w != first_width
                all_same = false
                break
            end
        end
    end
    if first_width == 0.0
        return Float64[]
    end
    return all_same ? [first_width] : unique([round(Float64(getIsolationWidthMz(spectra, i)), digits=1)
        for i in 1:length(spectra) if getMsOrder(spectra, i) == 2 && !ismissing(getIsolationWidthMz(spectra, i))])
end

#==========================================================
PSM scoring and processing
==========================================================#
function get_nearest_adjacent_scans(scan_idx::UInt32,
                            centerMz::AbstractArray{Union{Missing, T}},
                            isolationWidthMz::AbstractArray{Union{Missing, T}};
                            scans_to_search::Int64 = 500
        ) where {T<:AbstractFloat}
    
    upperBoundMz = centerMz[scan_idx] + isolationWidthMz[scan_idx]/T(2.0)
    min_diff, min_diff_idx = typemax(Float32), -1
    for near_scan_idx in range(scan_idx, min(scan_idx + scans_to_search, length(centerMz)))
        if ismissing(centerMz[near_scan_idx])
            continue
        end
        lowerBoundMz = centerMz[near_scan_idx] - isolationWidthMz[near_scan_idx]/T(2.0)
        if max(abs(upperBoundMz - lowerBoundMz), 0.1) < min_diff
            min_diff_idx = near_scan_idx
            min_diff = abs(upperBoundMz - lowerBoundMz) 
        end
    end
    next_scan_idx = sign(min_diff_idx)==-1 ? scan_idx : min_diff_idx

    min_diff, min_diff_idx = typemax(Float32), -1
    lowerBoundMz = centerMz[scan_idx] - isolationWidthMz[scan_idx]/T(2.0)
    for near_scan_idx in (scan_idx):-1:max(scan_idx - scans_to_search, 1)
        if ismissing(centerMz[near_scan_idx])
            continue
        end
        upperBoundMz = centerMz[near_scan_idx] + isolationWidthMz[near_scan_idx]/T(2.0)
        if max(abs(upperBoundMz - lowerBoundMz), 0.1) < min_diff
            min_diff_idx = near_scan_idx
            min_diff = abs(upperBoundMz - lowerBoundMz) 
        end
    end
    prev_scan_idx = sign(min_diff_idx)==-1 ? scan_idx : min_diff_idx



    return prev_scan_idx, next_scan_idx
end

function get_scan_to_prec_idx(
    scan_idxs::AbstractVector{UInt32},
    prec_idxs::AbstractVector{UInt32},
    centerMz::AbstractVector{Union{Missing, Float32}},
    isolationWidthMz::AbstractVector{Union{Missing, Float32}}
    )
    N = length(scan_idxs)
    scan_idx_to_prec_idx = Dictionary{UInt32, Vector{UInt32}}()
    adjacent_cache = Dict{UInt32, Tuple{Int, Int}}()

    for i in range(1, N)
        scan_idx = scan_idxs[i]
        prec_idx = prec_idxs[i]

        # Get or compute adjacent scans (cached per scan_idx)
        prev_scan_idx, next_scan_idx = get!(adjacent_cache, scan_idx) do
            get_nearest_adjacent_scans(scan_idx, centerMz, isolationWidthMz)
        end

        # Add precursor to its own scan and both adjacent scans
        for sid in (scan_idx, UInt32(prev_scan_idx), UInt32(next_scan_idx))
            if haskey(scan_idx_to_prec_idx, sid)
                push!(scan_idx_to_prec_idx[sid], prec_idx)
            else
                insert!(scan_idx_to_prec_idx, sid, [prec_idx])
            end
        end
    end
    return scan_idx_to_prec_idx
end

function add_columns!(
    precursor_idx::AbstractVector{UInt32},
    lib_precursor_mz::AbstractVector{Float32},
    lib_prec_charge::AbstractVector{UInt8},
    lib_sulfur_count::AbstractVector{UInt8},
    iso_idx::AbstractVector{UInt8},
    center_mz::AbstractVector{Float32},
    iso_splines::IsotopeSplineModel)

    mono_mz = [lib_precursor_mz[pid] for pid in precursor_idx]
    prec_charge = [lib_prec_charge[pid] for pid in precursor_idx]
    sulfur_count = [lib_sulfur_count[pid] for pid in precursor_idx]
    #The iso_idx is 1 indexed. So the M0 has an index of 1, the M+1 had an index of 2, etc.
    iso_mz = Float32.(mono_mz .+ C13_C12_MASS_DIFF.*(iso_idx.-1.0f0)./prec_charge)
    mz_offset = iso_mz .- center_mz
    δ = zeros(Float32, length(precursor_idx))
    for i in range(1, length(precursor_idx))
        s_count = min(Int64(sulfur_count[i]), 5)
        mono_mass = Float32((mono_mz[i] - PROTON)*prec_charge[i])
        δ[i] = iso_splines(s_count, 0, mono_mass)/iso_splines(s_count, 1, mono_mass)
    end
    return DataFrame((mono_mz = mono_mz, prec_charge = prec_charge, sulfur_count = sulfur_count, iso_mz = iso_mz, mz_offset = mz_offset, δ = δ))
end

"""
Filter quad PSMs based on criteria:
- M0/M1 isotopes only
- Minimum number of matches
- Non-zero abundance
"""
function filter_quad_psms(
    iso_idx::AbstractVector{UInt8},
    n_matches::AbstractVector{UInt8},
    weight::AbstractVector{Float32},
    charge::AbstractVector{UInt8},
    params::QuadTuningSearchParameters
)
    n = length(iso_idx)
    mask = Vector{Bool}(undef, n)
    @inbounds for i in 1:n
        mask[i] = (
            (iso_idx[i] < 3) &&  # M0/M1 only
            (n_matches[i] >= params.min_quad_tuning_fragments) &&  # Min matches
            (weight[i] > 0)  &&
            (charge[i] == 2)# Non-zero abundance
        )
    end
    return mask
end

function getCharges(prec_charges::AbstractVector{UInt8}, precursor_idx::AbstractVector{UInt32})
    charges = zeros(UInt8, length(precursor_idx))
    for i in range(1, length(precursor_idx))
        charges[i] = prec_charges[precursor_idx[i]]
    end
    return charges
end

#==========================================================
Quad Scan Priority Index (moved from separate file to ensure loading)
==========================================================#

"""
    QuadScanPriorityIndex

Contains prioritized scan indices for memory-efficient quadrupole tuning.
Only stores metadata, never loads actual scan data until needed.
"""
struct QuadScanPriorityIndex
    scan_indices::Vector{Int32}        # Ordered scan indices by priority
    center_mz_values::Vector{Float32}  # Center m/z values for each scan
    tic_values::Vector{Float32}        # TIC values for each scan
    ms_orders::Vector{UInt8}           # MS order for verification
    total_ms2_count::Int32             # Total MS2 scans
    n_mz_bins::Int32                   # Number of m/z bins used
end

"""
    build_quad_scan_priority_index(spectra::MassSpecData;
                                   n_mz_bins::Int = 20,
                                   target_ms_order::UInt8 = UInt8(2))::QuadScanPriorityIndex

Build prioritized scan index for quadrupole tuning without loading scan data.
"""
function build_quad_scan_priority_index(
    spectra::MassSpecData;
    n_mz_bins::Int = 20,
    target_ms_order::UInt8 = UInt8(2)
)::QuadScanPriorityIndex

    # Step 1: Extract metadata (NO peak data loaded)
    center_mz_values = getCenterMzs(spectra)
    tic_values = getTICs(spectra)
    ms_orders = getMsOrders(spectra)

    # Step 2: Filter for MS2 scans

    ms2_mask = ms_orders .== target_ms_order
    ms2_indices = findall(ms2_mask)
    n_ms2 = length(ms2_indices)

    if n_ms2 == 0
        return QuadScanPriorityIndex(
            Int32[], Float32[], Float32[], UInt8[],
            Int32(0), Int32(n_mz_bins)
        )
    end

    # Step 3: Create center m/z bins

    ms2_center_mz = center_mz_values[ms2_indices]
    ms2_tic = tic_values[ms2_indices]

    # Filter out missing values
    valid_mask = .!ismissing.(ms2_center_mz)
    if !all(valid_mask)
        ms2_indices = ms2_indices[valid_mask]
        ms2_center_mz = ms2_center_mz[valid_mask]
        ms2_tic = ms2_tic[valid_mask]
        n_ms2 = length(ms2_indices)
    end

    if n_ms2 == 0
        return QuadScanPriorityIndex(
            Int32[], Float32[], Float32[], UInt8[],
            Int32(0), Int32(n_mz_bins)
        )
    end

    mz_min, mz_max = extrema(ms2_center_mz)

    # Handle edge case of single m/z value
    if mz_min == mz_max
        priority_order = Int32.(ms2_indices[sortperm(ms2_tic, rev=true)])
        return QuadScanPriorityIndex(
            priority_order,
            center_mz_values[priority_order],
            tic_values[priority_order],
            ms_orders[priority_order],
            Int32(n_ms2),
            Int32(1)
        )
    end

    bin_width = (mz_max - mz_min) / n_mz_bins

    # Assign each MS2 scan to a bin
    bin_assignments = zeros(Int, n_ms2)
    for (i, center_mz) in enumerate(ms2_center_mz)
        bin_idx = min(max(1, ceil(Int, (center_mz - mz_min) / bin_width)), n_mz_bins)
        bin_assignments[i] = bin_idx
    end

    # Step 4: Sort within bins by TIC and create priority order

    # Group scans by bin
    bins = [Int32[] for _ in 1:n_mz_bins]
    for (i, bin_idx) in enumerate(bin_assignments)
        push!(bins[bin_idx], ms2_indices[i])
    end

    # Sort each bin by TIC (descending - highest TIC first)
    for bin in bins
        if !isempty(bin)
            sort!(bin, by=idx -> tic_values[idx], rev=true)
        end
    end

    # Round-robin selection from bins
    priority_order = Int32[]
    max_bin_size = maximum(length.(bins))

    for round in 1:max_bin_size
        for bin_idx in 1:n_mz_bins
            if round <= length(bins[bin_idx])
                push!(priority_order, bins[bin_idx][round])
            end
        end
    end


    return QuadScanPriorityIndex(
        priority_order,
        center_mz_values[priority_order],
        tic_values[priority_order],
        ms_orders[priority_order],
        Int32(n_ms2),
        Int32(n_mz_bins)
    )
end

"""
    progressive_quad_psm_collection!(scan_index::QuadScanPriorityIndex,
                                    spectra::MassSpecData,
                                    search_context::SearchContext,
                                    params::QuadTuningSearchParameters,
                                    ms_file_idx::Int64;
                                    initial_percent::Float64 = 10.0,
                                    min_psms_required::Int = 1000,
                                    verbose::Bool = true)

Progressively collect PSMs using sampling without replacement for quadrupole tuning.
"""
function progressive_quad_psm_collection!(
    scan_index::QuadScanPriorityIndex,
    spectra::MassSpecData,
    search_context::SearchContext,
    params::QuadTuningSearchParameters,
    ms_file_idx::Int64;
    target_psms::Int,
    fallback_min_psms::Int = target_psms,
    score_tiers::Tuple = TUNING_SCORE_TIERS,
    n_required_top::Int = TUNING_N_REQUIRED_TOP,
    fdr_threshold::Float64 = 0.01
)
    all_scan_indices = scan_index.scan_indices
    max_scans = length(all_scan_indices)
    fdr_scale = getLibraryFdrScaleFactor(search_context)
    precursors = getPrecursors(getSpecLib(search_context))

    if max_scans == 0
        return DataFrame(), false, 0
    end

    scored_psms = DataFrame()
    n_passing = 0
    total_scans_used = 0

    for (tier_idx, score) in enumerate(score_tiers)
        lut = make_top_n_required_lut(n_required_top, Int(score))
        setBitVecFilter!(search_context, ms_file_idx, lut)

        raw_psms = DataFrame()
        prev_scan_end = 0
        n_passing = 0

        # First tier: progressive scan schedule (start small, double until 80%,
        # then finish at 100%). Subsequent tiers: all scans immediately.
        scan_schedule = tier_idx == 1 ? [0.025, 0.05, 0.10, 0.20, 0.40, 0.80, 1.00] : [1.00]

        for frac in scan_schedule
            scan_end = min(ceil(Int, max_scans * frac), max_scans)
            scan_end <= prev_scan_end && continue

            chunk_indices = all_scan_indices[(prev_scan_end+1):scan_end]
            chunk_psms = library_search(spectra, search_context, params, ms_file_idx;
                                        scan_indices = Int.(chunk_indices))

            if !isempty(chunk_psms)
                add_tuning_search_columns!(
                    chunk_psms, spectra,
                    getIsDecoy(precursors), getIrt(precursors),
                    getCharge(precursors), getRetentionTimes(spectra),
                    getTICs(spectra))
                if isempty(raw_psms)
                    raw_psms = chunk_psms
                else
                    append!(raw_psms, chunk_psms)
                end
            end
            prev_scan_end = scan_end

            n_raw = nrow(raw_psms)
            if n_raw >= 50
                scored_tmp = copy(raw_psms)
                sort!(scored_tmp, [:rt, :precursor_idx])
                score_presearch!(scored_tmp)
                scored_tmp[!, :q_value] = zeros(Float16, nrow(scored_tmp))
                get_qvalues!(scored_tmp[!,:prob], scored_tmp[!,:target],
                             scored_tmp[!,:q_value]; fdr_scale_factor=fdr_scale)
                n_passing = count(row -> row.q_value::Float16 <= Float16(fdr_threshold) && row.target::Bool, eachrow(scored_tmp))
                scored_psms = scored_tmp
            end

            n_passing >= target_psms && break
        end

        total_scans_used = prev_scan_end
        n_passing >= target_psms && break
    end

    delete!(search_context.bitvec_filter, ms_file_idx)

    # Always filter to 1% FDR + target + best-per-precursor + charge-2 so downstream
    # code (including fallback-path plotting) gets clean PSMs. `converged` reflects
    # whether the count is sufficient; the caller may still use the PSMs otherwise.
    if !isempty(scored_psms)
        filter!(row -> row.q_value::Float16 <= Float16(fdr_threshold) && row.target::Bool, scored_psms)
        if !isempty(scored_psms)
            sort!(scored_psms, :prob, rev=true)
            scored_psms = combine(groupby(scored_psms, :precursor_idx), first)
            charges = getCharges(getCharge(precursors), UInt32.(scored_psms[!, :precursor_idx]))
            scored_psms = scored_psms[charges .== UInt8(2), :]
        end
    end

    converged = n_passing >= fallback_min_psms
    return scored_psms, converged, total_scans_used
end



"""
    summarize_precursor(iso_idx::AbstractVector{UInt8}, center_mz::AbstractVector{Float32},
                      iso_mz::AbstractVector{Float32}, prec_charge::AbstractVector{UInt8},
                      weight::AbstractVector{Float32}, δ::AbstractVector{Float32}) 
                      -> NamedTuple

Summarize isotope pair measurements for quadrupole transmission analysis.

# Arguments
- `iso_idx`: Indices of isotope peaks (1=M0, 2=M1)
- `center_mz`: Center m/z values of isolation windows
- `iso_mz`: m/z values of isotope peaks
- `prec_charge`: Precursor charge states
- `weight`: Fitted weights for each isotope peak
- `δ`: Theoretical isotope ratios

# Returns
NamedTuple containing:
- `center_mz`: Center m/z of isolation window
- `δ`: Theoretical isotope ratio
- `yt`: Log ratio of observed to theoretical isotope ratio
- `x0`: m/z offset of M0 peak from window center
- `x1`: m/z offset of M1 peak from window center
- `prec_charge`: Precursor charge state

Returns tuple of missing values if:
- More than 2 isotope peaks present
- Isotope peaks aren't M0 and M1
- Invalid isotope pattern

Note: Only processes M0/M1 isotope pairs for transmission modeling.
"""
function summarize_precursor(
    iso_idx::AbstractVector{UInt8},
    center_mz::AbstractVector{Float32},
    iso_mz::AbstractVector{Float32},
    prec_charge::AbstractVector{UInt8},
    weight::AbstractVector{Float32},
    δ::AbstractVector{Float32})
    
    # Handle empty groups gracefully
    if isempty(iso_idx)
        return (center_mz = missing, 
                δ = missing, 
                yt = missing, 
                x0 = missing, 
                x1 = missing, 
                prec_charge = missing)
    end
    
    if (length(iso_idx) == 2)
        if ((iso_idx[1] == 1) && (iso_idx[2] == 2))
            m0_idx, m1_idx = 0, 0
            if iso_idx[1] == 1
                m0_idx, m1_idx = 1, 2
            else
                m0_idx, m1_idx = 2, 1
            end
            return (center_mz = center_mz[m0_idx],
                    δ = δ[m0_idx],
                    yt = log(weight[m0_idx]/(weight[m1_idx]*δ[m0_idx])), 
                    x0 = iso_mz[m0_idx]-center_mz[m0_idx], 
                    x1 = iso_mz[m1_idx]-center_mz[m1_idx], 
                    prec_charge = prec_charge[m0_idx])
        end
    end
    
    #If we only got the M0
    if (length(iso_idx) == 1) && (iso_idx[1] == 1)
        m0_idx = 1
        m1_mz = iso_mz[m0_idx] + (C13_C12_MASS_DIFF/prec_charge[m0_idx])

        return (center_mz = center_mz[m0_idx], 
            δ = δ[m0_idx], 
            yt = 10.0f0, # assume this is the most extreme ratio we could see
            x0 = iso_mz[m0_idx]-center_mz[m0_idx],
            x1 = m1_mz - center_mz[m0_idx],
            prec_charge = prec_charge[m0_idx]) 
    end

    #If we only got the M1
    if (length(iso_idx) == 1) && (iso_idx[1] == 2)
        m1_idx = 1
        mo_mz = iso_mz[m1_idx] - (C13_C12_MASS_DIFF/prec_charge[m1_idx])

        return (center_mz = center_mz[m1_idx], 
               δ = δ[m1_idx], 
               yt = -10.0f0, # assume this is the most extreme ratio we could see
               x0 = mo_mz - center_mz[m1_idx],
               x1 = iso_mz[m1_idx]-center_mz[m1_idx],
               prec_charge = prec_charge[m1_idx])
    end

    return (center_mz = missing, 
            δ = missing, 
            yt = missing, 
            x0 = missing, 
            x1 = missing, 
            prec_charge = missing)
end

"""
    process_quad_results(psms::DataFrame, precursors::LibraryPrecursors,
                        iso_splines::IsotopeSplineModel) -> DataFrame

Process quadrupole transmission search results.

# Arguments
- `psms`: PSMs from quad transmission search
- `precursors`: Library precursor information
- `iso_splines`: Isotope spline models

# Process
1. Sorts PSMs by scan, precursor, and isotope index
2. Adds columns for precursor properties
3. Combines results by scan and precursor
4. Performs post-processing

# Returns
Processed DataFrame ready for quad model fitting.
"""
function process_quad_results(
    psms::DataFrame,
    precursors::LibraryPrecursors,
    iso_splines::IsotopeSplineModel
)
    sort!(psms, [:scan_idx, :precursor_idx, :iso_idx])
    
    processed = hcat(psms, add_columns!(
        psms[!, :precursor_idx],
        getMz(precursors),
        getCharge(precursors),#[:prec_charge],
        getSulfurCount(precursors),#[:sulfur_count],
        psms[!, :iso_idx],
        psms[!, :center_mz],
        iso_splines
    ))

    sort!(processed, [:scan_idx, :precursor_idx, :iso_idx])
    combined = combine(groupby(processed, [:scan_idx, :precursor_idx])) do group
        summarize_precursor(
            group[!, :iso_idx],
            group[!, :center_mz],
            group[!, :iso_mz],
            group[!, :prec_charge],
            group[!, :weight],
            group[!, :δ]
        )
    end
    
    # Filter out rows where summarize_precursor returned all missing values
    # Since summarize_precursor returns either all values or all missing,
    # we only need to check one column
    filter!(row -> !ismissing(row.center_mz), combined)
    
    postprocess_combined_results!(combined)
    return combined
end

"""
    postprocess_combined_results!(combined::DataFrame) -> DataFrame

Additional processing of combined quad search results.

# Arguments
- `combined`: DataFrame of combined results

# Actions
1. Removes rows with missing transmission values
2. Converts column types to appropriate numeric types
3. Returns processed DataFrame
"""
function postprocess_combined_results!(combined::DataFrame)
    filter!(row -> !ismissing(row.yt), combined)
    combined[!, :prec_charge] = UInt8.(combined[!, :prec_charge])
    combined[!, :x0] = Float32.(combined[!, :x0])
    combined[!, :yt] = Float32.(combined[!, :yt])
    combined[!, :x1] = Float32.(combined[!, :x1])
    #filter!(row -> row.prec_charge < 3, combined)
    return combined
end

#==========================================================
Quad Transmission Search
=========================================================#
"""
    perform_quad_transmission_search(spectra::MassSpecData,
                                  results::QuadTuningSearchResults,
                                  scan_idx_to_prec_idx::Dictionary{UInt32, Vector{UInt32}},
                                  search_context::SearchContext,
                                  params::QuadTuningSearchParameters,
                                  ms_file_idx::Int64) -> DataFrame

Execute quadrupole transmission search across MS data.

# Arguments
- `spectra`: MS/MS spectral data
- `results`: Current quad tuning results
- `scan_idx_to_prec_idx`: Mapping of scan indices to precursor indices
- `search_context`: Search context
- `params`: Search parameters
- `ms_file_idx`: MS file index

# Process
1. Partitions work across threads
2. For each scan:
   - Selects transitions for quad estimation
   - Matches peaks
   - Performs deconvolution
   - Records results
3. Combines results across threads

# Returns
DataFrame containing quad transmission search results.
"""
function perform_quad_transmission_search(
    spectra::MassSpecData,
    results::QuadTuningSearchResults,
    scan_idx_to_prec_idx::Dictionary{UInt32, Vector{UInt32}},
    search_context::SearchContext,
    params::QuadTuningSearchParameters,
    ms_file_idx::Int64
)
    """
    Record results for a single scan.
    """
    function record_scan_results!(
        search_data::SearchDataStructures,
        tuning_results::Vector{@NamedTuple{
            precursor_idx::UInt32,
            scan_idx::UInt32,
            weight::Float32,
            iso_idx::UInt8,
            center_mz::Float32,
            n_matches::UInt8
        }},
        weights::Vector{Float32},
        Hs::AbstractSparseDesignMatrix,
        scan_idx::Int,
        center_mz::Float32
    )
        precursor_weights = getPrecursorWeights(search_data)
        for (id, colid) in active_keys(getIdToCol(search_data))
            # Update precursor weights
            precursor_weights[id] = weights[colid]

            # Decode artificial column key: id = (pid - 1) * 3 + iso_pass
            isotope_idx = UInt8(((id - 1) % 3) + 1)
            pid = UInt32(((id - 1) ÷ 3) + 1)

            # Count matches via the abstract matched_at accessor — works
            # for both SparseArray and SparseArrayFused layouts.
            n_matches = 0
            @inbounds for j in Hs.colptr[colid]:(Hs.colptr[colid+1] - 1)
                n_matches += matched_at(Hs, j) ? 1 : 0
            end

            # Record result
            push!(tuning_results, (
                precursor_idx = pid,
                scan_idx = scan_idx,
                weight = weights[colid],
                iso_idx = isotope_idx,
                center_mz = center_mz,
                n_matches = n_matches
            ))
        end
        #Arrow.write("/Users/n.t.wamsley/Desktop/test.arrow", DataFrame(tuning_results))
        #ttable = DataFrame(Tables.columntable(Arrow.Table("/Users/n.t.wamsley/Desktop/test.arrow")))

    end

    """
    Process a single scan for quad transmission search (fused kernel).

    Uses `run_fused!(FusedQuadEst(), ...)` to produce `Hs_fused` directly
    in CSC order. No `selectTransitions!` / `matchPeaks!` / `sort!` /
    `buildDesignMatrix!` / `sortSparse!` / `ScoreFragmentMatches!` — all
    replaced by the single fused pass.
    """
    function process_scan!(
        scan_idx::Int,
        scan_idxs::Set{UInt32},
        spectra::MassSpecData,
        tuning_results::Vector{@NamedTuple{
            precursor_idx::UInt32,
            scan_idx::UInt32,
            weight::Float32,
            iso_idx::UInt8,
            center_mz::Float32,
            n_matches::UInt8
        }},
        search_context::SearchContext,
        search_data::SearchDataStructures,
        params::QuadTuningSearchParameters,
        ms_file_idx::Int64,
        Hs::SparseArrayFused{UInt32, Float32},
        weights::Vector{Float32},
        precursor_weights::AbstractPrecursorMap{Float32},
        residuals::Vector{Float32},
        scan_idx_to_prec_idx::Dictionary{UInt32, Vector{UInt32}}
    )
        scan_idx ∉ scan_idxs && return

        msn = getMsOrder(spectra, scan_idx)
        msn ∉ params.spec_order && return

        nce_model = getNceModel(search_context, ms_file_idx)
        mem       = getMassErrorModel(search_context, ms_file_idx)
        spec_lib  = getSpecLib(search_context)
        ion_list  = getFragmentLookupTable(spec_lib)
        precursors = getPrecursors(spec_lib)

        prec_mzs     = getMz(precursors)
        prec_charges = getCharge(precursors)
        prec_sulfs   = getSulfurCount(precursors)
        prec_irts    = getIrt(precursors)

        # Thread-local SoA peak table + scratch.
        corr_mz        = getScanCorrectedMz(search_data)
        obs_low        = getScanObsLow(search_data)
        obs_high       = getScanObsHigh(search_data)
        isotopes_buf   = getIsotopes(search_data)
        prec_trans_buf = getPrecursorTransmission(search_data)
        fused_scratch  = getFusedScratch(search_data)
        id_to_col      = getIdToCol(search_data)

        scan_mz  = getMzArray(spectra, scan_idx)
        scan_int = getIntensityArray(spectra, scan_idx)
        scan_rt  = Float32(getRetentionTime(spectra, scan_idx))
        peak_mz_len = prepare_scan_peaks!(corr_mz, obs_low, obs_high,
                                          mem, scan_mz, scan_int, scan_rt)

        quad_fn = getQuadTransmissionFunction(
            getQuadTransmissionModel(search_context, ms_file_idx),
            getCenterMz(spectra, scan_idx),
            getIsolationWidthMz(spectra, scan_idx))

        precs_vec = scan_idx_to_prec_idx[scan_idx]

        # FusedQuadEst: 3 outer iso-passes per precursor, one-hot
        # transmission, artificial column IDs, no inline scoring.
        # iRT/isotope_err_bounds filters skipped by dispatch.
        nmatches, nmisses = run_fused!(
            FusedQuadEst(),
            Hs, getTuningUnscoredPsms(search_data), id_to_col, fused_scratch,
            corr_mz, obs_low, obs_high, peak_mz_len,
            isotopes_buf, prec_trans_buf,
            ion_list, nce_model,
            precs_vec, 1:length(precs_vec),
            prec_mzs, prec_charges, prec_sulfs, prec_irts,
            getIsoSplines(search_data), quad_fn, mem,
            scan_int, 0f0, Float32(Inf),           # scan_irt / irt_tol — skipped by kind
            (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx)),
            4,                                      # n_frag_isotopes — overridden to 0:3 via max_frag_iso_idx(FusedQuadEst)
            (UInt8(0), UInt8(0))                    # isotope_err_bounds — skipped by kind
        )

        nmatches ≤ 2 && (reset!(id_to_col); reset!(Hs); return)

        # Deconvolve: grow weight buffers to fit new columns, seed from
        # precursor_weights, run solver. (Replaces classic
        # perform_deconvolution!'s buildDesignMatrix + solve path.)
        n_active_cols = n_active(id_to_col)
        if n_active_cols > length(weights)
            new_entries = n_active_cols - length(weights) + 1000
            resize!(weights, length(weights) + new_entries)
            resize!(getColNorm2(search_data), length(getColNorm2(search_data)) + new_entries)
        end

        @inbounds for (k, col) in active_keys(id_to_col)
            weights[col] = precursor_weights[k]
        end

        solve_deconvolution!(
            params.deconvolution_solver,
            Hs,
            residuals,
            weights,
            getColNorm2(search_data),
            getMu(search_data),
            getObserved(search_data),
            params.max_iter_outer,
            params.max_diff)

        # Record results (reads Hs.colptr + matched_at(Hs, j) per column).
        record_scan_results!(
            search_data,
            tuning_results,
            weights,
            Hs,
            scan_idx,
            getCenterMz(spectra, scan_idx))

        reset!(id_to_col)
        reset!(Hs)
    end

    thread_tasks = partition_scans(spectra, Threads.nthreads())
    scan_idxs = Set(keys(scan_idx_to_prec_idx))

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            # Get thread-specific data structures
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]
            tuning_results = results.tuning_results[thread_id]
            # Get working arrays (fused kernel writes to Hs_fused).
            Hs = getHsFused(search_data)
            weights = getTempWeights(search_data)
            precursor_weights = getPrecursorWeights(search_data)
            residuals = getResiduals(search_data)

            # Process each scan
            for scan_idx in last(thread_task)
                process_scan!(
                    scan_idx,
                    scan_idxs,
                    spectra,
                    tuning_results,
                    search_context,
                    search_data,
                    params,
                    ms_file_idx,
                    Hs,
                    weights,
                    precursor_weights,
                    residuals,
                    scan_idx_to_prec_idx
                )
            end
            
            return DataFrame(tuning_results)
        end
    end
    
    return vcat(fetch.(tasks)...)
end

"""
    adjust_precursor_arrays!(search_context::SearchContext) -> SearchContext

Adjust array sizes to accommodate isotope variants.

# Arguments
- `search_context`: Search context to modify

Expands precursor arrays to handle additional isotope variants.
Returns modified search context.
"""
function adjust_precursor_arrays!(search_context::SearchContext)
    # Sparse-backed: artificial keys (prec_id - 1) * 3 + iso_pass go up to
    # ~3*n_precursors. Just empty existing maps so they start clean for
    # quad tuning; sparse storage grows on demand.
    for search_data in getSearchData(search_context)
        reset!(search_data.id_to_col)
        reset!(search_data.precursor_weights)
    end
    return search_context
end

"""
    reset_precursor_arrays!(search_context::SearchContext) -> SearchContext

Reset arrays to original size after quad tuning.

# Arguments
- `search_context`: Search context to reset

Restores precursor arrays to original size.
Returns modified search context.
"""
function reset_precursor_arrays!(search_context::SearchContext)
    # Sparse-backed: simply empty the maps so the next search method starts
    # fresh. No size-shrinking step needed; the artificial-id keys from
    # quad tuning just go away.
    for search_data in getSearchData(search_context)
        reset!(search_data.id_to_col)
        reset!(search_data.precursor_weights)
    end
    return search_context
end


#==========================================================
Plotting and Model fitting
=========================================================#
"""
    plot_charge_distributions(psms::DataFrame, results::QuadTuningSearchResults, fname::String)

Generate visualization of charge state distributions.

# Arguments
- `psms`: Processed PSMs containing charge state data
- `results`: Quad tuning results containing plot directory
- `fname`: Filename for plot

Creates scatter plot of m/z offset vs transmission for different charge states.
Saves plot to quad_plot_dir/quad_data directory.
"""

function plot_charge_distributions(psms::DataFrame, results::QuadTuningSearchResults, fname::String;
                                    quad_model::Union{Nothing, QuadTransmissionModel} = nothing,
                                    window_width::Float64 = 2.0,
                                    pad::Float64 = 0.5)
    charges = 2:2  # collection filters to charge 2 only; charge 3 adds clutter
    n_by_charge = Dict(c => count(psms[!, :prec_charge] .== c) for c in charges)
    title_str = "Quad Model Data for $fname (n=$(sum(values(n_by_charge))))"
    x_half = Float32(window_width/2 + pad)
    scatter_colors = Dict(2 => :dodgerblue)
    fit_colors     = Dict(2 => :navy)

    dense_x0 = collect(LinRange(-x_half, x_half, 400))

    # Top panel: log-ratio scatter + model F curves.
    p_top = plot(title = title_str,
                 legend = :outertopright, legendfontsize = 6,
                 xlim = (-x_half, x_half),
                 ylim = (-3.5, 3.5),
                 ylabel = L"log(\delta_i\frac{x_0}{x_1})")
    for charge in charges
        mask = psms[!, :prec_charge] .== charge
        plot!(p_top,
            psms[mask, :x0],
            # Match the clamp applied in fit_quad_model so the scatter shows the
            # same range the fits actually see.
            clamp.(psms[mask, :yt], Float32(-3), Float32(3)),
            seriestype=:scatter,
            alpha=0.1,
            markersize=2,
            color=scatter_colors[charge],
            label="Charge $charge (n=$(n_by_charge[charge]))")
        δ = Float32(1.0 / charge)
        dense_x1 = dense_x0 .+ δ
        if quad_model isa RazoQuadModel
            preds = [F(quad_model.params, dense_x0[i], dense_x1[i]) for i in eachindex(dense_x0)]
            plot!(p_top, dense_x0, preds, lw=2, color=fit_colors[charge], alpha=0.9,
                  label="Charge $charge Razo (LM)")
        end
    end

    # Bottom panel: actual transmission f(x) ∈ [0, 1] on the same x-axis.
    p_bot = plot(legend = :outertopright, legendfontsize = 6,
                 xlim = (-x_half, x_half),
                 ylim = (0, 1.05),
                 xlabel = "m/z offset of M0",
                 ylabel = "transmission")
    if quad_model isa RazoQuadModel
        tf = getQuadTransmissionFunction(quad_model, 0.0f0, 2.0f0)
        plot!(p_bot, dense_x0, [Float64(tf(Float32(x))) for x in dense_x0],
              lw=2, color=:navy, alpha=0.9, label="Razo (LM)")
    end
    return plot(p_top, p_bot, layout = grid(2, 1, heights=[0.6, 0.4]),
                size = (900, 700), link = :x)
end


"""
    plot_quad_model(quad_model::QuadTransmissionModel, window_width::Float32, results::QuadTuningSearchResults, fname::String)

Generate visualization of quadrupole transmission fit.

# Arguments
- `quad_model`: fitted quad model
- `window_width`: isolation width used to model the quad
- `results`: Quad tuning results containing plot directory
- `fname`: Filename for plot

Creates line plot of m/z offset vs transmission for the fitted model.
Saves plot to quad_plot_dir/quad_models directory.
"""

function plot_quad_model(quad_model::QuadTransmissionModel, window_width::Float64,
                          results::QuadTuningSearchResults, fname::String;
                          initial_model::Union{Nothing, RazoQuadModel} = nothing)
    padding = 2
    half_width = padding + window_width/2
    plot_bins = LinRange(-half_width, half_width, 200)

    quad_func = getQuadTransmissionFunction(quad_model, 0.0f0, 2.0f0)
    fit_label = quad_model isa RazoQuadModel ? "Razo (LM)" :
        "$(nameof(typeof(quad_model))) (fallback)"

    title_str = if quad_model isa RazoQuadModel
        p_ = quad_model.params
        "$fname\nRazo: al=$(round(p_.al,digits=2)) ar=$(round(p_.ar,digits=2)) bl=$(round(p_.bl,digits=1)) br=$(round(p_.br,digits=1))"
    else
        "$fname (fallback)"
    end
    p = plot(plot_bins, quad_func.(plot_bins), lw=2, alpha=0.85, color=:navy,
             title=title_str, titlefontsize=8, label=fit_label, xlabel="m/z offset",
             ylabel="transmission", top_margin=8*Plots.mm)
    if initial_model !== nothing
        init_func = getQuadTransmissionFunction(initial_model, 0.0f0, 2.0f0)
        plot!(p, plot_bins, init_func.(plot_bins), lw=1.5, alpha=0.7,
              color=:steelblue, label="Razo initial")
    end
    return p
end

"""
    plot_sliding_median_smoother(psms, fname; window=10)

Diagnostic plot: scatter of per-charge decimated medians (sort by x0, take
median of each non-overlapping window of `window` points). Optionally overlays
the fitted Razo quad model as a dashed curve per charge.
"""
function plot_sliding_median_smoother(psms::DataFrame, fname::String;
                                       window_width::Float64 = 2.0,
                                       bin_resolution::Float64 = 0.0625,
                                       pad::Float64 = 0.5,
                                       quad_model::Union{Nothing, QuadTransmissionModel} = nothing,
                                       initial_model::Union{Nothing, RazoQuadModel} = nothing)
    n_bins = max(8, ceil(Int, (window_width + 2 * pad) / bin_resolution))

    _fmt_params(m) = isnothing(m) ? "n/a" :
        "al=$(round(m.params.al, digits=2)) ar=$(round(m.params.ar, digits=2)) bl=$(round(m.params.bl, digits=1)) br=$(round(m.params.br, digits=1))"
    fit_str  = quad_model isa RazoQuadModel ? _fmt_params(quad_model) : "n/a"
    init_str = _fmt_params(initial_model)
    title_str = "Equal-count median (bin=$bin_resolution m/z, n_bins=$n_bins) — $fname\ninit:    $init_str\nrazo-lm: $fit_str"

    p = plot(title = title_str, titlefontsize=8,
             xlabel = "m/z offset of M0",
             ylabel = L"log(\delta_i\frac{x_0}{x_1})",
             ylim = (-3.5, 3.5),
             legend = :outertopright, legendfontsize = 6,
             size = (900, 500),
             top_margin = 10*Plots.mm)
    scatter_colors = Dict(2 => :dodgerblue)
    l2_colors      = Dict(2 => :steelblue)
    fit_colors     = Dict(2 => :navy)

    # Dense x0 grid for smooth overlay curves.
    half_plot = Float32(window_width/2 + pad)
    dense_x0 = collect(LinRange(-half_plot, half_plot, 400))

    for charge in 2:2
        mask = psms[!, :prec_charge] .== charge
        n_raw = count(mask)
        n_raw >= n_bins || continue
        charge_psms = psms[mask, :]
        sm_x0, _, sm_yt = equal_count_median(charge_psms; n_bins=n_bins)
        # Match the clamp applied in fit_quad_model so the plot shows the same
        # values the fits actually see.
        sm_yt = clamp.(sm_yt, Float32(-3), Float32(3))
        n_out = length(sm_x0)
        plot!(p, sm_x0, sm_yt, seriestype=:scatter, alpha=0.5,
              markersize=2, color=scatter_colors[charge],
              label="Charge $charge (n=$n_out)")
        δ = Float32(1.0 / charge)
        dense_x1 = dense_x0 .+ δ
        if initial_model !== nothing
            init_preds = [F(initial_model.params, dense_x0[i], dense_x1[i]) for i in eachindex(dense_x0)]
            plot!(p, dense_x0, init_preds, lw=1.5, color=l2_colors[charge], alpha=0.8,
                  label="Charge $charge grid L2")
        end
        if quad_model isa RazoQuadModel
            preds = [F(quad_model.params, dense_x0[i], dense_x1[i]) for i in eachindex(dense_x0)]
            plot!(p, dense_x0, preds, lw=2, color=fit_colors[charge], alpha=0.9,
                  label="Charge $charge Razo (LM)")
        end
    end
    return p
end


"""
    process_quad_pipeline(initial_psms::DataFrame,
                         spectra::MassSpecData,
                         search_context::SearchContext,
                         results::QuadTuningSearchResults,
                         params::QuadTuningSearchParameters,
                         ms_file_idx::Int64,
                         window_width::Float64)::DataFrame

Process initial PSMs through the complete quadrupole analysis pipeline.

# Arguments
- `initial_psms`: Initial PSMs from library search
- `spectra`: MassSpecData object
- `search_context`: SearchContext for analysis
- `results`: QuadTuningSearchResults for storing intermediate data
- `params`: Quad tuning parameters
- `ms_file_idx`: File index
- `window_width`: Isolation window width

# Process
1. Get scan mapping for quadrupole transmission search
2. Perform quad transmission search with deconvolution
3. Filter and process results for quad modeling
4. Return processed DataFrame ready for model fitting

# Returns
Processed DataFrame containing quad transmission data
"""
function process_quad_pipeline(
    initial_psms::DataFrame,
    spectra::MassSpecData,
    search_context::SearchContext,
    results::QuadTuningSearchResults,
    params::QuadTuningSearchParameters,
    ms_file_idx::Int64,
    window_width::Float64
)::DataFrame

    scan_idx_to_prec_idx = get_scan_to_prec_idx(
        UInt32.(initial_psms[!, :scan_idx]),
        UInt32.(initial_psms[!, :precursor_idx]),
        getCenterMzs(spectra),
        getIsolationWidthMzs(spectra)
    )

    quad_psms = perform_quad_transmission_search(
        spectra,
        results,
        scan_idx_to_prec_idx,
        search_context,
        params,
        ms_file_idx
    )

    if isempty(quad_psms)
        return DataFrame(
            scan_idx = Int64[],
            precursor_idx = UInt32[],
            center_mz = Union{Float32, Missing}[],
            δ = Union{Float32, Missing}[],
            yt = Union{Float32, Missing}[],
            x0 = Union{Float32, Missing}[],
            x1 = Union{Float32, Missing}[],
            prec_charge = Union{UInt8, Missing}[],
            half_width_mz = Float32[]
        )
    end

    quad_psms[!,:charge] = getCharges(getCharge(getPrecursors(getSpecLib(search_context))), quad_psms[!,:precursor_idx])
    quad_psms = quad_psms[
        filter_quad_psms(
            quad_psms[!,:iso_idx],
            quad_psms[!,:n_matches],
            quad_psms[!,:weight],
            quad_psms[!,:charge],
            params
        ),:]

    processed_psms = process_quad_results(
        quad_psms,
        getPrecursors(getSpecLib(search_context)),
        getIsoSplines(first(getSearchData(search_context)))
    )

    processed_psms[!,:half_width_mz] = zeros(Float32, size(processed_psms, 1))
    for (i, scan_idx) in enumerate(processed_psms[!,:scan_idx])
        processed_psms[i,:half_width_mz] = getIsolationWidthMz(spectra, scan_idx)/2
    end
    keep_data = zeros(Bool, size(processed_psms, 1))
    for i in range(1, size(processed_psms, 1))
        x0 = processed_psms[i,:x0]::Float32
        hw = processed_psms[i,:half_width_mz]::Float32
        if x0 > zero(Float32)
            if (x0 - hw) < (C13_C12_MASS_DIFF/4 + 0.1)
                keep_data[i] = true
            end
        else
            if (abs(x0) - hw) < (C13_C12_MASS_DIFF/2 + 0.1)
                keep_data[i] = true
            end
        end
    end
    processed_psms = processed_psms[keep_data,:]

    return processed_psms
end

"""
    equal_count_median(psms; n_bins=40) -> (bm_x0, bm_x1, bm_yt)

Sort raw (x0, x1, yt) by x0 and split into `n_bins` equal-count chunks;
return the IQR-trimmed mean of each column (mean of points between the 25th
and 75th percentiles, per column).

Using trimmed mean instead of median keeps more of the signal (uses ~50% of
the data rather than a single middle point) while still discarding outliers
and clipped values in the tails.
"""
function equal_count_median(psms::DataFrame; n_bins::Int = 40)
    x0 = psms[!, :x0]::Vector{Float32}
    x1 = psms[!, :x1]::Vector{Float32}
    yt = psms[!, :yt]::Vector{Float32}
    n = length(x0)
    n < n_bins && return (Float32[], Float32[], Float32[])

    perm = sortperm(x0)
    x0s = x0[perm]; x1s = x1[perm]; yts = yt[perm]

    # Float-valued chunk boundaries so counts are spread evenly.
    chunk = Float64(n) / Float64(n_bins)
    out_x0 = Vector{Float32}(undef, n_bins)
    out_x1 = Vector{Float32}(undef, n_bins)
    out_yt = Vector{Float32}(undef, n_bins)
    @inbounds for k in 1:n_bins
        lo = floor(Int, (k - 1) * chunk) + 1
        hi = floor(Int, k * chunk)
        hi = max(hi, lo)
        out_x0[k] = _iqr_trimmed_mean(@view x0s[lo:hi])
        out_x1[k] = _iqr_trimmed_mean(@view x1s[lo:hi])
        out_yt[k] = _iqr_trimmed_mean(@view yts[lo:hi])
    end
    out_x0, out_x1, out_yt
end

"""
    _iqr_trimmed_mean(v) -> Float32

Mean of values between Q1 and Q3 (25th/75th percentiles). For very small
vectors (≤3 elements) falls back to `median` since IQR trimming leaves
too few points.
"""
function _iqr_trimmed_mean(v::AbstractVector{T}) where {T<:AbstractFloat}
    n = length(v)
    n == 0 && return zero(Float32)
    n <= 3 && return Float32(median(v))
    sorted = sort(collect(v))
    q1_idx = max(1,  ceil(Int, 0.25 * n))
    q3_idx = min(n, floor(Int, 0.75 * n))
    q3_idx < q1_idx && return Float32(median(sorted))
    s = zero(Float64)
    @inbounds for i in q1_idx:q3_idx
        s += sorted[i]
    end
    return Float32(s / (q3_idx - q1_idx + 1))
end

function fit_quad_model(psms::DataFrame, window_width::Float64;
                         bin_resolution::Float64 = 0.0625,
                         pad::Float64 = 0.5)
    # One bin per `bin_resolution` Da of m/z coverage (window_width + 2*pad),
    # floored at 8. Bins are equal-count chunks (sorted by x0) inside that total.
    n_bins = max(8, ceil(Int, (window_width + 2 * pad) / bin_resolution))
    sm_x0, sm_x1, sm_yt = equal_count_median(psms; n_bins=n_bins)
    # Clamp bin summaries to |yt| ≤ 3 (natural log; ~20:1 transmission ratio)
    # so extreme tail bins don't drag the fit into unphysical parameters.
    sm_yt = clamp.(sm_yt, Float32(-3), Float32(3))
    n_pts = length(sm_x0)
    if n_pts > 0
        @debug_l1 "  QuadTuning fit_quad_model: $(nrow(psms)) raw → $n_pts equal-count bins (bin_resolution=$bin_resolution), x0 ∈ [$(round(minimum(sm_x0), digits=2)), $(round(maximum(sm_x0), digits=2))]"
    end

    if n_pts < 4
        error("Razo quad fit needs ≥4 points, got $n_pts (from $(nrow(psms)) raw points). Increase bin_resolution or check PSM count.")
    end

    # L2 grid search → LM warm-start for the 4-param Razo model.
    initial_params = grid_search_razo_init(sm_x0, sm_x1, sm_yt, window_width;
                                            a_sweep_half=1.0, a_step=0.025, n_b=60)

    fitted_params = fitRazoQuadModel(
        initial_params.al, initial_params.ar, initial_params.bl, initial_params.br,
        (Float32(0.2), Float32(window_width)),
        (Float32(0.2), Float32(window_width)),
        (Float32(1e-3), Float32(24.0)),
        (Float32(1e-3), Float32(24.0)),
        sm_x0, sm_x1, sm_yt;
        weights=nothing,
        λ0=1e-2,
        ϵ1=1e-5,
        ϵ2=1e-4,
        ϵ3=1e-5,
        max_iter=500
    )
    return fitted_params, initial_params
end
