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
    window_widths = Set{String}()
    for i in 1:length(spectra)
        if getMsOrder(spectra, i) == 2 && !ismissing(getIsolationWidthMz(spectra, i))
            push!(window_widths, string(getIsolationWidthMz(spectra, i)))
        end
    end
    return window_widths
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
    #Maps scans to the list of precursors in the scan 
    scan_idx_to_prec_idx = Dictionary{UInt32, Vector{UInt32}}()
    for i in range(1, N)
        scan_idx = scan_idxs[i]
        prec_idx = prec_idxs[i]
        #Have encountered scan 
        if haskey(scan_idx_to_prec_idx, scan_idx)
            push!(scan_idx_to_prec_idx[scan_idx], prec_idx)#(zip(psms[!,:scan_idx], psms[!,:precursor_idx]))
            prev_scan_idx, next_scan_idx = get_nearest_adjacent_scans(
                scan_idx, centerMz, isolationWidthMz
            )
            #Have encountered nearest scan 
            if haskey(scan_idx_to_prec_idx, next_scan_idx)
                push!(scan_idx_to_prec_idx[next_scan_idx], prec_idx)#(zip(psms[!,:scan_idx], psms[!,:precursor_idx]))
            else
                insert!(scan_idx_to_prec_idx, next_scan_idx,[prec_idx])
            end
            if haskey(scan_idx_to_prec_idx, prev_scan_idx)
                push!(scan_idx_to_prec_idx[prev_scan_idx], prec_idx)#(zip(psms[!,:scan_idx], psms[!,:precursor_idx]))
            else
                insert!(scan_idx_to_prec_idx, prev_scan_idx,[prec_idx])
            end
        else
            insert!(scan_idx_to_prec_idx, scan_idx,[prec_idx])
            prev_scan_idx, next_scan_idx = get_nearest_adjacent_scans(
                scan_idx, centerMz, isolationWidthMz
            )
            if haskey(scan_idx_to_prec_idx, next_scan_idx)
                push!(scan_idx_to_prec_idx[next_scan_idx], prec_idx)#(zip(psms[!,:scan_idx], psms[!,:precursor_idx]))
            else
                insert!(scan_idx_to_prec_idx, next_scan_idx,[prec_idx])
            end
            if haskey(scan_idx_to_prec_idx, prev_scan_idx)
                push!(scan_idx_to_prec_idx[prev_scan_idx], prec_idx)#(zip(psms[!,:scan_idx], psms[!,:precursor_idx]))
            else
                insert!(scan_idx_to_prec_idx, prev_scan_idx,[prec_idx])
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
    iso_mz = Float32.(mono_mz .+ NEUTRON.*(iso_idx.-1.0f0)./prec_charge)
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
    collect_psms(spectra::MassSpecData, search_context::SearchContext,
                results::QuadTuningSearchResults, params::QuadTuningSearchParameters,
                ms_file_idx::Int64) -> DataFrame

Iteratively collect and process PSMs for quadrupole tuning analysis.

# Arguments
- `spectra`: MS/MS spectral data
- `search_context`: Search context containing spectral library
- `results`: Current quad tuning results
- `params`: Search parameters
- `ms_file_idx`: Index of MS file being processed

# Process
1. Performs library search to get initial PSMs
2. Processes initial PSMs with basic filters
3. Performs quad transmission search
4. Filters results for M0/M1 isotopes and quality criteria
5. Repeats until sufficient high-quality PSMs collected

# Returns
DataFrame containing processed PSMs suitable for quad model fitting.
"""
function collect_psms(
    spectra::MassSpecData,
    search_context::SearchContext,
    results::QuadTuningSearchResults,
    params::QuadTuningSearchParameters,
    ms_file_idx::Int64
)
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

    total_psms = DataFrame()
    n = 0
    unique_precursors = Set{UInt32}()
    function getCharges(prec_charges::AbstractVector{UInt8}, precursor_idx::AbstractVector{UInt32})
        charges = zeros(UInt8, length(precursor_idx))
        for i in range(1, length(precursor_idx))
            charges[i] = prec_charges[precursor_idx[i]]
        end
        return charges
    end
    for i in range(1, 10)
        # Get initial PSMs
        psms = library_search(spectra, search_context, params, ms_file_idx)
        isempty(psms) && return total_psms
        psms = psms[[pid∉unique_precursors for pid in psms[!,:precursor_idx]],:]
        unique_precursors = union(unique_precursors, Set(psms[!,:precursor_idx]))

        # Process PSMs
        processed_psms = process_initial_psms(psms, spectra, search_context)
        
        # Get scan mapping and perform quad search
        scan_idx_to_prec_idx = get_scan_to_prec_idx(
            processed_psms[!, :scan_idx],
            processed_psms[!, :precursor_idx],
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

        
        # Filter and process results
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
                if (x0 - hw) < (NEUTRON/4 + 0.1)
                    keep_data[i] = true
                end
            else
                if (abs(x0) - hw) < (NEUTRON/2 + 0.1)
                    keep_data[i] = true
                end
            end
        end
        processed_psms = processed_psms[keep_data,:]

        append!(total_psms, processed_psms)
        if size(total_psms, 1) > params.min_quad_tuning_psms
            break
        end
    end
    return total_psms
end


"""
    process_initial_psms(psms::DataFrame, spectra::MassSpecData,
                        search_context::SearchContext) -> DataFrame

Process initial PSMs from library search for quad tuning.

# Arguments
- `psms`: Raw PSMs from library search
- `spectra`: MS/MS spectral data
- `search_context`: Search context with spectral library

# Process
1. Adds pre-search columns (target/decoy, RT, charge, etc.)
2. Scores PSMs using presearch scoring
3. Filters by q-value (≤ 0.01) and target status
4. Selects best PSM per precursor

# Returns
Filtered DataFrame containing high-confidence PSMs.
"""
function process_initial_psms(
    psms::DataFrame,
    spectra::MassSpecData,
    search_context::SearchContext
)
    add_tuning_search_columns!(
        psms,
        spectra,
        getIsDecoy(getPrecursors(getSpecLib(search_context))),#[:is_decoy],
        getIrt(getPrecursors(getSpecLib(search_context))),#[:irt],
        getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge],
        getRetentionTimes(spectra),
        getTICs(spectra)
    )
    
    score_presearch!(psms)
    get_qvalues!(psms[!,:prob], psms[!,:target], psms[!,:q_value])

    filter!(:q_value => x -> x <= 0.01, psms)
    filter!(:target => identity, psms)
    
    psms[!, :best_psms] .= false
    for group in groupby(psms, :precursor_idx)
        best_idx = argmax(group.prob)
        group[best_idx, :best_psms] = true
    end
    
    filter!(row -> row.best_psms, psms)
    return psms
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
        if ((iso_idx[1] == 1) & (iso_idx[2] == 2))
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
        else
            println("iso_idx $iso_idx")
        end
    end
    if length(iso_idx) == 3
        println("bro... $iso_idx")
    end
    
    #If we only got the M0
    if (length(iso_idx) == 1) & (iso_idx[1] == 1)
        m0_idx = 1
        m1_mz = iso_mz[m0_idx] + (NEUTRON/prec_charge[m0_idx])

        return (center_mz = center_mz[m0_idx], 
            δ = δ[m0_idx], 
            yt = 10.0f0, # assume this is the most extreme ratio we could see
            x0 = iso_mz[m0_idx]-center_mz[m0_idx],
            x1 = m1_mz - center_mz[m0_idx],
            prec_charge = prec_charge[m0_idx]) 
    end

    #If we only got the M1
    if (length(iso_idx) == 1) & (iso_idx[1] == 2)
        m1_idx = 1
        mo_mz = iso_mz[m1_idx] - (NEUTRON/prec_charge[m1_idx])

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
    Perform deconvolution for a single scan.
    """
    function perform_deconvolution!(
        Hs::SparseArray,
        weights::Vector{Float32},
        precursor_weights::Vector{Float32},
        residuals::Vector{Float32},
        search_data::SearchDataStructures,
        nmatches::Int,
        nmisses::Int,
        stop_tolerance::Real,
        params::QuadTuningSearchParameters)
        buildDesignMatrix!(
            Hs,
            getIonMatches(search_data),
            getIonMisses(search_data),
            nmatches,
            nmisses,
            getIdToCol(search_data)
        )
        
        # Resize arrays if needed
        if getIdToCol(search_data).size > length(weights)
            new_entries = getIdToCol(search_data).size - length(weights) + 1000
            resize!(weights, length(weights) + new_entries)
            resize!(getSpectralScores(search_data), length(getSpectralScores(search_data)) + new_entries)
            append!(getUnscoredPsms(search_data), [eltype(getUnscoredPsms(search_data))() for _ in 1:new_entries])
        end
        
        # Initialize weights
        for i in 1:getIdToCol(search_data).size
            weights[getIdToCol(search_data)[getIdToCol(search_data).keys[i]]] = 
                precursor_weights[getIdToCol(search_data).keys[i]]
        end
        
        # Solve deconvolution problem
        initResiduals!(residuals, Hs, weights)
        solveHuber!(
            Hs,
            residuals,
            weights,
            Float32(100000.0f0),
            0.0f0,
            params.max_iter_newton,
            params.max_iter_bisection,
            params.max_iter_outer,
            params.accuracy_newton,
            params.accuracy_bisection,
            params.max_diff,
            NoNorm()
        )
    end

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
        Hs::SparseArray,
        scan_idx::Int,
        center_mz::Float32
    )
        #tuning_results = getTuningResults(search_data)
        
        for i in 1:getIdToCol(search_data).size
            id = getIdToCol(search_data).keys[i]
            colid = getIdToCol(search_data)[id]
            
            # Update precursor weights
            search_data.precursor_weights[id] = weights[colid]
            
            # Calculate indices
            isotope_idx = UInt8(((id - 1) % 3) + 1)
            pid = UInt32(((id - 1) ÷ 3) + 1)
            
            # Count matches
            n_matches = sum(Hs.matched[j] 
                for j in Hs.colptr[colid]:(Hs.colptr[colid+1] - 1))
            
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
    Process a single scan for quad transmission search.
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
        Hs::SparseArray,
        weights::Vector{Float32},
        precursor_weights::Vector{Float32},
        residuals::Vector{Float32},
        scan_idx_to_prec_idx::Dictionary{UInt32, Vector{UInt32}}  # Added this parameter
    )
        scan_idx ∉ scan_idxs && return
        
        msn = getMsOrder(spectra, scan_idx)
        msn ∉ params.spec_order && return
        
        # Select transitions
        ion_idx, _ = selectTransitions!(
            getIonTemplates(search_data),
            QuadEstimationTransitionSelection(),
            PartialPrecCapture(),
            getFragmentLookupTable(getSpecLib(search_context)),
            scan_idx_to_prec_idx[scan_idx],
            getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
            getCharge(getPrecursors(getSpecLib(search_context))),#[:prec_charge],
            getSulfurCount(getPrecursors(getSpecLib(search_context))),#[:sulfur_count],
            getIsoSplines(search_data),
            getPrecursorTransmission(search_data),
            getIsotopes(search_data),
            (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx));
            block_size = 10000
        )
        
        # Match peaks
        nmatches, nmisses = matchPeaks!(
            getIonMatches(search_data),
            getIonMisses(search_data),
            getIonTemplates(search_data),
            ion_idx,
            getMzArray(spectra, scan_idx),
            getIntensityArray(spectra, scan_idx),
            getMassErrorModel(search_context, ms_file_idx),
            getHighMz(spectra, scan_idx),
            UInt32(scan_idx),
            UInt32(ms_file_idx)
        )
        
        nmatches ≤ 2 && return
        
        # Sort matches
        sort!(@view(getIonMatches(search_data)[1:nmatches]), 
            by = x->(x.peak_ind, x.prec_id),
            alg=QuickSort)
        
        # Build design matrix and deconvolve
        perform_deconvolution!(
            Hs,
            weights,
            precursor_weights,
            residuals,
            search_data,
            nmatches,
            nmisses,
            search_context.deconvolution_stop_tolerance[],
            params
        )
        
        # Record results
        record_scan_results!(
            search_data,
            tuning_results,
            weights,
            Hs,
            scan_idx,
            getCenterMz(spectra, scan_idx)
        )
        
        # Reset for next scan
        reset!(getIdToCol(search_data))
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
            # Get working arrays
            Hs = getHs(search_data)
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
    target_size = length(getPrecursors(getSpecLib(search_context))) * 3 + 1

    # Only adjust if not already at target size
    if length(first(getSearchData(search_context)).precursor_weights) == length(getPrecursors(getSpecLib(search_context)))
        for search_data in getSearchData(search_context)
            search_data.id_to_col = ArrayDict(UInt32, UInt16, target_size)
            
            if length(search_data.precursor_weights) == length(getPrecursors(getSpecLib(search_context)))
                resize!(search_data.precursor_weights, target_size)
                search_data.precursor_weights[length(search_data.precursor_weights):end] .= zero(Float32)
            end
        end
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
    original_size = length(getPrecursors(getSpecLib(search_context)))

    if length(first(getSearchData(search_context)).precursor_weights) != original_size
        for search_data in getSearchData(search_context)
            search_data.id_to_col = ArrayDict(UInt32, UInt16, original_size)
            resize!(search_data.precursor_weights, original_size)
        end
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

function plot_charge_distributions(psms::DataFrame, results::QuadTuningSearchResults, fname::String)
    p = plot(title = "Quad Model Data for $fname")
    for charge in 2:3
        mask = psms[!, :prec_charge] .== charge
        plot!(p, 
            psms[mask, :x0],
            psms[mask, :yt],
            seriestype=:scatter,
            alpha=0.1,
            label="Charge $charge",
            xlabel = "m/z offset of M0",
            ylabel = L"log(\delta_i\frac{x_0}{x_1})",
            ylim = (-10.5, 10.5)
        )
    end
   return p
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

function plot_quad_model(quad_model::QuadTransmissionModel, window_width::Float64, results::QuadTuningSearchResults, fname::String)
    padding = 2
    half_width = padding + window_width/2
    plot_bins = LinRange(-half_width, half_width, 100)

    quad_func = getQuadTransmissionFunction(quad_model, 0.0f0, 2.0f0)
    p = plot(plot_bins, quad_func.(plot_bins), lw = 2, alpha = 0.5, title = "$fname")
    return p
end


"""
    fit_quad_model(psms::DataFrame, window_width::Float64) -> QuadTransmissionModel

Fit quadrupole transmission model to binned PSM data.

# Arguments
- `psms`: Processed PSMs with transmission data
- `window_width`: Isolation window width

# Process
1. Bins PSMs by m/z offset
2. Fits Razo quad model to binned data
3. Uses regularization parameters for stability

# Returns
Fitted quadrupole transmission model.
"""
function fit_quad_model(psms::DataFrame, window_width::Float64)
    binned_psms = MergeBins(
        psms,
        (-(window_width + 1.0), window_width + 1.0),
        min_bin_size=20,
        min_bin_width=0.1
    )

    return fitRazoQuadModel(
        window_width,
        binned_psms[!, :median_x0],
        binned_psms[!, :median_x1],
        binned_psms[!, :median_yt],
        λ0=1e-2,
        ϵ1=1e-5,
        ϵ2=1e-4,
        ϵ3=1e-5
    )
end



