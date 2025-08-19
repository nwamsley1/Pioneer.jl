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
    add_tuning_search_columns!(psms::DataFrame, MS_TABLE::Arrow.Table,
                             prec_is_decoy::Arrow.BoolVector{Bool},
                             prec_irt::Arrow.Primitive{T, Vector{T}},
                             prec_charge::Arrow.Primitive{UInt8, Vector{UInt8}},
                             scan_retention_time::AbstractVector{Float32},
                             tic::AbstractVector{Float32}) where {T<:AbstractFloat}

Adds essential columns to PSM DataFrame for parameter tuning analysis.

# Arguments
- `psms`: DataFrame containing PSMs to modify
- `MS_TABLE`: Mass spectrometry data table
- `prec_is_decoy`: Boolean vector indicating decoy status
- `prec_irt`: Vector of iRT values
- `prec_charge`: Vector of precursor charges
- `scan_retention_time`: Vector of scan retention times
- `tic`: Vector of total ion currents

# Added Columns
- Basic metrics: RT, iRT predicted, charge, TIC
- Analysis columns: target, decoy, matched_ratio
- Scoring columns: q_value, prob, intercept

Uses parallel processing for efficiency through data chunking.
"""
function add_tuning_search_columns!(psms::DataFrame, 
                                MS_TABLE::MassSpecData, 
                                prec_is_decoy::Arrow.BoolVector{Bool},
                                prec_irt::Arrow.Primitive{T, Vector{T}},
                                prec_charge::Arrow.Primitive{UInt8, Vector{UInt8}},
                                scan_retention_time::AbstractVector{Float32},
                                tic::AbstractVector{Float32}) where {T<:AbstractFloat}
    
    N = size(psms, 1)
    decoys = zeros(Bool, N);
    targets = zeros(Bool, N);
    TIC = zeros(Float16, N);
    charge = zeros(UInt8, N);
    spectrum_peak_count = zeros(UInt32, N);
    irt_pred = zeros(Float32, N);
    rt = zeros(Float32, N);
    scan_idx::Vector{UInt32} = psms[!,:scan_idx]
    precursor_idx::Vector{UInt32} = psms[!,:precursor_idx]
    matched_ratio::Vector{Float16} = psms[!,:matched_ratio]
    
    tasks_per_thread = 10
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
            decoys[i] = prec_is_decoy[precursor_idx[i]];
            targets[i] = decoys[i] == false
            irt_pred[i] = Float32(prec_irt[precursor_idx[i]]);
            rt[i] = Float32(scan_retention_time[scan_idx[i]]);
            TIC[i] = Float16(log2(tic[scan_idx[i]]));
            charge[i] = UInt8(prec_charge[precursor_idx[i]]);
            matched_ratio[i] = Float16(min(matched_ratio[i], 6e4))
            end
        end
    end
    fetch.(tasks)
    psms[!,:matched_ratio] = matched_ratio
    psms[!,:decoy] = decoys
    psms[!,:irt_predicted] = irt_pred
    psms[!,:rt] = rt
    psms[!,:TIC] = TIC
    psms[!,:target] = targets
    psms[!,:charge] = charge
    psms[!,:spectrum_peak_count] = spectrum_peak_count
    psms[!,:intercept] = ones(Float16, size(psms, 1))
    psms[!,:q_value] = zeros(Float16, N);
    psms[!,:prob] = zeros(Float16, N);
end


#==========================================================
PSM scoring and filtering
==========================================================#
"""
    filter_and_score_psms!(psms::DataFrame, params::P) where {P<:ParameterTuningSearchParameters}

Filters PSMs based on score and selects best matches per precursor.

# Arguments
- `psms`: DataFrame containing PSMs
- `params`: Parameter tuning search parameters

# Process
1. Scores PSMs using presearch scoring
2. Calculates q-values
3. Filters by q-value threshold
4. Selects best PSM per precursor based on probability score

# Returns
- Number of PSMs after all filtering (q-value, best per precursor, targets only) if successful
- -1 if insufficient PSMs pass q-value threshold (DataFrame is NOT filtered in this case)

Modifies psms DataFrame in place by filtering to best matches ONLY if successful.
If returns -1, DataFrame contains scored but unfiltered PSMs.
"""
function filter_and_score_psms!(
    psms::DataFrame,
    params::P,
    search_context::SearchContext
) where {P<:ParameterTuningSearchParameters}
    
    score_presearch!(psms)
    # Get FDR scale factor from search context
    fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
    get_qvalues!(psms[!,:prob], psms[!,:target], psms[!,:q_value]; fdr_scale_factor = fdr_scale_factor)
    
    max_q_val, min_psms = getMaxQVal(params), getMinPsms(params)
    n_passing_psms = sum(psms[!,:q_value] .<= max_q_val)
    
    if n_passing_psms >= min_psms
        filter!(row -> row.q_value::Float16 <= max_q_val, psms)
        psms[!,:best_psms] = zeros(Bool, size(psms, 1))
        
        # Select best PSM per precursor
        for sub_psms in groupby(psms, :precursor_idx)
            best_idx = argmax(sub_psms.prob::AbstractVector{Float32})
            sub_psms[best_idx,:best_psms] = true
        end
        
        filter!(x -> x.best_psms::Bool, psms)
        filter!(x->x.target::Bool, psms) #Otherwise fitting rt/irt and mass tolerance partly on decoys.
        
        # Return actual count after all filtering (not n_passing_psms which was pre-filtering)
        return size(psms, 1)
    end
    
    # Filtering failed - DataFrame remains with scored but unfiltered PSMs
    return -1
end


"""
    score_presearch!(psms::DataFrame)

Performs initial scoring of PSMs using probit regression.

# Arguments
- `psms`: DataFrame containing PSMs to score

# Process
1. Uses features: entropy_score, city_block, scribe, spectral_contrast,
   y_count, error, TIC, intercept
2. Performs probit regression with 20 max iterations
3. Calculates probability scores for each PSM

Adds 'prob' column to psms DataFrame containing probability scores.
Uses parallel processing for efficiency.
"""
function score_presearch!(psms::DataFrame)
    features = [:entropy_score,:city_block,
                    :scribe,:spectral_contrast,
                    :y_count,:error,
                    :TIC,:intercept]
    psms[!,:prob] = zeros(Float32, size(psms, 1))

    M = size(psms, 1)
    if M > 10
        tasks_per_thread = 10
        chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
        data_chunks = partition(1:M, chunk_size) # partition your data into chunks that
        β = zeros(Float64, length(features))
        β = ProbitRegression(β, psms[!,features], psms[!,:target], data_chunks, max_iter = 20)
        ModelPredictProbs!(psms[!,:prob], psms[!,features], β, data_chunks)
    else
        psms[!,:prob] = ones(Float32, size(psms, 1))
    end
end


# Note: get_qvalues! and get_PEP! functions have been moved to src/utils/ML/fdrUtilities.jl
# They are imported through the module loading system and are available in the Pioneer module namespace

#==========================================================
RT modeling
==========================================================#
"""
    fit_irt_model(params::P, psms::DataFrame) where {P<:ParameterTuningSearchParameters}

Fits retention time alignment model between library and empirical retention times.

# Arguments
- `params`: Parameter tuning search parameters
- `psms`: DataFrame containing PSMs with RT information

# Process
1. Performs initial spline fit between predicted and observed RTs
2. Calculates residuals and median absolute deviation
3. Removes outliers based on MAD threshold
4. Refits spline model on cleaned data

# Returns
Tuple containing:
- Final RT conversion model
- Valid RT values
- Valid iRT values
- iRT median absolute deviation
"""
function fit_irt_model(
    params::P,
    psms::DataFrame
) where {P<:ParameterTuningSearchParameters}
    
    # Initial spline fit
    rt_to_irt_map = UniformSpline(
        psms[!,:irt_predicted],
        psms[!,:rt],
        getSplineDegree(params),
        getSplineNKnots(params)
    )
    
    # Calculate residuals
    psms[!,:irt_observed] = rt_to_irt_map.(psms.rt::Vector{Float32})
    residuals = psms[!,:irt_observed] .- psms[!,:irt_predicted]
    irt_mad = mad(residuals, normalize=false)::Float32
    
    # Remove outliers and refit
    valid_psms = psms[abs.(residuals) .< (irt_mad * getOutlierThreshold(params)), :]
    
    final_model = SplineRtConversionModel(UniformSpline(
        valid_psms[!,:irt_predicted],
        valid_psms[!,:rt],
        getSplineDegree(params),
        getSplineNKnots(params)
    ))
    
    return (final_model, valid_psms[!,:rt], valid_psms[!,:irt_predicted], irt_mad)
end

"""
    fit_mass_err_model(params::P, fragments::Vector{FragmentMatch{Float32}}) where {P<:FragmentIndexSearchParameters}

Fits mass error model from fragment matches.

# Arguments
- `params`: Search parameters
- `fragments`: Vector of fragment matches

# Process
1. Calculates PPM errors for all matches
2. Determines systematic mass error offset
3. Calculates error bounds based on quantiles

# Returns
Tuple containing:
- MassErrorModel with offset and tolerance bounds
- Vector of PPM errors
"""

#==========================================================
Mass Error Modeling
==========================================================#

"""
    adjustMatchMz(match::FragmentMatch{T}, ppm_offset::Float32) where T<:AbstractFloat

Adjusts the match_mz field of a FragmentMatch by a given PPM offset.

# Arguments
- `match`: FragmentMatch to adjust
- `ppm_offset`: Parts-per-million offset to apply

# Returns
- New FragmentMatch with adjusted match_mz field
- Formula: adjusted_mz = match_mz * (1 + ppm_offset/1e6)
"""
function adjustMatchMz(match::FragmentMatch{T}, ppm_offset::Float32) where T<:AbstractFloat
    adjusted_mz = getMatchMZ(match) * (1.0f0 + ppm_offset/(Float32(1.0e6)))
    return FragmentMatch(
        getPredictedIntensity(match),
        getIntensity(match),
        getMZ(match),  # theoretical_mz stays the same
        adjusted_mz,   # adjusted match_mz
        getPeakInd(match),
        getFragInd(match),
        getCharge(match),
        getIsotope(match),
        getIonType(match),
        isIsotope(match),
        getPrecID(match),
        getCount(match),
        getScanID(match),
        getMSFileID(match),
        getRank(match)
    )
end

"""
    adjustMatchMz(match::PrecursorMatch{T}, ppm_offset::Float32) where T<:AbstractFloat

Adjusts the observed_mz field of a PrecursorMatch by a given PPM offset.

# Arguments
- `match`: PrecursorMatch to adjust
- `ppm_offset`: Parts-per-million offset to apply

# Returns
- New PrecursorMatch with adjusted observed_mz field
- Formula: adjusted_mz = observed_mz * (1 + ppm_offset/1e6)
"""
function adjustMatchMz(match::PrecursorMatch{T}, ppm_offset::Float32) where T<:AbstractFloat
    adjusted_mz = getMatchMz(match) * (1.0f0 + ppm_offset/1.0e6f0)
    return PrecursorMatch(
        getPredictedIntensity(match),
        getIntensity(match),
        getMZ(match),  # theoretical_mz stays the same
        adjusted_mz,   # adjusted observed_mz
        getIsoIdx(match),
        getPeakInd(match),
        getPrecID(match)
    )
end

function getScanToPrecIdx(scan_idxs::Vector{UInt32}, n_scans::Int64)
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, n_scans)
    start_idx = stop_idx = 1
    
    for i in 1:length(scan_idxs)
        stop_idx = i
        if scan_idxs[start_idx] == scan_idxs[stop_idx]
            scan_to_prec_idx[scan_idxs[i]] = start_idx:stop_idx
        else
            scan_to_prec_idx[scan_idxs[i]] = i:i
            start_idx = i
        end
    end
    scan_to_prec_idx
end

"""
    collectFragErrs(fmatches::Vector{M}, new_fmatches::Vector{M}, nmatches::Int, n::Int; 
                   ppm_offset::Float32 = 0.0f0) where {M<:MatchIon{Float32}}

Collects fragment matches from new_fmatches into fmatches, optionally adjusting match_mz by a PPM offset.

# Arguments
- `fmatches`: Vector to store accumulated fragment matches
- `new_fmatches`: Vector of new fragment matches to add
- `nmatches`: Number of matches to copy from new_fmatches
- `n`: Current index in fmatches
- `ppm_offset`: Optional PPM offset to apply to match_mz (default: 0.0)

# Returns
- Updated index n after adding matches

# Notes
- When ppm_offset is 0, performs direct copy (original behavior)
- When ppm_offset is non-zero, adjusts match_mz before copying
- Grows fmatches vector as needed
"""
function collectFragErrs(fmatches::Vector{M}, new_fmatches::Vector{M}, nmatches::Int, n::Int, 
                        ppm_offset::Float32) where {M<:MatchIon{Float32}}
    for match in range(1, nmatches)
        if n < length(fmatches)
            n += 1
            if iszero(ppm_offset)# == 0.0f0
                # Direct copy when no adjustment needed (original behavior)
                fmatches[n] = new_fmatches[match]
            else
                # Apply PPM adjustment to match_mz
                fmatches[n] = adjustMatchMz(new_fmatches[match],ppm_offset)
            end
        else
            fmatches = append!(fmatches, [M() for x in range(1, length(fmatches))])
        end
    end
    return n
end

"""
    mass_error_search(spectra::MassSpecData, scan_idxs::Vector{UInt32},
                     precursor_idxs::Vector{UInt32}, ms_file_idx::UInt32,
                     spec_lib::SpectralLibrary, search_data::AbstractVector{S},
                     mem::M, params::P) where {M<:MassErrorModel, S<:SearchDataStructures, P<:SearchParameters}

Performs mass error-focused fragment matching.

# Arguments
- `spectra`: MS spectra data
- `scan_idxs`: Scan indices to process
- `precursor_idxs`: Precursor indices
- `ms_file_idx`: MS file index
- `spec_lib`: Spectral library
- `search_data`: Search data structures
- `mem`: Mass error model
- `params`: Search parameters

# Process
1. Maps scans to precursor indices
2. Selects transitions for mass error estimation
3. Matches peaks and collects errors
4. Processes matches in parallel across threads

Returns matched fragments for mass error analysis.
"""
# Default method without chromatogram type - calls MS2CHROM version
function mass_error_search(
    spectra::MassSpecData,
    scan_idxs::Vector{UInt32},
    precursor_idxs::Vector{UInt32},
    ms_file_idx::UInt32,
    spec_lib::SpectralLibrary,
    search_data::AbstractVector{S},
    mem::M,
    params::P
) where {
        M<:MassErrorModel, 
        S<:SearchDataStructures, 
        P<:SearchParameters
        }
    return mass_error_search(spectra, scan_idxs, precursor_idxs, ms_file_idx, 
                            spec_lib, search_data, mem, params, MS2CHROM())
end

function mass_error_search(
    spectra::MassSpecData,
    scan_idxs::Vector{UInt32},
    precursor_idxs::Vector{UInt32},
    ms_file_idx::UInt32,
    spec_lib::SpectralLibrary,
    search_data::AbstractVector{S},
    mem::M,
    params::P,
    ::MS2CHROM
) where {
        M<:MassErrorModel, 
        S<:SearchDataStructures, 
        P<:SearchParameters
        }

    # Sort scans and setup thread tasks
    sorted_indices = sortperm(scan_idxs)
    scan_idxs = scan_idxs[sorted_indices]
    precursor_idxs = precursor_idxs[sorted_indices]
    
    scan_to_prec_idx = getScanToPrecIdx(scan_idxs, length(spectra))
    thread_tasks = partition_scans(spectra, Threads.nthreads())

    # Process mass errors in parallel
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            frag_err_idx = 0
            
            for scan_idx in last(thread_task)
                (scan_idx == 0 || scan_idx > length(spectra)) && continue
                ismissing(scan_to_prec_idx[scan_idx]) && continue

                # Select transitions for mass error estimation
                # Uses top N fragments per precursor for mass tolerance calibration
                ion_idx, _ = selectTransitions!(
                    getIonTemplates(search_data[thread_id]),
                    MassErrEstimationStrategy(),
                    getPrecEstimation(params),
                    getFragmentLookupTable(spec_lib),
                    scan_to_prec_idx[scan_idx],
                    precursor_idxs,
                    max_rank = Int64(getMaxFragsForMassErrEstimation(params))  # Convert UInt8 to Int64 for compatibility
                )

                # Match peaks and collect errors
                nmatches, nmisses = matchPeaks!(
                    getIonMatches(search_data[thread_id]),
                    getIonMisses(search_data[thread_id]),
                    getIonTemplates(search_data[thread_id]),
                    ion_idx,
                    getMzArray(spectra, scan_idx),
                    getIntensityArray(spectra, scan_idx),
                    mem,
                    getHighMz(spectra, scan_idx),
                    UInt32(scan_idx),
                    ms_file_idx
                )

                frag_err_idx = collectFragErrs(
                    getMassErrMatches(search_data[thread_id]),
                    getIonMatches(search_data[thread_id]),
                    nmatches,
                    frag_err_idx,
                    getMassOffset(mem)
                )
            end
            
            @view(getMassErrMatches(search_data[thread_id])[1:frag_err_idx])
        end
    end
    fetch.(tasks)
end

function mass_error_search(
    spectra::MassSpecData,
    scan_idxs::Vector{UInt32},
    precursor_idxs::Vector{UInt32},
    ms_file_idx::UInt32,
    spec_lib::SpectralLibrary,
    search_data::AbstractVector{S},
    mem::M,
    params::P,
    ::MS1CHROM
) where {
        M<:MassErrorModel, 
        S<:SearchDataStructures, 
        P<:SearchParameters
        }

    # Sort scans and setup thread tasks
    sorted_indices = sortperm(scan_idxs)
    scan_idxs = scan_idxs[sorted_indices]
    precursor_idxs = precursor_idxs[sorted_indices]
    prec_mzs = getMz(getPrecursors(spec_lib))
    prec_charges = getCharge(getPrecursors(spec_lib))
    #Convert MS2 scan idxs to nearest MS1 scan idxs
    prev_ms1_scan_idx = 1
    next_ms1_scan_idx = 1
    ms1_scan_idx = 1
    for (i, scan_idx) in enumerate(scan_idxs)
        if scan_idx > next_ms1_scan_idx
            while (ms1_scan_idx < length(spectra))
                if getMsOrder(spectra, ms1_scan_idx) == 1
                    prev_ms1_scan_idx = next_ms1_scan_idx
                    next_ms1_scan_idx = ms1_scan_idx 
                    if ms1_scan_idx  > scan_idx 
                        break
                    end
                end
                ms1_scan_idx += 1
            end
        end
        if abs(next_ms1_scan_idx - scan_idx) < abs(prev_ms1_scan_idx - scan_idx)
            scan_idxs[i] = next_ms1_scan_idx
        else
            scan_idxs[i] =  prev_ms1_scan_idx
        end
    end
    scan_to_prec_idx = getScanToPrecIdx(scan_idxs, length(spectra))
    thread_tasks = partition_scans(spectra, Threads.nthreads(), ms_order_select = 1)
    #println("thread_tasks $thread_tasks")
    # Process mass errors in parallel
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            mass_errs = Vector{Float32}(undef, 1000)
            peak_intensities = Vector{Float32}(undef, 1000)
            ion_templates = Vector{Isotope{Float32}}(undef, 10000)
            ion_matches = [PrecursorMatch{Float32}() for _ in range(1, 1000)]
            ion_misses = [PrecursorMatch{Float32}() for _ in range(1, 1000)]
            ion_idx = 0
            frag_err_idx = 0
            ion_idx = 0


            prec_to_match_dict = Dictionary{
                UInt32, #predc_id
                UInt8 #Number of matches to spectrum 
            }() 
            for scan_idx in last(thread_task)
                empty!( prec_to_match_dict)
                ion_idx = 0
                (scan_idx == 0 || scan_idx > length(spectra)) && continue
                ismissing(scan_to_prec_idx[scan_idx]) && continue
                prec_idx_range = nothing
                prec_idx_range = scan_to_prec_idx[scan_idx]

                #println("prec_idx_range $prec_idx_range")
                for i in prec_idx_range
                    prec_idx = precursor_idxs[i]
                    prec_mz = prec_mzs[prec_idx]
                    prec_charge = prec_charges[prec_idx]
                    for iso in range(0, 3)
                        ion_idx += 1
                        if ion_idx > length(ion_templates)
                            append!(ion_templates, Vector{Isotope{Float32}}(undef, length(ion_templates)))
                        end
                        mz = Float32(prec_mz + iso*NEUTRON/prec_charge)
                        ion_templates[ion_idx] = Isotope( #precursor monoisotope
                            mz,
                            0.0f0,
                            UInt8(iso+1),
                            prec_idx
                        )
                    end
                end

                #Sort the precursor isotopes by m/z
                sort!(@view(ion_templates[1:ion_idx]), by = x->(getMZ(x)), alg=PartialQuickSort(1:ion_idx))

                # Match peaks and collect errors
                nmatches, nmisses = matchPeaks!(
                    ion_matches,
                    ion_misses,
                    ion_templates,
                    ion_idx,
                    getMzArray(spectra, scan_idx),
                    getIntensityArray(spectra, scan_idx),
                    mem,
                    getHighMz(spectra, scan_idx),
                    UInt32(scan_idx),
                    UInt32(ms_file_idx)
                )

                #Which precursors matched isotopes
                for i in range(1, nmatches)
                    prec_idx = getPrecID(ion_matches[i])
                    if !haskey( prec_to_match_dict, prec_idx)
                        insert!(prec_to_match_dict, prec_idx, zero(UInt8))
                    end
                    n_match = prec_to_match_dict[prec_idx]
                    if getIsoIdx(ion_matches[i]) <= 4
                        n_match += 1
                    end
                    prec_to_match_dict[prec_idx] = n_match 
                end

                #Removed matched fragments for precursors that did not match sufficiently many isotopes 
                new_nmatches = 0
                for i in range(1, nmatches)
                    prec_idx = getPrecID(ion_matches[i])
                    n_match = prec_to_match_dict[prec_idx]
                    if n_match >= 4
                        new_nmatches += 1
                        ion_matches[new_nmatches] = ion_matches[i]
                    end
                end
                
                nmatches = new_nmatches

                for match_idx in range(1, nmatches)
                    match = ion_matches[match_idx]
                    if getIsoIdx(match)>3
                        continue
                    end
                    frag_err_idx += 1
                    if frag_err_idx > length(mass_errs)
                        append!(mass_errs, Vector{Float32}(undef, length(mass_errs)))
                        append!(peak_intensities, Vector{Float32}(undef, length(peak_intensities)))
                    end
                    mass_errs[frag_err_idx] = Float32((getMatchMz(match) - getMZ(match))/(getMZ(match)/1e6))
                    peak_intensities[frag_err_idx] = Float32(getIntensity(match))
                end
            end

            ########
            #Intensity filter to remove potentially erroneous matches 
            mass_errs = mass_errs[1:frag_err_idx]
            if frag_err_idx > 1
                peak_intensities = peak_intensities[1:frag_err_idx]
                med_intensity = quantile(peak_intensities, 0.75)
                return @view(mass_errs[peak_intensities.>med_intensity])
            else
                return mass_errs
            end
        end
    end
    fetch.(tasks)
end

function fit_mass_err_model(
    params::P,
    fragments::Vector{FragmentMatch{Float32}}
) where {P<:FragmentIndexSearchParameters}
    
    # Filter fragments by intensity to remove low-quality matches
    if length(fragments) > 0
        intensity_threshold = getIntensityFilterQuantile(params)
        if intensity_threshold > 0.0
            intensities = [fragment.intensity for fragment in fragments]
            min_intensity = quantile(intensities, intensity_threshold)
            fragments = filter(f -> f.intensity > min_intensity, fragments)
            # @info "Filtered to $(length(fragments)) fragments above $(round(intensity_threshold*100, digits=1))th percentile intensity"
        end
    end
    
    # Check if we have enough fragments after filtering
    if length(fragments) == 0
        @warn "No fragments remaining after intensity filtering"
        return MassErrorModel(0.0f0, (30.0f0, 30.0f0)), Float32[]
    end
    
    # Calculate PPM errors
    ppm_errs = [(f.match_mz - f.theoretical_mz)/(f.theoretical_mz/1e6) for f in fragments]
    mass_err = median(ppm_errs)
    ppm_errs .-= mass_err
    
    # Calculate error bounds using configured quantile
    frag_err_quantile = getFragErrQuantile(params)

    _mad = mad(ppm_errs)
    #r_mad = mad(ppm_errs[ppm_errs .> 0.0f0])
    l_bound = min(abs(quantile(ppm_errs, frag_err_quantile)), 5*_mad)
    r_bound = min(abs(quantile(ppm_errs, 1 - frag_err_quantile)), 5*_mad)
    
    # Diagnostic logging
    # @info "Mass error model: offset=$(round(mass_err, digits=2)) ppm, " *
    #       "tolerance=(-$(round(abs(l_bound), digits=1)), +$(round(abs(r_bound), digits=1))) ppm, " *
    #       "MAD=$(round(mad(ppm_errs, normalize=false)*4, digits=1)) ppm, " *
    #       "error_quantile=$frag_err_quantile"
    
    return MassErrorModel(
        Float32(mass_err),
        (Float32(abs(l_bound)), Float32(abs(r_bound)))
    ), ppm_errs
end

function mass_err_ms1(ppm_errs::Vector{Float32},
    params::P) where {P<:FragmentIndexSearchParameters}
    mass_err = median(ppm_errs)
    ppm_errs .-= mass_err
    # Calculate error bounds
    frag_err_quantile = getFragErrQuantile(params)
    l_bound = quantile(ppm_errs, frag_err_quantile)
    r_bound = quantile(ppm_errs, 1 - frag_err_quantile)

    return MassErrorModel(
        Float32(mass_err),
        (Float32(abs(l_bound)), Float32(abs(r_bound)))
    ), ppm_errs
end
"""
    get_matched_fragments(spectra::MassSpecData, psms::DataFrame,
                         results::ParameterTuningSearchResults,
                         search_context::SearchContext,
                         params::P, ms_file_idx::Int64) where {P<:FragmentIndexSearchParameters}

Collects matched fragments for mass error estimation.

# Arguments
- `spectra`: MS spectra data
- `psms`: Identified PSMs
- `results`: Current parameter tuning results
- `search_context`: Search context
- `params`: Search parameters
- `ms_file_idx`: Index of MS file being processed

Returns vector of fragment matches using mass error search strategy.
"""
function get_matched_fragments(
    spectra::MassSpecData,
    psms::DataFrame,
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64
) where {P<:FragmentIndexSearchParameters}
    
    return vcat(mass_error_search(
        spectra,
        psms[!,:scan_idx],
        psms[!,:precursor_idx],
        UInt32(ms_file_idx),
        getSpecLib(search_context),
        getSearchData(search_context),
        getMassErrorModel(search_context, ms_file_idx),#getMassErrorModel(results),
        params,
        MS2CHROM()
    )...)
end

function get_matched_precursors(
    spectra::MassSpecData,
    psms::DataFrame,
    results::ParameterTuningSearchResults,
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64
) where {P<:SearchParameters}
    
    return vcat(mass_error_search(
        spectra,
        psms[!,:scan_idx],
        psms[!,:precursor_idx],
        UInt32(ms_file_idx),
        getSpecLib(search_context),
        getSearchData(search_context),
        getMs1MassErrorModel(results),
        params,
        MS1CHROM()
    )...)
end

#==========================================================
Plotting Helpers
==========================================================#

"""
    generate_rt_plot(results::ParameterTuningSearchResults, plot_path::String, title::String)

Generates retention time alignment visualization plot.

# Arguments
- `results`: Parameter tuning results containing RT data
- `plot_path`: Path to save plot
- `title`: Plot title

Creates scatter plot of RT vs iRT with fitted spline curve.
"""
function generate_rt_plot(
    results::ParameterTuningSearchResults,
    title::String
)
    n = length(results.rt)
    p = plot(
        results.rt,
        results.irt,
        seriestype=:scatter,
        title = title*"\n n = $n",
        xlabel = "Retention Time RT (min)",
        ylabel = "Indexed Retention Time iRT (min)",
        label = nothing,
        alpha = 0.1,
        size = 100*[13.3, 7.5]
    )
    
    pbins = LinRange(minimum(results.rt), maximum(results.rt), 100)
    plot!(pbins, getRtToIrtModel(results).(pbins), lw=3, label=nothing)
    return p
end


"""
    generate_mass_error_plot(results::ParameterTuningSearchResults, fname::String, plot_path::String)

Generates mass error distribution visualization plot.

# Arguments
- `results`: Parameter tuning results containing mass error data
- `fname`: Filename for plot title
- `plot_path`: Path to save plot

Creates histogram of mass errors with marked boundaries and offset.
"""
function generate_mass_error_plot(
    results::ParameterTuningSearchResults,
    fname::String
)
    #p = histogram(results.ppm_errs)
    #savefig(p, plot_path)
    mem = results.mass_err_model[]
    errs = results.ppm_errs .+ getMassOffset(mem)
    n = length(errs)
    plot_title = fname
    mass_err = getMassOffset(mem)
    bins = LinRange(mass_err - 2*getLeftTol(mem), mass_err + 2*getRightTol(mem), 50)
    p = Plots.histogram(errs,
    orientation = :h, 
    yflip = true,
    #seriestype=:scatter,
    title = plot_title*"\n n = $n",
    xlabel = "Count",
    ylabel = "Mass Error (ppm)",
    label = nothing,
    bins = bins,
    ylim = (mass_err - 2*getLeftTol(mem), mass_err + 2*getRightTol(mem)),
    topmargin =15mm,
    #bottommargin = 10mm,
    )

    mass_err = getMassOffset(mem)
    Plots.hline!([mass_err], label = nothing, color = :black, lw = 2)
    Plots.annotate!(last(xlims(p)), mass_err, text("$mass_err", :black, :right, :bottom, 12))


    l_err = mass_err - getLeftTol(mem)
    Plots.hline!([l_err], label = nothing, color = :black, lw = 2)
    Plots.annotate!(last(xlims(p)), l_err, text("$l_err", :black, :right, :bottom, 12))
    r_err = mass_err + getRightTol(mem)
    Plots.hline!([r_err], label = nothing, color = :black, lw = 2)
    Plots.annotate!(last(xlims(p)), r_err, text("$r_err", :black, :right, :bottom, 12))

    return p

end

"""
    generate_fallback_rt_plot_in_memory(results, fname, search_context, ms_file_idx)

Generates a diagnostic RT plot when no data is available (fallback/borrowed parameters).
Returns the plot object without saving to disk.
"""
function generate_fallback_rt_plot_in_memory(
    results::ParameterTuningSearchResults,
    fname::String,
    search_context::SearchContext,
    ms_file_idx::Int64
)
    # Create an informative plot showing the fallback/borrowed status
    p = Plots.plot(
        title = fname * "\n⚠️ Using Fallback/Borrowed Parameters",
        xlabel = "Retention Time RT (min)",
        ylabel = "Indexed Retention Time iRT (min)",
        size = 100*[13.3, 7.5],
        grid = true
    )
    
    # Add identity line if using identity model
    rt_model = getRtToIrtModel(results)
    if isa(rt_model, IdentityModel)
        rt_range = [0, 120]  # Default RT range for visualization
        Plots.plot!(rt_range, rt_range, 
                   lw=2, ls=:dash, color=:red,
                   label="Identity Model (Fallback)")
    else
        # If borrowed, show the borrowed model
        rt_range = LinRange(0, 120, 100)
        Plots.plot!(rt_range, rt_model.(rt_range),
                   lw=2, color=:blue,
                   label="Borrowed RT Model")
    end
    
    # Add text annotation explaining the situation
    Plots.annotate!(60, 20, 
                   text("Insufficient PSMs for RT calibration\n" *
                        "Using conservative parameters", 
                        :center, 10))
    
    return p  # Return plot for collection
end

"""
    generate_fallback_mass_error_plot_in_memory(results, fname, search_context, ms_file_idx)

Generates a diagnostic mass error plot when no data is available (fallback/borrowed parameters).
Returns the plot object without saving to disk.
"""
function generate_fallback_mass_error_plot_in_memory(
    results::ParameterTuningSearchResults,
    fname::String,
    search_context::SearchContext,
    ms_file_idx::Int64
)
    mem = getMassErrorModel(results)
    mass_offset = getMassOffset(mem)
    left_tol = getLeftTol(mem)
    right_tol = getRightTol(mem)
    
    # Create an informative plot showing the tolerances being used
    p = Plots.plot(
        title = fname * "\n⚠️ Using Fallback/Borrowed Parameters",
        xlabel = "Count",
        ylabel = "Mass Error (ppm)",
        size = 100*[13.3, 7.5],
        ylim = (-max(left_tol, right_tol) - 10, max(left_tol, right_tol) + 10),
        grid = true,
        yflip = true,
        orientation = :h
    )
    
    # Draw tolerance boundaries
    Plots.hline!([mass_offset], color=:black, lw=2, label="Mass Offset: $(round(mass_offset, digits=2)) ppm")
    Plots.hline!([mass_offset - left_tol], color=:red, lw=1, ls=:dash, label="Lower Bound: -$(round(left_tol, digits=1)) ppm")
    Plots.hline!([mass_offset + right_tol], color=:red, lw=1, ls=:dash, label="Upper Bound: +$(round(right_tol, digits=1)) ppm")
    
    # Add text annotation
    Plots.annotate!(0, 0,
                   text("Insufficient fragment matches\n" *
                        "Tolerance: ±$(round((left_tol+right_tol)/2, digits=1)) ppm",
                        :center, 10))
    
    return p  # Return plot for collection
end

# Legacy functions that save to disk (kept for backward compatibility)
function generate_fallback_rt_plot(
    results::ParameterTuningSearchResults,
    plot_path::String,
    fname::String,
    search_context::SearchContext,
    ms_file_idx::Int64
)
    p = generate_fallback_rt_plot_in_memory(results, fname, search_context, ms_file_idx)
    savefig(p, plot_path)
    return p
end

function generate_fallback_mass_error_plot(
    results::ParameterTuningSearchResults,
    plot_path::String,
    fname::String,
    search_context::SearchContext,
    ms_file_idx::Int64
)
    p = generate_fallback_mass_error_plot_in_memory(results, fname, search_context, ms_file_idx)
    savefig(p, plot_path)
    return p
end

function generate_ms1_mass_error_plot(
    results::R,
    fname::String
) where {R<:SearchResults}
    #p = histogram(results.ppm_errs)
    #savefig(p, plot_path)
    mem = results.ms1_mass_err_model[]
    errs = results.ms1_ppm_errs .+ getMassOffset(mem)
    n = length(errs)
    plot_title = fname
    mass_err = getMassOffset(mem)
    bins = LinRange(mass_err - 2*getLeftTol(mem), mass_err + 2*getRightTol(mem), 50)
    p = Plots.histogram(errs,
    orientation = :h, 
    yflip = true,
    #seriestype=:scatter,
    title = plot_title*"\n n = $n",
    xlabel = "Count",
    ylabel = "Mass Error (ppm)",
    label = nothing,
    bins = bins,
    ylim = (mass_err - 2*getLeftTol(mem), mass_err + 2*getRightTol(mem)),
    topmargin =15mm,
    #bottommargin = 10mm,
    )

    mass_err = getMassOffset(mem)
    Plots.hline!([mass_err], label = nothing, color = :black, lw = 2)
    Plots.annotate!(last(xlims(p)), mass_err, text("$mass_err", :black, :right, :bottom, 12))


    l_err = mass_err - getLeftTol(mem)
    Plots.hline!([l_err], label = nothing, color = :black, lw = 2)
    Plots.annotate!(last(xlims(p)), l_err, text("$l_err", :black, :right, :bottom, 12))
    r_err = mass_err + getRightTol(mem)
    Plots.hline!([r_err], label = nothing, color = :black, lw = 2)
    Plots.annotate!(last(xlims(p)), r_err, text("$r_err", :black, :right, :bottom, 12))

    return p

end
#==========================================================
Utility Functions
==========================================================#

"""
    cap_mass_tolerance(tol::Float32, max_tol::Float32 = 50.0f0)

Caps a mass tolerance value to a maximum limit.

# Arguments
- `tol`: Tolerance value to cap
- `max_tol`: Maximum allowed tolerance (default: 50.0 ppm)

# Returns
- Capped tolerance value
"""
function cap_mass_tolerance(tol::Float32, max_tol::Float32 = 50.0f0)
    return min(tol, max_tol)
end

"""
    create_capped_mass_model(offset, left_tol, right_tol, max_tol)

Creates a MassErrorModel with tolerance values capped to a maximum.

# Arguments
- `offset`: Mass offset in ppm
- `left_tol`: Left tolerance value
- `right_tol`: Right tolerance value  
- `max_tol`: Maximum allowed tolerance (default: 50.0 ppm)

# Returns
- MassErrorModel with capped tolerances
"""
function create_capped_mass_model(
    offset::Float32, 
    left_tol::Float32, 
    right_tol::Float32, 
    max_tol::Float32 = 50.0f0
)
    return MassErrorModel(
        offset,
        (cap_mass_tolerance(left_tol, max_tol), 
         cap_mass_tolerance(right_tol, max_tol))
    )
end

function check_convergence(
    psms,
    new_mass_err_model::MassErrorModel,
    old_mass_err_model::MassErrorModel,
    ppms,
    min_psms::Int64
)
    if size(psms, 1) < min_psms
        @info "size(psms, 1) < $min_psms, skipping convergence check"
        return false 
    end
    if (getLeftTol(old_mass_err_model) - getLeftTol(new_mass_err_model))/(getLeftTol(new_mass_err_model)) < 0.05
        @info "Left tolerance failed to converge getLeftTol(old_mass_err_model) $(getLeftTol(old_mass_err_model)) getLeftTol(new_mass_err_model) $(getLeftTol(new_mass_err_model)) describe(ppms) $(describe(ppms))"
        return false 
    end
    if (getRightTol(old_mass_err_model) - getRightTol(new_mass_err_model))/(getRightTol(new_mass_err_model)) < 0.05
        @info "Right tolerance failed to converge getRightTol(old_mass_err_model) $(getRightTol(old_mass_err_model)) getRightTol(new_mass_err_model) $(getRightTol(new_mass_err_model)) describe(ppms) $(describe(ppms))"
        return false 
    end
    return true
end

"""
    test_tolerance_expansion!(search_context, params, ms_file_idx,
                             current_psms, current_model, current_ppm_errs,
                             collection_tolerance, filtered_spectra, spectra)

After convergence, test if expanding the collection tolerance yields more PSMs.
If successful, returns the expanded PSM set and refitted model.

# Returns
- `(final_psms, final_model, final_ppm_errs, was_expanded::Bool)`
"""
function test_tolerance_expansion!(
    search_context::SearchContext,
    params::ParameterTuningSearchParameters,
    ms_file_idx::Int64,
    current_psms::DataFrame,
    current_model::MassErrorModel,
    current_ppm_errs::Vector{<:AbstractFloat},
    collection_tolerance::Float32,
    filtered_spectra::FilteredMassSpecData,
    spectra::MassSpecData
)
    # Hardcoded expansion factor
    expansion_factor = 1.5f0
    
    current_psm_count = size(current_psms, 1)
    
    # Calculate expanded tolerance
    expanded_tolerance = collection_tolerance * expansion_factor
    
    # Create expanded model for collection
    # Keep the same bias, just expand the window
    expanded_model = MassErrorModel(
        getMassOffset(current_model),
        (expanded_tolerance, expanded_tolerance)
    )
    
     @info "Testing tolerance expansion for PSM collection:" *
           "\n  Original collection tolerance: ±$(round(collection_tolerance, digits=1)) ppm" *
           "\n  Expanded collection tolerance: ±$(round(expanded_tolerance, digits=1)) ppm" *
           "\n  Current PSM count: $current_psm_count"
    
           # Store original model
    original_model = getMassErrorModel(search_context, ms_file_idx)
    
    # Set expanded model for collection
    setMassErrorModel!(search_context, ms_file_idx, expanded_model)
    
    # Collect PSMs with expanded tolerance
    expanded_psms = collect_psms(filtered_spectra, spectra, search_context, params, ms_file_idx)
    expanded_psm_count = size(expanded_psms, 1)
    
    # Calculate improvement
    psm_increase = expanded_psm_count - current_psm_count
    improvement_ratio = current_psm_count > 0 ? psm_increase / current_psm_count : 0.0
    
    # @info "Expansion results:" *
    #       "\n  Expanded PSM count: $expanded_psm_count" *
    #       "\n  PSM increase: $psm_increase ($(round(100*improvement_ratio, digits=1))%)"
    
    # Check if expansion was beneficial (any improvement)
    if expanded_psm_count <= current_psm_count
        # No improvement, restore original and return
        setMassErrorModel!(search_context, ms_file_idx, original_model)
        # @info "No improvement found, keeping original results"
        return current_psms, current_model, current_ppm_errs, false
    end
    
    # Significant improvement found - refit model with expanded PSM set
    # @info "Improvement found, refitting model with expanded PSM set"
    
    # Get matched fragments for the expanded PSM set
    fragments = get_matched_fragments(spectra, expanded_psms, search_context, params, ms_file_idx)
    
    if length(fragments) == 0
        # Failed to get fragments, restore original
        setMassErrorModel!(search_context, ms_file_idx, original_model)
        @warn "Failed to extract fragments from expanded PSM set, keeping original"
        return current_psms, current_model, current_ppm_errs, false
    end
    
    # Fit new mass error model from expanded PSM set
    refitted_model, refitted_ppm_errs = fit_mass_err_model(params, fragments)
    
    if refitted_model === nothing
        # Failed to fit model, restore original
        setMassErrorModel!(search_context, ms_file_idx, original_model)
        @warn "Failed to fit mass error model from expanded PSM set, keeping original"
        return current_psms, current_model, current_ppm_errs, false
    end
    
    # Success! Use the expanded results
    # @info "Successfully expanded tolerance:" *
    #       "\n  Original fitted: ±$(round((getLeftTol(current_model) + getRightTol(current_model))/2, digits=1)) ppm" *
    #       "\n  Expanded fitted: ±$(round((getLeftTol(refitted_model) + getRightTol(refitted_model))/2, digits=1)) ppm" *
    #       "\n  PSM improvement: $psm_increase PSMs ($(round(100*improvement_ratio, digits=1))%)"
    
    # Update the model in search context
    setMassErrorModel!(search_context, ms_file_idx, refitted_model)
    
    return expanded_psms, refitted_model, refitted_ppm_errs, true
end

"""
    get_fallback_parameters(search_context, params, ms_file_idx, warnings)

Get fallback parameters from another file or use defaults.
Returns (mass_err_model, rt_model, borrowed_from_idx)
"""
function get_fallback_parameters(
    search_context::SearchContext,
    params::ParameterTuningSearchParameters,
    ms_file_idx::Int64,
    warnings::Vector{String}
)
    # Try to borrow from another successful file
    n_files = length(getMSData(search_context).file_paths)
    borrowed_from = nothing
    
    for other_idx in 1:n_files
        if other_idx != ms_file_idx && 
           haskey(search_context.mass_error_model, other_idx) &&
           haskey(search_context.rt_conversion_model, other_idx)
            mass_err = getMassErrorModel(search_context, other_idx)
            rt_model = getRtConversionModel(search_context, other_idx)
            push!(warnings, "Borrowed parameters from file $other_idx")
            return mass_err, rt_model, other_idx
        end
    end
    
    # Use default fallback
    @warn "No successful files to borrow from, using default parameters"
    push!(warnings, "Using default fallback parameters")
    
    # Default mass error model
    mass_err = MassErrorModel(0.0f0, (20.0f0, 20.0f0))
    
    # Default RT model (identity)
    rt_model = LinearRtConversionModel(1.0f0, 0.0f0)
    
    return mass_err, rt_model, nothing
end

"""
    record_file_status!(diagnostics, ms_file_idx, file_name, converged, used_fallback, warnings, iteration_state)

Record the parameter tuning status for a file.
"""
function record_file_status!(
    diagnostics::ParameterTuningDiagnostics,
    ms_file_idx::Int64,
    file_name::String,
    converged::Bool,
    used_fallback::Bool,
    warnings::Vector{String},
    iteration_state::IterationState
)
    status = ParameterTuningStatus(
        ms_file_idx,
        file_name,
        converged,
        used_fallback,
        used_fallback ? "Failed convergence" : "",
        iteration_state.total_iterations,
        0,  # PSM count - could be tracked if needed
        0.0f0,  # Mass offset - could be extracted if needed
        (0.0f0, 0.0f0),  # Tolerances - could be extracted if needed
        warnings
    )
    
    diagnostics.file_statuses[ms_file_idx] = status
    
    if converged && !used_fallback
        diagnostics.n_successful += 1
    elseif used_fallback
        diagnostics.n_fallback += 1
    else
        diagnostics.n_failed += 1
    end
end


