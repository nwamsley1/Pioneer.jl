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
- Number of passing PSMs if above minimum threshold
- -1 if insufficient PSMs pass threshold

Modifies psms DataFrame in place by filtering to best matches.
"""
function filter_and_score_psms!(
    psms::DataFrame,
    params::P
) where {P<:ParameterTuningSearchParameters}
    
    score_presearch!(psms)
    get_qvalues!(psms[!,:prob], psms[!,:target], psms[!,:q_value])
    
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
        return n_passing_psms
    end
    
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
    tasks_per_thread = 10
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:M, chunk_size) # partition your data into chunks that
    β = zeros(Float64, length(features))
    β = ProbitRegression(β, psms[!,features], psms[!,:target], data_chunks, max_iter = 20)
    ModelPredictProbs!(psms[!,:prob], psms[!,features], β, data_chunks)
end


"""
    get_qvalues!(probs::Vector{U}, labels::Vector{Bool}, qvals::Vector{T}) where {T,U<:AbstractFloat}
    get_qvalues!(PSMs::DataFrame, probs::Vector{Float64}, labels::Vector{Bool})

Calculates q-values (false discovery rate estimates) for PSMs.

# Arguments
- `probs`: Vector of probability scores
- `labels`: Vector of target/decoy labels
- `qvals`: Vector to store calculated q-values

# Process
1. Sorts PSMs by probability score
2. Calculates running ratio of decoys to targets
3. Assigns q-values based on decoy/target ratio

Implements target-decoy approach for FDR estimation.
"""
function get_qvalues!(probs::Vector{U}, labels::Vector{Bool}, qvals::Vector{T}
) where {T,U<:AbstractFloat}

    order = sortperm(probs, rev = true,alg=QuickSort) #Sort class probabilities
    targets = 0
    decoys = 0
    @inbounds @fastmath for i in order
            targets += labels[i]
            decoys += (1 - labels[i])
            qvals[i] = decoys/(targets)
    end
end
get_qvalues!(PSMs::DataFrame, probs::Vector{Float64}, labels::Vector{Bool}) = get_qvalues!(PSMs, allowmissing(probs), allowmissing(labels))


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

function collectFragErrs(fmatches::Vector{M}, new_fmatches::Vector{M}, nmatches::Int, n::Int) where {M<:MatchIon{Float32}}
    for match in range(1, nmatches)
        if n < length(fmatches)
            n += 1
            fmatches[n] = new_fmatches[match]
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
                ion_idx, _ = selectTransitions!(
                    getIonTemplates(search_data[thread_id]),
                    MassErrEstimationStrategy(),
                    getPrecEstimation(params),
                    getFragmentLookupTable(spec_lib),
                    scan_to_prec_idx[scan_idx],
                    precursor_idxs,
                    max_rank = 5#getMaxBestRank(params)
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
                    frag_err_idx
                )
            end
            
            @view(getMassErrMatches(search_data[thread_id])[1:frag_err_idx])
        end
    end
    fetch.(tasks)
end

function fit_mass_err_model(
    params::P,
    fragments::Vector{FragmentMatch{Float32}}
) where {P<:FragmentIndexSearchParameters}
    """
        calc_ppm_error(theoretical::T, observed::T) where T<:AbstractFloat

    Calculates parts-per-million (PPM) mass error.

    # Arguments
    - `theoretical`: Theoretical mass
    - `observed`: Observed mass

    Returns PPM error between theoretical and observed masses.
    """
    function calc_ppm_error(theoretical::T, observed::T) where T<:AbstractFloat
        return Float32((observed - theoretical)/(theoretical/1e6))
    end
    # Calculate PPM errors
    ppm_errs = [calc_ppm_error(match.theoretical_mz, match.match_mz) for match in fragments]
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
    results::ParameterTuningSearchResults,
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
        getMassErrorModel(results),
        params
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
    plot_path::String,
    title::String
)
    n = length(results.rt)
    p = Plots.plot(
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
    Plots.plot!(pbins, getRtToIrtModel(results).(pbins), lw=3, label=nothing)
    savefig(p, plot_path)
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
    fname::String,
    plot_path::String
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

    savefig(p, plot_path)

end

#==========================================================
Utility Functions
==========================================================#



