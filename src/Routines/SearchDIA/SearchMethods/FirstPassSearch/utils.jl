
"""
    add_main_search_columns!(psms::DataFrame, rt_irt::UniformSpline,
                           structural_mods::AbstractVector{Union{Missing, String}},
                           prec_missed_cleavages::Arrow.Primitive{UInt8, Vector{UInt8}},
                           prec_is_decoy::Arrow.BoolVector{Bool},
                           prec_irt::Arrow.Primitive{T, Vector{T}},
                           prec_charge::Arrow.Primitive{UInt8, Vector{UInt8}},
                           scan_retention_time::AbstractVector{Float32},
                           tic::AbstractVector{Float32},
                           masses::AbstractArray) where {T<:AbstractFloat}

Adds essential columns to PSM DataFrame for scoring and analysis.

# Arguments
- `psms`: DataFrame to modify
- `rt_irt`: Spline for RT to iRT conversion
- `structural_mods`: Vector of structural modifications
- `prec_missed_cleavages`: Vector of missed cleavage counts
- `prec_is_decoy`: Vector indicating decoy status
- `prec_irt`: Vector of iRT values
- `prec_charge`: Vector of precursor charges
- `scan_retention_time`: Vector of scan retention times
- `tic`: Vector of total ion currents
- `masses`: Array of mass spectra

# Added Columns
- Basic metrics: RT, iRT, charge, TIC
- Modification info: missed_cleavage, Mox
- Scoring columns: score, q_value
- Analysis columns: target, spectrum_peak_count, err_norm
"""
function add_main_search_columns!(psms::DataFrame, 
                                rt_irt::UniformSpline,
                                structural_mods::AbstractVector{Union{Missing, String}},
                                prec_missed_cleavages::Arrow.Primitive{UInt8, Vector{UInt8}},
                                prec_is_decoy::Arrow.BoolVector{Bool},
                                prec_irt::Arrow.Primitive{T, Vector{T}},
                                prec_charge::Arrow.Primitive{UInt8, Vector{UInt8}},
                                scan_retention_time::AbstractVector{Float32},
                                tic::AbstractVector{Float32},
                                masses::AbstractArray;
) where {T<:AbstractFloat}
    
    ###########################
    #Allocate new columns
    N = size(psms, 1)
    missed_cleavage = zeros(UInt8, N);
    Mox = zeros(UInt8, N);
    irt_pred = zeros(Float32, N);
    rt = zeros(Float32, N);
    irt = zeros(Float32, N);
    TIC = zeros(Float16, N);
    charge = zeros(UInt8, N);
    err_norm = zeros(Float16, N);
    targets = zeros(Bool, N);
    spectrum_peak_count = zeros(UInt32, N);
    scan_idx::Vector{UInt32} = psms[!,:scan_idx]
    precursor_idx::Vector{UInt32} = psms[!,:precursor_idx] 
    error::Vector{Float32} = psms[!,:error]
    total_ions::Vector{UInt8} = psms[!,:y_count] .+ psms[!,:b_count]
    function countMOX(seq::String)
        return UInt8(count("Unimod:35", seq))
    end

    #Split data into chunk ranges
    tasks_per_thread = 10
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                prec_idx = precursor_idx[i]

                targets[i] = prec_is_decoy[prec_idx] == false;
                missed_cleavage[i] = prec_missed_cleavages[prec_idx]
                Mox[i] = countMOX(coalesce(structural_mods[prec_idx], ""))::UInt8 #UInt8(length(collect(eachmatch(r"ox",  precursors[precursor_idx[i]].sequence))))
                irt_pred[i] = Float32(prec_irt[prec_idx]);
                rt[i] = Float32(scan_retention_time[scan_idx[i]]);
                irt[i] = rt_irt(rt[i])
                TIC[i] = Float16(log2(tic[scan_idx[i]]));
                charge[i] = UInt8(prec_charge[prec_idx]);
                err_norm[i] = min(Float16.(abs(error[i]/total_ions[i])), Float16(6e4))#Float16(min(abs((error[i])/(total_ions[i])), 6e4))
                spectrum_peak_count[i] = UInt32(length(masses[scan_idx[i]]))
            end
        end
    end
    fetch.(tasks)
    psms[!,:irt_predicted] = irt_pred
    psms[!,:rt] = rt
    psms[!,:irt] = irt
    psms[!,:TIC] = TIC
    psms[!,:total_ions] = total_ions
    psms[!,:err_norm] = err_norm
    psms[!,:target] = targets
    psms[!,:missed_cleavage] = missed_cleavage
    psms[!,:Mox] = Mox
    psms[!,:charge] = charge
    psms[!,:spectrum_peak_count] = spectrum_peak_count
    psms[!,:score] = zeros(Float32, N);
    psms[!,:q_value] = zeros(Float16, N);
    psms[!,:intercept] = ones(Float16, N)
end

"""
    score_main_search_psms!(psms::DataFrame, column_names::Vector{Symbol};
                           n_train_rounds::Int64=2,
                           max_iter_per_round::Int64=50,
                           max_q_value::Float64=0.01)

Performs iterative probit regression scoring of PSMs.

# Arguments
- `psms`: DataFrame containing PSMs to score
- `column_names`: Column names to use for scoring
- `n_train_rounds`: Number of training iterations
- `max_iter_per_round`: Maximum iterations per training round
- `max_q_value`: Maximum q-value threshold for filtering

# Process
1. First round: Trains on all data
2. Subsequent rounds: Trains on decoys and high-scoring targets
3. Updates scores and q-values after each round
"""
function score_main_search_psms!(psms::DataFrame, column_names::Vector{Symbol};
                             n_train_rounds::Int64 = 2,
                             max_iter_per_round::Int64 = 50,
                             max_q_value::Float64 = 0.01
)

    β = zeros(Float64, length(column_names));
    best_psms = nothing#ones(Bool, size(psms, 1))

    tasks_per_thread = 10
    M = size(psms, 1)
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:M, chunk_size) # partition your data into chunks

    for i in range(1, n_train_rounds)
        if i < 2 #Train on top scribe data during first round
            get_qvalues!(psms[!,:scribe], psms[!,:target], psms[!,:q_value])
        end

        if i < n_train_rounds #Get Data to train on
            best_psms = ((psms[!,:q_value].<=max_q_value).&(psms[!,:target])) .| (psms[!,:target].==false);
        end

        psms_targets = psms[best_psms,:target]
        M = size(psms_targets, 1)
        sub_chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
        sub_data_chunks = partition(1:M, sub_chunk_size) # partition your data into chunks that
        β = ProbitRegression(β, psms[best_psms,column_names], psms_targets, sub_data_chunks, max_iter = max_iter_per_round);
        ModelPredict!(psms[!,:score], psms[!,column_names], β, data_chunks); #Get Z-scores 
    end

    return 
end

"""
    get_probs!(psms::DataFrame, psms_scores::Vector{Float32})

Calculates probability scores for PSMs based on their Z-scores using error function.

# Arguments
- `psms`: DataFrame to modify
- `psms_scores`: Vector of PSM Z-scores

Adds 'prob' column to psms DataFrame containing probability scores.
Uses parallel processing for efficiency.
"""
function get_probs!(psms::DataFrame,
                    psms_scores::Vector{Float32}
)

    psms_probs = zeros(Float32, size(psms, 1))
    tasks_per_thread = 10
    M = size(psms, 1)
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:M, chunk_size) # partition your data into chunks that
    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                psms_probs[i] = (1 + SpecialFunctions.erf(psms_scores[i]/sqrt(2)))/2
            end
        end
    end
    fetch.(tasks)
    psms[!,:prob] = psms_probs
    return nothing
end


"""
    get_best_psms!(psms::DataFrame, prec_mz::Arrow.Primitive{T, Vector{T}}; 
                   max_q_val::Float32=0.01f0, max_psms::Int64=250000) where {T<:AbstractFloat}

Processes PSMs to identify the best matches and calculate peak characteristics.

# Arguments
- `psms`: DataFrame containing peptide-spectrum matches
- `prec_mz`: Vector of precursor m/z values
- `max_q_val`: Maximum q-value threshold for filtering PSMs
- `max_psms`: Maximum number of PSMs to retain

# Modifies PSMs DataFrame to add:
- `best_psm`: Boolean indicating highest scoring PSM for each precursor
- `fwhm`: Full width at half maximum of chromatographic peak
- `scan_count`: Number of scans below q-value threshold
- `prec_mz`: Precursor m/z values

Note: Assumes psms is sorted by retention time in ascending order.
"""
function get_best_psms!(psms::DataFrame,
                        prec_mz::Arrow.Primitive{T, Vector{T}};
                        max_local_fdr::Float32 = 1.00f0
) where {T<:AbstractFloat}

    #highest scoring psm for a given precursor
    psms[!,:local_fdr] = zeros(Float16, size(psms, 1))
    #highest scoring psm for a given precursor
    psms[!,:best_psm] = zeros(Bool, size(psms, 1))
    #fwhm estimate of the precursor
    psms[!,:fwhm] = zeros(Union{Missing, Float32}, size(psms, 1))
    #number of scans below the q value threshold for hte precursor
    psms[!,:scan_count] = zeros(UInt16,size(psms, 1))

    # Will use local FDR for final filter
    get_local_FDR!(psms[!,:score], psms[!,:target], psms[!,:local_fdr]);
    #Get best psm for each precursor 
    #ASSUMES psms IS SOrtED BY rt IN ASCENDING ORDER
    gpsms = groupby(psms,:precursor_idx)
    for (precursor_idx, prec_psms) in pairs(gpsms)

        #Get the best scoring psm, and the its row index. 
        #Get the maximum intensity psm (under the q_value threshold) and its row index. 
        max_irt, min_irt = missing, missing
        max_log2_intensity, max_idx = missing, one(Int64)
        best_psm_score, best_psm_idx = zero(Float32), one(Int64)
        scan_count = one(UInt8)
        for i in range(1, size(prec_psms, 1))
            if coalesce(max_log2_intensity, zero(Float32)) < prec_psms[i,:log2_summed_intensity]
                max_log2_intensity = prec_psms[i,:log2_summed_intensity]
                max_idx = i
            end
            if prec_psms[i,:score]>best_psm_score
                best_psm_idx = i
                best_psm_score = prec_psms[i,:score]
            end
        end
        #Mark the best psm 
        prec_psms[best_psm_idx,:best_psm] = true

        #Try to estimate the fwhm. 
        i = max_idx - 1
        while i > 0
            #Is the i'th psm above half the maximum 
            if (prec_psms[i,:log2_summed_intensity] > (max_log2_intensity - 1.0))
                scan_count += 1
                min_irt = prec_psms[i,:irt]
            else
                break
            end
            i -= 1
        end
        
        i = max_idx + 1
        while i <= size(prec_psms, 1)
            #Is the i'th psm above half the maximum 
            if (prec_psms[i,:log2_summed_intensity] > (max_log2_intensity - 1.0))
                scan_count += 1
                max_irt = prec_psms[i,:irt]
            else
                break
            end
            i += 1
        end

        prec_psms[best_psm_idx,:fwhm] = max_irt - min_irt
        prec_psms[best_psm_idx,:scan_count] = scan_count
    end

    filter!(x->x.best_psm, psms);
    sort!(psms,:score, rev = true)
    n = size(psms, 1)
    select!(psms, [:precursor_idx,:rt,:irt_predicted,:q_value,:score,:prob,:fwhm,:scan_count,:scan_idx,:local_fdr])
    #Instead of max_psms, 2x the number at 10% fdr. 
    psms_passing = 0
    local_fdrs = psms[!,:local_fdr]::Vector{Float16}
    @inbounds for i in range(1, n)
        psms_passing += 1
        if local_fdrs[i]>=max_local_fdr
            break
        end
    end
    #delete!(psms, min(n, max_psms + 1):n)
    deleteat!(psms, min(n, round(Int64, psms_passing) + 1):n)

    mz = zeros(T, size(psms, 1));
    precursor_idx = psms[!,:precursor_idx]::Vector{UInt32}
    Threads.@threads for i in range(1, size(psms, 1))
        mz[i] = prec_mz[precursor_idx[i]];
    end
    psms[!,:prec_mz] = mz

    return
end

"""
    map_retention_times!(search_context::SearchContext, results::FirstPassSearchResults, params::FirstPassSearchParameters)

Maps retention times between library and empirical scales for each MS file. For files with sufficient high-confidence PSMs 
(probability > min_prob_for_irt_mapping), creates spline-based conversion models between retention time (RT) and indexed 
retention time (iRT) scales.

# Arguments
- `search_context`: Contains MS data and mapping information
- `results`: FirstPassSearch results container
- `params`: Search parameters including minimum probability threshold

Creates and stores:
- RT to iRT mapping splines
- iRT to RT mapping splines

Throws an error if insufficient high-confidence PSMs are found for mapping.
"""
function map_retention_times!(
    search_context::SearchContext,
    results::FirstPassSearchResults,
    params::FirstPassSearchParameters
)
    @info "Mapping library to empirical retention times..."
    
    for (ms_file_idx, psms_path) in enumerate(getFirstPassPsms(getMSData(search_context)))
        #if getFailedIndicator(getMSData(search_context), ms_file_idx)==true
        #    continue
        #end
        psms = Arrow.Table(psms_path)
        best_hits = psms[:prob].>params.min_prob_for_irt_mapping#Map rts using only the best psms
        try#if sum(best_hits) > 100
            best_rts = psms[:rt][best_hits]
            best_irts = psms[:irt_predicted][best_hits]
            irt_to_rt_spline = UniformSpline(
                                        best_rts,
                                        best_irts,
                                        3, 
                                        5
            )
            rt_to_irt_spline = UniformSpline(
                best_irts,
                best_rts,
                3, 
                5
            )
            #Build rt=>irt and irt=> rt mappings for the file and add to the dictionaries 
            setRtIrtMap!(search_context, SplineRtConversionModel(rt_to_irt_spline), ms_file_idx)
            setIrtRtMap!(search_context, SplineRtConversionModel(irt_to_rt_spline), ms_file_idx)
        catch
            throw("add a default option here...")
            #sensible default here?
            continue
        end
    end
    return nothing
end


PrecToIrtType = Dictionary{UInt32, 
    NamedTuple{
        (:best_prob, :best_ms_file_idx, :best_scan_idx, :best_irt, :mean_irt, :var_irt, :n, :mz), 
        Tuple{Float32, UInt32, UInt32, Float32, Union{Missing, Float32}, Union{Missing, Float32}, Union{Missing, UInt16}, Float32}
    }
}

"""
    create_rt_indices!(search_context::SearchContext, results::FirstPassSearchResults,
                      precursor_dict::PrecToIrtType, params::FirstPassSearchParameters)

Creates retention time indices for improved search efficiency.

# Arguments
- `search_context`: Search context containing MS data
- `results`: FirstPassSearch results
- `precursor_dict`: Dictionary mapping precursors to iRT data
- `params`: Search parameters

# Process
1. Calculates iRT errors using peak width and variation
2. Creates precursor to iRT mapping
3. Generates RT indices for each MS file
4. Stores indices in search context
"""
function create_rt_indices!(
    search_context::SearchContext,
    results::FirstPassSearchResults,
    precursor_dict::PrecToIrtType,
    params::FirstPassSearchParameters
)


    # Calculate iRT errors
    @info "Calculating iRT errors..."
    irt_errs = get_irt_errs(results.fwhms, precursor_dict, params)


    setIrtErrors!(search_context, irt_errs)

    @info "Creating RT indices..."
    # Create precursor to iRT mapping
    prec_to_irt = map(x -> (irt=x[:best_irt], mz=x[:mz]), 
                      precursor_dict)

    # Set up indices folder
    rt_indices_folder = joinpath(getDataOutDir(search_context), "temp_data", "rt_indices")
    !isdir(rt_indices_folder) && mkdir(rt_indices_folder)

    # Make RT indices
    rt_index_paths = makeRTIndices(
        rt_indices_folder,
        getFirstPassPsms(getMSData(search_context)),
        prec_to_irt,
        getRtIrtMap(search_context),
        min_prob=params.max_prob_to_impute
    )
    for (ms_file_idx, rt_index_path) in pairs(rt_index_paths)
        setRtIndex!(getMSData(search_context), ms_file_idx, rt_index_path)
    end
end


"""
    get_irt_errs(fwhms::Dictionary, prec_to_irt::Dictionary, params::FirstPassSearchParameters)

Calculates iRT error tolerances based on peak widths and cross-run variation.

# Arguments
- `fwhms`: Dictionary of FWHM statistics per file
- `prec_to_irt`: Dictionary of precursor iRT data
- `params`: Parameters including FWHM and iRT standard deviation multipliers

# Returns
Dictionary mapping file indices to iRT tolerances, combining:
- Peak width variation (FWHM + n*MAD)
- Cross-run iRT variation
"""
function get_irt_errs(
    fwhms::Dictionary{Int64, 
                        @NamedTuple{
                            median_fwhm::Float32,
                            mad_fwhm::Float32
                        }},
    prec_to_irt::Dictionary{UInt32, 
    @NamedTuple{best_prob::Float32, 
                best_ms_file_idx::UInt32, 
                best_scan_idx::UInt32, 
                best_irt::Float32, 
                mean_irt::Union{Missing, Float32}, 
                var_irt::Union{Missing, Float32}, 
                n::Union{Missing, UInt16}, 
                mz::Float32}}
    ,
    params::FirstPassSearchParameters
)
    #Get upper bound on peak fwhm. Use median + n*standard_deviation
    #estimate standard deviation by the median absolute deviation. 
    #n is a user-defined paramter. 
    fwhms = map(x->x[:median_fwhm] + params.fwhm_nstd*x[:mad_fwhm],
    fwhms)
    #Get variance in irt of apex accross runs. Only consider precursor identified below q-value threshold
    #in more than two runs .
    irt_std = nothing
    variance_  = collect(skipmissing(map(x-> (x[:n] > 2) ? sqrt(x[:var_irt]/(x[:n] - 1)) : missing, prec_to_irt)))
    if !iszero(length(variance_))
        irt_std = median(variance_)
    else
        #This could happen if only two files are being searched 
        variance_  = collect(skipmissing(map(x-> (x[:n] == 2) ? sqrt(x[:var_irt]) : missing, prec_to_irt)))
        if iszero(length(variance_)) #only searching one file so 
            irt_std = 0.0f0
        else
            irt_std = median(variance_)
        end
    end
    #Number of standard deviations to cover 
    irt_std *= params.irt_nstd
    #dictionary maping file name to irt tolerance. 
    return map(x->Float32((x+irt_std))::Float32, fwhms)::Dictionary{Int64, Float32}
end




