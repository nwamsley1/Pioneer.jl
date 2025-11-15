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
                                rt_irt::T,
                                structural_mods::AbstractVector{Union{Missing, String}},
                                prec_missed_cleavages::Arrow.Primitive{UInt8, Vector{UInt8}},
                                prec_is_decoy::Arrow.BoolVector{Bool},
                                prec_irt::Arrow.Primitive{U, Vector{U}},
                                prec_charge::Arrow.Primitive{UInt8, Vector{UInt8}},
                                scan_retention_time::AbstractVector{Float32},
                                tic::AbstractVector{Float32},
                                masses::AbstractArray;
) where {T,U<:AbstractFloat}
    
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
                           max_q_value::Float64=0.01,
                           fdr_scale_factor::Float32=1.0f0)

Performs iterative probit regression scoring of PSMs.

# Arguments
- `psms`: DataFrame containing PSMs to score
- `column_names`: Column names to use for scoring
- `n_train_rounds`: Number of training iterations
- `max_iter_per_round`: Maximum iterations per training round
- `max_q_value`: Maximum q-value threshold for filtering
- `fdr_scale_factor`: Scale factor to correct for library target/decoy ratio

# Process
1. First round: Trains on all data
2. Subsequent rounds: Trains on decoys and high-scoring targets
3. Updates scores and q-values after each round
"""
function score_main_search_psms!(psms::DataFrame, column_names::Vector{Symbol};
                             n_train_rounds::Int64 = 2,
                             max_iter_per_round::Int64 = 50,
                             max_q_value::Float64 = 0.01,
                             fdr_scale_factor::Float32 = 1.0f0
)

    β = zeros(Float64, length(column_names));
    best_psms = nothing#ones(Bool, size(psms, 1))

    tasks_per_thread = 10
    M = size(psms, 1)
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:M, chunk_size) # partition your data into chunks

    for i in range(1, n_train_rounds)
        if i < 2 #Train on top scribe data during first round
            get_qvalues!(psms[!,:scribe], psms[!,:target], psms[!,:q_value]; fdr_scale_factor = fdr_scale_factor)
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
                   max_PEP::Float32=0.9f0, fdr_scale_factor::Float32=1.0f0) where {T<:AbstractFloat}

Processes PSMs to identify the best matches and calculate peak characteristics.

# Arguments
- `psms`: DataFrame containing peptide-spectrum matches
- `prec_mz`: Vector of precursor m/z values
- `max_PEP`: Maximum local FDR threshold for filtering PSMs
- `fdr_scale_factor`: Scale factor to correct for library target/decoy ratio

# Modifies PSMs DataFrame to add:
- `best_psm`: Boolean indicating highest scoring PSM for each precursor
- `fwhm`: Full width at half maximum of chromatographic peak
- `scan_count`: Number of scans below q-value threshold
- `prec_mz`: Precursor m/z values

Note: Assumes psms is sorted by retention time in ascending order.
"""
function get_best_psms!(psms::DataFrame,
                        prec_mz::Arrow.Primitive{T, Vector{T}};
                        max_PEP::Float32 = 0.9f0,
                        fdr_scale_factor::Float32 = 1.0f0
) where {T<:AbstractFloat}

    #highest scoring psm for a given precursor
    psms[!,:PEP] = zeros(Float16, size(psms, 1))
    #highest scoring psm for a given precursor
    psms[!,:best_psm] = zeros(Bool, size(psms, 1))
    #fwhm estimate of the precursor
    psms[!,:fwhm] = zeros(Union{Missing, Float32}, size(psms, 1))
    #number of scans below the q value threshold for hte precursor
    psms[!,:scan_count] = zeros(UInt16,size(psms, 1))

    
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
    # Will use PEP for final filter
    get_PEP!(psms[!,:score], psms[!,:target], psms[!,:PEP]; doSort=false, fdr_scale_factor=fdr_scale_factor);

    n = size(psms, 1)
    select!(psms, [:precursor_idx,:log2_summed_intensity,:rt,:irt_predicted,:q_value,:score,:prob,:fwhm,:scan_count,:scan_idx,:PEP,:target])

    first_fail = searchsortedfirst(psms[!,:PEP], Float16(max_PEP))
    if first_fail <= n
        deleteat!(psms, first_fail:n)
    end
    #println("unique IDs prefilter: ", n, " ", first_fail, "\n\n")

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

    # Only process valid (non-failed) files
    valid_files = get_valid_file_indices(search_context)
    all_psms_paths = getFirstPassPsms(getMSData(search_context))

    for ms_file_idx in valid_files
        # Skip files that have been marked as failed
        if is_file_failed(search_context, ms_file_idx)
            continue
        end
        psms_path = all_psms_paths[ms_file_idx]
        psms = Arrow.Table(psms_path)
        best_hits = psms[:prob].>params.min_prob_for_irt_mapping#Map rts using only the best psms
        try#if sum(best_hits) > 100
            if params.use_robust_fitting
                # Use robust fitting with RANSAC/penalty
                best_psms_df = DataFrame(
                    rt = psms[:rt][best_hits],
                    irt_predicted = psms[:irt_predicted][best_hits]
                )

                # Fit PRELIMINARY RT to iRT model (for refinement training)
                # Will be replaced with refined model if refinement is successful
                rt_model_prelim, valid_rt_prelim, valid_irt_prelim, irt_mad = Pioneer.fit_irt_model(
                    best_psms_df;
                    lambda_penalty = Float32(0.1),  # Default penalty
                    ransac_threshold = 1000,        # Default RANSAC threshold
                    min_psms = 10,                  # Default minimum PSMs
                    spline_degree = 3,              # Default spline degree
                    max_knots = 7,                  # Default max knots
                    outlier_threshold = Float32(5.0)  # Default outlier threshold
                )

                # Don't store preliminary models yet if refinement is enabled
                # They will be replaced with refined models if refinement succeeds

                # Fit PRELIMINARY inverse model (iRT to RT)
                irt_to_rt_df_prelim = DataFrame(
                    rt = valid_irt_prelim,           # Use iRT as input
                    irt_predicted = valid_rt_prelim  # Use RT as output
                )
                irt_model_prelim, _, _, _ = Pioneer.fit_irt_model(
                    irt_to_rt_df_prelim;
                    lambda_penalty = Float32(0.1),
                    ransac_threshold = 1000,
                    min_psms = 10,
                    spline_degree = 3,
                    max_knots = 7,
                    outlier_threshold = Float32(5.0)
                )

                # iRT refinement (if enabled)
                if params.enable_irt_refinement
                    try
                        # Get sequences for high-quality PSMs
                        precursor_indices = psms[:precursor_idx][best_hits]
                        precursors = getPrecursors(getSpecLib(search_context))
                        sequences = String[getSequence(precursors)[idx] for idx in precursor_indices]

                        # Convert observed RTs to iRT using PRELIMINARY model
                        observed_rts = psms[:rt][best_hits]
                        observed_irts = Float32[rt_model_prelim(Float32(rt)) for rt in observed_rts]

                        # Step 1: Fit file-specific refinement model
                        refinement_model = fit_irt_refinement_model(
                            sequences,
                            Float32.(psms[:irt_predicted][best_hits]),  # library iRTs
                            observed_irts,                               # RT converted to iRT
                            ms_file_idx=ms_file_idx,
                            min_psms=20,
                            train_fraction=0.67
                        )

                        # Store model in SearchContext for SecondPass (via getPredIrt)
                        setIrtRefinementModel!(search_context, ms_file_idx, refinement_model)

                        # Step 2: Refit RT models with refined iRT (if refinement successful)
                        if !isnothing(refinement_model) && refinement_model.use_refinement
                            @user_info "File $ms_file_idx: Refitting RT conversion models with refined iRT..."

                            # Calculate refined iRT for training PSMs using callable model
                            refined_irts_training = Vector{Float32}(undef, length(sequences))

                            for (i, seq) in enumerate(sequences)
                                library_irt = psms[:irt_predicted][best_hits][i]
                                # Use callable model interface (zero allocations!)
                                refined_irts_training[i] = refinement_model(seq, library_irt)
                            end

                            # Refit RT → iRT model using REFINED iRT
                            refined_psms_df = DataFrame(
                                rt = psms[:rt][best_hits],
                                irt_predicted = refined_irts_training
                            )

                            rt_model_refined, valid_rt_refined, valid_irt_refined, _ = Pioneer.fit_irt_model(
                                refined_psms_df;
                                lambda_penalty = Float32(0.1),
                                ransac_threshold = 1000,
                                min_psms = 10,
                                spline_degree = 3,
                                max_knots = 7,
                                outlier_threshold = Float32(5.0)
                            )

                            # Refit iRT → RT model using REFINED iRT
                            irt_to_rt_df_refined = DataFrame(
                                rt = valid_irt_refined,           # Refined iRT as input
                                irt_predicted = valid_rt_refined  # RT as output
                            )

                            irt_model_refined, _, _, _ = Pioneer.fit_irt_model(
                                irt_to_rt_df_refined;
                                lambda_penalty = Float32(0.1),
                                ransac_threshold = 1000,
                                min_psms = 10,
                                spline_degree = 3,
                                max_knots = 7,
                                outlier_threshold = Float32(5.0)
                            )

                            # Store REFINED RT models (replacing preliminary models)
                            setRtIrtMap!(search_context, rt_model_refined, ms_file_idx)
                            setIrtRtMap!(search_context, irt_model_refined, ms_file_idx)

                            @user_info "File $ms_file_idx: RT conversion models refitted with refined iRT"
                        else
                            # Refinement not successful - store preliminary models
                            setRtIrtMap!(search_context, rt_model_prelim, ms_file_idx)
                            setIrtRtMap!(search_context, irt_model_prelim, ms_file_idx)
                        end

                        # Step 3: Add :irt_refined column to PSMs file
                        # This column is used by:
                        # 1. plot_irt_comparison() to visualize refinement performance
                        # 2. get_best_precursors_across_runs() for precursor_dict building
                        # Note: SecondPass uses getPredIrt() which computes from model (not this column)
                        psms_path = all_psms_paths[ms_file_idx]
                        add_irt_refined_column!(psms_path, refinement_model, search_context)
                    catch e
                        @user_warn "iRT refinement failed for file $ms_file_idx: $e"
                        # Store preliminary models (refinement failed)
                        setRtIrtMap!(search_context, rt_model_prelim, ms_file_idx)
                        setIrtRtMap!(search_context, irt_model_prelim, ms_file_idx)
                    end
                else
                    # Refinement not enabled - store preliminary models
                    setRtIrtMap!(search_context, rt_model_prelim, ms_file_idx)
                    setIrtRtMap!(search_context, irt_model_prelim, ms_file_idx)
                end

                # Optionally generate plots
                # Use preliminary models for plotting (before refinement)
                if params.plot_rt_alignment
                    plot_rt_alignment_firstpass(
                        valid_rt_prelim,
                        valid_irt_prelim,
                        rt_model_prelim,
                        ms_file_idx,
                        getDataOutDir(search_context)
                    )
                end
            else
                # Use simple UniformSpline (legacy behavior)
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
            end
        catch e
            # Get file name for debugging
            file_name = try
                getFileIdToName(getMSData(search_context), ms_file_idx)
            catch
                "file_$ms_file_idx"
            end
            
            # Safely compute PSM count to avoid excessive output
            n_good_psms = try
                sum(best_hits)
            catch
                0  # Default to 0 if calculation fails
            end
            
            @user_warn "RT mapping failed for MS data file: $file_name ($n_good_psms good PSMs found, need >100 for spline). Using identity RT model."

            # Use identity mapping as fallback - no RT to iRT conversion
            identity_model = IdentityModel()
            setRtIrtMap!(search_context, identity_model, ms_file_idx)
            setIrtRtMap!(search_context, identity_model, ms_file_idx)
        end
    end
    
    # Set identity models for failed files
    ms_data = getMassSpecData(search_context)
    for failed_idx in 1:length(ms_data.file_paths)
        if getFailedIndicator(ms_data, failed_idx)
            file_name = try
                getFileIdToName(getMSData(search_context), failed_idx)
            catch
                "file_$failed_idx"
            end
            @user_warn "Setting identity RT models for failed file: $file_name"
            setRtIrtMap!(search_context, IdentityModel(), failed_idx)
            setIrtRtMap!(search_context, IdentityModel(), failed_idx)
        end
    end
    
    return nothing
end

"""
    plot_rt_alignment_firstpass(rt, irt, rt_model, ms_file_idx, output_dir)

Generate RT alignment diagnostic plot for FirstPassSearch.

# Arguments
- `rt`: Observed retention times
- `irt`: Predicted indexed retention times
- `rt_model`: Fitted RT conversion model
- `ms_file_idx`: MS file index
- `output_dir`: Output directory path

Creates a scatter plot of RT vs iRT with the fitted model curve overlay.
Saves to FirstPass RT alignment folder in QC plots directory.
"""
function plot_rt_alignment_firstpass(
    rt::Vector{Float32},
    irt::Vector{Float32},
    rt_model::RtConversionModel,
    ms_file_idx::Int,
    output_dir::String
)
    n = length(rt)

    # Create scatter plot
    p = plot(
        rt,
        irt,
        seriestype = :scatter,
        title = "FirstPass RT Alignment (File $ms_file_idx)\nn = $n",
        xlabel = "Observed RT (min)",
        ylabel = "Predicted iRT",
        label = nothing,
        alpha = 0.3,
        markersize = 2,
        size = (800, 600)
    )

    # Add fitted model line
    rt_range = LinRange(minimum(rt), maximum(rt), 100)
    fitted_irt = [rt_model(r) for r in rt_range]
    plot!(p, rt_range, fitted_irt, color = :red, linewidth = 2, label = nothing)

    # Create output directory if needed
    rt_plot_folder = joinpath(output_dir, "qc_plots", "rt_alignment_plots", "firstpass")
    !isdir(rt_plot_folder) && mkpath(rt_plot_folder)

    # Save plot
    plot_path = joinpath(rt_plot_folder, "file_$(ms_file_idx)_rt_alignment.pdf")
    savefig(p, plot_path)

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
    irt_errs = get_irt_errs(results.fwhms, precursor_dict, params)


    setIrtErrors!(search_context, irt_errs)

    # Create precursor to iRT mapping
    prec_to_irt = map(x -> (irt=x[:best_irt], mz=x[:mz]), 
                      precursor_dict)

    # Set up indices folder
    rt_indices_folder = joinpath(getDataOutDir(search_context), "temp_data", "rt_indices")
    !isdir(rt_indices_folder) && mkdir(rt_indices_folder)

    # Get PSM paths and RT models
    all_psm_paths = getFirstPassPsms(getMSData(search_context))
    all_rt_models = getRtIrtMap(search_context)

    # Filter to get valid file indices and their paths
    valid_indexed_data = []
    for (idx, path) in enumerate(all_psm_paths)
        if !isempty(path) && haskey(all_rt_models, idx)
            push!(valid_indexed_data, (idx, path))
        end
    end

    # If no valid paths, skip RT index creation
    if isempty(valid_indexed_data)
        @user_warn "No valid PSM files for RT index creation - all files may have failed"
        return nothing
    end

    # Extract just the paths and indices
    valid_indices = [idx for (idx, _) in valid_indexed_data]
    valid_psm_paths = [path for (_, path) in valid_indexed_data]
    valid_rt_models = Dict(i => all_rt_models[idx] for (i, idx) in enumerate(valid_indices))

    # Make RT indices only for valid files
    # Note: prec_to_irt contains refined iRT values from precursor_dict (not library iRT!)
    rt_index_paths = makeRTIndices(
        rt_indices_folder,
        valid_psm_paths,
        prec_to_irt,
        valid_rt_models,
        min_prob=params.max_prob_to_impute
    )

    # Map the results back to original indices
    for (new_idx, orig_idx) in enumerate(valid_indices)
        if new_idx <= length(rt_index_paths)
            setRtIndex!(getMSData(search_context), orig_idx, rt_index_paths[new_idx])
        end
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

"""
    mass_err_ms1(ppm_errs::Vector{Float32}, params::P) where {P<:FragmentIndexSearchParameters}

Estimate MS1 mass error model using exponential distribution fitting (consistent with MS2 approach).

Fits separate exponential distributions to positive and negative mass error tails,
then calculates tolerance bounds at the specified quantile level.

# Arguments
- `ppm_errs`: Vector of PPM mass errors from MS1 precursor matches
- `params`: Search parameters containing fragment error quantile

# Returns
- `(MassErrorModel, Vector{Float32})`: Fitted model and bias-corrected errors

# Details
Uses the same exponential distribution fitting approach as MS2 fragment tolerance estimation:
1. Calculates median bias and removes it from errors
2. Separates positive and negative errors for asymmetric modeling
3. Fits exponential distributions to each tail: Exponential(mean_of_tail)
4. Calculates tolerance bounds at 1-frag_err_quantile level
5. Returns zero tolerances if exponential fitting fails (consistent with MS2)
"""
function mass_err_ms1(ppm_errs::Vector{Float32},
    params::P) where {P<:FragmentIndexSearchParameters}
    mass_err = median(ppm_errs)
    pmp_errs_corrected = ppm_errs .- mass_err

    # Calculate error bounds using exponential distribution fitting (matching MS2 approach)
    frag_err_quantile = getFragErrQuantile(params)

    # Separate positive and negative errors for asymmetric modeling
    r_errs = abs.(pmp_errs_corrected[pmp_errs_corrected .> 0.0f0])
    l_errs = abs.(pmp_errs_corrected[pmp_errs_corrected .< 0.0f0])

    r_bound, l_bound = nothing, nothing
    try
        # Fit exponential distribution to right tail (positive errors)
        # Exponential(θ) where θ is the scale parameter (mean)
        r_bound = quantile(Exponential(
            sum(r_errs)/length(r_errs)  # Mean of positive errors
        ), 1 - frag_err_quantile)

        # Fit exponential distribution to left tail (negative errors)
        l_bound = quantile(Exponential(
            sum(l_errs)/length(l_errs)  # Mean of negative errors
        ), 1 - frag_err_quantile)
    catch e
        @debug_l1 "MS1 exponential fitting failed: length(r_errs)=$(length(r_errs)) length(l_errs)=$(length(l_errs)) sum(r_errs)=$(sum(r_errs)) sum(l_errs)=$(sum(l_errs))"
        # Return zero tolerances on fitting failure (consistent with MS2 approach)
        return MassErrorModel(
            Float32(mass_err),
            (zero(Float32), zero(Float32))
        ), pmp_errs_corrected
    end

    return MassErrorModel(
        Float32(mass_err),
        (Float32(abs(l_bound)), Float32(abs(r_bound)))
    ), pmp_errs_corrected
end

"""
    calculate_irt_error_stats(irt_predicted::Vector{Float64},
                               irt_refined::Vector{Float64},
                               irt_observed::Vector{Float64})

Calculate error statistics for library and refined iRT predictions.

# Returns
NamedTuple with:
- library_errors, refined_errors: Error vectors
- library_mae, library_std, library_mean, library_median, library_correlation
- refined_mae, refined_std, refined_mean, refined_median, refined_correlation
- mae_improvement, mae_improvement_pct
"""
function calculate_irt_error_stats(irt_predicted::Vector{Float64},
                                   irt_refined::Vector{Float64},
                                   irt_observed::Vector{Float64})
    # Calculate error distributions (prediction - observed)
    library_errors = irt_predicted .- irt_observed
    refined_errors = irt_refined .- irt_observed

    # Error statistics
    library_mae = mean(abs.(library_errors))
    library_std = std(library_errors)
    library_median = median(library_errors)
    library_mean = mean(library_errors)

    refined_mae = mean(abs.(refined_errors))
    refined_std = std(refined_errors)
    refined_median = median(refined_errors)
    refined_mean = mean(refined_errors)

    # Improvement metrics
    mae_improvement = library_mae - refined_mae
    mae_improvement_pct = (mae_improvement / library_mae) * 100.0

    # Correlation between predictions and observed
    library_correlation = cor(irt_predicted, irt_observed)
    refined_correlation = cor(irt_refined, irt_observed)

    return (
        library_errors = library_errors,
        refined_errors = refined_errors,
        library_mae = library_mae,
        library_std = library_std,
        library_mean = library_mean,
        library_median = library_median,
        library_correlation = library_correlation,
        refined_mae = refined_mae,
        refined_std = refined_std,
        refined_mean = refined_mean,
        refined_median = refined_median,
        refined_correlation = refined_correlation,
        mae_improvement = mae_improvement,
        mae_improvement_pct = mae_improvement_pct
    )
end

"""
    plot_irt_comparison(psms_paths::Vector{String},
                        rt_irt_models::Dict{Int64, RtConversionModel},
                        output_folder::String,
                        min_prob_for_irt_mapping::Float32)

Create comparison plots and statistics for refined iRT (from column) vs observed iRT (from RT conversion).

Generates two sets of overlaid histograms for each file:
1. High-quality PSMs (prob > min_prob_for_irt_mapping) - matches training data
2. All PSMs - shows performance on full dataset

# Arguments
- `psms_paths`: Paths to PSM Arrow files
- `rt_irt_models`: Dictionary mapping file index → RT→iRT conversion model
- `output_folder`: Base output directory for plots
- `min_prob_for_irt_mapping`: Probability threshold for filtering (matches training)

# Output
- PDFs saved to: `{output_folder}/irt_comparison/file_{idx}_irt_comparison_{filtered|all}.pdf`
- Console output: Statistical summary table for both filtered and all PSMs
"""
function plot_irt_comparison(
    psms_paths::Vector{String},
    rt_irt_models::Dict{Int64, RtConversionModel},
    output_folder::String,
    min_prob_for_irt_mapping::Float32
)
    # Create output directory
    comparison_dir = joinpath(output_folder, "irt_comparison")
    mkpath(comparison_dir)

    @user_info "Creating iRT comparison plots for $(length(psms_paths)) files..."

    # Summary statistics table
    println()
    println("="^100)
    println("iRT Comparison: Refined (Column) vs Observed (RT Conversion)")
    println("="^100)
    println()

    for (file_num, psms_path) in enumerate(psms_paths)
        # Load PSM data
        if !isfile(psms_path)
            @user_warn "PSM file not found: $psms_path, skipping"
            continue
        end

        psms = Arrow.Table(psms_path)

        # Get file index
        if isempty(psms[:ms_file_idx])
            @user_warn "Empty PSM file: $psms_path, skipping"
            continue
        end
        file_idx = first(psms[:ms_file_idx])

        # Get RT→iRT model
        if !haskey(rt_irt_models, file_idx)
            @user_warn "No RT→iRT model for file $file_idx, skipping"
            continue
        end
        rt_to_irt = rt_irt_models[file_idx]

        @user_info "  Reading iRT values from PSM file $(basename(psms_path)):"
        @user_info "    - :irt_refined (refined iRT from model)"
        @user_info "    - :irt_predicted (library iRT)"
        @user_info "    - RT → observed iRT (via RT conversion model)"

        # Extract ALL data
        irt_refined_all = Float64.(psms[:irt_refined])
        irt_predicted_all = Float64.(psms[:irt_predicted])  # Library iRT
        rts_all = Float64.(psms[:rt])
        probs_all = Float64.(psms[:prob])

        # Calculate observed iRT for all
        irt_observed_all = [rt_to_irt(rt) for rt in rts_all]

        # Filter for high-quality PSMs (matching training data)
        high_quality_mask = probs_all .> min_prob_for_irt_mapping
        n_high_quality = sum(high_quality_mask)

        irt_refined_filtered = irt_refined_all[high_quality_mask]
        irt_predicted_filtered = irt_predicted_all[high_quality_mask]
        irt_observed_filtered = irt_observed_all[high_quality_mask]

        # Calculate statistics for both datasets
        stats_all = calculate_irt_error_stats(irt_predicted_all, irt_refined_all, irt_observed_all)
        stats_filtered = calculate_irt_error_stats(irt_predicted_filtered, irt_refined_filtered, irt_observed_filtered)

        # Print statistics
        println("File $file_num (ms_file_idx=$file_idx): $(basename(psms_path))")
        println()
        println("  HIGH-QUALITY PSMs (prob > $(round(min_prob_for_irt_mapping, digits=3))) - N = $n_high_quality")
        println("  " * "-"^90)
        println("    Library iRT Error (before refinement):")
        println("      MAE:    $(round(stats_filtered.library_mae, digits=4))")
        println("      Mean:   $(round(stats_filtered.library_mean, digits=4))")
        println("      Std:    $(round(stats_filtered.library_std, digits=4))")
        println("      Median: $(round(stats_filtered.library_median, digits=4))")
        println("      Correlation: $(round(stats_filtered.library_correlation, digits=4))")
        println()
        println("    Refined iRT Error (after refinement):")
        println("      MAE:    $(round(stats_filtered.refined_mae, digits=4))")
        println("      Mean:   $(round(stats_filtered.refined_mean, digits=4))")
        println("      Std:    $(round(stats_filtered.refined_std, digits=4))")
        println("      Median: $(round(stats_filtered.refined_median, digits=4))")
        println("      Correlation: $(round(stats_filtered.refined_correlation, digits=4))")
        println()
        println("    Improvement:")
        println("      ΔMAE:   $(round(stats_filtered.mae_improvement, digits=4)) ($(round(stats_filtered.mae_improvement_pct, digits=2))%)")
        println("      Δσ:     $(round(stats_filtered.library_std - stats_filtered.refined_std, digits=4))")
        println()
        println("  ALL PSMs - N = $(length(irt_refined_all))")
        println("  " * "-"^90)
        println("    Library iRT Error (before refinement):")
        println("      MAE:    $(round(stats_all.library_mae, digits=4))")
        println("      Mean:   $(round(stats_all.library_mean, digits=4))")
        println("      Std:    $(round(stats_all.library_std, digits=4))")
        println("      Median: $(round(stats_all.library_median, digits=4))")
        println("      Correlation: $(round(stats_all.library_correlation, digits=4))")
        println()
        println("    Refined iRT Error (after refinement):")
        println("      MAE:    $(round(stats_all.refined_mae, digits=4))")
        println("      Mean:   $(round(stats_all.refined_mean, digits=4))")
        println("      Std:    $(round(stats_all.refined_std, digits=4))")
        println("      Median: $(round(stats_all.refined_median, digits=4))")
        println("      Correlation: $(round(stats_all.refined_correlation, digits=4))")
        println()
        println("    Improvement:")
        println("      ΔMAE:   $(round(stats_all.mae_improvement, digits=4)) ($(round(stats_all.mae_improvement_pct, digits=2))%)")
        println("      Δσ:     $(round(stats_all.library_std - stats_all.refined_std, digits=4))")
        println()
        println("-"^100)
        println()

        # Create plots for both filtered and all PSMs
        #using Plots

        # Plot 1: HIGH-QUALITY PSMs (filtered)
        all_errors_filtered = vcat(stats_filtered.library_errors, stats_filtered.refined_errors)
        max_abs_error_filtered = maximum(abs.(all_errors_filtered))
        bins_range_filtered = (-max_abs_error_filtered, max_abs_error_filtered)
        n_bins = 50
        bin_edges_filtered = range(bins_range_filtered[1], bins_range_filtered[2], length=n_bins+1)

        p_filtered = histogram(
            stats_filtered.library_errors,
            bins=bin_edges_filtered,
            alpha=0.6,
            label="Library iRT Error (before refinement)",
            color=:red,
            xlabel="iRT Error (Predicted - Observed)",
            ylabel="Count",
            title="File $file_num: HIGH-QUALITY PSMs (prob > $(round(min_prob_for_irt_mapping, digits=3)))",
            legend=:topright,
            size=(800, 600)
        )

        histogram!(
            p_filtered,
            stats_filtered.refined_errors,
            bins=bin_edges_filtered,
            alpha=0.6,
            label="Refined iRT Error (after refinement)",
            color=:blue
        )

        vline!(p_filtered, [0], color=:black, linestyle=:dash, linewidth=2, label="Perfect prediction")

        stats_text_filtered = "N = $n_high_quality\n" *
                             "Library:  MAE=$(round(stats_filtered.library_mae, digits=3)), σ=$(round(stats_filtered.library_std, digits=3))\n" *
                             "Refined:  MAE=$(round(stats_filtered.refined_mae, digits=3)), σ=$(round(stats_filtered.refined_std, digits=3))\n" *
                             "Improvement: ΔMAE=$(round(stats_filtered.mae_improvement, digits=3)) ($(round(stats_filtered.mae_improvement_pct, digits=1))%)"

        annotate!(
            p_filtered,
            bins_range_filtered[1] + 0.02 * (bins_range_filtered[2] - bins_range_filtered[1]),
            ylims(p_filtered)[2] * 0.95,
            text(stats_text_filtered, :left, 8, :monospace)
        )

        output_path_filtered = joinpath(comparison_dir, "file_$(file_num)_irt_comparison_filtered.pdf")
        savefig(p_filtered, output_path_filtered)
        @user_info "  Saved filtered plot: $output_path_filtered"

        # Plot 2: ALL PSMs
        all_errors_all = vcat(stats_all.library_errors, stats_all.refined_errors)
        max_abs_error_all = maximum(abs.(all_errors_all))
        bins_range_all = (-max_abs_error_all, max_abs_error_all)
        bin_edges_all = range(bins_range_all[1], bins_range_all[2], length=n_bins+1)

        p_all = histogram(
            stats_all.library_errors,
            bins=bin_edges_all,
            alpha=0.6,
            label="Library iRT Error (before refinement)",
            color=:red,
            xlabel="iRT Error (Predicted - Observed)",
            ylabel="Count",
            title="File $file_num: ALL PSMs",
            legend=:topright,
            size=(800, 600)
        )

        histogram!(
            p_all,
            stats_all.refined_errors,
            bins=bin_edges_all,
            alpha=0.6,
            label="Refined iRT Error (after refinement)",
            color=:blue
        )

        vline!(p_all, [0], color=:black, linestyle=:dash, linewidth=2, label="Perfect prediction")

        stats_text_all = "N = $(length(irt_refined_all))\n" *
                        "Library:  MAE=$(round(stats_all.library_mae, digits=3)), σ=$(round(stats_all.library_std, digits=3))\n" *
                        "Refined:  MAE=$(round(stats_all.refined_mae, digits=3)), σ=$(round(stats_all.refined_std, digits=3))\n" *
                        "Improvement: ΔMAE=$(round(stats_all.mae_improvement, digits=3)) ($(round(stats_all.mae_improvement_pct, digits=1))%)"

        annotate!(
            p_all,
            bins_range_all[1] + 0.02 * (bins_range_all[2] - bins_range_all[1]),
            ylims(p_all)[2] * 0.95,
            text(stats_text_all, :left, 8, :monospace)
        )

        output_path_all = joinpath(comparison_dir, "file_$(file_num)_irt_comparison_all.pdf")
        savefig(p_all, output_path_all)
        @user_info "  Saved all PSMs plot: $output_path_all"
    end

    println("="^100)
    println()
    @user_info "✓ iRT comparison plots complete - saved to: $comparison_dir"

    return nothing
end


