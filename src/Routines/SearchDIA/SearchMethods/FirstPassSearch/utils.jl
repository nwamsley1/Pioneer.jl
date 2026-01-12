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
    map_retention_times!(search_context::SearchContext, params::FirstPassSearchParameters, all_psms_paths::Vector{String})

Maps retention times between library and empirical scales for each MS file. For files with sufficient high-confidence PSMs
(probability > min_prob_for_irt_mapping), creates spline-based conversion models between retention time (RT) and indexed
retention time (iRT) scales. Includes edge trimming to remove sparse data at RT boundaries.

# Arguments
- `search_context`: Contains MS data and mapping information
- `params`: Search parameters including minimum probability threshold
- `all_psms_paths`: Paths to PSM Arrow files for each MS file

Creates and stores:
- RT to iRT mapping splines
- iRT to RT mapping splines

Uses identity model as fallback if insufficient PSMs or fitting fails.
"""
function map_retention_times!(
    search_context::SearchContext,
    params::FirstPassSearchParameters,
    all_psms_paths::Vector{String}
)
    for ms_file_idx in eachindex(all_psms_paths)
        file_name = try
            getFileIdToName(getMSData(search_context), ms_file_idx)
        catch
            "file_$ms_file_idx"
        end
        psms_path = all_psms_paths[ms_file_idx]
        psms = Arrow.Table(psms_path)
        best_hits = psms[:prob].>params.min_prob_for_irt_mapping

        try
            # Check if we have enough PSMs for RT alignment
            n_good_psms = sum(best_hits)

            if n_good_psms < 100
                # Get file name for warning
                file_name = try
                    getFileIdToName(getMSData(search_context), ms_file_idx)
                catch
                    "file_$ms_file_idx"
                end

                @user_warn "Too few PSMs ($n_good_psms) for RT alignment in file: $file_name (need ≥100). Using identity RT model."

                # Use identity mapping as fallback
                identity_model = IdentityModel()
                setRtIrtMap!(search_context, identity_model, ms_file_idx)
                setIrtRtMap!(search_context, identity_model, ms_file_idx)
                continue
            end

            # Extract best PSM data
            best_rts = psms[:rt][best_hits]
            best_irts = psms[:irt_predicted][best_hits]

            # Calculate adaptive knots: 1 per 100 PSMs, min 5, max 100
            n_knots = clamp(n_good_psms ÷ 100, 5, 100)

            # Trim edges where data is sparse to avoid unstable fitting
            min_psms_per_bin = 100  # Minimum PSMs required per bin

            if n_good_psms >= 1000 && n_knots >= 10  # Only trim if enough data and knots
                # Create uniform bin edges based on initial RT range
                rt_min_initial, rt_max_initial = extrema(best_rts)
                bin_edges = collect(LinRange(rt_min_initial, rt_max_initial, n_knots + 1))

                # Count PSMs per bin using optimized histogram
                hist = StatsBase.fit(Histogram, best_rts, bin_edges)
                bin_counts = hist.weights

                # Find first bin from left with sufficient data
                left_bin_idx = findfirst(count -> count >= min_psms_per_bin, bin_counts)
                # Find first bin from right with sufficient data
                right_bin_idx = findlast(count -> count >= min_psms_per_bin, bin_counts)

                if !isnothing(left_bin_idx) && !isnothing(right_bin_idx) && left_bin_idx <= right_bin_idx
                    # Trim to dense region
                    rt_min_trimmed = bin_edges[left_bin_idx]
                    rt_max_trimmed = bin_edges[right_bin_idx + 1]

                    # Cap trimming at 2.5% per side (5% total max)
                    max_trim_per_side = 0.025
                    left_trim_pct = sum(best_rts .< rt_min_trimmed) / n_good_psms
                    right_trim_pct = sum(best_rts .> rt_max_trimmed) / n_good_psms

                    if left_trim_pct > max_trim_per_side || right_trim_pct > max_trim_per_side
                        # Bin-based trimming too aggressive, fall back to percentile-based
                        rt_min_trimmed = Float32(quantile(best_rts, max_trim_per_side))
                        rt_max_trimmed = Float32(quantile(best_rts, 1.0 - max_trim_per_side))
                        @user_info "Edge trimming capped at 2.5% per side (bin-based would trim $(round(100*left_trim_pct, digits=1))% left, $(round(100*right_trim_pct, digits=1))% right)"
                    end

                    # Apply trimming
                    keep_mask = (best_rts .>= rt_min_trimmed) .& (best_rts .<= rt_max_trimmed)
                    n_excluded = n_good_psms - sum(keep_mask)

                    if n_excluded > 0
                        best_rts = best_rts[keep_mask]
                        best_irts = best_irts[keep_mask]

                        pct_excluded = round(100 * n_excluded / n_good_psms, digits=1)
                        #@user_info "Trimmed $n_excluded PSMs ($pct_excluded%) from sparse edge bins (RT: $(round(rt_min_trimmed, digits=2))-$(round(rt_max_trimmed, digits=2)) min)"
                    end
                end
            end

            # Fit adaptive UniformSpline for RT → iRT conversion
            rt_to_irt_spline = UniformSpline(
                best_irts,    # y values (iRT)
                best_rts,     # x values (RT)
                3,            # degree (cubic)
                n_knots       # adaptive based on PSM count
            )

            # Fit adaptive UniformSpline for iRT → RT conversion (inverse)
            irt_to_rt_spline = UniformSpline(
                best_rts,     # y values (RT)
                best_irts,    # x values (iRT)
                3,            # degree (cubic)
                n_knots       # adaptive based on PSM count
            )

            # Store models in SearchContext
            setRtIrtMap!(search_context, SplineRtConversionModel(rt_to_irt_spline), ms_file_idx)
            setIrtRtMap!(search_context, SplineRtConversionModel(irt_to_rt_spline), ms_file_idx)

            # Optionally generate diagnostic plots
            if params.plot_rt_alignment
                plot_rt_alignment_firstpass(
                    best_rts,
                    best_irts,
                    SplineRtConversionModel(rt_to_irt_spline),
                    ms_file_idx,
                    getDataOutDir(search_context)
                )
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

            @user_warn "RT mapping failed for MS data file: $file_name ($n_good_psms PSMs). Error: $e. Using identity RT model."

            # Use identity mapping as fallback - no RT to iRT conversion
            identity_model = IdentityModel()
            setRtIrtMap!(search_context, identity_model, ms_file_idx)
            setIrtRtMap!(search_context, identity_model, ms_file_idx)
        end
    end
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

    # Log RT tolerances for each file
   # @user_info "FirstPassSearch RT tolerances (per file):"
    ms_data = getMSData(search_context)
    for (file_idx, tol) in pairs(irt_errs)
        file_name = getFileIdToName(ms_data, file_idx)
        @debug_l1 "  File $(file_idx) ($file_name): RT tol = $(round(tol, digits=3)) min"
    end

    # Log summary statistics
    if !isempty(irt_errs)
        tol_values = collect(values(irt_errs))
        @debug_l1 "  Summary: min=$(round(minimum(tol_values), digits=3)), max=$(round(maximum(tol_values), digits=3)), median=$(round(median(tol_values), digits=3)) min"
    end

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



