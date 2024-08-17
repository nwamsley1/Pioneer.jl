function addPreSearchColumns!(PSMs::DataFrame, 
                                MS_TABLE::Arrow.Table, 
                                prec_is_decoy::Arrow.BoolVector{Bool},
                                prec_irt::Arrow.Primitive{T, Vector{T}},
                                prec_charge::Arrow.Primitive{UInt8, Vector{UInt8}}) where {T<:AbstractFloat}
    
    N = size(PSMs, 1)
    decoys = zeros(Bool, N);
    targets = zeros(Bool, N);
    TIC = zeros(Float16, N);
    charge = zeros(UInt8, N);
    spectrum_peak_count = zeros(UInt32, N);
    iRT_pred = zeros(Float32, N);
    RT = zeros(Float32, N);
    scan_idx::Vector{UInt32} = PSMs[!,:scan_idx]
    precursor_idx::Vector{UInt32} = PSMs[!,:precursor_idx]
    matched_ratio::Vector{Float16} = PSMs[!,:matched_ratio]
    tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    scan_retention_time = MS_TABLE[:retentionTime]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}

    tasks_per_thread = 10
    chunk_size = max(1, size(PSMs, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(PSMs, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
            decoys[i] = prec_is_decoy[precursor_idx[i]];
            targets[i] = decoys[i] == false
            iRT_pred[i] = Float32(prec_irt[precursor_idx[i]]);
            RT[i] = Float32(scan_retention_time[scan_idx[i]]);
            TIC[i] = Float16(log2(tic[scan_idx[i]]));
            charge[i] = UInt8(prec_charge[precursor_idx[i]]);
            matched_ratio[i] = Float16(min(matched_ratio[i], 6e4))
            end
        end
    end
    fetch.(tasks)
    PSMs[!,:matched_ratio] = matched_ratio
    PSMs[!,:decoy] = decoys
    PSMs[!,:iRT_predicted] = iRT_pred
    PSMs[!,:RT] = RT
    PSMs[!,:TIC] = TIC
    PSMs[!,:target] = targets
    PSMs[!,:charge] = charge
    PSMs[!,:spectrum_peak_count] = spectrum_peak_count
    PSMs[!,:intercept] = ones(Float16, size(PSMs, 1))
    PSMs[!,:q_value] = zeros(Float16, N);
    PSMs[!,:prob] = zeros(Float16, N);
    #=
    features = [:entropy_score,:city_block,:scribe,:spectral_contrast,:y_count,:error,
    #:topn,
    :TIC,:intercept]
    PSMs[!,:prob] = zeros(Float32, size(PSMs, 1))
    
    M = size(PSMs, 1)
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:M, chunk_size) # partition your data into chunks that
    β = zeros(Float64, length(features))
    β = ProbitRegression(β, PSMs[!,features], PSMs[!,:target], data_chunks, max_iter = 20)
    ModelPredictProbs!(PSMs[!,:prob], PSMs[!,features], β, data_chunks)
    filter!(:prob => x -> x>=min_prob, PSMs)

    PSMs[!,:best_psms] .= false
    grouped_psms = groupby(PSMs,:precursor_idx)
    for psms in grouped_psms
        best_idx = argmax(psms.prob)
        psms[best_idx,:best_psms] = true
    end
    filter!(x->x.best_psms, PSMs)
    =#
end

function scorePresearch!(PSMs::DataFrame)
    
    features = [:entropy_score,:city_block,
                    :scribe,:spectral_contrast,
                    :y_count,:error,
                    :TIC,:intercept]
    PSMs[!,:prob] = zeros(Float32, size(PSMs, 1))

    M = size(PSMs, 1)
    tasks_per_thread = 10
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:M, chunk_size) # partition your data into chunks that
    β = zeros(Float64, length(features))
    β = ProbitRegression(β, PSMs[!,features], PSMs[!,:target], data_chunks, max_iter = 20)
    ModelPredictProbs!(PSMs[!,:prob], PSMs[!,features], β, data_chunks)
end

function addMainSearchColumns!(PSMs::DataFrame, 
                                MS_TABLE::Arrow.Table, 
                                RT_iRT::UniformSpline,
                                structural_mods::Arrow.List{String, Int32, Vector{UInt8}},
                                prec_missed_cleavages::Arrow.Primitive{UInt8, Vector{UInt8}},
                                prec_is_decoy::Arrow.BoolVector{Bool},
                                prec_irt::Arrow.Primitive{T, Vector{T}},
                                prec_charge::Arrow.Primitive{UInt8, Vector{UInt8}};
                                ) where {T<:AbstractFloat}
    
    ###########################
    #Allocate new columns
    N = size(PSMs, 1)
    missed_cleavage = zeros(UInt8, N);
    Mox = zeros(UInt8, N);
    iRT_pred = zeros(Float32, N);
    RT = zeros(Float32, N);
    iRT = zeros(Float32, N);
    TIC = zeros(Float16, N);
    charge = zeros(UInt8, N);
    err_norm = zeros(Float16, N);
    targets = zeros(Bool, N);
    spectrum_peak_count = zeros(UInt32, N);
    scan_idx::Vector{UInt32} = PSMs[!,:scan_idx]
    precursor_idx::Vector{UInt32} = PSMs[!,:precursor_idx] 
    error::Vector{Float32} = PSMs[!,:error]
    total_ions::Vector{UInt8} = PSMs[!,:y_count] .+ PSMs[!,:b_count]
    tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    scan_retention_time = MS_TABLE[:retentionTime]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    masses = MS_TABLE[:masses]::Arrow.List{Union{Missing, SubArray{Union{Missing, Float32}, 1, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}, Tuple{UnitRange{Int64}}, true}}, Int64, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}}
    function countMOX(seq::String)
        return UInt8(count("Unimod:35", seq))
    end

    #Split data into chunk ranges
    tasks_per_thread = 10
    chunk_size = max(1, size(PSMs, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(PSMs, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                prec_idx = precursor_idx[i]

                targets[i] = prec_is_decoy[prec_idx] == false;
                missed_cleavage[i] = prec_missed_cleavages[prec_idx]
                Mox[i] = countMOX(structural_mods[prec_idx])::UInt8 #UInt8(length(collect(eachmatch(r"ox",  precursors[precursor_idx[i]].sequence))))
                iRT_pred[i] = Float32(prec_irt[prec_idx]);
                RT[i] = Float32(scan_retention_time[scan_idx[i]]);
                iRT[i] = RT_iRT(RT[i])
                TIC[i] = Float16(log2(tic[scan_idx[i]]));
                charge[i] = UInt8(prec_charge[prec_idx]);
                err_norm[i] = min(Float16.(abs(error[i]/total_ions[i])), Float16(6e4))#Float16(min(abs((error[i])/(total_ions[i])), 6e4))
                spectrum_peak_count[i] = UInt32(length(masses[scan_idx[i]]))
            end
        end
    end
    fetch.(tasks)
    PSMs[!,:iRT_predicted] = iRT_pred
    PSMs[!,:RT] = RT
    PSMs[!,:iRT] = iRT
    PSMs[!,:TIC] = TIC
    PSMs[!,:total_ions] = total_ions
    PSMs[!,:err_norm] = err_norm
    PSMs[!,:target] = targets
    PSMs[!,:missed_cleavage] = missed_cleavage
    PSMs[!,:Mox] = Mox
    PSMs[!,:charge] = charge
    PSMs[!,:spectrum_peak_count] = spectrum_peak_count
    PSMs[!,:score] = zeros(Float32, N);
    PSMs[!,:q_value] = zeros(Float16, N);
    PSMs[!,:intercept] = ones(Float16, N)
end

#=
function getRTErrs!(psms::DataFrame; 
                    n_bins::Int = 200, 
                    bandwidth::Float64 = 1.0,
                    min_prob::Float64 = 0.95)

    best_psms_bool = psms[!,:prob].>= min_prob

    linear_spline = KDEmapping(
                                psms[best_psms_bool,:RT],
                                psms[best_psms_bool,:iRT_predicted],
                                n = n_bins,
                                bandwidth = bandwidth,
                            )

    psms[!,:iRT_observed] = Float16.(linear_spline(psms[!,:RT]))
    psms[!,:iRT_error] = Float16.(abs.(psms[!,:iRT_observed] .- psms[!,:iRT_predicted]))
end
=#

function scoreMainSearchPSMs!(psms::DataFrame, column_names::Vector{Symbol};
                             n_train_rounds::Int64 = 2,
                             max_iter_per_round::Int64 = 50,
                             max_q_value::Float64 = 0.01)

    β = zeros(Float64, length(column_names));
    best_psms = nothing#ones(Bool, size(psms, 1))

    tasks_per_thread = 10
    M = size(psms, 1)
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:M, chunk_size) # partition your data into chunks that

    for i in range(1, n_train_rounds)
        if i < 2 #Train on all data during first round
            β = ProbitRegression(β, psms[!,column_names], psms[!,:target], data_chunks, max_iter = max_iter_per_round);
            ModelPredict!(psms[!,:score], psms[!,column_names], β, data_chunks); #Get Z-scores 
        else #Train using only decoys and high scoring targets from the previous round
            psms_targets = psms[best_psms,:target]
            M = size(psms_targets, 1)
            sub_chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
            sub_data_chunks = partition(1:M, sub_chunk_size) # partition your data into chunks that
            β = ProbitRegression(β, psms[best_psms,column_names], psms_targets, sub_data_chunks, max_iter = max_iter_per_round);
            ModelPredict!(psms[!,:score], psms[!,column_names], β, data_chunks); #Get Z-scores 
        end
      
        getQvalues!(psms[!,:score],psms[!,:target],psms[!,:q_value]);
        if i < n_train_rounds #Get Data to train on during subsequent round
            best_psms = ((psms[!,:q_value].<=max_q_value).&(psms[!,:target])) .| (psms[!,:target].==false);
        end
    end
    return 
end

function getProbs!(psms::DataFrame)
    psms[!,:prob] = zeros(Float32, size(psms, 1))
    tasks_per_thread = 10
    M = size(psms, 1)
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:M, chunk_size) # partition your data into chunks that
    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                psms[i,:prob] = (1 + SpecialFunctions.erf(psms[i,:score]/sqrt(2)))/2
            end
        end
    end
    fetch.(tasks)
    return nothing
end


function addSecondSearchColumns!(psms::DataFrame, 
                        MS_TABLE::Arrow.Table, 
                        prec_mz::AbstractVector{T},
                        prec_charge::AbstractVector{UInt8},
                        prec_is_decoy::AbstractVector{Bool},
                        prec_id_to_cv_fold::Dictionary{UInt32, UInt8}) where {T<:AbstractFloat}
    
    ###########################
    #Correct Weights by base-peak intensity
    filter!(x->x.weight>0.0, psms);
    ###########################
    #Allocate new columns
   
    #Threads.@threads for i in ProgressBar(range(1, size(PSMs)[1]))
    N = size(psms, 1);
    decoys = zeros(Bool, N);
    RT = zeros(Float32, N);
    #TIC = zeros(Float16, N);
    total_ions = zeros(UInt16, N);
    err_norm = zeros(Float16, N);
    targets = zeros(Bool, N);
    prec_charges = zeros(UInt8, N);
    #prec_mzs = zeros(Float32, N);
    cv_fold = zeros(UInt8, N);
    scan_idxs::Vector{UInt32} = psms[!,:scan_idx]
    prec_idxs::Vector{UInt32} = psms[!,:precursor_idx]
    y_count::Vector{UInt8} = psms[!,:y_count]
    b_count::Vector{UInt8} = psms[!,:b_count]
    isotope_count::Vector{UInt8} = psms[!,:isotope_count]
    error::Vector{Float32} = psms[!,:error]
    #PSMs[!,:total_ions]
    #tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    scan_retention_time = MS_TABLE[:retentionTime]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    matched_ratio::Vector{Float16} = psms[!,:matched_ratio]

    tasks_per_thread = 5
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                prec_idx = prec_idxs[i]
                scan_idx = scan_idxs[i]

                decoys[i] = prec_is_decoy[prec_idx];
                targets[i] = decoys[i] == false
                RT[i] = Float32(scan_retention_time[scan_idx]);
                prec_charges[i] = prec_charge[prec_idx];
                #prec_mzs[i] = prec_mz[prec_idx];
                #TIC[i] = Float16(log2(tic[scan_idx]));
                total_ions[i] = UInt16(y_count[i] + b_count[i] + isotope_count[i]);
                err_norm[i] = min(Float16((error[i])/(total_ions[i])), 6e4)
                if isinf(matched_ratio[i])
                    matched_ratio[i] = Float16(60000)*sign(matched_ratio[i])
                end
                cv_fold[i] = prec_id_to_cv_fold[prec_idx]
            end
        end
    end
    fetch.(tasks)
    psms[!,:matched_ratio] = matched_ratio
    psms[!,:decoy] = decoys
    psms[!,:RT] = RT
    #psms[!,:TIC] = TIC
    psms[!,:total_ions] = total_ions
    psms[!,:err_norm] = err_norm
    psms[!,:target] = targets
    psms[!,:charge] = prec_charges
    #psms[!,:prec_mz] = prec_mzs
    psms[!,:cv_fold] = cv_fold

    #psms[!,:best_scan] = zeros(Bool, N)
    #psms[!,:score] = zeros(Float32, size(psms, 1))
    psms[!,:intercept] = ones(Float16, size(psms, 1))
    #######
    sort!(psms,:RT); #Sorting before grouping is critical. 
    return nothing
end

function addIntegrationFeatures!(psms::DataFrame)
    ###########################
    #Add columns
    ###########################
    #Add columns
    new_cols = [
    (:peak_area,                   Union{Float32, Missing})
    (:GOF,                      Union{Float16, Missing})
    (:FWHM,                     Union{Float16, Missing})
    (:FWHM_01,                  Union{Float16, Missing})
    (:assymetry,                 Union{Float16, Missing})
    (:points_above_FWHM,        Union{UInt16, Missing})
    (:points_above_FWHM_01,     Union{UInt16, Missing})
    (:σ,                        Union{Float32, Missing})
    (:tᵣ,                       Union{Float16, Missing})
    (:τ,                        Union{Float32, Missing})
    (:H,                        Union{Float32, Missing})
    (:max_score,           Union{Float16, Missing})
    (:mean_score,           Union{Float16, Missing})
    (:max_entropy,              Union{Float16, Missing})
    (:max_scribe_score,     Union{Float16, Missing})
    (:max_city_fitted,     Union{Float16, Missing})
    (:mean_city_fitted,           Union{Float16, Missing})
    (:y_ions_sum,                 Union{UInt16, Missing})
    (:max_y_ions,                 Union{UInt16, Missing})
    (:data_points,              Union{UInt32, Missing})
    (:fraction_censored,              Union{Float16, Missing})
    (:max_matched_ratio,       Union{Float16, Missing})
    (:base_width_min,           Union{Float16, Missing})
    (:best_scan, Union{Bool, Missing})
    ];

    for column in new_cols
        col_type = last(column);
        col_name = first(column)
        psms[!,col_name] = zeros(col_type, size(psms, 1))
    end

end

function getIsotopesCaptured!(chroms::DataFrame, 
                        prec_charge::AbstractArray{UInt8},
                        prec_mz::AbstractArray{Float32},
                        MS_TABLE::Arrow.Table)
    #sum(MS2_CHROMS.weight.!=0.0)
    chroms[!,:isotopes_captured] = Vector{Tuple{Int8, Int8}}(undef, size(chroms, 1))
    
    tasks_per_thread = 5
    chunk_size = max(1, size(chroms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(chroms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        #last_mz = 0.0f0
        #last_prec_id = zero(UInt32)
        #isotopes = (zero(UInt8), zero(UInt8))
        Threads.@spawn begin
            for i in chunk
                
                prec_id = chroms[i,:precursor_idx]
                mz = prec_mz[prec_id]
                #if (abs(mz - last_mz)<1e-6) & (last_prec_id == prec_id)
                #    chroms[i,:isotopes_captured] = isotopes
                #    continue
                #end
                charge = prec_charge[prec_id]

                scan_id = chroms[i,:scan_idx]
                scan_mz = MS_TABLE[:centerMass][scan_id]
                window_width = MS_TABLE[:isolationWidth][scan_id]

                isotopes = getPrecursorIsotopeSet(mz, 
                                                    charge, 
                                                    Float32(scan_mz-window_width/2),
                                                    Float32(scan_mz+window_width/2) 
                                                    )                
                chroms[i,:isotopes_captured] = isotopes
                #last_mz = mz
                #last_prec_id = prec_id
            end
        end
    end
    fetch.(tasks)
    return nothing
end


function initSummaryColumns!(
    psms::DataFrame,
    )

    new_cols = [
        (:max_entropy,              Float16)
        (:max_gof,         Float16)
        (:max_fitted_manhattan_distance,          Float16)
        (:max_fitted_spectral_contrast,         Float16)
        (:y_ions_sum,               UInt16)
        (:max_y_ions,               UInt16)
        (:max_matched_ratio,        Float16)
        ];

        N = size(psms, 1)
        for column in new_cols
            col_type = last(column);
            col_name = first(column)
            psms[!,col_name] = zeros(col_type, N)
        end
        return psms
end

function scoreSecondSearchPSMs!(psms::DataFrame, 
                                features::Vector{Symbol},
                                max_iter::Int64 = 20,
                                k_cv_folds::Int64 = 2,
                                tasks_per_thread::Int = 3)
    cv_folds = psms[!,:cv_fold]::Vector{UInt8}
    M = size(psms, 1)
    allcv_chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    allcv_data_chunks = partition(1:M, allcv_chunk_size) # partition your data into chunks that
    β = zeros(Float64, length(features));
    fill!(psms[!,:score], zero(eltype(psms[!,:score])))
    for k in range(0, k_cv_folds - 1)
        cv_fold = cv_folds.==UInt8(k)
        M = sum(cv_fold)#size(PSMS, 1)
        #println("k $k, M $M")
        chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
        data_chunks = partition(1:M, chunk_size) # partition your data into chunks that
        β = zeros(Float64, length(features));
        β = ProbitRegression(β, psms[cv_fold,features], psms[cv_fold,:target], data_chunks, max_iter = max_iter);
        ModelPredictCVFold!(psms[!,:prob], psms[!,:cv_fold], UInt8(k), psms[!,features], β, allcv_data_chunks)  
        #println("sum(PSMs[cv_fold,:score]) ", sum(psms[cv_fold,:score]))
    end
end

function addPostIntegrationFeatures!(psms::DataFrame, 
                                    MS_TABLE::Arrow.Table, 
                                    precursor_sequence::AbstractArray{String},
                                    structural_mods::AbstractArray{String},
                                    prec_mz::AbstractArray{T},
                                    prec_irt::AbstractArray{T},
                                    prec_charge::AbstractArray{UInt8},
                                    precursor_missed_cleavage::AbstractArray{UInt8},
                                    ms_file_idx::Integer,
                                    ms_id_to_file_path::Dict{Int64, String},
                                    rt_to_irt_interp::Any,
                                    prec_id_to_irt::Dictionary{UInt32, @NamedTuple{irt::Float32, mz::Float32}}) where {T<:AbstractFloat}

    #filter!(x -> x.best_scan, psms);
    filter!(x->x.weight>0, psms);
    #filter!(x->x.data_points>0, psms)
    ###########################
    #Allocate new columns
    #println("TEST")
    N = size(psms, 1)
    irt_diff = zeros(Float32, N)
    irt_obs = zeros(Float32, N)
    irt_pred = zeros(Float32, N)
    irt_error = zeros(Float32, N)

    #PSMs[!,:missed_cleavage] .= zero(UInt8);
    #PSMs[!,:sequence] .= "";
    missed_cleavage = zeros(UInt8, N);
    #sequence = Vector{String}(undef, N);
    #stripped_sequence = Vector{String}(undef, N);
    adjusted_intensity_explained = zeros(Float16, N);
    #log2_base_peak_intensity = zeros(Float16, N);
    prec_charges = zeros(UInt8, N)
    sequence_length = zeros(UInt8, N);
    #b_y_overlap = zeros(Bool, N);
    spectrum_peak_count = zeros(Float16, N);
    prec_mzs = zeros(Float32, size(psms, 1));
    Mox = zeros(UInt8, N);
    TIC = zeros(Float16, N);

    tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    precursor_idx::Vector{UInt32} = psms[!,:precursor_idx] 
    scan_idx::Vector{UInt32} = psms[!,:scan_idx]
    masses = MS_TABLE[:masses]::Arrow.List{Union{Missing, SubArray{Union{Missing, Float32}, 1, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}, Tuple{UnitRange{Int64}}, true}}, Int64, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}}
    #longest_y::Vector{UInt8} = psms[!,:longest_y]
    #longest_b::Vector{UInt8} = psms[!,:longest_b]
    rt::Vector{Float32} = psms[!,:RT]
    tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    log2_intensity_explained = psms[!,:log2_intensity_explained]::Vector{Float16}
    #precursor_idx = PSMs[!,:precursor_idx]::Vector{UInt32}
    function countMOX(seq::String)
        return UInt8(count("Unimod:35", seq))
    end


    tasks_per_thread = 5
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin 
            for i in chunk
                prec_idx = precursor_idx[i]

                irt_obs[i] = rt_to_irt_interp[ms_id_to_file_path[ms_file_idx]](rt[i])
                irt_pred[i] = prec_irt[prec_idx]
                irt_diff[i] = abs(irt_obs[i] - first(prec_id_to_irt[prec_idx]))
                irt_error[i] = abs(irt_obs[i] - irt_pred[i])

                missed_cleavage[i] = precursor_missed_cleavage[prec_idx]
                #sequence[i] = precursor_sequence[prec_idx]
                sequence_length[i] = length(replace(precursor_sequence[prec_idx], r"\(.*?\)" => ""))#replace.(sequence[i], "M(ox)" => "M");
                Mox[i] = countMOX(structural_mods[prec_idx])::UInt8
                #sequence_length[i] = length(stripped_sequence[i])


                TIC[i] = Float16(log2(tic[scan_idx[i]]))
                adjusted_intensity_explained[i] = Float16(log2(TIC[i]) + log2_intensity_explained[i]);
                prec_charges[i] = prec_charge[prec_idx]
                #b_y_overlap[i] = ((sequence_length[i] - longest_y[i])>longest_b[i]) &  (longest_b[i] > 0) & (longest_y[i] > 0);

                spectrum_peak_count[i] = length(masses[scan_idx[i]])
         
                prec_mzs[i] = prec_mz[prec_idx];
            end
        end
    end
    fetch.(tasks)
    psms[!,:irt_obs] = irt_obs
    psms[!,:irt_pred] = irt_pred
    psms[!,:irt_diff] = irt_diff
    psms[!,:irt_error] = irt_error

    psms[!,:missed_cleavage] = missed_cleavage
    #psms[!,:sequence] = sequence
    #psms[!,:stripped_sequence] = stripped_sequence
    psms[!,:Mox] = Mox
    psms[!,:sequence_length] = sequence_length

    psms[!,:tic] = TIC
    psms[!,:adjusted_intensity_explained] = adjusted_intensity_explained
    psms[!,:charge] = prec_charges
    
    #psms[!,:b_y_overlap] = b_y_overlap
    psms[!,:spectrum_peak_count] = spectrum_peak_count
    psms[!,:prec_mz] = prec_mzs
    psms[!,:ms_file_idx] .= ms_file_idx
    return nothing
end