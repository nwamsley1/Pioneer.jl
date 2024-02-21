function addPreSearchColumns!(PSMs::DataFrame, MS_TABLE::Arrow.Table, precursors::Vector{LibraryPrecursorIon{T}}; 
                        max_rt_error::Float64 = 20.0,  
                        min_prob::Float64 = 0.9, 
                        n_bins::Int = 200, granularity::Int = 50) where {T<:AbstractFloat}
    
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
            decoys[i] = isDecoy(precursors[precursor_idx[i]]);
            targets[i] = decoys[i] == false
            iRT_pred[i] = Float32(getIRT(precursors[precursor_idx[i]]));
            RT[i] = Float32(scan_retention_time[scan_idx[i]]);
            TIC[i] = Float16(log2(tic[scan_idx[i]]));
            charge[i] = UInt8(getPrecCharge(precursors[precursor_idx[i]]));
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

    features = [:entropy_score,:city_block,:scribe,:spectral_contrast,:y_count,:error,:topn,:TIC,:intercept]
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

end

function addMainSearchColumns!(PSMs::DataFrame, 
                                MS_TABLE::Arrow.Table, 
                                precursors::Vector{LibraryPrecursorIon{T}}
                                ) where {T<:AbstractFloat}
    
    ###########################
    #Allocate new columns
    N = size(PSMs, 1)
    missed_cleavage = zeros(UInt8, N);
    Mox = zeros(UInt8, N);
    iRT_pred = zeros(Float32, N);
    RT = zeros(Float32, N);
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
        mox = zero(UInt8)
        in_mox = false
        for aa in seq
            if in_mox
                if aa == 'x'
                    mox += one(UInt8)
                end
                in_mox = false
            end
            if aa == 'o'
                in_mox = true
            end
        end
        return mox
    end

    #Split data into chunk ranges
    tasks_per_thread = 10
    chunk_size = max(1, size(PSMs, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(PSMs, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                prec = precursors[precursor_idx[i]]
                targets[i] = isDecoy(prec) == false;
                missed_cleavage[i] = prec.missed_cleavages
                Mox[i] = countMOX(prec.sequence)::UInt8 #UInt8(length(collect(eachmatch(r"ox",  precursors[precursor_idx[i]].sequence))))
                iRT_pred[i] = Float32(getIRT(prec));
                RT[i] = Float32(scan_retention_time[scan_idx[i]]);
                TIC[i] = Float16(log2(tic[scan_idx[i]]));
                charge[i] = UInt8(getPrecCharge(prec));
                err_norm[i] = min(Float16.(abs(error[i]/total_ions[i])), Float16(6e4))#Float16(min(abs((error[i])/(total_ions[i])), 6e4))
                spectrum_peak_count[i] = UInt32(length(masses[scan_idx[i]]))
            end
        end
    end
    fetch.(tasks)
    PSMs[!,:iRT_predicted] = iRT_pred
    PSMs[!,:RT] = RT
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

function getRTErrs!(psms::DataFrame; 
                    n_bins::Int = 200, bandwidth::Float64 = 1.0, w::Int = 11,
                    min_spectral_contrast::Float64 = 0.9,
                    min_entropy_score::Float64 = 0.9,
                    min_total_ions::Int64 = 6)

    best_psms_bool = (psms[!,:spectral_contrast].>min_spectral_contrast
    ) .& (psms[!,:target]
    ) .& (psms[!,:entropy_score].>min_entropy_score
    ) .& (psms[!,:total_ions].>min_total_ions)

    linear_spline = KDEmapping(
                                psms[best_psms_bool,:RT],
                                psms[best_psms_bool,:iRT_predicted],
                                n = n_bins,
                                bandwidth = bandwidth,
                                w = w
                            )

    psms[:,:iRT_observed] = Float16.(linear_spline(psms[:,:RT]))
    psms[!,:iRT_error] = Float16.(abs.(psms[!,:iRT_observed] .- psms[!,:iRT_predicted]))
end

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

function getBestPSMs!(psms::DataFrame,
                        precursors::Vector{LibraryPrecursorIon{T}}; 
                        max_q_value::Float64 = 0.10) where {T<:AbstractFloat}

    #Remove psms below q_value threshold
    #psms = psms[!,[:precursor_idx,:RT,:iRT_predicted,:q_value,:score]]
    #could this be optimized to make type stable?
    filter!(x->x.q_value<=max_q_value, psms);

    psms[!,:best_psm] = zeros(Bool, size(psms, 1))
    gpsms = groupby(psms,:precursor_idx)
    for (precursor_idx, prec_psms) in pairs(gpsms)
        best_psm_idx = argmax(prec_psms[!,:score])
        prec_psms[best_psm_idx,:best_psm] = true
    end

    filter!(x->x.best_psm, psms);
    select!(psms, [:precursor_idx,:RT,:iRT_predicted,:q_value,:score,:prob])

    prec_mz = zeros(Float32, size(psms, 1));
    precursor_idx = psms[!,:precursor_idx]::Vector{UInt32}
    Threads.@threads for i in range(1, size(psms, 1))
        prec_mz[i] = precursors[precursor_idx[i]].mz::Float32;
    end
    psms[!,:prec_mz] = prec_mz

    return
end

function addSecondSearchColumns!(psms::DataFrame, 
                        MS_TABLE::Arrow.Table, 
                        precursors::Vector{LibraryPrecursorIon{T}},
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
    TIC = zeros(Float16, N);
    total_ions = zeros(UInt16, N);
    err_norm = zeros(Float16, N);
    targets = zeros(Bool, N);
    charge = zeros(UInt8, N);
    prec_mz = zeros(Float32, N);
    cv_fold = zeros(UInt8, N);
    scan_idx::Vector{UInt32} = psms[!,:scan_idx]
    precursor_idx::Vector{UInt32} = psms[!,:precursor_idx]
    y_count::Vector{UInt8} = psms[!,:y_count]
    b_count::Vector{UInt8} = psms[!,:b_count]
    error::Vector{Float32} = psms[!,:error]
    #PSMs[!,:total_ions]
    tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    scan_retention_time = MS_TABLE[:retentionTime]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    matched_ratio::Vector{Float16} = psms[!,:matched_ratio]

    tasks_per_thread = 5
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                decoys[i] = isDecoy(precursors[precursor_idx[i]]);
                targets[i] = decoys[i] == false
                RT[i] = Float32(scan_retention_time[scan_idx[i]]);
                charge[i] = UInt8(precursors[precursor_idx[i]].prec_charge);
                prec_mz[i] = Float32(getMZ(precursors[precursor_idx[i]]));
                TIC[i] = Float16(log2(tic[scan_idx[i]]));
                total_ions[i] = UInt16(y_count[i] + b_count[i]);
                err_norm[i] = min(Float16((error[i])/(total_ions[i])), 6e4)
                if isinf(matched_ratio[i])
                    matched_ratio[i] = Float16(60000)*sign(matched_ratio[i])
                end
                cv_fold[i] = prec_id_to_cv_fold[precursor_idx[i]]
            end
        end
    end
    fetch.(tasks)
    psms[!,:matched_ratio] = matched_ratio
    psms[!,:decoy] = decoys
    psms[!,:RT] = RT
    psms[!,:TIC] = TIC
    psms[!,:total_ions] = total_ions
    psms[!,:err_norm] = err_norm
    psms[!,:target] = targets
    psms[!,:charge] = charge
    psms[!,:prec_mz] = prec_mz
    psms[!,:cv_fold] = cv_fold

    psms[!,:best_scan] = zeros(Bool, N)
    psms[!,:score] = zeros(Float32, size(psms, 1))
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

function getIsoRanks!(psms::DataFrame, 
                        MS_TABLE::Arrow.Table,
                        window_width::Float64)
    #sum(MS2_CHROMS.weight.!=0.0)
    psms[!,:iso_rank] = zeros(UInt8, size(psms, 1))

    
    tasks_per_thread = 5
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                charge = psms[i,:charge]
                mz = psms[i,:prec_mz]
                scan_id = psms[i,:scan_idx]
                scan_mz = MS_TABLE[:precursorMZ][scan_id]

                window = (Float32(scan_mz-window_width/2), 
                        Float32(scan_mz+window_width/2)

                        )
                isotopes = getPrecursorIsotopeSet(mz, charge, window)

                rank = zero(UInt8)
                if iszero(first(isotopes))
                    if last(isotopes) > 1
                        rank = UInt8(1)
                    elseif last(isotopes) == 1
                        rank = UInt8(2)
                    else
                        rank = UInt8(3)
                    end
                else
                    rank = UInt8(4)
                end
                psms[i,:iso_rank] = rank
            end
        end
    end
    fetch.(tasks)
    return nothing
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
        println("k $k, M $M")
        chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
        data_chunks = partition(1:M, chunk_size) # partition your data into chunks that
        β = zeros(Float64, length(features));
        β = ProbitRegression(β, psms[cv_fold,features], psms[cv_fold,:target], data_chunks, max_iter = max_iter);
        ModelPredictCVFold!(psms[!,:prob], psms[!,:cv_fold], UInt8(k), psms[!,features], β, allcv_data_chunks)  
        println("sum(PSMs[cv_fold,:score]) ", sum(psms[cv_fold,:score]))
    end
end

function addPostIntegrationFeatures!(psms::DataFrame, 
                                    MS_TABLE::Arrow.Table, 
                                    precursors::Vector{LibraryPrecursorIon{T}},
                                    ms_file_idx::Integer,
                                    ms_id_to_file_path::Dictionary{UInt32, String},
                                    rt_to_irt_interp::Any,
                                    prec_id_to_irt::Dictionary{UInt32, Tuple{Float64,Float32}}) where {T<:AbstractFloat}
    println("ms_file_idx test $ms_file_idx")
    filter!(x -> x.best_scan, psms);
    filter!(x->x.weight>0, psms);
    filter!(x->x.data_points>0, psms)
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
    sequence = Vector{String}(undef, N);
    stripped_sequence = Vector{String}(undef, N);
    adjusted_intensity_explained = zeros(Float16, N);
    #log2_base_peak_intensity = zeros(Float16, N);
    charge = zeros(UInt8, N)
    sequence_length = zeros(UInt8, N);
    b_y_overlap = zeros(Bool, N);
    spectrum_peak_count = zeros(Float16, N);
    prec_mz = zeros(Float32, size(psms, 1));
    Mox = zeros(UInt8, N);
    TIC = zeros(Float16, N);

    tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    precursor_idx::Vector{UInt32} = psms[!,:precursor_idx] 
    scan_idx::Vector{UInt32} = psms[!,:scan_idx]
    masses = MS_TABLE[:masses]::Arrow.List{Union{Missing, SubArray{Union{Missing, Float32}, 1, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}, Tuple{UnitRange{Int64}}, true}}, Int64, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}}
    longest_y::Vector{UInt8} = psms[!,:longest_y]
    longest_b::Vector{UInt8} = psms[!,:longest_b]
    rt::Vector{Float32} = psms[!,:RT]
    tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    log2_intensity_explained = psms[!,:log2_intensity_explained]::Vector{Float16}
    #precursor_idx = PSMs[!,:precursor_idx]::Vector{UInt32}
    function countMOX(seq::String)::UInt8
        mox = zero(UInt8)
        in_mox = false
        for aa in seq
            if in_mox
                if aa == 'x'
                    mox += one(UInt8)
                end
                in_mox = false
            end
            if aa == 'o'
                in_mox = true
            end
        end
        return mox
    end

    tasks_per_thread = 5
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin 
            for i in chunk
                prec = precursors[precursor_idx[i]]

                irt_obs[i] = rt_to_irt_interp[ms_id_to_file_path[ms_file_idx]](rt[i])
                irt_pred[i] = getIRT(prec)
                irt_diff[i] = abs(irt_obs[i] - first(prec_id_to_irt[precursor_idx[i]]))
                irt_error[i] = abs(irt_obs[i] - irt_pred[i])

                missed_cleavage[i] = prec.missed_cleavages
                sequence[i] = prec.sequence
                stripped_sequence[i] = replace(sequence[i], r"\(.*?\)" => "")#replace.(sequence[i], "M(ox)" => "M");
                Mox[i] = countMOX(prec.sequence)
                sequence_length[i] = length(stripped_sequence[i])


                TIC[i] = Float16(log2(tic[scan_idx[i]]))
                adjusted_intensity_explained[i] = Float16(log2(TIC[i]) + log2_intensity_explained[i]);
                charge[i] = getPrecCharge(prec)
                b_y_overlap[i] = ((sequence_length[i] - longest_y[i])>longest_b[i]) &  (longest_b[i] > 0) & (longest_y[i] > 0);

                spectrum_peak_count[i] = length(masses[scan_idx[i]])
         
                prec_mz[i] = precursors[precursor_idx[i]].mz::Float32;
            end
        end
    end
    fetch.(tasks)
    psms[!,:irt_obs] = irt_obs
    psms[!,:irt_pred] = irt_pred
    psms[!,:irt_diff] = irt_diff
    psms[!,:irt_error] = irt_error

    psms[!,:missed_cleavage] = missed_cleavage
    psms[!,:sequence] = sequence
    psms[!,:stripped_sequence] = stripped_sequence
    psms[!,:Mox] = Mox
    psms[!,:sequence_length] = sequence_length

    psms[!,:tic] = TIC
    psms[!,:adjusted_intensity_explained] = adjusted_intensity_explained
    psms[!,:charge] = charge
    psms[!,:b_y_overlap] = b_y_overlap
    psms[!,:spectrum_peak_count] = spectrum_peak_count
    psms[!,:prec_mz] = prec_mz
    psms[!,:ms_file_idx] .= ms_file_idx
    return nothing
end