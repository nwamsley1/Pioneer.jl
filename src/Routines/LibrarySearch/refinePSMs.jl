function refinePSMs!(PSMs::DataFrame, MS_TABLE::Arrow.Table, precursors::Vector{LibraryPrecursor{T}}; max_rt_error::Float64,  n_bins::Int = 200, granularity::Int = 50) where {T<:AbstractFloat}
    
    ###########################
    #Correct Weights by base-peak intensity
    #=
    for i in size(PSMs)[1]
        prec_id = PSMs[i,:precursor_idx]
        base_peak_intensity = precursors[prec_id].base_peak_intensity
        PSMs[i,:weight] = Float32(PSMs[i,:weight]/base_peak_intensity)
    end
    =#
    ###########################
    #Allocate new columns
    @time begin

    N = size(PSMs, 1)
    decoys = zeros(Bool, N);
    missed_cleavage = zeros(UInt8, N);
    Mox = zeros(UInt8, N);
    iRT_pred = zeros(Float32, N);
    RT = zeros(Float32, N);
    TIC = zeros(Float16, N);
    adjusted_intensity_explained = zeros(Float16, N);
    charge = zeros(UInt8, N);
    total_ions = zeros(UInt16, N);
    err_norm = zeros(Float16, N);
    targets = zeros(Bool, N);
    spectrum_peak_count = zeros(UInt32, N);
    scan_idx::Vector{UInt32} = PSMs[!,:scan_idx]
    precursor_idx::Vector{UInt32} = PSMs[!,:precursor_idx]
    log2_intensity_explained::Vector{Float16} = PSMs[!,:log2_intensity_explained]
    y_count::Vector{UInt8} = PSMs[!,:y_count]
    b_count::Vector{UInt8} = PSMs[!,:b_count]
    error::Vector{Float32} = PSMs[!,:error]
    matched_ratio::Vector{Float16} = PSMs[!,:matched_ratio]
    tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    scan_retention_time = MS_TABLE[:retentionTime]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    masses = MS_TABLE[:masses]::Arrow.List{Union{Missing, SubArray{Union{Missing, Float32}, 1, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}, Tuple{UnitRange{Int64}}, true}}, Int64, Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}}
    #PSMs[!,:total_ions]
    #SMs[!,:sequence_length] .= false
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

    Threads.@threads for i in ProgressBar(range(1, size(PSMs)[1]))
        decoys[i] = isDecoy(precursors[precursor_idx[i]]);
        targets[i] = decoys[i] == false
        missed_cleavage[i] = precursors[precursor_idx[i]].missed_cleavages
        #sequence[i] = precursors[precursor_idx[i]].sequence;
        Mox[i] = countMOX(precursors[precursor_idx[i]].sequence)::UInt8 #UInt8(length(collect(eachmatch(r"ox",  precursors[precursor_idx[i]].sequence))))
        iRT_pred[i] = Float32(getIRT(precursors[precursor_idx[i]]));
        RT[i] = Float32(scan_retention_time[scan_idx[i]]);
        #PSMs[i,:iRT_obdserved] = RT_iRT[PSMs[i,:file_path]](PSMs[i,:RT])
        TIC[i] = Float16(log2(tic[scan_idx[i]]));
        adjusted_intensity_explained[i] = Float16(min(log2(tic[scan_idx[i]]*(2^log2_intensity_explained[i])), 
                                                        6e4));
        charge[i] = UInt8(getCharge(precursors[precursor_idx[i]]));
        total_ions[i] = UInt16(y_count[i] + b_count[i]);
        err_norm[i] = min(Float16((error[i])/(total_ions[i])), 6e4)
        spectrum_peak_count[i] = UInt32(length(masses[scan_idx[i]]))
        matched_ratio[i] = Float16(min(matched_ratio[i], 6e4))
        #stripped_sequence[i] = replace.(sequence[i], "M(ox)" => "M");
    end
    PSMs[!,:matched_ratio] = matched_ratio
    PSMs[!,:decoy] = decoys
    PSMs[!,:iRT_predicted] = iRT_pred
    PSMs[!,:RT] = RT
    PSMs[!,:TIC] = TIC
    PSMs[!,:total_ions] = total_ions
    PSMs[!,:err_norm] = err_norm
    PSMs[!,:target] = targets
    PSMs[!,:missed_cleavage] = missed_cleavage
    PSMs[!,:Mox] = Mox
    PSMs[!,:adjusted_intensity_explained] = adjusted_intensity_explained
    PSMs[!,:charge] = charge
    PSMs[!,:spectrum_peak_count] = spectrum_peak_count
    end
    ###########################

    ###########################
    #Estimate RT Prediction Error
    @time begin
    #best_psms = combine(sdf -> sdf[argmax(sdf.matched_ratio), :], groupby(PSMs[(PSMs[!,spectral_contrast].>min_spectral_contrast) .& (PSMs[!,:decoy].==false),:], [:scan_idx]))
    linear_spline = KDEmapping(PSMs[(PSMs[!,:spectral_contrast].>0.9) .& (PSMs[!,:decoy].==false) .& (PSMs[!,:entropy_score].>0.9) .& (PSMs[!,:total_ions].>6),:RT],
                                        PSMs[(PSMs[!,:spectral_contrast].>0.9) .& (PSMs[!,:decoy].==false).& (PSMs[!,:entropy_score].>0.9) .& (PSMs[!,:total_ions].>6),:iRT_predicted]
                                    )
    #PSMs[!,:RT_pred] = 
    PSMs[:,:iRT_observed] = Float16.(linear_spline(PSMs[:,:RT]))
    PSMs[!,:iRT_error] = Float16.(abs.(PSMs[!,:iRT_observed] .- PSMs[!,:iRT_predicted]))
    end
    ############################
    #return 
    ###########################
    #Filter on Rank and Topn
    #filter!(x->x.RT_error<max_rt_error, PSMs);
    #filter!(x->isnan(x.entropy_score)==false, PSMs);
    #filter!(:best_rank => x -> x<2, PSMs);
    #filter!(:topn => x -> x>1, PSMs);
    ############################

    #######################
    #Clean Features
    #######################
    @time begin
    PSMs[:,:q_value] .= zero(Float16);
    
    FORM = FormulaTerm(
        (Term(:target),),
        (
         #Term(:scribe_corrected),
         #Term(:spectral_contrast_corrected),
         Term(:spectral_contrast),
         Term(:scribe),
         Term(:iRT_error),
         Term(:missed_cleavage),
         Term(:Mox),
         Term(:TIC),
         Term(:city_block),
         Term(:entropy_score),
         Term(:total_ions),
         Term(:err_norm),
         Term(:charge),
         Term(:spectrum_peak_count),
         #Term(:poisson)
         #Term(:adjusted_intensity_explained)
    ))

    model_fit = glm(FORM, PSMs[!,[:target,:spectral_contrast,:scribe,:iRT_error,:missed_cleavage,:Mox,:TIC,:city_block,:entropy_score,:total_ions,:err_norm,:charge,:spectrum_peak_count]], 
                            Binomial(), 
                            ProbitLink())
    Y′ = Float16.(GLM.predict(model_fit, PSMs[!,[:target,:spectral_contrast,:scribe,:iRT_error,:missed_cleavage,:Mox,:TIC,:city_block,:entropy_score,:total_ions,:err_norm,:charge,:spectrum_peak_count]]));
    getQvalues!(PSMs, allowmissing(Y′),  allowmissing(PSMs[:,:decoy]));

    println("Target PSMs at 25% FDR: ", sum((PSMs.q_value.<=0.25).&(PSMs.decoy.==false)))
    println("Target PSMs at 10% FDR: ", sum((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)))
    println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))
    end
    PSMs[:,:prob] = allowmissing(Y′);
    
    ########################
    #Add columns 
    PSMs[!,:best_scan] .= false
    
    return 
end

function _refinePSMs!(PSMs::DataFrame, MS_TABLE::Arrow.Table, precursors::Vector{LibraryPrecursor{T}}, RT_iRT::Any; min_spectral_contrast::AbstractFloat = 0.9,  n_bins::Int = 200, granularity::Int = 50) where {T<:AbstractFloat}
    
    ###########################
    #Correct Weights by base-peak intensity
    #=
    for i in size(PSMs)[1]
        prec_id = PSMs[i,:precursor_idx]
        base_peak_intensity = precursors[prec_id].base_peak_intensity
        PSMs[i,:weight] = Float32(PSMs[i,:weight]/base_peak_intensity)
    end
    =#
    ###########################
    #Allocate new columns
   
    test = @time begin
    #Threads.@threads for i in ProgressBar(range(1, size(PSMs)[1]))
    N = size(PSMs)[1]
    decoys = zeros(Bool, N);
    iRT_pred = zeros(Float32, N);
    iRT_obs = zeros(Float32, N);
    RT = zeros(Float32, N);
    TIC = zeros(Float16, N);
    total_ions = zeros(UInt16, N);
    err_norm = zeros(Float16, N);
    targets = zeros(Bool, N);
    charge = zeros(UInt8, N);
    prec_mz = zeros(Float32, N);
    scan_idx::Vector{UInt32} = PSMs[!,:scan_idx]
    precursor_idx::Vector{UInt32} = PSMs[!,:precursor_idx]
    y_count::Vector{UInt8} = PSMs[!,:y_count]
    b_count::Vector{UInt8} = PSMs[!,:b_count]
    error::Vector{Float32} = PSMs[!,:error]
    #PSMs[!,:total_ions]
    tic = MS_TABLE[:TIC]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    scan_retention_time = MS_TABLE[:retentionTime]::Arrow.Primitive{Union{Missing, Float32}, Vector{Float32}}
    matched_ratio::Vector{Float16} = PSMs[!,:matched_ratio]
    #file_path::Vector{String} = PSMs[!,:file_path]
    #err_norm = PSMs[!,:err_norm]
    Threads.@threads for i in range(1, size(PSMs)[1])
        decoys[i] = isDecoy(precursors[precursor_idx[i]]);
        targets[i] = decoys[i] == false
        #iRT_pred[i] = Float32(getIRT(precursors[precursor_idx[i]]));
        RT[i] = Float32(scan_retention_time[scan_idx[i]]);
        #iRT_obs = RT_iRT[PSMs[i,:file_path]](PSMs[i,:RT])
        charge[i] = UInt8(getCharge(precursors[precursor_idx[i]]));
        prec_mz[i] = Float32(getMz(precursors[precursor_idx[i]]));
        TIC[i] = Float16(log2(tic[scan_idx[i]]));
        total_ions[i] = UInt16(y_count[i] + b_count[i]);
        err_norm[i] = min(Float16((error[i])/(total_ions[i])), 6e4)
        if isinf(matched_ratio[i])
            matched_ratio[i] = Float16(60000)*sign(matched_ratio[i])
        end
        #if isinf(PSMs[i,:err_norm])
        #    PSMs[i,:err_norm] = Float16(60000)*sign(PSMs[i,:err_norm])
        #end
    end
    PSMs[!,:matched_ratio] = matched_ratio
    PSMs[!,:decoy] = decoys
    PSMs[!,:iRT_predicted] = iRT_pred
    PSMs[!,:RT] = RT
    PSMs[!,:TIC] = TIC
    PSMs[!,:total_ions] = total_ions
    PSMs[!,:err_norm] = err_norm
    PSMs[!,:target] = targets
    PSMs[!,:charge] = charge
    PSMs[!,:prec_mz] = prec_mz
    end
    #println("test $test")
    ###########################

    ###########################
    #Estimate RT Prediction Error
    #best_psms = combine(sdf -> sdf[argmax(sdf.matched_ratio), :], groupby(PSMs[(PSMs[!,spectral_contrast].>min_spectral_contrast) .& (PSMs[!,:decoy].==false),:], [:scan_idx]))
    #linear_spline = KDEmapping(PSMs[(PSMs[!,:spectral_contrast].>0.9) .& (PSMs[!,:decoy].==false) .& (PSMs[!,:entropy_score].>0.9) .& (PSMs[!,:total_ions].>6),:iRT_predicted],
    #                                   PSMs[(PSMs[!,:spectral_contrast].>0.9) .& (PSMs[!,:decoy].==false).& (PSMs[!,:entropy_score].>0.9) .& (PSMs[!,:total_ions].>6),:RT]
    #                               )
    #PSMs[!,:RT_pred] = Float16.(linear_spline(PSMs[!,:iRT_predicted]))
    #PSMs[!,:RT_error] = Float16.(abs.(PSMs[!,:RT_pred] .- PSMs[!,:RT]))
    ############################
     #######################
    #Clean Features
    #######################
    #model_time = @time begin
    #@time begin
    #=
    PSMs[:,:q_value] .= zero(Float16);
    FORM = FormulaTerm(
        (Term(:target),),
        (
         #Term(:scribe_corrected),
         #Term(:spectral_contrast_corrected),
         Term(:spectral_contrast),
         Term(:scribe),
         #Term(:RT_error),
         Term(:TIC),
         Term(:city_block),
         Term(:entropy_score),
         Term(:total_ions),
         Term(:err_norm),
         Term(:city_block_fitted),
         #Term(:scribe_fitted)
    ));
    end
    #filter!(x->isnan(x.entropy_score)==false, PSMs);
    @time model_fit = glm(FORM, PSMs, 
                            Binomial(), 
                            ProbitLink());
    @time Y′ = Float16.(GLM.predict(model_fit, PSMs));
    @time getQvalues!(PSMs, allowmissing(Y′),  allowmissing(PSMs[:,:decoy]));
    #println("Target PSMs at 25% FDR: ", sum((PSMs.q_value.<=0.25).&(PSMs.decoy.==false)));
    #println("Target PSMs at 10% FDR: ", sum((PSMs.q_value.<=0.1).&(PSMs.decoy.==false)));
    #println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)));
    PSMs[!,:prob] = allowmissing(Y′);
    =#
    #end
    #println("model_time $model_time")
    ########################
    #Add columns 
    new_features = @time begin
    new_cols = [(:peak_area,                   Union{Float32, Missing})
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
    (:max_weight,           Union{Float16, Missing})
    (:max_spectral_contrast,   Union{Float16, Missing})
    (:max_entropy,              Union{Float16, Missing})
    (:max_scribe_score,     Union{Float16, Missing})
    (:max_scribe_fitted,     Union{Float16, Missing})
    (:max_city_fitted,     Union{Float16, Missing})
    (:mean_city_fitted,     Union{Float16, Missing})
    (:ions_sum,                 Union{UInt32, Missing})
    (:max_ions,                 Union{UInt16, Missing})
    (:data_points,              Union{UInt32, Missing})
    (:fraction_censored,              Union{Float16, Missing})
    (:max_matched_ratio,       Union{Float32, Missing})
    (:base_width_min,           Union{Float16, Missing})
    (:best_scan, Union{Bool, Missing})];

    for column in new_cols
        col_type = last(column);
        col_name = first(column)
        PSMs[!,col_name] = zeros(col_type, size(PSMs)[1])
    end
    PSMs[!,:best_scan] .= false;
    end
    #println("new_features $new_features")
    return 
end

function addFeatures!(PSMs::DataFrame, MS_TABLE::Arrow.Table, precursors::Vector{LibraryPrecursor{T}}) where {T<:AbstractFloat}
    
    ###########################
    #Allocate new columns
    #println("TEST")
    PSMs[!,:missed_cleavage] .= zero(UInt8);
    PSMs[!,:sequence] .= "";
    PSMs[!,:log2_base_peak_intensity] .= zero(Float16);
    PSMs[!,:adjusted_intensity_explained] .= zero(Float16);
    PSMs[!,:charge] .= zero(UInt8);
    PSMs[!,:sequence_length] .= zero(UInt8);
    PSMs[!,:b_y_overlap] .= false;
    PSMs[!,:weight_log2] .= zero(Float16)
    PSMs[!,:max_weight] .= zero(Float32)
    PSMs[!,:stripped_sequence] .= "";
    PSMs[!,:spectrum_peak_count] .= zero(UInt16);
    PSMs[!,:sequence_length] .= zero(UInt8);
    PSMs[!,:Mox] .= zero(UInt8);
    @Threads.threads for i in ProgressBar(range(1, size(PSMs)[1]))
        PSMs[i,:missed_cleavage] = precursors[PSMs[i,:precursor_idx]].missed_cleavages
    #transform!(PSMs, AsTable(:) => ByRow(psm -> UInt8(length(collect(eachmatch(r"ox", psm[:sequence]))))) => [:Mox])
        PSMs[i,:sequence] = precursors[PSMs[i,:precursor_idx]].sequence;
        PSMs[i,:Mox] = UInt8(length(collect(eachmatch(r"ox",  PSMs[i,:sequence]))))
        PSMs[i,:adjusted_intensity_explained] = Float16(log2(MS_TABLE[:TIC][PSMs[i,:scan_idx]]*(2^PSMs[i,:log2_intensity_explained])));
        PSMs[i,:charge] = UInt8(getCharge(precursors[PSMs[i,:precursor_idx]]));
        #println(precursors[PSMs[i,:precursor_idx]].length)
        PSMs[i,:sequence_length] = UInt8(precursors[PSMs[i,:precursor_idx]].length);
        PSMs[i,:b_y_overlap] = ((PSMs[i,:sequence_length] - PSMs[i,:longest_y])>PSMs[i,:longest_b]) &  (PSMs[i,:longest_b] > 0) & (PSMs[i,:longest_y] > 0);
        PSMs[i,:weight_log2] = max(Float16(0.0), Float16(log2(PSMs[i,:weight])))
        PSMs[i,:spectrum_peak_count] = length(MS_TABLE[:masses][PSMs[i,:scan_idx]])
        if isinf(PSMs[i,:adjusted_intensity_explained])
            PSMs[i,:adjusted_intensity_explained] = Float16(6000)
        end
        #if isinf(coalesce(PSMs[i,:mean_log_probability], 0.0))
        #    PSMs[i,:mean_log_probability] = Float16(6000)*sign(PSMs[i,:mean_log_probability])
        #end

        PSMs[i,:stripped_sequence] = replace.(PSMs[i,:sequence], "M(ox)" => "M");
    end
    #######################
    #Clean Features
    #######################    
    return 
end



#=
function refinePSMs!(PSMs::DataFrame, MS_TABLE::Arrow.Table, precursors::Vector{LibraryPrecursor{T}}; min_spectral_contrast::AbstractFloat = 0.9,  n_bins::Int = 200, granularity::Int = 50) where {T<:AbstractFloat}
    
    ###########################
    #Correct Weights by base-peak intensity
    for i in size(PSMs)[1]
        prec_id = PSMs[i,:precursor_idx]
        base_peak_intensity = precursors[prec_id].base_peak_intensity
        PSMs[i,:weight] = Float32(PSMs[i,:weight]/base_peak_intensity)
    end
    ###########################
    #Get Precursor Features
    transform!(PSMs, AsTable(:) => ByRow(psm -> isDecoy(precursors[psm[:precursor_idx]])) => :decoy)
    transform!(PSMs, AsTable(:) => ByRow(psm -> precursors[psm[:precursor_idx]].missed_cleavages) => :missed_cleavage)
    transform!(PSMs, AsTable(:) => ByRow(psm -> precursors[psm[:precursor_idx]].sequence) => :sequence)
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float32(getIRT(precursors[psm[:precursor_idx]]))) => :iRT)
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float32(MS_TABLE[:retentionTime][psm[:scan_idx]])) => :RT)
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float16(log2(maximum(MS_TABLE[:intensities][psm[:scan_idx]])))) => :log2_base_peak_intensity)
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float16(log2(sum(MS_TABLE[:intensities][psm[:scan_idx]])))) => :TIC)
    PSMs[:,:adjusted_intensity_explained] = Float16.(log2.(PSMs[!,:intensity_explained] .* (2 .^ Float32.(PSMs[!,:TIC]))))
    transform!(PSMs, AsTable(:) => ByRow(psm -> UInt8(getCharge(precursors[psm[:precursor_idx]]))) => :charge)

    ###########################

    ###########################
    #Estimate RT Prediction Error
    #best_psms = combine(sdf -> sdf[argmax(sdf.matched_ratio), :], groupby(PSMs[(PSMs[!,spectral_contrast].>min_spectral_contrast) .& (PSMs[!,:decoy].==false),:], [:scan_idx]))
    println("TEST")
    linear_spline = KDEmapping(PSMs[(PSMs[!,:spectral_contrast].>0.9) .& (PSMs[!,:decoy].==false).& (PSMs[!,:entropy_sim].>0.9),:iRT],
                                        PSMs[(PSMs[!,:spectral_contrast].>0.9) .& (PSMs[!,:decoy].==false).& (PSMs[!,:entropy_sim].>0.9),:RT]
                                    )
                                    println("TEST2")
    PSMs[:,:RT_pred] = Float16.(linear_spline(PSMs[:,:iRT]))
    PSMs[:,:RT_error] = Float16.(abs.(PSMs[!,:RT_pred] .- PSMs[!,:RT]))
    ############################

    ###########################
    #Filter on Rank and Topn
    #filter!(:best_rank => x -> x<2, PSMs)
    #filter!(:topn => x -> x>1, PSMs)
    ############################


    # Number of PSMs occuring for each precursor 
    PSMs[:,:n_obs] .= zero(UInt16)
    #sort!(PSMs, [:precursor_idx]);
    #grouped_df = groupby(PSMs, :precursor_idx);
    #PSMs[:,:n_obs] = (combine(grouped_df) do sub_df
    #    repeat([size(sub_df)[1]], size(sub_df)[1])
    #end)[:,:x1]
    #grouped_df = nothing

    #######################
    #Clean Features
    PSMs[isnan.(PSMs[:,:matched_ratio]),:matched_ratio] .= Inf;
    PSMs[(PSMs[:,:matched_ratio]).==Inf,:matched_ratio] .= maximum(PSMs[(PSMs[:,:matched_ratio]).!=Inf,:matched_ratio]);
    #min_city_block = minimum(PSMs[isinf.(PSMs[!,:city_block]),:city_block])
    #replace!(PSMs[:,:city_block], -Inf => min_city_block);

    #min_scribe = minimum(PSMs[isinf.(PSMs[!,:scribe_score]),:scribe_score])
    #replace!(PSMs[:,:scribe_score], Inf => min_scribe);

    #transform!(PSMs, AsTable(:) => ByRow(psm -> isinf(psm[:city_block]) ? zero(Float16) : psm[:city_block] => [:weight_log2]));
    #transform!(PSMs, AsTable(:) => ByRow(psm -> isinf(psm[:matched_ratio]) ? zero(Float16) : psm[:matched_ratio] => [:matched_ratio]));


    transform!(PSMs, AsTable(:) => ByRow(psm -> UInt8(length(collect(eachmatch(r"ox", psm[:sequence]))))) => [:Mox]);

    #######################
    #Transform Features
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float16(log2(abs(psm[:error]/psm[:total_ions])))) => [:err_norm_log2]);
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float16(log2(abs(psm[:error])))) => [:error]);
    transform!(PSMs, AsTable(:) => ByRow(psm -> psm[:decoy]==false) => [:target]);
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float16(log2(psm[:weight]))) => [:weight_log2]);
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float16(log2(psm[:matched_ratio]))) => [:matched_ratio]);

    transform!(PSMs, AsTable(:) => ByRow(psm -> isinf(psm[:weight_log2]) ? zero(Float16) : psm[:weight_log2]) => [:weight_log2]);
    transform!(PSMs, AsTable(:) => ByRow(psm -> isinf(psm[:matched_ratio]) ? zero(Float16) : psm[:matched_ratio]) => [:matched_ratio]);

    transform!(PSMs, AsTable(:) => ByRow(psm -> isnan(psm[:matched_ratio]) ? zero(Float16) : psm[:matched_ratio]) => [:matched_ratio]);


    #min_log2_weight = minimum(PSMs[isinf.(PSMs[!,:weight_log2]),:weight_log2])
    #replace!(PSMs[:,:weight_log2], -Inf => min_log2_weight);
    #replace!(PSMs[:,:scribe_score], Inf => minimum(PSMs[PSMs[:,:scribe_score].!=Inf,:scribe_score]));

    #PSMs[isinf.(PSMs[:,:weight_log2]),:weight_log2] .= zero(Float16);
    #PSMs[:,:matched_ratio_log2] = Float16.(log2.(PSMs[:,:matched_ratio]))
    #PSMs[isinf.(PSMs[:,:matched_ratio_log2]),:matched_ratio_log2] .= zero(Float16);
    #filter!(:entropy_sim => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), PSMs);
    ########################
    #Rough Target-Decoy discrimination
    PSMs[:,:q_value] .= zero(Float16)
    #model_fit = glm(@formula(target ~ entropy_sim + poisson + hyperscore +
    #                            scribe_score + weight_log2 + topn + spectral_contrast + 
    #                            RT_error + missed_cleavage + Mox + intensity_explained + error + total_ions), PSMs, 
    #                            Binomial(), 
    #                           ProbitLink())
    model_fit = glm(@formula(target ~ entropy_sim +
                            scribe_score + weight_log2 + spectral_contrast + RT_error + missed_cleavage + Mox + TIC), PSMs, 
                            Binomial(), 
                            ProbitLink())
    Y′ = Float16.(GLM.predict(model_fit, PSMs));
    getQvalues!(PSMs, allowmissing(Y′),  allowmissing(PSMs[:,:decoy]));
    println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))
    PSMs[:,:prob] = allowmissing(Y′)

    ########################
    #Add columns 
    new_cols = [(:peak_area,                   Union{Float32, Missing})
    (:GOF,                      Union{Float16, Missing})
    (:FWHM,                     Union{Float16, Missing})
    (:FWHM_01,                  Union{Float16, Missing})
    (:points_above_FWHM,        Union{UInt16, Missing})
    (:points_above_FWHM_01,     Union{UInt16, Missing})
    (:σ,                        Union{Float32, Missing})
    (:tᵣ,                       Union{Float16, Missing})
    (:τ,                        Union{Float32, Missing})
    (:H,                        Union{Float32, Missing})
    (:log_sum_of_weights,           Union{Float16, Missing})
    (:mean_log_spectral_contrast,   Union{Float16, Missing})
    (:mean_log_entropy,              Union{Float16, Missing})
    (:mean_log_probability,     Union{Float16, Missing})
    (:mean_scribe_score,     Union{Float16, Missing})
    (:ions_sum,                 Union{UInt32, Missing})
    (:data_points,              Union{UInt32, Missing})
    (:mean_matched_ratio,       Union{Float32, Missing})
    (:base_width_min,           Union{Float16, Missing})
    (:best_scan, Union{Bool, Missing})]
    for column in new_cols
        col_type = last(column)
        PSMs[!,first(column)] = Vector{col_type}(undef, size(PSMs)[1])
    end
    PSMs[!,:best_scan] .= false
    return 
end
=#

#sum(psms_counts[:,:nrow].>2)
#=
function rtSpline(X::Vector{T}, Y::Vector{T}; n_bins::Int = 200, granularity::Int = 50) where {T<:AbstractFloat}
    sort_order = sortperm(X)

    #Divide RT space into estimation bins
    est_bins = [Int(bin÷1) for bin in range(1, length = n_bins, stop = length(sort_order))]

    #x and y values for each RT estimation bin
    xs = Vector{T}(undef, length(est_bins) - 1)
    ys = Vector{T}(undef, length(est_bins) - 1)
    for i in 1:(length(est_bins) - 1)

        #RTs for the i'th estimation bin
        obs = X[sort_order[est_bins[i]:est_bins[i + 1]]]
        x = Vector(LinRange(minimum(obs), maximum(obs), granularity))
        kde = KDEUniv(ContinuousDim(), 3.0, obs, MultiKDE.gaussian)
        y = [MultiKDE.pdf(kde, _x, keep_all=false) for _x in x]
        xs[i] = x[argmax(y)]
        ys[i] = mean(Y[sort_order[est_bins[i]]])#x[sort_order[bins[i]:bins[i + 1]]])
    end
    return LinearInterpolation(xs, ys, extrapolation_bc = Line() )
end

function refinePSMs!(PSMs::DataFrame, precursors::Vector{LibraryPrecursor{T}}; min_spectral_contrast::AbstractFloat = 0.9,  n_bins::Int = 200, granularity::Int = 50) where {T<:AbstractFloat}
    transform!(PSMs, AsTable(:) => ByRow(psm -> isDecoy(precursors[psm[:precursor_idx]])) => :decoy)
    transform!(PSMs, AsTable(:) => ByRow(psm -> precursors[psm[:precursor_idx]].missed_cleavages) => :missed_cleavage)
    transform!(PSMs, AsTable(:) => ByRow(psm -> precursors[psm[:precursor_idx]].sequence) => :sequence)
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float64(getIRT(precursors[psm[:precursor_idx]]))) => :iRT)
    transform!(PSMs, AsTable(:) => ByRow(psm -> Float64(MS_TABLE[:retentionTime][psm[:scan_idx]])) => :RT)
    transform!(PSMs, AsTable(:) => ByRow(psm -> psm[:weight] < 10.0) => :nmf)
    transform!(PSMs, AsTable(:) => ByRow(psm -> getCharge(precursors[psm[:precursor_idx]])) => :charge)

    best_psms = combine(sdf -> sdf[argmax(sdf.matched_ratio), :], groupby(PSMs[(PSMs[:,:spectral_contrast_all].>min_spectral_contrast) .& (PSMs[:,:decoy].==false),:], [:scan_idx]))
    @time linear_spline = KDEmapping(best_psms[:,:iRT], best_psms[:,:RT])
    PSMs[:,:RT_pred] = linear_spline(PSMs[:,:iRT])
    PSMs[:,:RT_error] = abs.(PSMs[:,:RT_pred] .- PSMs[:,:RT])

    sort!(PSMs, [:scan_idx, :total_ions]);
    # Group DataFrame by "day" column
    grouped_df = groupby(PSMs, :scan_idx);

    #PSMs[:,:next_best] = Vector{Union{Missing, UInt32}}(undef, size(PSMs)[1])
    PSMs[:,:next_best] = (combine(grouped_df) do sub_df
        pushfirst!(diff(sub_df.total_ions), zero(UInt32))
    end)[:,:x1]

    PSMs[:,:diff_hyper] = (combine(grouped_df) do sub_df
        sort!(sub_df, :hyperscore)
        pushfirst!(diff(sub_df.hyperscore), zero(Float64))
    end)[:,:x1]

    PSMs[:,:rank_hyper] = (combine(grouped_df) do sub_df
        StatsBase.ordinalrank(sub_df.hyperscore)
    end)[:,:x1]

    PSMs[:,:rank_scribe] = (combine(grouped_df) do sub_df
        StatsBase.ordinalrank(sub_df.scribe_score)
    end)[:,:x1]

    PSMs[:,:rank_poisson] = (combine(grouped_df) do sub_df
        StatsBase.ordinalrank(sub_df.poisson)
    end)[:,:x1]

    PSMs[:,:rank_total] = (combine(grouped_df) do sub_df
        StatsBase.ordinalrank(sub_df.total_ions)
    end)[:,:x1]

    PSMs[:,:diff_scribe] = (combine(grouped_df) do sub_df
        sort!(sub_df, :scribe_score)
        pushfirst!(diff(sub_df.scribe_score), zero(Float64))
    end)[:,:x1]

    PSMs[:,:median_ions] = (combine(grouped_df) do sub_df
        repeat([median(sub_df.total_ions)], size(sub_df)[1])
    end)[:,:x1]

    grouped_df = groupby(PSMs, :precursor_idx);

    PSMs[:,:n_obs] = (combine(grouped_df) do sub_df
        repeat([size(sub_df)[1]], size(sub_df)[1])
    end)[:,:x1]

    #return linear_spline
end

=#