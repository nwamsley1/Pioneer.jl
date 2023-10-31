function refinePSMs!(PSMs::DataFrame, MS_TABLE::Arrow.Table, precursors::Vector{LibraryPrecursor{T}}; min_spectral_contrast::AbstractFloat = 0.9,  n_bins::Int = 200, granularity::Int = 50) where {T<:AbstractFloat}
    
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
    PSMs[!,:decoy] .= false;
    PSMs[!,:Mox] .= zero(UInt8);
    PSMs[!,:missed_cleavage] .= zero(UInt8);
    PSMs[!,:sequence] .= "";
    PSMs[!,:iRT] .= zero(Float32);
    PSMs[!,:RT] .= zero(Float32);
    PSMs[!,:log2_base_peak_intensity] .= zero(Float16);
    PSMs[!,:TIC] .= zero(Float16);
    PSMs[!,:adjusted_intensity_explained] .= zero(Float16);
    PSMs[!,:charge] .= zero(UInt8);
    PSMs[!,:sequence_length] .= zero(UInt8);
    PSMs[!,:total_ions] .= zero(UInt16);
    PSMs[!,:b_y_overlap] .= false;
    PSMs[!,:err_norm] .= zero(Float16);
    PSMs[!,:target] .= false
    PSMs[!,:weight_log2] .= zero(Float16)
    PSMs[!,:max_weight] .= zero(Float32)
    #SMs[!,:sequence_length] .= false
    Threads.@threads for i in ProgressBar(range(1, size(PSMs)[1]))
        PSMs[i,:decoy] = isDecoy(precursors[PSMs[i,:precursor_idx]]);
        PSMs[i,:missed_cleavage] = precursors[PSMs[i,:precursor_idx]].missed_cleavages
    #transform!(PSMs, AsTable(:) => ByRow(psm -> UInt8(length(collect(eachmatch(r"ox", psm[:sequence]))))) => [:Mox])
        PSMs[i,:sequence] = precursors[PSMs[i,:precursor_idx]].sequence;
        PSMs[i,:Mox] = UInt8(length(collect(eachmatch(r"ox",  PSMs[i,:sequence]))))
        PSMs[i,:iRT] = Float32(getIRT(precursors[PSMs[i,:precursor_idx]]));
        PSMs[i,:RT] = Float32(MS_TABLE[:retentionTime][PSMs[i,:scan_idx]]);
        PSMs[i,:TIC] = Float16(log2(MS_TABLE[:TIC][PSMs[i,:scan_idx]]));
        PSMs[i,:adjusted_intensity_explained] = Float16(log2(MS_TABLE[:TIC][PSMs[i,:scan_idx]]*(2^PSMs[i,:log2_intensity_explained])));
        PSMs[i,:charge] = UInt8(getCharge(precursors[PSMs[i,:precursor_idx]]));
        PSMs[i,:sequence_length] = UInt8(precursors[PSMs[i,:precursor_idx]].length);
        PSMs[i,:total_ions] = UInt16(PSMs[i,:y_count])+UInt16(PSMs[i,:b_count]);
        PSMs[i,:b_y_overlap] = ((PSMs[i,:sequence_length] - PSMs[i,:longest_y])>PSMs[i,:longest_b]) &  (PSMs[i,:longest_b] > 0) & (PSMs[i,:longest_y] > 0);
        PSMs[i,:err_norm] = PSMs[i,:error]/PSMs[i,:total_ions]
        PSMs[i,:target] = PSMs[i,:decoy] == false
        PSMs[i,:weight_log2] = max(Float16(0.0), Float16(log2(PSMs[i,:weight])))
        if isinf(PSMs[i,:matched_ratio])
            PSMs[i,:matched_ratio] = Float16(60000)
        end

    end
    ###########################

    ###########################
    #Estimate RT Prediction Error
    #best_psms = combine(sdf -> sdf[argmax(sdf.matched_ratio), :], groupby(PSMs[(PSMs[!,spectral_contrast].>min_spectral_contrast) .& (PSMs[!,:decoy].==false),:], [:scan_idx]))
    linear_spline = KDEmapping(PSMs[(PSMs[!,:spectral_contrast].>0.9) .& (PSMs[!,:decoy].==false) .& (PSMs[!,:entropy_score].>0.9) .& (PSMs[!,:total_ions].>6),:iRT],
                                        PSMs[(PSMs[!,:spectral_contrast].>0.9) .& (PSMs[!,:decoy].==false).& (PSMs[!,:entropy_score].>0.9) .& (PSMs[!,:total_ions].>6),:RT]
                                    )
    PSMs[:,:RT_pred] = Float16.(linear_spline(PSMs[:,:iRT]))
    PSMs[:,:RT_error] = Float16.(abs.(PSMs[!,:RT_pred] .- PSMs[!,:RT]))
    ############################

    ###########################
    #Filter on Rank and Topn
    #filter!(:best_rank => x -> x<2, PSMs)
    #filter!(:topn => x -> x>1, PSMs)
    ############################

    #isnan.(PSMs[:,:matched_ratio])
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
    #######################
    #Transform Features
    #transform!(PSMs, AsTable(:) => ByRow(psm -> Float16(log2(abs(psm[:error]/psm[:total_ions])))) => [:err_norm]);
    #transform!(PSMs, AsTable(:) => ByRow(psm -> psm[:decoy]==false) => [:target]);
    #transform!(PSMs, AsTable(:) => ByRow(psm -> Float16(log2(psm[:weight]))) => [:weight_log2]);
    #transform!(PSMs, AsTable(:) => ByRow(psm -> Float16(log2(psm[:matched_ratio]))) => [:matched_ratio]);

    #transform!(PSMs, AsTable(:) => ByRow(psm -> isinf(psm[:weight_log2]) ? zero(Float16) : psm[:weight_log2]) => [:weight_log2]);
    #transform!(PSMs, AsTable(:) => ByRow(psm -> isinf(psm[:matched_ratio]) ? zero(Float16) : psm[:matched_ratio]) => [:matched_ratio]);

    #transform!(PSMs, AsTable(:) => ByRow(psm -> isnan(psm[:matched_ratio]) ? zero(Float16) : psm[:matched_ratio]) => [:matched_ratio]);


    #min_log2_weight = minimum(PSMs[isinf.(PSMs[!,:weight_log2]),:weight_log2])
    #replace!(PSMs[:,:weight_log2], -Inf => min_log2_weight);
    #replace!(PSMs[:,:scribe_score], Inf => minimum(PSMs[PSMs[:,:scribe_score].!=Inf,:scribe_score]));

    #PSMs[isinf.(PSMs[:,:weight_log2]),:weight_log2] .= zero(Float16);
    #PSMs[:,:matched_ratio_log2] = Float16.(log2.(PSMs[:,:matched_ratio]))
    #PSMs[isinf.(PSMs[:,:matched_ratio_log2]),:matched_ratio_log2] .= zero(Float16);
    #filter!(:entropy_sim => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), PSMs);
    ########################
    #Rough Target-Decoy discrimination
    
    #model_fit = glm(@formula(target ~ entropy_sim + poisson + hyperscore +
    #                            scribe_score + weight_log2 + topn + spectral_contrast + 
    #                            RT_error + missed_cleavage + Mox + intensity_explained + error + total_ions), PSMs, 
    #                            Binomial(), 
    #                           ProbitLink())
    PSMs[:,:q_value] .= zero(Float16);
    FORM = FormulaTerm(
        (Term(:target),),
        (Term(:entropy_score),
         Term(:scribe_corrected),
         Term(:weight_log2),
         Term(:spectral_contrast_corrected),
         Term(:RT_error),
         Term(:missed_cleavage),
         Term(:Mox),
         Term(:TIC),
         Term(:total_ions),
         Term(:poisson),
         Term(:error),
         Term(:b_y_overlap),
         Term(:charge)
    ))
    model_fit = glm(FORM, PSMs, 
                            Binomial(), 
                            ProbitLink())
    Y′ = Float16.(GLM.predict(model_fit, PSMs));
    getQvalues!(PSMs, allowmissing(Y′),  allowmissing(PSMs[:,:decoy]));
    println("Target PSMs at 1% FDR: ", sum((PSMs.q_value.<=0.01).&(PSMs.decoy.==false)))
    PSMs[:,:prob] = allowmissing(Y′);

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