function makeRTIndices(psms_dict::Dictionary{String, DataFrame}, 
                       precID_to_iRT::Dictionary{UInt32, Tuple{Float64, Float32}},
                       iRT_RT::Any;
                       min_prob::AbstractFloat = 0.75)
    #Maps filepath to a retentionTimeIndex (see buildRTIndex.jl)
    rt_indices = Dictionary{String, retentionTimeIndex{Float32, Float32}}()
    #iRT dict
    #iRT_dict = Dictionary{String, UnorderedDictionary{UInt32, Float32}}()
    for (file_path, psms) in ProgressBar(pairs(psms_dict)) #For each file in the experiment
        #insert!(iRT_dict, file_path, UnorderedDictionary{UInt32, Float32}())
        #Impute empirical iRT value for psms with probability lower than the threshold
        filter!(x->x.prob>=min_prob, psms); 
        RTs, mzs, prec_ids = zeros(Float32, length(precID_to_iRT)), zeros(Float32, length(precID_to_iRT)), zeros(UInt32, length(precID_to_iRT))

        prec_set = Dictionary(psms[!,:precursor_idx]::Vector{UInt32},
                            zip(psms[!,:RT]::Vector{Float32}, 
                                #psms[!,:iRT]::Vector{Float32},
                                psms[!,:prec_mz]::Vector{Float32})
                            )
        #Set of precursors where empirical iRTs can be used and do not need to be imputed 
        pset = Set(psms[!,:precursor_idx]::Vector{UInt32})
        #Loop through all precursors in the seed 
        Threads.@threads for (i, (prec_id, irt_mz_imputed)) in collect(enumerate(pairs(precID_to_iRT)))
            prec_ids[i] = prec_id
            if prec_id ∈ pset::Set{UInt32} #Use empirical iRT/RT
                rt, mz = prec_set[prec_id]::Tuple{Float32, Float32}
                RTs[i] = rt
                mzs[i] = mz
            else #Impute empirical iRT from the best observed psm for the precursor accross the experiment 
                irt, mz = irt_mz_imputed::Tuple{Float64, Float32}
                rt = iRT_RT[file_path](irt)::Float64 #Convert iRT to RT 
                RTs[i] = rt
                mzs[i] = mz
                #insert!(iRT_dict[file_path], prec_id, abs(irt - prec_set[prec_id][2]))
            end
        end
        
        #Build RT index 
        rt_df = DataFrame(Dict(:RT => RTs,
                                :prec_mz => mzs,
                                :precursor_idx => prec_ids))
        sort!(rt_df, :RT)
        rt_index = buildRTIndex(rt_df);
        insert!(rt_indices, file_path, rt_index)
    end
    return rt_indices
end

function mapRTandiRT(psms_dict::Dictionary{String, DataFrame}; min_prob::AbstractFloat = 0.9)

    #Dictionaries mapping fild_id names to data
    prec_ids = Dictionary{String, Set{UInt32}}() #File name => precursor id
    iRT_RT = Dictionary{String, Any}() #File name => KDEmapping from iRT to RT 
    RT_iRT = Dictionary{String, Any}() #File name => KDEmapping from iRT to RT 
    for (key, psms) in ProgressBar(pairs(psms_dict)) #For each data frame 
        psms[!,:file_path] .= key #Add column for file name

        insert!(prec_ids,key,Set(psms[!,:precursor_idx])) #Set of precursors ids in the dataframe

        best_hits = psms[!,:prob].>min_prob#Map RTs using only the best psms
        
        #Build RT=>iRT and iRT=> RT mappings for the file and add to the dictionaries 
        insert!(iRT_RT, key, KDEmapping(psms[best_hits,:iRT_predicted],
        psms[best_hits,:RT]
        ))
        insert!(RT_iRT, key, KDEmapping(psms[best_hits,:RT],
        psms[best_hits,:iRT_predicted]
        ))
    end
    return iRT_RT, RT_iRT
end

function getPrecIDtoiRT(psms_dict::Dictionary{String, DataFrame}, RT_iRT::Any; max_q_value::AbstractFloat = 0.1)
    psms = vcat(values(psms_dict)...); #Combine data from all files in the experiment
    #Only consider precursors with a q_value passing the threshold
    # at least once accross the entire experiment
    filter!(x->x.q_value<=max_q_value,psms); 

    #For each precursor, what was the best scan, maximum q-value, and observed iRT for the best scan
    psms[!,:best_scan] .= false #Best scan for a given precursor_idx accross the experiment
    psms[!,:iRT_observed] .= zero(Float64) #Observed iRT for each psm
    #psms[!,:iRT_observed_best] .= zero(Float64) #Best iRT accross the experiment for a given precursor
    grouped_psms = groupby(psms, :precursor_idx); #Groupby precursor idx

    for i in ProgressBar(range(1, length(grouped_psms))) #For each precursor 

        #Best scan for the precursor accross the experiment
        best_scan = argmax(grouped_psms[i].prob) 
        grouped_psms[i].best_scan[best_scan] = true #Makr columsn containing the best psm for a precursor
        #irt = RT_iRT[grouped_psms[i].file_path[best_scan]](grouped_psms[i].RT[best_scan])
        #grouped_psms[i].iRT_observed[best_scan] = irt

        #Could alternatively estimate the iRT using the best-n scans 
        for j in range(1, size(grouped_psms[i])[1])
            #Get empirical iRT for the psm given the empirical RT
            irt = RT_iRT[grouped_psms[i].file_path[j]](grouped_psms[i].RT[j])
            grouped_psms[i].iRT_observed[j] = irt
        end

        grouped_psms[i].iRT_observed[best_scan] = median(grouped_psms[i].iRT_observed)
    end

    filter!(x->x.best_scan, psms); #Filter out non-best scans 
    #Map the precursor idx to the best empirical iRT extimate and the precursor m/z
    return Dictionary(psms[!,:precursor_idx], zip(psms[!,:iRT_observed], psms[!,:prec_mz]))
end
