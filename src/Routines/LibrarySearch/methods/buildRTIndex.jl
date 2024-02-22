struct rtIndexBin{T,U<:AbstractFloat}
    lb::T
    ub::T
    prec::Vector{Tuple{UInt32, U}}
end

function compare_lb(rb::rtIndexBin{T,U}) where {T,U<:AbstractFloat}
    return rb.lb
end

getLB(rb::rtIndexBin{T,U}) where {T,U<:AbstractFloat} = rb.lb
getMZ(rb::rtIndexBin{T, U}) where {T,U<:AbstractFloat} = last(rb.prec)
getPrecID(rb::rtIndexBin{T, U}) where {T,U<:AbstractFloat} = first(rb.prec)

struct retentionTimeIndex{T,U<:AbstractFloat}
    rt_bins::Vector{rtIndexBin{T, U}}
end

function retentionTimeIndex(T::DataType, U::DataType) 
    return retentionTimeIndex(Vector{rtIndexBin{T, U}}())
end

function buildRTIndex(RTs::Vector{T}, prec_mzs::Vector{U}, prec_ids::Vector{I}, bin_rt_size::AbstractFloat) where {T,U<:AbstractFloat,I<:Integer}
    
    start_idx = 1
    start_RT =  RTs[start_idx]
    rt_index = retentionTimeIndex(T, U) #Initialize retention time index
    i = 1
    while i < length(RTs) + 1
        if ((RTs[min(i + 1, length(RTs))] - start_RT) > bin_rt_size) | (i == length(RTs))
            push!(rt_index.rt_bins, 
                    rtIndexBin(RTs[start_idx], #Retention time for first precursor in the bin
                          RTs[i],     #Retention time for last precursor in the bin
                        [(zero(UInt32), zero(Float32)) for _ in 1:(i - start_idx + 1)] #Pre-allocate precursors 
                        )
                )

            n = 1 #n'th precursor 
            for idx in start_idx:(min(i, length(RTs))) 
                rt_index.rt_bins[end].prec[n] = (prec_ids[idx], prec_mzs[idx]) #Add n'th precursor
                n += 1
            end

            sort!(rt_index.rt_bins[end].prec, by = x->last(x)) #Sort precursors by m/z
            i += 1
            start_idx = i
            start_RT = RTs[min(start_idx, length(RTs))]
            continue
        else
            i += 1
        end
    end


    function sortRTBins!(rt_index::retentionTimeIndex{T, U})
        for i in 1:length(rt_index.rt_bins)
            sort!(rt_index.rt_bins[i].prec, by = x->last(x));
        end
        return nothing
    end
    sortRTBins!(rt_index)
    return rt_index
end

buildRTIndex(PSMs::DataFrame, bin_rt_size::AbstractFloat = 0.1) = buildRTIndex(PSMs[:,:RT], PSMs[:,:prec_mz], PSMs[:,:precursor_idx], bin_rt_size)

function makeRTIndices(psms_dict::Dictionary{String, DataFrame}, 
                       precID_to_iRT::Dictionary{UInt32, Tuple{Float64, Float32}},
                       iRT_RT::Any;
                       min_prob::AbstractFloat = 0.5)
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
    #prec_ids = Dictionary{String, Set{UInt32}}() #File name => precursor id
    iRT_RT = Dictionary{String, Any}() #File name => KDEmapping from iRT to RT 
    RT_iRT = Dictionary{String, Any}() #File name => KDEmapping from iRT to RT 
    for (key, psms) in ProgressBar(pairs(psms_dict)) #For each data frame 
        psms[!,:file_path] .= key #Add column for file name

        #insert!(prec_ids,key,Set(psms[!,:precursor_idx])) #Set of precursors ids in the dataframe

        best_hits = psms[!,:q_value].<0.01#.>min_prob#Map RTs using only the best psms
        
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

function getPrecIDtoiRT(psms_dict::Dictionary{String, DataFrame}, RT_iRT::Any; 
                        max_q_value::AbstractFloat = 0.1,
                        max_precursors::Int = 250000)

    psms = vcat(values(psms_dict)...); #Combine data from all files in the experiment
    #Only consider precursors with a q_value passing the threshold
    # at least once accross the entire experiment
    filter!(x->x.q_value<=max_q_value,psms); 

    #For each precursor, what was the best scan, maximum q-value, and observed iRT for the best scan
    psms[!,:best_scan] .= false #Best scan for a given precursor_idx accross the experiment
    psms[!,:iRT_observed] .= zero(Float64) #Observed iRT for each psm
    #psms[!,:iRT_observed_best] .= zero(Float64) #Best iRT accross the experiment for a given precursor
    grouped_psms = groupby(psms, :precursor_idx); #Groupby precursor id
    prec_to_count = Dictionary{UInt32, Float32}()
    for (prec_idx, psms) in pairs(grouped_psms)
        insert!(prec_to_count, 
                prec_idx.precursor_idx,
                maximum(psms.prob))
    end
    sort!(prec_to_count, rev = true)
    best_precs = collect(keys(prec_to_count))[1:min(max_precursors, length(keys(prec_to_count)))]
    for prec_idx in ProgressBar(best_precs) #For each precursor 
        prec_psms = grouped_psms[(precursor_idx = prec_idx,)]
        #Best scan for the precursor accross the experiment
        best_scan = argmax(prec_psms.prob) 
        prec_psms.best_scan[best_scan] = true #Makr columsn containing the best psm for a precursor
        #Could alternatively estimate the iRT using the best-n scans 
        for j in range(1, size(prec_psms)[1])
            #Get empirical iRT for the psm given the empirical RT
            irt = RT_iRT[prec_psms.file_path[j]](prec_psms.RT[j])
            prec_psms.iRT_observed[j] = irt
        end
        prec_psms.iRT_observed[best_scan] = median(prec_psms.iRT_observed)
    end

    filter!(x->x.best_scan, psms); #Filter out non-best scans 
    #Map the precursor idx to the best empirical iRT extimate and the precursor m/z
    return Dictionary(psms[!,:precursor_idx], zip(psms[!,:iRT_observed], psms[!,:prec_mz]))
end

function getBestTrace!(psms::DataFrame)
    grouped_psms = groupby(psms, :precursor_idx)

    #For each precursor
    for i in ProgressBar(range(1, length(grouped_psms)))

        #Which traces had psms?
        iso_ranks = unique(grouped_psms[i][!,:iso_rank])

        length(iso_ranks)<=1 ? continue : nothing

        score, best_rank, best_score = 0, minimum(iso_ranks), 0
        
        for rank in iso_ranks
            for j in range(1, size(grouped_psms[i], 1))
                if grouped_psms[i][j,:iso_rank] == rank
                    if grouped_psms[i][j,:q_value] <=0.01
                        score += grouped_psms[i][j,:weight]
                    end
                end
            end
            if score == best_score
                if rank < best_rank
                    best_score = score
                    best_rank = rank
                end
                score = 0
            elseif score > best_score
                best_score = score
                best_rank = rank
                score = 0
            else
                score = 0
            end
        end
        
        for j in range(1, size(grouped_psms[i], 1))
            if grouped_psms[i][j,:iso_rank] != best_rank
                grouped_psms[i][j,:best_scan] = false
            end
        end

    end
    filter!(x->x.best_scan, psms);
end

function getCVFolds(precID_to_iRT::Dictionary{UInt32, Tuple{Float64, Float32}})
    precID_to_cv_fold = Dictionary{UInt32, UInt8}()
    for (prec_id, irt) in pairs(precID_to_iRT)
        insert!(precID_to_cv_fold,
        prec_id,
        rand(UInt8[0, 1]))
    end
    return precID_to_cv_fold
end
#=
N = 100000

good_precs = unique(PSMS[PSMS[!,:prob].>0.75,:precursor_idx])

N = 10000
test = MS2_CHROMS[(precursor_idx = good_precs[N], iso_rank = 1)][!,[:precursor_idx,:topn,:best_rank,:y_count,:b_count,:target,:scribe,
:city_block_fitted,:entropy_score,:matched_ratio,:RT,:weight,:prob]]
size(PSMS[PSMS[!,:precursor_idx].== good_precs[N],:])
plot(test[!,:RT],
    test[!,:weight], seriestype=:scatter)
N += 1
=#