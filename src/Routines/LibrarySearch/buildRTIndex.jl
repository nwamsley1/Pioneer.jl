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
                        iRT_RT::Any)
    rt_indices = Dictionary{String, retentionTimeIndex{Float32, Float32}}()
    for (file_path, psms) in ProgressBar(pairs(psms_dict))
        RTs = zeros(Float32, length(precID_to_iRT))
        mzs = zeros(Float32, length(precID_to_iRT))
        prec_ids = zeros(UInt32, length(precID_to_iRT))
        prec_set = Dictionary(psms[!,:precursor_idx]::Vector{UInt32},
                            zip(psms[!,:RT]::Vector{Float32}, 
                                psms[!,:prec_mz]::Vector{Float32})
                            )
        pset = Set(psms[!,:precursor_idx]::Vector{UInt32})
        #i = 1
        Threads.@threads for (i, (prec_id, value)) in collect(enumerate(pairs(precID_to_iRT)))
            prec_ids[i] = prec_id
            if prec_id ∈ pset::Set{UInt32}
                rt, mz = prec_set[prec_id]::Tuple{Float32, Float32}
                RTs[i] = rt
                mzs[i] = mz
            else
                irt, mz = precID_to_iRT[prec_id]::Tuple{Float64, Float32}
                rt = iRT_RT[file_path](irt)::Float64
                RTs[i] = rt
                mzs[i] = mz
            end
            #i += 1
        end
        rt_df = DataFrame(Dict(:RT => RTs,
        :prec_mz => mzs,
        :precursor_idx => prec_ids))
        sort!(rt_df, :RT)
        rt_index = buildRTIndex(rt_df);
        insert!(rt_indices, file_path, rt_index)
    end
    return rt_indices
end


function mapRTandiRT(psms_dict::Dictionary{String, DataFrame})

    #Dictionaries mapping fild_id names to data

    prec_ids = Dictionary{String, Set{UInt32}}() #File name => precursor id
    iRT_RT = Dictionary{String, Any}() #File name => KDEmapping from iRT to RT 
    RT_iRT = Dictionary{String, Any}() #File name => KDEmapping from iRT to RT 
    for (key, psms) in ProgressBar(pairs(psms_dict)) #For each data frame 
        psms[!,:file_path] .= key #Add column for file name

        insert!(prec_ids,key,Set(psms[!,:precursor_idx])) #Set of precursors ids in the dataframe

        best_hits = (psms[!,:spectral_contrast].>0.9 #Map RTs using only the best psms 
                    ) .& (psms[!,:decoy].==false
                    ) .& (psms[!,:entropy_score].>0.9
                    ) .& (psms[!,:total_ions].>6)
        
        #Build RT=>iRT and iRT=> RT mappings for the file and add to the dictionaries 
        insert!(iRT_RT, key, KDEmapping(psms[best_hits,:iRT],
        psms[best_hits,:RT]
        ))
        insert!(RT_iRT, key, KDEmapping(psms[best_hits,:RT],
        psms[best_hits,:iRT]
        ))
    end
    return (:iRT_RT => iRT_RT, :RT_iRT => RT_iRT)
end

function getPrecIDtoiRT(psms_dict::Dictionary{String, DataFrame}, RT_iRT::Any)

    psms = vcat(values(psms_dict)...); #Combine data from all files

    #For each precursor, what was the best scan, maximum q-value, and observed iRT for the best scan
    psms[!,:best_scan] .= false #Best scan for a given precursor_idx
    psms[!,:max_q] .= Float64(Inf) #Maximum q-value for a given precursor 
    psms[!,:iRT_observed] .= zero(Float64) #Observed iRT for the best scan
    grouped_psms = groupby(psms, :precursor_idx); #Groupby precursor idx
    for i in range(1, length(grouped_psms)) #For each precursor 
        best_scan = argmin(grouped_psms[i].q_value) #Best scan for teh precursor 
        grouped_psms[i].best_scan[best_scan] = true
        grouped_psms[i].max_q[best_scan] = maximum(grouped_psms[i].q_value)
        irt = RT_iRT[grouped_PSMs[i].file_path[best_scan]](grouped_psms[i].RT[best_scan])
        grouped_psms[i].iRT_observed[best_scan] = irt
    end
    filter!(x->x.best_scan, psms)
    return Dictionary(psms[!,:precursor_idx], zip(psms[!,:iRT_observed], psms[!,:prec_mz]))
end