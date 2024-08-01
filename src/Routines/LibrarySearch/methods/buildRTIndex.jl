struct rtIndexBin{T,U<:AbstractFloat}
    lb::T
    ub::T
    prec::Vector{Tuple{UInt32, U}}
end
getLow(r::rtIndexBin) = r.lb
getHigh(r::rtIndexBin) = r.ub
function compare_lb(rb::rtIndexBin{T,U}) where {T,U<:AbstractFloat}
    return rb.lb
end
getLB(rb::rtIndexBin{T,U}) where {T,U<:AbstractFloat} = rb.lb
getMZ(rb::rtIndexBin{T, U}) where {T,U<:AbstractFloat} = last(rb.prec)
getPrecID(rb::rtIndexBin{T, U}) where {T,U<:AbstractFloat} = first(rb.prec)

struct retentionTimeIndex{T,U<:AbstractFloat}
    rt_bins::Vector{rtIndexBin{T, U}}
end
getRTBins(rti::retentionTimeIndex) = rti.rt_bins
getRTBin(rti::retentionTimeIndex, rt_bin::Int) = rti.rt_bins[rt_bin]
function retentionTimeIndex(T::DataType, U::DataType) 
    return retentionTimeIndex(Vector{rtIndexBin{T, U}}())
end

function mapRTandiRT(psms_dict::Dictionary{String, DataFrame}; 
                        min_prob::AbstractFloat = 0.9
                    )

    #Dictionaries mapping fild_id names to data
    #prec_ids = Dictionary{String, Set{UInt32}}() #File name => precursor id
    iRT_RT = Dictionary{String, Any}() #File name => KDEmapping from iRT to RT 
    RT_iRT = Dictionary{String, Any}() #File name => KDEmapping from iRT to RT 
    for (key, psms) in pairs(psms_dict) #For each data frame 
        psms[!,:file_path] .= key #Add column for file name

        #insert!(prec_ids,key,Set(psms[!,:precursor_idx])) #Set of precursors ids in the dataframe

        best_hits = psms[!,:prob].>min_prob#Map RTs using only the best psms
        irt_to_rt_spline = UniformSpline(
                                    psms[best_hits,:RT],
                                    psms[best_hits,:iRT_predicted],
                                    
                                    3, 
                                    5
        )
        rt_to_irt_spline = UniformSpline(
            psms[best_hits,:iRT_predicted],
            psms[best_hits,:RT],
            3, 
            5
        )
        #Build RT=>iRT and iRT=> RT mappings for the file and add to the dictionaries 
        insert!(iRT_RT, key, irt_to_rt_spline)
        insert!(RT_iRT, key, rt_to_irt_spline)

        #Update iRT error and iRT observed based on first search results
        psms[!,:iRT_observed] = Float16.(rt_to_irt_spline.(psms[!,:RT]))
        psms[!,:iRT_error] = Float16.(abs.(psms[!,:iRT_observed] .- psms[!,:iRT_predicted]))

    end
    return iRT_RT, RT_iRT
end

function getPrecIDtoiRT(psms_dict::Dictionary{String, DataFrame}, 
                        RT_iRT::Any; 
                        max_precursors::Int = 250000)

    psms = vcat(values(psms_dict)...); #Combine data from all files in the experiment
    #Only consider precursors with a q_value passing the threshold
    # at least once accross the entire experiment

    #filter!(x->x.q_value<=max_q_value,psms); TRYING THIS 03/04/24

    #For each precursor, what was the best scan, maximum q-value, and observed iRT for the best scan
    psms[!,:best_scan] .= false #Best scan for a given precursor_idx accross the experiment
    psms[!,:iRT_observed] .= zero(Float64) #Observed iRT for each psm
    grouped_psms = groupby(psms, :precursor_idx); #Groupby precursor id
    prec_to_count = Dictionary{UInt32, Float32}()
    for (prec_idx, psms) in pairs(grouped_psms)
        insert!(prec_to_count, 
                prec_idx.precursor_idx,
                maximum(psms.prob))
    end
    sort!(prec_to_count, rev = true)
    best_precs = collect(keys(prec_to_count))[1:min(max_precursors, length(keys(prec_to_count)))]
    for prec_idx in best_precs #For each precursor 
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

buildRTIndex(PSMs::DataFrame; bin_rt_size::AbstractFloat = 0.1) = buildRTIndex(PSMs[:,:RT], PSMs[:,:prec_mz], PSMs[:,:precursor_idx], bin_rt_size)

buildRTIndex(PSMs::SubDataFrame; bin_rt_size::AbstractFloat = 0.1) = buildRTIndex(PSMs[:,:RT], PSMs[:,:prec_mz], PSMs[:,:precursor_idx], bin_rt_size)


function makeRTIndices(psms_dict::Dictionary{String, DataFrame}, 
                       precID_to_iRT::Dictionary{UInt32, Tuple{Float64, Float32}},
                       RT_to_iRT::Any;
                       bin_rt_size = 0.1,
                       min_prob::AbstractFloat = 0.5)
    #Maps filepath to a retentionTimeIndex (see buildRTIndex.jl)
    rt_indices = Dictionary{String, retentionTimeIndex{Float32, Float32}}()
    #iRT dict
    #iRT_dict = Dictionary{String, UnorderedDictionary{UInt32, Float32}}()
    for (file_path, psms) in pairs(psms_dict) #For each file in the experiment
        #Impute empirical iRT value for psms with probability lower than the threshold
        #filter!(x->x.prob>=min_prob, psms); 
        iRTs, mzs, prec_ids = zeros(Float32, length(precID_to_iRT)), zeros(Float32, length(precID_to_iRT)), zeros(UInt32, length(precID_to_iRT))

        prec_set = Dictionary(psms[!,:precursor_idx]::Vector{UInt32},
                            zip(psms[!,:RT]::Vector{Float32}, 
                                psms[!,:prec_mz]::Vector{Float32},
                                psms[!,:prob]::Vector{Float32})
                            )
        #Set of precursors where empirical iRTs can be used and do not need to be imputed 
        #pset = Set(psms[!,:precursor_idx]::Vector{UInt32})
        #Loop through all precursors in the seed 
        Threads.@threads for (i, (prec_id, irt_mz_imputed)) in collect(enumerate(pairs(precID_to_iRT)))
            prec_ids[i] = prec_id
            if haskey(prec_set, prec_id) #prec_id ∈ pset::Set{UInt32} 
                if last(prec_set[prec_id]) >= min_prob #Use empirical iRT/RT
                    rt, mz, _ = prec_set[prec_id]::Tuple{Float32, Float32, Float32}
                    iRTs[i] = RT_to_iRT[file_path](rt)
                    mzs[i] = mz
                    continue
                end
            end
            #Impute empirical iRT from the best observed psm for the precursor accross the experiment 
            irt, mz = irt_mz_imputed::Tuple{Float64, Float32}
            #rt = iRT_RT[file_path](irt)::Float64 #Convert iRT to RT 
            iRTs[i] = irt
            mzs[i] = mz
        end
        
        #Build RT index 
        rt_df = DataFrame(Dict(:RT => iRTs,
                                :prec_mz => mzs,
                                :precursor_idx => prec_ids))
        sort!(rt_df, :RT)
        rt_index = buildRTIndex(rt_df, bin_rt_size=bin_rt_size);
        insert!(rt_indices, file_path, rt_index)
    end
    return rt_indices
end

function getCVFolds(
                prec_ids::AbstractVector{UInt32},
                protein_groups::AbstractVector{String}
                )
    uniq_pgs = unique(protein_groups) #Unique protein groups
    pg_to_cv_fold = Dictionary{String, UInt8}() #Protein group to cross-validation fold
    #Map protein groups to cv folds 
    for pg in uniq_pgs
        insert!(pg_to_cv_fold,
        pg,
        rand(UInt8[0, 1])
        )
    end
    #Now map pid's to cv folds 
    pid_to_cv_fold = Dictionary{UInt32, UInt8}()
    for pid in prec_ids
        insert!(
        pid_to_cv_fold,
        pid,
        pg_to_cv_fold[protein_groups[pid]]
        )
    end
    return pid_to_cv_fold
end

function getRTErr(psms_dict::Dictionary{String, DataFrame},
                  RT_iRT::Dictionary{String, Any})
        
    for (key, psms) in pairs(psms_dict)
        psms[!,:iRT_observed] = RT_iRT[key].(psms[!,:RT])
    end

    combined_psms = vcat(values(psms_dict)...);
    filter!(x->x.q_value<=0.01, combined_psms)
    psms_by_prec = groupby(combined_psms, :precursor_idx)
    irt_mads = Vector{Union{Missing, Float32}}(undef, length(psms_by_prec))
    j = 1
    for (prec, psms) in pairs(psms_by_prec)
        if size(psms, 1) > 1
            irt_mads[j] = maximum(psms[!,:iRT_observed]) -  minimum(psms[!,:iRT_observed])
        else
            irt_mads[j] = missing
        end
        j += 1
    end
    return quantile(skipmissing(irt_mads), 0.99)#median(skipmissing(rt_mads))*4
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