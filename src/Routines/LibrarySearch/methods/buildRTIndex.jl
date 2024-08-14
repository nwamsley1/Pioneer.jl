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

#=
function getRTErr(Dictionary{UInt32, 
                            @NamedTuple{
                            best_prob::Float32, 
                            best_irt::Float32, 
                            min_irt::Union{Missing, Float32}, 
                            max_irt::Union{Missing, Float32}, 
                            mz::Float32}})
        
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
=#
#=
precs_passing = collect(keys(precID_to_iRT))
precs_irtest = Dictionary(
    precs_passing,
    [(mean_irt = zero(Float32), std_irt = zero(Float32), n = zero(UInt16)) for _ in range(1, length(precs_passing))]
    
)

for (key, psms_path) in pairs(psms_paths) #For each data frame 
    psms = Arrow.Table(psms_path)
    for i in ProgressBar(eachindex(psms[:precursor_idx]))
        precursor_idx = psms[:precursor_idx][i]
        if psms[:q_value][i]>0.01
            continue
        end
        if haskey(precs_irtest, precursor_idx)
            mean_irt, std_irt, n = precs_irtest[precursor_idx]
            n += one(UInt16)
            mean_irt += RT_iRT[key](psms[:RT][i])
            precs_irtest[precursor_idx]=(mean_irt = mean_irt, std_irt = std_irt, n = n)
        else
            insert!(
                precs_irtest,
                precursor_idx,
                (mean_irt = RT_iRT[key](psms[:RT][i]), std_irt = zero(Float32), n = one(UInt16))
            )
        end
    end
end
for (key, val) in pairs(precs_irtest)
    mean_irt, std_irt, n  = precs_irtest[key]
    mean_irt = mean_irt/n
    precs_irtest[key]= (mean_irt = mean_irt, std_irt = std_irt, n = n)
end
for (key, psms_path) in pairs(psms_paths) #For each data frame 
    psms = Arrow.Table(psms_path)
    for i in ProgressBar(eachindex(psms[:precursor_idx]))
        precursor_idx = psms[:precursor_idx][i]
        if haskey(precs_irtest, precursor_idx)
            mean_irt, std_irt, n = precs_irtest[precursor_idx]
            std_irt += (RT_iRT[key](psms[:RT][i]) - mean_irt)^2
            precs_irtest[precursor_idx]=(mean_irt = mean_irt, std_irt = std_irt, n = n)
        end
    end
end
for (key, val) in pairs(precs_irtest)
    mean_irt, std_irt, n  = precs_irtest[key]
    std_irt = sqrt(std_irt/(n-1))
    #std_irt = n/std_irt
    precs_irtest[key]= (mean_irt = mean_irt, std_irt = std_irt, n = n)
end
filter!(x->x[:n]>2, precs_irtest)
stds = [x[:std_irt] for x in values(precs_irtest)]

    zerso(@NamedTuple{
        best_prob::Float32, 
        best_irt::Float32, 
        min_irt::Union{Missing, Float32}, 
        =#