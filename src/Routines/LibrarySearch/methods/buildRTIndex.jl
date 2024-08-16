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

buildRTIndex(PSMs::DataFrame; bin_rt_size::AbstractFloat = 0.1) = buildRTIndex(PSMs[:,:irt], PSMs[:,:prec_mz], PSMs[:,:precursor_idx], bin_rt_size)

buildRTIndex(PSMs::SubDataFrame; bin_rt_size::AbstractFloat = 0.1) = buildRTIndex(PSMs[:,:irt], PSMs[:,:prec_mz], PSMs[:,:precursor_idx], bin_rt_size)


function makeRTIndices(temp_folder::String,

                       psms_paths::Dictionary{String, String}, 
                       prec_to_irt::Dictionary{UInt32, @NamedTuple{irt::Float32, mz::Float32}},
                       rt_to_irt_splines::Any;
                       min_prob::AbstractFloat = 0.5)

    #Maps filepath to a retentionTimeIndex (see buildRTIndex.jl)
    rt_index_paths = Dictionary{String, String}()
    #Fill retention time index for each file. 
    for (key, psms_path) in pairs(psms_paths)
        psms = Arrow.Table(psms_path)
        rt_to_irt = rt_to_irt_splines[key]
        #Impute empirical iRT value for psms with probability lower than the threshold
        irts = zeros(Float32, length(prec_to_irt))
        mzs = zeros(Float32, length(prec_to_irt))
        prec_ids = zeros(UInt32, length(prec_to_irt))
        #map observec precursors to irt and probability score
        prec_set = Dict(zip(
            psms[:precursor_idx],
            map(x->(irt=first(x),prob=last(x)), zip(rt_to_irt.(psms[:RT]), psms[:prob]))
        ))

        Threads.@threads for (i, (prec_id, irt_mz)) in collect(enumerate(pairs(prec_to_irt)))
            prec_ids[i] = prec_id
            irt, mz = irt_mz::@NamedTuple{irt::Float32, mz::Float32}
            #Don't impute irt, use empirical
            if haskey(prec_set, prec_id)
                irt, prob = prec_set[prec_id]
                if prob >= min_prob #Use empirical iRT/RT
                    irts[i] = irt
                    mzs[i] = mz
                    continue
                end
            end
            #Impute irt from the best observed psm for the precursor accross the experiment 
            irts[i] = irt
            mzs[i] = mz
        end

        #Build RT index 
        rt_df = DataFrame(Dict(:irt => irts,
                                :prec_mz => mzs,
                                :precursor_idx => prec_ids))
        sort!(rt_df, :irt)
        temp_path =joinpath(temp_folder, key*"rt_indices.arrow")
        Arrow.write(
            temp_path,
            rt_df,
            )
        insert!(
            rt_index_paths,
            key,
            temp_path
        )
        #rt_index = buildRTIndex(rt_df, bin_rt_size=bin_rt_size);
        #insert!(rt_indices, key, rt_index)
    end
    return rt_index_paths
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