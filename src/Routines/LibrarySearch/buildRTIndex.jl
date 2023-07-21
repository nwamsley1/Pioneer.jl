struct rtBin{T,U<:AbstractFloat}
    lb::T
    ub::T
    prec::Vector{Tuple{UInt32, U}}
end

function compare_lb(rb::rtBin{T,U}) where {T,U<:AbstractFloat}
    return rb.lb
end

getLB(rb::rtBin{T,U}) where {T,U<:AbstractFloat} = rb.lb
getMZ(rb::rtBin{T, U}) where {T,U<:AbstractFloat} = last(rb.prec)
getPrecID(rb::rtBin{T, U}) where {T,U<:AbstractFloat} = first(rb.prec)

struct retentionTimeIndex{T,U<:AbstractFloat}
    rt_bins::Vector{rtBin{T, U}}
end

function retentionTimeIndex(T::DataType, U::DataType) 
    return retentionTimeIndex(Vector{rtBin{T, U}}())
end

function buildRTIndex(RTs::Vector{T}, prec_mzs::Vector{U}, prec_ids::Vector{I}, bin_rt_size::AbstractFloat) where {T,U<:AbstractFloat,I<:Integer}
    
    start_idx = 1
    start_RT =  RTs[start_idx]
    rt_index = retentionTimeIndex(T, U) #Initialize retention time index

    for (i, RT) in enumerate(RTs)
        if (RT - start_RT) > bin_rt_size
            push!(rt_index.rt_bins, 
                    rtBin(RTs[start_idx], #Retention time for first precursor in the bin
                          RTs[i - 1],     #Retention time for last precursor in the bin
                        [(zero(UInt32), zero(Float32)) for _ in 1:(i - start_idx)] #Pre-allocate precursors 
                        )
                )

            n = 1 #n'th precursor 
            for idx in start_idx:(i - 1)
                rt_index.rt_bins[end].prec[n] = (prec_ids[idx], prec_mzs[idx]) #Add n'th precursor
                n += 1
            end

            sort!(rt_index.rt_bins[end].prec, by = x->last(x)) #Sort precursors by m/z

            start_idx = i
            start_RT = RTs[start_idx]
            continue
        end
    end

    return rt_index
end

buildRTIndex(PSMs::DataFrame, bin_rt_size::AbstractFloat = 0.1) = buildRTIndex(PSMs[:,:RT], PSMs[:,:prec_mz], PSMs[:,:precursor_idx], bin_rt_size)
