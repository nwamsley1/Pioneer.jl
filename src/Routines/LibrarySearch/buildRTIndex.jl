struct rtBin{T,U<:AbstractFloat}
    lb::T
    ub::T
    prec::Vector{Tuple{UInt32, U}}
end

function compare_lb(rb::rtBin{T,U}) where {T,U<:AbstractFloat}
    return rb.lb
end

getLB(rb::rtBin{T,U}) where {T,U<:AbstractFloat} = rb.lb

struct retentionTimeIndex{T,U<:AbstractFloat}
    rt_bins::Vector{rtBin{T, U}}
end

function retentionTimeIndex(T::DataType, U::DataType) 
    return retentionTimeIndex(Vector{rtBin{T, U}}())
end

function buildRTIndex(RTs::Vector{T}, prec_mzs::Vector{U}, prec_ids::Vector{UInt32}, bin_rt_size::AbstractFloat) where {T,U<:AbstractFloat}
    start_idx = 1
    start_RT =  RTs[start_idx]
    bin = one(UInt32)
    rt_index = retentionTimeIndex(Float64, Float32)
    for (i, RT) in enumerate(RTs)
        if (RT - start_RT) > bin_rt_size
            push!(rt_index.rt_bins, 
                    rtBin(RT, 
                                     RTs[i - 1], 
                                     [(zero(UInt32), zero(Float32)) for _ in 1:(i - start_idx)])
                )
            n = 1
            for idx in start_idx:(i - 1)
                rt_index.rt_bins[end].prec[n] = (prec_ids[idx], prec_mzs[idx])
                n += 1
            end
            sort!(rt_index.rt_bins[end].prec, by = x->last(x))
            start_idx = i
            start_RT = RTs[start_idx]
            bin += one(UInt32)
            continue
        end
    end
    return rt_index
end

buildRTIndex(PSMs::DataFrame, bin_rt_size::AbstractFloat = 0.1) = buildRTIndex(PSMs[:,:RT], PSMs[:,:prec_mz], PSMs[:,:precursor_idx], bin_rt_size)
