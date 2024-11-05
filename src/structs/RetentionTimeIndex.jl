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
getrtBins(rti::retentionTimeIndex) = rti.rt_bins
getrtBin(rti::retentionTimeIndex, rt_bin::Int) = rti.rt_bins[rt_bin]
function retentionTimeIndex(T::DataType, U::DataType) 
    return retentionTimeIndex(Vector{rtIndexBin{T, U}}())
end
