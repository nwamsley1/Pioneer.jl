mutable struct Counter{I<:Unsigned,T<:AbstractFloat}
    ids::Vector{I}
    counts::Vector{T}
    size::Int64
    matches::Int64
    function Counter(I::DataType, C::DataType, T::DataType, size::Int) #where {I,C<:Unsigned}
        new{I, T}(zeros(I, size), zeros(T, size), 1, 0)
    end
end

getCount(c::Counter{I,T}, id::I) where {I<:Unsigned,T<:AbstractFloat} = c.counts[id]
getSize(c::Counter{I,T}) where {I<:Unsigned,T<:AbstractFloat} = c.size
incSize!(c::Counter{I,T}) where {I<:Unsigned,T<:AbstractFloat} = c.size += 1
getID(c::Counter{I,T}, idx::Int) where {I<:Unsigned,T<:AbstractFloat} = c.ids[idx]

function update!(c::Counter{I,C,T}, id::I, pred_intensity::T) where {C,I<:Unsigned,T<:AbstractFloat}
    @inbounds @fastmath c.counts[id] += pred_intensity;
    return nothing
end

function reset!(c::Counter{I,C,T}, id::I) where {C,I<:Unsigned,T<:AbstractFloat}
    c.counts[id] = (zero(C), zero(T));
    return nothing
end

function inc!(c::Counter{I,C,T}, id::I, pred_intensity::T) where {I,C<:Unsigned,T<:AbstractFloat} 
    if iszero(c.counts[id])#c.counts[id]<1#iszero(c.counts[id])# == zero(C)
        c.ids[getSize(c)] = id;
        update!(c,id,pred_intensity);
        incSize!(c);
    else
        update!(c,id,pred_intensity);
    end
    return nothing
end


import Base.sort!
function sort!(counter::Counter{I,C,T}, topN::Int) where {I,C<:Unsigned,T<:AbstractFloat}
    sort!(
                @view(counter.ids[1:counter.matches]), 
                by = id -> counter.counts[id],
                rev = true,
                alg=PartialQuickSort(topN)
        )
    return nothing
end

function reset!(c::Counter{I,T}) where {I<:Unsigned,T<:AbstractFloat} 
    #for i in 1:(getSize(c) - 1)
    @turbo for i in 1:(getSize(c) - 1)
        c.counts[c.ids[i]] = zero(T);
        c.ids[i] = zero(I)
    end
    c.size, c.matches = 1, 0
    return nothing
end
