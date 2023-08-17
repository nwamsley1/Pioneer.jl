mutable struct Counter{I,C<:Unsigned,T<:AbstractFloat}
    ids::Vector{I}
    counts::Vector{Tuple{C, T}}
    size::Int64
    matches::Int64
    function Counter(I::DataType, C::DataType, T::DataType, size::Int) #where {I,C<:Unsigned}
        new{I, C, T}(zeros(I, size), [(zero(C), zero(T)) for x in 1:size], 1, 0)
    end
end

getCount(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = first(c.counts[id])
getSums(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = last(c.sums[id])
getSize(c::Counter{I,C,T}) where {I,C<:Unsigned,T<:AbstractFloat} = c.size
incSize!(c::Counter{I,C,T}) where {I,C<:Unsigned,T<:AbstractFloat} = c.size += 1
getID(c::Counter{I,C,T}, idx::Int) where {I,C<:Unsigned,T<:AbstractFloat} = c.ids[idx]

function update!(c::Counter{I,C,T}, id::I, pred_intensity::T) where {C,I<:Unsigned,T<:AbstractFloat}
    c.counts[id] = (first(c.counts[id]) + one(C), last(c.counts[id]) + pred_intensity);
    return nothing
end

function reset!(c::Counter{I,C,T}, id::I) where {C,I<:Unsigned,T<:AbstractFloat}
    c.counts[id] = (zero(C), zero(T));
    return nothing
end

import DataStructures.inc!
function inc!(c::Counter{I,C,T}, id::I, pred_intensity::T) where {I,C<:Unsigned,T<:AbstractFloat} 
    if iszero(first(c.counts[id]))#c.counts[id]<1#iszero(c.counts[id])# == zero(C)
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
                by = id -> last(counter.counts[id]),
                #by = id -> first(counter.counts[id]),
                rev = true,
                alg=PartialQuickSort(1:topN)
        )
    return nothing
end

function reset!(c::Counter{I,C,T}) where {I,C<:Unsigned,T<:AbstractFloat} 
    #for i in 1:(getSize(c) - 1)
    for i in 1:(getSize(c) - 1)
    #for i in 1:length(c.ids)
    #for i in 1:(getSize(c))
        #if iszero(c.ids[i])
        #    continue
        #end
        c.counts[c.ids[i]] = (zero(UInt16), zero(Float32));
        c.ids[i] = zero(UInt32)
    end
    c.size, c.matches = 1, 0
    return nothing
end

function countFragMatches(c::Counter{I,C,T}, min_count::Int, min_ratio::T) where {I,C<:Unsigned,T<:AbstractFloat} 
    matched_frags = 0
    for i in 1:(getSize(c) - 1)
        id = c.ids[i]
        frag_count = getCount(c, id)
        matched_frags += frag_count
        if frag_count >= min_count
            if last(c.counts[id])>=min_ratio
                    c.ids[c.matches + 1] = c.ids[i]
                    
                    c.matches += 1
            end
        end
        c.counts[id] = (zero(UInt16), zero(Float32));
    end
    return matched_frags
end