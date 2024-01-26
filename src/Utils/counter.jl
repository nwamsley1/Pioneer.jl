mutable struct ArrayDict{I<:Unsigned,C<:Real}
    keys::Vector{I}
    vals::Vector{C}
    size::Int64
    function ArrayDict(I::DataType, C::DataType, size::Int) #where {I,C<:Unsigned}
        new{I, C}(zeros(I, size), zeros(C, size), 0)
    end
end

function update!(c::ArrayDict{I,C}, key::I, val::C) where {I,C<:Unsigned}
    c.size += 1
    c.vals[key] = val
    c.keys[c.size] = key
end

function reset!(c::ArrayDict{I,C}) where {I,C<:Unsigned}
    @turbo for i in range(1, c.size)
        c.vals[c.keys[i]] = zero(I)
        c.keys[i] = zero(C)
    end
    c.size = 0
end

function Base.getindex(c::ArrayDict{I,C}, i::Ti) where {I,C<:Unsigned, Ti<:Integer}
    c.vals[i]
end

mutable struct Counter{I<:Unsigned,T<:AbstractFloat}
    ids::Vector{I}
    counts::Vector{T}
    size::Int64
    matches::Int64
    function Counter(I::DataType, T::DataType, size::Int) #where {I,C<:Unsigned}
        new{I, T}(zeros(I, size), zeros(T, size), 1, 0)
    end
end

getCount(c::Counter{I,T}, id::I) where {I<:Unsigned,T<:AbstractFloat} = c.counts[id]
getSize(c::Counter{I,T}) where {I<:Unsigned,T<:AbstractFloat} = c.size
incSize!(c::Counter{I,T}) where {I<:Unsigned,T<:AbstractFloat} = c.size += 1
getID(c::Counter{I,T}, idx::Int) where {I<:Unsigned,T<:AbstractFloat} = c.ids[idx]

function update!(c::Counter{I,T}, id::I, pred_intensity::T) where {I<:Unsigned,T<:AbstractFloat}
    @inbounds @fastmath c.counts[id] += pred_intensity;
    return nothing
end

function reset!(c::Counter{I,T}, id::I) where {I<:Unsigned,T<:AbstractFloat}
    c.counts[id] = zero(T);
    return nothing
end

function inc!(c::Counter{I,T}, id::I, pred_intensity::T) where {I<:Unsigned,T<:AbstractFloat} 
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
function sort!(counter::Counter{I,T}, topN::Int) where {I<:Unsigned,T<:AbstractFloat}
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

function countFragMatches(c::Counter{I,T}, min_count::T) where {I<:Unsigned,T<:AbstractFloat} 
    matched_frags = 0
    @inbounds for i in 1:(getSize(c) - 1)
        id = c.ids[i]
        if getCount(c, id)>=min_count
                c.ids[c.matches + 1] = c.ids[i]
                c.matches += 1
        end
        c.counts[id] = zero(Float32);
    end
    return matched_frags
end

#=
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
    @fastmath c.counts[id] = (first(c.counts[id]) + one(C), last(c.counts[id]) + pred_intensity);
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
    @inbounds for i in 1:(getSize(c) - 1)
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
    @inbounds for i in 1:(getSize(c) - 1)
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
=#