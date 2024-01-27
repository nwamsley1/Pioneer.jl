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

mutable struct Counter{I,C<:Unsigned}
    ids::Vector{I}
    counts::Vector{C}
    size::Int64
    matches::Int64
    function Counter(I::DataType, C::DataType, size::Int) #where {I,C<:Unsigned}
        new{I, C}(zeros(I, size), zeros(C, size), 1, 0)
    end
end

getCount(c::Counter{I,C}, id::I) where {I,C<:Unsigned} = c.counts[id]
getSize(c::Counter{I,C}) where {I,C<:Unsigned} = c.size
incSize!(c::Counter{I,C}) where {I,C<:Unsigned} = c.size += 1
getID(c::Counter{I,C}, idx::Int) where {I,C<:Unsigned} = c.ids[idx]

function update!(c::Counter{I,C}, id::I, pred_intensity::C) where {I,C<:Unsigned}
    @inbounds @fastmath c.counts[id] += pred_intensity;
    return nothing
end

function reset!(c::Counter{I,C}, id::I) where {I,C<:Unsigned}
    c.counts[id] = zero(T);
    return nothing
end

function inc!(c::Counter{I,C}, id::I, pred_intensity::C) where {I,C<:Unsigned} 
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
function sort!(counter::Counter{I,C}, topN::Int) where {I,C<:Unsigned}
    sort!(
                @view(counter.ids[1:counter.matches]), 
                by = id -> counter.counts[id],
                rev = true,
                alg=PartialQuickSort(topN)
        )
    return nothing
end

function reset!(c::Counter{I,C}) where {I,C<:Unsigned}
    #for i in 1:(getSize(c) - 1)
    @turbo for i in 1:(getSize(c) - 1)
        c.counts[c.ids[i]] = zero(T);
        c.ids[i] = zero(I)
    end
    c.size, c.matches = 1, 0
    return nothing
end

function countFragMatches(c::Counter{I,C}, min_count::C) where {I,C<:Unsigned}
    matched_frags = 0
    @inbounds for i in 1:(getSize(c) - 1)
        id = c.ids[i]
        if getCount(c, id)>min_count
                c.ids[c.matches + 1] = c.ids[i]
                c.matches += 1
        end
        c.counts[id] = zero(Float32);
    end
    return matched_frags
end
