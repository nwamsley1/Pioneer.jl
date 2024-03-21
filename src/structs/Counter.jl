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
    @inbounds c.counts[id] += pred_intensity;
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

function inc!(c::Counter{I,C}, id::I, pred_intensity::C) where {I,C<:Unsigned} 
    @inbounds @fastmath begin 
        no_previous_encounter = iszero(c.counts[id])
        c.ids[c.size] += id*no_previous_encounter;
        c.size += (no_previous_encounter)
        c.counts[id] += pred_intensity;
    end
    return nothing
end

function inc!(c::Counter{I,C}, id::I, pred_intensity::C) where {I,C<:Unsigned} 
    @inbounds @fastmath begin 
        no_previous_encounter = iszero(c.counts[id])
        c.ids[c.size] += no_previous_encounter ? id : zero(I)
        c.size += no_previous_encounter ? one(Int64) : zero(Int64)
        c.counts[id] += pred_intensity;
    end
    return nothing
end

import Base.sort!
function sortCounter!(counter::Counter{I,C}) where {I,C<:Unsigned}
    sort!(
                @view(counter.ids[1:counter.matches]), 
                by = id -> counter.counts[id],
                rev = true,
                alg=PartialQuickSort(counter.matches)
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
    @inbounds for i in 1:(getSize(c) - 1)
        id = c.ids[i]
        if getCount(c, id)>=min_count
                c.ids[c.matches + 1] = c.ids[i]
                c.matches += 1
        end
        c.counts[id] = zero(Float32);
    end
    return 0#c.matches
end
