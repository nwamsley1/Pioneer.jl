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
function getMatchedRatio(c::Counter{I,C,T}, totals::Vector{Float32}, id::I) where {C,I<:Unsigned,T<:AbstractFloat}
    #matched_intensity = getSums(c, id)
    return first(c.counts[id])/(totals[id] - first(c.counts[id]))
end

getSize(c::Counter{I,C,T}) where {I,C<:Unsigned,T<:AbstractFloat} = c.size
incSize!(c::Counter{I,C,T}) where {I,C<:Unsigned,T<:AbstractFloat} = c.size += 1
getID(c::Counter{I,C,T}, idx::Int) where {I,C<:Unsigned,T<:AbstractFloat} = c.ids[idx]

function update!(c::Counter{I,C,T}, id::I, pred_intensity::T) where {C,I<:Unsigned,T<:AbstractFloat}
    #count = c.counts[id]
    #c.counts[id] = (first(count) + one(C), last(count) + pred_intensity);
    c.counts[id] = (first(c.counts[id]) + one(C), last(c.counts[id]) + pred_intensity);
    return nothing
end

function reset!(c::Counter{I,C,T}, id::I) where {C,I<:Unsigned,T<:AbstractFloat}
    c.counts[id] = (zero(C), zero(T));
    return nothing
end
#getObs(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getObs(c.dotp[id])
#getPred(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getPred(c.dotp[id])
#getObsPred(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getObsPred(c.dotp[id])
#getCount(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getCount(c.dotp[id])
#getDP(c::Counter{I,C,T}, id::I) where {I,C<:Unsigned,T<:AbstractFloat} = getDP(c.[id])

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
function sort!(counter::Counter{I,C,T}, prosit_totals::Vector{Float32}, topN::Int) where {I,C<:Unsigned,T<:AbstractFloat}
    sort!(
                @view(counter.ids[1:counter.matches]), 
                by = id -> getMatchedRatio(counter, prosit_totals, id),
                rev = true,
                alg=PartialQuickSort(1:topN)
             )
    return nothing
end

function reset!(c::Counter{I,C,T}) where {I,C<:Unsigned,T<:AbstractFloat} 
    #@inbounds @simd for i in 1:(getSize(c) - 1)
    for i in 1:(getSize(c) - 1)
        c.counts[c.ids[i]] = (zero(UInt8), zero(Float32));
    end
    c.size = 1
    c.matches = 0
    return nothing
end

function countFragMatches(c::Counter{I,C,T}, prosit_totals::Vector{Float32}, min_count::Int) where {I,C<:Unsigned,T<:AbstractFloat} 
    frag_counts = 0
    for i in 1:(getSize(c) - 1)
        id = c.ids[i]
        frag_count = getCount(c, id)
        frag_counts += frag_count
        if frag_count >= min_count
            if getMatchedRatio(c, prosit_totals, id)>=-1000.0#0.8
                    c.ids[c.matches + 1] = c.ids[i]
                    c.matches += 1
            end
        end
    end
    return frag_counts
end

precs = Counter(UInt32, UInt8, Float32, 9387261) #Prec counter

prec_ids = [UInt32(x) for x in rand(1:1000, 1000000)]
prec_intensity = [Float32(x) for x in randn(Float32, 1000000)]
frag_matches = [(id, int) for (id, int) in zip(prec_ids, prec_intensity)]
prosit_totals = [randn(Float32) for x in 1:9387261]
function benchmark(precs::Counter{I,C,T}, frag_matches::Vector{Tuple{UInt32, Float32}}, prosit_totals::Vector{Float32}) where {I,C<:Unsigned,T<:AbstractFloat} 
    for match in frag_matches
        inc!(precs, first(match), last(match));
    end
    counts = countFragMatches(precs, prosit_totals, 1);
    sort!(precs,  prosit_totals, 20);
    reset!(precs);
    return counts
    #return 0
end

function benchmark(precs::Counter{I,C,T}, frag_matches::Vector{Tuple{UInt32, Float32}}, prosit_totals::Vector{Float32}, rep::Int) where {I,C<:Unsigned,T<:AbstractFloat} 
    for i in 1:rep
        benchmark(precs, frag_matches, prosit_totals)
    end
    return nothing
end