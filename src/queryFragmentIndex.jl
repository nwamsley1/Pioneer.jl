function QueryFragmentIndex(frag_index::Vector{FragBin{T}}, query::T) where {T<:AbstractFloat}
    #Binary Search
    lo, hi = 1, length(frag_index)
    while lo <= hi
        mid = (lo + hi) ÷ 2
        if getUB(frag_index[mid]) < (query)
             lo = mid + 1
        elseif getLB(frag_index[mid]) > (query)
            hi = mid - 1
        else
            return lo
        end
    end
end