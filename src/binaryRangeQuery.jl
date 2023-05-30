"""
    binaryGetNearest(arr::Vector{Union{Missing, U}}, query::T, low_tol::T, high_tol::T) where {U,T <:Real}

Finds and returns the element in arr nearest to the query value so long as it is within the interval [query - low_tol, query + high_tol].
Returns 0.0 if nothing in arr is in that interval. `arr` is assumed to be sorted in ascending order

- `arr::Vector{Union{Missing, U}}` -- Vector of real values. Assumed to be sorted in ascending order. 
- `query::T` -- Real valued number. Find the element in arr nearest this value
- `low_tol::T` Real valued number. Don't consider matches below query - low_tol
- `high_tol::T` Real valued number. Don't consider matches above query + high_tol

"""
function binaryGetNearest(arr::Vector{Union{Missing, U}}, query::T, low_tol::T, high_tol::T) where {U,T <:Real}

    #Check special cases (is the answer on the boundary or is the array empty?)
    n = length(arr)
    #Don't bother if there can't be any matches 
    if n == 0 return 0 end
    if query <= arr[1] - high_tol return 0 end
    if query >= arr[n] + low_tol return 0 end

    #Given a sorted vector `arr` and indices lo, and hi, use a linear search to
    #find the element in arr[lo:hi] with the smallesst absolute difference from `query`. 
    function getNearest(arr::Vector{Union{Missing, T}}, lo::Int, hi::Int, query::U) where {T,U <: Real}
        if hi - lo>1
            #Smallest distance observed so far
            smallest_distance = abs(query - arr[lo])
            best_idx = 1
            for (i, mass) in enumerate(@view(arr[lo:hi]))
                if abs(query - mass) < smallest_distance
                    smallest_distance = abs(query - mass)
                    best_idx = i
                end
            end
            return best_idx
        else
            return 1
        end
    end

    #Binary search. Start with indices lo = 1 and hi = n
    lo, hi = 1, n
    while lo <= hi
        mid = (lo + hi) ÷ 2
        if arr[mid] < (query - low_tol)
             lo = mid + 1
        elseif arr[mid] > (query + high_tol)
            hi = mid - 1
        else
            #Any index >=lo or <=hi is in the interval  [query - low_tol, query + high_tol].
            #Now find the element of arr[lo:hi] nearest to query. 
            return getNearest(arr,lo, hi, query) + lo - 1
        end
    end

    return 0

end
export binaryGetNearest

"""
    getPrecursors(window_center::T, precursorList::Vector{Precursor{T}}, params) where {T<:AbstractFloat}

Given a list of `Precursor{T}` sorted in ascending order by m/z, returns the subset where each precursor is
within the m/z tolerance of the window center. 

- `window_center::T` -- Center of the isolation window
- `precursorList::Vector{Precursor{T}}` -- Must be sorted in ascenging order by m/z. 
- `params` -- Must include params[:lower_tol] and params[:upper_tol] which are Real valued. 
"""
function getPrecursors(window_center::T, precursorList::Vector{Precursor{T}}, params) where {T<:AbstractFloat}
    l_bnd, u_bnd = window_center - params[:lower_tol], window_center + params[:upper_tol]
    start, stop = searchsortedfirst(precursorList, l_bnd,lt=(t,x)->getMZ(t)<x), searchsortedlast(precursorList, u_bnd,lt=(x,t)->getMZ(t)>x)
    return @view(precursorList[start:stop])
end
export getPrecursors