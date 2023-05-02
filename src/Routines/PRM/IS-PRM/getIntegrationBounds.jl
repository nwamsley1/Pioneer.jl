
"""
    getScanPairs(H::Vector{T}, L::Vector{T}; max_diff::T = T(3)) where T <: Number

Given sorted lists of retention times of MSN scans for the heavy and light precursor respectively, aligns
the two list. Returns two list of indices, `heavy_scans` and `light_scans`, for `H` and `L` respectively of equal length. 
heavy_scans[i] is paired with light_scans[i]. Each heavy scan is paired with its nearest light scan by retention time
given the difference bewteen the two is less than `max_diff` (in minutes). 

### Input

- `H::Vector{T}` -- Sorted (ascending) retention time for MSn scans targeting the heavy precursor
- `L::Vector{T}` -- Sorted (ascending) retention time for MSn scans targeting the light precursor
- `max_diff::T = T(3)` -- Maximum allowable difference in retention time (minutes) bewteen a paired light and heavy scan


### Output
- Returns a tuple (heavy_scans, light_scans). Tuple{Vector{Int64}, Vector{Int64}}. These are of equal length and have
1:1 correspondance. heavy_scans[i] corresponds to light_scans[i]. The vectors contain scan indices. The i'th heavy
scan was nearest the i'th light scan, so the two are "paired". 

### Notes
    

### Examples 
julia> heavy = [1.0, 2.0, 3.0, 4.0]
4-element Vector{Float64}:
 1.0
 2.0
 3.0
 4.0

julia> light = [1.4, 2.6, 3.3, 3.8, 4.1]
5-element Vector{Float64}:
 1.4
 2.6
 3.3
 3.8
 4.1

julia> getScanPairs(heavy, light)
([1, 2, 3, 4, 5], [1, 3, 3, 4, 4])
"""
function getScanPairs(H::Vector{T}, L::Vector{T}; max_diff::T = T(3)) where T <: Number
    heavy_scans = Vector{Int64}()
    light_scans = Vector{Int64}()

    if ((length(L) < 3) | (length(H) < 3)) #Not enough scans
        return heavy_scans, light_scans
    end

    j = 1
    while j <= length(L) #Advance first light index until the retention time is greater than
        if L[j] > H[1]   #that for the first heavy scan. 
            break
        else
            j += 1
        end
    end 

    for i in eachindex(@view(H[1:end - 1]))

        midpoint = H[i] + (H[i+1] - H[i])/2.0 

        while (L[j] <= H[i+1]) & (L[j] >= H[i])
            if (L[j] < midpoint) & (abs(L[j] - H[i]) < max_diff) #L[j] is clossest to H[i]
                push!(light_scans, j), push!(heavy_scans, i) #Pair j with i
            elseif (abs(L[j] - H[i+1]) < max_diff) #L[j] is clossest to H[i]+1
                push!(light_scans, j), push!(heavy_scans, i+1) #Pair j with i + 1
            end
            j+=1
            if j > length(L) #Bounds check, all light scans exhausted
                return heavy_scans, light_scans
            end
        end
    end

    if (j>1)
        #The first light scan after the last heavy scan is a special case
        if (abs(L[j] - H[end]) < abs(L[j-1] - H[end])) & (abs(L[j] - H[end]) < max_diff)
            push!(light_scans, j), push!(heavy_scans, length(H))
        end
    end

    return heavy_scans, light_scans
end