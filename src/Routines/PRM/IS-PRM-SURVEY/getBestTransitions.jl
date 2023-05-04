"""
    getBestTransitions(best_psm::NamedTuple{(:rt, :scan_idx, :name, :mz, :intensity), Tuple{Int, Int, Vector{String}, Vector{T}, Vector{T}}},
        maximum_fragment_count::I) where {T <: AbstractFloat, I <: Int}

Identifies the best subset of transitions (of size `maximum_fragment_count`) given transition names and intensities. 
Will not choose two of the same transiion at a different charge start (for example y3+1 and y3+2 will not both be selected)`. 
Returns a Vector{Int64} which are the indices of the best fragments.  

### Input

- `best_psm::NamedTuple{(:rt, :scan_idx, :name, :mz, :intensity), Tuple{Int, Int, Vector{String}, Vector{T}, Vector{T}}}` -- 
Contains the transition names, masses, and their intensities for the `best_psm`. 
- `maximum_fragment_count::I` -- Get this many of the best fragments from the `best_psm`. If there are not at least this many fragments, 
then all will be selected. 
### Output
- Fills the `id_to_prec` and `sorted_prec_ids` fields in a `PrecursorDatabase` 

### Examples 
```
test_best_psm = (rt = 10,
                scan_idx = 10,
                name = ["y3+2","y3+1", "y4+1"],
                mz = Float64[200, 100, 200],
                intensity = Float64[1000, 2000, 500]
                )
getBestTransitions(test_best_psm, 2)
```

In this example the fragments/transition y3+2, y3+1, and y4+1 have intensities 1000, 2000, and 500 respectively. 
Evemn though the y3+2 and y3+1 have the highest intensities, the best 2 transitions are y3+1 and y4+1. This is 
because we will not select two of the same transition at a different charge state (must pick either y3+1 or y3+2.)


"""
function getBestTransitions(best_psm::NamedTuple{(:rt, :scan_idx, :name, :mz, :intensity), Tuple{T, Int, Vector{String}, Vector{T}, Vector{T}}};
                            maximum_fragment_count::UInt8 = UInt8(5)) where T <: AbstractFloat

    if length(best_psm[:name]) == 0
        return Int64[]
    end

    #Sorts the names and returns the indices of the sorted peptide_list
    #sortperm("B","C","A") == [3, 1, 2]
    sorted_by_name = sortperm(best_psm[:name])
    ion_indices = Int64[sorted_by_name[1]]

    for i in 2:length(sorted_by_name)
        # If best_psm[:name][sorted_by_name[i]] == "y3+1" then best_psm[:name][sorted_by_name[i]][1:2] == "y3"
        # So this excludes the charge state information. 
        # Is the current ion the same as the previos ion excepting a difference in charge state?
        if (best_psm[:name][sorted_by_name[i]][1:2] == best_psm[:name][sorted_by_name[i - 1]][1:2])
            #If so, which has a greater intensity. 
            if best_psm[:intensity][sorted_by_name[i]] > best_psm[:intensity][sorted_by_name[i - 1]]
                ion_indices[end] = sorted_by_name[i]
            else
                ion_indices[end] = sorted_by_name[i-1]
            end
        else
            push!(ion_indices, sorted_by_name[i])
        end
    end
    #Ion indices 
    best_sorted = sortperm(best_psm[:intensity][ion_indices], rev = true)[1:min(maximum_fragment_count, end)]
    return ion_indices[best_sorted]
end