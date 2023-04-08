function getBestTransitions(best_psm::NamedTuple{(:rt, :scan_idx, :name, :mz, :intensity), Tuple{Float32, Int64, Vector{String}, Vector{Float32}, Vector{Float32}}},
                            maximum_fragment_count::UInt8)
    #Restrictions 
    #Can't choose two charge states for the same transition 
    sorted_by_name = sortperm(best_psm[:name])
    #names = best_psm[:name][sorted_by_name]
    #intensities = est_psm[:intensity][sorted_by_name]
    ion_indices = Int64[]

    for i in eachindex(sorted_by_name)
        if i>1
            if best_psm[:name][sorted_by_name[i]][1:2] == best_psm[:name][sorted_by_name[i - 1]][1:2]
                if best_psm[:intensity][sorted_by_name[i]] > best_psm[:intensity][sorted_by_name[i - 1]]
                    ion_indices[end] = sorted_by_name[i]
                else
                    ion_indices[end] = sorted_by_name[i-1]
                end
            end
            push!(ion_indices, sorted_by_name[i])
            continue
        end
        push!(ion_indices, sorted_by_name[i])
    end
    best_sorted = sortperm(best_psm[:intensity][ion_indices], rev = true)
    return best_sorted[1:min(maximum_fragment_count, end)]
end