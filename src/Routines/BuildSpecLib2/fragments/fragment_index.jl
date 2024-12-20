
"""
Builds an index of fragments for efficient searching
"""
struct FragmentIndex
    fragments::Vector{PioneerFrag}
    mz_bins::Vector{Tuple{Float32,Float32}}
    rt_bins::Vector{Tuple{Float32,Float32}}
    bin_indices::Vector{UInt32}
end

"""
Create a new fragment index from a set of fragments
"""
function build_fragment_index(fragments::Vector{T},
                            mz_tolerance::Float32,
                            rt_tolerance::Float32) where T <: AbstractFragment
    
    # Sort fragments by m/z
    sorted_frags = sort(fragments, by = f -> f.mz)
    
    # Create m/z bins
    mz_bins = Vector{Tuple{Float32,Float32}}()
    current_min = sorted_frags[1].mz
    current_max = current_min
    
    for frag in sorted_frags
        if (frag.mz - current_min)/current_min * 1e6 > mz_tolerance
            push!(mz_bins, (current_min, current_max))
            current_min = frag.mz
        end
        current_max = frag.mz
    end
    push!(mz_bins, (current_min, current_max))
    
    # Sort by retention time within m/z bins
    rt_bins = Vector{Tuple{Float32,Float32}}()
    bin_indices = Vector{UInt32}()
    
    start_idx = 1
    for (min_mz, max_mz) in mz_bins
        bin_frags = filter(f -> min_mz ≤ f.mz ≤ max_mz, sorted_frags)
        sort!(bin_frags, by = f -> f.retention_time)
        
        # Create RT bins within m/z bin
        current_rt = bin_frags[1].retention_time
        rt_start_idx = start_idx
        
        for (i, frag) in enumerate(bin_frags)
            if frag.retention_time - current_rt > rt_tolerance
                push!(rt_bins, (current_rt, bin_frags[i-1].retention_time))
                push!(bin_indices, rt_start_idx)
                current_rt = frag.retention_time
                rt_start_idx = start_idx + i - 1
            end
        end
        
        push!(rt_bins, (current_rt, bin_frags[end].retention_time))
        push!(bin_indices, rt_start_idx)
        start_idx += length(bin_frags)
    end
    
    return FragmentIndex(sorted_frags, mz_bins, rt_bins, bin_indices)
end