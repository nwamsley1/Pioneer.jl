# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

function findLocalMinima(arr::AbstractVector{T})::Vector{Int64} where {T<:Real}
    # Returns indices of local minima in arr 
    # A value x[i] is in bin j if x_bins[j] ≤ x[i] < x_bins[j+1]
    # Values less than x_bins[1] get bin 0
    # Values ≥ x_bins[end] get bin length(x_bins)
    N = length(arr)
    if N < 3
        if (arr[2] < arr[3]) & (arr[2] < arr[1])
            return Int64[2]
        else
            return  Int64[]
        end
    end
    
    minima_indices = Int64[]
    
    for i in 2:(N-1)
        if (arr[i] <= arr[i-1] && arr[i] < arr[i+1]) ||  # Handle regular minima
           (arr[i] < arr[i-1] && arr[i] <= arr[i+1])     # Handle plateau cases
            push!(minima_indices, i)
        end
    end

    return minima_indices
end

function findLocalMaxima(arr::AbstractVector{T})::Vector{Int64} where {T<:Real}
    # Returns indices of local minima in arr 
    # A value x[i] is in bin j if x_bins[j] ≤ x[i] < x_bins[j+1]
    # Values less than x_bins[1] get bin 0
    # Values ≥ x_bins[end] get bin length(x_bins)
    N = length(arr)
    if N < 3
        if (arr[2] > arr[3]) & (arr[2] > arr[1])
            return Int64[2]
        else
            return  Int64[]
        end
    end
    
    maxima_indices = Int64[]
    
    for i in 2:(N-1)
        if (arr[i] >= arr[i-1] && arr[i] > arr[i+1]) ||  # Handle regular minima
           (arr[i] > arr[i-1] && arr[i] >= arr[i+1])     # Handle plateau cases
            push!(maxima_indices, i)
        end
    end

    return maxima_indices
end

function getBinIdx(x::AbstractVector{T}, x_bins::AbstractVector{U}) where {T,U<:Real}
    # Returns bin indices for each value in x
    # A value x[i] is in bin j if x_bins[j] ≤ x[i] < x_bins[j+1]
    # Values less than x_bins[1] get bin 0
    # Values ≥ x_bins[end] get bin length(x_bins)
    
    bins = zeros(Int, length(x))
    
    for (i, val) in enumerate(x)
        # searchsorted returns range of indices where value could be inserted
        # while maintaining sort order. We want the lower bound.
        bins[i] = searchsorted(x_bins, val).start - 1
    end
    
    return bins
end

function mergeBins(
    bin_vals::AbstractVector{<:Real}, 
    bin_idxs::AbstractVector{Int64},
    n::AbstractVector{Int64};
    min_bin_size::Int64 = 50,
    min_bin_width::Real = 0.1)

    #=
        Given N bins, each bin has values (bin_values), 
        an integer id (bin_idxs) and a size (n). Merge
        bins that are too small or not sufficiently seperated in value.
        Then return a dictionary mapping the old bin ids to new, merged bin ids. 
    =#
    bin_idx_to_merged_bin_idx = Dict{Int64, Int64}()
    N = length(bin_vals)

    # Input validation
    @assert length(bin_idxs) == N "Length mismatch in input arrays"
    @assert length(n) == N "Length mismatch in input arrays"


    # First pass: identify bins that need merging
    needs_merging = falses(N)
    for i in 1:N
        bin_width = (i == N) ? typemax(Float32) : bin_vals[i + 1] - bin_vals[i]
        needs_merging[i] = (n[i] < min_bin_size) || (bin_width < min_bin_width)
    end

    # Second pass: assign new bin indices
    new_bin_idx = 1
    for i in 1:N
        if needs_merging[i]
            if i < N  # Merge with next bin
                bin_idx_to_merged_bin_idx[bin_idxs[i]] = new_bin_idx
            else     # Last bin merges with previous
                bin_idx_to_merged_bin_idx[bin_idxs[i]] = new_bin_idx - 1
            end
        else
            bin_idx_to_merged_bin_idx[bin_idxs[i]] = new_bin_idx
            new_bin_idx += 1
        end
    end

    return bin_idx_to_merged_bin_idx

end

function estimateKdeBins(
    x::AbstractVector{<:AbstractFloat})
    #=
        1) Kernel density estimate for x using Silverman's results_folder
        2) Evalueate kde at 1000 evenly spaced grid
        3) Find local minima in the kde grid. These are the bin boundaries
    =#

    #Kernel density estimate for x
    x_kde = KernelDensity.kde(x, bandwidth = 0.01)
    #Evaluate kernel density on grid 
    kde_eval_points = LinRange(minimum(x), maximum(x), 1000)

    density_atx = KernelDensity.pdf(x_kde, (kde_eval_points))
    #Get Local Minima/x bin boundaries
    x_bin_boundaries = kde_eval_points[findLocalMinima(density_atx)]
    #p = plot(kde_eval_points, density_atx)
    #vline!(p, x_bin_boundaries)
    #display(p)
    #Pad boundaires 
    pushfirst!(x_bin_boundaries, typemin(Float32))
    push!(x_bin_boundaries, typemax(Float32))

    #Get bin_idx for each element in the x array
    return getBinIdx(x, x_bin_boundaries), length(x_bin_boundaries) - 2 + 1
end

# Method overload to handle vectors with missing values
function estimateKdeBins(x::AbstractVector{Union{Missing, T}}) where T<:AbstractFloat
    # Filter out missing values
    x_clean = collect(skipmissing(x))
    
    # Check if we have enough data for meaningful KDE
    if length(x_clean) < 10
        throw(ArgumentError("Insufficient non-missing data for KDE estimation: $(length(x_clean)) points"))
    end
    
    # Call the original method with clean data
    return estimateKdeBins(x_clean)
end

function MergeBins(isotopes_ratio_data::SubDataFrame, x0_lim::Tuple{<:Real,<:Real}; min_bin_size=300, min_bin_width=0.1, max_iterations=100)
    
    isotopes_ratio_data = copy(isotopes_ratio_data)
    x0_min, x0_max = max(first(x0_lim), minimum(isotopes_ratio_data[!,:x0])), min(last(x0_lim), maximum(isotopes_ratio_data[!,:x0]))
    
    filter!(x->(x.x0>x0_min)&(x.x0<x0_max), isotopes_ratio_data)
    min_bin_count = max(
        (ceil(Int64, abs(x0_max - x0_min)) - 1)*2-2, #Should have two bins per Thompson (1 m/z)
        3 #but at least 6 bins 
    )
    #Estimate KDE Bins 
    bin_idxs, n_bins = estimateKdeBins(isotopes_ratio_data[!,:x0])
    isotopes_ratio_data[!,:kde_bin_idx] = bin_idxs
    if n_bins < 4
        println("<4")
        #plot(isotopes_ratio_data[!,:x0], isotopes_ratio_data[!,:yt], seriestype=:scatter, alpha = 0.1, show = true)
    end
    #Use KDE bins or uniform bins? Base this decision on the median intra-bin standard deviation 
    #A better bin selection means less intra-bin variance. 
    if n_bins < min_bin_count #KDE estimation did not give enough bins. Default to uniform binning 
        uniform_x_bounds = collect(LinRange{Float32, Int64}(x0_min, x0_max, min_bin_count))
        push!(uniform_x_bounds, typemax(Float32))
        pushfirst!(uniform_x_bounds, typemin(Float32))
        isotopes_ratio_data[!,:bin_idx] = getBinIdx(isotopes_ratio_data[!,:x0], uniform_x_bounds)
        @user_warn "KDE bin estimation failed... Default to uniform binning $n_bins $min_bin_count"
    else
        #Estimate uniform bins 
        uniform_x_bounds = collect(LinRange{Float32, Int64}(x0_min, x0_max, n_bins+2))
        push!(uniform_x_bounds, typemax(Float32))
        pushfirst!(uniform_x_bounds, typemin(Float32))
        isotopes_ratio_data[!,:uniform_bin_idx] = getBinIdx(isotopes_ratio_data[!,:x0], uniform_x_bounds)
        #Use KDE bins or uniform bins? Base this decision on the mean standard deviation of a bin. 
        #A better bin selection means less within-bin variance. 
        uniform_bin_stds = filter!(!isnan, values(combine(x->std(x.x0), groupby(isotopes_ratio_data,:uniform_bin_idx))[:,2]))
        kde_bin_stds = filter!(!isnan, values(combine(x->std(x.x0), groupby(isotopes_ratio_data,:kde_bin_idx))[:,2]))
        mean_uniform_std = median(skipmissing(uniform_bin_stds))
        mean_kde_std = median(skipmissing(kde_bin_stds))
        if mean_uniform_std > mean_kde_std*0.9
            isotopes_ratio_data[!,:bin_idx] = isotopes_ratio_data[!,:uniform_bin_idx]
            @user_warn "KDE bin estimation failed... Default to uniform binning"
        else
            isotopes_ratio_data[!,:bin_idx] = isotopes_ratio_data[!,:kde_bin_idx]
        end
    end
    
    
    last_n_bins = 0
    for iter in 1:1
        # Compute bin statistics
        grouped_ratio_data = groupby(isotopes_ratio_data, :bin_idx)
        bin_centers = combine(grouped_ratio_data) do k_mean_bin
            (
                median_x0 = median(k_mean_bin[!,:x0]),
                median_yt = median(k_mean_bin[!,:yt]),
                median_x1 = median(k_mean_bin[!,:x1]),
                std_x0 = std(k_mean_bin[!,:x0]),
                n = size(k_mean_bin, 1),
                n_bins =length(grouped_ratio_data)
            )
        end
        #current_n_bins = size(bin_centers, 1)
        filter!(x->x.n>min_bin_size,bin_centers)
        return bin_centers
        # Check for convergence
        #=
        if current_n_bins == last_n_bins
            return bin_centers
            break
        end
        
        if iter == max_iterations
            @user_warn "Reached maximum iterations without convergence"
        end
        
        last_n_bins = current_n_bins
        
        # Sort and merge bins
        sort!(bin_centers, :bin_idx)
        bin_idx_to_merged_bin_idx = mergeBins(
            bin_centers.median_x0,
            bin_centers.bin_idx,
            bin_centers.n,
            min_bin_size=min_bin_size,
            min_bin_width=min_bin_width
        )
        
        # Update bin indices
        isotopes_ratio_data[!,:bin_idx] = [bin_idx_to_merged_bin_idx[idx] for idx in isotopes_ratio_data[!,:bin_idx]]
        =#
    end

    return nothing
end

function MergeBins(isotopes_ratio_data::DataFrame, x0_lim::Tuple{<:Real,<:Real}; min_bin_size=50, min_bin_width=0.1, max_iterations=100)
    ##########
    #Checks 
    missing_cols = filter(col -> !(col in Symbol.(names(isotopes_ratio_data))), [:x0, :x1, :yt, :prec_charge])
    if !isempty(missing_cols)
        error("Missing required columns: $(missing_cols)")
    end
    
    bined_ratio_data = combine(groupby(isotopes_ratio_data,:prec_charge)) do subdf
        mb = MergeBins(subdf,x0_lim, min_bin_size = min_bin_size, min_bin_width = min_bin_width, max_iterations = max_iterations)
        mb
    end
    filter!(x->x.n_bins>=2, bined_ratio_data)
    filter!(x->x.n>=min_bin_size, bined_ratio_data)
    return bined_ratio_data
end
