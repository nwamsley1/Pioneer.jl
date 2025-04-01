base_path = "/Users/nathanwamsley/Data/Mar_2025/Kevin_DE_Tag_Pioneer/KMD_DE_out_L2Test_lambda1en4/qc_plots"
using JLD2, Plots
A = load(joinpath(base_path, "A_mat.jld2"))["H"]
y = load(joinpath(base_path, "y_mat.jld2"))["y"]
residuals = load(joinpath(base_path, "residuals.jld2"))["residuals"]
weights = load(joinpath(base_path, "weights.jld2"))["weights"]
id_to_col = load(joinpath(base_path, "getIdToCol.jld2"))["id_to_col"]
_matches_ = load(joinpath(base_path, "ion_matches.jld2"))["_matches_"]
_misses_ = load(joinpath(base_path, "ion_misses.jld2"))["_misses_"]
using Plots
using Plots

using Plots

using Plots

using Plots

"""
    visualize_precursor_fit_ms_style(A, y, weights, id_to_col, precursor_ids, mz_array; kwargs...)
    
Visualize mass spectrometry data in the traditional spectral plot style with vertical lines.

# Arguments
- `A::Matrix{T}`: Template matrix where each column is a template spectrum
- `y::Vector{T}`: Observed spectrum (intensity values)
- `weights::Vector{T}`: Fitted weights for each column in A
- `id_to_col`: Mapping from precursor ID to column index in A
- `precursor_ids::Vector{<:Integer}`: List of precursor IDs to include in the visualization
- `mz_array`: Array of m/z values corresponding to the observed spectrum

# Optional Arguments
- `matches::Vector{Pioneer.PrecursorMatch{T}}`: List of matched peaks
- `misses::Vector{Pioneer.PrecursorMatch{T}}`: List of missed peaks
- `mz_range::Union{Nothing, Tuple{<:Real, <:Real}}=nothing`: Option to limit m/z range
- `normalize::Bool=false`: Whether to normalize intensities
- `show_individual::Bool=true`: Whether to show individual precursor contributions
- `title::String="Precursor Fit Visualization"`: Main plot title
- `annotation_size::Int=8`: Font size for peak annotations
- `linewidth::Float64=1.0`: Width of the vertical lines
- `observed_color::Symbol=:black`: Color for observed spectrum
- `model_color::Symbol=:red`: Color for model fit
- `annotate_matches::Bool=true`: Whether to annotate matched peaks
- `flip_model::Bool=true`: Whether to flip model fit to negative values
- `max_annotations::Int=50`: Maximum number of annotations to show

# Returns
- `p::Plots.Plot`: The generated visualization
"""
function visualize_precursor_fit_ms_style(
    A::Matrix{T},                      # Template matrix
    y::Vector{T},                      # Observed spectrum 
    weights::Vector{T},                # Fitted weights
    id_to_col,                         # Mapping from precursor ID to column index
    precursor_ids::Vector{<:Integer},  # List of precursor IDs to plot
    mz_array::AbstractVector;          # M/Z values for x-axis
    matches::Union{Nothing, Vector}=nothing, # Matched peaks (optional)
    misses::Union{Nothing, Vector}=nothing,  # Missed peaks (optional)
    mz_range::Union{Nothing, Tuple{<:Real, <:Real}}=nothing, # Option to limit m/z range
    normalize::Bool=false,             # Option to normalize intensities
    show_individual::Bool=false,       # Whether to show individual precursor contributions
    title::String="Precursor Fit Visualization", # Plot title
    annotation_size::Int=8,            # Font size for peak annotations
    linewidth::Float64=1.0,            # Width of the vertical lines
    observed_color::Symbol=:black,     # Color for observed spectrum
    model_color::Symbol=:red,          # Color for model fit
    annotate_matches::Bool=true,       # Whether to annotate matched peaks
    flip_model::Bool=true,             # Whether to flip model fit to negative values
    max_annotations::Int=50            # Maximum number of annotations to show
) where T <: AbstractFloat
    
    # Verify mz_array isn't empty
    if isempty(mz_array)
        error("mz_array cannot be empty")
    end
    
    # Determine the valid range for plotting (indices that correspond to valid m/z values)
    valid_mz_range = 1:min(length(mz_array), length(y))
    
    # Process each precursor ID
    valid_precursors = UInt32[]
    precursor_contributions = Dict{UInt32, Vector{T}}()
    model_prediction = zeros(T, length(y))
    
    # Dictionary to map precursor IDs to colors for consistency
    precursor_colors = Dict{UInt32, Symbol}()
    available_colors = [:green, :orange, :purple, :brown, :pink, :cyan, :magenta, :yellow]
    
    for (i, id) in enumerate(precursor_ids)
        try
            # Convert to UInt32 if needed
            prec_id = convert(UInt32, id)
            
            # Get column index for this precursor
            col = id_to_col[prec_id]
            if col > 0 && col <= size(A, 2)  # Ensure valid column index
                # Calculate this precursor's contribution
                precursor_weight = weights[col]
                contribution = A[:, col] * precursor_weight
                precursor_contributions[prec_id] = contribution
                model_prediction .+= contribution
                push!(valid_precursors, prec_id)
                
                # Assign color
                color_idx = ((i-1) % length(available_colors)) + 1
                precursor_colors[prec_id] = available_colors[color_idx]
                
                @info "Precursor ID 0x$(string(prec_id, base=16)) (column $col) with weight $precursor_weight added"
            else
                @warn "Precursor ID 0x$(string(prec_id, base=16)) maps to invalid column $col"
            end
        catch e
            @warn "Error processing precursor ID $id: $e"
        end
    end
    
    if isempty(valid_precursors)
        error("None of the provided precursor IDs were found in the data or mapped to valid columns")
    end
    
    # Filter matches and misses for the selected precursor IDs
    selected_matches = []
    if matches !== nothing
        selected_matches = filter(m -> 
            m.prec_id in valid_precursors && 
            m.intensity > 0 && 
            m.peak_ind > 0 && 
            m.peak_ind <= length(valid_mz_range), 
            matches)
    end
    
    # Filter observed spectrum to valid m/z range and determine plotting range
    valid_indices = intersect(valid_mz_range, 1:length(y))
    
    # If m/z range provided, use it directly
    if mz_range !== nothing
        min_mz, max_mz = mz_range
        # Find indices that fall within the m/z range
        indices = findall(mz -> min_mz <= mz <= max_mz, mz_array[valid_indices])
        if isempty(indices)
            error("No m/z values found in the specified range ($min_mz, $max_mz)")
        end
        # Map back to original indices
        plot_indices = valid_indices[indices]
    else
        # Use all valid indices
        plot_indices = valid_indices
    end
    
    # Sort selected matches by intensity for annotation prioritization
    if annotate_matches && !isempty(selected_matches)
        sort!(selected_matches, by = m -> m.intensity, rev = true)
    end
    
    # Prepare data for plotting
    y_observed = y[plot_indices]
    mz_plot = mz_array[plot_indices]
    model_values = model_prediction[plot_indices]
    
    # Normalize if requested
    if normalize
        max_y = maximum(y_observed)
        if max_y > 0
            y_observed = y_observed ./ max_y
            model_values = model_values ./ max_y
            
            for id in valid_precursors
                precursor_contributions[id] = precursor_contributions[id][plot_indices] ./ max_y
            end
        end
    end
    
    # If flipping model values, negate them
    if flip_model
        model_values = -1.0 .* model_values
        for id in valid_precursors
            precursor_contributions[id] = -1.0 .* precursor_contributions[id][plot_indices]
        end
    end
    
    # Determine the y-axis range
    y_max = maximum(y_observed) * 1.1
    y_min = flip_model ? minimum(model_values) * 1.1 : 0
    
    # Create the plot
    p = plot(fontfamily="helvetica", size=(1000, 600), title=title, 
             xlabel="m/z", ylabel="Intensity", legend=:topright,
             xlims=(minimum(mz_plot), maximum(mz_plot)), 
             ylims=(y_min, y_max))
    
    # Plot observed spectrum (vertical lines)
    for i in 1:length(mz_plot)
        if y_observed[i] > 0
            plot!(p, [mz_plot[i], mz_plot[i]], [0, y_observed[i]], 
                  linewidth=linewidth, color=observed_color, label=(i==1 ? "Observed" : nothing))
        end
    end
    
    # Plot model prediction (vertical lines)
    for i in 1:length(mz_plot)
        if flip_model
            if model_values[i] < 0  # When flipped, model values are negative
                plot!(p, [mz_plot[i], mz_plot[i]], [0, model_values[i]], 
                      linewidth=linewidth, color=model_color, label=(i==1 ? "Model Fit" : nothing))
            end
        else
            if model_values[i] > 0
                plot!(p, [mz_plot[i], mz_plot[i]], [0, model_values[i]], 
                      linewidth=linewidth, color=model_color, label=(i==1 ? "Model Fit" : nothing))
            end
        end
    end
    
    # Plot individual precursor contributions if requested
    if show_individual
        for (idx, id) in enumerate(valid_precursors)
            color = precursor_colors[id]
            contribution = precursor_contributions[id]
            
            for i in 1:length(mz_plot)
                value = contribution[i]
                if (flip_model && value < 0) || (!flip_model && value > 0)
                    plot!(p, [mz_plot[i], mz_plot[i]], [0, value], 
                          linewidth=linewidth, color=color, 
                          label=(i==1 ? "Precursor 0x$(string(id, base=16))" : nothing))
                end
            end
        end
    end
    
    # Add annotations for matched peaks
    if annotate_matches && !isempty(selected_matches)
        # Limit number of annotations to avoid cluttering
        annotation_count = min(length(selected_matches), max_annotations)
        
        for i in 1:annotation_count
            match = selected_matches[i]
            peak_idx = match.peak_ind
            
            # Skip if outside plot range
            if !(peak_idx in plot_indices) || peak_idx > length(mz_array)
                continue
            end
            
            mz = mz_array[peak_idx]
            intensity = match.intensity
            
            if normalize && maximum(y_observed) > 0
                intensity = intensity / maximum(y_observed)
            end
            
            # Add annotation with fragment name or precursor ID if name not available
            annotation_text = "0x$(string(match.prec_id, base=16))"
            
            # Place annotation above the peak
            annotate!(p, mz, intensity + 0.02*y_max, 
                      text(annotation_text, precursor_colors[match.prec_id], :top, :center, annotation_size))
        end
    end
    
    return p
end

"""
    visualize_chromatogram(
        matched_precursors, 
        precursor_ids,
        transition_names=nothing;
        title="Precursor Chromatogram",
        flip_plot=false
    )
    
Create a chromatogram visualization for the specified precursors.

# Arguments
- `matched_precursors`: Dictionary mapping precursor IDs to chromatogram data
- `precursor_ids::Vector{<:Integer}`: List of precursor IDs to include
- `transition_names=nothing`: Optional list of transition names to include (if not specified, all are shown)

# Optional Arguments
- `title::String="Precursor Chromatogram"`: Plot title
- `flip_plot::Bool=false`: If true, will create a paired plot with negative values
- `legend_position=:topright`: Position of the legend

# Returns
- `p::Plots.Plot`: The generated visualization
"""
function visualize_chromatogram(
    matched_precursors,            # Dictionary with chromatogram data
    precursor_ids::Vector{<:Integer}, # Precursor IDs to visualize
    transition_names=nothing;      # Optional transition names to include
    title::String="Precursor Chromatogram",
    flip_plot::Bool=false,         # Create paired plot if true
    legend_position=:topright
)
    # Initialize plot
    p = plot(fontfamily="helvetica", size=(1000, 600), title=title, 
             xlabel="Retention Time", ylabel="Intensity", legend=legend_position)
    
    # Process each precursor ID
    for (idx, id) in enumerate(precursor_ids)
        # Convert if needed and check if precursor exists
        prec_id = convert(UInt32, id)
        if !haskey(matched_precursors, prec_id)
            @warn "Precursor ID 0x$(string(prec_id, base=16)) not found in matched_precursors"
            continue
        end
        
        # Extract chromatogram data
        precursor_data = matched_precursors[prec_id]
        
        # Check if this precursor has transitions and retention times
        if !hasfield(typeof(precursor_data), :transitions) || !hasfield(typeof(precursor_data), :rts)
            @warn "Precursor ID 0x$(string(prec_id, base=16)) lacks chromatogram data"
            continue
        end
        
        transitions = precursor_data.transitions
        rts = precursor_data.rts
        
        if isempty(transitions) || isempty(rts)
            @warn "Precursor ID 0x$(string(prec_id, base=16)) has empty chromatogram data"
            continue
        end
        
        # Determine which transitions to plot
        transitions_to_plot = isnothing(transition_names) ? keys(transitions) : 
                           filter(t -> t in transition_names, keys(transitions))
        
        if isempty(transitions_to_plot)
            @warn "No matching transitions found for precursor ID 0x$(string(prec_id, base=16))"
            continue
        end
        
        # Generate colors
        available_colors = [:green, :orange, :purple, :brown, :pink, :cyan, :magenta, :yellow]
        
        # Find max intensity for scaling if flipping
        max_intensity = 0.0
        if flip_plot
            for t in transitions_to_plot
                if maximum(transitions[t]) > max_intensity
                    max_intensity = maximum(transitions[t])
                end
            end
        end
        
        # Plot each transition
        for (color_idx, t) in enumerate(transitions_to_plot)
            color = available_colors[(color_idx-1) % length(available_colors) + 1]
            
            # Skip transitions with no data
            if isempty(transitions[t]) || all(iszero, transitions[t])
                continue
            end
            
            # Get label
            label = "ID 0x$(string(prec_id, base=16)) - $t"
            
            # Plot with or without flipping
            if flip_plot && idx == 2  # Second precursor gets flipped
                plot!(p, rts, -1.0 * transitions[t], color=color, linewidth=1.5, label=label)
                scatter!(p, rts, -1.0 * transitions[t], color=color, markersize=3, label=nothing)
            else
                plot!(p, rts, transitions[t], color=color, linewidth=1.5, label=label)
                scatter!(p, rts, transitions[t], color=color, markersize=3, label=nothing)
            end
        end
    end
    
    return p
end
# Example usage:
# p = visualize_precursor_fit(A, y, weights, id_to_col, [306817, 306889], 
#                            mz_array=ms_table[:mz_array][scan_idx],
#                            matches=_matches_, misses=_misses_)
# display(p)
ms_table = Arrow.Table("/Users/nathanwamsley/Data/Mar_2025/Kevin_DE_Tag_Pioneer/DE_benchmark_016.arrow")
function plot_fit(matches, misses, y_intensity, y_mz, weights, id_to_col, precursor_ids; kwargs...)
    p = plot(title = "test plot", fontfamily="helvetica", size=(1000, 600), xlabel="m/z", ylabel="Intensity", legend=:topright)
    min_mz, max_mz = typemax(Float32), typemin(Float32)
    for (color_idx, id) in enumerate(precursor_ids)
        plot!(label = id)
        for match in matches
            if match.prec_id == id
                mz = match.theoretical_mz
                if mz<min_mz
                    min_mz = mz
                end
                if mz>max_mz
                    max_mz = mz
                end
                theoretical_intensity = match.predicted_intensity*weights[id_to_col[id]]
                observed_intensity = match.intensity
                plot!(p, [mz, mz], [0, observed_intensity], color=:grey, label=nothing, linewidth = 5, alpha = 0.5)
                plot!(p, [mz, mz], [0, theoretical_intensity], color=color_idx, label=nothing, linewidth = 5, alpha = 0.5)
            end
        end
        for match in misses
            if match.prec_id == id
                mz = match.theoretical_mz
                if mz<min_mz
                    min_mz = mz
                end
                if mz>max_mz
                    max_mz = mz
                end
                theoretical_intensity = match.predicted_intensity*weights[id_to_col[id]]
                plot!(p, [mz, mz], [0, theoretical_intensity], color=color_idx, label=nothing, linewidth = 5, alpha = 0.5)
                #plot!(p, [mz, mz], [0, -observed_intensity], color=color_idx, label=nothing)
            end
        end
        for i in 1:length(y_mz)
            if (y_mz[i] > min_mz - 2.0) && (y_mz[i] < max_mz + 2.0)
                plot!(p, [y_mz[i], y_mz[i]], [0, -y_intensity[i]], color=:grey, label=nothing, linewidth = 1, alpha = 0.5)
            end
        end
    end
    return p
end 
p = plot_fit(_matches_, _misses_, ms_table[:intensity_array][8241], ms_table[:mz_array][8241], weights, id_to_col, [0x0004ae81
0x0004aee1
0x0004aef4
0x0004ae8d
0x0004aea3
0x0004aeb7
0x0004aec9])                           
display(p)


p = plot_fit(_matches_, _misses_, ms_table[:intensity_array][6625], ms_table[:mz_array][6625], weights, id_to_col, 
[
    0x00028547,
    0x00028579,
    0x00028555,
    0x00028561,
    0x00028567,
    0x0002856f
])                           
display(p)


sub_matches = [match for match in _matches_ if match.prec_id ∈ [
        0x0004ae81
       0x0004aee1
       0x0004aef4
       0x0004ae8d
       0x0004aea3
       0x0004aeb7
       0x0004aec9]];
sub_misses = [match for match in _misses_ if match.prec_id ∈ [
        0x0004ae81
       0x0004aee1
       0x0004aef4
       0x0004ae8d
       0x0004aea3
       0x0004aeb7
       0x0004aec9]];
a = sub_matches[3].theoretical_mz
b = sub_misses[1].theoretical_mz
(a-b)/(a/1e6)
#How were the masses calcualted to begin with. Do we need higher precision?

#Difference more than 6ppm, so makes sense that one missed
#however, the theoretical mz should actually be the same. 

# Example usage:
#scan_idx = 8241
#ms_table = Arrow.Table("/Users/nathanwamsley/Data/Mar_2025/Kevin_DE_Tag_Pioneer/DE_benchmark_016.arrow")
#p = visualize_precursor_fit(A, y, weights, id_to_col, [306817, 306889], 
#                           mz_array=ms_table[:mz_array][8241],
#                           matches=_matches_, misses=_misses_,
#                           highlight_matches=true, show_misses=true)
#function plot_fit(A, y, weights, id_to_col, precursor_ids; kwargs...)
#    )
#end
# display(p)
# Example usage:
# p = visualize_precursor_fit(A, y, weights, id_to_col, [0x0004c443, 0x0004fac9], matches=_matches_)
# display(p)
#Example usage:
#ms_table = Arrow.Table("/Users/nathanwamsley/Data/Mar_2025/Kevin_DE_Tag_Pioneer/DE_benchmark_016.arrow")
#ms_table[:mz_array][8241]
#p = visualize_precursor_fit(A, y, weights, id_to_col, [306817,306889], 
#                          matches=_matches_, misses=_misses_,
#                          highlight_matches=true, show_misses=true)

#p = visualize_precursor_fit(A, y, weights, id_to_col, [306817,306889], matches=_matches_)
#display(p)d