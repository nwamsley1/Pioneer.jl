using DataFrames

# Base type for EFDR methods
abstract type EFDRMethod end

# EFDR method structs with data fields
struct CombinedEFDR{T<:Real, I<:Integer} <: EFDRMethod
    score::Vector{T}
    original_target_score::Vector{T}
    entrapment_label::Vector{I}
    qval::Vector{T}
    r::T
    
    # Inner constructor to ensure vectors have same length
    function CombinedEFDR(score::Vector{T}, original_target_score::Vector{T}, 
                         entrapment_label::Vector{I}, qval::Vector{T}, r::T) where {T<:Real, I<:Integer}
        n = length(score)
        if length(original_target_score) != n || length(entrapment_label) != n || length(qval) != n
            error("All input vectors must have the same length")
        end
        new{T,I}(score, original_target_score, entrapment_label, qval, r)
    end
end

struct PairedEFDR{T<:Real, I<:Integer} <: EFDRMethod
    score::Vector{T}
    original_target_score::Vector{T}
    entrapment_label::Vector{I}
    qval::Vector{T}
    r::T
    
    # Inner constructor to ensure vectors have same length
    function PairedEFDR(score::Vector{T}, original_target_score::Vector{T}, 
                       entrapment_label::Vector{I}, qval::Vector{T}, r::T) where {T<:Real, I<:Integer}
        n = length(score)
        if length(original_target_score) != n || length(entrapment_label) != n || length(qval) != n
            error("All input vectors must have the same length")
        end
        new{T,I}(score, original_target_score, entrapment_label, qval, r)
    end
end

# Convenience constructors with default r value
CombinedEFDR(score, original_target_score, entrapment_label, qval; r=1.0) = 
    CombinedEFDR(convert(Vector{Float64}, score), convert(Vector{Float64}, original_target_score), 
                 convert(Vector{Int}, entrapment_label), convert(Vector{Float64}, qval), Float64(r))

PairedEFDR(score, original_target_score, entrapment_label, qval; r=1.0) = 
    PairedEFDR(convert(Vector{Float64}, score), convert(Vector{Float64}, original_target_score), 
               convert(Vector{Int}, entrapment_label), convert(Vector{Float64}, qval), Float64(r))

"""
    calculate_efdr(method::EFDRMethod)

Calculate empirical FDR using the data stored in the method struct.

# Methods
- `CombinedEFDR`: Standard combined empirical FDR calculation
- `PairedEFDR`: Paired empirical FDR that considers score relationships

# Returns
- Vector of empirical FDR values
"""
function calculate_efdr(method::CombinedEFDR)
    # Sort by q-value (ascending) first, then by score (descending) to break ties
    sort_indices = sortperm(collect(zip(-method.score, method.qval)))
    
    entrapment_fdr = zeros(eltype(method.qval), length(method.qval))
    Nτ = 0
    Nϵ = 0
    for i in sort_indices
        # Determine if this is a target or entrapment
        is_original_target = iszero(method.entrapment_label[i])
        
        # Calculate empirical FDR based on the sorted order
        if is_original_target
            Nτ += 1
        else
            Nϵ += 1
        end
        # Avoid division by zero and cap at 1.0
        if Nϵ + Nτ > 0
            entrapment_fdr[i] = min(1.0, (Nϵ*(1 + 1/method.r)) / (Nϵ + Nτ))
        else
            entrapment_fdr[i] = 0.0
        end
    end

    return entrapment_fdr 
end

function calculate_efdr(method::PairedEFDR)
    # Sort by q-value (ascending) first, then by score (descending) to break ties
    sort_indices = sortperm(collect(zip(-method.score, method.qval)))
    
    entrapment_fdr = zeros(eltype(method.qval), length(method.qval))

    total_ops = sum(1:length(sort_indices)) # Total operations for progress tracking
    pb = ProgressBar(total=total_ops)
    completed_ops = 0
    last_update = 0
    #Seems quadratic 
    for i in range(1,length(sort_indices))
        Nτ = 0
        Nϵ = 0
        Nϵsτ = 0
        Nϵτs = 0
        s = method.score[sort_indices[i]]   
        for j in 1:i
            # Determine if this is a target or entrapment
            is_original_target = method.entrapment_label[sort_indices[j]] == 0
            if is_original_target
                Nτ += 1
            else
                t = method.original_target_score[sort_indices[j]]
                e = method.score[sort_indices[j]]
                Nϵ += 1
                if (e >= s) & (t < s)
                    Nϵsτ += 1
                elseif  (e > t) & (t >= s)
                    Nϵτs += 1
                end
            end
            completed_ops += 1
            if completed_ops % 1000 == 0
                update(pb, completed_ops-last_update)
                last_update = completed_ops
            end
        end
        # Determine if this is a target or entrapment
        # Avoid division by zero and cap at 1.0
        if Nϵ + Nτ > 0
            entrapment_fdr[sort_indices[i]] = min(1.0, (Nϵ + Nϵsτ + 2*Nϵτs) / (Nϵ + Nτ))
        else
            entrapment_fdr[sort_indices[i]] = 0.0
        end
    end

    return entrapment_fdr 
end

# Backward compatibility functions
"""
    get_combined_efdr(score, original_target_score, entrapment_label, qval)

Calculate empirical FDR using entrapment pairs, sorting by q-value first and score second.

# Arguments
- `score`: Vector of scores for each peptide
- `original_target_score`: Vector of original target scores (from entrapment pairs)
- `entrapment_label`: Vector indicating entrapment group (0 = target, >0 = entrapment)
- `qval`: Vector of q-values for sorting

# Returns
- Vector of empirical FDR values
"""
function get_combined_efdr(
    score::AbstractVector{<:Real},
    original_target_score::AbstractVector{<:Real},
    entrapment_label::AbstractVector{<:Integer},
    qval::AbstractVector{<:Real},
    r::Real = 1.0
)
    method = CombinedEFDR(score, original_target_score, entrapment_label, qval; r=r)
    return calculate_efdr(method)
end

"""
    get_paired_efdr(score, original_target_score, entrapment_label, qval)

Calculate paired empirical FDR that considers score relationships between entrapment pairs.

# Arguments
- `score`: Vector of scores for each peptide  
- `original_target_score`: Vector of original target scores (from entrapment pairs)
- `entrapment_label`: Vector indicating entrapment group (0 = target, >0 = entrapment)
- `qval`: Vector of q-values for sorting
- `r`: Ratio parameter (default: 1.0)

# Returns
- Vector of empirical FDR values
"""
function get_paired_efdr(
    score::AbstractVector{<:Real},
    original_target_score::AbstractVector{<:Real},
    entrapment_label::AbstractVector{<:Integer},
    qval::AbstractVector{<:Real},
    r::Real = 1.0
)
    method = PairedEFDR(score, original_target_score, entrapment_label, qval; r=r)
    return calculate_efdr(method)
end

"""
    add_efdr_columns!(df::DataFrame, library_precursors::DataFrame;
                     method_types=[CombinedEFDR, PairedEFDR],
                     score_qval_pairs=[(:score, :qval)],
                     r=1.0)

Add empirical FDR columns to a DataFrame for multiple score/q-value pairs and methods.

# Arguments
- `df`: DataFrame containing precursor results with scores and q-values
- `library_precursors`: DataFrame with entrapment group information
- `method_types`: Vector of EFDR method types to apply (default: [CombinedEFDR, PairedEFDR])
- `score_qval_pairs`: Vector of (score_col, qval_col) tuples to process
- `r`: Ratio parameter for EFDR calculation (default: 1.0)

# Effects
Adds columns to `df` named as: score_col_methodname_efdr
For example: global_prob_combined_efdr, global_prob_paired_efdr
"""
function add_efdr_columns!(df::DataFrame, library_precursors::DataFrame;
                          method_types::Vector=[CombinedEFDR, PairedEFDR],
                          score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:score, :qval)],
                          r::Float64=1.0)
    
    # Check for required columns
    if !hasproperty(df, :precursor_idx)
        error("DataFrame must have :precursor_idx column")
    end
    if !hasproperty(library_precursors, :entrapment_group_id)
        error("library_precursors must have :entrapment_group_id column")
    end
    
    # Get entrapment labels for each precursor in results
    entrap_labels = [library_precursors.entrapment_group_id[pid] for pid in df.precursor_idx]
    
    # Process each score/qval pair
    for (score_col, qval_col) in score_qval_pairs
        # Check columns exist
        if !hasproperty(df, score_col)
            @warn "Column $score_col not found in DataFrame, skipping..."
            continue
        end
        if !hasproperty(df, qval_col)
            @warn "Column $qval_col not found in DataFrame, skipping..."
            continue
        end
        
        # Get original target score column name
        original_target_col = Symbol(String(score_col) * "_original_target")
        if !hasproperty(df, original_target_col)
            @warn "Column $original_target_col not found. Make sure to run add_original_target_scores! first."
            continue
        end
        
        # Extract vectors
        scores = Float64.(df[!, score_col])
        original_target_scores = Float64.(df[!, original_target_col])
        qvals = Float64.(df[!, qval_col])
        
        # Apply each method type
        for method_type in method_types
            # Determine method name for column
            method_name = if method_type == CombinedEFDR
                "combined"
            elseif method_type == PairedEFDR
                "paired"
            else
                lowercase(string(method_type))
            end
            
            # Create column name
            efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
            
            # Create method instance and calculate EFDR
            method = method_type(scores, original_target_scores, entrap_labels, qvals, r)
            efdr_values = calculate_efdr(method)
            
            # Add column to dataframe
            df[!, efdr_col] = efdr_values
        end
    end
    
    return nothing
end