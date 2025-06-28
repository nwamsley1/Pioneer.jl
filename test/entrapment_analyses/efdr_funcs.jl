# Base type for EFDR methods
abstract type EFDRMethod end

# EFDR method structs with data fields
struct CombinedEFDR{T<:Real} <: EFDRMethod
    score::Vector{T}
    original_target_score::Vector{T}
    entrapment_label::Vector{Int}
    qval::Vector{T}
    r::T
    
    # Inner constructor to ensure vectors have same length
    function CombinedEFDR(score::Vector{T}, original_target_score::Vector{T}, 
                         entrapment_label::Vector{Int}, qval::Vector{T}, r::T) where T<:Real
        n = length(score)
        if length(original_target_score) != n || length(entrapment_label) != n || length(qval) != n
            error("All input vectors must have the same length")
        end
        new{T}(score, original_target_score, entrapment_label, qval, r)
    end
end

struct PairedEFDR{T<:Real} <: EFDRMethod
    score::Vector{T}
    original_target_score::Vector{T}
    entrapment_label::Vector{Int}
    qval::Vector{T}
    r::T
    
    # Inner constructor to ensure vectors have same length
    function PairedEFDR(score::Vector{T}, original_target_score::Vector{T}, 
                       entrapment_label::Vector{Int}, qval::Vector{T}, r::T) where T<:Real
        n = length(score)
        if length(original_target_score) != n || length(entrapment_label) != n || length(qval) != n
            error("All input vectors must have the same length")
        end
        new{T}(score, original_target_score, entrapment_label, qval, r)
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
    sort_indices = sortperm(collect(zip(method.qval, -method.score)))
    
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
    sort_indices = sortperm(collect(zip(method.qval, -method.score)))
    
    entrapment_fdr = zeros(eltype(method.qval), length(method.qval))
    Nτ = 0
    Nϵ = 0
    Nϵsτ = 0
    Nϵτs = 0
    for i in sort_indices
        # Determine if this is a target or entrapment
        is_original_target = method.entrapment_label[i] == 0
        # Calculate empirical FDR based on the sorted order
        if is_original_target
            Nτ += 1
        else
            if method.original_target_score[i] < method.score[i]
                Nϵsτ += 1
            else 
                Nϵτs += 1
            end
            Nϵ += 1
        end
        # Avoid division by zero and cap at 1.0
        if Nϵ + Nτ > 0
            entrapment_fdr[i] = min(1.0, (Nϵ + Nϵsτ + 2*Nϵτs) / (Nϵ + Nτ))
        else
            entrapment_fdr[i] = 0.0
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
    add_efdr_columns!(prec_results::DataFrame, library_precursors::DataFrame; 
                      method_types=[CombinedEFDR, PairedEFDR],
                      score_qval_pairs=[(:global_prob, :global_qval), (:prec_prob, :qval)],
                      r=1.0)

Add empirical FDR columns for specified score/qval pairs using specified methods.

Creates columns named like:
- global_prob_combined_efdr
- global_prob_paired_efdr  
- prec_prob_combined_efdr
- prec_prob_paired_efdr

# Arguments
- `prec_results`: DataFrame with scoring results
- `library_precursors`: DataFrame with library information including entrapment_group_id
- `method_types`: Vector of EFDR method types to use (default: both CombinedEFDR and PairedEFDR)
- `score_qval_pairs`: Vector of (score_col, qval_col) tuples to process
- `r`: Ratio parameter for EFDR calculation (default: 1.0)
"""
function add_efdr_columns!(prec_results::DataFrame, library_precursors::DataFrame; 
                          method_types::Vector = [CombinedEFDR, PairedEFDR],
                          score_qval_pairs::Vector{Tuple{Symbol,Symbol}} = [(:global_prob, :global_qval), (:prec_prob, :qval)],
                          r::Real = 1.0)
    # Check for required columns
    if !hasproperty(prec_results, :precursor_idx)
        error("prec_results must have :precursor_idx column.")
    end
    if !hasproperty(library_precursors, :entrapment_group_id)
        error("library_precursors must have :entrapment_group_id column.")
    end
    
    # Add entrap_pair_ids if not already present
    if !hasproperty(prec_results, :entrap_pair_id)
        add_entrap_pair_ids!(prec_results, library_precursors)
    end
    
    # Add entrapment labels
    entrapment_labels = [library_precursors.entrapment_group_id[pid] for pid in prec_results.precursor_idx]
    
    # Process each score/qval pair
    score_cols = Symbol[]
    for (score_col, qval_col) in score_qval_pairs
        # Check columns exist
        if !hasproperty(prec_results, score_col)
            error("prec_results must have :$score_col column.")
        end
        if !hasproperty(prec_results, qval_col)
            error("prec_results must have :$qval_col column.")
        end
        push!(score_cols, score_col)
    end
    
    # Add original target scores for all score columns
    add_original_target_scores!(prec_results, library_precursors, score_cols)
    
    # Calculate EFDR for each method and score/qval pair
    for (score_col, qval_col) in score_qval_pairs
        original_target_col = Symbol(String(score_col) * "_original_target")
        
        for method_type in method_types
            # Determine column name based on method type
            method_name = if method_type == CombinedEFDR
                "combined"
            elseif method_type == PairedEFDR
                "paired"
            else
                error("Unknown EFDR method type: $method_type")
            end
            
            efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
            
            # Create method instance with data
            method = method_type(prec_results[!, score_col],
                               prec_results[!, original_target_col],
                               entrapment_labels,
                               prec_results[!, qval_col];
                               r=r)
            
            # Calculate EFDR
            efdr_values = calculate_efdr(method)
            
            # Add to dataframe
            prec_results[!, efdr_col] = efdr_values
        end
    end
    
    return nothing
end
