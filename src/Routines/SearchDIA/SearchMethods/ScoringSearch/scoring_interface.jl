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

"""
Interface functions for ScoringSearch that work exclusively with file references.

These functions provide a clean abstraction layer between ScoringSearch and file operations,
ensuring all file access goes through the reference system.
"""
#==========================================================
Reference-Based Wrappers for Direct File Operations
==========================================================#

"""
    calculate_and_add_global_scores!(pg_refs::Vector{ProteinGroupFileReference})
    
Calculate global protein scores and add them to files via references.
Returns the score dictionary for downstream use.
"""
function calculate_and_add_global_scores!(pg_refs::Vector{ProteinGroupFileReference})
    sqrt_n_runs = max(1, floor(Int, sqrt(length(pg_refs))))
    acc_to_scores = Dict{ProteinKey, Vector{Float32}}()

    # First pass: collect scores per protein across all files
    for ref in pg_refs
        process_with_memory_limit(ref,
            batch -> begin
                for row in eachrow(batch)
                    key = ProteinKey(row.protein_name, row.target, row.entrap_id)
                    push!(get!(acc_to_scores, key, Float32[]), row.pg_score)
                end
            end
        )
    end

    # Compute global score using log-odds combination
    acc_to_global_score = Dict{ProteinKey, Float32}()
    for (key, scores) in acc_to_scores
        acc_to_global_score[key] = logodds(scores, sqrt_n_runs)
    end
    
    # Second pass: add global_pg_score column and sort
    for ref in pg_refs
        add_column_and_sort!(ref, :global_pg_score, 
            batch -> begin
                scores = Vector{Float32}(undef, nrow(batch))
                for i in 1:nrow(batch)
                    key = ProteinKey(
                        batch.protein_name[i],
                        batch.target[i],
                        batch.entrap_id[i]
                    )
                    scores[i] = get(acc_to_global_score, key, batch.pg_score[i])
                end
                scores
            end,
            :global_pg_score, :target;  # sort keys
            reverse=true
        )
    end
    
    return acc_to_global_score
end

"""
    add_trace_qvalues(fdr_scale_factor::Float32)

Add a column `:trace_qval` based on the `:prob` column using target/decoy q-values.
"""
function add_trace_qvalues(fdr_scale_factor::Float32)
    op = function(df)
        qvals = Vector{Float32}(undef, nrow(df))
        get_qvalues!(df.prob, df.target, qvals; fdr_scale_factor=fdr_scale_factor)
        df[!, :trace_qval] = qvals
        return df
    end
    return "add_trace_qvalues" => op
end



"""
    add_prec_prob(prob_col::Symbol)

Compute run-specific precursor probabilities from the given probability column.
"""
function add_prec_prob(prob_col::Symbol)
    op = function(df)
        transform!(groupby(df, [:precursor_idx, :ms_file_idx]),
                   prob_col => (p -> 1.0f0 - 0.000001f0 - exp(sum(log1p.(-p)))) => :prec_prob)
        return df
    end
    return "add_prec_prob" => op
end

"""
    logodds(probs::AbstractVector{T}, top_n::Int) where {T<:AbstractFloat}

Combine probabilities using a log-odds average. 
The final value is converted back to a probability via the logistic function.
"""
function logodds(probs::AbstractVector{T}, top_n::Int) where {T<:AbstractFloat}
    isempty(probs) && return 0.0f0
    n = min(length(probs), top_n)
    # Sort descending and select the top n probabilities
    sorted = sort(probs; rev=true)
    selected = sorted[1:n]
    eps = 1f-6
    # Convert to log-odds, clip to avoid Inf or negative contribution
    logodds = log.(clamp.(selected, 0.1f0, 1 - eps) ./ (1 .- clamp.(selected, 0.1f0, 1 - eps)))
    avg = sum(logodds) / n
    return 1.0f0 / (1 + exp(-avg))
end


#==========================================================
Reference-based Sort and Filter Functions
==========================================================#





"""
    apply_probit_scores!(pg_refs::Vector{ProteinGroupFileReference}, 
                        β_fitted::Vector{Float64}, feature_names::Vector{Symbol})
    
Apply probit regression scores to protein group files.
Note: This function is called from utils.jl and needs access to calculate_probit_scores.
"""
function apply_probit_scores!(pg_refs::Vector{ProteinGroupFileReference},
                             β_fitted::Vector{Float64},
                             feature_names::Vector{Symbol})
    for ref in pg_refs
        transform_and_write!(ref) do df
            # Calculate probit scores (function from utils.jl)
            X_file = Matrix{Float64}(df[:, feature_names])
            prob_scores = calculate_probit_scores(X_file, β_fitted)
            
            # Update scores
            df[!, :old_pg_score] = copy(df.pg_score)
            df[!, :pg_score] = Float32.(prob_scores)
            
            # Sort by new scores
            sort!(df, [:pg_score, :target], rev = [true, true])
            
            return df
        end
    end
end

#==========================================================
Scoring-Specific Pipeline Operations
==========================================================#

function apply_mbr_filter!(
    merged_df::DataFrame,
    params,
    fdr_scale_factor::Float32,
)
    # 1) identify transfer candidates
    candidate_mask = merged_df.MBR_transfer_candidate

    # 2) identify bad transfers
    is_bad_transfer = candidate_mask .& (
         (merged_df.target .& coalesce.(merged_df.MBR_is_best_decoy, false)) .| # T->D
         (merged_df.decoy .& .!coalesce.(merged_df.MBR_is_best_decoy, false)) # D->T
        #merged_df.decoy # D->D or T->D
    )

    # 3) compute threshold using the local bad_mask
    τ = get_ftr_threshold(
        merged_df.prob,
        merged_df.target,
        is_bad_transfer,
        params.max_MBR_false_transfer_rate;
        mask = candidate_mask,
    )

    @user_warn "τ is $τ at α = $(params.max_MBR_false_transfer_rate)"
    # 4) one fused pass to clamp probs
    merged_df._filtered_prob = ifelse.(
        candidate_mask .& (merged_df.prob .< τ),
        0.0f0,
        merged_df.prob,
    )

    # if downstream code expects a Symbol for the prob-column
    return :_filtered_prob
end


"""
    add_best_trace_indicator(isotope_type::IsotopeTraceType, best_traces::Set)

Add best trace indicator based on isotope trace type.
"""
function add_best_trace_indicator(isotope_type::IsotopeTraceType, best_traces::Set)
    op = function(df)  # df is passed by transform_and_write!
        if seperateTraces(isotope_type)
            # Extract columns with type assertions for performance
            precursor_idx_col = df.precursor_idx::AbstractVector{UInt32}
            isotopes_captured_col = df.isotopes_captured::AbstractVector{Tuple{Int8, Int8}}   
            
            # Efficient vectorized operation for separate traces
            df[!,:best_trace] = [
                (precursor_idx=precursor_idx_col[i], 
                 isotopes_captured=isotopes_captured_col[i]) ∈ best_traces
                for i in eachindex(precursor_idx_col)
            ]
        else
            # Group-based operation for combined traces
            transform!(groupby(df, :precursor_idx),
                      :prob => (p -> begin
                          best_idx = argmax(p)
                          result = falses(length(p))
                          result[best_idx] = true
                          result
                      end) => :best_trace)
        end
        return df
    end
    return "" => op
end

"""
    get_quant_necessary_columns() -> Vector{Symbol}

Get the standard columns needed for quantification analysis.
"""
function get_quant_necessary_columns()
    return [
        :precursor_idx,
        :global_prob,
        :prec_prob,
        :trace_prob,
        :global_qval,
        :run_specific_qval,
        :prec_mz,
        :pep,
        :weight,
        :target,
        :rt,
        :irt_obs,
        :missed_cleavage,
        :Mox,
        :isotopes_captured,
        :scan_idx,
        :entrapment_group_id,
        :ms_file_idx
    ]
end
