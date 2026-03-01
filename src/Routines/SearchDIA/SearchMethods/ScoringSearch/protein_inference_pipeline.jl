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

#==========================================================
Pipeline Operations for Protein Inference
==========================================================#

"""
    add_peptide_metadata(precursors::LibraryPrecursors)

Add peptide sequence and protein information from precursor library.
Enriches PSM data with metadata needed for inference.
"""
function add_peptide_metadata(precursors::LibraryPrecursors)
    desc = "add_peptide_metadata"
    
    op = function(df)
        # Get library data
        all_sequences = getSequence(precursors)::AbstractVector{String}
        all_accessions = getAccessionNumbers(precursors)::AbstractVector{String}
        all_is_decoys = getIsDecoy(precursors)::AbstractVector{Bool}
        all_entrap_ids = getEntrapmentGroupId(precursors)::AbstractVector{UInt8}
        
        # Extract precursor_idx column with type assertion for performance
        precursor_idx = df.precursor_idx::AbstractVector{UInt32}
        n_rows = length(precursor_idx)
        
        # Add columns if they don't exist
        if !hasproperty(df, :sequence)
            sequences = Vector{String}(undef, n_rows)
            for i in 1:n_rows
                sequences[i] = all_sequences[precursor_idx[i]]
            end
            df.sequence = sequences
        end
        
        if !hasproperty(df, :accession_numbers)
            accessions = Vector{String}(undef, n_rows)
            for i in 1:n_rows
                accessions[i] = all_accessions[precursor_idx[i]]
            end
            df.accession_numbers = accessions
        end
        
        # Add target/decoy status if needed
        if !hasproperty(df, :is_decoy)
            is_decoy_vec = Vector{Bool}(undef, n_rows)
            for i in 1:n_rows
                is_decoy_vec[i] = all_is_decoys[precursor_idx[i]]
            end
            df.is_decoy = is_decoy_vec
        end
        
        # Add entrapment group if needed
        if !hasproperty(df, :entrap_id)
            entrap_vec = Vector{UInt8}(undef, n_rows)
            for i in 1:n_rows
                entrap_vec[i] = all_entrap_ids[precursor_idx[i]]
            end
            df.entrap_id = entrap_vec
        end
        
        return df
    end
    
    return desc => op
end

"""
    validate_peptide_data()

Ensure all required columns exist for protein inference.
"""
function validate_peptide_data()
    desc = "validate_peptide_data"
    
    op = function(df)
        required_cols = [:sequence, :accession_numbers, :is_decoy, 
                        :entrap_id, :precursor_idx]
        
        missing_cols = setdiff(required_cols, Symbol.(names(df)))
        if !isempty(missing_cols)
            error("Missing required columns for protein inference: $missing_cols")
        end
        
        return df
    end
    
    return desc => op
end

#==========================================================
Inference Wrapper Functions
==========================================================#

"""
    apply_inference_to_dataframe(df::DataFrame, precursors::LibraryPrecursors)

Apply the core infer_proteins algorithm to a prepared DataFrame.
Returns structured inference results using type-safe protein inference.
"""
function apply_inference_to_dataframe(df::DataFrame, precursors::LibraryPrecursors)
    if nrow(df) == 0
        # Return empty result for empty input
        return InferenceResult(
            Dictionary{PeptideKey, ProteinKey}()
        )
    end
    
    # Extract unique peptide-protein pairs
    unique_pairs = unique(df, [:sequence, :accession_numbers, :is_decoy, :entrap_id])
    
    # Convert to format expected by infer_proteins
    proteins_vec = Vector{ProteinKey}(undef, nrow(unique_pairs))
    peptides_vec = Vector{PeptideKey}(undef, nrow(unique_pairs))
    
    for (i, row) in enumerate(eachrow(unique_pairs))
        proteins_vec[i] = ProteinKey(
            row.accession_numbers,
            !row.is_decoy,  # Convert is_decoy to is_target
            row.entrap_id
        )
        
        peptides_vec[i] = PeptideKey(
            row.sequence,
            !row.is_decoy,  # Convert is_decoy to is_target
            row.entrap_id
        )
    end
    
    # Call inference algorithm
    return infer_proteins(proteins_vec, peptides_vec)
end

#==========================================================
Post-Inference Pipeline Operations
==========================================================#

"""
    add_inferred_protein_column(inference_result::InferenceResult)

Add inferred protein group assignments to PSMs based on inference results.
"""
function add_inferred_protein_column(inference_result::InferenceResult)
    desc = "add_inferred_protein_column"
    
    op = function(df)
        # Extract columns with type assertions for performance
        sequences = df.sequence::AbstractVector{String}
        is_decoy = df.is_decoy::AbstractVector{Bool}
        entrap_ids = df.entrap_id::AbstractVector{UInt8}
        accession_numbers = df.accession_numbers::AbstractVector{String}
        
        n_rows = length(sequences)
        inferred_proteins = Vector{Union{Missing, String}}(undef, n_rows)
        
        for i in 1:n_rows
            pep_key = PeptideKey(
                sequences[i],
                !is_decoy[i],
                entrap_ids[i]
            )
            
            if haskey(inference_result.peptide_to_protein, pep_key)
                protein_key = inference_result.peptide_to_protein[pep_key]
                inferred_proteins[i] = protein_key.name
            else
                # Fallback to original protein
                inferred_proteins[i] = missing#accession_numbers[i]
            end
        end
        
        df.inferred_protein_group = inferred_proteins
        return df
    end
    
    return desc => op
end

"""
    add_quantification_flag(inference_result::InferenceResult)

Add flag indicating whether peptide should be used for quantification.

Peptides present in the inference result (unique peptides assigned to minimal protein set)
are marked as usable for quantification. Peptides not in the result (shared peptides or
filtered out) are marked as NOT usable for quantification.
"""
function add_quantification_flag(inference_result::InferenceResult)
    desc = "add_quantification_flag"

    op = function(df)
        # Extract columns with type assertions for performance
        sequences = df.sequence::AbstractVector{String}
        is_decoy = df.is_decoy::AbstractVector{Bool}
        entrap_ids = df.entrap_id::AbstractVector{UInt8}

        n_rows = length(sequences)
        use_for_quant = Vector{Bool}(undef, n_rows)

        for i in 1:n_rows
            pep_key = PeptideKey(
                sequences[i],
                !is_decoy[i],
                entrap_ids[i]
            )

            # Peptide is usable for quantification if it's in the inference results
            # (i.e., it's a unique peptide assigned to a protein)
            use_for_quant[i] = haskey(inference_result.peptide_to_protein, pep_key)
        end

        df.use_for_protein_quant = use_for_quant
        return df
    end

    return desc => op
end

function _median_abs_deviation(values::Vector{Float64})::Float64
    isempty(values) && return 0.0
    med = Statistics.median(values)
    abs_devs = Vector{Float64}(undef, length(values))
    @inbounds for i in eachindex(values)
        abs_devs[i] = abs(values[i] - med)
    end
    return Statistics.median(abs_devs)
end

"""
    estimate_weight_detection_model(df::DataFrame)

Estimate per-file weight calibration for protein coverage surprise features.
Uses unique inferred peptides with valid positive weight values.
"""
function estimate_weight_detection_model(df::DataFrame)
    default_model = (
        log_threshold = 0.0f0,
        sigma_log = 1.0f0,
        n_unique_peptides = 0,
        used_fallback = true,
        coverage_overdispersion_rho = 0.0f0,
        n_dispersion_proteins = 0,
        used_dispersion_fallback = true
    )

    required_cols = (:use_for_protein_quant, :sequence, :inferred_protein_group, :weight)
    if any(col -> !hasproperty(df, col), required_cols)
        return default_model
    end

    n_rows = nrow(df)
    if n_rows == 0
        return default_model
    end

    best_weight_by_protein_peptide = Dict{Tuple{String, String}, Float64}()

    for i in 1:n_rows
        if df.use_for_protein_quant[i] != true
            continue
        end

        protein_name = df.inferred_protein_group[i]
        if ismissing(protein_name)
            continue
        end

        weight = df.weight[i]
        if ismissing(weight)
            continue
        end

        weight_val = Float64(weight)
        if !isfinite(weight_val) || weight_val <= 0.0
            continue
        end

        key = (String(protein_name), String(df.sequence[i]))
        if haskey(best_weight_by_protein_peptide, key)
            if weight_val > best_weight_by_protein_peptide[key]
                best_weight_by_protein_peptide[key] = weight_val
            end
        else
            best_weight_by_protein_peptide[key] = weight_val
        end
    end

    if isempty(best_weight_by_protein_peptide)
        return default_model
    end

    log_weights = Float64[log(weight) for weight in values(best_weight_by_protein_peptide)]
    log_threshold = Float64(Statistics.quantile(log_weights, 0.05))

    protein_to_log_weights = Dict{String, Vector{Float64}}()
    for ((protein_name, _), weight) in best_weight_by_protein_peptide
        push!(get!(protein_to_log_weights, protein_name, Float64[]), log(weight))
    end

    pooled_residuals = Float64[]
    for log_vals in values(protein_to_log_weights)
        if length(log_vals) < 2
            continue
        end
        med = Statistics.median(log_vals)
        for v in log_vals
            push!(pooled_residuals, v - med)
        end
    end

    sigma_log = 1.0
    used_fallback = true
    if !isempty(pooled_residuals)
        mad_val = _median_abs_deviation(pooled_residuals)
        if isfinite(mad_val) && mad_val > 0.0
            sigma_log = clamp(mad_val / 0.6744897501960817, 0.25, 2.5)
            used_fallback = false
        end
    end

    return (
        log_threshold = Float32(log_threshold),
        sigma_log = Float32(sigma_log),
        n_unique_peptides = length(log_weights),
        used_fallback = used_fallback,
        coverage_overdispersion_rho = 0.0f0,
        n_dispersion_proteins = 0,
        used_dispersion_fallback = true
    )
end

"""
    estimate_weight_coverage_overdispersion(df::DataFrame, calibration::NamedTuple)

Estimate a single file-level beta-binomial overdispersion parameter for the
weight-based coverage surprise features. Uses target protein groups weighted by
their raw initial `pg_score`.
"""
function estimate_weight_coverage_overdispersion(df::DataFrame, calibration::NamedTuple)
    default_result = (
        coverage_overdispersion_rho = 0.0f0,
        n_dispersion_proteins = 0,
        used_dispersion_fallback = true
    )

    required_cols = (:target, :pg_score, :top_pep_weight, :n_possible_peptides, :n_peptides)
    if any(col -> !hasproperty(df, col), required_cols) || nrow(df) == 0
        return default_result
    end

    log_threshold = Float64(hasproperty(calibration, :log_threshold) ? calibration.log_threshold : 0.0f0)
    sigma_log = Float64(hasproperty(calibration, :sigma_log) ? calibration.sigma_log : 1.0f0)
    sigma_log = (isfinite(sigma_log) && sigma_log > 0.0) ? sigma_log : 1.0
    std_normal = Distributions.Normal()

    numerator = 0.0
    denominator = 0.0
    n_dispersion_proteins = 0

    for i in 1:nrow(df)
        if df.target[i] != true
            continue
        end

        score_val = Float64(df.pg_score[i])
        if !isfinite(score_val) || score_val <= 0.0
            continue
        end

        top_weight = Float64(df.top_pep_weight[i])
        if !isfinite(top_weight) || top_weight <= 0.0
            continue
        end

        N_total = max(Int(df.n_possible_peptides[i]), 0)
        k_obs = max(Int(df.n_peptides[i]), 0)
        N_add = max(N_total - 1, 0)
        k_add = max(k_obs - 1, 0)
        if N_add < 2
            continue
        end

        z_norm = (log(top_weight) - log_threshold) / sigma_log
        p_other = clamp(Distributions.cdf(std_normal, z_norm / sqrt(2.0)), 1e-6, 1.0 - 1e-6)
        binomial_variance = N_add * p_other * (1.0 - p_other)
        if !isfinite(binomial_variance) || binomial_variance <= 1e-6
            continue
        end

        expected = N_add * p_other
        residual_sq = (k_add - expected)^2

        numerator += score_val * (residual_sq - binomial_variance)
        denominator += score_val * binomial_variance * (N_add - 1)
        n_dispersion_proteins += 1
    end

    if n_dispersion_proteins < 5 || !isfinite(denominator) || denominator <= 0.0
        return default_result
    end

    rho = numerator / denominator
    if !isfinite(rho)
        return default_result
    end

    return (
        coverage_overdispersion_rho = Float32(clamp(rho, 0.0, 0.2)),
        n_dispersion_proteins = n_dispersion_proteins,
        used_dispersion_fallback = false
    )
end

"""
    add_weight_observation_features(calibration)

Add weight-based protein coverage surprise features using a log-normal observation model.
"""
function add_weight_observation_features(calibration::NamedTuple)
    desc = "add_weight_observation_features"

    op = function(df)
        if !hasproperty(df, :protein_name)
            return df
        end

        n_rows = nrow(df)
        coverage_miss_pval = Vector{Float32}(undef, n_rows)
        coverage_miss_surprisal = Vector{Float32}(undef, n_rows)
        coverage_deficit_z = Vector{Float32}(undef, n_rows)
        top_weight_vs_threshold_z = Vector{Float32}(undef, n_rows)

        log_threshold = Float64(hasproperty(calibration, :log_threshold) ? calibration.log_threshold : 0.0f0)
        sigma_log = Float64(hasproperty(calibration, :sigma_log) ? calibration.sigma_log : 1.0f0)
        sigma_log = (isfinite(sigma_log) && sigma_log > 0.0) ? sigma_log : 1.0
        rho = Float64(hasproperty(calibration, :coverage_overdispersion_rho) ? calibration.coverage_overdispersion_rho : 0.0f0)
        rho = clamp((isfinite(rho) && rho >= 0.0) ? rho : 0.0, 0.0, 0.2)
        std_normal = Distributions.Normal()

        has_top_weight_col = hasproperty(df, :top_pep_weight)
        has_n_possible_col = hasproperty(df, :n_possible_peptides)
        has_n_peptides_col = hasproperty(df, :n_peptides)

        for i in 1:n_rows
            valid_top_weight = false
            top_weight = 0.0
            if has_top_weight_col
                weight_value = df.top_pep_weight[i]
                if !ismissing(weight_value)
                    top_weight = Float64(weight_value)
                    valid_top_weight = isfinite(top_weight) && top_weight > 0.0
                end
            end

            N_total = if has_n_possible_col
                max(Int(df.n_possible_peptides[i]), 0)
            else
                0
            end
            k_obs = if has_n_peptides_col
                max(Int(df.n_peptides[i]), 0)
            else
                0
            end

            N_add = max(N_total - 1, 0)
            k_add = max(k_obs - 1, 0)

            if !(valid_top_weight && isfinite(top_weight) && top_weight > 0.0) || N_add == 0
                coverage_miss_pval[i] = 1.0f0
                coverage_miss_surprisal[i] = 0.0f0
                coverage_deficit_z[i] = 0.0f0
                top_weight_vs_threshold_z[i] = 0.0f0
                continue
            end

            z_norm = (log(top_weight) - log_threshold) / sigma_log
            z_pair = z_norm / sqrt(2.0)
            p_other = clamp(Distributions.cdf(std_normal, z_pair), 1e-6, 1.0 - 1e-6)

            expected = N_add * p_other
            binomial_variance = N_add * p_other * (1.0 - p_other)

            pval = if rho > 1e-8 && N_add > 1
                concentration = (1.0 - rho) / rho
                alpha = max(p_other * concentration, 1e-6)
                beta = max((1.0 - p_other) * concentration, 1e-6)
                Distributions.cdf(Distributions.BetaBinomial(N_add, alpha, beta), min(k_add, N_add))
            else
                Distributions.cdf(Distributions.Binomial(N_add, p_other), min(k_add, N_add))
            end
            variance = (binomial_variance * (1.0 + (N_add - 1) * rho)) + 1e-6

            coverage_miss_pval[i] = Float32(pval)
            coverage_miss_surprisal[i] = Float32(-log10(max(pval, 1e-12)))
            coverage_deficit_z[i] = Float32((k_add - expected) / sqrt(variance))
            top_weight_vs_threshold_z[i] = Float32(z_norm)
        end

        df.coverage_miss_pval = coverage_miss_pval
        df.coverage_miss_surprisal = coverage_miss_surprisal
        df.coverage_deficit_z = coverage_deficit_z
        df.top_weight_vs_threshold_z = top_weight_vs_threshold_z
        return df
    end

    return desc => op
end

"""
    add_pg_score_interaction_features()

Add pg_score interaction terms for the current coverage-surprise feature family.
"""
function add_pg_score_interaction_features()
    desc = "add_pg_score_interaction_features"

    op = function(df)
        if !hasproperty(df, :pg_score)
            return df
        end

        n_rows = nrow(df)
        pg_score_x_coverage_miss_surprisal = Vector{Float32}(undef, n_rows)
        pg_score_x_coverage_deficit_z = Vector{Float32}(undef, n_rows)
        pg_score_x_top_weight_vs_threshold_z = Vector{Float32}(undef, n_rows)

        has_surprisal = hasproperty(df, :coverage_miss_surprisal)
        has_deficit = hasproperty(df, :coverage_deficit_z)
        has_threshold = hasproperty(df, :top_weight_vs_threshold_z)

        for i in 1:n_rows
            score = df.pg_score[i]
            pg_score_val = if ismissing(score)
                0.0f0
            else
                score_val = Float32(score)
                isfinite(score_val) ? score_val : 0.0f0
            end

            pg_score_x_coverage_miss_surprisal[i] = has_surprisal ? pg_score_val * Float32(df.coverage_miss_surprisal[i]) : 0.0f0
            pg_score_x_coverage_deficit_z[i] = has_deficit ? pg_score_val * Float32(df.coverage_deficit_z[i]) : 0.0f0
            pg_score_x_top_weight_vs_threshold_z[i] = has_threshold ? pg_score_val * Float32(df.top_weight_vs_threshold_z[i]) : 0.0f0
        end

        df.pg_score_x_coverage_miss_surprisal = pg_score_x_coverage_miss_surprisal
        df.pg_score_x_coverage_deficit_z = pg_score_x_coverage_deficit_z
        df.pg_score_x_top_weight_vs_threshold_z = pg_score_x_top_weight_vs_threshold_z
        return df
    end

    return desc => op
end

"""
    group_psms_by_protein(df::DataFrame)

Transform PSMs into protein groups by aggregating peptides.
Returns a DataFrame with one row per protein group.
"""
function group_psms_by_protein(df::DataFrame)
    if nrow(df) == 0
        # Return empty protein groups DataFrame with expected schema
        return DataFrame(
            protein_name = String[],
            target = Bool[],
            entrap_id = UInt8[],
            n_peptides = Int64[],
            peptide_list = String[],
            pg_score = Float32[],
            any_common_peps = Bool[],
            top_pep_weight = Float32[]
        )
    end

    # Determine which probability column to use for protein scoring
    if hasproperty(df, :MBR_boosted_prec_prob)
        prob_col = :MBR_boosted_prec_prob
    else
        prob_col = :prec_prob
    end
    has_weight = hasproperty(df, :weight)

    # Group by protein
    grouped = groupby(df, [:inferred_protein_group, :target, :entrap_id])
    
    # Aggregate to protein groups
    protein_groups = combine(grouped) do gdf
        quant_mask = gdf.use_for_protein_quant .== true
        # Get unique peptides that are used for quantification
        quant_peptides = unique(gdf[quant_mask, :sequence])
        n_peptides = length(quant_peptides)
        
        # Calculate initial protein score (log-sum)
        peptide_probs = gdf[quant_mask, prob_col]
        if isempty(peptide_probs)
            pg_score = 0.0f0
        else
            # Use best probability per peptide
            unique_pep_probs = Float32[]
            for pep in quant_peptides
                pep_mask = (gdf.sequence .== pep) .& quant_mask
                if any(pep_mask)
                    push!(unique_pep_probs, maximum(gdf[pep_mask, prob_col]))
                end
            end
            pg_score = -sum(log.(1.0f0 .- unique_pep_probs))
        end        

        top_pep_weight = 0.0f0
        if has_weight
            best_weight_by_peptide = Dict{String, Float32}()
            for i in 1:nrow(gdf)
                if quant_mask[i] != true
                    continue
                end
                weight = gdf.weight[i]
                if ismissing(weight)
                    continue
                end
                weight_val = Float32(weight)
                if !isfinite(weight_val) || weight_val <= 0.0f0
                    continue
                end
                pep = gdf.sequence[i]
                if haskey(best_weight_by_peptide, pep)
                    if weight_val > best_weight_by_peptide[pep]
                        best_weight_by_peptide[pep] = weight_val
                    end
                else
                    best_weight_by_peptide[pep] = weight_val
                end
            end

            if !isempty(best_weight_by_peptide)
                top_pep_weight = maximum(values(best_weight_by_peptide))
            end
        end

        has_common = any(
            quant_mask .&
            (gdf.missed_cleavage .== 0) .&
            (gdf.Mox .== 0)
        )
        
        DataFrame(
            n_peptides = n_peptides,
            peptide_list = join(quant_peptides, ";"),
            pg_score = pg_score,
            any_common_peps = has_common,
            top_pep_weight = top_pep_weight
        )
    end
    
    # Rename the grouping column
    rename!(protein_groups, :inferred_protein_group => :protein_name)
    
    return protein_groups
end

"""
    filter_by_min_peptides(min_peptides::Int)

Filter protein groups that don't meet minimum peptide requirement.
"""
function filter_by_min_peptides(min_peptides::Int)
    desc = "filter_by_min_peptides(min=$min_peptides)"
    
    op = function(df)
        if hasproperty(df, :n_peptides)
            # For protein groups
            filter!(row -> row.n_peptides >= min_peptides, df)
        end
        return df
    end
    
    return desc => op
end

"""
    add_protein_features(protein_catalog::Dict)

Add protein-level features like peptide coverage.
"""
function add_protein_features(protein_catalog::Dict)
    desc = "add_protein_features"
    
    op = function(df)
        if !hasproperty(df, :protein_name)
            return df
        end
        
        n_rows = nrow(df)
        n_possible = Vector{Int64}(undef, n_rows)
        peptide_coverage = Vector{Float32}(undef, n_rows)
        grouped_catalog_cache = Dict{
            @NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8},
            Set{String}
        }()
        
        for i in 1:n_rows
            key = (
                protein_name = df.protein_name[i],
                target = df.target[i],
                entrap_id = df.entrap_id[i]
            )

            possible_peptides = if haskey(protein_catalog, key)
                protein_catalog[key]
            elseif haskey(grouped_catalog_cache, key)
                grouped_catalog_cache[key]
            elseif occursin(';', key.protein_name)
                merged_peptides = Set{String}()
                for member in split(key.protein_name, ';')
                    member_key = (
                        protein_name = strip(member),
                        target = key.target,
                        entrap_id = key.entrap_id
                    )
                    if haskey(protein_catalog, member_key)
                        union!(merged_peptides, protein_catalog[member_key])
                    end
                end
                grouped_catalog_cache[key] = merged_peptides
                merged_peptides
            else
                Set{String}()
            end

            n_possible[i] = length(possible_peptides)

            if n_possible[i] > 0
                peptide_coverage[i] = Float32(df.n_peptides[i]) / Float32(n_possible[i])
            else
                peptide_coverage[i] = 0.0f0
            end
        end
        
        df.n_possible_peptides = n_possible
        df.peptide_coverage = peptide_coverage
        
        return df
    end
    
    return desc => op
end

#==========================================================
High-Level Interface
==========================================================#

"""
    perform_protein_inference_pipeline(psm_refs, output_folder, precursors, protein_catalog; kwargs...)

Perform protein inference using a composable pipeline approach.

# Arguments
- `psm_refs`: Vector of PSM file references
- `output_folder`: Directory for protein group output
- `precursors`: Library precursor information
- `protein_catalog`: Pre-computed protein-to-peptide mappings
- `min_peptides`: Minimum peptides per protein group (default: 2)

# Returns
- `pg_refs`: Vector of protein group file references
- `psm_to_pg_mapping`: Dictionary mapping PSM paths to protein group paths
"""
function perform_protein_inference_pipeline(
    psm_refs::Vector{PSMFileReference},
    output_folder::String,
    precursors::LibraryPrecursors,
    protein_catalog::Dict;
    min_peptides::Int = 2
)
    # Ensure output folder exists
    !isdir(output_folder) && mkpath(output_folder)
    
    # Build pre-inference pipeline
    pre_inference_pipeline = TransformPipeline() |>
        add_peptide_metadata(precursors) |>
        validate_peptide_data()
    
    # Process each file
    pg_refs = ProteinGroupFileReference[]
    psm_to_pg_mapping = Dict{String, String}()
    
    for (idx, psm_ref) in ProgressBar(collect(enumerate(psm_refs)))
        if !exists(psm_ref)
            continue
        end
        
        # Step 1: Apply pre-inference pipeline to add necessary columns
        apply_pipeline!(psm_ref, pre_inference_pipeline)
        
        # Step 2: Run inference on the prepared PSMs
        prepared_df = load_dataframe(psm_ref)
        inference_result = apply_inference_to_dataframe(prepared_df, precursors)
        
        # Step 3: Update PSMs with inference results
        psm_update_pipeline = TransformPipeline() |>
            add_inferred_protein_column(inference_result) |> #Protein group assigned to each peptide
            add_quantification_flag(inference_result) #Whether or not the peptide should be used for protein quant and inference (non ambiguous)
        
        # Apply updates to the same PSM file (which now has necessary columns)
        apply_pipeline!(psm_ref, psm_update_pipeline)
        
        # Step 4: Create protein groups
        # Reload updated PSMs
        updated_psms = load_dataframe(psm_ref)
        weight_calibration = estimate_weight_detection_model(updated_psms)
        
        # Group by protein
        protein_groups_df = group_psms_by_protein(updated_psms)
        (_, protein_feature_op) = add_protein_features(protein_catalog)
        protein_groups_df = protein_feature_op(protein_groups_df)
        coverage_overdispersion = estimate_weight_coverage_overdispersion(protein_groups_df, weight_calibration)
        weight_calibration = merge(weight_calibration, coverage_overdispersion)

        # Build post-inference pipeline
        post_inference_pipeline = TransformPipeline() |>
            filter_by_min_peptides(min_peptides) |>
            add_weight_observation_features(weight_calibration) |>
            add_pg_score_interaction_features()

        @user_info "Weight protein coverage calibration file_idx=$(idx) n_unique_peptides=$(weight_calibration.n_unique_peptides) log_threshold=$(weight_calibration.log_threshold) sigma_log=$(weight_calibration.sigma_log) rho=$(weight_calibration.coverage_overdispersion_rho) n_dispersion_proteins=$(weight_calibration.n_dispersion_proteins) used_fallback=$(weight_calibration.used_fallback) used_dispersion_fallback=$(weight_calibration.used_dispersion_fallback)"
        
        # Apply post-processing
        initial_rows = nrow(protein_groups_df)
        for (desc, op) in post_inference_pipeline.operations
            protein_groups_df = op(protein_groups_df)
            #@debug_l1 "Pipeline operation on protein groups" operation=desc rows_before=initial_rows rows_after=nrow(protein_groups_df)
            initial_rows = nrow(protein_groups_df)
        end
        
        # Write protein groups
        pg_filename = "protein_groups_$(lpad(idx, 3, '0')).arrow"
        pg_path = joinpath(output_folder, pg_filename)

        if nrow(protein_groups_df) > 0
            # Add file_idx for later merge/split operations
            protein_groups_df[!, :file_idx] = fill(Int64(idx), nrow(protein_groups_df))
            writeArrow(pg_path, protein_groups_df)
            pg_ref = ProteinGroupFileReference(pg_path)
            push!(pg_refs, pg_ref)
            psm_to_pg_mapping[file_path(psm_ref)] = pg_path
        else
            @user_warn "No protein groups to write after filtering" file_idx=idx pg_path=pg_path
        end
    end
    
    return pg_refs, psm_to_pg_mapping
end
