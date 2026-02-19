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
Trace Selection  
==========================================================#

"""
    get_best_traces(second_pass_psms_paths::Vector{String}, min_prob::Float32=0.75f0)
    -> Set{@NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}}

Identify best scoring isotope traces for each precursor.

# Process
1. Accumulates scores across files using MBR_boosted_trace_prob if available, otherwise trace_prob
2. Selects highest scoring trace per precursor
3. Returns set of best precursor-isotope combinations
"""
function get_best_traces(
    second_pass_psms_paths::Vector{String},
    min_prob::Float32 = 0.75f0
)

    #The sum of scores for a given precursor trace (precursor_idx and isotopes_captured) accross the 
    #entire experiment (all runs)!
    psms_trace_scores = Dictionary{
            @NamedTuple{precursor_idx::UInt32, isotopes_captured::Tuple{Int8, Int8}}, Float32}()

    # Track aggregated stats
    total_rows_processed = 0
    files_processed = 0
    
    for (file_idx, file_path) in enumerate(second_pass_psms_paths)

        if splitext(file_path)[end] != ".arrow"
            continue
        end
        
        row_score = zero(Float32)
        psms_table = Arrow.Table(file_path)
        n_rows = length(psms_table[1])
        total_rows_processed += n_rows
        files_processed += 1

        # Use MBR_boosted_trace_prob if available, otherwise trace_prob
        use_mbr_column = hasproperty(psms_table, :MBR_boosted_trace_prob)
        score_column = use_mbr_column ? :MBR_boosted_trace_prob : :trace_prob

        for i in range(1, n_rows)
            psms_key = (precursor_idx = psms_table[:precursor_idx][i],  isotopes_captured = psms_table[:isotopes_captured][i])

            row_score = psms_table[score_column][i]
            if haskey(psms_trace_scores, psms_key)
                psms_trace_scores[psms_key] = psms_trace_scores[psms_key] + row_score
            else
                insert!(
                    psms_trace_scores,
                    psms_key,
                    row_score
                )
            end
        end
    end

    #Convert the dictionary to a DataFrame 
    psms_trace_df = DataFrame(
        (precursor_idx = [key[:precursor_idx] for key in keys(psms_trace_scores)],
        isotopes_captured = [key[:isotopes_captured] for key in keys(psms_trace_scores)],
        score = [val for val in values(psms_trace_scores)])
        );
    #Now retain only the very best trace!
    psms_trace_df[!,:best_trace] .= false;
    gpsms = groupby(psms_trace_df,:precursor_idx)
    for (precursor_idx, psms) in pairs(gpsms)
        psms[argmax(psms[!,:score]),:best_trace] = true
    end
    filter!(x->x.best_trace, psms_trace_df);
    traces_passing = Set([(precursor_idx = x.precursor_idx, isotopes_captured = x.isotopes_captured) for x in eachrow(psms_trace_df)]);
    return traces_passing
end


"""
    adjust_any_common_peps!(feature_names::Vector{Symbol}, df::AbstractDataFrame)

Remove :any_common_peps from `feature_names` if the column is constant.
This avoids singular matrices in probit regression when the feature has
the same value for all rows.
"""
function adjust_any_common_peps!(feature_names::Vector{Symbol}, df::AbstractDataFrame)
    if :any_common_peps in feature_names && hasproperty(df, :any_common_peps)
        col = df.any_common_peps
        if all(col) || all(.!col)
            filter!(x -> x != :any_common_peps, feature_names)
            @user_warn "Removed constant :any_common_peps feature" constant_value = first(col)
        end
    end
    return feature_names
end

"""
    remove_zero_variance_columns!(feature_names::Vector{Symbol}, df::AbstractDataFrame)

Remove columns with zero variance from `feature_names` to prevent singular matrices
in probit regression. This includes constant columns, columns with all missing values,
and columns containing Inf/NaN values.

# Arguments
- `feature_names`: Vector of feature symbols to filter
- `df`: DataFrame containing the feature data

# Returns
- `feature_names`: Filtered vector with problematic features removed
"""
function remove_zero_variance_columns!(feature_names::Vector{Symbol}, df::AbstractDataFrame)
    removed_features = Symbol[]

    # Filter out problematic columns
    filter!(feature_names) do feature
        if !hasproperty(df, feature)
            push!(removed_features, feature)
            return false
        end

        col_data = df[!, feature]

        # Check for columns with all missing values
        if all(ismissing, col_data)
            push!(removed_features, feature)
            return false
        end

        # Check for Inf or NaN values
        if any(x -> !ismissing(x) && (isinf(x) || isnan(x)), col_data)
            push!(removed_features, feature)
            return false
        end

        # Check for zero variance (constant columns)
        non_missing_data = collect(skipmissing(col_data))
        if isempty(non_missing_data) || length(unique(non_missing_data)) <= 1 || var(non_missing_data) ≈ 0.0
            push!(removed_features, feature)
            return false
        end

        return true
    end

    # Log removed features if any
    if !isempty(removed_features)
        @user_warn "Removed $(length(removed_features)) problematic features to prevent singular matrix" removed_features = removed_features
    end

    return feature_names
end

"""
    filter_ms1_features_if_disabled!(feature_names::Vector{Symbol}, ms1_scoring::Bool)

Remove MS1 features from `feature_names` if MS1 scoring is disabled.
This prevents including MS1 features with zero variance when ms1_scoring=false.

# Arguments
- `feature_names`: Vector of feature symbols to filter
- `ms1_scoring`: Whether MS1 scoring is enabled

# Returns
- `feature_names`: Filtered vector with MS1 features removed if ms1_scoring=false
"""
function filter_ms1_features_if_disabled!(feature_names::Vector{Symbol}, ms1_scoring::Bool)
    if !ms1_scoring
        # Define MS1 features that should be excluded when ms1_scoring=false
        ms1_features = Set([
            :ms1_irt_diff, :weight_ms1, :gof_ms1, :max_matched_residual_ms1,
            :max_unmatched_residual_ms1, :fitted_spectral_contrast_ms1, :error_ms1,
            :m0_error_ms1, :n_iso_ms1, :big_iso_ms1, :rt_max_intensity_ms1,
            :rt_diff_max_intensity_ms1, :ms1_features_missing
        ])

        original_count = length(feature_names)
        filter!(feature -> !(feature in ms1_features), feature_names)
        removed_count = original_count - length(feature_names)

        if removed_count > 0
            @user_info "Excluded $removed_count MS1 features (ms1_scoring=false)"
        end
    end

    return feature_names
end





#==========================================================
PSM score merging and processing 
==========================================================#


"""
    get_pep_spline(merged_psms_path::String, score_col::Symbol;
                   min_pep_points_per_bin=5000, n_spline_bins=20) -> UniformSpline

Create posterior error probability spline from merged scores.

Returns spline for PEP calculation based on target/decoy distributions.
"""
function get_pep_spline(
                            merged_psms_path::String,
                            score_col::Symbol;
                            min_pep_points_per_bin = 5000,
                            n_spline_bins = 20
)

    psms_scores = Arrow.Table(merged_psms_path)
    Q = length(psms_scores[1])
    M = ceil(Int, Q / min_pep_points_per_bin)
    bin_target_fraction, bin_mean_prob = Vector{Float32}(undef, M), Vector{Float32}(undef, M)
    bin_size = 0
    bin_idx = 0
    mean_prob, targets = 0.0f0, 0
    for i in range(1, Q)
        bin_size += 1
        targets += psms_scores[:target][i]
        mean_prob += psms_scores[score_col][i]
        if bin_size == min_pep_points_per_bin
            bin_idx += 1
            bin_target_fraction[bin_idx] = targets/bin_size
            bin_mean_prob[bin_idx] = mean_prob/bin_size
            bin_size, targets, mean_prob = zero(Int64), zero(Int64), zero(Float32)
        end
    end
    bin_target_fraction[end] = targets/max(bin_size, 1)
    bin_mean_prob[end] = mean_prob/max(bin_size, 1)
    try 
        if length(bin_target_fraction)<20
            @user_warn "Less than 20 bins to estimate PEP. PEP results suspect..."
        end
        return UniformSpline(bin_target_fraction, bin_mean_prob, 3, 3)
    catch
        @user_warn "Failed to estimate PEP spline"
        return UniformSpline(SVector{4, Float32}([0, 0, 0, 0]), 3, 0.0f0, 1.0f0, 100.0f0)
    end
end

"""
    get_pep_interpolation(merged_psms_path::String, score_col::Symbol;
                          fdr_scale_factor=1.0f0)

Create an interpolation function mapping scores to posterior error probabilities.

PEP values are estimated across all runs using `get_PEP!` and then fit with a
linear interpolation.  The returned object can be used with
`add_interpolated_column` to populate per-file PEP columns.
"""
function get_pep_interpolation(
    merged_psms_path::String,
    score_col::Symbol;
    fdr_scale_factor::Float32 = 1.0f0,
)
    df = DataFrame(Arrow.Table(merged_psms_path))
    scores  = Vector{Float32}(df[!, score_col])
    targets = Vector{Bool}(df[!, :target])

    pep_vals = Vector{Float32}(undef, length(scores))
    get_PEP!(scores, targets, pep_vals; doSort=true,
             fdr_scale_factor=fdr_scale_factor)

    # Collapse redundant PEP values by taking the minimum score for each PEP
    pep_to_score = Dict{Float32, Float32}()
    order = sortperm(scores)
    for idx in order
        p = pep_vals[idx]
        s = scores[idx]
        if !haskey(pep_to_score, p)
            pep_to_score[p] = s
        end
    end

    # Sort the resulting pairs by score for interpolation
    ordered = sort(collect(pep_to_score), by = last)
    xs_tmp = Float32[last(pair) for pair in ordered]
    ys_tmp = Float32[first(pair) for pair in ordered]

    # Remove duplicate knot values explicitly
    xs = Float32[]
    ys = Float32[]
    prev_x = NaN32
    for (x, y) in zip(xs_tmp, ys_tmp)
        if x != prev_x && !isnan(x)
            push!(xs, x)
            push!(ys, clamp(y, 0.0f0, 1.0f0))
            prev_x = x
        end
    end

    if length(xs) < 2
        @user_warn "Insufficient unique points for PEP interpolation, using default"
        xs = Float32[0.0, 1.0]
        ys = Float32[1.0, 0.0]
    end

    Interpolations.deduplicate_knots!(xs)

    return linear_interpolation(xs, ys; extrapolation_bc=Interpolations.Flat())
end

"""
    get_qvalue_spline(merged_psms_path::String, score_col::Symbol;
                      min_pep_points_per_bin=1000, fdr_scale_factor=1.0f0) -> Interpolation

Create q-value interpolation function from merged scores.

Returns interpolation for mapping scores to q-values.
"""
function get_qvalue_spline(
                            merged_psms_path::String,
                            score_col::Symbol,
                            use_unique::Bool;
                            min_pep_points_per_bin = 1000,
                            fdr_scale_factor::Float32 = 1.0f0
)


    psms_scores = DataFrame(Arrow.Table(merged_psms_path))
    if use_unique
        # select the columns needed to identify globally unique precursors or proteins
        if score_col == :global_prob # precursors
            select!(psms_scores, [:precursor_idx, :target, score_col])
        elseif score_col == :global_pg_score # proteins
            select!(psms_scores, [:protein_name, :target, :entrap_id, score_col])
        end
        psms_scores = unique(psms_scores)
    end

    Q = size(psms_scores, 1)
    M = ceil(Int, (Q - 1) / min_pep_points_per_bin) + 1
    bin_qval, bin_mean_prob = Vector{Float32}(undef, M), Vector{Float32}(undef, M)
    bin_size = 0
    bin_idx = 0
    mean_prob, targets, decoys = 0.0f0, 0, 0
    for i in range(1, Q)
        targets += psms_scores[!, :target][i]
        decoys += (1 - psms_scores[!, :target][i])
    end
    min_q_val = typemax(Float32)
    for i in reverse(range(1, Q))
        bin_size += 1
        targets -= psms_scores[!, :target][i]
        decoys -= (1 - psms_scores[!, :target][i])
        mean_prob += psms_scores[!, score_col][i]
        if bin_idx == 0 || bin_size == min_pep_points_per_bin
            bin_idx += 1
            # Apply FDR scale factor to correct for library target/decoy ratio
            qval = (decoys * fdr_scale_factor) / max(targets, 1)
            if qval > min_q_val
                bin_qval[bin_idx] = min_q_val
            else
                min_q_val = qval
                bin_qval[bin_idx] = qval
            end
            bin_mean_prob[bin_idx] = mean_prob/bin_size
            bin_size, mean_prob = zero(Int64), zero(Float32)
        end
    end
    # Apply FDR scale factor to final bin calculation
    bin_qval[end] = (decoys * fdr_scale_factor) / max(targets, 1)
    bin_mean_prob[end] = mean_prob/bin_size
    prepend!(bin_qval, 1.0f0)
    prepend!(bin_mean_prob, 0.0f0)
    bin_qval = bin_qval[isnan.(bin_mean_prob).==false]
    bin_mean_prob = bin_mean_prob[isnan.(bin_mean_prob).==false]
    
    # Sort and deduplicate knot vectors
    # First, create paired array for sorting
    paired = collect(zip(bin_mean_prob, bin_qval))
    # Sort by probability (x-values)
    sort!(paired, by = x -> x[1])
    
    # Manual deduplication keeping the first occurrence of each x value
    xs = Float32[]
    ys = Float32[]
    prev_x = NaN32
    for (x, y) in paired
        if x != prev_x && !isnan(x)
            push!(xs, x)
            push!(ys, y)
            prev_x = x
        end
    end
    
    # Ensure we have at least 2 points for interpolation
    if length(xs) < 2
        @user_warn "Insufficient unique points for q-value interpolation, using default"
        xs = Float32[0.0, 1.0]
        ys = Float32[1.0, 0.0]
    end
    
    # Final check for sorted and unique
    if length(xs) > 1
        for i in 2:length(xs)
            if xs[i] <= xs[i-1]
                error("Failed to properly deduplicate knot vectors. Debug info: xs=$xs")
            end
        end
    end
    #=
    xs = copy(bin_mean_prob)
    ys = copy(bin_qval)
    # remove duplicates and keep the **first** occurrence
    mask = [i == 1 || xs[i] != xs[i-1] for i in eachindex(xs)]
    xs   = xs[mask]
    ys   = ys[mask]
    =#
    return linear_interpolation(xs, ys; extrapolation_bc = Line())
end
#==========================================================
Protein group analysis
==========================================================#
#=
"""
    getProteinGroupsDict(protein_inference_dict, psm_precursor_idx, psm_score, 
                        psm_is_target, psm_entrapment_id, precursors; min_peptides=2)

Create protein groups from PSMs and calculate group scores.

# Arguments
- `protein_inference_dict`: Maps peptides to inferred protein groups
- `psm_precursor_idx`: Precursor indices from PSMs
- `psm_score`: PSM probability scores
- `psm_is_target`: Boolean array indicating targets
- `psm_entrapment_id`: Entrapment group IDs
- `precursors`: Library precursor information
- `min_peptides`: Minimum peptides required per group

# Returns
- `pg_score`: Protein group scores for each PSM
- `inferred_protein_group_names`: Protein names for each PSM
- `protein_groups`: Dictionary of protein groups with scores and peptide sets
"""
function getProteinGroupsDict(
    protein_inference_dict::Dictionary{NamedTuple{(:peptide, :decoy, :entrap_id), Tuple{String, Bool, UInt8}}, NamedTuple{(:protein_name, :decoy, :entrap_id, :retain), Tuple{String, Bool, UInt8, Bool}}},
    psm_precursor_idx::AbstractVector{UInt32},
    psm_score::AbstractVector{Float32},
    psm_is_target::AbstractVector{Bool},
    psm_entrapment_id::AbstractVector{UInt8},
    precursors::LibraryPrecursors;
    min_peptides::Int64 = 1)

    precursor_sequence = getSequence(precursors)
    
    # Use new builder pattern
    group_builders = Dictionary{ProteinKey, ProteinGroupBuilder}()

    for i in range(1, length(psm_precursor_idx))
        precursor_idx = psm_precursor_idx[i]
        sequence = precursor_sequence[precursor_idx]
        
        # Create key for protein_inference_dict lookup
        peptide_key = (peptide = sequence, decoy = !psm_is_target[i], entrap_id = psm_entrapment_id[i])
        
        # Check if this peptide exists in our protein inference dictionary
        if !haskey(protein_inference_dict, peptide_key)
            throw("Peptide key not found error!")
            continue
        end
        
        # Skip if not retained for quantification
        if protein_inference_dict[peptide_key][:retain] == false
            continue
        end
        
        score = psm_score[i]
        protein_name = protein_inference_dict[peptide_key][:protein_name]
        protein_key = ProteinKey(protein_name, psm_is_target[i], psm_entrapment_id[i])
        
        # Get or create builder
        if !haskey(group_builders, protein_key)
            insert!(group_builders, protein_key, ProteinGroupBuilder(protein_key))
        end
        
        # Add peptide to group
        add_peptide!(group_builders[protein_key], sequence, score)
    end
    
    # Filter by minimum peptides
    filter!(x -> length(x.peptides) >= min_peptides, group_builders)

    # Convert back to old format for compatibility
    protein_groups = Dictionary{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8},
                               @NamedTuple{pg_score::Float32, peptides::Set{String}}}()
    
    for (key, builder) in pairs(group_builders)
        old_key = to_namedtuple(key)
        insert!(protein_groups, old_key, (pg_score = -builder.score, peptides = builder.peptides))
    end
    
    return protein_groups
end
=#

"""
    calculate_protein_features(builder::ProteinGroupBuilder, catalog::Dictionary{ProteinKey, Set{String}}) -> ProteinFeatures

Calculate statistical features for a protein group.

# Arguments
- `builder`: ProteinGroupBuilder with observed peptides
- `catalog`: Dictionary mapping protein keys to all possible peptides

# Returns
- `ProteinFeatures`: Statistical features for the protein group
"""
function calculate_protein_features(builder::ProteinGroupBuilder, catalog::Dictionary{ProteinKey, Set{String}})::ProteinFeatures
    n_peptides = UInt16(length(builder.peptides))
    
    # Get possible peptides for this protein
    possible_peptides = get(catalog, builder.key, Set{String}())
    n_possible_peptides = UInt16(length(possible_peptides))
    
    # Calculate coverage
    peptide_coverage = n_possible_peptides > 0 ? Float32(n_peptides) / Float32(n_possible_peptides) : 0.0f0
    
    # Calculate total peptide length
    total_peptide_length = UInt16(sum(length(peptide) for peptide in builder.peptides))
    
    # Calculate log features for regression
    log_n_possible_peptides = n_possible_peptides > 0 ? log(Float32(n_possible_peptides)) : 0.0f0
    
    # Calculate log binomial coefficient: log(n_possible choose n_peptides)
    log_binom_coeff = if n_possible_peptides > 0 && n_peptides > 0 && n_peptides <= n_possible_peptides
        # Use the more numerically stable logbinom function if available, otherwise approximate
        try
            # Simple approximation for log binomial coefficient
            lgamma(n_possible_peptides + 1) - lgamma(n_peptides + 1) - lgamma(n_possible_peptides - n_peptides + 1)
        catch
            0.0f0
        end
    else
        0.0f0
    end
    
    return ProteinFeatures(
        n_peptides,
        n_possible_peptides,
        peptide_coverage,
        total_peptide_length,
        log_n_possible_peptides,
        Float32(log_binom_coeff)
    )
end

"""
    write_protein_groups_arrow(protein_groups::Dictionary{ProteinKey, ProteinGroup}, output_path::String)

Write protein groups to Arrow file.

# Arguments
- `protein_groups`: Dictionary of protein groups
- `output_path`: Path to output Arrow file
"""
function write_protein_groups_arrow(protein_groups::Dictionary{ProteinKey, ProteinGroup}, output_path::String)
    if isempty(protein_groups)
        @user_warn "No protein groups to write to $output_path"
        return
    end
    
    # Extract data from protein groups
    n_groups = length(protein_groups)
    
    # Initialize arrays
    protein_names = Vector{String}(undef, n_groups)
    targets = Vector{Bool}(undef, n_groups)
    entrap_ids = Vector{UInt8}(undef, n_groups)
    pg_scores = Vector{Float32}(undef, n_groups)
    n_peptides = Vector{UInt16}(undef, n_groups)
    n_possible_peptides = Vector{UInt16}(undef, n_groups)
    peptide_coverages = Vector{Float32}(undef, n_groups)
    total_peptide_lengths = Vector{UInt16}(undef, n_groups)
    log_n_possible_peptides = Vector{Float32}(undef, n_groups)
    log_binom_coeffs = Vector{Float32}(undef, n_groups)
    
    # Fill arrays
    i = 1
    for (key, group) in pairs(protein_groups)
        protein_names[i] = key.name
        targets[i] = key.is_target
        entrap_ids[i] = key.entrap_id
        pg_scores[i] = group.score
        n_peptides[i] = group.features.n_peptides
        n_possible_peptides[i] = group.features.n_possible_peptides
        peptide_coverages[i] = group.features.peptide_coverage
        total_peptide_lengths[i] = group.features.total_peptide_length
        log_n_possible_peptides[i] = group.features.log_n_possible_peptides
        log_binom_coeffs[i] = group.features.log_binom_coeff
        i += 1
    end
    
    # Create DataFrame
    df = DataFrame(
        protein_name = protein_names,
        target = targets,
        entrap_id = entrap_ids,
        pg_score = pg_scores,
        n_peptides = n_peptides,
        n_possible_peptides = n_possible_peptides,
        peptide_coverage = peptide_coverages,
        total_peptide_length = total_peptide_lengths,
        log_n_possible_peptides = log_n_possible_peptides,
        log_binom_coeff = log_binom_coeffs
    )
    
    # Sort by pg_score and target in descending order
    sort!(df, [:pg_score, :target], rev = [true, true])
    
    # Write to Arrow file
    writeArrow(output_path, df)
end

"""
    writeProteinGroups(acc_to_max_pg_score, protein_groups, 
                        protein_to_possible_peptides, protein_groups_path)

Write protein groups with features to Arrow file.

# Arguments
- `acc_to_max_pg_score`: Maximum scores across runs for each protein
- `protein_groups`: Dictionary of protein groups with scores and peptides
- `protein_to_possible_peptides`: All possible peptides for each protein
- `protein_groups_path`: Output file path

# Returns
- Number of protein groups written

# Output columns
- Basic: protein_name, target, entrap_id, pg_score, global_pg_score
- Features: n_peptides, total_peptide_length, n_possible_peptides, peptide_coverage
"""
function writeProteinGroups(
                                protein_groups::Dictionary{
                                    @NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8},
                                    @NamedTuple{pg_score::Float32,  peptides::Set{String}}
                                },
                                protein_accession_to_possible_peptides::Dict{
                                    @NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8},
                                    Set{String}
                                },
                                protein_groups_path::String)
    # Convert to new type system
    catalog = Dictionary{ProteinKey, Set{String}}()
    for (old_key, peptide_set) in protein_accession_to_possible_peptides
        new_key = ProteinKey(old_key.protein_name, old_key.target, old_key.entrap_id)
        insert!(catalog, new_key, peptide_set)
    end
    
    # Build ProteinGroup objects with features
    full_protein_groups = Dictionary{ProteinKey, ProteinGroup}()
    
    for (old_key, old_value) in pairs(protein_groups)
        key = ProteinKey(old_key.protein_name, old_key.target, old_key.entrap_id)
        
        # Create a builder and calculate features
        builder = ProteinGroupBuilder(key)
        builder.score = -old_value.pg_score  # Convert back to log space temporarily
        builder.peptides = old_value.peptides
        
        features = calculate_protein_features(builder, catalog)
        protein_group = finalize(builder, features)
        
        insert!(full_protein_groups, key, protein_group)
    end
    
    # Write using the new function
    write_protein_groups_arrow(full_protein_groups, protein_groups_path)
end

"""
    count_protein_peptides(precursors::LibraryPrecursors)

Count all possible peptides for each protein in the library.

# Arguments
- `precursors`: Library precursors

# Returns
- Dictionary mapping protein keys to sets of peptide sequences
"""
function count_protein_peptides(precursors::LibraryPrecursors)
    peptide_count_start = time()
    protein_to_possible_peptides = Dict{@NamedTuple{protein_name::String, target::Bool, entrap_id::UInt8}, Set{String}}()
    
    # Count all peptides in the library for each protein
    all_accession_numbers = getAccessionNumbers(precursors)
    all_sequences = getSequence(precursors)
    all_decoys = getIsDecoy(precursors)
    all_entrap_ids = getEntrapmentGroupId(precursors)
    n_precursors = length(all_accession_numbers)

    for i in 1:n_precursors
        protein_names = split(all_accession_numbers[i], ';')  # Handle shared peptides
        is_decoy = all_decoys[i]
        entrap_id = all_entrap_ids[i]
        
        for protein_name in protein_names
            key = (protein_name = String(protein_name), target = !is_decoy, entrap_id = entrap_id)
            if !haskey(protein_to_possible_peptides, key)
                protein_to_possible_peptides[key] = Set{String}()
            end
            push!(protein_to_possible_peptides[key], all_sequences[i])
        end
    end

    return protein_to_possible_peptides
end


"""
    perform_protein_probit_regression(pg_refs::Vector{ProteinGroupFileReference},
                                    max_psms_in_memory::Int64,
                                    qc_folder::String,
                                    precursors::LibraryPrecursors;
                                    protein_to_cv_fold::Union{Nothing, Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}} = nothing,
                                    ms1_scoring::Bool = true)

Perform probit regression on protein groups.

# Arguments
- `pg_refs`: Vector of protein group file references
- `max_psms_in_memory`: Memory limit for in-memory vs OOM processing
- `qc_folder`: Folder for QC plots
- `precursors`: Library precursors
- `protein_to_cv_fold`: Optional pre-built mapping of proteins to CV folds
- `ms1_scoring`: Whether MS1 scoring is enabled (affects feature selection)
"""
function perform_protein_probit_regression(
    pg_refs::Vector{ProteinGroupFileReference},
    max_psms_in_memory::Int64,
    qc_folder::String,
    precursors::LibraryPrecursors;
    protein_to_cv_fold::Union{Nothing, Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}} = nothing,
    ms1_scoring::Bool = true
)
    # Extract paths for compatibility with existing code
    passing_pg_paths = [file_path(ref) for ref in pg_refs]
    
    # First, count the total number of protein groups across all files
    total_protein_groups = 0
    for ref in pg_refs
        if exists(ref)
            total_protein_groups += row_count(ref)
        end
    end
    
    # Set protein group limit to 5x the precursor limit
    max_protein_groups_in_memory_limit = 5 * max_psms_in_memory

    # HARDCODED: Always use in-memory processing (OOM path disabled)
    if false  # total_protein_groups > max_protein_groups_in_memory_limit
        #Need to implement safety checks for minimal number of targets/decoys in each split (DISABLED) 
        
        # Check if we should skip scoring in OOM path
        # We need to load a sample to check targets/decoys
        sample_df = DataFrame()
        for ref in pg_refs[1:min(3, length(pg_refs))]  # Sample first few files
            if exists(ref)
                append!(sample_df, DataFrame(Tables.columntable(Arrow.Table(file_path(ref)))))
            end
        end
        n_sample_targets = sum(sample_df.target)
        n_sample_decoys = sum(.!sample_df.target)
        
        # Estimate total based on sample
        est_targets = n_sample_targets * (total_protein_groups / nrow(sample_df))
        est_decoys = n_sample_decoys * (total_protein_groups / nrow(sample_df))
        
        skip_scoring_oom = !(est_targets > 10 && est_decoys > 10 && total_protein_groups > 1000)
        
        perform_probit_analysis_oom(pg_refs, total_protein_groups, max_protein_groups_in_memory_limit, qc_folder;
                                   skip_scoring = skip_scoring_oom, ms1_scoring = ms1_scoring)
    else
        # Load all protein group tables into a single DataFrame
        all_protein_groups = DataFrame()
        for pg_path in passing_pg_paths
            if isfile(pg_path) && endswith(pg_path, ".arrow")
                append!(all_protein_groups, DataFrame(Tables.columntable(Arrow.Table(pg_path))))
            end
        end
        
        # Only perform analysis if we have both targets and decoys
        n_targets = sum(all_protein_groups.target)
        n_decoys = sum(.!all_protein_groups.target)
        
        # Always run through the same probit workflow, just skip scoring if insufficient data
        skip_scoring = !(n_targets > 50 && n_decoys > 50 && nrow(all_protein_groups) > 1000)
        
        # Always use the same function, just with different skip_scoring parameter
        perform_probit_analysis_multifold(
            all_protein_groups,
            qc_folder,
            pg_refs,
            precursors;
            protein_to_cv_fold = protein_to_cv_fold,
            skip_scoring = skip_scoring,
            ms1_scoring = ms1_scoring
        )
    end
end

"""
    update_psms_with_probit_scores_refs(paired_refs::Vector{PairedSearchFiles},
                                       pg_name_to_global_pg_score::Dict{ProteinKey,Float32},
                                       pg_score_to_qval::Interpolations.Extrapolation,
                                       global_pg_score_to_qval_dict::Dict{Tuple{String,Bool,UInt8}, Float32})

Update PSMs with probit-scored pg_score values and q-values using references.

# Arguments
- `paired_refs`: Paired PSM/protein group file references
- `pg_name_to_global_pg_score`: Dictionary mapping protein keys to global scores
- `pg_score_to_qval`: Interpolation function for pg_score to q-value
- `global_pg_score_to_qval_dict`: Dictionary mapping (protein_name, target, entrap_id) to global q-value
"""
function update_psms_with_probit_scores_refs(
    paired_refs::Vector{PairedSearchFiles},
    pg_name_to_global_pg_score::Dict{ProteinKey,Float32},
    pg_score_to_qval::Interpolations.Extrapolation,
    global_pg_score_to_qval_dict::Dict{Tuple{String,Bool,UInt8}, Float32}
)

    total_psms_updated = 0
    files_processed = 0
    
    for paired_ref in paired_refs
        psm_ref = paired_ref.psm_ref
        pg_ref = paired_ref.protein_ref
        
        # Verify both files exist
        if !exists(psm_ref) || !exists(pg_ref)
            @user_warn "Missing file" psm_path=file_path(paired_ref.psm_ref) pg_path=file_path(paired_ref.protein_ref)
            continue
        end
        
        # Load protein groups to get probit-scored pg_score values
        pg_table = Arrow.Table(file_path(pg_ref))
        
        # Create lookup dictionary: ProteinKey -> (probit pg_score, PEP)
        pg_score_lookup = Dict{ProteinKey, Tuple{Float32, Float32}}()
        n_pg_rows = length(pg_table[:protein_name])
        
        for i in 1:n_pg_rows
            key = ProteinKey(
                pg_table[:protein_name][i],
                pg_table[:target][i],
                pg_table[:entrap_id][i]
            )
            pep_val = hasproperty(pg_table, :pg_pep) ? pg_table[:pg_pep][i] : 1.0f0
            pg_score_lookup[key] = (pg_table[:pg_score][i], pep_val) # pg_score is now the probit score
        end
        
        # Transform PSM file
        transform_and_write!(psm_ref) do psms_df
            n_psms = nrow(psms_df)
            n_missing_pg = 0
            n_not_for_quant = 0
            
            # Update pg_score column with probit scores
            probit_pg_scores = Vector{Union{Missing, Float32}}(undef, n_psms)
            global_pg_scores = Vector{Union{Missing, Float32}}(undef, n_psms)
            pg_qvals = Vector{Union{Missing, Float32}}(undef, n_psms)
            global_pg_qvals = Vector{Union{Missing, Float32}}(undef, n_psms)
            pg_peps = Vector{Union{Missing, Float32}}(undef, n_psms)
            
            for i in 1:n_psms
                # Skip if missing inferred protein group (shared peptides excluded from inference)
                if ismissing(psms_df[i, :inferred_protein_group])
                    n_missing_pg += 1
                    probit_pg_scores[i] = missing
                    global_pg_scores[i] = missing
                    pg_qvals[i] = missing
                    global_pg_qvals[i] = missing
                    pg_peps[i] = missing
                    continue
                end

                # Skip if peptide didn't match to a distinct protein group
                if psms_df[i,:use_for_protein_quant] == false
                    n_not_for_quant += 1
                    #Should be able to make this 'missing' since that is more clear
                    probit_pg_scores[i] = missing
                    global_pg_scores[i] =  missing
                    pg_qvals[i] =  missing
                    global_pg_qvals[i] =  missing
                    pg_peps[i] =  missing
                    continue
                end
                
                # Create key for lookup
                key = ProteinKey(
                    psms_df[i, :inferred_protein_group],
                    psms_df[i, :target],
                    psms_df[i, :entrapment_group_id]
                )
                
                # Get scores and PEP
                if !haskey(pg_score_lookup, key)
                    n_missing_pg += 1
                    # Add detailed logging to understand why keys don't match
                    #if n_missing_pg <= 5  # Only log first few to avoid spam
                    #    @user_warn "Protein group not found in lookup" missing_key=key psm_idx=i available_keys_sample=collect(keys(pg_score_lookup))[1:min(5, length(pg_score_lookup))]
                    #end
                    #Should be able to make this 'missing' since that is more clear 
                    probit_pg_scores[i] = missing
                    global_pg_scores[i] = missing
                    pg_qvals[i] =  missing
                    global_pg_qvals[i] =  missing
                    pg_peps[i] = missing
                    continue
                end
                scores_tuple = pg_score_lookup[key]
                probit_pg_scores[i] = scores_tuple[1]
                pg_peps[i] = scores_tuple[2]
                
                if !haskey(pg_name_to_global_pg_score, key)
                    throw("Missing global pg score lookup key!!!")
                end
                global_pg_scores[i] = pg_name_to_global_pg_score[key]

                # Calculate q-values
                pg_qvals[i] = pg_score_to_qval(probit_pg_scores[i])
                # Look up global q-value from dictionary using protein group key
                dict_key = (key.name, key.is_target, key.entrap_id)
                global_pg_qvals[i] = get(global_pg_score_to_qval_dict, dict_key, missing)
            end
            
            # Update columns
            psms_df[!, :pg_score] = probit_pg_scores
            psms_df[!, :global_pg_score] = global_pg_scores
            psms_df[!, :pg_qval] = pg_qvals
            psms_df[!, :qlobal_pg_qval] = global_pg_qvals
            psms_df[!, :pg_pep] = pg_peps
            
            total_psms_updated += n_psms
            
            #if n_missing_pg > 0 || n_not_for_quant > 0
            #    @user_warn "PSMs with missing protein scores" file=file_path(psm_ref) total_psms=n_psms missing_pg=n_missing_pg not_for_quant=n_not_for_quant psms_with_scores=(n_psms - n_missing_pg - n_not_for_quant)
            #end
            
            return psms_df
        end
        
        files_processed += 1
    end
end




# DISABLED: OOM protein probit function - only used when total_protein_groups > max_protein_groups_in_memory_limit
# This function is commented out in favor of always using in-memory processing.
# Preserved for potential future use if needed for extremely large datasets.
#=
"""
    perform_probit_analysis_oom(pg_paths::Vector{String}, total_protein_groups::Int, 
                               max_protein_groups_in_memory::Int, qc_folder::String)

Perform out-of-memory probit regression analysis on protein groups.

# Arguments
- `pg_paths`: Vector of paths to protein group Arrow files
- `total_protein_groups`: Total number of protein groups across all files
- `max_protein_groups_in_memory`: Maximum number of protein groups to hold in memory
- `qc_folder`: Folder for QC plots

# Process
1. Calculate sampling ratio for each file
2. Sample protein groups from each file proportionally
3. Fit probit model on sampled data
4. Apply model to all protein groups file by file
5. Calculate and report performance metrics
"""
function perform_probit_analysis_oom(pg_refs::Vector{ProteinGroupFileReference}, total_protein_groups::Int,
                                    max_protein_groups_in_memory::Int, qc_folder::String;
                                    skip_scoring = false, ms1_scoring::Bool = true)
    
    # Calculate sampling ratio
    sampling_ratio = max_protein_groups_in_memory / total_protein_groups
    sampled_targets, sampled_decoys = 0, 0
    # Sample protein groups from each file
    sampled_protein_groups = DataFrame()
    for ref in pg_refs
        if exists(ref)
            n_rows = row_count(ref)
            n_sample = ceil(Int, n_rows * sampling_ratio)
            
            if n_sample > 0
                # Sample indices
                sample_indices = sort(sample(1:n_rows, n_sample, replace=false))
                
                # Load sampled rows
                table = Arrow.Table(file_path(ref))
                df = DataFrame(Tables.columntable(table))
                # Count statistics
                sampled_df = df[sample_indices, :]
                _df_targets = sum(sampled_df.target)
                _df_decoys = size(sampled_df, 1) - _df_targets
                sampled_targets += _df_targets
                sampled_decoys += _df_decoys
                append!(sampled_protein_groups, sampled_df)
            end
        end
    end
    
    n_sampled = sampled_targets + sampled_decoys
    if (n_sampled < 1000) || (sampled_targets < 50) || (sampled_decoys < 50)
        @user_warn "Insufficient sampled protein groups for OOM probit regression: targets=$sampled_targets, decoys=$sampled_decoys"
        return
    end
    
    # Define features to use
    feature_names = [:pg_score, :peptide_coverage, :n_possible_peptides, :log_binom_coeff, :any_common_peps]

    # Apply feature filtering
    adjust_any_common_peps!(feature_names, sampled_protein_groups)
    remove_zero_variance_columns!(feature_names, sampled_protein_groups)

    if isempty(feature_names)
        @user_warn "No valid features remaining for OOM probit regression after filtering"
        return
    end

    X = Matrix{Float64}(sampled_protein_groups[:, feature_names])
    y = sampled_protein_groups.target
    
    # Fit probit model on sampled data (skip if skip_scoring = true)
    if !skip_scoring
        β_fitted = fit_probit_model(X, y)
    else
        β_fitted = Float64[]  # Empty model when skipping
    end
    
    total_targets, total_decoys = 0, 0
    # Process each file
    for ref in pg_refs
        if exists(ref)
            # Load file
            df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
            
            if !skip_scoring
                # Calculate probit scores
                X_file = Matrix{Float64}(df[:, feature_names])
                prob_scores = calculate_probit_scores(X_file, β_fitted)
                
                # Overwrite pg_score with probit scores
                df[!, :pg_score] = Float32.(prob_scores)
            end
            # If skip_scoring, keep original pg_score values
            
            # Sort by pg_score and target in descending order
            sort!(df, [:pg_score, :target], rev = [true, true])
            
            total_targets += sum(df.target)
            total_decoys += sum(.!df.target)
            # Write back the file
            writeArrow(file_path(ref), df)
        end
    end
end
=#

"""
    perform_probit_analysis(all_protein_groups::DataFrame, qc_folder::String, 
                          pg_paths::Union{Vector{String}, Nothing} = nothing)

Perform probit regression analysis on protein groups with comparison to baseline.

# Arguments
- `all_protein_groups::DataFrame`: Protein group data with features
- `qc_folder::String`: Folder for QC plots
- `pg_paths::Union{Vector{String}, Nothing}`: Paths to protein group files for re-processing

# Process
1. Fits probit model with multiple features
2. Calculates performance metrics
3. Compares to pg_score-only baseline
4. Generates decision boundary plots
5. Re-processes individual files with probit scores if pg_paths provided
"""
function perform_probit_analysis(all_protein_groups::DataFrame, qc_folder::String,
                               pg_refs::Vector{ProteinGroupFileReference};
                               show_improvement = true)
    n_targets = sum(all_protein_groups.target)
    n_decoys = sum(.!all_protein_groups.target)
    # Define features to use
    feature_names = [:pg_score, :peptide_coverage, :n_possible_peptides, :any_common_peps, :all_precursors_mbr] # :log_binom_coeff]

    # Apply feature filtering
    adjust_any_common_peps!(feature_names, all_protein_groups)
    remove_zero_variance_columns!(feature_names, all_protein_groups)

    if isempty(feature_names)
        @user_warn "No valid features remaining for probit regression after filtering"
        return
    end

    X = Matrix{Float64}(all_protein_groups[:, feature_names])
    y = all_protein_groups.target

    # Fit probit model
    #β_fitted, X_mean, X_std = fit_probit_model(X, y)
    β_fitted = fit_probit_model(X, y)
    # Re-process individual files if references are provided
    if !isempty(pg_refs)
        # Use the new apply_probit_scores! function with references
        apply_probit_scores!(pg_refs, β_fitted, feature_names)
    end
end


"""
    fit_probit_model(X::Matrix{Float64}, y::Vector{Bool})

Fit a probit regression model for protein group classification.

# Arguments
- `X::Matrix{Float64}`: Feature matrix
- `y::Vector{Bool}`: Target labels (true for targets, false for decoys)

# Returns
- `β_fitted`: Fitted coefficients
- `X_mean`: Feature means for standardization
- `X_std`: Feature standard deviations for standardization
"""
function fit_probit_model(X::Matrix{Float64}, y::Vector{Bool})
    # Check for problematic columns in the feature matrix
    X_df_temp = DataFrame(X, Symbol.("feature_", 1:size(X, 2)))
    feature_names = Symbol.("feature_", 1:size(X, 2))

    # Remove zero-variance columns to prevent singular matrices
    valid_features = remove_zero_variance_columns!(copy(feature_names), X_df_temp)

    if isempty(valid_features)
        throw(ArgumentError("No valid features remaining for probit regression after variance filtering"))
    end

    # Use only valid features
    valid_indices = [findfirst(==(f), feature_names) for f in valid_features]
    X_filtered = X[:, valid_indices]

    # Standardize features
    #X_mean = mean(X_filtered, dims=1)
    #X_std = std(X_filtered, dims=1)
    #X_std[X_std .== 0] .= 1.0  # Should not happen after filtering, but defensive
    X_standardized = X_filtered #(X_filtered .- X_mean) ./ X_std

    # Add intercept column
    X_with_intercept = hcat(ones(size(X_standardized, 1)), X_standardized)
    X_df = DataFrame(X_with_intercept, [:intercept; valid_features])

    # Initialize coefficients
    β = zeros(Float64, size(X_with_intercept, 2))

    # Create data chunks for parallel processing
    n_chunks = max(1, Threads.nthreads())
    chunk_size = max(1, ceil(Int, length(y) / n_chunks))
    data_chunks = Iterators.partition(1:length(y), chunk_size)

    # Fit probit model
    β_fitted = Pioneer.ProbitRegression(β, X_df, y, data_chunks, max_iter=30)

    return β_fitted#, vec(X_mean), vec(X_std)
end

"""
    calculate_probit_scores(X::Matrix{Float64}, β::Vector{Float64}, X_mean::Vector{Float64}, X_std::Vector{Float64})

Calculate probit probability scores for new data.

# Arguments
- `X::Matrix{Float64}`: Feature matrix
- `β::Vector{Float64}`: Fitted coefficients
- `X_mean::Vector{Float64}`: Feature means from training
- `X_std::Vector{Float64}`: Feature standard deviations from training

# Returns
- `Vector{Float64}`: Probability scores
"""
function calculate_probit_scores(X::Matrix{Float64}, β::Vector{Float64}
    #, X_mean::Vector{Float64}, X_std::Vector{Float64}
    )
    # Standardize using training statistics
    X_standardized = X#(X .- X_mean') ./ X_std'
    
    # Add intercept
    X_with_intercept = hcat(ones(size(X_standardized, 1)), X_standardized)
    X_df = DataFrame(X_with_intercept, [:intercept; Symbol.("feature_", 1:size(X, 2))])
    
    # Create data chunks
    n_chunks = max(1, Threads.nthreads())
    chunk_size = max(1, ceil(Int, size(X, 1) / n_chunks))
    data_chunks = Iterators.partition(1:size(X, 1), chunk_size)
    
    # Calculate probabilities
    prob_scores = zeros(Float64, size(X, 1))
    Pioneer.ModelPredictProbs!(prob_scores, X_df, β, data_chunks)
    
    return prob_scores
end

"""
    add_pep_column(new_col::Symbol, score_col::Symbol, target_col::Symbol;
                   doSort::Bool=false, fdr_scale_factor::Float32=1.0f0)

Add a column of posterior error probabilities using `get_PEP!`.
The input DataFrame must already be sorted by `score_col` if `doSort` is false.
"""
function add_pep_column(new_col::Symbol, score_col::Symbol, target_col::Symbol;
                        doSort::Bool=false, fdr_scale_factor::Float32=1.0f0)
    desc = "add_pep_column($new_col from $score_col)"
    op = function(df)
        scores = df[!, score_col]::AbstractVector{Float32}
        targets = df[!, target_col]::AbstractVector{Bool}
        pep_vals = Vector{Float32}(undef, length(scores))
        get_PEP!(scores, targets, pep_vals; doSort=doSort, fdr_scale_factor=fdr_scale_factor)
        df[!, new_col] = pep_vals
        return df
    end
    return desc => op
end




#=
#
#=
Have N of these tables. Need to combine into one sorted Arrow table without loading all tables
into memory at once.
julia> DataFrame(Arrow.Table(readdir(second_quant_folder, join = true)[1]))
280488×11 DataFrame
    Row │ precursor_idx  trace_prob  weight         target  irt_obs    missed_cleavage  isotopes_captured  scan_idx  ms_file_idx  peak_area   new_best_scan
        │ UInt32         Float32     Float32        Bool    Float32    UInt8            Tuple{Int8, Int8}  UInt32    Int64        Float32     UInt32
────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
      1 │       1468734    0.808784    201.467        true   0.659221                1  (0, 3)                49270            1     6.03541          50180
      2 │        262434    0.989585   2696.17         true   0.659221                0  (0, 3)                76753            1   121.201            76753
=#

getColNames(at::Arrow.Table) = keys(at)
getColTypes(at::Arrow.Table) = [eltype(at[col]) for col in getColNames(at)]
function getEmptyDF(at::Arrow.Table, N::Int)
    df = DataFrame()
    [df[!,Symbol(col)] = Vector{coltype}(undef, N) for (col, coltype) in zip(getColNames(at), getColTypes(at))]
    return df
end
    
function addPrecursorToHeap!(
    precursor_heap::BinaryMinHeap{Tuple{UInt32, UInt32, Int64}},
    first_sort_key::AbstractVector{UInt32},
    second_sort_key::AbstractVector{UInt32},
    table_idx::Int64,
    row_idx::Int64)
    push!(
        precursor_heap,
        (
        first_sort_key[row_idx],
        second_sort_key[row_idx],
        table_idx
        )
    )
end

function addPrecursorToHeap!(
    precursor_heap::BinaryMinHeap{Tuple{S, UInt32, Int64}},
    first_sort_key::AbstractVector{S},
    second_sort_key::AbstractVector{UInt32},
    table_idx::Int64,
    row_idx::Int64) where {S<:AbstractString}
    push!(
        precursor_heap,
        (
        first_sort_key[row_idx],
        second_sort_key[row_idx],
        table_idx
        )
    )
end

function mergeSortedArrowTables(
    input_dir::String, 
    output_path::String,
    sort_keys::Tuple{Symbol, Symbol};
    N = 1000000
)

    function fillColumn!(
        peptide_batch_col::Vector{R},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    ) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::R
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{Union{Missing, R}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    ) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::Union{Missing, R}
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{String},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    )
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::String
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{Union{Missing, String}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    )
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::Union{String,Missing}
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{Tuple{R, R}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    ) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::Tuple{R, R}
        end
    end

    #Get all .arrow files in the input 
    input_paths = [path for path in readdir(input_dir, join=true) if endswith(path, ".arrow")]
    #Keep track of which tables have 
    tables = [Arrow.Table(path) for path in input_paths]
    table_idxs = ones(Int64, length(tables))

    peptide_batch = getEmptyDF(first(tables), N)
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, N)
    peptide_batch_names = Symbol.(names(peptide_batch))
    precursor_heap = BinaryMinHeap{Tuple{eltype(first(tables)[first(sort_keys)]), UInt32, Int64}}()
    for (i, table) in enumerate(tables)
        addPrecursorToHeap!(
            precursor_heap,
            table[first(sort_keys)],
            table[last(sort_keys)],
            i,
            1
        )
    end
    i = 1
    n_writes = 0
    while length(precursor_heap) > 0
        _, _, table_idx = pop!( precursor_heap)
        table = tables[table_idx]
        idx = table_idxs[table_idx]
        sorted_tuples[i] = (table_idx, idx)
        table_idxs[table_idx] += 1
        idx = table_idxs[table_idx]
        if (idx > length(table[1]))
            continue
        end
        addPrecursorToHeap!(
            precursor_heap,
            table[first(sort_keys)],
            table[last(sort_keys)],
            table_idx,
            idx
        )
        i += 1
        if i > N
            for col in peptide_batch_names
                fillColumn!(
                    peptide_batch[!,col],
                    col,
                    sorted_tuples,
                    tables,
                    i
                )
            end
            if iszero(n_writes)
                if isfile(output_path)
                    rm(output_path)
                end
                open(output_path, "w") do io
                    Arrow.write(io, peptide_batch; file=false)  # file=false creates stream format
                end
            else
                Arrow.append(
                    output_path,
                    peptide_batch
                )
            end
            n_writes += 1
            i = 1
        end
    end
    for col in peptide_batch_names
        fillColumn!(
            peptide_batch[!,col],
            col,
            sorted_tuples,
            tables,
            i
        )
    end
    Arrow.append(
        output_path,
        peptide_batch[range(1, max(1, i - 1)),:]
    )
    return nothing
end

function addPrecursorToHeap!(
    precursor_heap::BinaryMinHeap{Tuple{Float32, Int64}},
    first_sort_key::AbstractVector{Float32},
    table_idx::Int64,
    row_idx::Int64)
    push!(
        precursor_heap,
        (
        first_sort_key[row_idx],
        table_idx
        )
    )
end

function addPrecursorToHeap!(
    precursor_heap::BinaryMinHeap{Tuple{Float32, Int64}},
    first_sort_key::AbstractVector{Union{Missing, Float32}},
    table_idx::Int64,
    row_idx::Int64)
    push!(
        precursor_heap,
        (
        coalesce(first_sort_key[row_idx], 0.0f0),
        table_idx
        )
    )
end

function mergeSortedArrowTables(
    input_dir::String, 
    output_path::String,
    sort_key::Symbol;
    N = 1000000)

    function fillColumn!(
        peptide_batch_col::Vector{R},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    ) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::R
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{Union{Missing, R}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    ) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::Union{Missing, R}
        end
    end

    function fillColumn!(
        peptide_batch_col::Vector{String},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    )
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::String
        end
    end


    function fillColumn!(
        peptide_batch_col::Vector{Union{Missing, String}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n
    )
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::String
        end
    end


    function fillColumn!(
        peptide_batch_col::Vector{Tuple{R, R}},
        col::Symbol,
        sorted_tuples::Vector{Tuple{Int64, Int64}},
        tables::Vector{Arrow.Table},
        n) where {R<:Real}
        for i in range(1, max(min(length(peptide_batch_col), n - 1), 1))
            table_idx, idx = sorted_tuples[i]
            peptide_batch_col[i] = tables[table_idx][col][idx]::Tuple{R, R}
        end
    end

    #Get all .arrow files in the input 
    input_paths = [path for path in readdir(input_dir, join=true) if endswith(path, ".arrow")]
    #Keep track of which tables have 
    tables = [Arrow.Table(path) for path in input_paths]
    table_idxs = ones(Int64, length(tables))

    peptide_batch = getEmptyDF(first(tables), N)
    sorted_tuples = Vector{Tuple{Int64, Int64}}(undef, N)
    peptide_batch_names = Symbol.(names(peptide_batch))
    precursor_heap = BinaryMinHeap{Tuple{Float32, Int64}}()
    for (i, table) in enumerate(tables)
        addPrecursorToHeap!(
            precursor_heap,
            table[sort_key],
            i,
            1
        )
    end
    i = 1
    n_writes = 0
    while length(precursor_heap) > 0
        _, table_idx = pop!(precursor_heap)
        table = tables[table_idx]
        idx = table_idxs[table_idx]
        sorted_tuples[i] = (table_idx, idx)
        table_idxs[table_idx] += 1
        idx = table_idxs[table_idx]
        if (idx > length(table[1]))
            continue
        end
        addPrecursorToHeap!(
            precursor_heap,
            table[sort_key],
            table_idx,
            idx
        )
        i += 1
        if i > N
            for col in peptide_batch_names
                fillColumn!(
                    peptide_batch[!,col],
                    col,
                    sorted_tuples,
                    tables,
                    i
                )
            end

            if iszero(n_writes)
                if isfile(output_path)
                    rm(output_path)
                end
                open(output_path, "w") do io
                    Arrow.write(io, peptide_batch; file=false)  # file=false creates stream format
                end
            else
                Arrow.append(
                    output_path,
                    peptide_batch
                )
            end
            n_writes += 1
            i = 1
        end
    end
    for col in peptide_batch_names
        fillColumn!(
            peptide_batch[!,col],
            col,
            sorted_tuples,
            tables,
            i
        )
    end
    Arrow.append(
        output_path,
        peptide_batch[range(1, max(1, i - 1)),:]
    )
    return nothing
end


=#

#==========================================================
Multi-fold Cross-Validation Probit Analysis
==========================================================#

"""
    detect_unique_cv_folds(precursors::LibraryPrecursors) -> Vector{UInt8}

Detect unique CV fold values from the library precursors.

# Returns
- Sorted vector of unique CV fold values
"""
function detect_unique_cv_folds(precursors::LibraryPrecursors)
    return sort(unique(precursors.pid_to_cv_fold))
end

"""
    get_corresponding_psm_path(pg_ref::ProteinGroupFileReference) -> String

Map a protein group file reference to its corresponding PSM file path.

# Arguments
- `pg_ref`: Protein group file reference

# Returns
- Path to the corresponding PSM file
"""
function get_corresponding_psm_path(pg_ref::ProteinGroupFileReference)
    pg_path = file_path(pg_ref)
    # Replace "passing_proteins" with "scored_PSMs" in path
    return replace(pg_path, "passing_proteins" => "scored_PSMs")
end

"""
    build_protein_cv_fold_mapping(psm_paths::Vector{String}, precursors::LibraryPrecursors)
    -> Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}

Build a mapping from protein names to their CV fold assignments based on highest-scoring peptides.

# Arguments
- `psm_paths`: Vector of paths to PSM files
- `precursors`: Library precursors containing cv_fold information

# Returns
- Dictionary mapping protein_name to named tuple with best_score and cv_fold

# Process
1. Scans PSM files to find peptides for each protein
2. Determines cv_fold of highest-scoring peptide per protein
3. Returns the protein_to_cv_fold mapping
"""
function build_protein_cv_fold_mapping(
    psm_paths::Vector{String},
    precursors::LibraryPrecursors
)
    @user_info "Building protein to CV fold mapping from PSM files"
    # Create mapping: protein_name -> (best_score=score, cv_fold=fold)
    protein_to_cv_fold = Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}()

    # Process each PSM file
    for psm_path in psm_paths
        # Skip if PSM file doesn't exist
        if !isfile(psm_path)
            @user_warn "PSM file not found: $psm_path"
            continue
        end

        # Load PSM data
        psms = DataFrame(Arrow.Table(psm_path))

        # Verify required columns exist
        required_columns = [:inferred_protein_group, :prec_prob, :precursor_idx]
        missing_columns = [col for col in required_columns if !hasproperty(psms, col)]
        if !isempty(missing_columns)
            error("PSM file $psm_path is missing required columns: $missing_columns")
        end

        # Filter for valid PSMs with protein group assignments
        psms = filter(row -> !ismissing(row.inferred_protein_group), psms)

        # Skip if no valid PSMs
        if nrow(psms) == 0
            continue
        end

        # Group by inferred_protein_group
        for group in groupby(psms, :inferred_protein_group)
            protein_name = first(group.inferred_protein_group)

            # Find highest scoring PSM
            best_idx = argmax(group.prec_prob)
            best_score = group.prec_prob[best_idx]
            precursor_idx = group.precursor_idx[best_idx]

            # Get cv_fold from library (more reliable than PSM file)
            cv_fold = getCvFold(precursors, precursor_idx)

            # Update if this is the best score for this protein
            value = (best_score = best_score, cv_fold = cv_fold)
            if !haskey(protein_to_cv_fold, protein_name)
                insert!(protein_to_cv_fold, protein_name, value)
            elseif best_score > protein_to_cv_fold[protein_name].best_score
                protein_to_cv_fold[protein_name] = value
            end
        end
    end

    return protein_to_cv_fold
end

#=
"""
    build_protein_cv_fold_mapping(psm_paths::Vector{String}, precursors::LibraryPrecursors)
    -> Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}

Build a mapping from protein names to CV fold assignments using RANDOM assignment.

**TEMPORARY IMPLEMENTATION**: Replaces peptide-based assignment for testing.
Each unique protein is randomly assigned to fold 0 or 1 with equal probability.

# Arguments
- `psm_paths`: Vector of paths to PSM files (used to discover proteins)
- `precursors`: Library precursors (not used in random version, kept for signature compatibility)

# Returns
- Dictionary mapping protein_name to named tuple with:
  - `best_score`: Set to 1.0f0 (dummy value, not used in CV)
  - `cv_fold`: Randomly assigned to 0 or 1

# Process
1. Scans PSM files to collect all unique protein names
2. Randomly assigns each protein to fold 0 or 1
3. Uses fixed random seed (1234) for reproducibility

# Notes
- **This is a simplified replacement for the peptide-based assignment**
- All instances of the same protein get the same CV fold
- Random seed ensures reproducible results across runs
"""
function build_protein_cv_fold_mapping(
    psm_paths::Vector{String},
    precursors::LibraryPrecursors
)

    # Fixed seed for reproducibility
    rng = MersenneTwister(1776)

    # Create mapping: protein_name -> (best_score=1.0, cv_fold=random)
    protein_to_cv_fold = Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}()

    # Collect all unique protein names across all PSM files
    all_proteins = Set{String}()

    for psm_path in psm_paths
        # Skip if PSM file doesn't exist
        if !isfile(psm_path)
            @user_warn "PSM file not found: $psm_path"
            continue
        end

        # Load PSM data
        psms = DataFrame(Arrow.Table(psm_path))

        # Verify required column exists
        if !hasproperty(psms, :inferred_protein_group)
            @user_warn "PSM file missing :inferred_protein_group column: $psm_path"
            continue
        end

        # Filter for valid PSMs with protein group assignments
        psms = filter(row -> !ismissing(row.inferred_protein_group), psms)

        # Add unique proteins to the set
        for protein_name in psms.inferred_protein_group
            push!(all_proteins, protein_name)
        end
    end

    # Convert to sorted vector for deterministic iteration
    protein_names = sort(collect(all_proteins))

    # Randomly assign each protein to fold 0 or 1
    n_proteins = length(protein_names)

    @user_info "Random CV fold assignment: assigning $n_proteins unique proteins to 2 folds"

    for protein_name in protein_names
        # Randomly choose fold 0 or 1
        cv_fold = rand(rng, UInt8[0, 1])

        # Create entry with dummy best_score
        value = (best_score = 1.0f0, cv_fold = cv_fold)
        insert!(protein_to_cv_fold, protein_name, value)
    end

    # Report fold distribution
    fold_0_count = count(p -> p.cv_fold == 0, values(protein_to_cv_fold))
    fold_1_count = count(p -> p.cv_fold == 1, values(protein_to_cv_fold))

    @user_info "Random CV fold distribution: Fold 0: $fold_0_count proteins, Fold 1: $fold_1_count proteins"

    return protein_to_cv_fold
end
=#
"""
    assign_protein_group_cv_folds!(all_protein_groups::DataFrame, 
                                  protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}})

Assign CV fold to each protein group based on a pre-built mapping.

# Arguments
- `all_protein_groups`: DataFrame of protein groups to update (modified in-place)
- `protein_to_cv_fold`: Pre-built dictionary mapping protein names to cv_fold assignments

# Process
Adds cv_fold column to protein groups DataFrame based on the mapping
"""
function assign_protein_group_cv_folds!(
    all_protein_groups::DataFrame,
    protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}
)
    # Assign cv_fold to protein groups using the pre-built mapping
    cv_folds = Vector{UInt8}(undef, nrow(all_protein_groups))
    missing_count = 0
    
    for (i, protein_name) in enumerate(all_protein_groups.protein_name)
        if haskey(protein_to_cv_fold, protein_name)
            cv_folds[i] = protein_to_cv_fold[protein_name].cv_fold
        else
            missing_count += 1
            # Default to fold 0 if no matching peptide found
            cv_folds[i] = UInt8(0)
        end
    end
    
    if missing_count > 0
        @user_warn "There were $missing_count protein groups without matching peptides, assigned to fold 0"
    end

    all_protein_groups[!, :cv_fold] = cv_folds
end

"""
    apply_probit_scores_multifold!(pg_refs::Vector{ProteinGroupFileReference},
                                  protein_to_cv_fold::Dict{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}},
                                  models::Dict{UInt8, Vector{Float64}},
                                  feature_names::Vector{Symbol})

Apply probit models to protein group files based on their CV fold.

# Arguments
- `pg_refs`: Vector of protein group file references to update
- `protein_to_cv_fold`: Pre-built mapping of protein names to cv_fold
- `models`: Dictionary mapping CV fold to fitted model coefficients
- `feature_names`: Feature names used in the model
"""
function apply_probit_scores_multifold!(
    pg_refs::Vector{ProteinGroupFileReference},
    protein_to_cv_fold::Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}},
    models::Dict{UInt8, Vector{Float64}},
    feature_names::Vector{Symbol};
    skip_scoring = false
)
    for ref in pg_refs
        transform_and_write!(ref) do df
            # Assign CV folds using pre-built mapping
            cv_folds = Vector{UInt8}(undef, nrow(df))
            for (i, protein_name) in enumerate(df.protein_name)
                if haskey(protein_to_cv_fold, protein_name)
                    cv_folds[i] = protein_to_cv_fold[protein_name].cv_fold
                else
                    # Default to fold 0 if not found
                    cv_folds[i] = UInt8(0)
                end
            end
            df[!, :cv_fold] = cv_folds
            
            # Save original scores for comparison
            df[!, :old_pg_score] = copy(df.pg_score)
            
            # Apply appropriate model to each fold (skip if skip_scoring = true)
            if !skip_scoring
                for (fold, model) in models
                    mask = df.cv_fold .== fold
                    if sum(mask) > 0
                        X = Matrix{Float64}(df[mask, feature_names])
                        df[mask, :pg_score] = Float32.(calculate_probit_scores(X, model))
                    end
                end
            else
                # When skip_scoring is true, pg_score contains log-sum scores: -sum(log1p(-p))
                # Convert to probabilities: p = 1 - exp(-pg_score)
                df[!, :pg_score] = 1.0f0 .- exp.(-df.pg_score)
                # Clamp to avoid numerical issues with logodds function
                df[!, :pg_score] = clamp.(df.pg_score, 1f-6, 1f0 - 1f-6)
            end
            
            # Sort by pg_score and target in descending order
            sort!(df, [:pg_score, :target], rev = [true, true])
            
            # Remove temporary cv_fold column before returning
            select!(df, Not(:cv_fold))
            
            return df
        end
    end
end

"""
    perform_probit_analysis_multifold(all_protein_groups::DataFrame,
                                     qc_folder::String,
                                     pg_refs::Vector{ProteinGroupFileReference},
                                     precursors::LibraryPrecursors;
                                     protein_to_cv_fold::Union{Nothing, Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}} = nothing,
                                     show_improvement = true,
                                     skip_scoring = false)

Perform probit regression analysis with automatic CV fold detection from library.

# Arguments
- `all_protein_groups`: DataFrame with protein group data
- `qc_folder`: Folder for QC plots
- `pg_refs`: Protein group file references
- `precursors`: Library precursors containing CV fold information
- `protein_to_cv_fold`: Optional pre-built mapping of proteins to CV folds
- `show_improvement`: Whether to report improvement metrics

# Process
1. Detects unique CV folds from library
2. Assigns CV folds to protein groups based on mapping (or builds it if not provided)
3. Trains separate probit model for each fold
4. Applies models to held-out data
5. Updates protein group files if provided
"""
function perform_probit_analysis_multifold(
    all_protein_groups::DataFrame,
    qc_folder::String,
    pg_refs::Vector{ProteinGroupFileReference},
    precursors::LibraryPrecursors;
    protein_to_cv_fold::Union{Nothing, Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}} = nothing,
    show_improvement = true,
    skip_scoring = false,
    ms1_scoring::Bool = true
)

    #skip_scoring = true 
    #@user_info "Skipped scoring!!!"
    # 1. Detect unique CV folds from library
    unique_cv_folds = detect_unique_cv_folds(precursors)
    n_folds = length(unique_cv_folds)
    
    # 2. Use provided mapping or build it if not provided
    if protein_to_cv_fold === nothing
        # Build the mapping from PSM files
        psm_paths = [get_corresponding_psm_path(ref) for ref in pg_refs]
        protein_to_cv_fold = build_protein_cv_fold_mapping(psm_paths, precursors)
    end
    
    # 3. Assign CV folds to protein groups based on the mapping
    assign_protein_group_cv_folds!(all_protein_groups, protein_to_cv_fold)
    
    # 3. Check distribution
    fold_counts = Dict{UInt8, Int}()
    for fold in all_protein_groups.cv_fold
        fold_counts[fold] = get(fold_counts, fold, 0) + 1
    end

    # 4. Define features (same as original)
    #feature_names = [:pg_score, :peptide_coverage, :n_possible_peptides] #:any_common_peps]
    feature_names = [:pg_score,
    #:peptide_coverage, :n_possible_peptides,
    :any_common_peps, :all_precursors_mbr]
    # Apply feature filtering
    adjust_any_common_peps!(feature_names, all_protein_groups)
    remove_zero_variance_columns!(feature_names, all_protein_groups)

    if isempty(feature_names)
        @user_warn "No valid features remaining for multifold probit regression after filtering"
        return
    end
    
    # 5. Train probit model for each fold (skip if skip_scoring = true)
    models = Dict{UInt8, Vector{Float64}}()

    if !skip_scoring
        # Determine positive training examples using 1% FDR on pg_score
        n_proteins = nrow(all_protein_groups)
        qvals = Vector{Float32}(undef, n_proteins)
        get_qvalues!(all_protein_groups.pg_score, all_protein_groups.target, qvals)
        passing_mask = (qvals .<= 0.01f0) .& all_protein_groups.target
        train_mask_FDR = passing_mask .| .!all_protein_groups.target

        for test_fold in unique_cv_folds
            # Get training data (all folds except test_fold)
            train_mask = (all_protein_groups.cv_fold .!= test_fold) #.& train_mask_FDR
            
            # Check if we have sufficient data
            n_train_targets = sum(all_protein_groups[train_mask, :target])
            n_train_decoys = sum(.!all_protein_groups[train_mask, :target])
            
            if n_train_targets < 10 || n_train_decoys < 10
                @user_warn "Insufficient training data for fold $test_fold" targets=n_train_targets decoys=n_train_decoys
                continue
            end
            
            X_train = Matrix{Float64}(all_protein_groups[train_mask, feature_names])
            y_train = all_protein_groups[train_mask, :target]
            
            # Fit model
            β_fitted = fit_probit_model(X_train, y_train)
            models[test_fold] = β_fitted
        end
    end
    
    # 6. Apply models to their respective test folds (skip if skip_scoring = true)
    all_protein_groups[!, :old_pg_score] = copy(all_protein_groups[!, :pg_score])  # Save for comparison
    
    if !skip_scoring
        for test_fold in unique_cv_folds
            if !haskey(models, test_fold)
                @user_warn "No model available for fold $test_fold, keeping original scores"
                continue
            end
            
            test_mask = all_protein_groups.cv_fold .== test_fold
            n_test = sum(test_mask)
            
            if n_test > 0
                X_test = Matrix{Float64}(all_protein_groups[test_mask, feature_names])
                prob_scores = calculate_probit_scores(X_test, models[test_fold])
                all_protein_groups[test_mask, :pg_score] = Float32.(prob_scores)
            end
        end
    end
    
    # 7. Report improvement if requested (skip if skip_scoring = true)
    if show_improvement && !skip_scoring
        # Calculate improvement at 1% FDR
        old_qvalues = zeros(Float32, nrow(all_protein_groups))
        new_qvalues = zeros(Float32, nrow(all_protein_groups))
        
        get_qvalues!(all_protein_groups[!, :old_pg_score], all_protein_groups[!, :target], old_qvalues)
        get_qvalues!(all_protein_groups[!, :pg_score], all_protein_groups[!, :target], new_qvalues)
        
        old_passing = sum((old_qvalues .<= 0.01f0) .& all_protein_groups.target)
        new_passing = sum((new_qvalues .<= 0.01f0) .& all_protein_groups.target)
        
        percent_improvement = 0.0
        if old_passing > 0
            percent_improvement = round(100.0 * (new_passing - old_passing) / old_passing, digits=1)
        end
    end
    
    # 8. Update protein group files if provided
    if !isempty(pg_refs)
        apply_probit_scores_multifold!(pg_refs, protein_to_cv_fold, models, feature_names; skip_scoring = skip_scoring)
    end
    
    # Clean up temporary column
    select!(all_protein_groups, Not(:cv_fold))
end
