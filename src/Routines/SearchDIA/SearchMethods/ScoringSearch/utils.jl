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
1. Accumulates scores across files
2. Selects highest scoring trace per precursor
3. Returns set of best precursor-isotope combinations
"""
function get_best_traces(
    second_pass_psms_paths::Vector{String},
    min_prob::Float32 = 0.75f0
)
    start_time = time()
    initial_memory = Base.gc_live_bytes() / 1024^2  # MB
    
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
        
        for i in range(1, n_rows)
            psms_key = (precursor_idx = psms_table[:precursor_idx][i],  isotopes_captured = psms_table[:isotopes_captured][i])

            row_score = psms_table[:prob][i]
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

    psms_trace_df = DataFrame(
    (precursor_idx = [key[:precursor_idx] for key in keys(psms_trace_scores)],
    isotopes_captured = [key[:isotopes_captured] for key in keys(psms_trace_scores)],
    score = [val for val in values(psms_trace_scores)])
    );
    psms_trace_df[!,:best_trace] .= false;
    gpsms = groupby(psms_trace_df,:precursor_idx)
    for (precursor_idx, psms) in pairs(gpsms)
        psms[argmax(psms[!,:score]),:best_trace] = true
    end
    filter!(x->x.best_trace, psms_trace_df);
    traces_passing = Set([(precursor_idx = x.precursor_idx, isotopes_captured = x.isotopes_captured) for x in eachrow(psms_trace_df)]);
    return traces_passing
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
            @warn "Less than 20 bins to estimate PEP. PEP results suspect..."
        end
        return UniformSpline(bin_target_fraction, bin_mean_prob, 3, 3)
    catch
        @warn "Failed to estimate PEP spline"
        return UniformSpline(SVector{4, Float32}([0, 0, 0, 0]), 3, 0.0f0, 1.0f0, 100.0f0)
    end
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
        psms_scores = unique(psms_scores)
    end

    Q = size(psms_scores, 1)
    M = ceil(Int, Q / min_pep_points_per_bin)
    bin_qval, bin_mean_prob = Vector{Float32}(undef, M), Vector{Float32}(undef, M)
    bin_size = 0
    bin_idx = 0
    mean_prob, targets, decoys = 0.0f0, 0, 0
    targets, decoys = 0, 0
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
        if bin_size == min_pep_points_per_bin
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
        @warn "Insufficient unique points for q-value interpolation, using default"
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
        @warn "No protein groups to write to $output_path"
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
    
    @info "Wrote $(n_groups) protein groups to $output_path"
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
                                    protein_to_cv_fold::Union{Nothing, Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}} = nothing)

Perform probit regression on protein groups.

# Arguments
- `pg_refs`: Vector of protein group file references
- `max_psms_in_memory`: Memory limit for in-memory vs OOM processing
- `qc_folder`: Folder for QC plots
- `precursors`: Library precursors
- `protein_to_cv_fold`: Optional pre-built mapping of proteins to CV folds
"""
function perform_protein_probit_regression(
    pg_refs::Vector{ProteinGroupFileReference},
    max_psms_in_memory::Int64,
    qc_folder::String,
    precursors::LibraryPrecursors;
    protein_to_cv_fold::Union{Nothing, Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}} = nothing
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
    
    if total_protein_groups > max_protein_groups_in_memory_limit
        @info "Using out-of-memory probit regression (exceeds limit of $max_protein_groups_in_memory_limit)"
        perform_probit_analysis_oom(pg_refs, total_protein_groups, max_protein_groups_in_memory_limit, qc_folder)
    else
        @info "Using in-memory probit regression"
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
        
        if n_targets > 0 && n_decoys > 0 && nrow(all_protein_groups) > 10
            # Perform probit regression analysis
            @info "Performing Probit Analysis (targets: $n_targets, decoys: $n_decoys)"
            #perform_probit_analysis(all_protein_groups, qc_folder, pg_refs)
            perform_probit_analysis_multifold(
                all_protein_groups,
                qc_folder,
                pg_refs,
                precursors;
                protein_to_cv_fold = protein_to_cv_fold
            )
        else
            @info "Skipping Probit analysis: insufficient data (targets: $n_targets, decoys: $n_decoys)"
        end
    end
end

"""
    update_psms_with_probit_scores_refs(paired_refs::Vector{PairedSearchFiles},
                                       acc_to_max_pg_score::Dict{ProteinKey,Float32},
                                       pg_score_to_qval::Interpolations.Extrapolation,
                                       global_pg_score_to_qval::Interpolations.Extrapolation)

Update PSMs with probit-scored pg_score values and q-values using references.

# Arguments
- `paired_refs`: Paired PSM/protein group file references
- `acc_to_max_pg_score`: Dictionary mapping protein keys to global scores
- `pg_score_to_qval`: Interpolation function for pg_score to q-value
- `global_pg_score_to_qval`: Interpolation function for global_pg_score to q-value
"""
function update_psms_with_probit_scores_refs(
    paired_refs::Vector{PairedSearchFiles},
    acc_to_max_pg_score::Dict{ProteinKey,Float32},
    pg_score_to_qval::Interpolations.Extrapolation,
    global_pg_score_to_qval::Interpolations.Extrapolation
)

    total_psms_updated = 0
    files_processed = 0
    
    for paired_ref in paired_refs
        psm_ref = paired_ref.psm_ref
        pg_ref = paired_ref.protein_ref
        
        # Verify both files exist
        if !exists(psm_ref) || !exists(pg_ref)
            @warn "Missing file" psm_path=file_path(paired_ref.psm_ref) pg_path=file_path(paired_ref.protein_ref)
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
            
            # Update pg_score column with probit scores
            probit_pg_scores = Vector{Union{Missing, Float32}}(undef, n_psms)
            global_pg_scores = Vector{Union{Missing, Float32}}(undef, n_psms)
            pg_qvals = Vector{Union{Missing, Float32}}(undef, n_psms)
            global_pg_qvals = Vector{Union{Missing, Float32}}(undef, n_psms)
            pg_peps = Vector{Union{Missing, Float32}}(undef, n_psms)
            
            for i in 1:n_psms
                # Skip if missing inferred protein group
                if ismissing(psms_df[i, :inferred_protein_group])
                    throw("Missing Inferred Protein Group!!!")
                end
                
                # Skip if peptide didn't match to a distinct protein group
                if psms_df[i,:use_for_protein_quant] == false
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
                
                if !haskey(acc_to_max_pg_score, key)
                    throw("Missing global pg score lookup key!!!")
                end
                global_pg_scores[i] = acc_to_max_pg_score[key]
                
                # Calculate q-values
                pg_qvals[i] = pg_score_to_qval(probit_pg_scores[i])
                global_pg_qvals[i] = global_pg_score_to_qval(global_pg_scores[i])
            end
            
            # Update columns
            psms_df[!, :pg_score] = probit_pg_scores
            psms_df[!, :global_pg_score] = global_pg_scores
            psms_df[!, :pg_qval] = pg_qvals
            psms_df[!, :global_qval_pg] = global_pg_qvals
            psms_df[!, :pg_pep] = pg_peps
            
            total_psms_updated += n_psms
            
            return psms_df
        end
        
        files_processed += 1
    end
end





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
                                    max_protein_groups_in_memory::Int, qc_folder::String)
    
    # Calculate sampling ratio
    sampling_ratio = max_protein_groups_in_memory / total_protein_groups
    @info "Sampling ratio: $(round(sampling_ratio * 100, digits=2))%"
    
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
                append!(sampled_protein_groups, df[sample_indices, :])
            end
        end
    end
    
    @info "Sampled $(nrow(sampled_protein_groups)) protein groups for training"
    
    # Define features to use
    feature_names = [:pg_score, :peptide_coverage, :n_possible_peptides, :log_binom_coeff]
    X = Matrix{Float64}(sampled_protein_groups[:, feature_names])
    y = sampled_protein_groups.target
    
    # Fit probit model on sampled data
    #β_fitted, X_mean, X_std = fit_probit_model(X, y)
    β_fitted = fit_probit_model(X, y)
    @info "Fitted probit model coefficients: $β_fitted"
    
    # Now apply the model to all protein groups file by file
    # and collect statistics
    total_targets = 0
    total_decoys = 0
    
    # Process each file
    for ref in pg_refs
        if exists(ref)
            # Load file
            df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
            
            # Calculate probit scores
            X_file = Matrix{Float64}(df[:, feature_names])
            prob_scores = calculate_probit_scores(X_file, β_fitted)#, X_mean, X_std)
            
            # Overwrite pg_score with probit scores
            df[!, :pg_score] = Float32.(prob_scores)
            
            # Sort by pg_score and target in descending order
            sort!(df, [:pg_score, :target], rev = [true, true])
            
            # Count statistics
            total_targets += sum(df.target)
            total_decoys += sum(.!df.target)
            
            # Write back the file with probit scores
            writeArrow(file_path(ref), df)
        end
    end
    
    # Report results
    @info "Out-of-Memory Probit Regression Results:"
    @info "  Total protein groups: $(total_targets + total_decoys) (targets: $total_targets, decoys: $total_decoys)"
    @info "  Training sample size: $(nrow(sampled_protein_groups))"
    @info "  Model features: $(join(feature_names, ", "))"
    @info "  Model coefficients: $(round.(β_fitted, digits=3))"
end

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
    @info "In memory probit regression analysis" n_targets=n_targets n_decoys=n_decoys total_protein_groups=nrow(all_protein_groups)
    # Define features to use
    feature_names = [:pg_score, :peptide_coverage, :n_possible_peptides] # :log_binom_coeff] 
    X = Matrix{Float64}(all_protein_groups[:, feature_names])
    y = all_protein_groups.target
    
    # Fit probit model
    #β_fitted, X_mean, X_std = fit_probit_model(X, y)
    β_fitted = fit_probit_model(X, y)
    # Report basic model statistics
    @info "Probit Regression completed:"
    @info "  Total protein groups: $(n_targets + n_decoys) (targets: $n_targets, decoys: $n_decoys)"
    @info "  Model features: $(join(feature_names, ", "))"
    @info "  Model coefficients: $(round.(β_fitted, digits=3))"
    
    # Create decision boundary plots
    # TODO: Fix plotting type error
    # plot_probit_decision_boundary(all_protein_groups, β_fitted, X_mean, X_std, feature_names, qc_folder)
    
    # Re-process individual files if references are provided
    if !isempty(pg_refs)
        @info "Re-processing individual protein group files with probit scores"
        # Use the new apply_probit_scores! function with references
        apply_probit_scores!(pg_refs, β_fitted, feature_names)
    end
    #=
    if show_improvement && !isempty(pg_refs)
        # Merge all protein group files to calculate improvement
        pg_paths = [file_path(ref) for ref in pg_refs]
        all_pgs = vcat([DataFrame(Arrow.Table(path)) for path in pg_paths]...)
        old_qvalues = zeros(Float32, size(all_pgs, 1))
        get_qvalues!(all_pgs[!,:old_pg_score], all_pgs[!,:target], old_qvalues)
        old_passing = sum((old_qvalues .<= 0.01f0) .& all_pgs.target)
        @info "Old passing qvals " sum((old_qvalues .<= 0.01f0) .& all_pgs.target) # Count targets with qval < 0.01
        new_qvalues = zeros(Float32, size(all_pgs, 1))
        get_qvalues!(all_pgs[!,:pg_score], all_pgs[!,:target], new_qvalues)
        new_passing = sum((new_qvalues .<= 0.01f0) .& all_pgs.target)
        @info "New passing qvals " sum((new_qvalues .<= 0.01f0) .& all_pgs.target) # Count targets with qval < 0.01
        percent_improv = 100.0*(new_passing - old_passing)/old_passing |> round 
        @info "Probit regression improved passing targets by $(percent_improv)%"
    end
    =#
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
    # Standardize features
    #X_mean = mean(X, dims=1)
    #X_std = std(X, dims=1)
    #X_std[X_std .== 0] .= 1.0  # Avoid division by zero
    X_standardized =X#(X .- X_mean) ./ X_std
    
    # Add intercept column
    X_with_intercept = hcat(ones(size(X_standardized, 1)), X_standardized)
    X_df = DataFrame(X_with_intercept, [:intercept; Symbol.("feature_", 1:size(X, 2))])
    
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
    Row │ precursor_idx  prob      weight         target  irt_obs    missed_cleavage  isotopes_captured  scan_idx  ms_file_idx  peak_area   new_best_scan 
        │ UInt32         Float32   Float32        Bool    Float32    UInt8            Tuple{Int8, Int8}  UInt32    Int64        Float32     UInt32        
────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
      1 │       1468734  0.808784    201.467        true   0.659221                1  (0, 3)                49270            1     6.03541          50180
      2 │        262434  0.989585   2696.17         true   0.659221                0  (0, 3)                76753            1   121.201            76753
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
    # Create mapping: protein_name -> (best_score=score, cv_fold=fold)
    protein_to_cv_fold = Dictionary{String, @NamedTuple{best_score::Float32, cv_fold::UInt8}}()
    
    # Process each PSM file
    for psm_path in psm_paths
        # Skip if PSM file doesn't exist
        if !isfile(psm_path)
            @warn "PSM file not found: $psm_path"
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
        @warn "There were $missing_count protein groups without matching peptides, assigned to fold 0"
    end
    @info "sum(cv_folds.==0) = $(sum(cv_folds .== 0))"
    @info "sum(cv_folds.==1) = $(sum(cv_folds .== 1))"
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
    feature_names::Vector{Symbol}
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
            
            # Apply appropriate model to each fold
            for (fold, model) in models
                mask = df.cv_fold .== fold
                if sum(mask) > 0
                    X = Matrix{Float64}(df[mask, feature_names])
                    df[mask, :pg_score] = Float32.(calculate_probit_scores(X, model))
                end
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
                                     show_improvement = true)

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
    show_improvement = true
)
    # 1. Detect unique CV folds from library
    unique_cv_folds = detect_unique_cv_folds(precursors)
    n_folds = length(unique_cv_folds)
    
    @info "Multi-fold probit regression analysis" n_folds=n_folds cv_folds=unique_cv_folds
    
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
    @info "Protein group distribution across folds" fold_counts
    
    # 4. Define features (same as original)
    feature_names = [:pg_score, :peptide_coverage, :n_possible_peptides]
    
    # 5. Train probit model for each fold
    models = Dict{UInt8, Vector{Float64}}()
    
    for test_fold in unique_cv_folds
        # Get training data (all folds except test_fold)
        train_mask = all_protein_groups.cv_fold .!= test_fold
        
        # Check if we have sufficient data
        n_train_targets = sum(all_protein_groups[train_mask, :target])
        n_train_decoys = sum(.!all_protein_groups[train_mask, :target])
        
        if n_train_targets < 10 || n_train_decoys < 10
            @warn "Insufficient training data for fold $test_fold" targets=n_train_targets decoys=n_train_decoys
            continue
        end
        
        X_train = Matrix{Float64}(all_protein_groups[train_mask, feature_names])
        y_train = all_protein_groups[train_mask, :target]
        
        # Fit model
        β_fitted = fit_probit_model(X_train, y_train)
        models[test_fold] = β_fitted
        
        @info "Fitted model for fold $test_fold" n_train=sum(train_mask) coefficients=round.(β_fitted, digits=3)
    end
    
    # 6. Apply models to their respective test folds
    all_protein_groups[!, :old_pg_score] = all_protein_groups[!, :pg_score]  # Save for comparison
    
    for test_fold in unique_cv_folds
        if !haskey(models, test_fold)
            @warn "No model available for fold $test_fold, keeping original scores"
            continue
        end
        
        test_mask = all_protein_groups.cv_fold .== test_fold
        n_test = sum(test_mask)
        
        if n_test > 0
            X_test = Matrix{Float64}(all_protein_groups[test_mask, feature_names])
            prob_scores = calculate_probit_scores(X_test, models[test_fold])
            all_protein_groups[test_mask, :pg_score] = Float32.(prob_scores)
            
            @info "Applied model to fold $test_fold" n_test=n_test
        end
    end
    
    # 7. Report improvement if requested
    if show_improvement
        @info "Show improvement "
        # Calculate improvement at 1% FDR
        old_qvalues = zeros(Float32, nrow(all_protein_groups))
        new_qvalues = zeros(Float32, nrow(all_protein_groups))
        
        get_qvalues!(all_protein_groups[!, :old_pg_score], all_protein_groups[!, :target], old_qvalues)
        get_qvalues!(all_protein_groups[!, :pg_score], all_protein_groups[!, :target], new_qvalues)
        
        old_passing = sum((old_qvalues .<= 0.01f0) .& all_protein_groups.target)
        new_passing = sum((new_qvalues .<= 0.01f0) .& all_protein_groups.target)
        
        if old_passing > 0
            percent_improvement = round(100.0 * (new_passing - old_passing) / old_passing, digits=1)
            @info "Multi-fold probit regression results" old_passing=old_passing new_passing=new_passing improvement="$percent_improvement%"
        else
            @info "Multi-fold probit regression results" old_passing=old_passing new_passing=new_passing
        end
    end
    
    # 8. Update protein group files if provided
    if !isempty(pg_refs)
        @info "Updating individual protein group files with probit scores"
        apply_probit_scores_multifold!(pg_refs, protein_to_cv_fold, models, feature_names)
    end
    
    # Clean up temporary column
    select!(all_protein_groups, Not(:cv_fold))
end
