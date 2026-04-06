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

# ==========================================================
# Protein group analysis
# ==========================================================

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
    peptide_coverage_logits = Vector{Float32}(undef, n_groups)
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
        peptide_coverage_logits[i] = smoothed_coverage_logit(group.features.n_peptides, group.features.n_possible_peptides)
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
        peptide_coverage_logit = peptide_coverage_logits,
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
