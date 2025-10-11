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
    infer_proteins(proteins::Vector{ProteinKey}, peptides::Vector{PeptideKey})::InferenceResult

Infer minimal set of proteins that explain observed peptides, following the parsimony principle.

This function implements a protein inference algorithm for shotgun proteomics using:
1. Connected component analysis to partition the problem into independent protein-peptide clusters
2. Greedy set cover algorithm to find the minimal set of proteins explaining all peptides
3. Systematic handling of ambiguous peptide assignments

The algorithm handles the following cases (based on Nesvizhskii & Aebersold, 2005):

- Case A: Distinct proteins (proteins with non-overlapping peptide sets)
- Case B: Differentiable proteins (proteins with shared peptides but also unique peptides)
- Case C: Indistinguishable proteins (proteins with identical peptide sets)
- Case D: Subset proteins (when one protein's peptides are a subset of another's)
- Case E: Subsumable proteins (when multiple proteins together can explain all peptides)
- Case F: Protein groups with shared peptides only (no unique peptides for any protein)
- Case G: Complex interconnected proteins (with various shared peptide patterns)
- Case H: Concatenation of cases A-G

# Arguments
- `proteins::Vector{ProteinKey}`: Protein identifiers for each peptide. Each ProteinKey can contain 
  multiple protein identifiers separated by semicolons in the name field, indicating that the peptide 
  maps to multiple proteins.
- `peptides::Vector{PeptideKey}`: Peptide identifiers corresponding to the proteins vector.

# Returns
- `InferenceResult`: Contains two dictionaries:
  - `peptide_to_protein`: Maps each PeptideKey to its assigned ProteinKey
  - `use_for_quant`: Maps each PeptideKey to whether it should be used for quantification

# Examples
```julia
using Pioneer

# Case A: Distinct proteins
protein_keys = [
    ProteinKey("A", true, UInt8(1)),
    ProteinKey("A", true, UInt8(1)),
    ProteinKey("B", true, UInt8(1)),
    ProteinKey("B", true, UInt8(1))
]
peptide_keys = [
    PeptideKey("pep1", true, UInt8(1)),
    PeptideKey("pep2", true, UInt8(1)),
    PeptideKey("pep3", true, UInt8(1)),
    PeptideKey("pep4", true, UInt8(1))
]

result = infer_proteins(protein_keys, peptide_keys)
# Result maps:
# pep1 => A (use_for_quant: true)
# pep2 => A (use_for_quant: true)
# pep3 => B (use_for_quant: true)
# pep4 => B (use_for_quant: true)

# Case B: Differentiable proteins with shared peptides
protein_keys = [
    ProteinKey("A", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("B", true, UInt8(1))
]
peptide_keys = [
    PeptideKey("pep1", true, UInt8(1)),
    PeptideKey("pep2", true, UInt8(1)),
    PeptideKey("pep3", true, UInt8(1)),
    PeptideKey("pep4", true, UInt8(1))
]

result = infer_proteins(protein_keys, peptide_keys)
# Result maps:
# pep1 => A (use_for_quant: true)
# pep2 => A;B (use_for_quant: false)  # Ambiguous peptide
# pep3 => A;B (use_for_quant: false)  # Ambiguous peptide
# pep4 => B (use_for_quant: true)
```

# References
Nesvizhskii AI, Aebersold R. Interpretation of shotgun proteomic data: the protein inference problem.
Mol Cell Proteomics. 2005 Oct;4(10):1419-40. doi: 10.1074/mcp.R500012-MCP200.
Zhang B, Chambers MC, Tabb DL. Proteomic parsimony through bipartite graph analysis improves accuracy
and transparency. J Proteome Res. 2007 Sep;6(9):3549-57. doi: 10.1021/pr070230d.

Notes... Better input. might be Vector{Tuple{ProteinKey, PeptideKey}}. The input is really an edge list. 
This also enforces equal length because unequal length inputs should be invalid. 
""" 
function infer_proteins(
    proteins::Vector{ProteinKey}, 
    peptides::Vector{PeptideKey}
)::InferenceResult
    # Validate input lengths match
    if length(proteins) != length(peptides)
        throw(ArgumentError("proteins and peptides vectors must have the same length"))
    end
    
    # Build peptide-to-protein and protein-to-peptide mappings
    peptide_to_proteins = Dictionary{PeptideKey, Set{ProteinKey}}()
    original_groups = Dictionary{PeptideKey, ProteinKey}()
    
    for i in 1:length(peptides)
        peptide_key = peptides[i]
        protein_key = proteins[i]
        
        if !haskey(peptide_to_proteins, peptide_key)
            insert!(peptide_to_proteins, peptide_key, Set{ProteinKey}())
            # Store the original protein group
            insert!(original_groups, peptide_key, protein_key)
        end
        
        # Split the protein name string by ";" and treat each part as a protein
        for protein_part in split(protein_key.name, ";")
            individual_protein = ProteinKey(
                protein_part, 
                protein_key.is_target, 
                protein_key.entrap_id
            )
            push!(peptide_to_proteins[peptide_key], individual_protein)
        end
    end
    
    # Map from protein to peptides
    protein_to_peptides = Dictionary{ProteinKey, Set{PeptideKey}}()
    
    for (peptide_key, protein_set) in pairs(peptide_to_proteins)
        for protein_key in protein_set
            if !haskey(protein_to_peptides, protein_key)
                insert!(protein_to_peptides, protein_key, Set{PeptideKey}())
            end
            
            push!(protein_to_peptides[protein_key], peptide_key)
        end
    end
    
    # Find connected components (independent protein-peptide clusters)
    visited_peptides = Set{PeptideKey}()
    components = Vector{Tuple{Set{PeptideKey}, Set{ProteinKey}}}()  # (peptides, proteins)
    
    for i in 1:length(peptides)
        peptide_key = peptides[i]
        
        if peptide_key in visited_peptides
            continue
        end
        
        component_peptides = Set{PeptideKey}()
        component_proteins = Set{ProteinKey}()
        queue = [peptide_key]
        
        while !isempty(queue)
            current = pop!(queue)
            
            if current in visited_peptides
                continue
            end
            
            push!(component_peptides, current)
            push!(visited_peptides, current)
            
            # Add proteins connected to this peptide
            for protein_key in peptide_to_proteins[current]
                push!(component_proteins, protein_key)
                
                # Add peptides connected to this protein
                for p in protein_to_peptides[protein_key]
                    if !(p in visited_peptides)
                        push!(queue, p)
                    end
                end
            end
        end
        
        push!(components, (component_peptides, component_proteins))
    end
    
    # Initialize result dictionaries
    peptide_to_protein = Dictionary{PeptideKey, ProteinKey}()
    use_for_quant = Dictionary{PeptideKey, Bool}()
    
    # Process each component independently
    for (component_peptides, component_proteins) in components
        # Check if all proteins have identical peptide sets (Case C/F)
        if length(component_proteins) > 1
            identical_sets = true
            first_protein = first(component_proteins)
            first_peptides = intersect(protein_to_peptides[first_protein], component_peptides)
            
            for protein in component_proteins
                if protein != first_protein
                    protein_peptides = intersect(protein_to_peptides[protein], component_peptides)
                    if protein_peptides != first_peptides
                        identical_sets = false
                        break
                    end
                end
            end
            
            if identical_sets
                # Case where all proteins are indistinguishable
                # Group proteins by their target status and entrapment group
                protein_groups = Dictionary{Tuple{Bool, UInt8}, Vector{String}}()
                
                for protein in component_proteins
                    key = (protein.is_target, protein.entrap_id)
                    if !haskey(protein_groups, key)
                        insert!(protein_groups, key, String[])
                    end
                    push!(protein_groups[key], protein.name)
                end
                
                # Process each group and assign peptides
                for ((is_target, entrap_id), protein_names) in pairs(protein_groups)
                    if !isempty(protein_names)
                        # Get protein names sorted and joined
                        sorted_protein_names = sort(protein_names)
                        protein_group_name = join(sorted_protein_names, ";")
                        
                        # Create final protein key
                        final_protein = ProteinKey(protein_group_name, is_target, entrap_id)
                        
                        # Assign to peptides with matching target status and entrapment group
                        for peptide_key in component_peptides
                            if (is_target == peptide_key.is_target) && (entrap_id == peptide_key.entrap_id)
                                insert!(peptide_to_protein, peptide_key, final_protein)
                                insert!(use_for_quant, peptide_key, true)
                            end
                        end
                    end
                end
                continue
            end
        end
        
        # Apply greedy set cover for minimal protein list
        remaining_peptides = copy(component_peptides)
        necessary_proteins = Set{ProteinKey}()
        
        # Find peptides unique to a protein
        unique_peptide_to_protein = Dictionary{PeptideKey, ProteinKey}()
        
        for peptide_key in component_peptides
            proteins_for_peptide = intersect(peptide_to_proteins[peptide_key], component_proteins)
            if length(proteins_for_peptide) == 1
                insert!(unique_peptide_to_protein, peptide_key, first(proteins_for_peptide))
            end
        end
        
        # Case F handling: No protein has unique peptides
        if isempty(unique_peptide_to_protein) && !isempty(component_proteins)
            # Group proteins by their target status and entrapment group
            protein_groups = Dictionary{Tuple{Bool, UInt8}, Vector{String}}()
            
            for protein in component_proteins
                key = (protein.is_target, protein.entrap_id)
                if !haskey(protein_groups, key)
                    insert!(protein_groups, key, String[])
                end
                push!(protein_groups[key], protein.name)
            end
            
            # Process each group and assign peptides
            for ((is_target, entrap_id), protein_names) in pairs(protein_groups)
                if !isempty(protein_names)
                    # Get protein names sorted and joined
                    sorted_protein_names = sort(protein_names)
                    protein_group_name = join(sorted_protein_names, ";")
                    
                    # Create final protein key
                    final_protein = ProteinKey(protein_group_name, is_target, entrap_id)
                    
                    # Assign to peptides with matching target status and entrapment group
                    for peptide_key in component_peptides
                        if (is_target == peptide_key.is_target) && (entrap_id == peptide_key.entrap_id)
                            insert!(peptide_to_protein, peptide_key, final_protein)
                            insert!(use_for_quant, peptide_key, true)
                        end
                    end
                end
            end
            continue
        end
        
        # First include proteins with unique peptides
        unique_proteins = Set{ProteinKey}()
        for (_, protein) in pairs(unique_peptide_to_protein)
            push!(unique_proteins, protein)
        end
        
        for protein in unique_proteins
            push!(necessary_proteins, protein)
            for peptide_key in intersect(protein_to_peptides[protein], remaining_peptides)
                delete!(remaining_peptides, peptide_key)
            end
        end
        
        # Continue with greedy set cover for remaining peptides
        candidate_proteins = collect(setdiff(component_proteins, necessary_proteins))

        while !isempty(remaining_peptides) && !isempty(candidate_proteins)
            # Step 1: Merge indistinguishable proteins before greedy selection
            # Group proteins by (remaining peptide set, is_target, entrap_id)
            peptide_set_to_proteins = Dictionary{Tuple{UInt64, Bool, UInt8}, Tuple{Set{PeptideKey}, Vector{ProteinKey}}}()

            for protein in candidate_proteins
                remaining_peps = intersect(protein_to_peptides[protein], remaining_peptides)
                # Use hash of sorted peptide set as key for grouping
                pep_set_hash = hash(sort(collect(remaining_peps)))
                group_key = (pep_set_hash, protein.is_target, protein.entrap_id)

                if !haskey(peptide_set_to_proteins, group_key)
                    insert!(peptide_set_to_proteins, group_key, (remaining_peps, ProteinKey[]))
                end
                push!(peptide_set_to_proteins[group_key][2], protein)
            end

            # Create merged candidates
            merged_candidates = ProteinKey[]

            for ((pep_set_hash, is_target, entrap_id), (pep_set, proteins)) in pairs(peptide_set_to_proteins)
                if length(proteins) == 1
                    # No merge needed - single protein
                    push!(merged_candidates, proteins[1])
                else
                    # Merge proteins with identical remaining peptide sets
                    protein_names = sort([p.name for p in proteins])
                    merged_protein = ProteinKey(
                        join(protein_names, ";"),
                        is_target,
                        entrap_id
                    )
                    push!(merged_candidates, merged_protein)

                    # Update protein_to_peptides mapping for merged protein
                    # Use the full peptide set from the first protein (they should all be identical)
                    insert!(protein_to_peptides, merged_protein, protein_to_peptides[proteins[1]])
                end
            end

            candidate_proteins = merged_candidates

            # Step 2: Find protein/group that covers the most remaining peptides
            best_protein = nothing
            best_coverage = 0

            for protein in candidate_proteins
                coverage = length(intersect(protein_to_peptides[protein], remaining_peptides))
                if coverage > best_coverage
                    best_coverage = coverage
                    best_protein = protein
                end
            end

            if best_coverage == 0
                break  # No more peptides can be covered
            end

            # Step 3: Add best protein/group to necessary set
            push!(necessary_proteins, best_protein)
            filter!(p -> p != best_protein, candidate_proteins)

            # Step 4: Remove covered peptides
            for peptide_key in intersect(protein_to_peptides[best_protein], remaining_peptides)
                delete!(remaining_peptides, peptide_key)
            end
        end
        
        # Create a mapping to track peptides that can be uniquely attributed to a protein in the necessary set
        peptide_to_necessary_protein = Dictionary{PeptideKey, ProteinKey}()
        
        # Track peptides that can be attributed to multiple necessary proteins
        ambiguous_peptides = Set{PeptideKey}()
        
        for peptide_key in component_peptides
            # Get all necessary proteins that contain this peptide
            proteins_with_peptide = Set{ProteinKey}()
            for protein in necessary_proteins
                if peptide_key in protein_to_peptides[protein]
                    push!(proteins_with_peptide, protein)
                end
            end
            
            if length(proteins_with_peptide) == 1
                # This peptide is unique to one necessary protein
                insert!(peptide_to_necessary_protein, peptide_key, first(proteins_with_peptide))
            else
                # This peptide is shared among multiple necessary proteins
                push!(ambiguous_peptides, peptide_key)
            end
        end
        
        # Assign peptides to proteins
        for peptide_key in component_peptides
            if haskey(peptide_to_necessary_protein, peptide_key)
                # Peptide is unique to one necessary protein - assign with use_for_quant=true
                protein = peptide_to_necessary_protein[peptide_key]
                insert!(peptide_to_protein, peptide_key, protein)
                insert!(use_for_quant, peptide_key, true)
            else
                # Shared peptide - use original protein group and mark as use_for_quant=false
                original_group = original_groups[peptide_key]
                insert!(peptide_to_protein, peptide_key, original_group)
                insert!(use_for_quant, peptide_key, false)
            end
        end
    end

    return InferenceResult(peptide_to_protein, use_for_quant)
end

