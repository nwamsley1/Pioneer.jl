function infer_proteins(proteins::Vector{String}, peptides::Vector{String})::Dictionary{String, Tuple{String, Bool}}
    # Build peptide-to-protein and protein-to-peptide mappings
    peptide_to_proteins = Dictionary{String, Set{String}}()
    original_groups = Dictionary{String, String}()
    
    for i in 1:length(peptides)
        peptide = peptides[i]
        if !haskey(peptide_to_proteins, peptide)
            insert!(peptide_to_proteins, peptide, Set{String}())
            insert!(original_groups, peptide, proteins[i])
        end
        
        for protein in split(proteins[i], ";")
            push!(peptide_to_proteins[peptide], protein)
        end
    end
    
    protein_to_peptides = Dictionary{String, Set{String}}()
    
    for (peptide, protein_set) in pairs(peptide_to_proteins)
        for protein in protein_set
            if !haskey(protein_to_peptides, protein)
                insert!(protein_to_peptides, protein, Set{String}())
            end
            push!(protein_to_peptides[protein], peptide)
        end
    end
    
    # Find connected components (independent protein-peptide clusters)
    visited_peptides = Set{String}()
    components = Vector{Tuple{Set{String}, Set{String}}}()  # (peptides, proteins)
    
    for peptide in peptides
        if peptide in visited_peptides
            continue
        end
        
        component_peptides = Set{String}()
        component_proteins = Set{String}()
        queue = [peptide]
        
        while !isempty(queue)
            current = pop!(queue)
            
            if current in visited_peptides
                continue
            end
            
            push!(component_peptides, current)
            push!(visited_peptides, current)
            
            # Add proteins connected to this peptide
            for protein in peptide_to_proteins[current]
                push!(component_proteins, protein)
                
                # Add peptides connected to this protein
                for p in protein_to_peptides[protein]
                    if !(p in visited_peptides)
                        push!(queue, p)
                    end
                end
            end
        end
        
        push!(components, (component_peptides, component_proteins))
    end
    
    # Initialize result dictionary
    result = Dictionary{String, Tuple{String, Bool}}()
    println("components $components")
    println("result $result")
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
                all_proteins = sort(collect(component_proteins))
                protein_group = join(all_proteins, ";")
                
                for peptide in component_peptides
                    insert!(result, peptide, (protein_group, true))
                end
                continue
            end
        end
        
        # Find peptides unique to a protein
        unique_peptide_to_protein = Dictionary{String, String}()
        
        for peptide in component_peptides
            proteins_for_peptide = intersect(peptide_to_proteins[peptide], component_proteins)
            if length(proteins_for_peptide) == 1
                insert!(unique_peptide_to_protein, peptide, first(proteins_for_peptide))
            end
        end
        
        # Case F handling: No protein has unique peptides
        if isempty(unique_peptide_to_protein) && !isempty(component_proteins)
            # If no protein has unique peptides, merge all
            all_proteins = sort(collect(component_proteins))
            protein_group = join(all_proteins, ";")
            
            for peptide in component_peptides
                insert!(result, peptide, (protein_group, true))
            end
            continue
        end
        
        # Apply greedy set cover for minimal protein list
        remaining_peptides = copy(component_peptides)
        necessary_proteins = Set{String}()
        
        # First include proteins with unique peptides
        unique_proteins = Set{String}()
        for (_, protein) in pairs(unique_peptide_to_protein)
            push!(unique_proteins, protein)
        end
        
        for protein in unique_proteins
            push!(necessary_proteins, protein)
            for peptide in intersect(protein_to_peptides[protein], remaining_peptides)
                delete!(remaining_peptides, peptide)
            end
        end
        
        # Continue with greedy set cover for remaining peptides
        candidate_proteins = collect(setdiff(component_proteins, necessary_proteins))
        
        while !isempty(remaining_peptides) && !isempty(candidate_proteins)
            # Find protein that covers the most remaining peptides
            best_protein = ""
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
            
            push!(necessary_proteins, best_protein)
            filter!(p -> p != best_protein, candidate_proteins)
            
            for peptide in intersect(protein_to_peptides[best_protein], remaining_peptides)
                delete!(remaining_peptides, peptide)
            end
        end
        
        # Assign peptides to proteins
        for peptide in component_peptides
            proteins_for_peptide = intersect(peptide_to_proteins[peptide], component_proteins)
            
            if length(proteins_for_peptide) == 1
                # Unique peptide for a single protein - mark as true
                insert!(result, peptide, (first(proteins_for_peptide), true))
            else
                # Shared peptide - use original protein group and mark as false
                original_group = original_groups[peptide]
                insert!(result, peptide, (original_group, false))
            end
        end
    end
    
    return result
end


"""
    protein_inference_tests()
Test protein inference algorithm. Based on Figure 5 from Nesvishkii and Aebersold "Interpretation of Shotgun Proteomic Data",
Mol & Cell Proteomics 4:1419-1140, 2005
"""
