"""
    infer_proteins(proteins::Vector{String}, peptides::Vector{String})::Dictionary{String, Tuple{String, Bool}}

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
# Parameters
- `proteins::Vector{String}`: Protein IDs for each peptide. Each element can contain multiple protein 
  identifiers separated by semicolons (e.g., "A;B;C") indicating that the peptide maps to multiple proteins.
- `peptides::Vector{String}`: Peptide sequence identifiers corresponding to the proteins vector.

# Returns
- `Dictionary{String, Tuple{String, Bool}}`: Maps each peptide to a tuple of:
  - Protein group assignment (either a single protein ID or a semicolon-separated list)
  - Boolean indicating whether the peptide should be used for protein inference (true) or is ambiguous (false)

# Examples
```julia
# Case A: Distinct proteins
proteins = ["A", "A", "B", "B"]
peptides = ["pep1", "pep2", "pep3", "pep4"]
result = infer_proteins(proteins, peptides)
# Returns:
# "pep1" => ("A", true)
# "pep2" => ("A", true)
# "pep3" => ("B", true)
# "pep4" => ("B", true)

# Case B: Differentiable proteins
proteins = ["A", "A;B", "A;B", "B"]
peptides = ["pep1", "pep2", "pep3", "pep4"]
result = infer_proteins(proteins, peptides)
# Returns:
# "pep1" => ("A", true)
# "pep2" => ("A;B", false)
# "pep3" => ("A;B", false)
# "pep4" => ("B", true)

# Case C: Indistinguishable proteins
proteins = ["A;B", "A;B", "A;B", "A;B"]
peptides = ["pep1", "pep2", "pep3", "pep4"]
result = infer_proteins(proteins, peptides)
# Returns:
# "pep1" => ("A;B", true)
# "pep2" => ("A;B", true)
# "pep3" => ("A;B", true)
# "pep4" => ("A;B", true)

# Case D: Subset proteins
proteins = ["A", "A;B", "A;B", "A;B"]
peptides = ["pep1", "pep2", "pep3", "pep4"]
result = infer_proteins(proteins, peptides)
# Returns:
# "pep1" => ("A", true)
# "pep2" => ("A;B", false)
# "pep3" => ("A;B", false)
# "pep4" => ("A;B", false)

# Case E: Subsumable proteins
proteins = ["A", "A;B", "B;C", "C"]
peptides = ["pep1", "pep2", "pep3", "pep4"]
result = infer_proteins(proteins, peptides)
# Returns:
# "pep1" => ("A", true)
# "pep2" => ("A;B", false)
# "pep3" => ("B;C", false)
# "pep4" => ("C", true)

# Case F: Proteins with shared peptides only
proteins = ["A;B", "A;B;C", "A;B;C", "A;C"]
peptides = ["pep1", "pep2", "pep3", "pep4"]
result = infer_proteins(proteins, peptides)
# Returns:
# "pep1" => ("A;B;C", true)
# "pep2" => ("A;B;C", true)
# "pep3" => ("A;B;C", true)
# "pep4" => ("A;B;C", true)

# Case G: Complex case
proteins = ["A", "A;B", "B;C", "B;C"]
peptides = ["pep1", "pep2", "pep3", "pep4"]
result = infer_proteins(proteins, peptides)
# Returns:
# "pep1" => ("A", true)
# "pep2" => ("A;B", false)
# "pep3" => ("B;C", false)
# "pep4" => ("B;C", false)

# Test Case H:
#Combination of all previous test cases 
proteins = ["A","A","B","B",
"C","C;D","C;D","D",
"E;F","E;F","E;F","E;F",
"G","G;H","G;H","G;H",
"I","I;J","J;K","K",
"L;M","L;M;N","L;M;N","L;N",
"O","O;P","P;Q","P;Q"]
peptides = ["pep1","pep2","pep3","pep4",
            "pep5","pep6","pep7","pep8",
            "pep9","pep10","pep11","pep12",
            "pep13","pep14","pep15","pep16",
            "pep17","pep18","pep19","pep20",
            "pep21","pep22","pep23","pep24",
            "pep25","pep26","pep27","pep28"]
dict_h = Dictionary{String, Tuple{String, Bool}}()
insert!(dict_h, "pep1", ("A", true))
insert!(dict_h, "pep2", ("A", true))
insert!(dict_h, "pep3", ("B", true))
insert!(dict_h, "pep4", ("B", true))
insert!(dict_h, "pep5", ("C", true))
insert!(dict_h, "pep6", ("C;D", false))
insert!(dict_h, "pep7", ("C;D", false))
insert!(dict_h, "pep8", ("D", true))
insert!(dict_h, "pep9", ("E;F", true))
insert!(dict_h, "pep10", ("E;F", true))
insert!(dict_h, "pep11", ("E;F", true))
insert!(dict_h, "pep12", ("E;F", true))
insert!(dict_h, "pep13", ("G", true))
insert!(dict_h, "pep14", ("G;H", false))
insert!(dict_h, "pep15", ("G;H", false))
insert!(dict_h, "pep16", ("G;H", false))
insert!(dict_h, "pep17", ("I", true))
insert!(dict_h, "pep18", ("I;J", false))
insert!(dict_h, "pep19", ("J;K", false))
insert!(dict_h, "pep20", ("K", true))
insert!(dict_h, "pep21", ("L;M;N", true))
insert!(dict_h, "pep22", ("L;M;N", true))
insert!(dict_h, "pep23", ("L;M;N", true))
insert!(dict_h, "pep24", ("L;M;N", true))
insert!(dict_h, "pep25", ("O", true))
insert!(dict_h, "pep26", ("O;P", false))
insert!(dict_h, "pep27", ("P;Q", false))
insert!(dict_h, "pep28", ("P;Q", false))
solution = infer_proteins(proteins, peptides)
@test all([dict_h[key] == solution[key] for (key, val) in pairs(solution)])
```

See Also

Nesvizhskii AI, Aebersold R. Interpretation of shotgun proteomic data: the protein inference problem.
Mol Cell Proteomics. 2005 Oct;4(10):1419-40. doi: 10.1074/mcp.R500012-MCP200.
Zhang B, Chambers MC, Tabb DL. Proteomic parsimony through bipartite graph analysis improves accuracy
and transparency. J Proteome Res. 2007 Sep;6(9):3549-57. doi: 10.1021/pr070230d.
"""
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
