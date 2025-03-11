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
# "pep3" => ("C", true)
# "pep4" => ("C", true)

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
insert!(dict_h, "pep27", ("Q", true))
insert!(dict_h, "pep28", ("Q", ecoli_test_results))
solution = infer_proteins(proteins, peptides)
@test all([dict_h[key] == solution[key] for (key, val) in pairs(solution)])
```

See Also

Nesvizhskii AI, Aebersold R. Interpretation of shotgun proteomic data: the protein inference problem.
Mol Cell Proteomics. 2005 Oct;4(10):1419-40. doi: 10.1074/mcp.R500012-MCP200.
Zhang B, Chambers MC, Tabb DL. Proteomic parsimony through bipartite graph analysis improves accuracy
and transparency. J Proteome Res. 2007 Sep;6(9):3549-57. doi: 10.1021/pr070230d.
"""
function infer_proteins(proteins::Vector{NamedTuple{(:protein_name, :decoy), Tuple{String, Bool}}}, 
                        peptides::Vector{String})::Dictionary{NamedTuple{(:peptide, :decoy), Tuple{String, Bool}}, NamedTuple{(:protein_name, :decoy, :retain), Tuple{String, Bool, Bool}}}
    # Build peptide-to-protein and protein-to-peptide mappings
    peptide_to_proteins = Dictionary{NamedTuple{(:peptide, :decoy), Tuple{String, Bool}}, Set{NamedTuple{(:protein_name, :decoy), Tuple{String, Bool}}}}()
    original_groups = Dictionary{NamedTuple{(:peptide, :decoy), Tuple{String, Bool}}, NamedTuple{(:protein_name, :decoy), Tuple{String, Bool}}}()
    
    for i in 1:length(peptides)
        peptide = peptides[i]
        peptide_key = (peptide = peptide, decoy = proteins[i].decoy)
        
        if !haskey(peptide_to_proteins, peptide_key)
            insert!(peptide_to_proteins, peptide_key, Set{NamedTuple{(:protein_name, :decoy), Tuple{String, Bool}}}())
            # Store the original protein group
            insert!(original_groups, peptide_key, proteins[i])
        end
        
        # Split the protein name string by ";" and treat each part as a protein
        for protein_part in split(proteins[i].protein_name, ";")
            push!(peptide_to_proteins[peptide_key], (protein_name = protein_part, decoy = proteins[i].decoy))
        end
    end
    
    # Map from protein (name + decoy status) to peptides
    protein_to_peptides = Dictionary{NamedTuple{(:protein_name, :decoy), Tuple{String, Bool}}, Set{NamedTuple{(:peptide, :decoy), Tuple{String, Bool}}}}()
    
    for (peptide_key, protein_set) in pairs(peptide_to_proteins)
        for protein_tuple in protein_set
            if !haskey(protein_to_peptides, protein_tuple)
                insert!(protein_to_peptides, protein_tuple, Set{NamedTuple{(:peptide, :decoy), Tuple{String, Bool}}}())
            end
            
            push!(protein_to_peptides[protein_tuple], peptide_key)
        end
    end
    
    # Find connected components (independent protein-peptide clusters)
    visited_peptides = Set{NamedTuple{(:peptide, :decoy), Tuple{String, Bool}}}()
    components = Vector{Tuple{Set{NamedTuple{(:peptide, :decoy), Tuple{String, Bool}}}, Set{NamedTuple{(:protein_name, :decoy), Tuple{String, Bool}}}}}()  # (peptides, proteins)
    
    for i in 1:length(peptides)
        peptide_key = (peptide = peptides[i], decoy = proteins[i].decoy)
        
        if peptide_key in visited_peptides
            continue
        end
        
        component_peptides = Set{NamedTuple{(:peptide, :decoy), Tuple{String, Bool}}}()
        component_proteins = Set{NamedTuple{(:protein_name, :decoy), Tuple{String, Bool}}}()
        queue = [peptide_key]
        
        while !isempty(queue)
            current = pop!(queue)
            
            if current in visited_peptides
                continue
            end
            
            push!(component_peptides, current)
            push!(visited_peptides, current)
            
            # Add proteins connected to this peptide
            for protein_tuple in peptide_to_proteins[current]
                push!(component_proteins, protein_tuple)
                
                # Add peptides connected to this protein
                for p in protein_to_peptides[protein_tuple]
                    if !(p in visited_peptides)
                        push!(queue, p)
                    end
                end
            end
        end
        
        push!(components, (component_peptides, component_proteins))
    end
    
    # Initialize result dictionary
    result = Dictionary{NamedTuple{(:peptide, :decoy), Tuple{String, Bool}}, NamedTuple{(:protein_name, :decoy, :retain), Tuple{String, Bool, Bool}}}()
    
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
                # Group proteins by their decoy status (composite keys already ensure this grouping)
                protein_groups = Dictionary{Bool, Vector{String}}()
                
                for protein in component_proteins
                    if !haskey(protein_groups, protein.decoy)
                        insert!(protein_groups, protein.decoy, String[])
                    end
                    push!(protein_groups[protein.decoy], protein.protein_name)
                end
                
                # Process each group and assign peptides
                for (is_decoy, protein_names) in pairs(protein_groups)
                    if !isempty(protein_names)
                        # Get protein names sorted and joined
                        sorted_protein_names = sort(protein_names)
                        protein_group = join(sorted_protein_names, ";")
                        
                        # Assign to peptides with matching decoy status
                        for peptide_key in component_peptides
                            if peptide_key.decoy == is_decoy
                                insert!(result, peptide_key, (
                                    protein_name = protein_group,
                                    decoy = is_decoy,
                                    retain = true
                                ))
                            end
                        end
                    end
                end
                continue
            end
        end
        
        # Apply greedy set cover for minimal protein list
        remaining_peptides = copy(component_peptides)
        necessary_proteins = Set{NamedTuple{(:protein_name, :decoy), Tuple{String, Bool}}}()
        
        # Find peptides unique to a protein
        unique_peptide_to_protein = Dictionary{NamedTuple{(:peptide, :decoy), Tuple{String, Bool}}, NamedTuple{(:protein_name, :decoy), Tuple{String, Bool}}}()
        
        for peptide_key in component_peptides
            proteins_for_peptide = intersect(peptide_to_proteins[peptide_key], component_proteins)
            if length(proteins_for_peptide) == 1
                insert!(unique_peptide_to_protein, peptide_key, first(proteins_for_peptide))
            end
        end
        
        # Case F handling: No protein has unique peptides
        if isempty(unique_peptide_to_protein) && !isempty(component_proteins)
            # Group proteins by their decoy status (composite keys already ensure this grouping)
            protein_groups = Dictionary{Bool, Vector{String}}()
            
            for protein in component_proteins
                if !haskey(protein_groups, protein.decoy)
                    insert!(protein_groups, protein.decoy, String[])
                end
                push!(protein_groups[protein.decoy], protein.protein_name)
            end
            
            # Process each group and assign peptides
            for (is_decoy, protein_names) in pairs(protein_groups)
                if !isempty(protein_names)
                    # Get protein names sorted and joined
                    sorted_protein_names = sort(protein_names)
                    protein_group = join(sorted_protein_names, ";")
                    
                    # Assign to peptides with matching decoy status
                    for peptide_key in component_peptides
                        if peptide_key.decoy == is_decoy
                            insert!(result, peptide_key, (
                                protein_name = protein_group,
                                decoy = is_decoy,
                                retain = true
                            ))
                        end
                    end
                end
            end
            continue
        end
        
        # First include proteins with unique peptides
        unique_proteins = Set{NamedTuple{(:protein_name, :decoy), Tuple{String, Bool}}}()
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
            # Find protein that covers the most remaining peptides
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
            
            push!(necessary_proteins, best_protein)
            filter!(p -> p != best_protein, candidate_proteins)
            
            for peptide_key in intersect(protein_to_peptides[best_protein], remaining_peptides)
                delete!(remaining_peptides, peptide_key)
            end
        end
        
        # Create a mapping to track peptides that can be uniquely attributed to a protein in the necessary set
        peptide_to_necessary_protein = Dictionary{NamedTuple{(:peptide, :decoy), Tuple{String, Bool}}, NamedTuple{(:protein_name, :decoy), Tuple{String, Bool}}}()
        
        # Track peptides that can be attributed to multiple necessary proteins
        ambiguous_peptides = Set{NamedTuple{(:peptide, :decoy), Tuple{String, Bool}}}()
        
        for peptide_key in component_peptides
            # Get all necessary proteins that contain this peptide
            proteins_with_peptide = Set{NamedTuple{(:protein_name, :decoy), Tuple{String, Bool}}}()
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
                # Peptide is unique to one necessary protein - assign with retain=true
                protein = peptide_to_necessary_protein[peptide_key]
                insert!(result, peptide_key, (
                    protein_name = protein.protein_name,
                    decoy = protein.decoy,
                    retain = true
                ))
            else
                # Shared peptide - use original protein group and mark as retain=false
                original_group = original_groups[peptide_key]
                insert!(result, peptide_key, (
                    protein_name = original_group.protein_name,
                    decoy = original_group.decoy,
                    retain = false
                ))
            end
        end
    end

    return result
end