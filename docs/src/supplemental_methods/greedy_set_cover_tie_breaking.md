# Greedy Set Cover Tie-Breaking Bug in Protein Inference

## The Problem

**Issue**: When multiple proteins have identical coverage (tie for `best_coverage`), the current greedy set cover algorithm arbitrarily selects the first one encountered. This causes incorrect protein group assignments when proteins are indistinguishable.

### Concrete Example

```julia
proteins = [
    ProteinKey("A", true, UInt8(1)),        # pep1
    ProteinKey("A;B;C", true, UInt8(1)),    # pep2
    ProteinKey("B;C", true, UInt8(1)),      # pep3
    ProteinKey("B;C;D", true, UInt8(1)),    # pep4
    ProteinKey("D", true, UInt8(1)),        # pep5
]

peptides = [
    PeptideKey("1", true, UInt8(1)),
    PeptideKey("2", true, UInt8(1)),
    PeptideKey("3", true, UInt8(1)),
    PeptideKey("4", true, UInt8(1)),
    PeptideKey("5", true, UInt8(1)),
]
```

**Peptide-Protein Mapping:**
- Peptide 1 → {A}
- Peptide 2 → {A, B, C}
- Peptide 3 → {B, C}
- Peptide 4 → {B, C, D}
- Peptide 5 → {D}

**Algorithm Execution:**

1. **Phase 1: Select proteins with unique peptides**
   - Protein A has unique peptide 1 → Select A, remove pep1, pep2
   - Protein D has unique peptide 5 → Select D, remove pep5, pep4
   - Remaining peptides: {pep3}

2. **Phase 2: Greedy set cover for pep3**
   - Candidate proteins: {B, C}
   - Protein B covers: {pep3} → coverage = 1
   - Protein C covers: {pep3} → coverage = 1
   - **TIE!** But algorithm picks B arbitrarily (first in iteration order)
   - Result: pep3 → B (WRONG!)

**Expected Behavior:**
- Proteins B and C are **indistinguishable** (identical peptide sets)
- Peptide 3 should map to **"B;C"** protein group
- `use_for_quant` should be **true** (the group B;C is unambiguous)

**Actual Behavior:**
- Peptide 3 maps to **"B"**
- Protein C is ignored
- This violates the parsimony principle for indistinguishable proteins

---

## Root Cause Analysis

### Location in Code
File: `src/utils/proteinInference.jl`, lines 318-341

```julia
while !isempty(remaining_peptides) && !isempty(candidate_proteins)
    # Find protein that covers the most remaining peptides
    best_protein = nothing
    best_coverage = 0

    for protein in candidate_proteins
        coverage = length(intersect(protein_to_peptides[protein], remaining_peptides))
        if coverage > best_coverage        # <-- ISSUE: Only checks >, not ==
            best_coverage = coverage
            best_protein = protein
        end
    end

    if best_coverage == 0
        break
    end

    push!(necessary_proteins, best_protein)
    filter!(p -> p != best_protein, candidate_proteins)

    for peptide_key in intersect(protein_to_peptides[best_protein], remaining_peptides)
        delete!(remaining_peptides, peptide_key)
    end
end
```

### Why This Happens

1. When two proteins (B and C) have identical peptide sets, they both cover the same remaining peptides
2. The algorithm only checks `if coverage > best_coverage`, so ties are not detected
3. The first protein in iteration order wins the tie
4. The second protein is never selected, even though it's indistinguishable
5. Peptides get assigned to the single protein instead of the group

---

## Proposed Solutions

### Option 1: Detect Indistinguishable Proteins in Tie Situations

**Idea**: When there's a tie, check if the tied proteins have identical peptide sets. If so, merge them into a group.

**Implementation:**
```julia
while !isempty(remaining_peptides) && !isempty(candidate_proteins)
    # Find ALL proteins that cover the most remaining peptides
    best_coverage = 0
    tied_proteins = Set{ProteinKey}()

    for protein in candidate_proteins
        coverage = length(intersect(protein_to_peptides[protein], remaining_peptides))
        if coverage > best_coverage
            best_coverage = coverage
            tied_proteins = Set([protein])
        elseif coverage == best_coverage && coverage > 0
            push!(tied_proteins, protein)
        end
    end

    if best_coverage == 0
        break
    end

    # Check if tied proteins are indistinguishable (identical peptide sets)
    if length(tied_proteins) > 1
        indistinguishable = true
        first_protein = first(tied_proteins)
        first_peptides = protein_to_peptides[first_protein]

        for protein in tied_proteins
            if protein_to_peptides[protein] != first_peptides
                indistinguishable = false
                break
            end
        end

        if indistinguishable
            # Merge into a protein group
            protein_names = sort([p.name for p in tied_proteins])
            merged_protein = ProteinKey(
                join(protein_names, ";"),
                first_protein.is_target,
                first_protein.entrap_id
            )
            push!(necessary_proteins, merged_protein)

            # Remove all tied proteins from candidates
            filter!(p -> !(p in tied_proteins), candidate_proteins)

            # Update protein_to_peptides with merged group
            insert!(protein_to_peptides, merged_protein, copy(first_peptides))
        else
            # Tie but not indistinguishable - pick arbitrarily (first one)
            best_protein = first(tied_proteins)
            push!(necessary_proteins, best_protein)
            filter!(p -> p != best_protein, candidate_proteins)
        end
    else
        # No tie - add the single best protein
        best_protein = first(tied_proteins)
        push!(necessary_proteins, best_protein)
        filter!(p -> p != best_protein, candidate_proteins)
    end

    # Remove covered peptides
    for protein in (length(tied_proteins) > 1 && indistinguishable) ? [merged_protein] : [best_protein]
        for peptide_key in intersect(protein_to_peptides[protein], remaining_peptides)
            delete!(remaining_peptides, peptide_key)
        end
    end
end
```

**Pros:**
- Correctly handles indistinguishable proteins
- Maintains parsimony principle
- Ties are resolved by checking for identical peptide sets

**Cons:**
- More complex logic
- Need to handle merged protein groups carefully in subsequent steps
- Need to update `protein_to_peptides` mapping with merged groups

---

### Option 2: Pre-process Component to Merge Indistinguishable Proteins

**Idea**: Before running greedy set cover, scan the component and merge all indistinguishable proteins upfront.

**Implementation:**
```julia
# After finding connected components, before greedy set cover
# Merge indistinguishable proteins within each component

for component in components
    component_proteins = component[2]

    # Group proteins by their peptide sets
    peptide_set_to_proteins = Dict()
    for protein in component_proteins
        pep_set = protein_to_peptides[protein]
        pep_set_key = hash(sort(collect(pep_set)))  # Use hash as key

        if !haskey(peptide_set_to_proteins, pep_set_key)
            peptide_set_to_proteins[pep_set_key] = (pep_set, Set{ProteinKey}())
        end
        push!(peptide_set_to_proteins[pep_set_key][2], protein)
    end

    # For each group with >1 protein, merge them
    for (pep_set, protein_set) in values(peptide_set_to_proteins)
        if length(protein_set) > 1
            # Merge proteins with identical peptide sets
            protein_names = sort([p.name for p in protein_set])
            merged_protein = ProteinKey(
                join(protein_names, ";"),
                first(protein_set).is_target,
                first(protein_set).entrap_id
            )

            # Update mappings
            insert!(protein_to_peptides, merged_protein, pep_set)

            # Remove individual proteins
            for protein in protein_set
                delete!(protein_to_peptides, protein)
                # Update component_proteins
                filter!(p -> p != protein, component_proteins)
            end

            # Add merged protein to component
            push!(component_proteins, merged_protein)
        end
    end
end
```

**Pros:**
- Cleaner separation of concerns
- Greedy set cover remains simple
- All indistinguishable proteins handled upfront

**Cons:**
- More preprocessing overhead
- Need to carefully update all mappings
- Changes component structure before main algorithm

---

### Option 3: Post-process Assignments to Detect Singleton Group Members

**Idea**: After greedy set cover, check if any selected proteins should have been grouped together.

**Implementation:**
```julia
# After greedy set cover completes
# Check if any necessary proteins have identical peptide sets

protein_set_groups = Dict()
for protein in necessary_proteins
    pep_set = intersect(protein_to_peptides[protein], component_peptides)
    pep_set_key = hash(sort(collect(pep_set)))

    if !haskey(protein_set_groups, pep_set_key)
        protein_set_groups[pep_set_key] = Set{ProteinKey}()
    end
    push!(protein_set_groups[pep_set_key], protein)
end

# Merge proteins with identical peptide sets
for protein_set in values(protein_set_groups)
    if length(protein_set) > 1
        # These proteins are indistinguishable
        protein_names = sort([p.name for p in protein_set])
        merged_protein = ProteinKey(
            join(protein_names, ";"),
            first(protein_set).is_target,
            first(protein_set).entrap_id
        )

        # Replace all instances in necessary_proteins
        for protein in protein_set
            delete!(necessary_proteins, protein)
        end
        push!(necessary_proteins, merged_protein)

        # Update protein_to_peptides
        insert!(protein_to_peptides, merged_protein,
                protein_to_peptides[first(protein_set)])
    end
end
```

**Pros:**
- Minimal changes to existing algorithm
- Handles the issue after the fact
- Simpler to implement

**Cons:**
- Reactive rather than proactive
- May miss some edge cases
- Still requires careful mapping updates

---

## Recommended Solution

**I recommend Option 1: Detect and handle ties during greedy set cover**

### Rationale:

1. **Correctness at the source**: Handles the tie-breaking problem exactly where it occurs
2. **Maintains algorithm semantics**: The greedy set cover naturally identifies when proteins are equivalent
3. **Aligns with LaTeX description**: The algorithm description in `protein_inference.tex` should specify this tie-breaking behavior
4. **Minimal side effects**: Doesn't require restructuring components or extensive post-processing

### Additional Considerations:

1. **Peptide assignment after merging**: When a merged protein group is created during greedy set cover, we need to ensure:
   - The merged group is added to `necessary_proteins`
   - The `protein_to_peptides` mapping includes the merged group
   - Peptides assigned to the merged group get `use_for_quant = true` (since the group is unambiguous)

2. **Edge case: Ties with different peptide sets**
   - If two proteins tie on coverage but have different peptide sets, we should pick deterministically
   - Options: alphabetical order, smallest protein name, etc.
   - Current arbitrary selection is actually fine for this case

3. **Updating the algorithm description**
   - The `protein_inference.tex` should be updated to specify tie-breaking rules
   - Lines 47-55 of the algorithm should note: "If multiple proteins have equal coverage, check if they are indistinguishable (identical peptide sets). If so, merge into a protein group."

---

## Test Case to Validate Fix

```julia
# This should pass after fix
proteins = [
    ProteinKey("A", true, UInt8(1)),
    ProteinKey("A;B;C", true, UInt8(1)),
    ProteinKey("B;C", true, UInt8(1)),
    ProteinKey("B;C;D", true, UInt8(1)),
    ProteinKey("D", true, UInt8(1)),
]

peptides = [
    PeptideKey("1", true, UInt8(1)),
    PeptideKey("2", true, UInt8(1)),
    PeptideKey("3", true, UInt8(1)),
    PeptideKey("4", true, UInt8(1)),
    PeptideKey("5", true, UInt8(1)),
]

result = infer_proteins(proteins, peptides)

# Expected results:
@test result.peptide_to_protein[peptides[1]].name == "A"
@test result.peptide_to_protein[peptides[2]].name == "A"  # or "A;B;C" - need to decide
@test result.peptide_to_protein[peptides[3]].name == "B;C"  # KEY: Should be "B;C", not "B"
@test result.peptide_to_protein[peptides[4]].name == "D"  # or "B;C;D" - need to decide
@test result.peptide_to_protein[peptides[5]].name == "D"

@test result.use_for_quant[peptides[1]] == true
@test result.use_for_quant[peptides[3]] == true  # KEY: Should be true for "B;C" group
@test result.use_for_quant[peptides[5]] == true
```

---

## LaTeX Algorithm Update

The algorithm in `protein_inference.tex` (lines 318-341) should be updated to include tie-breaking logic:

```latex
\State \textit{// Phase 2: Greedy set cover for remaining peptides (lines 312-338)}
\State $K \gets A_c \setminus S$ \Comment{Candidate proteins}
\While{$R \neq \emptyset$ and $K \neq \emptyset$}
    \State $T \gets \{a \in K : |P[a] \cap R| = \max_{a' \in K} |P[a'] \cap R|\}$ \Comment{Proteins with max coverage}
    \If{$\max_{a \in K} |P[a] \cap R| = 0$}
        \State \textbf{break} \Comment{No more coverage possible}
    \EndIf
    \If{$|T| > 1$} \Comment{Tie detected}
        \State Check if $P[a] = P[a']$ for all $a, a' \in T$ \Comment{Indistinguishable?}
        \If{yes}
            \State $a^* \gets$ merged protein group from $T$ \Comment{Merge tied proteins}
            \State $P[a^*] \gets P[a]$ for any $a \in T$
        \Else
            \State $a^* \gets$ arbitrary choice from $T$ \Comment{Break tie deterministically}
        \EndIf
    \Else
        \State $a^* \gets$ the single protein in $T$
    \EndIf
    \State $S \gets S \cup \{a^*\}$ \Comment{Add to necessary proteins}
    \State $R \gets R \setminus P[a^*]$ \Comment{Remove covered peptides}
    \State $K \gets K \setminus T$ \Comment{Remove all tied proteins from candidates}
\EndWhile
```
