# Merge-First Strategy for Greedy Set Cover in Protein Inference

## The Proposal

**Key Idea**: At each iteration of greedy set cover, before selecting the protein with maximum coverage, first identify and merge all proteins with identical peptide sets into protein groups.

### Example Scenario

```
Proteins after Phase 1:
- Protein B: {pep7, pep8}
- Protein C: {pep7, pep8}
- Protein D: {pep9, pep10, pep11}
```

**Current behavior**: Greedy picks B (or C arbitrarily), then picks D
- Result: pep7, pep8 → B; pep9, pep10, pep11 → D
- **Problem**: C is ignored, peptides incorrectly assigned to B alone

**Proposed behavior**:
1. Detect B and C have identical peptide sets
2. Merge into "B;C" group
3. Greedy picks B;C, then picks D
- Result: pep7, pep8 → B;C; pep9, pep10, pep11 → D
- **Correct**: Indistinguishable proteins properly grouped

---

## Algorithm Modification

### Current Greedy Set Cover (Simplified)

```julia
while !isempty(remaining_peptides) && !isempty(candidate_proteins)
    # Find protein with max coverage
    best_protein = argmax(p -> coverage(p, remaining_peptides), candidate_proteins)

    # Add to necessary set
    push!(necessary_proteins, best_protein)
    remove!(candidate_proteins, best_protein)

    # Remove covered peptides
    for pep in covered_by(best_protein)
        delete!(remaining_peptides, pep)
    end
end
```

### Proposed: Merge-First Greedy Set Cover

```julia
while !isempty(remaining_peptides) && !isempty(candidate_proteins)
    # *** NEW STEP 1: Merge indistinguishable proteins ***
    candidate_proteins = merge_indistinguishable_proteins(candidate_proteins, protein_to_peptides)

    # STEP 2: Find protein/group with max coverage
    best_protein = argmax(p -> coverage(p, remaining_peptides), candidate_proteins)

    # STEP 3: Add to necessary set
    push!(necessary_proteins, best_protein)
    remove!(candidate_proteins, best_protein)

    # STEP 4: Remove covered peptides
    for pep in covered_by(best_protein)
        delete!(remaining_peptides, pep)
    end
end

function merge_indistinguishable_proteins(candidates, protein_to_peptides)
    # Group proteins by their peptide sets
    peptide_set_to_proteins = Dict()

    for protein in candidates
        pep_set = protein_to_peptides[protein]
        pep_set_hash = hash(sort(collect(pep_set)))

        if !haskey(peptide_set_to_proteins, pep_set_hash)
            peptide_set_to_proteins[pep_set_hash] = (pep_set, ProteinKey[])
        end
        push!(peptide_set_to_proteins[pep_set_hash][2], protein)
    end

    # Merge groups with >1 protein
    merged_candidates = ProteinKey[]

    for (pep_set, proteins) in values(peptide_set_to_proteins)
        if length(proteins) == 1
            # No merge needed
            push!(merged_candidates, proteins[1])
        else
            # Merge into protein group
            protein_names = sort([p.name for p in proteins])
            merged_protein = ProteinKey(
                join(protein_names, ";"),
                proteins[1].is_target,
                proteins[1].entrap_id
            )
            push!(merged_candidates, merged_protein)

            # Update protein_to_peptides mapping
            protein_to_peptides[merged_protein] = pep_set
        end
    end

    return merged_candidates
end
```

---

## Concrete Example Walkthrough

### Setup

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

**Protein peptide sets:**
- A: {pep1, pep2}
- B: {pep2, pep3, pep4}
- C: {pep2, pep3, pep4}
- D: {pep4, pep5}

### Phase 1: Select Unique Proteins

- A has unique pep1 → Select A, remove pep1, pep2
- D has unique pep5 → Select D, remove pep5, pep4
- Remaining peptides: {pep3}
- Candidate proteins: {B, C}

### Phase 2: Merge-First Greedy Set Cover

**Iteration 1:**

1. **Merge step:**
   - Check B: peptides = {pep3} (after removing pep2, pep4)
   - Check C: peptides = {pep3} (after removing pep2, pep4)
   - **B and C have identical remaining peptides!**
   - Merge: Create "B;C" with peptides {pep3}
   - Updated candidates: {"B;C"}

2. **Select step:**
   - "B;C" covers {pep3} → coverage = 1
   - Select "B;C"
   - Remove pep3

3. **Termination:**
   - Remaining peptides: {} (empty)
   - Done!

**Final necessary proteins:** {A, D, B;C}

**Peptide assignments:**
- pep1 → A (unique to A) → use_for_quant = true
- pep2 → A;B;C (shared) → use_for_quant = false
- pep3 → B;C (unique to B;C group) → use_for_quant = true ✓
- pep4 → D (after A, B;C, D selected, D explains pep4) → use_for_quant = ?
- pep5 → D (unique to D) → use_for_quant = true

---

## Advantages of Merge-First Approach

### 1. **Eliminates Tie-Breaking Ambiguity**
- No need to detect ties during selection
- Indistinguishable proteins are merged before comparison
- Greedy algorithm sees protein groups as single entities

### 2. **Simpler Logic**
- Separation of concerns: merge indistinguishable → greedy select
- No special case handling during selection
- Easier to understand and maintain

### 3. **Correctness by Construction**
- Impossible to select only one protein from an indistinguishable pair
- Peptide assignments automatically correct
- use_for_quant flags naturally correct

### 4. **Matches Parsimony Principle**
- Indistinguishable proteins should be treated as a single entity
- Merging upfront makes this explicit
- Aligns with biological interpretation

### 5. **Efficient**
- Merge check is O(k²) where k = number of candidates
- Only happens at start of each greedy iteration
- Total overhead is small compared to coverage calculations

---

## Implementation Details

### When to Check for Identical Peptide Sets

**Option A: Check full peptide sets**
```julia
# Compare complete peptide sets for each protein
B_peptides = protein_to_peptides[B]  # {pep2, pep3, pep4}
C_peptides = protein_to_peptides[C]  # {pep2, pep3, pep4}
are_identical = (B_peptides == C_peptides)
```

**Option B: Check only remaining peptides (RECOMMENDED)**
```julia
# Compare only peptides still in remaining_peptides
B_remaining = intersect(protein_to_peptides[B], remaining_peptides)  # {pep3}
C_remaining = intersect(protein_to_peptides[C], remaining_peptides)  # {pep3}
are_identical = (B_remaining == C_remaining)
```

**Recommendation**: Use Option B
- More efficient (smaller sets to compare)
- Semantically correct (we only care about coverage of remaining peptides)
- Handles dynamic merging as algorithm progresses

### Handling Target/Decoy and Entrapment Groups

**Important**: Only merge proteins with matching metadata

```julia
function can_merge(protein1::ProteinKey, protein2::ProteinKey)
    return (protein1.is_target == protein2.is_target) &&
           (protein1.entrap_id == protein2.entrap_id)
end

# During merge step:
for (pep_set, proteins) in values(peptide_set_to_proteins)
    if length(proteins) > 1
        # Group by metadata
        metadata_groups = Dict()
        for protein in proteins
            key = (protein.is_target, protein.entrap_id)
            if !haskey(metadata_groups, key)
                metadata_groups[key] = ProteinKey[]
            end
            push!(metadata_groups[key], protein)
        end

        # Merge within each metadata group
        for ((is_target, entrap_id), group) in metadata_groups
            if length(group) > 1
                # Merge this group
                names = sort([p.name for p in group])
                merged = ProteinKey(join(names, ";"), is_target, entrap_id)
                push!(merged_candidates, merged)
            else
                push!(merged_candidates, group[1])
            end
        end
    end
end
```

---

## Potential Issues and Solutions

### Issue 1: Merging Creates New Protein Keys

**Problem**: After merging B and C into "B;C", we need to update mappings

**Solution**:
```julia
# After creating merged protein
merged_protein = ProteinKey("B;C", true, UInt8(1))

# Update protein_to_peptides
protein_to_peptides[merged_protein] = protein_to_peptides[B]  # Same as C

# Remove old entries (optional, for cleanup)
delete!(protein_to_peptides, B)
delete!(protein_to_peptides, C)
```

### Issue 2: What if Merge Happens Mid-Algorithm?

**Scenario**: B and C start with different peptide sets, but after removing peptides, they become identical

**Example:**
- Initially: B = {pep2, pep3}, C = {pep2, pep4}
- After A removes pep2: B = {pep3}, C = {pep4}
- Still different, no merge
- After D removes pep4: B = {pep3}, C = {}
- C has no remaining peptides, removed from candidates

**This is fine!** The merge check at each iteration handles this dynamically.

### Issue 3: Performance Concerns

**Question**: Does checking for identical sets at each iteration slow things down?

**Analysis**:
- Number of candidates decreases each iteration
- Set comparison is O(n) where n = peptides per protein
- Greedy typically runs O(k) iterations where k = number of proteins
- Total: O(k² × n) worst case
- In practice: k is small after Phase 1 (unique proteins already selected)

**Verdict**: Not a bottleneck. The merge check is cheap compared to other operations.

---

## Comparison to Other Approaches

### Approach 1: Detect Ties During Selection (from previous analysis)

**Pros:**
- Only checks for indistinguishability when tie detected
- Minimal changes to greedy logic

**Cons:**
- Tie detection adds complexity
- Must handle merge during selection
- Easy to miss edge cases

### Approach 2: Pre-merge All Indistinguishable Proteins (before greedy)

**Pros:**
- One-time merge before algorithm starts
- Very clean separation

**Cons:**
- Misses proteins that become indistinguishable mid-algorithm
- Example: B = {pep3, pep4}, C = {pep3, pep5} become identical after D removes pep4 and E removes pep5
- Less flexible

### Approach 3: Merge-First at Each Iteration (RECOMMENDED)

**Pros:**
- ✓ Handles dynamic indistinguishability
- ✓ Clean separation of merge and select logic
- ✓ Correct by construction
- ✓ Easy to understand and verify

**Cons:**
- Slight performance overhead (negligible in practice)

---

## Recommendation

**Implement Approach 3: Merge-first at each iteration**

### Algorithm Pseudocode (Updated)

```
# Phase 2: Greedy Set Cover with Merge-First

remaining_peptides = peptides not covered by Phase 1
candidate_proteins = proteins without unique peptides

while remaining_peptides ≠ ∅ and candidate_proteins ≠ ∅:

    # Step 1: Merge indistinguishable proteins
    grouped_proteins = group by (remaining peptide set, is_target, entrap_id)
    merged_candidates = []

    for each group in grouped_proteins:
        if |group| = 1:
            add single protein to merged_candidates
        else:
            # Merge proteins in group
            protein_names = sort([p.name for p in group])
            merged_protein = ProteinKey(join(protein_names, ";"), ...)
            add merged_protein to merged_candidates
            update protein_to_peptides[merged_protein]

    candidate_proteins = merged_candidates

    # Step 2: Find protein with max coverage (standard greedy)
    best_protein = argmax(p -> |protein_to_peptides[p] ∩ remaining_peptides|)

    if coverage(best_protein) = 0:
        break

    # Step 3: Select protein
    add best_protein to necessary_proteins
    remove best_protein from candidate_proteins

    # Step 4: Remove covered peptides
    remaining_peptides -= protein_to_peptides[best_protein]

return necessary_proteins
```

---

## Test Cases to Validate

### Test 1: Basic Indistinguishable Case
```julia
# B and C identical after Phase 1
proteins = [ProteinKey("A"), ProteinKey("B;C"), ProteinKey("B;C"), ProteinKey("D")]
peptides = [pep1, pep2, pep3, pep4]
# Expect: pep2 and pep3 both assigned to "B;C"
```

### Test 2: Dynamic Indistinguishability
```julia
# B and C become identical mid-algorithm
proteins = [ProteinKey("A"), ProteinKey("B"), ProteinKey("C"), ProteinKey("D")]
peptides = [pep1, pep2, pep3, pep4]
# Initially: B = {pep2, pep3}, C = {pep2, pep4}
# After removing pep2: B = {pep3}, C = {pep4} (still different)
# No merge needed in this case
```

### Test 3: Multiple Indistinguishable Groups
```julia
# Both (B,C) and (D,E) are indistinguishable pairs
proteins = [ProteinKey("A"), ProteinKey("B;C"), ProteinKey("B;C"),
            ProteinKey("D;E"), ProteinKey("D;E")]
peptides = [pep1, pep2, pep3, pep4, pep5]
# Expect: two merged groups "B;C" and "D;E"
```

---

## LaTeX Algorithm Update

The algorithm in `protein_inference.tex` should be updated:

```latex
\State \textit{// Phase 2: Greedy set cover with merge-first (lines 312-338)}
\State $K \gets A_c \setminus S$ \Comment{Candidate proteins}
\While{$R \neq \emptyset$ and $K \neq \emptyset$}
    \State \textit{// Merge indistinguishable proteins}
    \State Group $K$ by $(P[a] \cap R, \text{is\_target}, \text{entrap\_id})$
    \State $K' \gets \emptyset$
    \For{each group $G$ with $|G| > 1$}
        \State Create merged protein $a^* = \text{join}(\text{sort}(G), ";")$
        \State $P[a^*] \gets P[a]$ for any $a \in G$
        \State $K' \gets K' \cup \{a^*\}$
    \EndFor
    \For{each group $G$ with $|G| = 1$}
        \State $K' \gets K' \cup G$
    \EndFor
    \State $K \gets K'$
    \State
    \State \textit{// Standard greedy selection}
    \State $a^* \gets \arg\max_{a \in K} |P[a] \cap R|$
    \If{$|P[a^*] \cap R| = 0$}
        \State \textbf{break}
    \EndIf
    \State $S \gets S \cup \{a^*\}$
    \State $R \gets R \setminus P[a^*]$
    \State $K \gets K \setminus \{a^*\}$
\EndWhile
```

---

## Conclusion

The **merge-first approach** is the cleanest and most correct solution to the tie-breaking problem in protein inference:

1. ✅ Eliminates arbitrary tie-breaking
2. ✅ Correctly handles indistinguishable proteins
3. ✅ Works dynamically as algorithm progresses
4. ✅ Simple to understand and implement
5. ✅ Aligns with biological interpretation
6. ✅ Negligible performance overhead

This should be the approach implemented in `src/utils/proteinInference.jl`.
