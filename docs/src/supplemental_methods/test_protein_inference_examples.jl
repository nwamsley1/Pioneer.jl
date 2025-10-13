# Test Protein Inference Examples
#
# This file demonstrates the protein inference algorithm described in protein_inference.tex
# using the Pioneer.jl implementation. It provides concrete examples for each case
# described in the supplemental methods.

using Pkg
Pkg.develop(".")  # Activate the Pioneer.jl environment

using Pioneer
using Pioneer: infer_proteins, ProteinKey, PeptideKey, InferenceResult
using Dictionaries

println("="^80)
println("Protein Inference Algorithm Examples")
println("Based on the two-phase parsimony algorithm in protein_inference.tex")
println("="^80)
println()

# Helper function to display results
function display_inference_result(result::InferenceResult, peptides::Vector{PeptideKey})
    println("\nInference Results:")
    println("-"^60)
    for pep in peptides
        prot = result.peptide_to_protein[pep]
        quant = result.use_for_quant[pep]
        println("  Peptide: $(pep.sequence) → Protein: $(prot.name)")
        println("    Use for quant: $(quant)")
    end
    println("-"^60)
end

#=============================================================================
Case 1: Distinct Proteins (Separate Connected Components)
==============================================================================#
println("\n" * "="^80)
println("Case 1: Distinct Proteins")
println("="^80)
println("""
Description:
  - Protein A has peptides 1, 2 (unique to A)
  - Protein B has peptides 3, 4 (unique to B)
  - These form two separate connected components in the bipartite graph

Expected Result:
  - Each peptide maps uniquely to its protein
  - All peptides used for quantification
""")

proteins_case1 = [
    ProteinKey("A", true, UInt8(1)),
    ProteinKey("A;B;C", true, UInt8(1)),
    ProteinKey("B;C", true, UInt8(1)),
    ProteinKey("B;C;D", true, UInt8(1)),
    ProteinKey("D", true, UInt8(1)),
];

peptides_case1 = [
    PeptideKey("1", true, UInt8(1)),
    PeptideKey("2", true, UInt8(1)),
    PeptideKey("3", true, UInt8(1)),
    PeptideKey("4", true, UInt8(1)),
    PeptideKey("5", true, UInt8(1)),
];

result_case1 = infer_proteins(proteins_case1, peptides_case1)
display_inference_result(result_case1, peptides_case1)


#=============================================================================
Case 1: Distinct Proteins (Separate Connected Components)
==============================================================================#
println("\n" * "="^80)
println("Case 1: Distinct Proteins")
println("="^80)
println("""
Description:
  - Protein A has peptides 1, 2 (unique to A)
  - Protein B has peptides 3, 4 (unique to B)
  - These form two separate connected components in the bipartite graph

Expected Result:
  - Each peptide maps uniquely to its protein
  - All peptides used for quantification
""")

proteins = [
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1))
]
peptides = [
    PeptideKey("pep1", true, UInt8(1)),
    PeptideKey("pep2", true, UInt8(1)),
    PeptideKey("pep3", true, UInt8(1)),
    PeptideKey("pep4", true, UInt8(1))
]
        

result = infer_proteins(proteins, peptides)
display_inference_result(result, peptides)



proteins = [
    ProteinKey("A", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("A;C", true, UInt8(1)),
    ProteinKey("B", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1))
]
peptides = [
    PeptideKey("1", true, UInt8(1)),
    PeptideKey("2", true, UInt8(1)),
    PeptideKey("3", true, UInt8(1)),
    PeptideKey("4", true, UInt8(1)),
    PeptideKey("5", true, UInt8(1)),
    PeptideKey("6", true, UInt8(1)),
    PeptideKey("7", true, UInt8(1)),
]
        

result = infer_proteins(proteins, peptides)
display_inference_result(result, peptides)


proteins = [
    ProteinKey("A", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("B;C", true, UInt8(1)),
    ProteinKey("C;D", true, UInt8(1)),
    ProteinKey("D", true, UInt8(1)),
    ProteinKey("D;E", true, UInt8(1)),
    ProteinKey("F;E", true, UInt8(1)),
    ProteinKey("A;E", true, UInt8(1))
]
peptides = [
    PeptideKey("1", true, UInt8(1)),
    PeptideKey("2", true, UInt8(1)),
    PeptideKey("3", true, UInt8(1)),
    PeptideKey("4", true, UInt8(1)),
    PeptideKey("5", true, UInt8(1)),
    PeptideKey("6", true, UInt8(1)),
    PeptideKey("7", true, UInt8(1)),
    PeptideKey("8", true, UInt8(1))
]
        

result = infer_proteins(proteins, peptides)
display_inference_result(result, peptides)


#=============================================================================
Case 2: Differentiable Proteins (Unique + Shared Peptides)
==============================================================================#
Pioneer.PeptideKey[
  Pioneer.PeptideKey("KESYSIYVYK", true, 0x00), 
  Pioneer.PeptideKey("ESYSIYVYK", true, 0x00), 
  Pioneer.PeptideKey("ETYSSYIYK", true, 0x00), 
  Pioneer.PeptideKey("KESYSVYVYK", true, 0x00), 
  Pioneer.PeptideKey("EIQTAVR", true, 0x00), 
  Pioneer.PeptideKey("QVHPDTGISSK", true, 0x00), 
  Pioneer.PeptideKey("ESYSIYIYK", true, 0x00), 
  Pioneer.PeptideKey("AMGIMNSFVNDIFER", true, 0x00), 
  Pioneer.PeptideKey("ESYSVYVYK", true, 0x00), 
  Pioneer.PeptideKey("SMSILNSFVNDIFER", true, 0x00), 
  Pioneer.PeptideKey("AMSIMNSFVTDIFER", true, 0x00)
  ]

Pioneer.ProteinKey[
  Pioneer.ProteinKey("O60814", true, 0x00), 
  Pioneer.ProteinKey("P33778", true, 0x00), 
  Pioneer.ProteinKey("P02293", true, 0x00), 
  Pioneer.ProteinKey("Q6DRA6", true, 0x00), 
  Pioneer.ProteinKey("Q93079", true, 0x00), 
  Pioneer.ProteinKey("Q99880", true, 0x00), 
  Pioneer.ProteinKey("P58876", true, 0x00), 
  Pioneer.ProteinKey("P23527", true, 0x00), 
  Pioneer.ProteinKey("Q99879", true, 0x00), 
  Pioneer.ProteinKey("Q5QNW6", true, 0x00), 
  Pioneer.ProteinKey("P62807", true, 0x00), 
  Pioneer.ProteinKey("P02294", true, 0x00), 
  Pioneer.ProteinKey("P06899", true, 0x00), 
  Pioneer.ProteinKey("Q99877", true, 0x00), 
  Pioneer.ProteinKey("Q96A08", true, 0x00), 
  Pioneer.ProteinKey("Q16778", true, 0x00), 
  Pioneer.ProteinKey("Q6DN03", true, 0x00), 
  Pioneer.ProteinKey("Q8N257", true, 0x00), 
  Pioneer.ProteinKey("P57053", true, 0x00)]

#=============================================================================
Case 2: Differentiable Proteins (Unique + Shared Peptides)
==============================================================================#
println("\n" * "="^80)
println("Case 2: Differentiable Proteins")
println("="^80)
println("""
Description:
  - Protein A has unique peptide 1
  - Proteins A and B share peptides 2, 3
  - Protein B has unique peptide 4
  - Single connected component, but proteins distinguishable by unique peptides

Algorithm Steps:
  1. Identify unique peptides: pep1 (A only), pep4 (B only)
  2. Phase 1: Select proteins A and B (both have unique peptides)
  3. Phase 2: Greedy set cover not needed (all peptides covered)
  4. Assign quantification flags: unique peptides = true, shared = false

Expected Result:
  - Unique peptides map to single proteins (use_for_quant = true)
  - Shared peptides map to "A;B" group (use_for_quant = false)
""")

proteins_case2 = [
    ProteinKey("A", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("B", true, UInt8(1))
]

peptides_case2 = [
    PeptideKey("PEPTIDE_A_UNIQUE", true, UInt8(1)),
    PeptideKey("SHARED_AB_1", true, UInt8(1)),
    PeptideKey("SHARED_AB_2", true, UInt8(1)),
    PeptideKey("PEPTIDE_B_UNIQUE", true, UInt8(1))
]

result_case2 = infer_proteins(proteins_case2, peptides_case2)
display_inference_result(result_case2, peptides_case2)

#=============================================================================
Case 3: Indistinguishable Proteins (Identical Peptide Sets)
==============================================================================#
println("\n" * "="^80)
println("Case 3: Indistinguishable Proteins")
println("="^80)
println("""
Description:
  - Proteins A and B have identical peptide sets
  - Cannot distinguish between them based on observed evidence
  - Algorithm detects: P[A] = P[B] for all proteins in component

Algorithm Steps:
  1. Check if all proteins have identical peptide sets
  2. If yes, assign all peptides to combined group "A;B"
  3. All peptides flagged for quantification (group is unambiguous)

Expected Result:
  - All peptides → "A;B"
  - All peptides use_for_quant = true (protein group is the minimal explanation)
""")

proteins_case3 = [
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1))
]

peptides_case3 = [
    PeptideKey("PEPTIDE_AB1", true, UInt8(1)),
    PeptideKey("PEPTIDE_AB2", true, UInt8(1)),
    PeptideKey("PEPTIDE_AB3", true, UInt8(1)),
    PeptideKey("PEPTIDE_AB4", true, UInt8(1))
]

result_case3 = infer_proteins(proteins_case3, peptides_case3)
display_inference_result(result_case3, peptides_case3)

#=============================================================================
Case 4: Subset Proteins (Parsimony Principle)
==============================================================================#
println("\n" * "="^80)
println("Case 4: Subset Proteins (Parsimony)")
println("="^80)
println("""
Description:
  - Protein A has unique peptide 1
  - Peptides 2, 3, 4 are shared by both A and B
  - Protein A alone can explain all peptides (parsimony)

Algorithm Steps:
  1. Identify unique peptides: pep1 (A only)
  2. Phase 1: Select protein A (has unique peptide)
  3. Check remaining peptides: all covered by A
  4. Phase 2: Not needed, B is redundant
  5. Assign all peptides to A with use_for_quant = true

Expected Result:
  - All peptides → "A" (minimal explanation)
  - All peptides use_for_quant = true
""")

proteins_case4 = [
    ProteinKey("A", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1))
]

peptides_case4 = [
    PeptideKey("PEPTIDE_A_UNIQUE", true, UInt8(1)),
    PeptideKey("SHARED_AB_1", true, UInt8(1)),
    PeptideKey("SHARED_AB_2", true, UInt8(1)),
    PeptideKey("SHARED_AB_3", true, UInt8(1))
]

result_case4 = infer_proteins(proteins_case4, peptides_case4)
display_inference_result(result_case4, peptides_case4)

#=============================================================================
Case 5: Subsumable Proteins (Greedy Set Cover)
==============================================================================#
println("\n" * "="^80)
println("Case 5: Subsumable Proteins (Greedy Set Cover)")
println("="^80)
println("""
Description:
  - Protein A has unique peptide 1, shares pep2 with B
  - Protein B shares pep2 with A, shares pep3 with C
  - Protein C has unique peptide 4, shares pep3 with B
  - Need both A and C to explain all peptides (B is redundant)

Algorithm Steps:
  1. Identify unique peptides: pep1 (A), pep4 (C)
  2. Phase 1: Select A and C (both have unique peptides)
  3. Check coverage: A covers {pep1, pep2}, C covers {pep3, pep4}
  4. All peptides covered, B not needed (subsumed by A and C)
  5. Assign peptides to proteins that can uniquely explain them

Expected Result:
  - pep1, pep2 → "A"
  - pep3, pep4 → "C"
  - All peptides use_for_quant = true (each uniquely assigned)
""")

proteins_case5 = [
    ProteinKey("A", true, UInt8(1)),
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("B;C", true, UInt8(1)),
    ProteinKey("C", true, UInt8(1))
]

peptides_case5 = [
    PeptideKey("PEPTIDE_A_UNIQUE", true, UInt8(1)),
    PeptideKey("SHARED_AB", true, UInt8(1)),
    PeptideKey("SHARED_BC", true, UInt8(1)),
    PeptideKey("PEPTIDE_C_UNIQUE", true, UInt8(1))
]

result_case5 = infer_proteins(proteins_case5, peptides_case5)
display_inference_result(result_case5, peptides_case5)

#=============================================================================
Case 6: No Unique Peptides (Greedy Set Cover Required)
==============================================================================#
println("\n" * "="^80)
println("Case 6: No Unique Peptides (Greedy Set Cover)")
println("="^80)
println("""
Description:
  - No protein has a unique peptide
  - All peptides are shared among multiple proteins
  - Must use greedy set cover in Phase 2

Algorithm Steps:
  1. Identify unique peptides: none found (U = ∅)
  2. Phase 1: No proteins selected (S = ∅)
  3. Phase 2: Greedy set cover
     - Find protein covering most remaining peptides
     - Add to S, update remaining peptides R
     - Repeat until all covered
  4. Assign peptides based on final selected proteins

Expected Result:
  - Minimal set of proteins selected
  - Peptides mapped to protein groups
  - use_for_quant depends on whether peptide maps to single protein
""")

proteins_case6 = [
    ProteinKey("A;B", true, UInt8(1)),
    ProteinKey("A;B;C", true, UInt8(1)),
    ProteinKey("A;B;C", true, UInt8(1)),
    ProteinKey("A;C", true, UInt8(1))
]

peptides_case6 = [
    PeptideKey("SHARED_AB", true, UInt8(1)),
    PeptideKey("SHARED_ABC_1", true, UInt8(1)),
    PeptideKey("SHARED_ABC_2", true, UInt8(1)),
    PeptideKey("SHARED_AC", true, UInt8(1))
]

result_case6 = infer_proteins(proteins_case6, peptides_case6)
display_inference_result(result_case6, peptides_case6)

#=============================================================================
Case 7: Target/Decoy Separation
==============================================================================#
println("\n" * "="^80)
println("Case 7: Target/Decoy Separation")
println("="^80)
println("""
Description:
  - Target and decoy proteins are processed separately
  - Target peptides only map to target proteins
  - Decoy peptides only map to decoy proteins
  - Maintains target-decoy separation for FDR calculation

Expected Result:
  - Target peptide → target protein
  - Decoy peptide → decoy protein
  - is_target flag preserved in assignments
""")

proteins_case7 = [
    ProteinKey("A", true, UInt8(1)),    # Target
    ProteinKey("B", false, UInt8(1))    # Decoy
]

peptides_case7 = [
    PeptideKey("TARGET_PEPTIDE", true, UInt8(1)),
    PeptideKey("DECOY_PEPTIDE", false, UInt8(1))
]

result_case7 = infer_proteins(proteins_case7, peptides_case7)
display_inference_result(result_case7, peptides_case7)

#=============================================================================
Case 8: Entrapment Group Separation
==============================================================================#
println("\n" * "="^80)
println("Case 8: Entrapment Group Separation")
println("="^80)
println("""
Description:
  - Same peptide sequence in different entrapment groups
  - Treated as distinct peptides for inference
  - Each maps to protein in matching entrapment group
  - Used for FDR calibration (Wen et al. 2024)

Expected Result:
  - Each peptide maps to protein with matching entrap_id
  - Entrapment separation maintained
""")

proteins_case8 = [
    ProteinKey("A", true, UInt8(1)),   # Entrapment group 1
    ProteinKey("A", true, UInt8(2))    # Entrapment group 2
]

peptides_case8 = [
    PeptideKey("SAME_SEQUENCE", true, UInt8(1)),
    PeptideKey("SAME_SEQUENCE", true, UInt8(2))
]

result_case8 = infer_proteins(proteins_case8, peptides_case8)
display_inference_result(result_case8, peptides_case8)

#=============================================================================
Case 9: Complex Multi-Component Graph
==============================================================================#
println("\n" * "="^80)
println("Case 9: Complex Multi-Component Graph")
println("="^80)
println("""
Description:
  - Multiple disconnected components in bipartite graph
  - Each component processed independently
  - Demonstrates depth-first search (DFS) for component discovery
  - Shows V (visited set) update after each component

Graph Structure:
  Component 1: A-pep1, A-pep2, B-pep3, B-pep4 (distinct proteins)
  Component 2: C-pep5, C;D-pep6, C;D-pep7, D-pep8 (differentiable)
  Component 3: E;F-pep9, E;F-pep10 (indistinguishable)

Algorithm Steps (per component):
  1. DFS discovers component (P_c, A_c)
  2. V ← V ∪ P_c (mark peptides as visited)
  3. Process component with two-phase algorithm
  4. Move to next unvisited peptide

Expected Result:
  - Component 1: All unique assignments
  - Component 2: Unique to proteins, shared peptides flagged
  - Component 3: All to "E;F" group
""")

proteins_case9 = [
    # Component 1: Distinct proteins
    ProteinKey("A", true, UInt8(1)),
    ProteinKey("A", true, UInt8(1)),
    ProteinKey("B", true, UInt8(1)),
    ProteinKey("B", true, UInt8(1)),

    # Component 2: Differentiable proteins
    ProteinKey("C", true, UInt8(1)),
    ProteinKey("C;D", true, UInt8(1)),
    ProteinKey("C;D", true, UInt8(1)),
    ProteinKey("D", true, UInt8(1)),

    # Component 3: Indistinguishable proteins
    ProteinKey("E;F", true, UInt8(1)),
    ProteinKey("E;F", true, UInt8(1))
]

peptides_case9 = [
    PeptideKey("A_PEP1", true, UInt8(1)),
    PeptideKey("A_PEP2", true, UInt8(1)),
    PeptideKey("B_PEP3", true, UInt8(1)),
    PeptideKey("B_PEP4", true, UInt8(1)),
    PeptideKey("C_UNIQUE", true, UInt8(1)),
    PeptideKey("CD_SHARED1", true, UInt8(1)),
    PeptideKey("CD_SHARED2", true, UInt8(1)),
    PeptideKey("D_UNIQUE", true, UInt8(1)),
    PeptideKey("EF_PEP1", true, UInt8(1)),
    PeptideKey("EF_PEP2", true, UInt8(1))
]

result_case9 = infer_proteins(proteins_case9, peptides_case9)
display_inference_result(result_case9, peptides_case9)

#=============================================================================
Summary
==============================================================================#
println("\n" * "="^80)
println("Summary")
println("="^80)
println("""
These examples demonstrate the two-phase parsimony-based protein inference
algorithm described in protein_inference.tex:

Phase 1: Automatic Selection
  - Select all proteins with unique peptide evidence
  - Satisfies parsimony principle (never exclude unique evidence)

Phase 2: Greedy Set Cover
  - For remaining uncovered peptides
  - Select protein covering most peptides
  - Repeat until all peptides covered

Key Properties:
  1. DFS decomposes graph into connected components
  2. Components processed independently
  3. Visited set (V) prevents reprocessing
  4. Target/decoy separation maintained
  5. Entrapment groups handled separately
  6. Quantification flags based on unique assignment

Corresponds to:
  - Algorithm 1 in protein_inference.tex (lines 1-69)
  - Implementation in src/utils/proteinInference.jl (lines 107-380)
""")

#=============================================================================
Case: Real-World Integration Test Failure
==============================================================================#
println("\n" * "="^80)
println("Real-World Case: Integration Test Failure (2025-10-13)")
println("="^80)
println("""
Description:
  - 19 proteins in a single connected component (simplified as A-S)
  - 11 peptides with complex sharing patterns (p1-p11)
  - Multiple proteins share identical peptide sets (requiring merge-first)
  - Tests the sort() operation that triggered the original isless error

Component Structure (Bipartite Graph):
  Group 1 (9 proteins, 5 peptides): A,C,D,E,F,G,H,I,K share {p1,p2,p3,p4,p5}
  Group 2 (4 proteins, 5 peptides): B,L,M,N share {p6,p2,p7,p4,p5}
  Group 3 (2 proteins, 3 peptides): O,P (histones) share {p8,p9,p5}
  Group 4 (2 proteins, 2 peptides): Q,R share {p6,p7}
  Group 5 (1 protein, 4 peptides): S has {p2,p10,p11,p5}
  Group 6 (1 protein, 4 peptides): J has unique pattern {p2,p12,p13,p5}

Key Observations:
  - p5 (EIQTAVR) is shared by most proteins → will be excluded
  - p2 (QVHPDTGISSK) is shared by many → likely excluded
  - Groups 1 and 2 will compete in greedy set cover
  - Proteins within each group with identical peptide sets will be merged

Expected Result:
  - Proteins with identical remaining peptide sets should be merged
  - Only unique peptides assigned to protein groups
  - Highly shared peptides (p5, p2) excluded from results
""")

# Data from failing_component_2025-10-13T12-00-28.902.jl
# Simplified protein/peptide names for readability:
# Original → Simple:
#   O60814→A, P33778→B, P02293→O, Q6DRA6→Q, Q93079→C, Q99880→D, P58876→E
#   P23527→L, Q99879→F, Q5QNW6→G, P62807→H, P02294→P, P06899→M, Q99877→I
#   Q96A08→S, Q16778→N, Q6DN03→R, Q8N257→J, P57053→K
#   ESYSVYVYK→p1, QVHPDTGISSK→p2, KESYSVYVYK→p3, AMGIMNSFVNDIFER→p4, EIQTAVR→p5
#   ESYSIYVYK→p6, KESYSIYVYK→p7, ETYSSYIYK→p8, SMSILNSFVNDIFER→p9
#   ESYSIYIYK→p10, AMSIMNSFVTDIFER→p11

#=
  Proteins: A-S (19 proteins organized into 6 groups)
  - Group 1: A,C,D,E,F,G,H,I,K → {p1,p2,p3,p4,p5}
  - Group 2: B,L,M,N → {p6,p2,p7,p4,p5}
  - Group 3: O,P (histones) → {p8,p9,p5}
  - Group 4: Q,R → {p6,p7}
  - Group 5: S → {p2,p10,p11,p5}
  - Group 6: J → {p6,p2,p7,p4} (no p5)
=#

proteins_realworld = [
    ProteinKey("A", true, UInt8(0)),  # A -> {p1,p2,p3,p4,p5} (Group 1)
    ProteinKey("A", true, UInt8(0)),
    ProteinKey("A", true, UInt8(0)),
    ProteinKey("A", true, UInt8(0)),
    ProteinKey("A", true, UInt8(0)),
    ProteinKey("B", true, UInt8(0)),  # B -> {p6,p2,p7,p4,p5} (Group 2)
    ProteinKey("B", true, UInt8(0)),
    ProteinKey("B", true, UInt8(0)),
    ProteinKey("B", true, UInt8(0)),
    ProteinKey("B", true, UInt8(0)),
    ProteinKey("O", true, UInt8(0)),  # O -> {p8,p9,p5} (Histone, Group 3)
    ProteinKey("O", true, UInt8(0)),
    ProteinKey("O", true, UInt8(0)),
    ProteinKey("Q", true, UInt8(0)),  # Q -> {p6,p7} (Group 4)
    ProteinKey("Q", true, UInt8(0)),
    ProteinKey("C", true, UInt8(0)),  # C -> {p1,p2,p3,p4,p5} (Group 1)
    ProteinKey("C", true, UInt8(0)),
    ProteinKey("C", true, UInt8(0)),
    ProteinKey("C", true, UInt8(0)),
    ProteinKey("C", true, UInt8(0)),
    ProteinKey("D", true, UInt8(0)),  # D -> {p1,p2,p3,p4,p5} (Group 1)
    ProteinKey("D", true, UInt8(0)),
    ProteinKey("D", true, UInt8(0)),
    ProteinKey("D", true, UInt8(0)),
    ProteinKey("D", true, UInt8(0)),
    ProteinKey("E", true, UInt8(0)),  # E -> {p1,p2,p3,p4,p5} (Group 1)
    ProteinKey("E", true, UInt8(0)),
    ProteinKey("E", true, UInt8(0)),
    ProteinKey("E", true, UInt8(0)),
    ProteinKey("E", true, UInt8(0)),
    ProteinKey("L", true, UInt8(0)),  # L -> {p6,p2,p7,p4,p5} (Group 2)
    ProteinKey("L", true, UInt8(0)),
    ProteinKey("L", true, UInt8(0)),
    ProteinKey("L", true, UInt8(0)),
    ProteinKey("L", true, UInt8(0)),
    ProteinKey("F", true, UInt8(0)),  # F -> {p1,p2,p3,p4,p5} (Group 1)
    ProteinKey("F", true, UInt8(0)),
    ProteinKey("F", true, UInt8(0)),
    ProteinKey("F", true, UInt8(0)),
    ProteinKey("F", true, UInt8(0)),
    ProteinKey("G", true, UInt8(0)),  # G -> {p1,p2,p3,p4,p5} (Group 1)
    ProteinKey("G", true, UInt8(0)),
    ProteinKey("G", true, UInt8(0)),
    ProteinKey("G", true, UInt8(0)),
    ProteinKey("G", true, UInt8(0)),
    ProteinKey("H", true, UInt8(0)),  # H -> {p1,p2,p3,p4,p5} (Group 1)
    ProteinKey("H", true, UInt8(0)),
    ProteinKey("H", true, UInt8(0)),
    ProteinKey("H", true, UInt8(0)),
    ProteinKey("H", true, UInt8(0)),
    ProteinKey("P", true, UInt8(0)),  # P -> {p8,p9,p5} (Histone, Group 3)
    ProteinKey("P", true, UInt8(0)),
    ProteinKey("P", true, UInt8(0)),
    ProteinKey("M", true, UInt8(0)),  # M -> {p6,p2,p7,p4,p5} (Group 2)
    ProteinKey("M", true, UInt8(0)),
    ProteinKey("M", true, UInt8(0)),
    ProteinKey("M", true, UInt8(0)),
    ProteinKey("M", true, UInt8(0)),
    ProteinKey("I", true, UInt8(0)),  # I -> {p1,p2,p3,p4,p5} (Group 1)
    ProteinKey("I", true, UInt8(0)),
    ProteinKey("I", true, UInt8(0)),
    ProteinKey("I", true, UInt8(0)),
    ProteinKey("I", true, UInt8(0)),
    ProteinKey("S", true, UInt8(0)),  # S -> {p2,p10,p11,p5} (Group 5)
    ProteinKey("S", true, UInt8(0)),
    ProteinKey("S", true, UInt8(0)),
    ProteinKey("S", true, UInt8(0)),
    ProteinKey("N", true, UInt8(0)),  # N -> {p6,p2,p7,p4,p5} (Group 2)
    ProteinKey("N", true, UInt8(0)),
    ProteinKey("N", true, UInt8(0)),
    ProteinKey("N", true, UInt8(0)),
    ProteinKey("N", true, UInt8(0)),
    ProteinKey("R", true, UInt8(0)),  # R -> {p6,p7} (Group 4)
    ProteinKey("R", true, UInt8(0)),
    ProteinKey("J", true, UInt8(0)),  # J -> {p6,p2,p7,p4} (Group 6, no p5!)
    ProteinKey("J", true, UInt8(0)),
    ProteinKey("J", true, UInt8(0)),
    ProteinKey("J", true, UInt8(0)),
    ProteinKey("K", true, UInt8(0)),  # K -> {p1,p2,p3,p4,p5} (Group 1)
    ProteinKey("K", true, UInt8(0)),
    ProteinKey("K", true, UInt8(0)),
    ProteinKey("K", true, UInt8(0)),
    ProteinKey("K", true, UInt8(0)),
]

peptides_realworld = [
    PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
    PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR (highly shared!)
    PeptideKey("p6", true, UInt8(0)),  # ESYSIYVYK
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p7", true, UInt8(0)),  # KESYSIYVYK
    PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
    PeptideKey("p8", true, UInt8(0)),  # ETYSSYIYK
    PeptideKey("p9", true, UInt8(0)),  # SMSILNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
    PeptideKey("p6", true, UInt8(0)),  # ESYSIYVYK
    PeptideKey("p7", true, UInt8(0)),  # KESYSIYVYK
    PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
    PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
    PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
    PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
    PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
    PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
    PeptideKey("p6", true, UInt8(0)),  # ESYSIYVYK
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p7", true, UInt8(0)),  # KESYSIYVYK
    PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
    PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
    PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
    PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
    PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
    PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
    PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
    PeptideKey("p8", true, UInt8(0)),  # ETYSSYIYK
    PeptideKey("p9", true, UInt8(0)),  # SMSILNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
    PeptideKey("p6", true, UInt8(0)),  # ESYSIYVYK
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p7", true, UInt8(0)),  # KESYSIYVYK
    PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
    PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
    PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p10", true, UInt8(0)),  # ESYSIYIYK
    PeptideKey("p11", true, UInt8(0)),  # AMSIMNSFVTDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
    PeptideKey("p6", true, UInt8(0)),  # ESYSIYVYK
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p7", true, UInt8(0)),  # KESYSIYVYK
    PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
    PeptideKey("p6", true, UInt8(0)),  # ESYSIYVYK
    PeptideKey("p7", true, UInt8(0)),  # KESYSIYVYK
    PeptideKey("p6", true, UInt8(0)),  # ESYSIYVYK
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p7", true, UInt8(0)),  # KESYSIYVYK
    PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
    PeptideKey("p1", true, UInt8(0)),  # ESYSVYVYK
    PeptideKey("p2", true, UInt8(0)),  # QVHPDTGISSK
    PeptideKey("p3", true, UInt8(0)),  # KESYSVYVYK
    PeptideKey("p4", true, UInt8(0)),  # AMGIMNSFVNDIFER
    PeptideKey("p5", true, UInt8(0)),  # EIQTAVR
]

println("\nRunning inference on real-world failing component...")
result_realworld = infer_proteins(proteins_realworld, peptides_realworld)

println("SUCCESS! The sort() operation now works with Base.isless defined.")
println("\nResults summary:")
println("  Total unique peptides in result: $(length(result_realworld.peptide_to_protein))")
println("  (Many peptides excluded as shared - this is correct behavior)")

# Show a few example assignments
println("\nSample peptide assignments:")
count = 0
for (pep, prot) in pairs(result_realworld.peptide_to_protein)
    #if count < 5
        println("  $(pep.sequence) → $(prot.name)")
        count += 1
    #end
end

println("\nAnalysis:")
println("  - Proteins in Group 1 (A,C,D,E,F,G,H,I,K) share identical peptide sets")
println("  - They should be merged into a single protein group during inference")
println("  - Peptide p5 (EIQTAVR) is shared by most proteins → excluded from result")
println("  - Peptide p2 (QVHPDTGISSK) is also highly shared → likely excluded")
println("  - The merge-first approach ensures proteins with identical peptides are combined")

println("\n" * "="^80)
println("All examples completed successfully!")
println("="^80)
