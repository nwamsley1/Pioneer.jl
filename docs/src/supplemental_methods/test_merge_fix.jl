# Test script for the merge-first bug fix

using Pkg
Pkg.activate(".")

using Pioneer
using Pioneer: infer_proteins, ProteinKey, PeptideKey

println("="^80)
println("Testing Merge-First Bug Fix")
println("="^80)
println()

# The failing example from the bug report
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

println("Test Case: Indistinguishable Proteins B and C")
println("-"^80)
println("Protein-Peptide Mapping:")
println("  Peptide 1 → {A}")
println("  Peptide 2 → {A, B, C}")
println("  Peptide 3 → {B, C}")
println("  Peptide 4 → {B, C, D}")
println("  Peptide 5 → {D}")
println()
println("Expected Behavior:")
println("  Phase 1: Select A (unique pep1), D (unique pep5)")
println("  Phase 2: B and C have identical remaining peptides {pep3}")
println("           → Merge into 'B;C' → Select 'B;C'")
println()
println("Expected Results:")
println("  pep1 → A (use_for_quant = true)")
println("  pep3 → B;C (use_for_quant = true)  ← KEY FIX")
println("  pep5 → D (use_for_quant = true)")
println()

result = infer_proteins(proteins, peptides)

println("Actual Results:")
println("-"^80)
for pep in peptides
    if haskey(result.peptide_to_protein, pep)
        prot = result.peptide_to_protein[pep]
        quant = result.use_for_quant[pep]
        println("  Peptide $(pep.sequence) → Protein $(prot.name)")
        println("    use_for_quant: $(quant)")
    end
end
println()

# Check if the bug is fixed
pep3 = peptides[3]
if haskey(result.peptide_to_protein, pep3)
    prot3 = result.peptide_to_protein[pep3]
    if prot3.name == "B;C"
        println("✓ BUG FIXED: Peptide 3 correctly assigned to 'B;C'")
        if result.use_for_quant[pep3]
            println("✓ use_for_quant correctly set to true")
        else
            println("✗ ERROR: use_for_quant should be true")
        end
    else
        println("✗ BUG NOT FIXED: Peptide 3 assigned to '$(prot3.name)' instead of 'B;C'")
    end
else
    println("✗ ERROR: Peptide 3 not found in results")
end
