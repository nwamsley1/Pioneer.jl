"""
Example: Protein-level EFDR analysis with entrapment sequences

This example shows how to run empirical FDR analysis on protein-level data
from proteomics experiments with entrapment sequences.
"""

push!(LOAD_PATH, dirname(dirname(@__DIR__)))
using EntrapmentAnalysis
using DataFrames
using Arrow

# Example 1: Basic protein EFDR analysis
println("=== Example 1: Basic Protein EFDR Analysis ===\n")

# Path to protein results file (update with your actual file path)
protein_results_path = "/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_test_1p/protein_groups_long.arrow"

# Check if file exists
if isfile(protein_results_path)
    # Run the analysis
    results = run_protein_efdr_analysis(
        protein_results_path;
        output_dir="protein_efdr_output",
        score_qval_pairs=[(:global_pg_score, :global_qval), (:pg_score, :qval)],
        verbose=true
    )
    
    println("\nAnalysis complete! Check the output directory for results.")
    println("Generated files:")
    for file in results.output_files
        println("  - $file")
    end
else
    println("File not found: $protein_results_path")
    println("Using mock data for demonstration...")
    
    # Create mock protein data
    mock_protein_data = DataFrame(
        file_name = repeat(["file1.raw", "file2.raw"], inner=6),
        target = fill(true, 12),
        entrap_id = repeat(UInt8[0, 1, 2, 0, 1, 2], 2),
        species = fill("YEAST", 12),
        protein = repeat(["sp|P00330|ADH1_YEAST", "sp|P00330|ADH1_YEAST", "sp|P00330|ADH1_YEAST",
                         "sp|P00331|ADH2_YEAST", "sp|P00331|ADH2_YEAST", "sp|P00331|ADH2_YEAST"], 2),
        n_peptides = UInt32[10, 10, 10, 8, 8, 8, 10, 10, 10, 8, 8, 8],
        global_qval = Float32[0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02],
        qval = Float32[0.01, 0.015, 0.018, 0.02, 0.025, 0.028, 0.01, 0.012, 0.016, 0.02, 0.022, 0.026],
        pg_pep = Float32[0.001, 0.001, 0.001, 0.002, 0.002, 0.002, 0.001, 0.001, 0.001, 0.002, 0.002, 0.002],
        pg_score = Float32[100.5, 95.3, 92.1, 80.2, 75.1, 72.5, 98.7, 96.1, 93.2, 82.3, 77.4, 74.8],
        global_pg_score = Float32[100.5, 95.3, 92.1, 80.2, 75.1, 72.5, 100.5, 95.3, 92.1, 82.3, 77.4, 74.8],
        abundance = Float32[1e6, 0.9e6, 0.85e6, 0.8e6, 0.75e6, 0.7e6, 1.1e6, 0.95e6, 0.88e6, 0.85e6, 0.78e6, 0.73e6]
    )
    
    # Save mock data
    mock_path = "mock_protein_data.arrow"
    Arrow.write(mock_path, mock_protein_data)
    
    # Run analysis on mock data
    results = run_protein_efdr_analysis(
        mock_path;
        output_dir="protein_efdr_demo",
        verbose=true
    )
    
    println("\nDemo analysis complete!")
end

# Example 2: Manual protein entrapment analysis
println("\n\n=== Example 2: Manual Protein Entrapment Analysis ===\n")

# Create example protein data
protein_df = DataFrame(
    protein = ["PROT_A", "PROT_A", "PROT_A", "PROT_B", "PROT_B", "PROT_B"],
    entrap_id = UInt8[0, 1, 2, 0, 1, 2],
    file_name = ["exp1", "exp1", "exp1", "exp1", "exp1", "exp1"],
    pg_score = Float32[150.0, 140.0, 135.0, 120.0, 110.0, 105.0],
    global_pg_score = Float32[150.0, 140.0, 135.0, 120.0, 110.0, 105.0]
)

println("Original data:")
println(protein_df)

# Step 1: Assign entrapment pairs
println("\n1. Assigning entrapment pairs...")
assign_protein_entrapment_pairs!(protein_df)
println("Pairs assigned:")
println(select(protein_df, :protein, :entrap_id, :entrap_pair_id))

# Step 2: Add original target scores
println("\n2. Adding original target scores...")
add_original_target_protein_scores!(protein_df; score_col=:pg_score)
println("With original target scores:")
println(select(protein_df, :protein, :entrap_id, :pg_score, :pg_score_original_target))

# Step 3: Create global results
println("\n3. Creating global protein results...")
global_df = create_global_protein_results_df(protein_df; score_col=:global_pg_score)
println("Global results (best score per protein):")
println(global_df)

# Example 3: Analyzing proteins with different entrapment patterns
println("\n\n=== Example 3: Different Entrapment Patterns ===\n")

# Scenario with uneven entrapment groups
complex_df = DataFrame(
    protein = ["PROT_X", "PROT_X", "PROT_X", "PROT_Y", "PROT_Y", "PROT_Z"],
    entrap_id = UInt8[0, 1, 0, 0, 1, 0],  # PROT_X has 2 originals, PROT_Z has no entrapment
    file_name = ["exp1", "exp1", "exp2", "exp1", "exp1", "exp1"],
    pg_score = Float32[200.0, 190.0, 195.0, 180.0, 170.0, 160.0]
)

println("Complex scenario:")
println(complex_df)

assign_protein_entrapment_pairs!(complex_df)
println("\nAfter pairing:")
println(select(complex_df, :protein, :entrap_id, :file_name, :entrap_pair_id))

# Note how:
# - PROT_X: Two originals are paired with the single entrapment (round-robin)
# - PROT_Y: Simple 1:1 pairing
# - PROT_Z: No pairing (no entrapment version)

println("\n=== Summary ===")
println("Protein-level EFDR analysis differs from precursor-level:")
println("1. Proteins pair by name only (no charge/modification considerations)")
println("2. Each protein can appear in multiple files")
println("3. Global analysis selects best-scoring instance across files")
println("4. Simpler than precursor analysis but follows same EFDR principles")