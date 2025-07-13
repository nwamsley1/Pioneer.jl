# Simple test to load and examine protein data
push!(LOAD_PATH, dirname(dirname(@__DIR__)))

using EntrapmentAnalysis
using Arrow
using DataFrames

# Based on the schema shown by the user, create a mock protein dataset
# Schema from user:
# :file_name        Union{Missing, String}
# :target           Union{Missing, Bool}
# :entrap_id        Union{Missing, UInt8}
# :species          Union{Missing, String}
# :protein          Union{Missing, String}
# :peptides         Union{Missing, SubArray{Union{Missing, UInt32}, 1, Arrow.Primitive{Union{Missing, UInt32}, Vector{UInt32}}, Tuple{UnitRange{Int64}}, true}}
# :n_peptides       Union{Missing, UInt32}
# :global_qval      Union{Missing, Float32}
# :qval             Union{Missing, Float32}
# :pg_pep           Union{Missing, Float32}
# :pg_score         Union{Missing, Float32}
# :global_pg_score  Union{Missing, Float32}
# :abundance        Float32

# Create mock protein data that mimics the structure
mock_protein_data = DataFrame(
    file_name = ["file1.raw", "file1.raw", "file1.raw", "file1.raw", 
                 "file2.raw", "file2.raw", "file2.raw", "file2.raw"],
    target = [true, true, true, true, true, true, true, true],
    entrap_id = UInt8[0, 1, 0, 2, 0, 1, 0, 2],  # 0 = original, 1,2 = entrapments
    species = ["YEAST", "YEAST", "YEAST", "YEAST", "YEAST", "YEAST", "YEAST", "YEAST"],
    protein = ["sp|P00330|ADH1_YEAST", "sp|P00330|ADH1_YEAST", 
               "sp|P00331|ADH2_YEAST", "sp|P00331|ADH2_YEAST",
               "sp|P00330|ADH1_YEAST", "sp|P00330|ADH1_YEAST", 
               "sp|P00331|ADH2_YEAST", "sp|P00331|ADH2_YEAST"],
    peptides = [missing, missing, missing, missing, missing, missing, missing, missing],  # Simplified
    n_peptides = UInt32[10, 10, 8, 8, 10, 10, 8, 8],
    global_qval = Float32[0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.02, 0.02],
    qval = Float32[0.01, 0.015, 0.02, 0.025, 0.01, 0.012, 0.02, 0.022],
    pg_pep = Float32[0.001, 0.001, 0.002, 0.002, 0.001, 0.001, 0.002, 0.002],
    pg_score = Float32[100.5, 95.3, 80.2, 75.1, 98.7, 96.1, 82.3, 77.4],
    global_pg_score = Float32[100.5, 95.3, 80.2, 75.1, 100.5, 95.3, 82.3, 77.4],
    abundance = Float32[1e6, 0.9e6, 0.8e6, 0.75e6, 1.1e6, 0.95e6, 0.85e6, 0.78e6]
)

println("=== Mock Protein Data Structure ===")
println("Columns: ", names(mock_protein_data))
println("Number of rows: ", nrow(mock_protein_data))
println("\nFirst few rows:")
show(first(mock_protein_data, 5), allcols=true)

println("\n\n=== Key Observations for Implementation ===")
println("1. File identification: Uses 'file_name' column (not ms_file_idx)")
println("2. Entrapment identification: 'entrap_id' column (0=original, >0=entrapment)")
println("3. Protein identification: 'protein' column contains full protein names")
println("4. Score columns: pg_score (per-file), global_pg_score (global)")
println("5. Q-value columns: qval (per-file), global_qval (global)")

println("\n=== Entrapment Pairing Logic ===")
# Group by protein name to see pairing needs
grouped = groupby(mock_protein_data, :protein)
for (key, group) in pairs(grouped)
    println("\nProtein: $(key.protein)")
    unique_entrap = sort(unique(group.entrap_id))
    println("  Entrap IDs present: $unique_entrap")
    println("  Files present: $(unique(group.file_name))")
end