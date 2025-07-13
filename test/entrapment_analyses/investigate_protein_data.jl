# Script to investigate protein data structure
using Arrow
using DataFrames

println("Loading protein data...")
protein_data_path = "/Users/nathanwamsley/Documents/PIONEER/RAW/YEAST_TEST/MtacYeastAlternating3M/yeast_test_1p/protein_groups_long.arrow"

# Load the data
df = DataFrame(Arrow.Table(protein_data_path))

println("\n=== Data Overview ===")
println("Number of rows: ", nrow(df))
println("Number of columns: ", ncol(df))

println("\n=== Column Information ===")
for col in names(df)
    col_type = eltype(df[!, col])
    println("$col: $col_type")
end

println("\n=== First 5 rows ===")
show(first(df, 5), allcols=true)

println("\n\n=== Unique entrap_id values ===")
println(sort(unique(skipmissing(df.entrap_id))))

println("\n=== Sample protein names ===")
println("First 10 unique proteins:")
unique_proteins = unique(skipmissing(df.protein))
for (i, prot) in enumerate(first(unique_proteins, 10))
    println("  $i. $prot")
end

println("\n=== File information ===")
if hasproperty(df, :ms_file_idx)
    println("Has ms_file_idx column")
    println("Unique ms_file_idx values: ", sort(unique(skipmissing(df.ms_file_idx))))
else
    println("No ms_file_idx column")
end

if hasproperty(df, :file_name)
    println("Has file_name column")
    println("Unique file names: ")
    for fname in unique(skipmissing(df.file_name))
        println("  - $fname")
    end
end

println("\n=== Score columns ===")
score_cols = filter(x -> occursin("score", String(x)) || occursin("qval", String(x)) || x == :pg_pep, names(df))
println("Score-related columns: ", score_cols)

println("\n=== Entrapment Analysis ===")
# Group by protein to see entrapment patterns
protein_groups = groupby(df, :protein)
println("Number of unique proteins: ", length(protein_groups))

# Find proteins with multiple entrap_id values
proteins_with_entrapments = []
for (key, group) in pairs(protein_groups)
    unique_entrap_ids = unique(skipmissing(group.entrap_id))
    if length(unique_entrap_ids) > 1
        push!(proteins_with_entrapments, (
            protein = key.protein,
            entrap_ids = sort(unique_entrap_ids),
            count = nrow(group)
        ))
    end
end

println("\nProteins with multiple entrap_id values: ", length(proteins_with_entrapments))
if length(proteins_with_entrapments) > 0
    println("First 5 examples:")
    for (i, example) in enumerate(first(proteins_with_entrapments, 5))
        println("  $i. $(example.protein)")
        println("     Entrap IDs: $(example.entrap_ids)")
        println("     Total rows: $(example.count)")
    end
end

# Check target column
println("\n=== Target column ===")
if hasproperty(df, :target)
    target_counts = combine(groupby(df, :target), nrow => :count)
    println("Target value counts:")
    for row in eachrow(target_counts)
        println("  target=$(row.target): $(row.count) rows")
    end
else
    println("No target column found")
end

# Check for missing values in key columns
println("\n=== Missing values in key columns ===")
key_cols = [:protein, :entrap_id, :pg_score, :global_pg_score]
for col in key_cols
    if hasproperty(df, col)
        missing_count = sum(ismissing.(df[!, col]))
        println("$col: $missing_count missing values ($(round(missing_count/nrow(df)*100, digits=2))%)")
    end
end

# Sample a protein with entrapments to understand structure
println("\n=== Example: Detailed view of one protein with entrapments ===")
if length(proteins_with_entrapments) > 0
    example_protein = proteins_with_entrapments[1].protein
    example_data = df[df.protein .== example_protein, :]
    println("Protein: $example_protein")
    show(example_data, allcols=true)
end