#
#=
Have N of these tables. Need to combine into one sorted Arrow table without loading all tables
into memory at once. 
julia> DataFrame(Arrow.Table(readdir(second_quant_folder, join = true)[1]))
280488×11 DataFrame
    Row │ precursor_idx  prob      weight         target  irt_obs    missed_cleavage  isotopes_captured  scan_idx  ms_file_idx  peak_area   new_best_scan 
        │ UInt32         Float32   Float32        Bool    Float32    UInt8            Tuple{Int8, Int8}  UInt32    Int64        Float32     UInt32        
────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
      1 │       1468734  0.808784    201.467        true   0.659221                1  (0, 3)                49270            1     6.03541          50180
      2 │        262434  0.989585   2696.17         true   0.659221                0  (0, 3)                76753            1   121.201            76753
=#


using Arrow
using DataFrames
using Heaps

function merge_sorted_arrow_tables(input_paths, output_path)
    tables = [Arrow.Table(path) for path in input_paths]
    iters = [Iterators.Stateful(table) for table in tables]
    heap = BinaryMinHeap{Tuple{Int,Int,Int}}()
    
    for (i, iter) in enumerate(iters)
        if !isempty(iter)
            row = first(iter)
            push!(heap, (row.precursor_idx, row.precursor_idx, i))
        end
    end
    
    Arrow.write(output_path, Tables.Schema(first(tables))) do writer
        while !isempty(heap)
            _, _, table_idx = pop!(heap)
            row = popfirst!(iters[table_idx])
            Arrow.write(writer, Tables.Row(row))
            
            if !isempty(iters[table_idx])
                next_row = first(iters[table_idx])
                push!(heap, (next_row.precursor_idx, next_row.precursor_idx, table_idx))
            end
        end
    end
end

function cascade_merge(input_paths, final_output_path)
    temp_dir = mktempdir()
    current_paths = copy(input_paths)
    
    while length(current_paths) > 1
        new_paths = String[]
        
        for chunk in Iterators.partition(current_paths, MAX_OPEN_FILES - 1)
            output_path = joinpath(temp_dir, "merged_$(hash(chunk)).arrow")
            merge_sorted_arrow_tables(chunk, output_path)
            push!(new_paths, output_path)
        end
        
        # Clean up intermediate files
        for path in current_paths
            if startswith(path, temp_dir)
                rm(path)
            end
        end
        
        current_paths = new_paths
    end
    
    # Rename the final merged file
    mv(first(current_paths), final_output_path, force=true)
    
    # Clean up temporary directory
    rm(temp_dir, recursive=true)
end

# Usage
input_paths = ["table1.arrow", "table2.arrow", "table3.arrow", "table4.arrow", "table5.arrow"]
sorted_paths = [path * ".sorted" for path in input_paths]

# Step 1: Sort each table individually
for (input, output) in zip(input_paths, sorted_paths)
    sort_arrow_table(input, output)
end

# Step 2: Cascade merge sorted tables
cascade_merge(sorted_paths, "final_combined_sorted.arrow")

# Clean up temporary sorted files
for path in sorted_paths
    rm(path)
end