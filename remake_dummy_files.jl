#!/usr/bin/env julia

using Arrow
using DataFrames
using Random

# First, let's examine the existing dummy files to understand their structure
dummy_dir = "/Users/nathanwamsley/Data/MS_DATA/ARROW/SingleCell/250pg/24ms/arrow_out"
existing_files = [
    "dummy_empty.arrow",
    "dummy_empty_scans.arrow", 
    "dummy_sparse.arrow",
    "dummy_minimal.arrow"
]

println("=== Examining existing dummy files ===")
for file in existing_files
    filepath = joinpath(dummy_dir, file)
    if isfile(filepath)
        println("\n--- $file ---")
        try
            table = Arrow.Table(filepath)
            println("Columns: ", propertynames(table))
            println("Number of rows: ", length(table))
            
            # Check the types of m/z and intensity arrays
            if :mz_array in propertynames(table)
                mz_type = eltype(table.mz_array)
                println("mz_array type: $mz_type")
                if length(table) > 0
                    first_mz = table.mz_array[1]
                    println("First mz_array element type: $(typeof(first_mz))")
                    if !isempty(first_mz)
                        println("First mz value type: $(typeof(first_mz[1]))")
                    end
                end
            end
            
            if :intensity_array in propertynames(table)
                intensity_type = eltype(table.intensity_array)
                println("intensity_array type: $intensity_type")
                if length(table) > 0
                    first_intensity = table.intensity_array[1]
                    println("First intensity_array element type: $(typeof(first_intensity))")
                    if !isempty(first_intensity)
                        println("First intensity value type: $(typeof(first_intensity[1]))")
                    end
                end
            end
        catch e
            println("Error reading $file: $e")
        end
    else
        println("File not found: $filepath")
    end
end

println("\n=== Creating new dummy files with correct types ===")

# Function to create a dummy MS data file with proper Union{Missing, Float32} types
function create_dummy_ms_file(output_path::String, n_scans::Int, scan_type::Symbol)
    Random.seed!(42)  # For reproducible results
    
    # Create basic scan metadata
    scan_numbers = collect(1:n_scans)
    scan_headers = ["scan_header_$i" for i in 1:n_scans]
    packet_types = fill(Int32(1), n_scans)  # MS data packet type
    
    # Create retention times (Float32)
    retention_times = Float32.(sort(rand(n_scans) * 60.0))  # 0-60 minutes
    
    # Create m/z ranges (Float32)
    low_mzs = fill(Float32(100.0), n_scans)
    high_mzs = fill(Float32(2000.0), n_scans)
    
    # Create TIC values (Float32)
    tics = Float32.(rand(n_scans) * 1e6)
    
    # Create center m/z and isolation width (Union{Missing, Float32})
    center_mzs = Union{Missing, Float32}[]
    isolation_widths = Union{Missing, Float32}[]
    ms_orders = UInt8[]
    
    # Create m/z and intensity arrays based on scan type
    mz_arrays = Vector{Union{Missing, Float32}}[]
    intensity_arrays = Vector{Union{Missing, Float32}}[]
    
    for i in 1:n_scans
        # Alternate between MS1 and MS2 scans
        if i % 4 == 1
            # MS1 scan
            push!(ms_orders, UInt8(1))
            push!(center_mzs, missing)
            push!(isolation_widths, missing)
        else
            # MS2 scan
            push!(ms_orders, UInt8(2))
            push!(center_mzs, Float32(500.0 + rand() * 1000.0))  # Random precursor m/z
            push!(isolation_widths, Float32(1.0))  # 1 Da isolation window
        end
        
        if scan_type == :empty
            # Empty arrays
            push!(mz_arrays, Union{Missing, Float32}[])
            push!(intensity_arrays, Union{Missing, Float32}[])
        elseif scan_type == :minimal
            # Minimal data - just a few peaks
            n_peaks = rand(1:5)
            mzs = Union{Missing, Float32}[Float32(200.0 + rand() * 1600.0) for _ in 1:n_peaks]
            intensities = Union{Missing, Float32}[Float32(rand() * 1000.0) for _ in 1:n_peaks]
            push!(mz_arrays, sort(mzs))
            push!(intensity_arrays, intensities)
        elseif scan_type == :sparse
            # Sparse data - reasonable number of peaks
            n_peaks = rand(10:50)
            mzs = Union{Missing, Float32}[Float32(200.0 + rand() * 1600.0) for _ in 1:n_peaks]
            intensities = Union{Missing, Float32}[Float32(rand() * 10000.0) for _ in 1:n_peaks]
            push!(mz_arrays, sort(mzs))
            push!(intensity_arrays, intensities)
        end
    end
    
    # Create the DataFrame with proper types
    df = DataFrame(
        scanNumber = scan_numbers,
        scanHeader = scan_headers,
        packetType = packet_types,
        retentionTime = retention_times,
        lowMz = low_mzs,
        highMz = high_mzs,
        TIC = tics,
        centerMz = center_mzs,
        isolationWidthMz = isolation_widths,
        msOrder = ms_orders,
        mz_array = mz_arrays,
        intensity_array = intensity_arrays
    )
    
    # Write to Arrow file
    Arrow.write(output_path, df)
    println("Created: $output_path")
    
    # Verify the types
    table = Arrow.Table(output_path)
    println("  mz_array type: $(eltype(table.mz_array))")
    println("  intensity_array type: $(eltype(table.intensity_array))")
    if length(table) > 0 && !isempty(table.mz_array[1])
        println("  First mz value type: $(typeof(table.mz_array[1][1]))")
        println("  First intensity value type: $(typeof(table.intensity_array[1][1]))")
    end
    
    return output_path
end

# Create new dummy files in the same directory
new_dummy_dir = "/Users/nathanwamsley/Data/MS_DATA/ARROW/SingleCell/250pg/24ms/arrow_out"

# Create the corrected dummy files
println("\nCreating corrected dummy files...")
create_dummy_ms_file(joinpath(new_dummy_dir, "dummy_empty_corrected.arrow"), 100, :empty)
create_dummy_ms_file(joinpath(new_dummy_dir, "dummy_minimal_corrected.arrow"), 10, :minimal)
create_dummy_ms_file(joinpath(new_dummy_dir, "dummy_sparse_corrected.arrow"), 100, :sparse)

# Create one file with only empty scans (different from completely empty)
println("\nCreating empty scans file...")
create_dummy_ms_file(joinpath(new_dummy_dir, "dummy_empty_scans_corrected.arrow"), 100, :empty)

# Keep one file with the old pure Float32 type for testing
println("\nCopying one original file for testing type mismatch...")
old_file = joinpath(new_dummy_dir, "dummy_minimal.arrow")
test_file = joinpath(new_dummy_dir, "dummy_minimal_pure_float32.arrow")
if isfile(old_file)
    cp(old_file, test_file, force=true)
    println("Created test file with pure Float32 types: $test_file")
else
    println("Original file not found: $old_file")
end

println("\n=== Summary ===")
println("Created corrected dummy files with Union{Missing, Float32} types:")
println("- dummy_empty_corrected.arrow")
println("- dummy_minimal_corrected.arrow") 
println("- dummy_sparse_corrected.arrow")
println("- dummy_empty_scans_corrected.arrow")
println("- dummy_minimal_pure_float32.arrow (for testing type mismatch)")