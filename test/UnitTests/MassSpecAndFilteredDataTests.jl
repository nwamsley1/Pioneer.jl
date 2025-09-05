using Test
using Arrow
using DataFrames
using Random

# Import only what we need from Pioneer (avoid Pioneer.func calls)
using Pioneer: BasicMassSpecData,
               getScanHeader, getScanNumber, getBasePeakMz, getBasePeakIntensity,
               getRetentionTime, getLowMz, getHighMz, getTIC, getMsOrder,
               getMzArray, getIntensityArray,
               FilteredMassSpecData,
               getOriginalScanIndex, getOriginalScanIndices,
               compute_rt_bins, sort_scans_by_peak_density, create_priority_order

function make_union_vec(xs::Vector{Float32})
    return Vector{Union{Missing, Float32}}(xs)
end

function write_basic_ms_arrow(path::String; n::Int=5)
    mz_arrays = Vector{Vector{Union{Missing,Float32}}}(undef, n)
    int_arrays = Vector{Vector{Union{Missing,Float32}}}(undef, n)
    headers = ["scan_$(i)" for i in 1:n]
    scan_nums = Int32.(1:n)
    base_mz = Float32.(100 .+ collect(1:n))
    base_int = Float32.(1000 .+ 10 .* collect(1:n))
    rt = Float32.(range(10, length=n, step=0.5))
    low_mz = Float32.(fill(100.0, n))
    high_mz = Float32.(fill(1000.0, n))
    tic = Float32.(10000 .+ 100 .* collect(1:n))
    center_mz = Vector{Union{Missing,Float32}}(Float32.(200 .+ collect(1:n)))
    isol_width = Vector{Union{Missing,Float32}}(Float32.(fill(1.0, n)))
    ms_order = UInt8.(fill(2, n))
    packet_type = Int32.(fill(0, n))

    for i in 1:n
        mz_arrays[i] = make_union_vec(Float32[100.0 + i, 150.0 + i, 200.0 + i])
        int_arrays[i] = make_union_vec(Float32[1000.0 + i, 900.0 + i, 800.0 + i])
    end

    df = DataFrame(
        mz_array = mz_arrays,
        intensity_array = int_arrays,
        scanHeader = headers,
        scanNumber = scan_nums,
        basePeakMz = base_mz,
        basePeakIntensity = base_int,
        retentionTime = rt,
        lowMz = low_mz,
        highMz = high_mz,
        TIC = tic,
        centerMz = center_mz,
        isolationWidthMz = isol_width,
        msOrder = ms_order,
        packetType = packet_type,
    )
    Arrow.write(path, df)
    return path
end

@testset "MassSpec + Filtered Data (basic + filtered)" begin
    # Create temp Arrow file
    temp_dir = mktempdir()
    file_path = joinpath(temp_dir, "ms.arrow")
    write_basic_ms_arrow(file_path, n=6)

    # Load as Basic
    original = BasicMassSpecData(file_path)
    @test length(original) == 6

    # Spot-check getters (Basic)
    @test getScanHeader(original, 1) == "scan_1"
    @test getScanNumber(original, 2) == Int32(2)
    @test getBasePeakMz(original, 3) ≈ Float32(103)
    @test getBasePeakIntensity(original, 4) ≈ Float32(1040)
    @test getRetentionTime(original, 5) ≈ Float32(12.0)
    @test getLowMz(original, 1) ≈ Float32(100.0)
    @test getHighMz(original, 1) ≈ Float32(1000.0)
    @test getTIC(original, 6) ≈ Float32(10600)
    @test getMsOrder(original, 1) == UInt8(2)

    # Array getters
    mz1 = getMzArray(original, 1)
    it1 = getIntensityArray(original, 1)
    @test length(mz1) == 3
    @test length(it1) == 3
    @test all(!ismissing, mz1)
    @test all(!ismissing, it1)

    # Build filtered data (MS2 only)
    filtered = FilteredMassSpecData(original; max_scans=4, topn=nothing, target_ms_order=UInt8(2), seed=42)
    @test 0 < length(filtered) <= 4

    # Index mapping and metadata consistency for a couple scans
    for i in 1:min(2, length(filtered))
        orig_idx = getOriginalScanIndex(filtered, i)
        @test getScanNumber(filtered, i) == getScanNumber(original, orig_idx)
        @test getRetentionTime(filtered, i) ≈ getRetentionTime(original, orig_idx)
        @test getMsOrder(filtered, i) == getMsOrder(original, orig_idx)
    end

    # Append a couple more scans
    prev_len = length(filtered)
    added = append!(filtered; max_additional_scans=2)
    @test length(filtered) == prev_len + added
    all_indices = getOriginalScanIndices(filtered)
    @test length(unique(all_indices)) == length(all_indices)

    # TopN filtering path
    filtered_topn = FilteredMassSpecData(original; max_scans=4, topn=2, target_ms_order=UInt8(2))
    @test length(filtered_topn) > 0
    for i in 1:length(filtered_topn)
        @test length(getMzArray(filtered_topn, i)) <= 2
    end

    # Edge cases
    empty_filtered = FilteredMassSpecData(original; max_scans=0, topn=nothing, target_ms_order=UInt8(2))
    @test length(empty_filtered) == 0

    # Priority selection helpers
    bins, rtmin, rtmax, binw = compute_rt_bins(original, 5)
    @test length(bins) == length(original)
    @test binw >= 0
    scans, starts, ends = sort_scans_by_peak_density(original, UInt8(2), bins, 5)
    prio = create_priority_order(scans, starts, ends)
    @test length(prio) == length(scans)
end
