using Test
using Arrow
using DataFrames
using Random

# Import only what we need from Pioneer (avoid Pioneer.func calls)
using Pioneer: BasicMassSpecData,
               getScanHeader, getScanNumber, getBasePeakMz, getBasePeakIntensity,
               getRetentionTime, getLowMz, getHighMz, getTIC, getMsOrder,
               getMzArray, getIntensityArray, getCenterMz, getIsolationWidthMz,
               getCenterMzs, getIsolationWidthMzs, getRetentionTimes, getTICs, getMsOrders,
               FilteredMassSpecData, IndexedMassSpecData, create_scan_mapping,
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

@testset "IndexedMassSpecData Tests" begin
    # Create temp Arrow file
    temp_dir = mktempdir()
    file_path = joinpath(temp_dir, "ms.arrow")
    write_basic_ms_arrow(file_path, n=10)

    # Load as Basic
    original = BasicMassSpecData(file_path)
    @test length(original) == 10

    # Test IndexedMassSpecData construction
    selected_indices = Int32[2, 4, 6, 8]
    indexed = IndexedMassSpecData(original, selected_indices)

    # Basic properties
    @test length(indexed) == 4
    @test indexed.n_scans == 4
    @test indexed.scan_indices == selected_indices

    # Test virtual-to-actual index mapping
    @test getOriginalScanIndex(indexed, 1) == 2  # Virtual index 1 → actual index 2
    @test getOriginalScanIndex(indexed, 2) == 4  # Virtual index 2 → actual index 4
    @test getOriginalScanIndex(indexed, 3) == 6  # Virtual index 3 → actual index 6
    @test getOriginalScanIndex(indexed, 4) == 8  # Virtual index 4 → actual index 8

    # Test batch mapping
    all_indices = getOriginalScanIndices(indexed)
    @test all_indices == selected_indices

    # Test metadata getters map correctly
    for i in 1:length(indexed)
        orig_idx = getOriginalScanIndex(indexed, i)
        @test getScanHeader(indexed, i) == getScanHeader(original, orig_idx)
        @test getScanNumber(indexed, i) == getScanNumber(original, orig_idx)
        @test getBasePeakMz(indexed, i) ≈ getBasePeakMz(original, orig_idx)
        @test getBasePeakIntensity(indexed, i) ≈ getBasePeakIntensity(original, orig_idx)
        @test getRetentionTime(indexed, i) ≈ getRetentionTime(original, orig_idx)
        @test getLowMz(indexed, i) ≈ getLowMz(original, orig_idx)
        @test getHighMz(indexed, i) ≈ getHighMz(original, orig_idx)
        @test getTIC(indexed, i) ≈ getTIC(original, orig_idx)
        @test getMsOrder(indexed, i) == getMsOrder(original, orig_idx)
    end

    # Test array getters
    for i in 1:length(indexed)
        orig_idx = getOriginalScanIndex(indexed, i)
        mz_indexed = getMzArray(indexed, i)
        mz_original = getMzArray(original, orig_idx)
        int_indexed = getIntensityArray(indexed, i)
        int_original = getIntensityArray(original, orig_idx)

        @test length(mz_indexed) == length(mz_original)
        @test length(int_indexed) == length(int_original)
        @test all(mz_indexed .≈ mz_original)
        @test all(int_indexed .≈ int_original)
    end

    # Test batch getters
    center_mzs = getCenterMzs(indexed)
    isolation_widths = getIsolationWidthMzs(indexed)
    retention_times = getRetentionTimes(indexed)
    tics = getTICs(indexed)
    ms_orders = getMsOrders(indexed)

    @test length(center_mzs) == length(indexed)
    @test length(isolation_widths) == length(indexed)
    @test length(retention_times) == length(indexed)
    @test length(tics) == length(indexed)
    @test length(ms_orders) == length(indexed)

    # Verify batch getters match individual getters
    for i in 1:length(indexed)
        @test getCenterMz(indexed, i) ≈ center_mzs[i]
        @test getIsolationWidthMz(indexed, i) ≈ isolation_widths[i]
        @test getRetentionTime(indexed, i) ≈ retention_times[i]
        @test getTIC(indexed, i) ≈ tics[i]
        @test getMsOrder(indexed, i) == ms_orders[i]
    end

    # Test create_scan_mapping function
    mapping = create_scan_mapping(indexed)
    @test length(mapping) == length(indexed)
    for i in 1:length(indexed)
        @test mapping[Int32(i)] == selected_indices[i]
    end

    # Edge cases

    # Empty IndexedMassSpecData
    empty_indices = Int32[]
    empty_indexed = IndexedMassSpecData(original, empty_indices)
    @test length(empty_indexed) == 0
    @test empty_indexed.n_scans == 0
    @test isempty(getOriginalScanIndices(empty_indexed))

    # Single scan
    single_indices = Int32[5]
    single_indexed = IndexedMassSpecData(original, single_indices)
    @test length(single_indexed) == 1
    @test getOriginalScanIndex(single_indexed, 1) == 5
    @test getScanNumber(single_indexed, 1) == getScanNumber(original, 5)

    # All scans (identity mapping)
    all_scan_indices = Int32.(1:length(original))
    all_indexed = IndexedMassSpecData(original, all_scan_indices)
    @test length(all_indexed) == length(original)
    for i in 1:length(all_indexed)
        @test getOriginalScanIndex(all_indexed, i) == i
        @test getScanNumber(all_indexed, i) == getScanNumber(original, i)
    end

    # Constructor validation tests
    @test_throws BoundsError IndexedMassSpecData(original, Int32[0])  # Invalid index
    @test_throws BoundsError IndexedMassSpecData(original, Int32[11]) # Out of bounds
    @test_throws BoundsError IndexedMassSpecData(original, Int32[-1]) # Negative index

    # Test with duplicated indices (should work - represents repeated sampling)
    duplicate_indices = Int32[3, 3, 7, 7]
    dup_indexed = IndexedMassSpecData(original, duplicate_indices)
    @test length(dup_indexed) == 4
    @test getOriginalScanIndex(dup_indexed, 1) == 3
    @test getOriginalScanIndex(dup_indexed, 2) == 3
    @test getOriginalScanIndex(dup_indexed, 3) == 7
    @test getOriginalScanIndex(dup_indexed, 4) == 7

    # Test with unsorted indices (should work - maintains order)
    unsorted_indices = Int32[9, 2, 7, 3]
    unsorted_indexed = IndexedMassSpecData(original, unsorted_indices)
    @test length(unsorted_indexed) == 4
    @test getOriginalScanIndex(unsorted_indexed, 1) == 9
    @test getOriginalScanIndex(unsorted_indexed, 2) == 2
    @test getOriginalScanIndex(unsorted_indexed, 3) == 7
    @test getOriginalScanIndex(unsorted_indexed, 4) == 3

    # Verify scan order preservation
    @test getScanNumber(unsorted_indexed, 1) == getScanNumber(original, 9)
    @test getScanNumber(unsorted_indexed, 2) == getScanNumber(original, 2)
    @test getScanNumber(unsorted_indexed, 3) == getScanNumber(original, 7)
    @test getScanNumber(unsorted_indexed, 4) == getScanNumber(original, 3)
end
