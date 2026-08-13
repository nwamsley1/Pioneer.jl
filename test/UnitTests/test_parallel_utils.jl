# Unit tests for parallel utilities:
#   - parallel_foreach!  (src/utils/parallelUtils.jl)
#   - _distribute_scans_to_threads, _sort_scans_by_mz_in_rt_bins!  (partitionThreadTasks.jl)
#
# Run standalone: julia --project=. test/UnitTests/test_parallel_utils.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
end

# ═══════════════════════════════════════════════════════════════════════════════
# parallel_foreach!
# ═══════════════════════════════════════════════════════════════════════════════

@testset "parallel_foreach!" begin

@testset "writes to every index" begin
    n = 1000
    results = zeros(Int, n)
    Pioneer.parallel_foreach!(n) do chunk
        for i in chunk
            results[i] = i * 2
        end
    end
    @test results == [2i for i in 1:n]
end

@testset "returns nothing" begin
    ret = Pioneer.parallel_foreach!(10) do chunk
        for _ in chunk; end
    end
    @test ret === nothing
end

@testset "n = 0 does nothing" begin
    called = Ref(false)
    Pioneer.parallel_foreach!(0) do chunk
        called[] = true
    end
    @test called[] == false
end

@testset "n = 1 processes single element" begin
    results = zeros(Int, 1)
    Pioneer.parallel_foreach!(1) do chunk
        for i in chunk
            results[i] = 42
        end
    end
    @test results == [42]
end

@testset "custom tasks_per_thread" begin
    n = 100
    results = zeros(Int, n)
    Pioneer.parallel_foreach!(n; tasks_per_thread=1) do chunk
        for i in chunk
            results[i] = i
        end
    end
    @test results == collect(1:n)
end

@testset "per-chunk local state (thread-local pattern)" begin
    n = 200
    results = zeros(Float64, n)
    Pioneer.parallel_foreach!(n) do chunk
        # Each chunk gets its own accumulator — simulates thread-local buffers
        local_buf = zeros(Float64, 3)
        for i in chunk
            local_buf[1] = Float64(i)
            results[i] = local_buf[1]^2
        end
    end
    @test results == [Float64(i)^2 for i in 1:n]
end

@testset "large n with default settings" begin
    n = 100_000
    results = zeros(Int, n)
    Pioneer.parallel_foreach!(n) do chunk
        for i in chunk
            results[i] = 1
        end
    end
    @test sum(results) == n
end

end # parallel_foreach!


# ═══════════════════════════════════════════════════════════════════════════════
# _distribute_scans_to_threads
# ═══════════════════════════════════════════════════════════════════════════════

@testset "_distribute_scans_to_threads" begin

@testset "empty input returns empty tasks for each thread" begin
    tasks = Pioneer._distribute_scans_to_threads(Int64[], 4)
    @test length(tasks) == 4
    for (i, task) in enumerate(tasks)
        @test task[1] == i
        @test isempty(task[2])
    end
end

@testset "all scan indices are distributed exactly once" begin
    scan_indices = Int64[10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
                         110, 120, 130, 140, 150]
    tasks = Pioneer._distribute_scans_to_threads(scan_indices, 3)

    @test length(tasks) == 3

    # Collect all distributed indices
    all_distributed = Int64[]
    for task in tasks
        append!(all_distributed, task[2])
    end
    sort!(all_distributed)
    @test all_distributed == sort(scan_indices)
end

@testset "single thread gets all scans" begin
    scan_indices = Int64[1, 2, 3, 4, 5]
    tasks = Pioneer._distribute_scans_to_threads(scan_indices, 1)
    @test length(tasks) == 1
    @test sort(tasks[1][2]) == scan_indices
end

@testset "more threads than scans" begin
    scan_indices = Int64[42, 99]
    tasks = Pioneer._distribute_scans_to_threads(scan_indices, 8)
    @test length(tasks) == 8

    all_distributed = Int64[]
    for task in tasks
        append!(all_distributed, task[2])
    end
    @test sort(all_distributed) == sort(scan_indices)
end

@testset "no trailing zeros in task arrays" begin
    scan_indices = Int64[1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    tasks = Pioneer._distribute_scans_to_threads(scan_indices, 4)
    for task in tasks
        @test all(x -> x != 0, task[2])
    end
end

@testset "thread IDs are 1:n_threads" begin
    tasks = Pioneer._distribute_scans_to_threads(Int64[1, 2, 3], 5)
    @test [task[1] for task in tasks] == [1, 2, 3, 4, 5]
end

@testset "large input distributes evenly-ish" begin
    n = 1000
    scan_indices = Int64.(1:n)
    n_threads = 4
    tasks = Pioneer._distribute_scans_to_threads(scan_indices, n_threads)

    counts = [length(task[2]) for task in tasks]
    # Each thread should get roughly n/n_threads scans (within 2x)
    expected = n ÷ n_threads
    for c in counts
        @test c > 0
        @test c <= 2 * expected
    end
    @test sum(counts) == n
end

end # _distribute_scans_to_threads


# ═══════════════════════════════════════════════════════════════════════════════
# _sort_scans_by_mz_in_rt_bins!
# ═══════════════════════════════════════════════════════════════════════════════

@testset "_sort_scans_by_mz_in_rt_bins!" begin

@testset "single RT bin sorts by m/z" begin
    # All scans within 1 minute of RT → single bin → sorted by m/z
    scan_indices = Int64[1, 2, 3, 4]
    rt = Float32[0.1, 0.2, 0.3, 0.4]
    prec_mz = Union{Missing, Float32}[400.0f0, 200.0f0, 300.0f0, 100.0f0]

    Pioneer._sort_scans_by_mz_in_rt_bins!(scan_indices, rt, prec_mz)
    # Should be sorted by prec_mz[scan_idx]: 100, 200, 300, 400
    @test scan_indices == [4, 2, 3, 1]
end

@testset "multiple RT bins sort independently" begin
    # Two RT bins: scans 1-2 (RT 0.0, 0.5) and scans 3-4 (RT 2.0, 2.5)
    scan_indices = Int64[1, 2, 3, 4]
    rt = Float32[0.0, 0.5, 2.0, 2.5]
    prec_mz = Union{Missing, Float32}[500.0f0, 100.0f0, 400.0f0, 200.0f0]

    Pioneer._sort_scans_by_mz_in_rt_bins!(scan_indices, rt, prec_mz)
    # Bin 1 (RT < 1.0): scans 1,2 sorted by mz → [2, 1] (100, 500)
    # Bin 2 (RT > 1.0): scans 3,4 sorted by mz → [4, 3] (200, 400)
    @test scan_indices == [2, 1, 4, 3]
end

@testset "empty and single-element inputs" begin
    empty = Int64[]
    rt = Float32[]
    mz = Union{Missing, Float32}[]
    Pioneer._sort_scans_by_mz_in_rt_bins!(empty, rt, mz)
    @test isempty(empty)

    single = Int64[1]
    rt1 = Float32[1.0]
    mz1 = Union{Missing, Float32}[500.0f0]
    Pioneer._sort_scans_by_mz_in_rt_bins!(single, rt1, mz1)
    @test single == [1]
end

@testset "missing m/z treated as zero" begin
    scan_indices = Int64[1, 2, 3]
    rt = Float32[0.1, 0.2, 0.3]
    prec_mz = Union{Missing, Float32}[300.0f0, missing, 100.0f0]

    Pioneer._sort_scans_by_mz_in_rt_bins!(scan_indices, rt, prec_mz)
    # missing → 0.0, so order by mz: 0(scan2), 100(scan3), 300(scan1)
    @test scan_indices == [2, 3, 1]
end

end # _sort_scans_by_mz_in_rt_bins!


# ═══════════════════════════════════════════════════════════════════════════════
# _round_float32
# ═══════════════════════════════════════════════════════════════════════════════

@testset "_round_float32" begin
    @test Pioneer._round_float32(1.23456f0, 2) ≈ 1.23f0
    @test Pioneer._round_float32(1.23456f0, 4) ≈ 1.2346f0
    @test Pioneer._round_float32(missing, 6) == 0.0f0
    @test Pioneer._round_float32(0.0f0, 3) == 0.0f0
end
