#!/usr/bin/env julia
# Microbenchmark: countingSortByCol! vs specialsort! on realistic SparseArray data
# Usage: julia --project=. scripts/benchmark_sparse_sort.jl

using Pioneer
using Random, Statistics, Printf

# Access internals
const SA = Pioneer.SparseArray
const specialsort! = Pioneer.specialsort!
const countingSortByCol! = Pioneer.countingSortByCol!

"""
Fill a SparseArray with realistic data mimicking buildDesignMatrix! output.
  n_entries: total non-zero entries (typically 100-500 per scan)
  n_cols: number of unique precursors (typically 5-50)
  n_rows: number of unique peaks (typically 50-200)
"""
function fill_sparse!(sa::SA{UInt32, Float32}, n_entries::Int, n_cols::Int, n_rows::Int)
    sa.n_vals = n_entries
    for i in 1:n_entries
        sa.colval[i] = UInt16(rand(1:n_cols))
        sa.rowval[i] = UInt32(rand(1:n_rows))
        sa.nzval[i] = rand(Float32)
        sa.x[i] = rand(Float32)
        sa.matched[i] = rand(Bool)
        sa.isotope[i] = UInt8(rand(0:2))
    end
end

function reset_sparse!(sa::SA{UInt32, Float32})
    for i in 1:sa.n_vals
        sa.colval[i] = zero(UInt16)
        sa.rowval[i] = zero(UInt32)
        sa.nzval[i] = zero(Float32)
        sa.x[i] = zero(Float32)
        sa.matched[i] = true
        sa.colptr[i] = zero(UInt32)
        sa.isotope[i] = zero(UInt8)
    end
    sa.n_vals = 0
    sa.m = 0
    sa.n = 0
end

function benchmark_old_sort!(sa, n_entries, n_cols, n_rows, n_iters)
    times = zeros(n_iters)
    for iter in 1:n_iters
        fill_sparse!(sa, n_entries, n_cols, n_rows)
        t = @elapsed begin
            specialsort!(sa, 1, sa.n_vals, Base.Order.Forward)
            # Build colptr (part of old sortSparse!)
            max_col = 1
            max_row = 1
            sa.colptr[1] = 1
            for i in 1:(sa.n_vals - 1)
                if sa.rowval[i + 1] > max_row; max_row = sa.rowval[i + 1]; end
                if sa.colval[i + 1] != sa.colval[i]
                    max_col += 1
                    sa.colptr[max_col] = i + 1
                end
            end
            sa.colptr[max_col + 1] = sa.n_vals
            sa.n = max_col
            sa.m = max_row
        end
        times[iter] = t
        reset_sparse!(sa)
    end
    return times
end

function benchmark_counting_sort!(sa, n_entries, n_cols, n_rows, n_iters)
    times = zeros(n_iters)
    for iter in 1:n_iters
        fill_sparse!(sa, n_entries, n_cols, n_rows)
        t = @elapsed countingSortByCol!(sa)
        times[iter] = t
        reset_sparse!(sa)
    end
    return times
end

function run_benchmark(n_entries, n_cols, n_rows; n_iters=5000, n_warmup=500)
    # Allocate with enough capacity for 2x (counting sort needs slack)
    sa = SA(UInt32(max(n_entries * 3, 5000)))

    # Warmup
    benchmark_old_sort!(sa, n_entries, n_cols, n_rows, n_warmup)
    benchmark_counting_sort!(sa, n_entries, n_cols, n_rows, n_warmup)

    # Benchmark
    old_times = benchmark_old_sort!(sa, n_entries, n_cols, n_rows, n_iters)
    new_times = benchmark_counting_sort!(sa, n_entries, n_cols, n_rows, n_iters)

    old_median = median(old_times) * 1e6  # μs
    new_median = median(new_times) * 1e6
    speedup = old_median / new_median

    Printf = Base.Libc  # just for formatting
    @printf("  n=%4d cols=%3d rows=%3d | old %7.1f μs | new %7.1f μs | speedup %.2fx\n",
            n_entries, n_cols, n_rows, old_median, new_median, speedup)
    return (old_median, new_median, speedup)
end

println("=" ^ 72)
println("SparseArray sort microbenchmark: specialsort! vs countingSortByCol!")
println("  5000 iterations per config, median times")
println("=" ^ 72)

# Realistic sizes from profiling (typical per-scan values)
configs = [
    (50,   5,   20),    # small scan
    (100,  10,  40),    # typical scan
    (200,  20,  80),    # medium scan
    (500,  30,  150),   # large scan
    (1000, 50,  300),   # very large scan
    (2000, 50,  500),   # extreme scan
]

results = []
for (n, c, r) in configs
    push!(results, run_benchmark(n, c, r))
end

println("=" ^ 72)
avg_speedup = mean([r[3] for r in results])
@printf("Average speedup: %.2fx\n", avg_speedup)
println("=" ^ 72)
