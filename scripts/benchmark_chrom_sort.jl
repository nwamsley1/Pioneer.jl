#!/usr/bin/env julia
#
# Benchmark chromatogram sorting strategies
#
# Usage: julia scripts/benchmark_chrom_sort.jl path/to/unsorted_chroms_ms2.arrow
#
# Tests 6 approaches for sorting a ~22M-row chromatogram DataFrame by
# [:precursor_idx, :rt] and reports wall time + allocations for each.

using Arrow, DataFrames, SortingAlgorithms

if isempty(ARGS)
    error("Usage: julia scripts/benchmark_chrom_sort.jl <path_to_unsorted_chroms_ms2.arrow>")
end

arrow_path = ARGS[1]
@assert isfile(arrow_path) "File not found: $arrow_path"

println("Loading Arrow table from $arrow_path ...")
df_orig = DataFrame(Arrow.Table(arrow_path))
println("  $(nrow(df_orig)) rows, $(ncol(df_orig)) columns: $(names(df_orig))")
println("  Column types: ", join(["$c::$(eltype(df_orig[!,c]))" for c in names(df_orig)], ", "))
println()

const N_TRIALS = 5

# ──────────────────────────────────────────────────────────────────────
# Baseline: sort a copy and keep for verification
# ──────────────────────────────────────────────────────────────────────
println("Computing baseline sort (MergeSort) for verification...")
baseline = copy(df_orig)
sort!(baseline, [:precursor_idx, :rt])
println("  Done.\n")

function verify_match(sorted_df::DataFrame, label::AbstractString)
    if sorted_df.precursor_idx == baseline.precursor_idx &&
       sorted_df.rt == baseline.rt &&
       sorted_df.intensity == baseline.intensity &&
       sorted_df.scan_idx == baseline.scan_idx
        println("  ✓ $label matches baseline")
        return true
    else
        mismatches = sum(
            (sorted_df.precursor_idx .!= baseline.precursor_idx) .|
            (sorted_df.rt .!= baseline.rt)
        )
        println("  ✗ $label MISMATCH — $mismatches rows differ")
        return false
    end
end

function fresh_copy()
    # Rebuild from column vectors to get a pristine unsorted DataFrame
    DataFrame(
        :precursor_idx => copy(df_orig.precursor_idx),
        :rt => copy(df_orig.rt),
        :intensity => copy(df_orig.intensity),
        :scan_idx => copy(df_orig.scan_idx),
    )
end

function run_benchmark(sort_fn!::Function, label::String; verify=true)
    println("=" ^ 70)
    println(label)
    println("=" ^ 70)

    # Warmup
    df_warmup = fresh_copy()
    sort_fn!(df_warmup)

    # Timed trials
    times = Float64[]
    allocs = Int[]
    local df_final
    for trial in 1:N_TRIALS
        df = fresh_copy()
        stats = @timed sort_fn!(df)
        push!(times, stats.time)
        push!(allocs, stats.bytes)
        df_final = df
    end

    median_t = sort(times)[div(N_TRIALS, 2) + 1]
    min_t = minimum(times)
    mean_alloc = round(Int, sum(allocs) / N_TRIALS)

    println("  Times:  ", join(["$(round(t, digits=3))s" for t in times], ", "))
    println("  Median: $(round(median_t, digits=3))s | Min: $(round(min_t, digits=3))s | Alloc: $(round(mean_alloc / 1e6, digits=1)) MB")

    if verify
        verify_match(df_final, label)
    end
    println()
    return (median=median_t, min=min_t, alloc_mb=mean_alloc / 1e6)
end

# ──────────────────────────────────────────────────────────────────────
# 1. Baseline: sort!(df, [:precursor_idx, :rt])  (MergeSort default)
# ──────────────────────────────────────────────────────────────────────
r1 = run_benchmark("1. sort!(df, [:precursor_idx, :rt])  — MergeSort (baseline)") do df
    sort!(df, [:precursor_idx, :rt])
end

# ──────────────────────────────────────────────────────────────────────
# 2. sort!(df, [:precursor_idx, :rt], alg=QuickSort)
# ──────────────────────────────────────────────────────────────────────
r2 = run_benchmark("2. sort!(df, [:precursor_idx, :rt], alg=QuickSort)") do df
    sort!(df, [:precursor_idx, :rt], alg=QuickSort)
end

# ──────────────────────────────────────────────────────────────────────
# 3. Manual sortperm on tuple column
# ──────────────────────────────────────────────────────────────────────
r3 = run_benchmark("3. Manual sortperm on zipped columns") do df
    keys = collect(zip(df.precursor_idx, df.rt))
    perm = sortperm(keys)
    for col in names(df)
        df[!, col] = df[!, col][perm]
    end
end

# ──────────────────────────────────────────────────────────────────────
# 4. RadixSort on packed UInt64 key (float → sortable uint)
# ──────────────────────────────────────────────────────────────────────
function float32_to_sortable_uint32(x::Float32)::UInt32
    u = reinterpret(UInt32, x)
    mask = ifelse(u & 0x80000000 != 0, 0xffffffff, 0x80000000)
    return xor(u, mask)
end

function pack_sort_key(precursor_idx::UInt32, rt::Float32)::UInt64
    return (UInt64(precursor_idx) << 32) | UInt64(float32_to_sortable_uint32(rt))
end

r4 = run_benchmark("4. RadixSort on packed UInt64 key (precursor_idx << 32 | float2sortable(rt))") do df
    n = nrow(df)
    keys = Vector{UInt64}(undef, n)
    @inbounds for i in 1:n
        keys[i] = pack_sort_key(df.precursor_idx[i], df.rt[i])
    end
    perm = sortperm(keys, alg=RadixSort)
    for col in names(df)
        df[!, col] = df[!, col][perm]
    end
end

# ──────────────────────────────────────────────────────────────────────
# 5. Counting sort on precursor_idx + within-group sort by rt
# ──────────────────────────────────────────────────────────────────────
r5 = run_benchmark("5. Counting sort precursor_idx + within-group sort by rt") do df
    n = nrow(df)
    pidx = df.precursor_idx
    max_pid = maximum(pidx)

    # Count occurrences of each precursor_idx
    counts = zeros(Int, max_pid + 1)
    @inbounds for i in 1:n
        counts[pidx[i] + 1] += 1
    end

    # Compute offsets (exclusive prefix sum)
    offsets = Vector{Int}(undef, length(counts) + 1)
    offsets[1] = 1
    @inbounds for i in eachindex(counts)
        offsets[i + 1] = offsets[i] + counts[i]
    end

    # Place row indices into bins by precursor_idx
    perm = Vector{Int}(undef, n)
    pos = copy(offsets)
    @inbounds for i in 1:n
        bin = pidx[i] + 1
        perm[pos[bin]] = i
        pos[bin] += 1
    end

    # Sort within each group by rt
    rt = df.rt
    @inbounds for bin in eachindex(counts)
        lo = offsets[bin]
        hi = offsets[bin + 1] - 1
        if hi > lo
            sort!(@view(perm[lo:hi]), by=j -> rt[j])
        end
    end

    # Apply permutation
    for col in names(df)
        df[!, col] = df[!, col][perm]
    end
end

# ──────────────────────────────────────────────────────────────────────
# 6. RadixSort on packed UInt64 key with integer RT (milliminutes)
# ──────────────────────────────────────────────────────────────────────
# If the file has rt_milliminutes, use it directly. Otherwise compute
# it on the fly from the Float32 rt column (simulating what the pipeline
# would do if we stored RT as integer).
has_int_rt = hasproperty(df_orig, :rt_milliminutes)
if !has_int_rt
    println("Note: rt_milliminutes column not in file — computing from rt * 1000\n")
end

r6 = run_benchmark("6. RadixSort on packed UInt64 key with integer RT (milliminutes)") do df
    n = nrow(df)
    keys = Vector{UInt64}(undef, n)
    if has_int_rt
        @inbounds for i in 1:n
            keys[i] = (UInt64(df.precursor_idx[i]) << 32) | UInt64(df.rt_milliminutes[i])
        end
    else
        @inbounds for i in 1:n
            rt_int = round(UInt32, df.rt[i] * 1000)
            keys[i] = (UInt64(df.precursor_idx[i]) << 32) | UInt64(rt_int)
        end
    end
    perm = sortperm(keys, alg=RadixSort)
    for col in names(df)
        df[!, col] = df[!, col][perm]
    end
end

# ──────────────────────────────────────────────────────────────────────
# Summary table
# ──────────────────────────────────────────────────────────────────────
println("=" ^ 70)
println("SUMMARY  ($(nrow(df_orig)) rows, $N_TRIALS trials each)")
println("=" ^ 70)
results = [
    ("1. MergeSort (baseline)", r1),
    ("2. QuickSort", r2),
    ("3. Zip sortperm", r3),
    ("4. RadixSort float key", r4),
    ("5. Counting + group sort", r5),
    ("6. RadixSort int key", r6),
]
println(rpad("Method", 30), rpad("Median", 10), rpad("Min", 10), "Alloc (MB)")
println("-" ^ 60)
for (name, r) in results
    println(
        rpad(name, 30),
        rpad("$(round(r.median, digits=3))s", 10),
        rpad("$(round(r.min, digits=3))s", 10),
        "$(round(r.alloc_mb, digits=1))"
    )
end
println()
best = argmin([r.median for (_, r) in results])
println("Winner: $(results[best][1])")
