#!/usr/bin/env julia
#
# K-Way Merge Benchmark for Transition Sorting
#
# Compares full QuickSort (current) vs k-way merge (proposed) for sorting
# transitions by m/z after filling from multiple precursors.
#
# Usage: julia --project=. scripts/kway_merge_benchmark.jl
#

using Pioneer
using Pioneer: DetailedFrag, getMZ, getFragments, getPrecFragRange,
               StandardFragmentLookup, load_detailed_frags, deserialize_from_jls
using BenchmarkTools
using DataStructures: MutableBinaryMinHeap
using Printf
using Random
using Statistics

const NEUTRON_F32 = Float32(1.00335)

# ============================================================================
# Library loading helper
# ============================================================================

"""
Load a StandardFragmentLookup directly from a .poin directory.
Returns (frags, prec_frag_ranges) or nothing if files are missing.
"""
function load_fragment_lookup(lib_dir::String)
    frags_path = joinpath(lib_dir, "detailed_fragments.jls")
    ranges_path = joinpath(lib_dir, "precursor_to_fragment_indices.jls")

    if !isfile(frags_path) || !isfile(ranges_path)
        return nothing
    end

    frags = load_detailed_frags(frags_path)
    ranges = deserialize_from_jls(ranges_path)
    return StandardFragmentLookup(frags, ranges)
end

# ============================================================================
# Part 1: Verify the assumption — library fragments are m/z-sorted per precursor
# ============================================================================

function verify_library_sorting(lib_dir::String)
    println("=" ^ 70)
    println("Part 1: Verify library fragment m/z ordering")
    println("=" ^ 70)

    lookup = load_fragment_lookup(lib_dir)
    if lookup === nothing
        println("  Skipped: library files not found in $lib_dir")
        println("  (Run the ecoli test first to build the library)")
        println()
        return nothing
    end

    frags = getFragments(lookup)
    n_precs = length(lookup.prec_frag_ranges) - 1

    n_sorted = 0
    n_unsorted = 0
    n_empty = 0
    violation_gaps = Float32[]

    for prec_i in 1:n_precs
        frag_range = getPrecFragRange(lookup, UInt32(prec_i))
        isempty(frag_range) && (n_empty += 1; continue)

        sorted = true
        for j in (first(frag_range)+1):last(frag_range)
            if frags[j].mz < frags[j-1].mz
                sorted = false
                push!(violation_gaps, frags[j-1].mz - frags[j].mz)
            end
        end
        sorted ? (n_sorted += 1) : (n_unsorted += 1)
    end

    pct_sorted = 100.0 * n_sorted / (n_sorted + n_unsorted)
    println("  Total precursors:    $(n_sorted + n_unsorted + n_empty)")
    println("  Non-empty:           $(n_sorted + n_unsorted)")
    println("  m/z-sorted:          $n_sorted ($(@sprintf("%.1f", pct_sorted))%)")
    println("  Unsorted:            $n_unsorted")
    if !isempty(violation_gaps)
        println("  Violation gaps:      median=$(@sprintf("%.4f", median(violation_gaps))) Da, " *
                "max=$(@sprintf("%.4f", maximum(violation_gaps))) Da")
    end
    println()

    # Now check: if we interleave isotopes per fragment, how often does m/z order break?
    println("  Isotope interleaving analysis (n_iso=2, charge=1 and 2):")
    for charge in UInt8[1, 2]
        n_broken = 0
        n_checked = 0
        for prec_i in 1:min(n_precs, 50000)
            frag_range = getPrecFragRange(lookup, UInt32(prec_i))
            length(frag_range) < 2 && continue
            for j in (first(frag_range)+1):last(frag_range)
                n_checked += 1
                frag_prev = frags[j-1]
                frag_curr = frags[j]
                # After interleaving: prev_M0, prev_M1, curr_M0, curr_M1
                # Violation if curr_M0 < prev_M1
                prev_m1 = frag_prev.mz + NEUTRON_F32 / charge
                if frag_curr.mz < prev_m1
                    n_broken += 1
                end
            end
        end
        pct = n_checked > 0 ? 100.0 * n_broken / n_checked : 0.0
        println("    charge=$charge: $n_broken / $n_checked consecutive pairs broken ($(@sprintf("%.1f", pct))%)")
    end
    println()

    return lookup
end

# ============================================================================
# Part 2: K-way merge implementation
# ============================================================================

"""
    HeapEntry

Entry for the min-heap used in k-way merge. Stores value + stream index.
"""
struct HeapEntry
    mz::Float32
    src_idx::Int32   # index into the flat source array
end
Base.isless(a::HeapEntry, b::HeapEntry) = a.mz < b.mz

"""
    kway_merge!(output, input, stream_ranges)

Merge K sorted streams (defined by `stream_ranges` into `input`) into `output`.
Both `input` and `output` are vectors of `DetailedFrag{Float32}`.
Returns the number of elements written to output.
"""
function kway_merge!(
    output::Vector{DetailedFrag{Float32}},
    input::Vector{DetailedFrag{Float32}},
    stream_ranges::Vector{UnitRange{Int}},
)
    K = length(stream_ranges)
    K == 0 && return 0

    # Fast path: single stream — just copy
    if K == 1
        r = stream_ranges[1]
        n = length(r)
        copyto!(output, 1, input, first(r), n)
        return n
    end

    # Initialize heap with first element of each non-empty stream
    # Also track current position within each stream
    cursors = Vector{Int}(undef, K)  # current position in each stream
    ends    = Vector{Int}(undef, K)  # end position of each stream

    # Build initial heap entries
    heap = MutableBinaryMinHeap{HeapEntry}()
    sizehint!(heap, K)

    for k in 1:K
        r = stream_ranges[k]
        if !isempty(r)
            cursors[k] = first(r)
            ends[k] = last(r)
            push!(heap, HeapEntry(input[first(r)].mz, Int32(k)))
        else
            cursors[k] = 0
            ends[k] = -1
        end
    end

    out_idx = 0
    while !isempty(heap)
        entry = pop!(heap)
        k = Int(entry.src_idx)
        pos = cursors[k]

        out_idx += 1
        @inbounds output[out_idx] = input[pos]

        # Advance cursor for this stream
        next_pos = pos + 1
        if next_pos <= ends[k]
            cursors[k] = next_pos
            push!(heap, HeapEntry(@inbounds(input[next_pos].mz), Int32(k)))
        end
    end

    return out_idx
end

"""
    kway_merge_inplace!(data, n, stream_ranges, buffer)

In-place variant: merges streams within `data[1:n]` using `buffer` as scratch,
then copies result back to `data[1:n]`.
"""
function kway_merge_inplace!(
    data::Vector{DetailedFrag{Float32}},
    n::Int,
    stream_ranges::Vector{UnitRange{Int}},
    buffer::Vector{DetailedFrag{Float32}},
)
    m = kway_merge!(buffer, data, stream_ranges)
    @assert m == n "Merge output $m != expected $n"
    copyto!(data, 1, buffer, 1, n)
    return n
end

# ============================================================================
# Part 3: Generate realistic test data & benchmark
# ============================================================================

"""
Create a synthetic transition array mimicking real fill patterns.
Returns (transitions, n_transitions, stream_ranges) where stream_ranges
describes the K sorted streams if we fill isotopes in separate passes.
"""
function generate_test_transitions(;
    n_precursors::Int = 200,
    frags_per_prec::Int = 10,
    n_isotopes::Int = 2,
    seed::Int = 42,
)
    N = n_precursors * frags_per_prec * n_isotopes
    transitions = Vector{DetailedFrag{Float32}}(undef, N + 1000)  # some slack

    # Pre-generate per-precursor properties using deterministic per-precursor RNGs
    # This ensures both fill strategies produce identical fragment data
    prec_charges = Vector{UInt8}(undef, n_precursors)
    prec_frag_mzs = Vector{Vector{Float32}}(undef, n_precursors)
    for p in 1:n_precursors
        prng = MersenneTwister(seed + p)
        prec_charges[p] = rand(prng, UInt8[1, 2, 3])
        prec_frag_mzs[p] = sort(rand(prng, Float32, frags_per_prec) .* 500.0f0 .+ 100.0f0)
    end

    # --- Strategy A: Current interleaved fill (per-fragment isotopes) ---
    # This is what the current code does: for each precursor, for each fragment,
    # write M+0, M+1, M+2, ... then next fragment.
    idx = 0
    for p in 1:n_precursors
        charge = prec_charges[p]
        frag_offsets = prec_frag_mzs[p]
        for f in 1:frags_per_prec
            frag_mz = frag_offsets[f]
            for iso in 0:(n_isotopes-1)
                idx += 1
                mz = Float32(frag_mz + iso * NEUTRON_F32 / charge)
                transitions[idx] = DetailedFrag{Float32}(
                    UInt32(p), mz, Float16(1.0),
                    UInt16(0), true, false, false, iso > 0,
                    charge, UInt8(f), UInt8(charge), UInt8(f), UInt8(0)
                )
            end
        end
    end
    n_interleaved = idx

    # --- Strategy B: Separate isotope passes ---
    # Fill all M+0 for all precursors, then all M+1, etc.
    # Each (precursor, isotope_level) block is m/z-sorted.
    transitions_b = Vector{DetailedFrag{Float32}}(undef, N + 1000)
    stream_ranges = UnitRange{Int}[]
    idx = 0
    for iso in 0:(n_isotopes-1)
        for p in 1:n_precursors
            charge = prec_charges[p]
            frag_offsets = prec_frag_mzs[p]
            stream_start = idx + 1
            for f in 1:frags_per_prec
                frag_mz = frag_offsets[f]
                mz = Float32(frag_mz + iso * NEUTRON_F32 / charge)
                idx += 1
                transitions_b[idx] = DetailedFrag{Float32}(
                    UInt32(p), mz, Float16(1.0),
                    UInt16(0), true, false, false, iso > 0,
                    charge, UInt8(f), UInt8(charge), UInt8(f), UInt8(0)
                )
            end
            push!(stream_ranges, stream_start:idx)
        end
    end
    n_streamed = idx

    @assert n_interleaved == n_streamed "Mismatch: $n_interleaved vs $n_streamed"

    return (
        interleaved = transitions,
        streamed = transitions_b,
        n = n_interleaved,
        stream_ranges = stream_ranges,
        n_precursors = n_precursors,
        n_isotopes = n_isotopes,
    )
end

"""
Benchmark current approach (full sort) vs k-way merge.
"""
function run_benchmarks()
    println("=" ^ 70)
    println("Part 3: Benchmark — QuickSort vs K-way Merge")
    println("=" ^ 70)

    configs = [
        (n_prec=50,  frags=10, iso=2, label="Small  (50×10×2=1000)"),
        (n_prec=200, frags=10, iso=2, label="Medium (200×10×2=4000)"),
        (n_prec=200, frags=15, iso=2, label="Large  (200×15×2=6000)"),
        (n_prec=500, frags=10, iso=2, label="XLarge (500×10×2=10000)"),
        (n_prec=200, frags=10, iso=3, label="Med+3iso (200×10×3=6000)"),
        (n_prec=1000, frags=10, iso=2, label="Huge (1000×10×2=20000)"),
    ]

    for cfg in configs
        println("\n  Config: $(cfg.label)")
        data = generate_test_transitions(
            n_precursors=cfg.n_prec,
            frags_per_prec=cfg.frags,
            n_isotopes=cfg.iso,
        )
        N = data.n
        K = length(data.stream_ranges)

        # --- Approach 1: Current — PartialQuickSort on interleaved data ---
        work_a = copy(data.interleaved)
        sort_bench = @benchmark begin
            sort!(@view($work_a[1:$N]), by=x -> x.mz, alg=PartialQuickSort(1:$N))
        end setup=begin
            copyto!($work_a, 1, $(data.interleaved), 1, $N)
        end evals=1 samples=200

        # --- Approach 1b: Plain QuickSort (not partial) ---
        work_a2 = copy(data.interleaved)
        qsort_bench = @benchmark begin
            sort!(@view($work_a2[1:$N]), by=x -> x.mz, alg=QuickSort)
        end setup=begin
            copyto!($work_a2, 1, $(data.interleaved), 1, $N)
        end evals=1 samples=200

        # --- Approach 2: K-way merge on streamed data ---
        buffer = Vector{DetailedFrag{Float32}}(undef, N + 100)
        work_b = copy(data.streamed)
        sr = data.stream_ranges
        merge_bench = @benchmark begin
            kway_merge!($buffer, $work_b, $sr)
        end evals=1 samples=200

        # Verify correctness: both should produce same sorted order
        copyto!(work_a, 1, data.interleaved, 1, N)
        sort!(@view(work_a[1:N]), by=x -> x.mz, alg=PartialQuickSort(1:N))
        m = kway_merge!(buffer, data.streamed, data.stream_ranges)
        @assert m == N

        sorted_ok = true
        for i in 1:N
            if work_a[i].mz != buffer[i].mz
                sorted_ok = false
                break
            end
        end

        t_partial = median(sort_bench).time / 1e3  # μs
        t_quick   = median(qsort_bench).time / 1e3
        t_merge   = median(merge_bench).time / 1e3
        alloc_partial = median(sort_bench).allocs
        alloc_quick   = median(qsort_bench).allocs
        alloc_merge   = median(merge_bench).allocs

        @printf("    N=%d, K=%d streams\n", N, K)
        @printf("    PartialQuickSort: %8.1f μs  (%d allocs)\n", t_partial, alloc_partial)
        @printf("    QuickSort:        %8.1f μs  (%d allocs)\n", t_quick, alloc_quick)
        @printf("    K-way merge:      %8.1f μs  (%d allocs)\n", t_merge, alloc_merge)
        @printf("    Speedup (vs PQS): %5.2fx\n", t_partial / t_merge)
        @printf("    Speedup (vs QS):  %5.2fx\n", t_quick / t_merge)
        @printf("    Correct:          %s\n", sorted_ok ? "YES ✓" : "NO ✗")
    end
    println()
end

# ============================================================================
# Part 3b: Simpler merge approaches — tournament tree / manual heap
# ============================================================================

"""
    merge_sorted_streams_manual!(output, input, stream_starts, stream_ends, K)

A manual tournament-tree k-way merge without DataStructures.jl heap.
Uses a simple array-based min-heap for lower overhead at small K.
"""
function merge_sorted_streams_manual!(
    output::Vector{DetailedFrag{Float32}},
    input::Vector{DetailedFrag{Float32}},
    stream_starts::Vector{Int},
    stream_ends::Vector{Int},
    K::Int,
)
    # Heap entries: (mz, stream_index)
    heap_mz = Vector{Float32}(undef, K)
    heap_k  = Vector{Int}(undef, K)
    cursors = copy(stream_starts)

    # Initialize heap with first element of each stream
    heap_size = 0
    for k in 1:K
        if cursors[k] <= stream_ends[k]
            heap_size += 1
            heap_mz[heap_size] = @inbounds input[cursors[k]].mz
            heap_k[heap_size] = k
        end
    end

    # Build min-heap
    for i in div(heap_size, 2):-1:1
        _siftdown!(heap_mz, heap_k, i, heap_size)
    end

    out_idx = 0
    while heap_size > 0
        # Extract min
        out_idx += 1
        best_k = @inbounds heap_k[1]
        @inbounds output[out_idx] = input[cursors[best_k]]
        cursors[best_k] += 1

        if cursors[best_k] <= stream_ends[best_k]
            @inbounds heap_mz[1] = input[cursors[best_k]].mz
            # heap_k[1] unchanged (same stream)
        else
            # Stream exhausted — replace root with last element
            @inbounds heap_mz[1] = heap_mz[heap_size]
            @inbounds heap_k[1] = heap_k[heap_size]
            heap_size -= 1
        end
        heap_size > 0 && _siftdown!(heap_mz, heap_k, 1, heap_size)
    end

    return out_idx
end

@inline function _siftdown!(mz::Vector{Float32}, ks::Vector{Int}, i::Int, n::Int)
    while true
        smallest = i
        l = 2i
        r = 2i + 1
        if l <= n && @inbounds(mz[l]) < @inbounds(mz[smallest])
            smallest = l
        end
        if r <= n && @inbounds(mz[r]) < @inbounds(mz[smallest])
            smallest = r
        end
        smallest == i && return
        @inbounds mz[i], mz[smallest] = mz[smallest], mz[i]
        @inbounds ks[i], ks[smallest] = ks[smallest], ks[i]
        i = smallest
    end
end

function run_manual_merge_benchmarks()
    println("=" ^ 70)
    println("Part 3b: Manual heap merge vs DataStructures heap vs sort")
    println("=" ^ 70)

    configs = [
        (n_prec=200, frags=10, iso=2, label="Medium (200×10×2=4000)"),
        (n_prec=500, frags=10, iso=2, label="XLarge (500×10×2=10000)"),
        (n_prec=1000, frags=10, iso=2, label="Huge (1000×10×2=20000)"),
    ]

    for cfg in configs
        println("\n  Config: $(cfg.label)")
        data = generate_test_transitions(
            n_precursors=cfg.n_prec,
            frags_per_prec=cfg.frags,
            n_isotopes=cfg.iso,
        )
        N = data.n
        sr = data.stream_ranges
        K = length(sr)

        buffer = Vector{DetailedFrag{Float32}}(undef, N + 100)

        # Pre-compute starts/ends for manual merge
        starts = [first(r) for r in sr]
        ends = [last(r) for r in sr]

        # Manual heap merge
        manual_bench = @benchmark begin
            merge_sorted_streams_manual!($buffer, $(data.streamed), $starts, $ends, $K)
        end evals=1 samples=200

        # DataStructures heap merge
        ds_bench = @benchmark begin
            kway_merge!($buffer, $(data.streamed), $sr)
        end evals=1 samples=200

        # PartialQuickSort
        work = copy(data.interleaved)
        pqs_bench = @benchmark begin
            sort!(@view($work[1:$N]), by=x -> x.mz, alg=PartialQuickSort(1:$N))
        end setup=begin
            copyto!($work, 1, $(data.interleaved), 1, $N)
        end evals=1 samples=200

        t_manual = median(manual_bench).time / 1e3
        t_ds     = median(ds_bench).time / 1e3
        t_pqs    = median(pqs_bench).time / 1e3

        @printf("    N=%d, K=%d\n", N, K)
        @printf("    PartialQuickSort:    %8.1f μs\n", t_pqs)
        @printf("    DS heap merge:       %8.1f μs\n", t_ds)
        @printf("    Manual heap merge:   %8.1f μs\n", t_manual)
        @printf("    Speedup manual/PQS:  %5.2fx\n", t_pqs / t_manual)

        # Correctness check
        copyto!(work, 1, data.interleaved, 1, N)
        sort!(@view(work[1:N]), by=x -> x.mz, alg=PartialQuickSort(1:N))
        merge_sorted_streams_manual!(buffer, data.streamed, starts, ends, K)
        ok = all(work[i].mz == buffer[i].mz for i in 1:N)
        @printf("    Correct:             %s\n", ok ? "YES ✓" : "NO ✗")
    end
    println()
end

# ============================================================================
# Part 4: End-to-end feasibility estimate
# ============================================================================

function estimate_savings()
    println("=" ^ 70)
    println("Part 4: End-to-end feasibility estimate")
    println("=" ^ 70)

    # Representative call: 200 precs × 10 frags × 2 isotopes = 4000 transitions
    data = generate_test_transitions(n_precursors=200, frags_per_prec=10, n_isotopes=2)
    N = data.n
    sr = data.stream_ranges
    K = length(sr)

    buffer = Vector{DetailedFrag{Float32}}(undef, N + 100)
    starts = [first(r) for r in sr]
    ends = [last(r) for r in sr]

    # Warm up
    work = copy(data.interleaved)
    sort!(@view(work[1:N]), by=x -> x.mz, alg=PartialQuickSort(1:N))
    merge_sorted_streams_manual!(buffer, data.streamed, starts, ends, K)

    # Measure
    sort_bench = @benchmark begin
        sort!(@view($work[1:$N]), by=x -> x.mz, alg=PartialQuickSort(1:$N))
    end setup=begin
        copyto!($work, 1, $(data.interleaved), 1, $N)
    end evals=1 samples=500

    merge_bench = @benchmark begin
        merge_sorted_streams_manual!($buffer, $(data.streamed), $starts, $ends, $K)
    end evals=1 samples=500

    t_sort  = median(sort_bench).time / 1e3  # μs
    t_merge = median(merge_bench).time / 1e3

    calls_per_file = 40_000
    n_files = 3
    total_calls = calls_per_file * n_files

    total_sort_s  = t_sort  * total_calls / 1e6  # seconds
    total_merge_s = t_merge * total_calls / 1e6
    savings_s = total_sort_s - total_merge_s

    println()
    @printf("  Per-call:  sort=%.1f μs, merge=%.1f μs (%.2fx speedup)\n",
            t_sort, t_merge, t_sort / t_merge)
    @printf("  Per-file (%dk calls):  sort=%.2f s, merge=%.2f s, savings=%.2f s\n",
            calls_per_file ÷ 1000, total_sort_s / n_files, total_merge_s / n_files, savings_s / n_files)
    @printf("  3-file run (%dk calls): sort=%.2f s, merge=%.2f s, savings=%.2f s\n",
            total_calls ÷ 1000, total_sort_s, total_merge_s, savings_s)
    println()

    # Note overhead of tracking stream boundaries
    println("  Overhead considerations:")
    println("    - Must track stream start/end during fill (one UnitRange per precursor×isotope)")
    println("    - Extra buffer of size N needed for merge output")
    println("    - Fill order change: loop over isotope levels first, then precursors")
    println("    - K = n_precursors × n_isotopes (typically $K for this config)")
    println()
end

# ============================================================================
# Part 5: Real library benchmark (if available)
# ============================================================================

function benchmark_with_real_library(lib_dir::String)
    println("=" ^ 70)
    println("Part 5: Real library fragment m/z analysis")
    println("=" ^ 70)

    lookup = load_fragment_lookup(lib_dir)
    if lookup === nothing
        println("  Skipped: library files not found in $lib_dir")
        println()
        return
    end

    frags = getFragments(lookup)
    n_precs = length(lookup.prec_frag_ranges) - 1

    # Gather real fragment counts per precursor
    frag_counts = Int[]
    for prec_i in 1:n_precs
        r = getPrecFragRange(lookup, UInt32(prec_i))
        !isempty(r) && push!(frag_counts, length(r))
    end

    println("  Library: $lib_dir")
    println("  Precursors (non-empty): $(length(frag_counts))")
    if !isempty(frag_counts)
        @printf("  Fragments per precursor: median=%d, mean=%.1f, max=%d\n",
                median(frag_counts), mean(frag_counts), maximum(frag_counts))
    end

    # Simulate a realistic selectTransitions! call:
    # Pick ~200 random precursors, generate their transitions
    n_sample = min(200, length(frag_counts))
    sample_precs = sort(randperm(n_precs)[1:n_sample])

    # Count total fragments
    total_frags = sum(length(getPrecFragRange(lookup, UInt32(p))) for p in sample_precs)
    n_iso = 2
    N_expected = total_frags * n_iso

    println("  Sample: $n_sample precursors, $total_frags fragments, ~$(N_expected) transitions (with $n_iso isotopes)")

    # Build synthetic transitions in both fill orders
    transitions_interleaved = Vector{DetailedFrag{Float32}}(undef, N_expected + 1000)
    transitions_streamed = Vector{DetailedFrag{Float32}}(undef, N_expected + 1000)
    stream_ranges = UnitRange{Int}[]

    # Interleaved fill (current)
    idx = 0
    for p in sample_precs
        r = getPrecFragRange(lookup, UInt32(p))
        for fi in r
            frag = frags[fi]
            for iso in 0:(n_iso-1)
                idx += 1
                mz = Float32(frag.mz + iso * NEUTRON_F32 / max(frag.frag_charge, UInt8(1)))
                transitions_interleaved[idx] = DetailedFrag{Float32}(
                    frag.prec_id, mz, frag.intensity,
                    frag.ion_type, frag.is_y, frag.is_b, frag.is_p, iso > 0,
                    frag.frag_charge, frag.ion_position, frag.prec_charge,
                    frag.rank, frag.sulfur_count
                )
            end
        end
    end
    N = idx

    # Streamed fill (proposed: all M+0, then all M+1, ...)
    idx = 0
    for iso in 0:(n_iso-1)
        for p in sample_precs
            r = getPrecFragRange(lookup, UInt32(p))
            stream_start = idx + 1
            for fi in r
                frag = frags[fi]
                idx += 1
                mz = Float32(frag.mz + iso * NEUTRON_F32 / max(frag.frag_charge, UInt8(1)))
                transitions_streamed[idx] = DetailedFrag{Float32}(
                    frag.prec_id, mz, frag.intensity,
                    frag.ion_type, frag.is_y, frag.is_b, frag.is_p, iso > 0,
                    frag.frag_charge, frag.ion_position, frag.prec_charge,
                    frag.rank, frag.sulfur_count
                )
            end
            push!(stream_ranges, stream_start:idx)
        end
    end
    @assert idx == N

    K = length(stream_ranges)
    buffer = Vector{DetailedFrag{Float32}}(undef, N + 100)
    starts = [first(r) for r in stream_ranges]
    ends   = [last(r) for r in stream_ranges]

    # Verify each stream is sorted
    n_unsorted_streams = 0
    for r in stream_ranges
        for j in (first(r)+1):last(r)
            if transitions_streamed[j].mz < transitions_streamed[j-1].mz
                n_unsorted_streams += 1
                break
            end
        end
    end
    println("  Unsorted streams: $n_unsorted_streams / $K")

    # Benchmark
    work = copy(transitions_interleaved)

    pqs_bench = @benchmark begin
        sort!(@view($work[1:$N]), by=x -> x.mz, alg=PartialQuickSort(1:$N))
    end setup=begin
        copyto!($work, 1, $transitions_interleaved, 1, $N)
    end evals=1 samples=200

    manual_bench = @benchmark begin
        merge_sorted_streams_manual!($buffer, $transitions_streamed, $starts, $ends, $K)
    end evals=1 samples=200

    t_pqs    = median(pqs_bench).time / 1e3
    t_manual = median(manual_bench).time / 1e3

    # Correctness
    copyto!(work, 1, transitions_interleaved, 1, N)
    sort!(@view(work[1:N]), by=x -> x.mz, alg=PartialQuickSort(1:N))
    merge_sorted_streams_manual!(buffer, transitions_streamed, starts, ends, K)
    ok = all(work[i].mz == buffer[i].mz for i in 1:N)

    @printf("\n  N=%d, K=%d\n", N, K)
    @printf("  PartialQuickSort:  %8.1f μs\n", t_pqs)
    @printf("  Manual heap merge: %8.1f μs\n", t_manual)
    @printf("  Speedup:           %5.2fx\n", t_pqs / t_manual)
    @printf("  Correct:           %s\n", ok ? "YES ✓" : "NO ✗")
    println()
end

# ============================================================================
# Part 6: Alternative Sort Strategies for Transition Sorting
# ============================================================================
#
# Context: selectTransitions! (selectTransitions.jl:90-92) sorts ~4000
# DetailedFrag{Float32} transitions by m/z after filling. Called ~40k times
# per file. Current code uses PartialQuickSort(1:N) on a view.
#
# matchPeaks! (matchPeaks.jl:350) does a single-pass merge scan that requires
# transitions to be strictly m/z-sorted.
#
# DetailedFrag{Float32} is 24 bytes (21 used + 3 padding). A UInt16 sort key
# fits in existing padding (struct stays 24 bytes).
#
# Strategies tested:
#   Tier 1 (no structural changes):
#     1. QuickSort on view + by=getMZ
#     2. PartialQuickSort(1:N) on FULL array (no view)
#     3. QuickSort on view + lt comparator
#   Tier 2 (UInt16 radix key in padding):
#     4. UInt16 radix + insertion sort cleanup
#     5. UInt16 keys in parallel array + radix + permute + insertion cleanup
#   Tier 3 (packed UInt64 radix, exact precision):
#     6. Packed UInt64 = (round(UInt32, mz*10000) << 32) | index
# ============================================================================

"""
Generate a simple vector of N random DetailedFrag{Float32} transitions with
m/z values in [100, 2000]. Uses a full-size buffer of `buf_size` elements
(default 100_000) to simulate the real pre-allocated transitions array.
Returns (buffer, N).
"""
function generate_sort_test_data(N::Int; buf_size::Int=100_000, seed::Int=42)
    rng = MersenneTwister(seed)
    buf = Vector{DetailedFrag{Float32}}(undef, buf_size)
    for i in 1:N
        mz = Float32(rand(rng) * 1900.0 + 100.0)
        buf[i] = DetailedFrag{Float32}(
            UInt32(i), mz, Float16(1.0),
            UInt16(0), true, false, false, false,
            UInt8(2), UInt8(1), UInt8(2), UInt8(1), UInt8(0)
        )
    end
    return buf, N
end

# --- Helper: insertion sort on a view by .mz (O(N) on almost-sorted data) ---
function insertion_sort_mz!(v::AbstractVector{DetailedFrag{Float32}}, lo::Int, hi::Int)
    @inbounds for i in (lo+1):hi
        key = v[i]
        key_mz = key.mz
        j = i - 1
        while j >= lo && v[j].mz > key_mz
            v[j+1] = v[j]
            j -= 1
        end
        v[j+1] = key
    end
end

# --- UInt16 sort key computation (linear quantization over [0, 2000] Da) ---
# Resolution = 2000 / 65535 ≈ 0.031 Da
@inline function mz_to_uint16_key(mz::Float32)::UInt16
    return round(UInt16, clamp(mz * Float32(65535.0 / 2000.0), 0.0f0, 65535.0f0))
end

# --- Zero-allocation counting sort on UInt16 keys ---
# Single-pass counting sort: O(N + 65536). Requires pre-allocated:
#   counts::Vector{Int} of length 65536 (reused across calls)
#   scratch::Vector{DetailedFrag{Float32}} of length ≥ N (reused across calls)
# After counting sort, data is sorted by UInt16 key. Insertion sort cleanup
# fixes within-bucket misordering (O(N) since data is almost-sorted).
function countsort_uint16_mz!(
    data::Vector{DetailedFrag{Float32}}, N::Int,
    scratch::Vector{DetailedFrag{Float32}},
    counts::Vector{Int},
)
    NBUCKETS = 65536

    # Step 1: Zero counts
    @inbounds for i in 1:NBUCKETS
        counts[i] = 0
    end

    # Step 2: Count occurrences
    @inbounds for i in 1:N
        k = Int(mz_to_uint16_key(data[i].mz)) + 1  # 1-indexed
        counts[k] += 1
    end

    # Step 3: Prefix sum (exclusive → gives write positions)
    running = 0
    @inbounds for i in 1:NBUCKETS
        c = counts[i]
        counts[i] = running
        running += c
    end

    # Step 4: Scatter into scratch in sorted order
    @inbounds for i in 1:N
        k = Int(mz_to_uint16_key(data[i].mz)) + 1
        counts[k] += 1
        scratch[counts[k]] = data[i]
    end

    # Step 5: Copy back
    @inbounds copyto!(data, 1, scratch, 1, N)

    # Step 6: Insertion sort cleanup on actual Float32 mz
    insertion_sort_mz!(data, 1, N)
end

# --- Zero-allocation 2-pass LSD radix sort on Float32 mz ---
# Reinterpret Float32 as UInt32 with sign-flip for correct ordering,
# then do 2-pass radix sort (low 16 bits, high 16 bits). Each pass
# uses a 65536-element count array and scatters into scratch buffer.
# No insertion cleanup needed — this is exact.
@inline function float32_to_sortable_uint32(f::Float32)::UInt32
    u = reinterpret(UInt32, f)
    # IEEE 754: if sign bit set, flip all bits; else flip sign bit only
    # This makes the UInt32 ordering match Float32 ordering
    mask = ifelse(u & 0x80000000 != 0, 0xFFFFFFFF, 0x80000000)
    return xor(u, mask)
end

function radixsort_float32_mz!(
    data::Vector{DetailedFrag{Float32}}, N::Int,
    scratch::Vector{DetailedFrag{Float32}},
    counts::Vector{Int},  # length ≥ 65536
)
    NBUCKETS = 65536

    # --- Pass 1: sort by low 16 bits of sortable key ---
    # Zero counts
    @inbounds for i in 1:NBUCKETS; counts[i] = 0; end

    # Count
    @inbounds for i in 1:N
        k = Int(float32_to_sortable_uint32(data[i].mz) & 0xFFFF) + 1
        counts[k] += 1
    end

    # Prefix sum
    running = 0
    @inbounds for i in 1:NBUCKETS
        c = counts[i]
        counts[i] = running
        running += c
    end

    # Scatter data → scratch
    @inbounds for i in 1:N
        k = Int(float32_to_sortable_uint32(data[i].mz) & 0xFFFF) + 1
        counts[k] += 1
        scratch[counts[k]] = data[i]
    end

    # --- Pass 2: sort by high 16 bits of sortable key ---
    # Zero counts
    @inbounds for i in 1:NBUCKETS; counts[i] = 0; end

    # Count
    @inbounds for i in 1:N
        k = Int(float32_to_sortable_uint32(scratch[i].mz) >> 16) + 1
        counts[k] += 1
    end

    # Prefix sum
    running = 0
    @inbounds for i in 1:NBUCKETS
        c = counts[i]
        counts[i] = running
        running += c
    end

    # Scatter scratch → data
    @inbounds for i in 1:N
        k = Int(float32_to_sortable_uint32(scratch[i].mz) >> 16) + 1
        counts[k] += 1
        data[counts[k]] = scratch[i]
    end
end

# --- Zero-allocation 4-pass LSD radix sort on Float32 mz (byte-wise) ---
# 4 passes, one per byte. Each pass uses a 256-element count array.
# Much smaller count array (256 vs 65536) — better cache behavior.
function radixsort_float32_mz_4pass!(
    data::Vector{DetailedFrag{Float32}}, N::Int,
    scratch::Vector{DetailedFrag{Float32}},
    counts::Vector{Int},  # length ≥ 256
)
    NBUCKETS = 256

    # Pass 1-4: sort by byte 0 (LSB) through byte 3 (MSB)
    src = data
    dst = scratch
    for pass in 0:3
        shift = pass * 8

        # Zero counts
        @inbounds for i in 1:NBUCKETS; counts[i] = 0; end

        # Count
        @inbounds for i in 1:N
            k = Int((float32_to_sortable_uint32(src[i].mz) >> shift) & 0xFF) + 1
            counts[k] += 1
        end

        # Prefix sum
        running = 0
        @inbounds for i in 1:NBUCKETS
            c = counts[i]
            counts[i] = running
            running += c
        end

        # Scatter src → dst
        @inbounds for i in 1:N
            k = Int((float32_to_sortable_uint32(src[i].mz) >> shift) & 0xFF) + 1
            counts[k] += 1
            dst[counts[k]] = src[i]
        end

        # Swap src/dst for next pass
        src, dst = dst, src
    end

    # After 4 passes (even number), result is back in `data` (src started as data,
    # swapped 4 times → src=data again). If odd passes, we'd need a final copy.
    # 4 passes: data→scratch→data→scratch→data. Actually let's verify:
    # pass 0: src=data, dst=scratch → scatter into scratch. swap: src=scratch, dst=data
    # pass 1: src=scratch, dst=data → scatter into data. swap: src=data, dst=scratch
    # pass 2: src=data, dst=scratch → scatter into scratch. swap: src=scratch, dst=data
    # pass 3: src=scratch, dst=data → scatter into data. swap: src=data, dst=scratch
    # After loop: src=data, dst=scratch. Result is in data. ✓
end

"""
    run_sort_strategy_benchmarks()

Part 6: Benchmark alternative sort strategies for transition sorting.
Compares approaches across multiple array sizes, including zero-allocation
custom radix sorts.
"""
function run_sort_strategy_benchmarks()
    println("=" ^ 70)
    println("Part 6: Alternative Sort Strategies for Transition Sorting")
    println("=" ^ 70)

    sizes = [1000, 4000, 10000, 20000]

    # --- Strategy functions ---

    # Baseline: current code (PartialQuickSort on view + by=getMZ)
    function sort_baseline!(buf, N)
        sort!(@view(buf[1:N]), by=x -> getMZ(x), alg=PartialQuickSort(1:N))
    end

    # Strategy 1: QuickSort on view + by=getMZ
    function sort_quicksort_by!(buf, N)
        sort!(@view(buf[1:N]), by=x -> getMZ(x), alg=QuickSort)
    end

    # Strategy 2: PartialQuickSort on FULL array (no view)
    # NOTE: Expected to fail correctness — PQS(1:N) on the full 100k buffer
    # selects the N smallest from ALL elements (including uninitialized slots),
    # not just the first N. Included to demonstrate this pitfall.
    function sort_pqs_full!(buf, N)
        sort!(buf, by=x -> getMZ(x), alg=PartialQuickSort(1:N))
    end

    # Strategy 3: QuickSort on view + lt comparator
    function sort_quicksort_lt!(buf, N)
        sort!(@view(buf[1:N]), lt=(a,b) -> a.mz < b.mz, alg=QuickSort)
    end

    # Strategy 4: UInt16 radix + insertion sort cleanup
    # Compute UInt16 keys into a parallel array, sortperm by key, permute, then
    # insertion sort on actual mz.
    function sort_uint16_radix_insertion!(buf, N, keys_buf, perm_buf, temp_buf)
        # Step 1: Extract UInt16 keys
        @inbounds for i in 1:N
            keys_buf[i] = mz_to_uint16_key(buf[i].mz)
        end
        # Step 2: Radix sort keys to get permutation
        # Julia sorts integers with radix sort automatically
        @inbounds for i in 1:N
            perm_buf[i] = i
        end
        sort!(@view(perm_buf[1:N]), by=i -> @inbounds(keys_buf[i]))
        # Step 3: Apply permutation
        @inbounds for i in 1:N
            temp_buf[i] = buf[perm_buf[i]]
        end
        @inbounds copyto!(buf, 1, temp_buf, 1, N)
        # Step 4: Insertion sort cleanup on actual mz (data is almost-sorted)
        insertion_sort_mz!(buf, 1, N)
    end

    # Strategy 5: UInt16 keys in parallel array + radix + permute + insertion cleanup
    # Same idea as 4, but use a direct key array sorted with sortperm! for less overhead
    function sort_uint16_parallel!(buf, N, keys_buf, perm_buf, temp_buf)
        # Step 1: Extract UInt16 keys
        @inbounds for i in 1:N
            keys_buf[i] = mz_to_uint16_key(buf[i].mz)
        end
        # Step 2: Sort keys directly and track permutation via packed UInt64
        # Pack (key << 16) | index won't work for UInt16 indices > 65535
        # Instead, pack (UInt32(key) << 16) | UInt32(index) — works for N ≤ 65535
        # For N > 65535, use UInt64 packing
        @inbounds for i in 1:N
            perm_buf[i] = (UInt64(keys_buf[i]) << 32) | UInt64(i)
        end
        sort!(@view(perm_buf[1:N]))  # Radix sort on UInt64
        # Step 3: Extract permutation and apply
        @inbounds for i in 1:N
            idx = Int(perm_buf[i] & 0xFFFFFFFF)
            temp_buf[i] = buf[idx]
        end
        @inbounds copyto!(buf, 1, temp_buf, 1, N)
        # Step 4: Insertion sort cleanup
        insertion_sort_mz!(buf, 1, N)
    end

    # Strategy 6: Packed UInt64 = (round(UInt32, mz*10000) << 32) | index
    # Near-exact precision (0.0001 Da). Julia radix-sorts UInt64 automatically.
    # Insertion cleanup handles Float32 ties within the same UInt32 bucket.
    function sort_packed_uint64!(buf, N, pack_buf, temp_buf)
        # Step 1: Pack (mz_key << 32) | index
        @inbounds for i in 1:N
            mz_key = round(UInt32, clamp(Float64(buf[i].mz) * 10000.0, 0.0, 4.294967295e9))
            pack_buf[i] = (UInt64(mz_key) << 32) | UInt64(i)
        end
        # Step 2: Radix sort
        sort!(@view(pack_buf[1:N]))
        # Step 3: Extract indices and permute
        @inbounds for i in 1:N
            idx = Int(pack_buf[i] & 0xFFFFFFFF)
            temp_buf[i] = buf[idx]
        end
        @inbounds copyto!(buf, 1, temp_buf, 1, N)
        # Step 4: Insertion sort cleanup for Float32 ties within same UInt32 bucket
        insertion_sort_mz!(buf, 1, N)
    end

    # Pre-allocate working buffers (worst case: 20000)
    max_n = maximum(sizes)
    buf_size = 100_000
    work      = Vector{DetailedFrag{Float32}}(undef, buf_size)
    keys_buf  = Vector{UInt16}(undef, max_n)
    perm_buf_int = Vector{Int}(undef, max_n)
    perm_buf_u64 = Vector{UInt64}(undef, max_n)
    temp_buf  = Vector{DetailedFrag{Float32}}(undef, max_n)
    pack_buf  = Vector{UInt64}(undef, max_n)
    radix_counts = Vector{Int}(undef, 65536)  # reused for counting/radix sorts
    radix_scratch = Vector{DetailedFrag{Float32}}(undef, max_n)

    for N in sizes
        println("\n  N = $N (buffer size = $buf_size)")
        println("  " * "-"^64)

        orig, _ = generate_sort_test_data(N; buf_size=buf_size)

        # --- Compute baseline sorted result for correctness verification ---
        copyto!(work, 1, orig, 1, N)
        sort_baseline!(work, N)
        baseline_mzs = [work[i].mz for i in 1:N]

        # Verify baseline is actually sorted
        for i in 2:N
            @assert baseline_mzs[i] >= baseline_mzs[i-1] "Baseline not sorted at index $i"
        end

        strategies = [
            ("Baseline (PQS+view+by)",
             () -> sort_baseline!(work, N),
             () -> copyto!(work, 1, orig, 1, N)),
            ("1. QS+view+by=getMZ",
             () -> sort_quicksort_by!(work, N),
             () -> copyto!(work, 1, orig, 1, N)),
            ("2. PQS full array",
             () -> sort_pqs_full!(work, N),
             () -> copyto!(work, 1, orig, 1, N)),
            ("3. QS+view+lt",
             () -> sort_quicksort_lt!(work, N),
             () -> copyto!(work, 1, orig, 1, N)),
            ("4. UInt16 radix+ins",
             () -> sort_uint16_radix_insertion!(work, N, keys_buf, perm_buf_int, temp_buf),
             () -> copyto!(work, 1, orig, 1, N)),
            ("5. UInt16 packed+ins",
             () -> sort_uint16_parallel!(work, N, keys_buf, perm_buf_u64, temp_buf),
             () -> copyto!(work, 1, orig, 1, N)),
            ("6. UInt64 packed radix",
             () -> sort_packed_uint64!(work, N, pack_buf, temp_buf),
             () -> copyto!(work, 1, orig, 1, N)),
            ("7. CountSort U16+ins",
             () -> countsort_uint16_mz!(work, N, radix_scratch, radix_counts),
             () -> copyto!(work, 1, orig, 1, N)),
            ("8. Radix F32 2-pass",
             () -> radixsort_float32_mz!(work, N, radix_scratch, radix_counts),
             () -> copyto!(work, 1, orig, 1, N)),
            ("9. Radix F32 4-pass",
             () -> radixsort_float32_mz_4pass!(work, N, radix_scratch, radix_counts),
             () -> copyto!(work, 1, orig, 1, N)),
        ]

        results = []
        for (name, sort_fn, setup_fn) in strategies
            bench = @benchmark $sort_fn() setup=($setup_fn()) evals=1 samples=200

            # Correctness check
            setup_fn()
            sort_fn()
            correct = true
            for i in 1:N
                if work[i].mz != baseline_mzs[i]
                    correct = false
                    break
                end
            end

            t_us = median(bench).time / 1e3
            allocs = median(bench).allocs
            push!(results, (name=name, time_us=t_us, allocs=allocs, correct=correct))
        end

        # Print results table
        baseline_time = results[1].time_us
        @printf("    %-25s %10s %8s %8s %8s\n", "Strategy", "Median μs", "Allocs", "Speedup", "Correct")
        @printf("    %-25s %10s %8s %8s %8s\n", "-"^25, "-"^10, "-"^8, "-"^8, "-"^8)
        for r in results
            speedup = baseline_time / r.time_us
            @printf("    %-25s %10.1f %8d %7.2fx %8s\n",
                    r.name, r.time_us, r.allocs, speedup,
                    r.correct ? "YES" : "NO")
        end
    end

    # --- End-to-end savings estimate ---
    println("\n  " * "="^64)
    println("  End-to-end savings estimate (40k calls × 3 files = 120k calls)")
    println("  " * "="^64)

    # Use N=4000 as representative
    N_rep = 4000
    orig_rep, _ = generate_sort_test_data(N_rep; buf_size=buf_size)

    # Measure baseline and best strategies with more samples
    # Inline the sort calls to avoid closure scoping issues with @benchmark
    baseline_bench = @benchmark begin
        sort!(@view($work[1:$N_rep]), by=x -> getMZ(x), alg=PartialQuickSort(1:$N_rep))
    end setup=(copyto!($work, 1, $orig_rep, 1, $N_rep)) evals=1 samples=500

    qs_lt_bench = @benchmark begin
        sort!(@view($work[1:$N_rep]), lt=(a,b) -> a.mz < b.mz, alg=QuickSort)
    end setup=(copyto!($work, 1, $orig_rep, 1, $N_rep)) evals=1 samples=500

    countsort_bench = @benchmark begin
        countsort_uint16_mz!($work, $N_rep, $radix_scratch, $radix_counts)
    end setup=(copyto!($work, 1, $orig_rep, 1, $N_rep)) evals=1 samples=500

    radix2_bench = @benchmark begin
        radixsort_float32_mz!($work, $N_rep, $radix_scratch, $radix_counts)
    end setup=(copyto!($work, 1, $orig_rep, 1, $N_rep)) evals=1 samples=500

    radix4_bench = @benchmark begin
        radixsort_float32_mz_4pass!($work, $N_rep, $radix_scratch, $radix_counts)
    end setup=(copyto!($work, 1, $orig_rep, 1, $N_rep)) evals=1 samples=500

    t_base   = median(baseline_bench).time / 1e3
    t_lt     = median(qs_lt_bench).time / 1e3
    t_csort  = median(countsort_bench).time / 1e3
    t_radix2 = median(radix2_bench).time / 1e3
    t_radix4 = median(radix4_bench).time / 1e3

    total_calls = 40_000 * 3
    for (name, t) in [("Baseline (PQS+view+by)", t_base),
                       ("3. QS+view+lt", t_lt),
                       ("7. CountSort U16+ins", t_csort),
                       ("8. Radix F32 2-pass", t_radix2),
                       ("9. Radix F32 4-pass", t_radix4)]
        total_s = t * total_calls / 1e6
        savings_s = (t_base - t) * total_calls / 1e6
        @printf("    %-25s  %6.1f μs/call  → %5.2f s total  (savings: %+.2f s)\n",
                name, t, total_s, savings_s)
    end

    println()
end


# ============================================================================
# Part 7: Sort-Free Fragment Matching via Binary Search
# ============================================================================
#
# Context: selectTransitions! (selectTransitions.jl:90-92) sorts ~4000
# DetailedFrag{Float32} transitions by m/z before every matchPeaks! call
# (~40k times per file). matchPeaks! (matchPeaks.jl:350-425) then does a
# single-pass O(T+P) merge scan.
#
# DIA-NN avoids sorting entirely by binary-searching each fragment into
# pre-sorted peaks. This benchmark tests whether eliminating the sort via
# binary search is faster end-to-end.
#
# Approaches:
#   A. QuickSort + merge-scan (current pipeline baseline)
#   B. Radix sort (4-pass) + merge-scan
#   C. Binary search per transition (unsorted, no sort needed)
#   D. Binary search + DIA-NN-style bound narrowing (sorted per precursor)
# ============================================================================

"""
    generate_matching_data(n_peaks, n_precursors, frags_per_prec; kwargs...)

Generate synthetic data for matching benchmarks.

Returns `(peaks, trans_unsorted, trans_sorted_prec, prec_bounds, N)`:
- `peaks`: sorted `Vector{Float32}` of empirical peak m/z in [100, 2000] Da
- `trans_unsorted`: `DetailedFrag{Float32}` interleaved across precursors (current behavior)
- `trans_sorted_prec`: same fragments, contiguous per precursor, m/z-sorted within each group
- `prec_bounds`: `Vector{Tuple{Int,Int}}` with `(start, stop)` for each precursor group
- `N`: total transition count

~`match_fraction` of transitions are placed within 80% of `ppm_tol` of a random peak;
the rest are uniformly random and almost certainly outside tolerance.
"""
function generate_matching_data(
    n_peaks::Int, n_precursors::Int, frags_per_prec::Int;
    match_fraction::Float64=0.3, ppm_tol::Float64=20.0, seed::Int=42
)
    rng = MersenneTwister(seed)
    N = n_precursors * frags_per_prec

    # Sorted peaks in [100, 2000] Da
    peaks = sort!(Float32[rand(rng) * 1900.0f0 + 100.0f0 for _ in 1:n_peaks])

    # Generate (prec_id, mz) pairs
    frag_data = Vector{Tuple{UInt32,Float32}}(undef, N)
    idx = 0
    for p in 1:n_precursors
        pid = UInt32(p)
        n_match = round(Int, frags_per_prec * match_fraction)
        n_rand = frags_per_prec - n_match
        for _ in 1:n_match
            peak_mz = peaks[rand(rng, 1:n_peaks)]
            max_off = Float32(ppm_tol * 1e-6 * peak_mz * 0.8)
            mz = peak_mz + Float32((rand(rng) * 2.0 - 1.0) * max_off)
            idx += 1
            frag_data[idx] = (pid, mz)
        end
        for _ in 1:n_rand
            mz = Float32(rand(rng) * 1900.0 + 100.0)
            idx += 1
            frag_data[idx] = (pid, mz)
        end
    end

    # Unsorted: shuffle across precursors (simulates current interleaved fill)
    shuffle!(rng, frag_data)
    trans_unsorted = Vector{DetailedFrag{Float32}}(undef, N)
    @inbounds for i in 1:N
        pid, mz = frag_data[i]
        trans_unsorted[i] = DetailedFrag{Float32}(
            pid, mz, Float16(1.0), UInt16(0),
            true, false, false, false,
            UInt8(2), UInt8(1), UInt8(2), UInt8(1), UInt8(0))
    end

    # Sorted per precursor: contiguous groups, m/z-sorted within each
    by_prec = Dict{UInt32,Vector{Float32}}()
    for (pid, mz) in frag_data
        push!(get!(Vector{Float32}, by_prec, pid), mz)
    end

    trans_sorted_prec = Vector{DetailedFrag{Float32}}(undef, N)
    prec_bounds = Vector{Tuple{Int,Int}}()
    sizehint!(prec_bounds, n_precursors)
    tidx = 0
    for pid in sort!(collect(keys(by_prec)))
        fmzs = sort!(by_prec[pid])
        start = tidx + 1
        for mz in fmzs
            tidx += 1
            trans_sorted_prec[tidx] = DetailedFrag{Float32}(
                pid, mz, Float16(1.0), UInt16(0),
                true, false, false, false,
                UInt8(2), UInt8(1), UInt8(2), UInt8(1), UInt8(0))
        end
        push!(prec_bounds, (start, tidx))
    end

    return peaks, trans_unsorted, trans_sorted_prec, prec_bounds, N
end

"""
    mock_merge_match!(sorted_trans, N, peaks, ppm_tol) -> (n_matched, n_missed)

Two-pointer merge scan matching sorted transitions against sorted peaks.
Mirrors matchPeaks! (matchPeaks.jl:350-425): advance ion on match or overshoot,
advance peak when below tolerance. Includes setNearest!-style forward scan
to find best peak within tolerance window.
"""
function mock_merge_match!(
    sorted_trans::AbstractVector{DetailedFrag{Float32}},
    N::Int, peaks::Vector{Float32}, ppm_tol::Float64
)
    n_matched = 0
    n_missed = 0
    peak_idx = 1
    ion_idx = 1
    n_peaks = length(peaks)
    (N < 1 || n_peaks < 1) && return (0, N)

    # Pre-compute bounds for first ion (mirrors matchPeaks! line 348)
    ion_mz = sorted_trans[1].mz
    tol = Float32(ppm_tol * 1e-6 * ion_mz)
    low = ion_mz - tol
    high = ion_mz + tol

    @inbounds @fastmath while peak_idx <= n_peaks && ion_idx <= N
        peak_mz = peaks[peak_idx]
        if peak_mz >= low
            if peak_mz <= high
                # Match — scan forward for nearest peak (setNearest! analog)
                best_diff = abs(peak_mz - ion_mz)
                j = peak_idx + 1
                while j <= n_peaks && peaks[j] <= high
                    d = abs(peaks[j] - ion_mz)
                    best_diff = ifelse(d < best_diff, d, best_diff)
                    j += 1
                end
                n_matched += 1
                # Advance ion, keep peak_idx (multiple ions can match same peak)
                ion_idx += 1
                (ion_idx > N) && return (n_matched, n_missed)
                ion_mz = sorted_trans[ion_idx].mz
                tol = Float32(ppm_tol * 1e-6 * ion_mz)
                low = ion_mz - tol
                high = ion_mz + tol
                continue
            end
            # Peak past high bound — no match for this ion
            n_missed += 1
            ion_idx += 1
            (ion_idx > N) && return (n_matched, n_missed)
            ion_mz = sorted_trans[ion_idx].mz
            tol = Float32(ppm_tol * 1e-6 * ion_mz)
            low = ion_mz - tol
            high = ion_mz + tol
            continue
        end
        # Peak below low bound — advance peak
        peak_idx += 1
    end

    # Remaining ions past all peaks are misses
    n_missed += N - ion_idx + 1
    return (n_matched, n_missed)
end

"""
    mock_bsearch_match!(trans, N, peaks, ppm_tol) -> (n_matched, n_missed)

Per-transition binary search into sorted peaks. No sort required.
For each transition: binary search for first peak >= low bound,
then linear scan for best match within tolerance.
"""
function mock_bsearch_match!(
    trans::AbstractVector{DetailedFrag{Float32}},
    N::Int, peaks::Vector{Float32}, ppm_tol::Float64
)
    n_matched = 0
    n_missed = 0
    n_peaks = length(peaks)

    @inbounds @fastmath for i in 1:N
        ion_mz = trans[i].mz
        tol = Float32(ppm_tol * 1e-6 * ion_mz)
        low = ion_mz - tol
        high = ion_mz + tol

        # Binary search: first index where peaks[lo] >= low
        lo = 1
        hi = n_peaks + 1
        while lo < hi
            mid = (lo + hi) >>> 1
            if peaks[mid] < low
                lo = mid + 1
            else
                hi = mid
            end
        end

        # Linear scan for best match within [low, high]
        found = false
        best_diff = typemax(Float32)
        j = lo
        while j <= n_peaks && peaks[j] <= high
            d = abs(peaks[j] - ion_mz)
            best_diff = ifelse(d < best_diff, d, best_diff)
            found = true
            j += 1
        end

        n_matched += found
        n_missed += !found
    end

    return (n_matched, n_missed)
end

"""
    mock_bsearch_bounds_match!(trans, N, peaks, ppm_tol, prec_bounds) -> (n_matched, n_missed)

Binary search with DIA-NN-style bound narrowing. Transitions must be contiguous
per precursor and m/z-sorted within each group.

First fragment per precursor group: full binary search into peaks.
Subsequent fragments: narrow forward from previous position (since m/z increases
within a precursor group, the search start can only move forward).

Amortized cost per precursor group: O(F + log M) instead of O(F log M).
"""
function mock_bsearch_bounds_match!(
    trans::AbstractVector{DetailedFrag{Float32}},
    N::Int, peaks::Vector{Float32}, ppm_tol::Float64,
    prec_bounds::Vector{Tuple{Int,Int}}
)
    n_matched = 0
    n_missed = 0
    n_peaks = length(peaks)

    @inbounds @fastmath for (start_idx, stop_idx) in prec_bounds
        search_lo = 1  # Reset per precursor group

        for i in start_idx:stop_idx
            ion_mz = trans[i].mz
            tol = Float32(ppm_tol * 1e-6 * ion_mz)
            low = ion_mz - tol
            high = ion_mz + tol

            if i == start_idx
                # Full binary search for first fragment in group
                lo = search_lo
                hi = n_peaks + 1
                while lo < hi
                    mid = (lo + hi) >>> 1
                    if peaks[mid] < low
                        lo = mid + 1
                    else
                        hi = mid
                    end
                end
                search_lo = lo
            else
                # Forward scan from previous position (m/z monotonically increases)
                while search_lo <= n_peaks && peaks[search_lo] < low
                    search_lo += 1
                end
            end

            # Linear scan for best match within [low, high]
            found = false
            best_diff = typemax(Float32)
            j = search_lo
            while j <= n_peaks && peaks[j] <= high
                d = abs(peaks[j] - ion_mz)
                best_diff = ifelse(d < best_diff, d, best_diff)
                found = true
                j += 1
            end

            n_matched += found
            n_missed += !found
        end
    end

    return (n_matched, n_missed)
end

"""
    bench_matching_cell(work, radix_scratch, radix_counts, peaks, trans_u, trans_s,
                        prec_bounds, N, ppm_tol; samples=200)

Run all four matching approaches on one data configuration.
Returns `(t_A, t_B, t_C, t_D, n_match, N)` with times in μs (minimum).
Uses `minimum` time — BenchmarkTools' recommended metric for microbenchmarks.

Must be a top-level function (not a closure) for BenchmarkTools' \$ interpolation.
"""
function bench_matching_cell(
    work::Vector{DetailedFrag{Float32}},
    radix_scratch::Vector{DetailedFrag{Float32}},
    radix_counts::Vector{Int},
    peaks::Vector{Float32},
    trans_u::Vector{DetailedFrag{Float32}},
    trans_s::Vector{DetailedFrag{Float32}},
    prec_bounds::Vector{Tuple{Int,Int}},
    N::Int,
    ppm_tol::Float64;
    samples::Int=200
)
    # Correctness verification
    copyto!(work, 1, trans_u, 1, N)
    sort!(@view(work[1:N]), lt=(a,b) -> a.mz < b.mz, alg=QuickSort)
    m_A, _ = mock_merge_match!(work, N, peaks, ppm_tol)
    m_C, _ = mock_bsearch_match!(trans_u, N, peaks, ppm_tol)
    m_D, _ = mock_bsearch_bounds_match!(trans_s, N, peaks, ppm_tol, prec_bounds)
    if m_A != m_C || m_A != m_D
        @printf("    !! MISMATCH: A=%d C=%d D=%d\n", m_A, m_C, m_D)
    end

    bA = @benchmark begin
        sort!(@view($work[1:$N]), lt=(a,b) -> a.mz < b.mz, alg=QuickSort)
        mock_merge_match!($work, $N, $peaks, $ppm_tol)
    end setup=(copyto!($work, 1, $trans_u, 1, $N)) evals=1 samples=samples

    bB = @benchmark begin
        radixsort_float32_mz_4pass!($work, $N, $radix_scratch, $radix_counts)
        mock_merge_match!($work, $N, $peaks, $ppm_tol)
    end setup=(copyto!($work, 1, $trans_u, 1, $N)) evals=1 samples=samples

    bC = @benchmark mock_bsearch_match!(
        $trans_u, $N, $peaks, $ppm_tol
    ) evals=1 samples=samples

    bD = @benchmark mock_bsearch_bounds_match!(
        $trans_s, $N, $peaks, $ppm_tol, $prec_bounds
    ) evals=1 samples=samples

    return (minimum(bA).time/1e3, minimum(bB).time/1e3,
            minimum(bC).time/1e3, minimum(bD).time/1e3,
            m_A, N)
end

function print_matching_header()
    @printf("    %-36s %9s %9s %9s %9s\n",
            "", "A:QS+mrg", "B:Rd+mrg", "C:BSrch", "D:BS+Nar")
    @printf("    %-36s %9s %9s %9s %9s\n",
            "", "(μs)", "(μs)", "(μs)", "(μs)")
    println("    " * "-"^72)
end

"""
    run_matching_strategy_benchmarks()

Part 7: Benchmark sort-free fragment matching via binary search.

Complexity analysis:
  A (QS + merge):      O(T·log(T) + M)     — comparison sort dominates
  B (radix + merge):   O(T + M)            — linear radix + linear merge
  C (binary search):   O(T·log(M))         — no sort, but log(M) per fragment
  D (BS + narrowing):  O(T + P·log(M))     — amortized: P = T/F precursors

Key question: For T=4000, log(T)≈12, log(M=1000)≈10. The sort cost T·log(T)
is similar to T·log(M). But D with large F gets O(T + P·log M) which can be
much less. Table 2 (F sweep) tests this directly.

Uses `minimum` time — BenchmarkTools' recommended metric for microbenchmarks.
"""
function run_matching_strategy_benchmarks()
    println("=" ^ 70)
    println("Part 7: Sort-Free Fragment Matching via Binary Search")
    println("=" ^ 70)
    println()
    println("  Approaches:")
    println("    A. QuickSort + merge-scan        O(T·log T + M)")
    println("    B. Radix sort (4-pass) + merge   O(T + M)  [linear radix + linear merge]")
    println("    C. Binary search (unsorted)       O(T·log M)")
    println("    D. BS + bound narrowing            O(T + P·log M)  [P = T/F precursors]")
    println()

    ppm_tol = 20.0

    # Pre-allocate sort working buffers
    max_T = 20000
    work = Vector{DetailedFrag{Float32}}(undef, max_T)
    radix_scratch = Vector{DetailedFrag{Float32}}(undef, max_T)
    radix_counts = Vector{Int}(undef, 65536)

    # --- JIT warmup (ensures compilation doesn't affect first benchmark) ---
    println("  Warming up JIT...")
    let
        pw, tu, ts, pb, nw = generate_matching_data(50, 10, 10)
        ww = Vector{DetailedFrag{Float32}}(undef, 100)
        sw = Vector{DetailedFrag{Float32}}(undef, 100)
        cw = Vector{Int}(undef, 65536)
        copyto!(ww, 1, tu, 1, nw)
        sort!(@view(ww[1:nw]), lt=(a,b) -> a.mz < b.mz, alg=QuickSort)
        mock_merge_match!(ww, nw, pw, ppm_tol)
        copyto!(ww, 1, tu, 1, nw)
        radixsort_float32_mz_4pass!(ww, nw, sw, cw)
        mock_merge_match!(ww, nw, pw, ppm_tol)
        mock_bsearch_match!(tu, nw, pw, ppm_tol)
        mock_bsearch_bounds_match!(ts, nw, pw, ppm_tol, pb)
        # Also warmup bench_matching_cell itself
        bench_matching_cell(ww, sw, cw, pw, tu, ts, pb, nw, ppm_tol; samples=3)
    end
    println()

    # ====================================================================
    # Table 1: T × M matrix, fixed F=20 frags/precursor
    # ====================================================================
    println("  Table 1: T (transitions) × M (peaks), fixed F=20 frags/precursor")
    println()
    print_matching_header()

    T_configs = [(1000, 50), (4000, 200), (10000, 500)]
    M_sizes = [100, 300, 1000]

    for (T, P) in T_configs
        for M in M_sizes
            peaks, tu, ts, pb, N = generate_matching_data(M, P, 20; ppm_tol=ppm_tol)
            tA, tB, tC, tD, nm, _ = bench_matching_cell(
                work, radix_scratch, radix_counts, peaks, tu, ts, pb, N, ppm_tol)
            mp = round(Int, 100 * nm / N)
            @printf("    T=%-5d M=%-5d (%2d%% match=%4d) %8.1f %9.1f %9.1f %9.1f\n",
                    T, M, mp, nm, tA, tB, tC, tD)
        end
        println()
    end

    # ====================================================================
    # Table 2: F sweep — how frags/precursor affects each approach
    # ====================================================================
    # This is the key table. D's bound narrowing amortizes log(M) across F
    # fragments per precursor group. With large F, D approaches O(T+M).
    # A/B don't depend on F (they sort globally). C doesn't either (no grouping).
    println("  Table 2: F (frags/precursor) sweep, fixed T=4000")
    println("  D's cost = O(T + P·log M) where P = T/F. Larger F = fewer binary searches.")
    println()
    print_matching_header()

    T_fixed = 4000
    F_values = [2, 5, 10, 20, 50, 100]

    for M in M_sizes
        for F in F_values
            P = T_fixed ÷ F
            peaks, tu, ts, pb, N = generate_matching_data(M, P, F; ppm_tol=ppm_tol)
            tA, tB, tC, tD, nm, _ = bench_matching_cell(
                work, radix_scratch, radix_counts, peaks, tu, ts, pb, N, ppm_tol)
            @printf("    F=%-3d P=%-4d M=%-4d (%2d%% m=%4d) %8.1f %9.1f %9.1f %9.1f\n",
                    F, P, M, round(Int, 100*nm/N), nm, tA, tB, tC, tD)
        end
        println()
    end

    # ====================================================================
    # Table 3: Large T regime (T=20000) — does sort cost dominate?
    # ====================================================================
    println("  Table 3: Large T=20000, F=20 — sort cost O(T·log T) should dominate A")
    println()
    print_matching_header()

    for M in M_sizes
        P = 1000  # 20000 / 20
        peaks, tu, ts, pb, N = generate_matching_data(M, P, 20; ppm_tol=ppm_tol)
        tA, tB, tC, tD, nm, _ = bench_matching_cell(
            work, radix_scratch, radix_counts, peaks, tu, ts, pb, N, ppm_tol)
        @printf("    T=20000 M=%-5d (%2d%% match=%5d) %7.1f %9.1f %9.1f %9.1f\n",
                M, round(Int, 100*nm/N), nm, tA, tB, tC, tD)
    end
    println()

    # ====================================================================
    # End-to-end savings estimate
    # ====================================================================
    println("  " * "="^70)
    println("  End-to-end estimate: T=4000, M=300, F=20")
    println("  40k calls/file × 3 files = 120k calls")
    println("  " * "="^70)
    println()

    peaks, tu, ts, pb, N = generate_matching_data(300, 200, 20; ppm_tol=ppm_tol)
    tA, tB, tC, tD, _, _ = bench_matching_cell(
        work, radix_scratch, radix_counts, peaks, tu, ts, pb, N, ppm_tol; samples=500)

    total_calls = 40_000 * 3
    for (name, t) in [("A. QuickSort + merge-scan", tA),
                       ("B. Radix sort + merge-scan", tB),
                       ("C. Binary search (unsorted)", tC),
                       ("D. Binary search + narrowing", tD)]
        total_s = t * total_calls / 1e6
        savings_s = (tA - t) * total_calls / 1e6
        @printf("    %-34s %6.1f μs/call → %5.2f s total (savings: %+.2f s vs A)\n",
                name, t, total_s, savings_s)
    end

    println()
    best = argmin([tA, tB, tC, tD])
    best_name = ["A", "B", "C", "D"][best]
    best_t = [tA, tB, tC, tD][best]
    @printf("    Winner: %s (%.1f μs). ", best_name, best_t)
    @printf("D vs C: %.2fx. D vs B: %.2fx.\n", tC/max(tD,0.001), tB/max(tD,0.001))
    println()
end


# ============================================================================
# Main
# ============================================================================

function main()
    println("\n" * "=" ^ 70)
    println("  K-Way Merge Benchmark for Transition Sorting")
    println("=" ^ 70 * "\n")

    # Part 1: Check real library if available
    lib_dir = "./data/ecoli_test/altimeter_ecoli.poin"
    try
        verify_library_sorting(lib_dir)
    catch e
        println("  Skipping Part 1 (library verification): $e")
        println()
    end

    # Parts 2-3: Synthetic benchmarks
    run_benchmarks()
    run_manual_merge_benchmarks()

    # Part 4: Feasibility estimate
    estimate_savings()

    # Part 5: Real library benchmark
    try
        benchmark_with_real_library(lib_dir)
    catch e
        println("  Skipping Part 5: $e")
        println()
    end

    # Part 6: Alternative sort strategies
    run_sort_strategy_benchmarks()

    # Part 7: Sort-free matching via binary search
    run_matching_strategy_benchmarks()

    println("Done.")
end

main()
