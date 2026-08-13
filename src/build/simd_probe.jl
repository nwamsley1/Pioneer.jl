# Cross-platform probe for Pioneer's hand-written SIMD primitives.
#
# Pioneer emits explicit LLVM vector IR in three places (fragment-index and fused scan):
#   src/structs/SpectralLibrary/PartitionedFragmentIndex/search.jl   (_vcmpge_mask)
#   src/Routines/SearchDIA/CommonSearchUtils/fusedScan.jl            (_fused_vcmpge_mask, _fused_vcmple_mask)
# all of the shape  `fcmp oge <8 x float>` -> `bitcast <8 x i1> to i8`, i.e. a FIXED 256-bit vector.
#
# That width is native on AVX/AVX2 (one instruction) but not on ARM NEON (128-bit, so LLVM emits
# two operations) and does not widen on AVX-512. This script answers, per machine:
#
#   1. CORRECTNESS  -- do the SIMD scanners agree with a scalar reference? (the alignment and
#      mask-lowering questions on AArch64 are checked here, not assumed)
#   2. SPEED        -- is the SIMD path actually beating scalar on this hardware?
#   3. WIDTH        -- would a different vector width do better here?
#
# Run:  julia --project=. src/build/simd_probe.jl
using Pioneer, Printf
const P = Pioneer

# ---------------------------------------------------------------- reference implementation
scalar_find_first_ge(a, lo, hi, t) = begin
    @inbounds for i in lo:hi
        a[i] >= t && return i
    end
    return hi + 1
end

# ---------------------------------------------------------------- 1. correctness
function check_correctness()
    ok = true
    for n in (1, 2, 7, 8, 9, 15, 16, 17, 1000, 4096, 10_007)
        a = sort!(rand(Float32, n))
        for t in (-1.0f0, 0.0f0, 0.25f0, 0.5f0, 0.75f0, 1.0f0, 2.0f0)
            got = P.simd_find_first_ge(a, 1, n, t)
            ref = scalar_find_first_ge(a, 1, n, t)
            if got != ref
                @printf("  MISMATCH n=%d t=%.3f: simd=%d scalar=%d\n", n, t, got, ref)
                ok = false
            end
        end
        # unaligned starts: the SIMD load is 32 bytes from an arbitrary element offset
        for lo in 2:min(n, 9)
            t = 0.5f0
            got = P.simd_find_first_ge(a, lo, n, t)
            ref = scalar_find_first_ge(a, lo, n, t)
            got == ref || (@printf("  MISMATCH unaligned lo=%d n=%d\n", lo, n); ok = false)
        end
    end
    println(ok ? "  correctness: PASS" : "  correctness: *** FAIL ***")
    return ok
end

# ---------------------------------------------------------------- 2. speed vs scalar
# The accumulator matters: without consuming the result, LLVM deletes the scalar loop entirely and
# it benchmarks at 0 ns. `sink` is a non-const global so the compiler cannot prove the value unused.
sink = 0

@noinline function drive_scalar(a, n, t, reps)
    acc = 0
    @inbounds for _ in 1:reps
        acc += scalar_find_first_ge(a, 1, n, t)
    end
    return acc
end

@noinline function drive_simd(a, n, t, reps)
    acc = 0
    @inbounds for _ in 1:reps
        acc += Int(P.simd_find_first_ge(a, 1, n, t))
    end
    return acc
end

function bench(f, reps)
    global sink
    sink += f(reps)                          # consume, so nothing is elided
    best = Inf
    for _ in 1:5
        t0 = time_ns()
        sink += f(reps)
        best = min(best, (time_ns() - t0) / reps)
    end
    return best
end

function check_speed()
    for n in (64, 256, 4096, 65_536)
        a = sort!(rand(Float32, n))
        t = 0.997f0                       # near the end: forces a full scan
        reps = max(200, 2_000_000 ÷ n)
        s = bench(r -> drive_scalar(a, n, t, r), reps)
        v = bench(r -> drive_simd(a, n, t, r), reps)
        @printf("  n=%6d   scalar %8.1f ns   simd %8.1f ns   speedup %5.2fx\n", n, s, v, s / v)
    end
end

function main()
    println("host: ", Sys.MACHINE, "  CPU: ", Sys.CPU_NAME, "  threads: ", Threads.nthreads())
    println("julia ", VERSION)
    println("\n1. correctness of the hand-written SIMD scanners")
    check_correctness()
    println("\n2. SIMD vs scalar (a speedup near 1.0 means the vector path is not paying off here)")
    check_speed()
    println("\nReport back: the host/CPU line, PASS/FAIL, and the speedup column.")
end

main()
