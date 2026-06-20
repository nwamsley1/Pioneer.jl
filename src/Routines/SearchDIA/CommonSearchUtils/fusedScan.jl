# Fused per-precursor scan helpers.
#
# Low-level building blocks for the fused match+build+score path used by
# MainSearch. Per-precursor invariants (m/z-sorted detailed fragments) are
# enforced once at search setup via `verify_mz_sorted`.

"""
    FusedScratch

Per-thread scratch for one precursor's matches + misses. Sized at init
for the largest single-precursor fragment list the library could produce;
reused across all precursors and scans.

Entries are written out-of-order within a precursor (one per matched or
missed fragment-isotope pair). `finalize_column!` sorts them by `row` if
non-monotonic, dedupes duplicates (same peak, multiple fragments), and
flushes to a `SparseArrayFused` in CSC order.

Miss rows are assigned unique synthetic row IDs (> n_peaks) by the caller
so they never collide with matched peak indices.
"""
mutable struct FusedScratch
    # Match entries: one per matched (fragment, iso) pair.
    row::Vector{UInt32}
    nzval::Vector{Float32}
    x::Vector{Float32}
    isotope::Vector{UInt8}
    rank::Vector{UInt8}
    # Miss entries: buffered here and flushed alongside matches in finalize.
    # Misses from a precursor with zero matches are discarded.
    miss_nzval::Vector{Float32}
    miss_isotope::Vector{UInt8}
    miss_rank::Vector{UInt8}
    n::Int
    miss_n::Int
end

function FusedScratch(max_per_prec::Int)
    FusedScratch(
        Vector{UInt32}(undef, max_per_prec),
        Vector{Float32}(undef, max_per_prec),
        Vector{Float32}(undef, max_per_prec),
        Vector{UInt8}(undef, max_per_prec),
        Vector{UInt8}(undef, max_per_prec),
        Vector{Float32}(undef, max_per_prec),
        Vector{UInt8}(undef, max_per_prec),
        Vector{UInt8}(undef, max_per_prec),
        0, 0,
    )
end

@inline reset_fused_scratch!(s::FusedScratch) = (s.n = 0; s.miss_n = 0; nothing)

"""
    grow_fused_scratch!(s, needed)

Expand all entry-indexed vectors so they can hold at least `needed`
entries in either the match or miss buffer. Called from `run_fused!`
when a precursor's fragment list exceeds the initial capacity.
"""
function grow_fused_scratch!(s::FusedScratch, needed::Int)
    cap = length(s.row)
    if needed > cap
        new_cap = max(needed, cap * 2)
        extra = new_cap - cap
        append!(s.row,          zeros(UInt32, extra))
        append!(s.nzval,        zeros(Float32, extra))
        append!(s.x,            zeros(Float32, extra))
        append!(s.isotope,      zeros(UInt8, extra))
        append!(s.rank,         zeros(UInt8, extra))
        append!(s.miss_nzval,   zeros(Float32, extra))
        append!(s.miss_isotope, zeros(UInt8, extra))
        append!(s.miss_rank,    zeros(UInt8, extra))
    end
    return
end

# ── SIMD primitives (mirror structs/SpectralLibrary/PartitionedFragmentIndex/search.jl) ──

const _FUSED_F32x8 = NTuple{8, Core.VecElement{Float32}}

@inline _fused_vbroadcast8(x::Float32)::_FUSED_F32x8 =
    ntuple(_ -> Core.VecElement(x), Val(8))

@inline _fused_vload8(arr::AbstractVector{Float32}, i::Int)::_FUSED_F32x8 =
    unsafe_load(Ptr{_FUSED_F32x8}(pointer(arr, i)))

@inline _fused_vcmpge_mask(a::_FUSED_F32x8, b::_FUSED_F32x8)::UInt8 =
    Core.Intrinsics.llvmcall("""
        %cmp  = fcmp oge <8 x float> %0, %1
        %mask = bitcast <8 x i1> %cmp to i8
        ret i8 %mask
    """, UInt8, Tuple{_FUSED_F32x8, _FUSED_F32x8}, a, b)

@inline _fused_vcmple_mask(a::_FUSED_F32x8, b::_FUSED_F32x8)::UInt8 =
    Core.Intrinsics.llvmcall("""
        %cmp  = fcmp ole <8 x float> %0, %1
        %mask = bitcast <8 x i1> %cmp to i8
        ret i8 %mask
    """, UInt8, Tuple{_FUSED_F32x8, _FUSED_F32x8}, a, b)

"""
    simd_find_first_ge(a, lo, hi, target) -> Int

SIMD forward scan: first index `i ∈ lo:hi` with `a[i] >= target`, else `hi+1`.
Checks 8 Float32 elements per iteration; scalar tail.
"""
@inline function simd_find_first_ge(a::AbstractVector{Float32}, lo::Int, hi::Int,
                                     target::Float32)
    thr = _fused_vbroadcast8(target)
    i = lo
    @inbounds while i + 7 <= hi
        mask = _fused_vcmpge_mask(_fused_vload8(a, i), thr)
        mask != 0x00 && return i + trailing_zeros(mask)
        i += 8
    end
    @inbounds while i <= hi
        a[i] >= target && return i
        i += 1
    end
    return hi + 1
end

"""
    bsearch_hybrid(a, target, lo, hi) -> Int

Branchless binary iterations until the remaining window ≤ `FUSED_SIMD_CUTOFF`,
then SIMD forward scan. Lower-bound semantics; returns `hi+1` on no match.
"""
const FUSED_SIMD_CUTOFF = 32

@inline function bsearch_hybrid(a::AbstractVector{Float32}, target::Float32,
                                 lo::Int, hi::Int)
    len = hi - lo + 1
    @inbounds while len > FUSED_SIMD_CUTOFF
        half = len >> 1
        mid  = lo + half - 1
        lo   = a[mid] < target ? mid + 1 : lo
        len -= half
    end
    return simd_find_first_ge(a, lo, lo + len - 1, target)
end

"""
    scan_for_nearest(mz, lo, hi, target, ppm_tol) -> (peak_idx, ppm_err)

Scalar forward scan for the peak nearest `target` within ppm tolerance,
starting at `lo`. Returns `(0, 0f0)` if no peak is within tolerance. Called
after `bsearch_hybrid` has narrowed the window to ~1-3 peaks.
"""
@inline function scan_for_nearest(mz::AbstractVector{Float32}, lo::Int, hi::Int,
                                   target::Float32, ppm_tol::Float32)
    tol = target * ppm_tol * 1f-6
    target_hi = target + tol
    best_peak = 0
    best_err = typemax(Float32)
    @inbounds while lo <= hi && mz[lo] <= target_hi
        err = abs(mz[lo] - target)
        if err <= tol && err < best_err
            best_err = err
            best_peak = lo
        end
        lo += 1
    end
    if best_peak > 0
        @inbounds ppm_err = (target - mz[best_peak]) / (target * 1f-6)
        return best_peak, ppm_err
    else
        return 0, 0f0
    end
end

"""
    insertion_sort_scratch!(s, n)

In-place insertion sort of the first `n` entries of a `FusedScratch` by
ascending `row`. Used by `finalize_column!` when the per-precursor match
sequence turns out to be non-monotonic (rare in practice).
"""
@inline function insertion_sort_scratch!(s::FusedScratch, n::Int)
    @inbounds for i in 2:n
        r = s.row[i]; nv = s.nzval[i]; xv = s.x[i]; iv = s.isotope[i]; rv = s.rank[i]
        j = i - 1
        while j >= 1 && s.row[j] > r
            s.row[j+1]     = s.row[j]
            s.nzval[j+1]   = s.nzval[j]
            s.x[j+1]       = s.x[j]
            s.isotope[j+1] = s.isotope[j]
            s.rank[j+1]    = s.rank[j]
            j -= 1
        end
        s.row[j+1] = r; s.nzval[j+1] = nv; s.x[j+1] = xv; s.isotope[j+1] = iv; s.rank[j+1] = rv
    end
    nothing
end

"""
    finalize_column!(H, col, s, entry, miss_row_start) -> (new_entry, new_miss_row)

Flush a per-precursor `FusedScratch` to `SparseArrayFused` `H` as column
`col`, in CSC order:
  1. Sort scratch by row (insertion sort when non-monotonic).
  2. Assign each miss a fresh synthetic row ID starting from
     `miss_row_start + 1` (so misses never collide with peak indices or
     with misses from other precursors), write them after matches.
  3. Dedupe consecutive equal rows (multiple fragments landing on the
     same peak) by accumulating `nzval`; isotope of last writer wins.
  4. Write `H.colptr[col]` to the first entry of this column.

Returns `(updated_entry, updated_miss_row)`.
"""
function finalize_column!(H::SparseArrayFused{Ti,T}, col::UInt16,
                          s::FusedScratch, entry::Int,
                          miss_row_start::UInt32) where {Ti, T}
    n_match = s.n
    n_miss  = s.miss_n
    n_total = n_match + n_miss
    n_total == 0 && return (entry, miss_row_start)

    grow_if_needed!(H, entry + n_total)
    grow_colptr_if_needed!(H, Int(col))

    # Monotonicity check on matches.
    monotonic = true
    @inbounds for i in 2:n_match
        if s.row[i] < s.row[i-1]
            monotonic = false; break
        end
    end
    monotonic || insertion_sort_scratch!(s, n_match)

    H.colptr[col] = Ti(entry + 1)

    # Flush matches with in-pass dedup.
    prev_row = UInt32(0)
    prev_entry = 0
    @inbounds for i in 1:n_match
        r = s.row[i]
        if r == prev_row && prev_entry > 0
            H.nzval[prev_entry] += s.nzval[i]
            H.meta[prev_entry]   = pack_meta(true, s.isotope[i])
            H.rank[prev_entry]   = s.rank[i]
        else
            entry += 1
            H.rowval[entry] = Ti(r)
            H.nzval[entry]  = s.nzval[i]
            H.x[entry]      = s.x[i]
            H.meta[entry]   = pack_meta(true, s.isotope[i])
            H.rank[entry]   = s.rank[i]
            prev_row = r
            prev_entry = entry
        end
    end

    # Flush misses with synthetic row IDs (> all peak indices, all unique).
    miss_row = miss_row_start
    @inbounds for i in 1:n_miss
        miss_row += UInt32(1)
        entry += 1
        H.rowval[entry] = Ti(miss_row)
        H.nzval[entry]  = s.miss_nzval[i]
        H.x[entry]      = zero(T)
        H.meta[entry]   = pack_meta(false, s.miss_isotope[i])
        H.rank[entry]   = s.miss_rank[i]
    end

    return (entry, miss_row)
end


"""
    verify_mz_sorted(fragments, prec_ranges)

Assert that each precursor's fragment block in `fragments` is sorted by
ascending m/z. Called once at MainSearch setup when `use_fused_scan=true`.
`prec_ranges[k]` is the first fragment index for precursor `k`; precursor
`k`'s range is `prec_ranges[k]:prec_ranges[k+1]-1`.

Throws when fragments aren't m/z-sorted, telling the user to rebuild
the library with current `BuildSpecLib(...)` (which sorts at build time).
"""
function verify_mz_sorted(fragments, prec_ranges::AbstractVector)
    n_prec = length(prec_ranges) - 1
    @inbounds for pid in 1:n_prec
        lo = Int(prec_ranges[pid])
        hi = Int(prec_ranges[pid + 1]) - 1
        lo > hi && continue
        prev = fragments[lo].mz
        for i in (lo + 1):hi
            cur = fragments[i].mz
            if cur < prev
                error("Library is not m/z-sorted per precursor (precursor $pid, " *
                      "fragment indices $lo:$hi). The search path requires each " *
                      "precursor's fragments to be ascending in m/z. Rebuild the " *
                      "library from FASTA with current Pioneer — `BuildSpecLib(...)` " *
                      "sorts detailed fragments by m/z at build time.")
            end
            prev = cur
        end
    end
    return nothing
end
