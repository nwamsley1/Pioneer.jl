# FirstPassSearch Deconvolution Profile Analysis

**Date:** 2026-03-11
**Dataset:** OlsenAstralThreeProteome200ng (1 file)
**Config:** 10 threads, ~26K scans/thread, 27.5M PSMs produced
**Total deconvolution time:** 27.4s
**Profiler:** Julia `Profile` at 100us sampling interval, ~25,844 application samples

---

## How to View the Interactive Flamegraph

The profile was saved as a PProf protobuf file on your Desktop. To open an interactive flamegraph in your browser:

```julia
# In your existing Julia session:
PProf.pprof("/Users/nathanwamsley/Desktop/pioneer_firstpass_profile.pb.gz"; web=true)

# Or from a fresh terminal:
julia -e 'using PProf; PProf.pprof("/Users/nathanwamsley/Desktop/pioneer_firstpass_profile.pb.gz"; web=true)'
```

This starts a local web server (typically `http://localhost:62261`) and opens the flamegraph in your browser. You can click to zoom into callstacks, search by function name, and toggle between flame/graph views.

The raw text profiles are also available:
- `/Users/nathanwamsley/Desktop/pioneer_firstpass_profile.txt` — tree-format call stack
- `/Users/nathanwamsley/Desktop/pioneer_firstpass_profile_flat.txt` — flat profile sorted by total count

---

## Executive Summary

The deconvolution pipeline (`perform_second_pass_search`) breaks down as follows:

| Component | Total Samples | Self Samples | % of Runtime | Key Files |
|-----------|--------------|-------------|-------------|-----------|
| **selectTransitions + fillTransitionList** | ~11,074 | ~4,000+ | **43%** | selectTransitions.jl, fillTransitionList.jl |
| **getDistanceMetrics (scoring)** | ~5,748 | ~2,400+ | **22%** | spectralDistanceMetrics.jl |
| **libraryBSpline B()** | ~3,630 | ~3,630 | **14%** | libraryBSpline.jl |
| **buildDesignMatrix + sortSparse** | ~2,396 | ~500+ | **9%** | buildDesignMatrix.jl |
| **partitionScansToThreads** | ~1,939 | ~1,400+ | **7.5%** | partitionThreadTasks.jl |
| **DataFrame construction** | ~1,794 | ~1,600+ | **7%** | Tables.jl fallbacks |
| **matchPeaks** | ~1,453 | ~600+ | **6%** | matchPeaks.jl |
| **solveOLS** | ~508 | ~435 | **2%** | spectralLinearRegression.jl |
| **Score! + ScoreFragmentMatches** | ~313 | ~300+ | **1.2%** | ScoredPSMs.jl, UnscoredPSMs.jl |
| JIT Compiler overhead | ~920 | ~920 | 3.5% | (background compilation) |

Note: Percentages sum to >100% because B-spline time is counted within selectTransitions, and some functions are called from multiple places. The "Self Samples" column shows time in the function's own code (not callees).

---

## Detailed Findings

### 1. selectTransitions Pipeline (43% total) — Largest Bottleneck

The transition selection pipeline is called once per scan (~260K total calls across threads) and consists of:

**Call chain:**
```
selectTransitions!  (selectTransitions.jl:76)
  -> _select_transitions_impl!  (standardTransitionSelection.jl:65)
       -> fillTransitionList!  (fillTransitionList.jl:70)
            -> fillTransitionListPrecomputed!  (fillTransitionList.jl:91)
                 -> getFragIsotopes!  (fillTransitionList.jl:174)
                      -> getFragAbundance!  (isotopeSplines.jl:326)
                           -> isotope()  (isotopeSplines.jl:136)
                                -> splevl()  (libraryBSpline.jl:51)
                                     -> B()  (libraryBSpline.jl:23)  ← HOTTEST LEAF
                 -> addTransitionIsotopes!  (fillTransitionList.jl:148)
  -> sort!(@view(...), by=getMZ)  (selectTransitions.jl:90)
```

**Top self-time functions in this chain:**

| Function | Self Samples | Description |
|----------|-------------|-------------|
| `B()` (libraryBSpline.jl:23) | 3,630 | Recursive B-spline basis evaluation |
| `sort!` / `partition!` / `smallsort!` | ~800 | Sorting transitions by m/z after selection |
| `addTransitionIsotopes!` (fillTransitionList.jl:155-167) | ~300 | DetailedFrag struct creation per isotope |
| `fillTransitionListPrecomputed!` inner loop | ~200 | Fragment iteration overhead |
| `ensureTransitionCapacity!` (selectTransitions.jl:51) | ~290 | Array resize checks |

#### B-spline: The single hottest leaf function

The recursive `B()` function at `src/utils/ML/libraryBSpline.jl:23-44` is the single largest self-time contributor at **3,630 samples (14% of total)**.

```julia
function B(x::T, k::Int, i::Int, t::NTuple{N,T}) where {N,T<:AbstractFloat}
    if k == 0
        return T(t[i] ≤ x < t[i+1])
    end
    c1 = if t[i+k] == t[i]
        zero(T)
    else
        ((x - t[i]) / (t[i+k] - t[i])) * B(x, k-1, i, t)    # recursive
    end
    c2 = if t[i+k+1] == t[i+1]
        zero(T)
    else
        ((t[i+k+1] - x) / (t[i+k+1] - t[i+1])) * B(x, k-1, i+1, t)  # recursive
    end
    return c1 + c2
end
```

For degree k=3 (cubic), each call to `B(x, 3, i, t)` generates 2^3 = 8 recursive calls down to k=0. The `splevl()` function calls `B()` for each of `n` basis functions, making the total cost O(n * 2^k) per spline evaluation. This is called via `getFragAbundance!` -> `isotope()` -> `splevl()` for every fragment of every precursor of every scan.

**Optimization opportunity:** Replace recursive evaluation with the iterative de Boor algorithm, which evaluates all basis functions in O(n * k) instead of O(n * 2^k). This alone could save ~10-14% of total runtime.

#### Sort after selection

After filling the transition list, `selectTransitions!` calls `sort!(@view(transitions[1:transition_idx]), by=getMZ)` (line 90). The sort uses PartialQuickSort and contributes ~800 samples. This is called once per scan.

---

### 2. getDistanceMetrics / computeFittedMetricsFor (22% total)

**File:** `src/Routines/SearchDIA/PSMs/spectralDistanceMetrics.jl`

Called once per precursor per scan to compute 9 spectral similarity metrics. The outer function `getDistanceMetrics` (line 216) iterates over all columns (precursors) in the sparse design matrix, and for each one calls `computeFittedMetricsFor` (line 302) in a while loop that iteratively drops the worst-matching peak.

**Self-time breakdown:**

| Function/Line | Self | What it does |
|---------------|------|-------------|
| `computeFittedMetricsFor` main loop (lines 347-400) | ~1,236 | 15+ accumulators, sqrt/log per fragment |
| `computeFittedMetricsFor` third loop (lines 403-407) | ~400 | Recomputes sqrt(fitted_peak)/h_sqrt_sum for scribe |
| `getDistanceMetrics` residual init (lines 231-247) | ~350 | Sparse matrix-vector product for residuals |
| `getDistanceMetrics` outer while loop (lines 269-283) | ~200 | Iterative peak dropping with `deleteat!` |

**Key observations:**
- Three separate passes over `included_indices`: (1) compute total_h/total_x (lines 314-322), (2) main accumulator loop (lines 347-400), (3) scribe sqrt loop (lines 403-407). The first and third could be folded into the second.
- `sqrt(fitted_peak)` and `sqrt(shadow_peak)` are computed in the main loop (line 365-366) but the scribe normalization (line 406) recomputes them.
- `deleteat!(incl, worst_pos)` on line 282 causes array shifting. Since the while loop typically runs only 1-3 iterations, this is minor.
- The `incl = Int[]` allocation on line 254 with `push!` in the inner loop creates garbage per precursor. A pre-allocated buffer would be better.

---

### 3. buildDesignMatrix + sortSparse (9% total)

**File:** `src/Routines/SearchDIA/CommonSearchUtils/buildDesignMatrix.jl`

Constructs the sparse design matrix mapping fragments to precursors for each scan.

| Function | Total | Self |
|----------|-------|------|
| `buildDesignMatrix!` (line 80) | 2,396 | ~500 |
| `sortSparse!` (line 170, called at end) | 2,011 | ~400 |

The `sortSparse!` at the end of `buildDesignMatrix!` sorts each column of the sparse matrix. The custom sort implementation uses partition/smallsort for small arrays. This is called once per scan.

---

### 4. partitionScansToThreads (7.5% total — one-time cost)

**File:** `src/Routines/SearchDIA/CommonSearchUtils/partitionThreadTasks.jl`

This runs **once** at the start of `perform_second_pass_search` to distribute scans across threads. Despite being one-time, it takes 1,939 samples (~7.5%).

**Root cause (line 48-59):**
```julia
spectra_ids = collect([x for x in range(1, length(spectra)) if ms_order[x]==2])
for i in range(1, length(spectra_ids))
    if rt[i] - rt[bin_start] > 1.0f0
        sort!(@view(spectra_ids[bin_start:bin_stop]), by = x->round_float32_alt(prec_mz[x],6))
        ...
    end
end
```

The `round_float32_alt` function (1,016 self samples on line 53) does `Float32(round(x; digits=decimals))` in the sort comparison, which is expensive. This anonymous closure is called O(n log n) times per RT bin for the sort. Precomputing the rounded values and sorting by a lookup would be much faster.

Additionally, line 67-76 distributes scans with a stride-10 interleaving pattern using nested loops, which itself contributes ~300 samples.

---

### 5. DataFrame construction (7% total — end-of-thread)

**File:** `SecondPassSearch/utils.jl:398`

```julia
return (psms = DataFrame(@view(getComplexScoredPsms(search_data)[1:last_val])), ...)
```

This converts a `SubArray{ComplexScoredPSM}` of ~2.7M structs into a DataFrame. The Tables.jl fallback `_buildcolumns` iterates the SubArray of structs and extracts each field into column vectors — 1,794 samples.

**Optimization opportunity:** Return the raw scored PSM vector and defer DataFrame construction, or use a struct-of-arrays representation from the start.

---

### 6. matchPeaks (6% total)

**File:** `src/Routines/SearchDIA/CommonSearchUtils/matchPeaks.jl`

The single-pass O(T+P) matching algorithm is well-optimized with `@inbounds @fastmath`. Self-time is distributed across:

| Function | Self | Line |
|----------|------|------|
| `setMatch!` (struct construction) | ~370 | lines 49-90 |
| `setNearest!` (peak scanning) | ~160 | lines 149-214 |
| `matchPeaks!` main loop | ~140 | lines 346-392 |

No obvious low-hanging fruit — the algorithm is already efficient.

---

### 7. solveOLS (2% total — surprisingly fast)

**File:** `src/utils/ML/spectralLinearRegression.jl:310-371`

The coordinate descent solver contributes only 435 self-time samples. The inner loops at lines 334-336 and 345-347 are tight `@inbounds @fastmath` loops over the sparse matrix. The previous `time_ns()` timing showed ~33% for "solveOLS" but that bucket included `buildDesignMatrix` + `initResiduals` overhead.

The solver itself is very efficient — no optimization needed here.

---

### 8. Score! + ScoreFragmentMatches! (1.2% total)

**Files:** `PSMs/ScoredPSMs.jl:207-289`, `PSMs/UnscoredPSMs.jl:81-200`

Minimal overhead. The `Score!` function constructs `ComplexScoredPSM` structs with 25+ fields, and `ScoreFragmentMatches!` / `ModifyFeatures!` extracts per-fragment features. Well-optimized.

---

## Optimization Recommendations (Ranked by Expected Impact)

### Priority 1: Iterative B-spline evaluation
**Expected speedup: 10-14% of total runtime**
**File:** `src/utils/ML/libraryBSpline.jl:23-64`

Replace the recursive `B()` with an iterative de Boor algorithm. The current implementation does O(n * 2^k) work per `splevl()` call; an iterative approach does O(n * k). For k=3, this is an 8x/3x = 2.7x speedup in spline evaluation, which at 14% of total runtime translates to ~10% overall.

The iterative approach builds up basis values from k=0 to k=3 using a small working array of size k+1, reusing intermediate values instead of recomputing them.

### Priority 2: Fuse loops in computeFittedMetricsFor
**Expected speedup: 5-8% of total runtime**
**File:** `src/Routines/SearchDIA/PSMs/spectralDistanceMetrics.jl:302-431`

Currently three passes over `included_indices`:
1. Lines 314-322: Compute `total_h`, `total_x`, `num_matching_peaks`
2. Lines 347-400: Main 15-accumulator loop (also computes `h_sqrt_sum`, `x_sqrt_sum`)
3. Lines 403-407: Recompute `sqrt(fitted_peak)/h_sqrt_sum` for scribe score

Fuse pass 1 into pass 2 (accumulate `total_h/total_x` alongside everything else). For the scribe score, cache `sqrt(fitted_peak)` values from pass 2 in a small buffer and compute scribe in a second tight pass without re-accessing H, or accept a two-pass approach but eliminate the redundant `w[col]*H.nzval[i]` computation.

Also: replace `incl = Int[]` + `push!` (line 254-258) with a pre-allocated buffer.

### Priority 3: Avoid per-thread DataFrame materialization
**Expected speedup: 3-7% of total runtime**
**File:** `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl:398`

Instead of `DataFrame(@view(getComplexScoredPsms(search_data)[1:last_val]))`, return the raw `ComplexScoredPSM` vector (or a copy of the view) and build the DataFrame once after all threads complete and merge. This avoids the expensive Tables.jl struct-of-arrays conversion per thread.

### Priority 4: Precompute sort keys in partitionScansToThreads
**Expected speedup: 2-4% of total runtime (one-time)**
**File:** `src/Routines/SearchDIA/CommonSearchUtils/partitionThreadTasks.jl:48-59`

The sort comparison `by = x->round_float32_alt(prec_mz[x], 6)` calls `Float32(round(x; digits=6))` on every comparison. Precompute the rounded m/z values into a vector and sort by direct lookup:

```julia
rounded_mz = [round_float32_alt(prec_mz[i], 6) for i in spectra_ids]
sort!(@view(spectra_ids[bin_start:bin_stop]), by = i -> rounded_mz[i])
```

### Priority 5: Reduce fillTransitionList overhead
**Expected speedup: 1-3% of total runtime**
**File:** `src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/fillTransitionList.jl:148-172`

In `addTransitionIsotopes!`, precompute `NEUTRON/frag.frag_charge` once before the isotope loop instead of computing it per isotope. The `ensureTransitionCapacity!` check on line 169 could be moved outside the loop or batched.

---

## What's NOT a Bottleneck

- **solveOLS** (2%): The coordinate descent solver is tight and efficient
- **Score! / ScoreFragmentMatches!** (1.2%): Minimal overhead
- **matchPeaks** (6%): Already well-optimized with `@inbounds @fastmath`
- **JIT compilation** (3.5%): Background noise from Julia compiler, will decrease on subsequent runs

---

## Comparison with Previous time_ns() Timing

The previous per-section timing showed different proportions because each bucket included overhead from its surrounding code:

| time_ns() bucket | Reported % | Profile finding |
|-----------------|-----------|-----------------|
| selectTransitions | ~30% | 43% (includes B-spline, isotope calc) |
| solveOLS | ~33% | 2% solver + 9% buildDesignMatrix + initResiduals |
| scoring | ~21% | 22% getDistanceMetrics + 1.2% Score! |
| matchPeaks | ~8% | 6% |
| buildMatrix | ~8% | 9% (now separated from solveOLS bucket) |

The profile reveals that `time_ns()` lumped `buildDesignMatrix` + `initResiduals` into the "solveOLS" bucket (they were timed together between `_t0 = time_ns()` and `t_solve_ols += time_ns() - _t0`). The actual solver is very fast; the matrix construction is the expensive part of that pipeline stage.
