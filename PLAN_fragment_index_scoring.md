# Plan: Likelihood-Ratio Fragment Scoring for Fragment Index

## Context

Replace rank-based fragment scoring with likelihood-ratio weights derived from predicted intensities. The UInt8 `score` field in `SimpleFrag`/`IndexFragment` is recomputed at library build time; search-time code (`Counter`, `inc!`, `filterPrecursorMatches!`) is unchanged.

**Focus**: Altimeter/SplineCoefficientModel path. Spline coefficients evaluated at reference NCE (default 27) to get representative intensities for LR scoring.

**Both scoring modes must remain available** via `"fragment_scoring": "rank_based"` or `"likelihood_ratio"` in the config JSON. This allows building libraries with either method for comparison.

## Current State (already committed on `feature/fragment-index-scoring`)

Done:
- `fragment_likelihood_weights_uint8()` in `build_poin_lib.jl` (uses OLD normalization — needs update)
- `getSimpleFrags` dual-path (rank-based vs LR) with pre-allocated buffers
- SplineCoefficientModel: spline evaluation via `splevl()` at reference NCE
- KoinaModelType: passes `fragments_table[:intensity]` directly
- Config plumbing: `BuildSpecLib.jl` reads LR params from JSON, passes to `buildPionLib`
- `defaultBuildLibParams.json` has LR config fields (default `"rank_based"`)

Still needed: update normalization in `fragment_likelihood_weights_uint8()`.

## Normalization Strategy: Per-Precursor + Max Clamp

### Approach

1. Compute raw Δ_i for each fragment using the Binomial detection model
2. Clamp negatives to 0
3. Scale so sum of all Δ_i = 255 (per-precursor normalization)
4. Clamp max individual score to `K × floor(255/M) - 1` where K is min required fragment matches
5. Fragments with Δ_raw = 0 get score **0** (undetectable = no evidence)
6. If all Δ_raw = 0, fall back to equal weights: each = floor(255/M)

### Threshold

- `threshold = K × floor(255/M)` where K=5, M=7 → **threshold = 180**
- Max individual fragment score = 179
- **No single fragment match can pass the threshold** — enforces multi-fragment evidence
- A precursor with only 1 detectable fragment (score ≤ 179) cannot pass the pre-screen

### Score examples (N=200, d=10, α=0.01, M=7, K=5, max_score=179)

| Distribution | UInt8 weights (post-clamp) | Sum | Pass 180? |
|---|---|---|---|
| Even [1/7 each] | [36, 36, 36, 36, 36, 36, 36] | 252 | 5+ frags ✓ |
| All-in-one [1,0,...,0] | [179, 0, 0, 0, 0, 0, 0] | 179 | Never ✗ |
| Slight gradient | [70, 52, 42, 34, 25, 18, 13] | 254 | 3+ frags ✓ |
| Typical [.40,.25,.15,.10,.05,.03,.02] | [114, 67, 35, 22, 12, 5, 0] | 255 | 2+ frags ✓ |
| Extreme skew [.70,.15,.08,.04,...] | [179, 40, 20, 10, 0, 0, 0] | 249 | 2+ frags ✓ |

### Why this works

- Uses full 0-255 UInt8 range for individual scores
- Good differentiation between strong and weak fragments
- Sum always ≤ 255 → no Counter overflow, no search-time changes
- The max clamp prevents false positives from a single lucky fragment match
- Precursors with very few detectable fragments are correctly filtered out

## Config Parameters

In `defaultBuildLibParams.json` under `library_params`:

| Parameter | Default | Description |
|---|---|---|
| `fragment_scoring` | `"rank_based"` | `"rank_based"` or `"likelihood_ratio"` |
| `likelihood_ratio_M` | `7` | Target fragments per precursor |
| `likelihood_ratio_N` | `200` | Precursor ion count |
| `likelihood_ratio_d` | `10` | Detection limit in ions |
| `likelihood_ratio_alpha` | `0.01` | Random m/z coincidence probability |
| `likelihood_ratio_ref_nce` | `27.0` | Reference NCE for Altimeter splines |
| `likelihood_ratio_min_fragments` | `5` | K: minimum fragment matches for threshold |

The threshold `K × floor(255/M)` and max individual score `K × floor(255/M) - 1` are derived from K and M — not stored separately.

## Files to Modify

| File | Change |
|------|--------|
| `src/Routines/BuildSpecLib/build/build_poin_lib.jl` | Update `fragment_likelihood_weights_uint8` normalization, add `lr_K` param |
| `assets/example_config/defaultBuildLibParams.json` | Add `likelihood_ratio_min_fragments`, update N/d/α defaults |
| `src/Routines/BuildSpecLib.jl` | Read `likelihood_ratio_min_fragments`, update N/d/α fallback defaults, pass K to `buildPionLib` |

## Implementation

### Step 1: Update `fragment_likelihood_weights_uint8` in `build_poin_lib.jl`

```julia
function fragment_likelihood_weights_uint8(
    intensities::AbstractVector{<:Real};
    M::Int, N::Int, d::Int, α::Float64, K::Int = 5
)::Vector{UInt8}
    n_frags = length(intensities)
    n_frags == 0 && return UInt8[]

    total = sum(Float64, intensities)
    if total ≤ 0
        return fill(UInt8(max(1, 255 ÷ M)), n_frags)
    end
    p = Float64.(intensities) ./ total

    ε_floor = 1e-300  # avoid log(0)
    Δ = Vector{Float64}(undef, n_frags)
    for i in 1:n_frags
        p_miss = clamp(cdf(Binomial(N, p[i]), d - 1), ε_floor, 1.0 - ε_floor)
        p_detect = 1.0 - p_miss
        Δ[i] = max(0.0, log(p_detect / α) - log(p_miss / (1.0 - α)))
    end

    Δ_sum = sum(Δ)
    if Δ_sum ≤ 0.0
        return fill(UInt8(max(1, 255 ÷ M)), n_frags)
    end

    # Per-precursor normalization: scale so sum = 255
    scale = 255.0 / Δ_sum
    # Max individual score: ensures no single match passes the threshold
    max_individual = K * (255 ÷ M) - 1  # e.g., 5*36-1 = 179 for M=7

    return UInt8[clamp(round(Int, Δ[i] * scale), 0, max_individual) for i in 1:n_frags]
end
```

### Step 2: Plumb K parameter

Add `lr_K::Int = 5` keyword to:
- Both `buildPionLib` overloads (KoinaModelType and SplineCoefficientModel)
- `getSimpleFrags` function
- Pass through to `fragment_likelihood_weights_uint8(...; K=lr_K)`

In `BuildSpecLib.jl`:
- Read: `_lr_K = Int(get(_lp, "likelihood_ratio_min_fragments", 5))`
- Pass: `lr_K = _lr_K` in `buildPionLib` calls

### Step 3: Update config defaults

In `defaultBuildLibParams.json`:
```json
"likelihood_ratio_N": 200,
"likelihood_ratio_d": 10,
"likelihood_ratio_alpha": 0.01,
"likelihood_ratio_min_fragments": 5
```

In `BuildSpecLib.jl` get-fallback defaults: update N→200, d→10, α→0.01.

In `build_poin_lib.jl` `buildPionLib` overloads: update default keyword values to match.

### Step 4: Search threshold

The threshold = K × floor(255/M) = 5 × 36 = 180.

This goes in `min_index_search_score` in search configs. For now, keep existing defaults (which are for rank-based scoring). The threshold will be tuned when testing with LR scoring enabled.

## Verification

1. Unit test: `fragment_likelihood_weights_uint8` with even, all-in-one, gradient, realistic, and skewed distributions — verify sums ≤ 255, max ≤ 179, score ordering matches intensity ordering
2. Build library: `BuildSpecLib` with `"fragment_scoring": "likelihood_ratio"`
3. Search: `SearchDIA` with `min_index_search_score` = 180 on ecoli test data
4. Verify that all-in-one precursors (1 detectable fragment) are correctly filtered out
5. Compare PSM counts: rank-based vs LR scoring
