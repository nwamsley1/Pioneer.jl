# Plan: Likelihood-Ratio Fragment Scoring for Fragment Index

## Context

The current fragment index assigns UInt8 scores to fragments by intensity rank via a static `rank_to_score` vector (default `[8, 4, 4, 2, 2, 1, 1]`). This ignores the actual predicted intensity distribution — a fragment with 80% of total intensity gets the same score as one with 30%, as long as both are rank 1. The likelihood-ratio approach uses predicted intensities + instrument physics (M, N, d, α) to compute per-fragment weights (Δ_i) that reflect how much *evidence* each detection provides. This is Neyman-Pearson optimal for the binary detect/miss decision.

**Key constraint**: The existing search-time architecture (`Counter{UInt32, UInt8}`, `inc!`, `filterPrecursorMatches!`) is performance-critical and must not change. We only change how the UInt8 `score` field in `SimpleFrag`/`IndexFragment` is computed at library build time.

**Simplification**: We drop the per-precursor baseline B (sum of w⁻ᵢ) and store only Δ_i values. The baseline shifts all scores for a precursor by a constant; without it, we lose the penalty for missing dominant fragments, but the relative weight of each detection is preserved. This fits the existing additive Counter architecture with zero search-time changes.

## Current State: How fragment ranking works now

**Standard models (Prosit, UniSpec)**: Koina API returns `:intensities` column (Float32). Pioneer sorts fragments by descending intensity per precursor via `sort_fragments!()` (`fragment_predict.jl:275`). In `getSimpleFrags`, fragments are iterated in this pre-sorted order and assigned `rank_to_score[rank]` (rank 1 = highest intensity). `fragments_table[:intensity]` has real predicted intensities.

**Spline model (Altimeter)**: Koina API returns `:coefficients` (spline basis), NOT `:intensities`. The spline path does NOT call `sort_fragments!()`. In `process_spline_batch!` (`fragment_parse.jl:813`), `frag_intensity = rank` — a sequential counter, not predicted intensity. So `fragments_table[:intensity]` contains rank numbers, not real intensities. To get actual intensities, we evaluate the spline coefficients at a reference NCE using `splevl()` from `libraryBSpline.jl`.

The knot vector needed for spline evaluation is saved to `spline_knots.jls` in the library directory during `BuildSpecLib.jl:252-258`.

## Quantization Strategy: Global ε-clamped Normalization

### Why not per-precursor normalization?

Per-precursor normalization (scaling each precursor's Δ sum to 255) is wrong because it erases the total evidence strength. A precursor with 1 marginally-detectable fragment would get score 255 just like a precursor with 7 highly-detectable fragments. The whole point of likelihood-ratio scoring is that some precursors produce more total evidence than others.

### The problem with raw Δ values

For any reasonable N (e.g., 10000), the Binomial CDF gives astronomically small P(miss) values for dominant fragments. For example, a fragment with p = 0.3 of total intensity gives P(miss) ≈ 10^{-1300}. The log-likelihood ratio Δ then becomes enormous (~3000), far exceeding UInt8 range.

### Solution: Derive ε from the UInt8 budget

We clamp P(miss) at a floor ε, which caps Δ_max. We derive ε so that `M × Δ_max = 255`, where M is a configurable parameter representing the target number of indexed fragments per precursor.

**Derivation**: When P(miss) = ε and P(detect) = 1 - ε ≈ 1:

```
Δ_max = log((1-ε)/α) - log(ε/(1-α))
      ≈ log(1/α) - log(ε/(1-α))        [since ε << 1]
      = log((1-α)/(α × ε))
```

Setting `M × Δ_max = 255`:

```
ε = (1-α) / (α × exp(255/M))
```

**For α = 0.002, M = 7**: ε ≈ 8.2×10⁻¹⁴, Δ_max ≈ 36.4

M is a config parameter, NOT `length(intensities)`. This ensures the UInt8 scale is consistent across all precursors. A precursor with fewer than M fragments can still max out if all its fragments are highly detectable.

### Score interpretation

- A precursor with M highly-detectable fragments (all clamped at ε) gets total score ≈ 255
- A precursor with mixed fragments gets total score < 255, proportional to total evidence
- A precursor with mostly weak fragments gets a much lower total score
- The relative ordering within a precursor is preserved, AND absolute scale differences between precursors are preserved

## Configurable Parameters

All in `defaultBuildLibParams.json` under `library_params`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fragment_scoring` | `"rank_based"` | Scoring mode: `"rank_based"` (original) or `"likelihood_ratio"` (new) |
| `likelihood_ratio_M` | `7` | Target fragments per precursor (sets UInt8 scale via ε derivation) |
| `likelihood_ratio_N` | `10000` | Total precursor ion count (~10k for Astral MS2) |
| `likelihood_ratio_d` | `1` | Detection limit in ions (1 for Astral, 3-10 for Orbitrap) |
| `likelihood_ratio_alpha` | `0.002` | Random m/z coincidence probability (~0.002 for ±10 ppm) |
| `likelihood_ratio_ref_nce` | `27.0` | Reference NCE for evaluating Altimeter spline coefficients |

When `fragment_scoring` is `"rank_based"` (or absent), the existing `rank_to_score` path is used unchanged.

## Files to Modify

1. `src/Routines/BuildSpecLib/build/build_poin_lib.jl` — New scoring function + modified `getSimpleFrags` + spline evaluation in `buildPionLib`
2. `src/Routines/BuildSpecLib.jl` — Pass new params to `buildPionLib`
3. `assets/example_config/defaultBuildLibParams.json` — Add LR config fields
4. `src/Routines/SearchDIA/ParseInputs/paramDefaults.jl` — Add defaults for new params
5. `data/ecoli_test/ecoli_test_params.json` — Update min_index_search_score for new scale
6. `assets/example_config/defaultSearchParams.json` — Update min_index_search_score defaults

## Key Existing Functions to Reuse

- **`splevl(x, knots, c, k)`** in `src/utils/ML/libraryBSpline.jl:51` — Evaluates B-spline at point x given knots, coefficients, and degree k. Use to get fragment intensities from Altimeter spline coefficients at reference NCE.
- **`deserialize_from_jls(path)`** — Load serialized Julia objects (for knot vector)
- **`Binomial`, `cdf`** from `Distributions.jl` (already a dependency)

## Implementation

### Step 1: New scoring function in `build_pion_lib.jl`

```julia
using Distributions: Binomial, cdf

"""
    fragment_likelihood_weights_uint8(intensities; M, N, d, α) → Vector{UInt8}

Compute likelihood-ratio Δ_i weights for fragments and quantize to UInt8.

Uses global ε-clamped normalization: derives ε = (1-α)/(α·exp(255/M)) so that
M × Δ_max = 255. M is a config parameter (target fragments per precursor), not
the actual number of fragments. This preserves absolute score differences between
precursors.
"""
function fragment_likelihood_weights_uint8(
    intensities::AbstractVector{<:Real};
    M::Int, N::Int, d::Int, α::Float64
)::Vector{UInt8}
    n_frags = length(intensities)
    n_frags == 0 && return UInt8[]

    total = sum(Float64, intensities)
    if total ≤ 0
        return fill(UInt8(max(1, 255 ÷ M)), n_frags)
    end
    p = Float64.(intensities) ./ total

    # Derive ε so M × Δ_max = 255 (M is config param, not n_frags)
    ε = (1.0 - α) / (α * exp(255.0 / M))

    Δ = Vector{Float64}(undef, n_frags)
    for i in 1:n_frags
        p_miss = clamp(cdf(Binomial(N, p[i]), d - 1), ε, 1.0 - 1e-15)
        p_detect = 1.0 - p_miss
        w_plus = log(p_detect / α)
        w_minus = log(p_miss / (1.0 - α))
        Δ[i] = max(0.0, w_plus - w_minus)
    end

    return UInt8[clamp(round(Int, Δ[i]), 1, 255) for i in 1:n_frags]
end
```

### Step 2: Add intensity parameter to `getSimpleFrags`

Add five new positional parameters after `rank_to_score`:
```julia
    frag_bounds, rank_to_score,
    # New: likelihood ratio scoring (optional)
    frag_intensity::Union{Nothing, AbstractVector{<:Real}} = nothing,
    lr_M::Int = 7,
    lr_N::Int = 10000,
    lr_d::Int = 1,
    lr_alpha::Float64 = 0.002
```

When `frag_intensity !== nothing`, use two-phase per-precursor loop:
1. Collect filtered fragment indices + intensities into pre-allocated buffers
2. Call `fragment_likelihood_weights_uint8` on the collected intensities (passing M from config)
3. Create `SimpleFrag`s with LR weights instead of `rank_to_score[rank]`

When `frag_intensity === nothing`, keep the existing rank-based path unchanged.

### Step 3: Evaluate spline coefficients in `buildPionLib`

Use the same `splevl` call pattern as SearchDIA (`LibraryIon.jl:332`):
```julia
# SearchDIA: splevl(getNCE(intensity_type), getKnots(intensity_type), pf.intensity, getDegree(intensity_type))
# We replicate: splevl(nce, knots, coefficients, degree) where degree=3
```

**SplineCoefficientModel overload (~line 261)**:
1. Load knot vector: `spl_knots = deserialize_from_jls(joinpath(spec_lib_path, "spline_knots.jls"))`
2. Convert to Tuple (same as `loadSpectralLibrary.jl:105`): `knots = Tuple(spl_knots)`
3. Evaluate each fragment's spline at reference NCE using `splevl`:
   ```julia
   coefficients_col = fragments_table[:coefficients]
   n_frags = length(coefficients_col)
   ref_nce_f32 = Float32(ref_nce)
   evaluated_intensities = Vector{Float32}(undef, n_frags)
   for i in 1:n_frags
       evaluated_intensities[i] = max(0f0, splevl(ref_nce_f32, knots, coefficients_col[i], 3))
   end
   ```
4. Pass `evaluated_intensities` to `getSimpleFrags` as `frag_intensity`

**KoinaModelType overload (~line 65)**: Pass `fragments_table[:intensity]` directly (already has real predicted intensities).

### Step 4: Config parameters

In `assets/example_config/defaultBuildLibParams.json`, inside `library_params`:
```json
{
    "library_params": {
        "rank_to_score": [8, 4, 4, 2, 2, 1, 1],
        "fragment_scoring": "likelihood_ratio",
        "likelihood_ratio_M": 7,
        "likelihood_ratio_N": 10000,
        "likelihood_ratio_d": 1,
        "likelihood_ratio_alpha": 0.002,
        "likelihood_ratio_ref_nce": 27.0,
        ...
    }
}
```

When `fragment_scoring` is `"rank_based"` (or absent), use existing `rank_to_score`. When `"likelihood_ratio"`, use the new path.

### Step 5: Adjust `min_index_search_score` thresholds

With global ε-clamped normalization:
- Max possible per-precursor score: M × Δ_max = 255 (when all fragments are highly detectable)
- Typical per-precursor score: varies based on actual fragment detectability
- Current `rank_to_score` sum = 22, threshold = 3 ≈ 14% of max

Equivalent thresholds on the new 0-255 scale:
- `ParameterTuningSearch`: `min_index_search_score` 3 → 35
- `FirstPassSearch`: `min_index_search_score` 3 → 35
- Presearch: 22 → ~200 (slightly less than 255 for tolerance)
- Test config: adjust accordingly

**Note**: These thresholds only matter when using `likelihood_ratio` scoring. The `rank_based` path produces the same scores as before.

**Interpreting a threshold**: Since Δ_max = 255/M, a threshold of `k × (255/M)` means "at least k max-evidence fragment detections' worth of evidence." For M=7, Δ_max ≈ 36.4:
- threshold 36 ≈ 1 perfect fragment match
- threshold 73 ≈ 2 perfect fragment matches
- threshold 109 ≈ 3 perfect fragment matches

A weaker fragment (Δ < Δ_max) counts as a fraction of a max-evidence fragment. The format remains raw UInt8, same as current `min_index_search_score`.

### Step 6: Plumb config through BuildSpecLib.jl

In `src/Routines/BuildSpecLib.jl`, the `buildPionLib` calls need the new params:
- Read `fragment_scoring`, `likelihood_ratio_M/N/d/alpha/ref_nce` from parsed JSON config
- Pass through to `buildPionLib` call sites

## How the scoring works (reference)

From the spec, for each fragment i of a precursor:

```
w⁺ᵢ = log(P(detect i | correct precursor) / P(detect i | random match))
w⁻ᵢ = log(P(miss i | correct precursor) / P(miss i | random match))
Δᵢ  = w⁺ᵢ - w⁻ᵢ  (score increment when fragment i is detected)
```

Where:
- `P(detect i | correct) = 1 - CDF_Binomial(d-1; N, p'_i)`
- `P(detect i | random) = α`
- `P(miss i | correct) = CDF_Binomial(d-1; N, p'_i)`
- `P(miss i | random) = 1 - α`
- `p'_i = r_i / sum(r)` (renormalized predicted intensity)

Physical parameters:
- `M` = target number of indexed fragments per precursor (sets UInt8 scale)
- `N` = total precursor ion count (~10,000 for Astral MS2)
- `d` = detection limit in ions (1 for Astral, 3-10 for Orbitrap)
- `α` = random m/z coincidence probability (~0.002 for ±10 ppm at density 0.2 peaks/Da)

Quantization:
- `ε = (1-α) / (α × exp(255/M))` — derived so M × Δ_max = 255
- `P(miss)` is clamped at ε from below
- Δ values are rounded to nearest integer for UInt8 storage
- For α=0.002, M=7: ε ≈ 8.2×10⁻¹⁴, Δ_max ≈ 36.4, sum at max ≈ 255

## Verification

1. **Unit test**: Call `fragment_likelihood_weights_uint8` with known intensities and verify:
   - Dominant fragment gets the largest weight
   - All weights ≥ 1
   - Total sum ≤ 255 when n_frags ≤ M
   - Edge cases: equal intensities, single fragment, very skewed distribution
   - Verify that different intensity distributions produce different total scores

2. **Build a library**: Run `BuildSpecLib` with a test FASTA to confirm the pipeline completes and produces valid `.poin` output with likelihood-ratio scores.

3. **Integration test**: Run `SearchDIA("./data/ecoli_test/ecoli_test_params.json")` with adjusted thresholds to confirm the search pipeline works end-to-end.

4. **Comparison**: Build the same library with `rank_based` and `likelihood_ratio` scoring, run the same search, compare PSM counts and FDR curves.
