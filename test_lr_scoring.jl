"""
Exploration script for likelihood-ratio fragment scoring.

Tests different parameter combinations (N, d, α) and intensity distributions
to understand score behavior and evaluate normalization strategies.

Run with: julia test_lr_scoring.jl
Or interactively: include("test_lr_scoring.jl")
"""

using Distributions: Binomial, cdf

#=============================================================================
  RAW SCORING FUNCTIONS
=============================================================================#

"""
Raw Δ_i values (no clamping, no normalization).
Returns the full likelihood-ratio score increment for each fragment.
Negative values mean detecting fragment i is evidence AGAINST the precursor.
"""
function raw_fragment_weights(r::Vector{Float64}; N::Int, d::Int, α::Float64)
    p = r ./ sum(r)
    ε = 1e-300  # just to avoid log(0)
    p_miss = [clamp(cdf(Binomial(N, p[i]), d - 1), ε, 1.0 - ε) for i in eachindex(p)]
    p_detect = 1.0 .- p_miss
    w_plus = log.(p_detect ./ α)
    w_minus = log.(p_miss ./ (1.0 - α))
    Δ = w_plus .- w_minus
    baseline = sum(w_minus)
    return (; Δ, baseline, w_plus, w_minus, p_miss, p_detect, p)
end

"""
Current Pioneer implementation: ε-clamped global normalization.
ε = (1-α)/(α·exp(255/M)) so M × Δ_max = 255.
Clamps P(miss) at ε, then clamps Δ to [0, ∞), then rounds to UInt8 [1, 255].
"""
function pioneer_current_weights(r::Vector{Float64}; M::Int, N::Int, d::Int, α::Float64)
    n = length(r)
    total = sum(r)
    total ≤ 0 && return fill(max(1, 255 ÷ M), n)
    p = r ./ total

    ε = (1.0 - α) / (α * exp(255.0 / M))

    Δ = Vector{Float64}(undef, n)
    for i in 1:n
        p_miss = clamp(cdf(Binomial(N, p[i]), d - 1), ε, 1.0 - 1e-15)
        p_detect = 1.0 - p_miss
        Δ[i] = max(0.0, log(p_detect / α) - log(p_miss / (1.0 - α)))
    end
    return (;
        weights_uint8 = UInt8[clamp(round(Int, Δ[i]), 1, 255) for i in 1:n],
        Δ_raw = Δ,
        ε = ε,
        Δ_max = 255.0 / M,
        p = p
    )
end

"""
Per-precursor normalization: compute raw Δ, clamp negatives to 0,
then scale so the sum = 255. Every precursor gets the same total budget.
"""
function perprecursor_normalized_weights(r::Vector{Float64}; N::Int, d::Int, α::Float64)
    n = length(r)
    total = sum(r)
    total ≤ 0 && return fill(UInt8(255 ÷ n), n)
    p = r ./ total

    ε_floor = 1e-300
    Δ_raw = Vector{Float64}(undef, n)
    for i in 1:n
        p_miss = clamp(cdf(Binomial(N, p[i]), d - 1), ε_floor, 1.0 - ε_floor)
        p_detect = 1.0 - p_miss
        Δ_raw[i] = max(0.0, log(p_detect / α) - log(p_miss / (1.0 - α)))
    end

    Δ_sum = sum(Δ_raw)
    if Δ_sum ≤ 0
        return (;
            weights_uint8 = fill(UInt8(max(1, 255 ÷ n)), n),
            Δ_raw = Δ_raw,
            fractions = fill(1.0/n, n),
            p = p
        )
    end

    fractions = Δ_raw ./ Δ_sum
    scaled = fractions .* 255.0
    weights_uint8 = UInt8[clamp(round(Int, scaled[i]), 1, 255) for i in 1:n]

    return (; weights_uint8, Δ_raw, fractions, p)
end

#=============================================================================
  DISPLAY HELPERS
=============================================================================#

function print_header(title)
    println()
    println("=" ^ 80)
    println("  ", title)
    println("=" ^ 80)
end

function print_comparison(label, r; N, d, α, M=7)
    println("\n--- $label ---")
    println("  r = $r")
    println("  N=$N, d=$d, α=$α, M=$M")

    # Raw weights
    raw = raw_fragment_weights(r; N, d, α)
    println("\n  Raw Δ values (unclamped):")
    for i in eachindex(r)
        p_i = round(raw.p[i], digits=4)
        pmiss = round(raw.p_miss[i], sigdigits=3)
        pdet = round(raw.p_detect[i], sigdigits=4)
        delta = round(raw.Δ[i], digits=2)
        println("    frag $i: p=$p_i  P(miss)=$pmiss  P(det)=$pdet  Δ=$delta")
    end
    println("  Raw sum(Δ): ", round(sum(raw.Δ), digits=2))
    println("  Baseline (sum w⁻): ", round(raw.baseline, digits=2))

    # Pioneer current (ε-clamped global)
    cur = pioneer_current_weights(r; M, N, d, α)
    println("\n  Pioneer current (ε-clamped global, M=$M):")
    println("    ε = ", round(cur.ε, sigdigits=3))
    println("    Δ_max = ", round(cur.Δ_max, digits=1))
    println("    UInt8 weights: ", Int.(cur.weights_uint8))
    println("    sum = ", sum(Int, cur.weights_uint8))

    # Per-precursor normalized
    norm = perprecursor_normalized_weights(r; N, d, α)
    println("\n  Per-precursor normalized (sum→255):")
    println("    UInt8 weights: ", Int.(norm.weights_uint8))
    println("    sum = ", sum(Int, norm.weights_uint8))
    println("    fractions: ", round.(norm.fractions, digits=3))
end

#=============================================================================
  TEST CASES
=============================================================================#

function run_all_tests()

    # ---- Typical fragment distributions ----

    print_header("1. INTENSITY DISTRIBUTION COMPARISON (N=100, d=1, α=0.002)")

    dists = [
        ("Typical (decreasing)",   [100.0, 60.0, 35.0, 20.0, 10.0, 5.0, 2.0]),
        ("Skewed (one dominant)",   [100.0, 10.0, 10.0, 5.0, 1.0, 1.0, 0.5]),
        ("Even (all equal)",        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
        ("Two strong, rest weak",   [50.0, 50.0, 2.0, 1.0, 0.5, 0.2, 0.1]),
        ("One fragment only",       [1.0, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]),
        ("Realistic Altimeter",     [0.35, 0.22, 0.15, 0.12, 0.08, 0.05, 0.03]),
        ("3 good + 4 noise",        [0.40, 0.30, 0.20, 0.03, 0.03, 0.02, 0.02]),
    ]

    for (label, r) in dists
        print_comparison(label, r; N=100, d=1, α=0.002)
    end

    # ---- Effect of d (detection limit) ----

    print_header("2. EFFECT OF DETECTION LIMIT d (N=100, α=0.002)")

    r_typical = [100.0, 60.0, 35.0, 20.0, 10.0, 5.0, 2.0]
    for d in [1, 2, 3, 5, 10]
        print_comparison("d=$d", r_typical; N=100, d=d, α=0.002)
    end

    # ---- Effect of N ----

    print_header("3. EFFECT OF N (d=1, α=0.002)")

    for N in [50, 100, 200, 500, 1000]
        print_comparison("N=$N", r_typical; N=N, d=1, α=0.002)
    end

    # ---- Effect of N with d=5 ----

    print_header("4. EFFECT OF N WITH d=5 (α=0.002)")

    for N in [50, 100, 200, 500, 1000]
        print_comparison("N=$N, d=5", r_typical; N=N, d=5, α=0.002)
    end

    # ---- Effect of α ----

    print_header("5. EFFECT OF α (N=100, d=1)")

    for α in [0.001, 0.002, 0.005, 0.01, 0.05]
        print_comparison("α=$α", r_typical; N=100, d=1, α=α)
    end

    # ---- Different fragment counts ----

    print_header("6. DIFFERENT FRAGMENT COUNTS (N=100, d=1, α=0.002)")

    frag_counts = [
        ("3 fragments", [0.50, 0.30, 0.20]),
        ("5 fragments", [0.35, 0.25, 0.20, 0.12, 0.08]),
        ("7 fragments", [0.30, 0.22, 0.15, 0.12, 0.10, 0.06, 0.05]),
        ("10 fragments", [0.25, 0.15, 0.12, 0.10, 0.08, 0.07, 0.06, 0.06, 0.06, 0.05]),
    ]

    for (label, r) in frag_counts
        print_comparison(label, r; N=100, d=1, α=0.002)
    end

    # ---- Threshold analysis ----

    print_header("7. THRESHOLD ANALYSIS: WHAT HAPPENS IF ONLY TOP-K FRAGMENTS MATCH?")

    r_test = [100.0, 60.0, 35.0, 20.0, 10.0, 5.0, 2.0]
    println("\n  r = $r_test")
    println("\n  Scenario: precursor is present, k of 7 fragments match (strongest first)")

    for (N, d) in [(100, 1), (100, 5), (200, 1), (200, 5)]
        println("\n  N=$N, d=$d, α=0.002:")

        cur = pioneer_current_weights(r_test; M=7, N=N, d=d, α=0.002)
        norm = perprecursor_normalized_weights(r_test; N=N, d=d, α=0.002)

        println("    Pioneer current weights: ", Int.(cur.weights_uint8))
        println("    Per-prec normalized:     ", Int.(norm.weights_uint8))

        for k in 1:7
            cur_score = sum(Int, cur.weights_uint8[1:k])
            norm_score = sum(Int, norm.weights_uint8[1:k])
            println("    Top-$k match: current=$cur_score, normalized=$norm_score")
        end
    end

    # ---- Cross-precursor comparison ----

    print_header("8. CROSS-PRECURSOR COMPARISON: SAME # MATCHES, DIFFERENT DISTRIBUTIONS")

    println("\n  All 7 fragments match. Which precursor scores higher?")

    for (N, d) in [(100, 1), (100, 5), (200, 5)]
        println("\n  N=$N, d=$d, α=0.002:")

        distributions = [
            ("Even",     [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
            ("Typical",  [100.0, 60.0, 35.0, 20.0, 10.0, 5.0, 2.0]),
            ("Skewed",   [100.0, 10.0, 5.0, 2.0, 1.0, 0.5, 0.1]),
            ("2 strong", [50.0, 50.0, 1.0, 0.5, 0.2, 0.1, 0.05]),
        ]

        for (label, r) in distributions
            cur = pioneer_current_weights(r; M=7, N=N, d=d, α=0.002)
            norm = perprecursor_normalized_weights(r; N=N, d=d, α=0.002)
            raw = raw_fragment_weights(r; N=N, d=d, α=0.002)

            println("    $label: current_sum=$(sum(Int, cur.weights_uint8)), norm_sum=$(sum(Int, norm.weights_uint8)), raw_sum=$(round(sum(raw.Δ), digits=1))")
        end
    end

    # ---- The key question: does per-precursor normalization lose important info? ----

    print_header("9. PER-PRECURSOR NORMALIZATION: WHAT DO WE LOSE?")

    println("""
    With per-precursor normalization (sum→255):
    ✓ Every precursor has the same max possible score (255)
    ✓ Threshold is simple: e.g., 128 = "half the fragments' evidence"
    ✓ Score reflects which fragments matched, weighted by importance
    ✗ Can't distinguish "strong evidence precursor" from "weak evidence precursor"

    But: this is a PRE-SCREENING step. The real scoring happens later with
    spectral contrast, matched ratio, RT alignment, etc. So all we need here
    is "did enough of the expected fragments match, weighted by importance?"

    The per-precursor normalization makes the threshold universal and interpretable.
    """)

    # ---- Recommended parameter exploration ----

    print_header("10. PARAMETER SWEEP: FINDING GOOD DEFAULTS")

    r_test = [0.35, 0.22, 0.15, 0.12, 0.08, 0.05, 0.03]
    println("\n  Realistic Altimeter distribution: $r_test")
    println("\n  Goal: meaningful differentiation between fragments,")
    println("  weak fragments should still get nonzero weight,")
    println("  no negative Δ values (or very few).\n")

    println("  " * rpad("N", 6) * rpad("d", 4) * rpad("α", 8) * "Raw Δ values" * " " ^ 30 * "# negative")
    println("  " * "-"^100)

    for N in [50, 100, 200, 500]
        for d in [1, 3, 5]
            for α in [0.002, 0.01]
                raw = raw_fragment_weights(r_test; N=N, d=d, α=α)
                n_neg = count(x -> x < 0, raw.Δ)
                deltas = join([lpad(round(raw.Δ[i], digits=1), 7) for i in eachindex(raw.Δ)], " ")
                println("  ", rpad(N, 6), rpad(d, 4), rpad(α, 8), deltas, "   ", n_neg, " neg")
            end
        end
    end
end

# Run everything
run_all_tests()
