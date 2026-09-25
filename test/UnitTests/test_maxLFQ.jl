# Unit tests for MaxLFQ algorithm (src/utils/maxLFQ.jl).
#
# Tests the core linear algebra functions (getS, getB, getA) and the
# end-to-end getProtAbundance protein quantification.
#
# Run standalone: julia --project=. test/UnitTests/test_maxLFQ.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
    using Statistics
end

using Random

@testset "MaxLFQ" begin

@testset "MaxLFQ matches independent dense equations" begin
    function dense_reference(X)
        n = size(X, 2)
        A = zeros(n + 1, n + 1)
        b = zeros(n + 1)
        for i in 1:n-1, j in i+1:n
            ratios = [Float64(X[k, j]) - Float64(X[k, i]) for k in axes(X, 1)
                      if !ismissing(X[k, i]) && !ismissing(X[k, j])]
            isempty(ratios) && continue
            A[i, i] += 2
            A[j, j] += 2
            A[i, j] -= 2
            A[j, i] -= 2
            ratio = 2 * median(ratios)
            b[i] -= ratio
            b[j] += ratio
        end
        A[1:n, end] .= 1
        A[end, 1:n] .= 1
        profile = (A \ b)[1:n]
        observed = collect(skipmissing(vec(X)))
        scale = maximum(observed) + log2(sum(exp2, observed .- maximum(observed))) -
                maximum(profile) - log2(sum(exp2, profile .- maximum(profile)))
        return profile .+ scale
    end

    rng = MersenneTwister(123)
    for T in (Float32, Float64), n in (2, 3, 17), p in (1, 2, 7, 8)
        X = Matrix{Union{Missing, T}}(20 .+ 3 .* randn(rng, T, p, n))
        # Keep one shared precursor, but vary the other pairwise overlap counts.
        for j in 1:n, i in 2:p
            rand(rng) < 0.4 && (X[i, j] = missing)
        end
        @test Pioneer.solve_maxlfq_component(X) ≈ dense_reference(X) rtol=1e-6
    end

    # Every pair overlaps, although no precursor is shared by all three runs.
    X = Union{Missing, Float64}[10 12 missing; missing 15 18; 8 missing 9]
    @test Pioneer.solve_maxlfq_component(X) ≈ dense_reference(X) rtol=1e-6
    @test Pioneer.solve_maxlfq_component(@view X[:, [3, 1, 2]]) ≈
          dense_reference(X)[[3, 1, 2]] rtol=1e-6

    # Connected through a chain, with incomplete pairwise overlap.
    X = Union{Missing, Float64}[10 12 missing missing; missing 15 18 missing;
                               missing missing 8 9]
    @test Pioneer.solve_maxlfq_component(X) ≈ dense_reference(X) rtol=1e-6
    @test all(ismissing, Pioneer.solve_maxlfq_component(
        Union{Missing, Float64}[10 missing; missing 12]))
    @test all(ismissing, Pioneer.solve_maxlfq_component(
        fill!(Matrix{Union{Missing, Float64}}(undef, 2, 3), missing)))
    @test all(ismissing, Pioneer.solve_maxlfq_component(Union{Missing, Float64}[10;;]))
    @test isempty(Pioneer.solve_maxlfq_component(Matrix{Union{Missing, Float64}}(undef, 2, 0)))

    # Pairwise ratio calculation reuses its buffer across all run pairs.
    # Inputs and compilation are excluded from the allocation measurement.
    X = Matrix{Union{Missing, Float64}}(randn(rng, 8, 256))
    Pioneer.getB(X)
    @test (@allocated Pioneer.getB(X)) < 200_000
end

# ═══════════════════════════════════════════════════════════════════════════════
# getS — build peptide × experiment intensity matrix
# ═══════════════════════════════════════════════════════════════════════════════

@testset "getS" begin
    # 2 peptides, 3 experiments, 4 observations (one peptide missing from one experiment)
    peptides  = UInt32[1, 1, 2, 2]
    pep_dict  = Dict(UInt32(1) => 1, UInt32(2) => 2)
    exps      = UInt16[1, 2, 1, 3]
    exp_dict  = Dict(UInt16(1) => 1, UInt16(2) => 2, UInt16(3) => 3)
    abundance = Union{Float32, Missing}[10.0f0, 20.0f0, 5.0f0, 15.0f0]

    S = Pioneer.getS(peptides, pep_dict, exps, exp_dict, abundance, 2, 3)

    @testset "dimensions" begin
        @test size(S) == (2, 3)
    end

    @testset "filled values" begin
        @test S[1, 1] == 10.0f0  # peptide 1, experiment 1
        @test S[1, 2] == 20.0f0  # peptide 1, experiment 2
        @test S[2, 1] == 5.0f0   # peptide 2, experiment 1
        @test S[2, 3] == 15.0f0  # peptide 2, experiment 3
    end

    @testset "missing entries" begin
        @test ismissing(S[1, 3])  # peptide 1 not seen in experiment 3
        @test ismissing(S[2, 2])  # peptide 2 not seen in experiment 2
    end

    @testset "zero abundance treated as missing" begin
        abundance_zero = Union{Float32, Missing}[0.0f0, 20.0f0, 5.0f0, 15.0f0]
        S2 = Pioneer.getS(peptides, pep_dict, exps, exp_dict, abundance_zero, 2, 3)
        @test ismissing(S2[1, 1])
    end

    @testset "missing abundance stays missing" begin
        abundance_miss = Union{Float32, Missing}[missing, 20.0f0, 5.0f0, 15.0f0]
        S3 = Pioneer.getS(peptides, pep_dict, exps, exp_dict, abundance_miss, 2, 3)
        @test ismissing(S3[1, 1])
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# getA — MaxLFQ normal-equation LHS (AtA, n_runs × n_runs).
#
# Old API: `getA(N::Int)` returned an (N+1)×(N+1) matrix encoding a
# constraint row. The algorithm was rewritten to assemble AtA directly
# from the log2 intensity matrix; the constraint row no longer exists.
# ═══════════════════════════════════════════════════════════════════════════════

@testset "getA" begin
    @testset "fully connected 3-run / 2-peptide matrix" begin
        # Every run pair shares both peptides → fully connected graph.
        X = Union{Float32, Missing}[
            log2(10.0f0)  log2(20.0f0)  log2(40.0f0);
            log2(5.0f0)   log2(10.0f0)  log2(20.0f0)
        ]
        A = Pioneer.getA(X)
        @test size(A) == (3, 3)
        # Diagonal = degree of each node (connections to other runs):
        # each run is connected to the other 2 → degree 2.
        @test A[1,1] == 2.0
        @test A[2,2] == 2.0
        @test A[3,3] == 2.0
        # Off-diagonal = -1 for every connected run pair, symmetric.
        @test A[1,2] == -1.0
        @test A[2,1] == -1.0
        @test A[1,3] == -1.0
        @test A[2,3] == -1.0
    end

    @testset "disconnected pair" begin
        # Run 1 has only peptide 1, run 2 has only peptide 2 — no shared
        # peptides → no edge → disconnected.
        X = Union{Float32, Missing}[
            10.0f0   missing;
            missing  20.0f0
        ]
        A = Pioneer.getA(X)
        @test size(A) == (2, 2)
        @test A[1,1] == 0.0  # no edges → degree 0
        @test A[2,2] == 0.0
        @test A[1,2] == 0.0  # no edge
        @test A[2,1] == 0.0
    end

    @testset "single run" begin
        X = Union{Float32, Missing}[10.0f0;;]  # 1×1
        A = Pioneer.getA(X)
        @test size(A) == (1, 1)
        @test A[1,1] == 0.0
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# getB — MaxLFQ normal-equation RHS (Atb, length n_runs) + valid_pairs count.
#
# Old API: `getB(S, N, n_peptides)` returned a length-(N+1) vector. The
# new API takes only the log2 intensity matrix and returns
# `(Atb, valid_pairs::Int)`.
# ═══════════════════════════════════════════════════════════════════════════════

@testset "getB" begin
    @testset "uniform 2× ratios" begin
        # 2 peptides, 3 experiments. Each successive run is 2× the prior →
        # log2 ratios 1.0 between adjacent runs. Pre-log the matrix because
        # solve_maxlfq_component now expects log2 inputs.
        X = Union{Float32, Missing}[
            log2(10.0f0)  log2(20.0f0)  log2(40.0f0);
            log2(5.0f0)   log2(10.0f0)  log2(20.0f0)
        ]
        Atb, valid_pairs = Pioneer.getB(X)
        @test length(Atb) == 3
        @test valid_pairs == 3  # all 3 pairs share both peptides
        # Run 1 sits below run 3 in intensity → Atb[1] < Atb[3].
        @test Atb[1] < Atb[2] < Atb[3]
    end

    @testset "single run" begin
        X = Union{Float32, Missing}[10.0f0; 20.0f0;;]  # 2×1
        Atb, valid_pairs = Pioneer.getB(X)
        @test length(Atb) == 1
        @test valid_pairs == 0  # no pairs
        @test Atb[1] == 0.0
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# solve_maxlfq_component — relative ratios and cumulative-intensity scale
# ═══════════════════════════════════════════════════════════════════════════════

@testset "solve_maxlfq_component" begin
    X = Union{Float64, Missing}[
        log2(10.0)  log2(20.0)  log2(40.0);
        log2(5.0)   log2(10.0)  log2(20.0)
    ]

    estimates = Pioneer.solve_maxlfq_component(X)
    linear_estimates = exp2.(Float64.(estimates))
    observed_cumulative_intensity = sum(exp2, skipmissing(vec(X)))

    @testset "preserves precursor ratios" begin
        @test isapprox(estimates[2] - estimates[1], 1.0; atol = 1.0e-6)
        @test isapprox(estimates[3] - estimates[2], 1.0; atol = 1.0e-6)
    end

    @testset "preserves cumulative precursor intensity" begin
        @test isapprox(
            sum(linear_estimates),
            observed_cumulative_intensity;
            rtol = 1.0e-6
        )
        @test all(isapprox.(linear_estimates, [15.0, 30.0, 60.0]; rtol = 1.0e-6))
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# getProtAbundance — end-to-end protein quantification
# ═══════════════════════════════════════════════════════════════════════════════

@testset "getProtAbundance" begin
    # Setup: 2 peptides seen in 3 experiments with 2× fold changes
    peptides  = UInt32[1, 1, 1, 2, 2, 2]
    exps      = UInt16[1, 2, 3, 1, 2, 3]
    use_quant = Bool[true, true, true, true, true, true]
    abundance = Union{Float32, Missing}[10.0f0, 20.0f0, 40.0f0, 5.0f0, 10.0f0, 20.0f0]
    qvals     = Union{Float32, Missing}[0.01f0, 0.01f0, 0.01f0, 0.02f0, 0.02f0, 0.02f0]
    gqvals    = Union{Float32, Missing}[0.005f0, 0.005f0, 0.005f0, 0.01f0, 0.01f0, 0.01f0]
    peps      = Union{Float32, Missing}[0.05f0, 0.05f0, 0.05f0, 0.1f0, 0.1f0, 0.1f0]
    scores    = Union{Float32, Missing}[5.0f0, 5.0f0, 5.0f0, 3.0f0, 3.0f0, 3.0f0]
    gscores   = Union{Float32, Missing}[4.0f0, 4.0f0, 4.0f0, 2.0f0, 2.0f0, 2.0f0]

    n_experiments = 3

    # Pre-allocate output vectors
    target_out       = Vector{Union{Missing, Bool}}(missing, n_experiments)
    entrap_out       = Vector{Union{Missing, UInt8}}(missing, n_experiments)
    species_out      = Vector{Union{Missing, String}}(missing, n_experiments)
    protein_out      = Vector{Union{Missing, String}}(missing, n_experiments)
    peptides_out     = Vector{Union{Missing, Vector{Union{Missing, UInt32}}}}(missing, n_experiments)
    experiments_out  = Vector{Union{Missing, UInt32}}(missing, n_experiments)
    log2_abund_out   = Vector{Union{Missing, Float32}}(missing, n_experiments)
    gqval_out        = Vector{Union{Missing, Float32}}(missing, n_experiments)
    qval_out         = Vector{Union{Missing, Float32}}(missing, n_experiments)
    pep_out          = Vector{Union{Missing, Float32}}(missing, n_experiments)
    score_out        = Vector{Union{Missing, Float32}}(missing, n_experiments)
    gscore_out       = Vector{Union{Missing, Float32}}(missing, n_experiments)
    peak_area_out    = Vector{Union{Missing, Float32}}(missing, n_experiments)
    n_quant_out      = Vector{Union{Missing, UInt32}}(missing, n_experiments)

    Pioneer.getProtAbundance(
        "ProteinA", 1, true, UInt8(0), "HUMAN",
        peptides, exps, use_quant, abundance,
        gqvals, qvals, peps, scores, gscores,
        target_out, entrap_out, species_out, protein_out,
        peptides_out, experiments_out, log2_abund_out,
        gqval_out, qval_out, pep_out, score_out, gscore_out,
        peak_area_out, n_quant_out,
    )

    @testset "output filled for all experiments" begin
        @test all(!ismissing, protein_out)
        @test all(!ismissing, experiments_out)
        @test all(!ismissing, log2_abund_out)
    end

    @testset "protein name correct" begin
        @test all(x -> x == "ProteinA", protein_out)
    end

    @testset "metadata propagated" begin
        @test all(x -> x == true, target_out)
        @test all(x -> x == UInt8(0), entrap_out)
        @test all(x -> x == "HUMAN", species_out)
    end

    @testset "log2 abundances have correct fold changes" begin
        # With uniform 2× ratios, differences should be ~1.0
        vals = collect(skipmissing(log2_abund_out))
        diffs = diff(vals)
        @test all(d -> isapprox(d, 1.0, atol=0.1), diffs)
    end

    @testset "protein profile preserves cumulative precursor intensity" begin
        protein_cumulative_intensity = sum(exp2, skipmissing(log2_abund_out))
        @test isapprox(protein_cumulative_intensity, sum(abundance); rtol = 1.0e-6)
    end

    @testset "peptide lists populated" begin
        for p in peptides_out
            @test !ismissing(p)
            # Both peptides seen in all experiments, so 2 non-missing entries
            @test count(!ismissing, p) == 2
        end
    end

    @testset "quality metrics propagated" begin
        @test all(!ismissing, qval_out)
        @test all(!ismissing, pep_out)
        @test all(!ismissing, score_out)
    end

    @testset "quantified precursor counts" begin
        # Every precursor has a positive abundance here, so the quantified count
        # matches the peptide list in every experiment.
        @test all(!ismissing, n_quant_out)
        @test all(x -> x == UInt32(2), n_quant_out)
    end
end

@testset "getProtAbundance accepts non-nullable normalized quant values" begin
    peptides  = UInt32[1, 1, 1, 2, 2, 2]
    exps      = UInt16[1, 2, 3, 1, 2, 3]
    use_quant = Bool[true, true, true, true, true, true]
    abundance = Float64[10.0, 20.0, 40.0, 5.0, 10.0, 20.0]
    qvals     = Float32[0.01f0, 0.01f0, 0.01f0, 0.02f0, 0.02f0, 0.02f0]
    gqvals    = Float32[0.005f0, 0.005f0, 0.005f0, 0.01f0, 0.01f0, 0.01f0]
    peps      = Float32[0.05f0, 0.05f0, 0.05f0, 0.1f0, 0.1f0, 0.1f0]
    scores    = Float32[5.0f0, 5.0f0, 5.0f0, 3.0f0, 3.0f0, 3.0f0]
    gscores   = Float32[4.0f0, 4.0f0, 4.0f0, 2.0f0, 2.0f0, 2.0f0]

    n_experiments = 3
    target_out       = Vector{Union{Missing, Bool}}(missing, n_experiments)
    entrap_out       = Vector{Union{Missing, UInt8}}(missing, n_experiments)
    species_out      = Vector{Union{Missing, String}}(missing, n_experiments)
    protein_out      = Vector{Union{Missing, String}}(missing, n_experiments)
    peptides_out     = Vector{Union{Missing, Vector{Union{Missing, UInt32}}}}(missing, n_experiments)
    experiments_out  = Vector{Union{Missing, UInt32}}(missing, n_experiments)
    log2_abund_out   = Vector{Union{Missing, Float32}}(missing, n_experiments)
    gqval_out        = Vector{Union{Missing, Float32}}(missing, n_experiments)
    qval_out         = Vector{Union{Missing, Float32}}(missing, n_experiments)
    pep_out          = Vector{Union{Missing, Float32}}(missing, n_experiments)
    score_out        = Vector{Union{Missing, Float32}}(missing, n_experiments)
    gscore_out       = Vector{Union{Missing, Float32}}(missing, n_experiments)
    peak_area_out    = Vector{Union{Missing, Float32}}(missing, n_experiments)
    n_quant_out      = Vector{Union{Missing, UInt32}}(missing, n_experiments)

    Pioneer.getProtAbundance(
        "ProteinA", 1, true, UInt8(0), "HUMAN",
        peptides, exps, use_quant, abundance,
        gqvals, qvals, peps, scores, gscores,
        target_out, entrap_out, species_out, protein_out,
        peptides_out, experiments_out, log2_abund_out,
        gqval_out, qval_out, pep_out, score_out, gscore_out,
        peak_area_out, n_quant_out,
    )

    @test all(!ismissing, log2_abund_out)
    @test all(!ismissing, peak_area_out)
    @test all(x -> x == "ProteinA", protein_out)
end

# ═══════════════════════════════════════════════════════════════════════════════
# getProtAbundance — missing data handling
# ═══════════════════════════════════════════════════════════════════════════════

@testset "getProtAbundance missing data" begin
    # Peptide 2 is missing from experiment 2
    peptides  = UInt32[1, 1, 2, 2]
    exps      = UInt16[1, 2, 1, 3]
    use_quant = Bool[true, true, true, true]
    abundance = Union{Float32, Missing}[10.0f0, 20.0f0, 5.0f0, 10.0f0]
    qvals     = Union{Float32, Missing}[0.01f0, 0.01f0, 0.02f0, 0.02f0]
    gqvals    = Union{Float32, Missing}[0.005f0, 0.005f0, 0.01f0, 0.01f0]
    peps      = Union{Float32, Missing}[0.05f0, 0.05f0, 0.1f0, 0.1f0]
    scores    = Union{Float32, Missing}[5.0f0, 5.0f0, 3.0f0, 3.0f0]
    gscores   = Union{Float32, Missing}[4.0f0, 4.0f0, 2.0f0, 2.0f0]

    n_experiments = 3

    target_out       = Vector{Union{Missing, Bool}}(missing, n_experiments)
    entrap_out       = Vector{Union{Missing, UInt8}}(missing, n_experiments)
    species_out      = Vector{Union{Missing, String}}(missing, n_experiments)
    protein_out      = Vector{Union{Missing, String}}(missing, n_experiments)
    peptides_out     = Vector{Union{Missing, Vector{Union{Missing, UInt32}}}}(missing, n_experiments)
    experiments_out  = Vector{Union{Missing, UInt32}}(missing, n_experiments)
    log2_abund_out   = Vector{Union{Missing, Float32}}(missing, n_experiments)
    gqval_out        = Vector{Union{Missing, Float32}}(missing, n_experiments)
    qval_out         = Vector{Union{Missing, Float32}}(missing, n_experiments)
    pep_out          = Vector{Union{Missing, Float32}}(missing, n_experiments)
    score_out        = Vector{Union{Missing, Float32}}(missing, n_experiments)
    gscore_out       = Vector{Union{Missing, Float32}}(missing, n_experiments)
    peak_area_out    = Vector{Union{Missing, Float32}}(missing, n_experiments)
    n_quant_out      = Vector{Union{Missing, UInt32}}(missing, n_experiments)

    Pioneer.getProtAbundance(
        "ProteinB", 1, true, UInt8(0), "HUMAN",
        peptides, exps, use_quant, abundance,
        gqvals, qvals, peps, scores, gscores,
        target_out, entrap_out, species_out, protein_out,
        peptides_out, experiments_out, log2_abund_out,
        gqval_out, qval_out, pep_out, score_out, gscore_out,
        peak_area_out, n_quant_out,
    )

    @testset "still produces finite abundances" begin
        vals = collect(skipmissing(log2_abund_out))
        @test all(isfinite, vals)
        @test length(vals) == 3
        @test isapprox(sum(exp2, vals), sum(abundance); rtol = 1.0e-6)
    end

    @testset "peptide list reflects missing data" begin
        # Find experiment 2's peptide list — should only have peptide 1
        for i in 1:n_experiments
            if experiments_out[i] == UInt32(2)
                @test count(!ismissing, peptides_out[i]) == 1
            end
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# getProtAbundance — the scoring set vs. the quantified set
# ═══════════════════════════════════════════════════════════════════════════════

@testset "getProtAbundance reports the scoring set, not the quantified set" begin
    # A precursor observed in a run but whose area was zeroed -- integration found
    # no peak, or QUANT_MIN_AREA_SURVIVING_RATIO withheld it -- still scored the
    # protein group; the roll-up mask has no peak_area term. `getS` drops it,
    # because MaxLFQ must not ingest a zero, so deriving the support columns from
    # `S` reported protein groups with no peptides beside a passing q-value.
    # `getPresence` is the independent occupancy matrix that fixes that, and
    # `n_precursors_quantified` keeps the quantified subset visible alongside it.
    peptides  = UInt32[1, 1, 1, 2, 2, 2]
    exps      = UInt16[1, 2, 3, 1, 2, 3]
    use_quant = Bool[true, true, true, true, true, true]
    # Peptide 2 is present in experiment 2 but carries a zeroed area.
    abundance = Union{Float32, Missing}[10.0f0, 20.0f0, 40.0f0, 5.0f0, 0.0f0, 20.0f0]
    qvals     = Union{Float32, Missing}[0.01f0, 0.01f0, 0.01f0, 0.02f0, 0.02f0, 0.02f0]
    gqvals    = Union{Float32, Missing}[0.005f0, 0.005f0, 0.005f0, 0.01f0, 0.01f0, 0.01f0]
    peps      = Union{Float32, Missing}[0.05f0, 0.05f0, 0.05f0, 0.1f0, 0.1f0, 0.1f0]
    scores    = Union{Float32, Missing}[5.0f0, 5.0f0, 5.0f0, 3.0f0, 3.0f0, 3.0f0]
    gscores   = Union{Float32, Missing}[4.0f0, 4.0f0, 4.0f0, 2.0f0, 2.0f0, 2.0f0]

    n_experiments = 3

    target_out       = Vector{Union{Missing, Bool}}(missing, n_experiments)
    entrap_out       = Vector{Union{Missing, UInt8}}(missing, n_experiments)
    species_out      = Vector{Union{Missing, String}}(missing, n_experiments)
    protein_out      = Vector{Union{Missing, String}}(missing, n_experiments)
    peptides_out     = Vector{Union{Missing, Vector{Union{Missing, UInt32}}}}(missing, n_experiments)
    experiments_out  = Vector{Union{Missing, UInt32}}(missing, n_experiments)
    log2_abund_out   = Vector{Union{Missing, Float32}}(missing, n_experiments)
    gqval_out        = Vector{Union{Missing, Float32}}(missing, n_experiments)
    qval_out         = Vector{Union{Missing, Float32}}(missing, n_experiments)
    pep_out          = Vector{Union{Missing, Float32}}(missing, n_experiments)
    score_out        = Vector{Union{Missing, Float32}}(missing, n_experiments)
    gscore_out       = Vector{Union{Missing, Float32}}(missing, n_experiments)
    peak_area_out    = Vector{Union{Missing, Float32}}(missing, n_experiments)
    n_quant_out      = Vector{Union{Missing, UInt32}}(missing, n_experiments)

    Pioneer.getProtAbundance(
        "ProteinC", 1, true, UInt8(0), "HUMAN",
        peptides, exps, use_quant, abundance,
        gqvals, qvals, peps, scores, gscores,
        target_out, entrap_out, species_out, protein_out,
        peptides_out, experiments_out, log2_abund_out,
        gqval_out, qval_out, pep_out, score_out, gscore_out,
        peak_area_out, n_quant_out,
    )

    exp2_row = findfirst(==(UInt32(2)), experiments_out)
    @test exp2_row !== nothing

    @testset "zeroed precursor still counts toward support" begin
        # Both peptides scored the group in experiment 2, so both must be listed.
        # Filling this from S is what produced n_peptides == 0.
        @test count(!ismissing, peptides_out[exp2_row]) == 2
        @test UInt32(2) in skipmissing(peptides_out[exp2_row])
    end

    @testset "quantified count excludes it" begin
        @test n_quant_out[exp2_row] == UInt32(1)
    end

    @testset "runs with no zeroed area are unaffected" begin
        for i in 1:n_experiments
            i == exp2_row && continue
            @test count(!ismissing, peptides_out[i]) == 2
            @test n_quant_out[i] == UInt32(2)
        end
    end

    @testset "quantification itself is unchanged" begin
        # getS is untouched, so the zeroed area stays out of MaxLFQ and the
        # protein profile still sums to the observed precursor intensity.
        @test all(!ismissing, log2_abund_out)
        @test isapprox(sum(exp2, skipmissing(log2_abund_out)), sum(abundance); rtol = 1.0e-6)
    end
end

# `appendPeptides!` was a private inner helper inside `getProtAbundance`
# that no longer has an exported wrapper at the Pioneer.* level — its
# behavior is exercised end-to-end by the `getProtAbundance` testsets
# above (the "peptide lists populated" / "peptide list reflects missing
# data" cases).

end # top-level testset
