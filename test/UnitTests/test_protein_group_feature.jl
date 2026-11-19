# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

using Test
using Pioneer
using Arrow
using DataFrames
using Random
using Statistics

const PG_TEST_HP = (
    num_iterations = 4, num_leaves = 3, max_depth = 2,
    learning_rate = 0.2, feature_fraction = 1.0,
    min_data_in_leaf = 1, min_gain_to_split = 0.0, num_threads = 1,
)

@testset "PG_MAX protein-group column" begin
    # Precursors 1-2 are two charges of peptide 10; 3 is peptide 11 (same group 1);
    # 4-5 are group 2; 6 is alone in group 3.
    group = UInt32[1, 1, 1, 2, 2, 3]
    peptide = UInt32[10, 10, 11, 12, 13, 14]
    # Two rows for precursor 1 (its best counts), one for the rest.
    pidx = UInt32[1, 1, 2, 3, 4, 5, 6]
    target = Bool[true, true, true, true, false, true, false]
    scores = Float32[0.2, 0.9, 0.5, 0.4, 0.3, 0.8, 0.1]
    floor = Float32(median([0.3, 0.1]))   # decoy precursors 4 and 6

    @testset "v1 leaves out the precursor's own peptide" begin
        keys = Pioneer.ProteinGroupKeys(:v1, group, peptide, Float32[])
        col = Pioneer.protein_group_columns(keys, pidx, target, scores)
        @test size(col) == (7, 1)
        # Peptide 10 (precursors 1, 2) sees only peptide 11 (0.4); peptide 11 sees peptide 10's best (0.9).
        @test col[:, 1] == Float32[0.4, 0.4, 0.4, 0.9, 0.8, 0.3, floor]
    end

    @testset "v0 includes the own peptide" begin
        keys = Pioneer.ProteinGroupKeys(:v0, group, peptide, Float32[])
        col = Pioneer.protein_group_columns(keys, pidx, target, scores)
        @test col[:, 1] == Float32[0.9, 0.9, 0.9, 0.9, 0.8, 0.8, 0.1]
    end

    @testset "v2 adds the group size" begin
        sizes = Pioneer._distinct_sequences_per_group(group, peptide)
        @test sizes == Float32[2, 2, 2, 2, 2, 1]
        keys = Pioneer.ProteinGroupKeys(:v2, group, peptide, sizes)
        col = Pioneer.protein_group_columns(keys, pidx, target, scores)
        @test size(col, 2) == 2
        @test col[:, 2] == Float32[2, 2, 2, 2, 2, 2, 1]
        @test Pioneer.pg_feature_names(keys) == [:pg_max, :pg_size]
    end

    @testset "negative control keeps group sizes, peptide units and folds" begin
        rng = MersenneTwister(7)
        n = 400
        g = UInt32.(rand(rng, 1:40, n))
        base = UInt32.((0:n - 1) .÷ 2)                # precursors (1,2), (3,4), ... share a base peptide
        g[2:2:end] .= g[1:2:end]                      # a base peptide has one accession string
        fold = UInt8.(g .% 2)
        p = Pioneer._permute_groups_within_fold(g, base, fold)
        counts(v) = (d = Dict{UInt32, Int}(); foreach(x -> d[x] = get(d, x, 0) + 1, v); d)
        @test sort(collect(values(counts(g)))) == sort(collect(values(counts(p))))
        @test all(p[2:2:end] .== p[1:2:end])          # units stay together
        new_fold = Dict(p[i] => fold[i] for i in 1:n)
        @test all(new_fold[p[i]] == fold[i] for i in 1:n)  # groups stay within one fold
        @test p != g
    end

    @testset "mode parsing" begin
        withenv("PIONEER_PG_MAX" => nothing) do
            @test Pioneer.pg_max_mode() == :off
        end
        withenv("PIONEER_PG_MAX" => "V1") do
            @test Pioneer.pg_max_mode() == :v1
        end
        withenv("PIONEER_PG_MAX" => "bogus") do
            @test_throws ArgumentError Pioneer.pg_max_mode()
        end
    end
end

@testset "PG_MAX inside the Pass-1 loop" begin
    mktempdir() do dir
        rng = MersenneTwister(11)
        n_prec = 240
        group = UInt32.((0:n_prec - 1) .÷ 4 .+ 1)          # groups of 4 precursors
        peptide = UInt32.(1:n_prec)
        fold_of = UInt8.(group .% 2)                        # a group sits in one fold
        is_target = [isodd(g) || rand(rng) < 0.5 for g in group]
        paths = String[]
        for f in 1:2
            ids = collect(1:n_prec)
            signal = Float32.(is_target[ids]) .+ 0.8f0 .* randn(rng, Float32, n_prec)
            path = joinpath(dir, "run$f.arrow")
            Arrow.write(path, (
                precursor_idx = UInt32.(ids), scan_idx = UInt32.(ids),
                cv_fold = fold_of[ids], target = is_target[ids],
                signal = signal, noise = randn(rng, Float32, n_prec),
            ))
            push!(paths, path)
        end
        for mode in (:v1, :vr)
            keys = Pioneer.ProteinGroupKeys(mode, group, peptide, Float32[])
            result = Pioneer.train_and_predict_pass1_oom!(
                paths; features = [:signal, :noise], compute_infold = true,
                lgbm_hp = PG_TEST_HP, k_per_fold = 1_000, semisupervised = true,
                semisupervised_train_q_threshold = 1f0, semisupervised_min_gain = 0.0f0,
                semisupervised_max_iterations = 3, pg_keys = keys,
            )
            @test result.n_total == 2 * n_prec
            # The kept models' feature list includes PG_MAX exactly when an iteration after the first was kept.
            @test (:pg_max in result.available_features) == (result.semisupervised_iter > 1)
            for path in paths
                sidecar = DataFrame(Arrow.Table(path * Pioneer.PASS1_SIDECAR_SUFFIX))
                @test nrow(sidecar) == n_prec
                @test all(0 .< sidecar.trace_prob_prepass .< 1)
                @test all(isfinite, sidecar.trace_prob_infold)
            end
        end
    end
end
