using Test
using Pioneer
using Arrow
using DataFrames
using LightGBM
using Random

const PASS1_POOL_TEST_HP = (
    num_iterations = 4, num_leaves = 3, max_depth = 2,
    learning_rate = 0.2, feature_fraction = 1.0,
    min_data_in_leaf = 1, min_gain_to_split = 0.0, num_threads = 1,
)

function pass1_pool_test_write(path, ids; folds = UInt8.((ids .- 1) .% 2),
                               targets = ((ids .- 1) .÷ 2) .% 3 .!= 0)
    Arrow.write(path, (
        precursor_idx = UInt32.(ids), scan_idx = UInt32.(ids .+ 100),
        cv_fold = folds, target = targets,
        signal = Float32.(targets), position = Float64.(ids) ./ 100,
    ))
    return path
end

# The original reservoir copies features at each accepted replacement.
# Compare final slots and RNG state with the deferred gather implementation.
function pass1_pool_reference(paths, features, cap, n0, n1, rng)
    X = (zeros(Float32, min(cap, n0), length(features)),
         zeros(Float32, min(cap, n1), length(features)))
    y = (fill(false, size(X[1], 1)), fill(false, size(X[2], 1)))
    seen = [0, 0]
    for path in paths
        table = Arrow.Table(path)
        for row in eachindex(table.cv_fold)
            fold = Int(table.cv_fold[row]) + 1
            seen[fold] += 1
            slot = seen[fold] <= length(y[fold]) ? seen[fold] : rand(rng, 1:seen[fold])
            slot <= length(y[fold]) || continue
            y[fold][slot] = table.target[row]
            for (j, feature) in enumerate(features)
                X[fold][slot, j] = Float32(getproperty(table, feature)[row])
            end
        end
    end
    return X[1], y[1], X[2], y[2]
end

pass1_pool_predictions(model, X) = clamp.(
    Float32.(vec(LightGBM.predict(model, X; num_threads = 1))), 1f-6, 1f0 - 1f-4,
)

@testset "Pass-1 representative training pool" begin
    mktempdir() do dir
        paths = [pass1_pool_test_write(joinpath(dir, "run$i.arrow"), collect(rows))
                 for (i, rows) in enumerate((1:31, 32:63, 64:96, 97:96))]
        # Empty files need no feature columns during final prediction.
        Arrow.write(paths[end], (precursor_idx = UInt32[], scan_idx = UInt32[],
                                 cv_fold = UInt8[], target = Bool[]))
        features = [:signal, :position]
        @testset "Final reservoir matches serial sampling" begin
            # Exercise replacements across files and a pool larger than the input.
            for cap in (12, 100)
                actual_rng, reference_rng = MersenneTwister(53), MersenneTwister(53)
                expected = pass1_pool_reference(paths, features, cap, 48, 48, reference_rng)
                actual = Pioneer._sample_both_folds(paths, features, cap, 48, 48, actual_rng)
                @test actual == expected
                @test rand(actual_rng, UInt64) == rand(reference_rng, UInt64)
            end
        end

        @testset "Iterations reuse the pool and score every source row" begin
            sample_rng, training_rng = MersenneTwister(53), MersenneTwister(53)
            X0, y0, X1, y1 = Pioneer._sample_both_folds(paths, features, 12, 48, 48, sample_rng)
            result = Pioneer.train_and_predict_pass1_oom!(
                paths; features, compute_infold = true, lgbm_hp = PASS1_POOL_TEST_HP,
                k_per_fold = 12, rng = training_rng, semisupervised = true,
                # All rows stay eligible; identical fits stop after two iterations.
                semisupervised_train_q_threshold = 1f0,
                semisupervised_min_gain = 0.01f0, semisupervised_max_iterations = 4,
            )
            @test result.n_total == 96
            @test result.n_fold0 == result.n_fold1 == 48
            @test result.n_pool_fold0 == result.n_pool_fold1 == 12
            @test result.metric_scope == :training_pool
            @test result.semisupervised_iter == 2
            @test rand(training_rng, UInt64) == rand(sample_rng, UInt64)
            scores = vcat(pass1_pool_predictions(result.cls_trained_on.fold1, X0),
                          pass1_pool_predictions(result.cls_trained_on.fold0, X1))
            targets, qvalues = vcat(y0, y1), similar(scores)
            Pioneer.get_qvalues!(scores, targets, qvalues)
            @test result.target_q01 == count(targets .& (qvalues .<= 0.01f0))
            @test result.decoy_q01 == count((.!targets) .& (qvalues .<= 0.01f0))
            @test 0 < result.target_q01 <= 24
            for path in paths
                source = DataFrame(Arrow.Table(path))
                sidecar = DataFrame(Arrow.Table(path * Pioneer.PASS1_SIDECAR_SUFFIX))
                @test sidecar.precursor_idx == source.precursor_idx
                @test sidecar.scan_idx == source.scan_idx
                for fold in (0, 1)
                    rows = findall(==(fold), source.cv_fold)
                    isempty(rows) && continue
                    X = Matrix{Float32}(source[rows, features])
                    opposite = fold == 0 ? result.cls_trained_on.fold1 : result.cls_trained_on.fold0
                    same = fold == 0 ? result.cls_trained_on.fold0 : result.cls_trained_on.fold1
                    @test sidecar.trace_prob_prepass[rows] == pass1_pool_predictions(opposite, X)
                    @test sidecar.trace_prob_infold[rows] == pass1_pool_predictions(same, X)
                end
            end
        end
    end

    @testset "Retained models have no training datasets" begin
        X = randn(MersenneTwister(41), Float32, 200, 3)
        y = collect(X[:, 1] .+ 0.2f0 .* X[:, 2] .> 0f0)
        baseline = Pioneer.build_lightgbm_classifier(; PASS1_POOL_TEST_HP...)
        LightGBM.fit!(baseline, X, Pioneer._prepare_labels(y); verbosity = -1)
        model = Pioneer._fit_pass1_booster(X, y, PASS1_POOL_TEST_HP)
        @test !isempty(baseline.booster.datasets)
        @test isempty(model.booster.datasets)
        @test LightGBM.predict(model, X; num_threads = 1) == LightGBM.predict(baseline, X; num_threads = 1)
        @test LightGBM.split_importance(model) == LightGBM.split_importance(baseline)
    end

    @testset "Filtering leaves excluded rows available for later iterations" begin
        X, y = reshape(Float32[2, 0, 1, 0], 4, 1), Bool[true, false, true, false]
        original_X, original_y = copy(X), copy(y)
        fit = Pioneer._fit_pass1_pool_fold(X, y, Bool[false, true, true, true], PASS1_POOL_TEST_HP)
        @test (fit.targets, fit.decoys) == (1, 2)
        @test X == original_X && y == original_y
        decoys_only = Pioneer._fit_pass1_pool_fold(X, y, .!y, PASS1_POOL_TEST_HP)
        @test (decoys_only.targets, decoys_only.decoys) == (0, 2)
        @test decoys_only.classifier.value == 0f0
        @test X == original_X && y == original_y
    end

    @testset "Constant folds preserve OOF direction and optional in-fold scores" begin
        mktempdir() do dir
            folds = UInt8.((0:9) .% 2)
            targets = folds .== 0
            path = pass1_pool_test_write(joinpath(dir, "constant.arrow"), collect(1:10); folds, targets)
            for compute_infold in (true, false)
                Pioneer.train_and_predict_pass1_oom!(
                    [path]; features = [:signal], compute_infold,
                    lgbm_hp = PASS1_POOL_TEST_HP, k_per_fold = 2,
                )
                sidecar = DataFrame(Arrow.Table(path * Pioneer.PASS1_SIDECAR_SUFFIX))
                @test sidecar.trace_prob_prepass[targets] == fill(1f-6, 5)
                @test sidecar.trace_prob_prepass[.!targets] == fill(1f0 - 1f-4, 5)
                if compute_infold
                    @test sidecar.trace_prob_infold[targets] == fill(1f0 - 1f-4, 5)
                    @test sidecar.trace_prob_infold[.!targets] == fill(1f-6, 5)
                else
                    @test all(isnan, sidecar.trace_prob_infold)
                end
            end
        end
    end
end
