using Test, Pioneer, LightGBM, DataFrames, Random

@testset "LightGBM prediction-only models" begin
    rng = MersenneTwister(23)
    x = rand(rng, Float32, 600, 3)
    y = Int.(x[:, 1] .+ x[:, 2] .> 1)
    model = Pioneer.build_lightgbm_classifier(num_iterations=5, min_data_in_leaf=10,
        min_gain_to_split=0.0, num_threads=2)
    LightGBM.fit!(model, x, y; verbosity=-1)
    predictions = LightGBM.predict(model, x)
    gains = LightGBM.gain_importance(model)
    trained = model.booster
    datasets = copy(trained.datasets)
    @test !isempty(datasets)
    @test Pioneer._detach_lightgbm_training_data!(model) === model
    @test isempty(model.booster.datasets)
    @test trained.handle == C_NULL
    @test all(ds -> ds.handle == C_NULL, datasets)
    GC.gc()
    @test LightGBM.predict(model, x) == predictions
    # LightGBM model text rounds split gains, but preserves prediction values.
    @test LightGBM.gain_importance(model) ≈ gains rtol=1e-5
    @test Pioneer._detach_lightgbm_training_data!(model) === model
    @test LightGBM.predict(model, x) == predictions

    wrapper = Pioneer.fit_lightgbm_model(
        Pioneer.build_lightgbm_classifier(num_iterations=5, min_data_in_leaf=10,
            min_gain_to_split=0.0, num_threads=2), DataFrame(x, :auto), y)
    @test isempty(wrapper.booster.booster.datasets)
    @test LightGBM.predict(wrapper.booster, x) == predictions
end

@testset "MBR detached OOF models preserve predictions" begin
    n = 600
    rng = MersenneTwister(41)
    frame = DataFrame(real=rand(rng, Float32, n), cf1=rand(rng, Float32, n),
        cf2=rand(rng, Float32, n), cf3=rand(rng, Float32, n))
    features = [:real]
    false_features = [[:cf1], [:cf2], [:cf3]]
    folds = UInt8.(mod.(1:n, 2))
    positives = BitVector(mod.(1:n, 3) .!= 0)
    decoys = .!positives
    present = trues(n, 3)
    test_rows = Pioneer._mbr_test_rows_by_fold(folds, present)
    scores, model = mktemp() do path, io
        old_file, old_level = Pioneer.DEBUG_FILE[], Pioneer.DEBUG_FILE_LEVEL[]
        try
            Pioneer.DEBUG_FILE[] = io
            Pioneer.DEBUG_FILE_LEVEL[] = 1
            result = Pioneer._mbr_fit_oof_iteration(
                frame, features, false_features, folds, positives, decoys, present, test_rows)
            log = read(path, String)
            for fold in 0:1, phase in ("training gather complete", "fit complete", "prediction complete")
                @test occursin("MBR transfer model OOF fold $fold $phase:", log)
            end
            @test occursin("threads=", log)
            return result
        finally
            Pioneer.DEBUG_FILE[] = old_file
            Pioneer.DEBUG_FILE_LEVEL[] = old_level
        end
    end
    @test isempty(model.booster.datasets)
    for (i, test_fold) in enumerate(UInt8[0, 1])
        training = Pioneer._mbr_training_rows(folds, positives, decoys, present, UInt8(1)-test_fold)
        reference = Pioneer.build_lightgbm_classifier(; Pioneer.SHARED_LGBM_HP...)
        LightGBM.fit!(reference,
            Pioneer._mbr_gather_feature_rows(frame, features, false_features, training.rows, n),
            Pioneer._prepare_labels(training.labels); verbosity=-1)
        expected = vec(LightGBM.predict(reference,
            Pioneer._mbr_gather_feature_rows(frame, features, false_features, test_rows[i], n)))
        @test scores[test_rows[i]] == Float32.(expected)
        Pioneer._detach_lightgbm_training_data!(reference)
    end
end

@testset "MainSearch retains only prediction models" begin
    n = 22000
    target = mod.(1:n, 3) .!= 0
    psms = DataFrame(target=target, cv_fold=UInt8.(mod.(1:n, 2)), signal=Float32.(target))
    scores, infold, model, info = Pioneer.train_psm_classifier_with_fallback(psms;
        features=[:signal], compute_infold=true,
        lgbm_hp=(num_iterations=5, num_threads=2, feature_fraction=1.0))
    @test model !== nothing
    @test isempty(model.booster.datasets)
    @test all(p -> p.kind == :lgbm && isempty(p.model.booster.datasets), info.predictor.fold_predictors)
    GC.gc()
    @test Pioneer.predict_psm_classifier_scores(psms, info.predictor) == scores
    @test all(isfinite, infold)
end
