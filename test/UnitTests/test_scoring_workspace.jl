# Unit tests for InMemoryScoringWorkspace
# (src/utils/ML/PSMScoring/workspace.jl).
#
# Tests CV fold setup, phase-dispatched accessors, and output array management.
#
# Run standalone: julia --project=. test/UnitTests/test_scoring_workspace.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
    using DataFrames
end

# Helper: build a PSM DataFrame and ScoringConfig suitable for workspace tests
function _make_workspace_psms(n::Int; n_folds::Int=2)
    df = DataFrame(
        target            = [i <= n÷2 for i in 1:n],
        cv_fold           = UInt8.(mod.(0:(n-1), n_folds)),
        precursor_idx     = UInt32.(1:n),
        ms_file_idx       = UInt32.(ones(Int, n)),
        irt_pred          = rand(Float32, n),
        isotopes_captured = [(Int8(0), Int8(1)) for _ in 1:n],
        score             = rand(Float32, n),
    )
    container = Pioneer.DataFramePSMContainer(df)
    config = Pioneer.ScoringConfig(
        Pioneer.LightGBMScorer(Dict{Symbol,Any}()),
        Pioneer.AllDataSelection(),
        Pioneer.StaticFeatureSelection([:score]),
        200
    )
    return container, config
end

@testset "ScoringWorkspace" begin

# ═══════════════════════════════════════════════════════════════════════════════
# setup_scoring_workspace — in-memory path
# ═══════════════════════════════════════════════════════════════════════════════

@testset "setup_scoring_workspace" begin
    container, config = _make_workspace_psms(20; n_folds=2)
    ws = Pioneer.setup_scoring_workspace(container, config)

    @testset "type" begin
        @test ws isa Pioneer.InMemoryScoringWorkspace
    end

    @testset "cv_folds" begin
        folds = Pioneer.get_cv_folds(ws)
        @test sort(folds) == UInt8[0, 1]
    end

    @testset "fold indices partition all rows" begin
        idx0 = Pioneer.get_test_indices(ws, UInt8(0))
        idx1 = Pioneer.get_test_indices(ws, UInt8(1))
        @test sort(vcat(idx0, idx1)) == collect(1:20)
        @test isempty(intersect(idx0, idx1))
    end

    @testset "train indices complement test indices" begin
        for fold in Pioneer.get_cv_folds(ws)
            test = Pioneer.get_test_indices(ws, fold)
            train = Pioneer.get_train_indices(ws, fold)
            @test sort(vcat(test, train)) == collect(1:20)
            @test isempty(intersect(test, train))
        end
    end

    @testset "output arrays initialized to zero" begin
        @test length(ws.prob_test) == 20
        @test length(ws.prob_train) == 20
        @test all(ws.prob_test .== 0.0f0)
        @test all(ws.prob_train .== 0.0f0)
    end

    @testset "models dict starts empty" begin
        @test isempty(ws.models)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Phase-dispatched accessors
# ═══════════════════════════════════════════════════════════════════════════════

@testset "phase accessors" begin
    container, config = _make_workspace_psms(20; n_folds=2)
    ws = Pioneer.setup_scoring_workspace(container, config)
    fold = UInt8(0)

    @testset "TrainingPhase view matches train indices" begin
        train_view = Pioneer.get_phase_view(ws, Pioneer.TrainingPhase(), fold)
        train_idx = Pioneer.get_phase_indices(ws, Pioneer.TrainingPhase(), fold)
        @test Pioneer.nrows(train_view) == length(train_idx)
    end

    @testset "PredictionPhase view matches test indices" begin
        test_view = Pioneer.get_phase_view(ws, Pioneer.PredictionPhase(), fold)
        test_idx = Pioneer.get_phase_indices(ws, Pioneer.PredictionPhase(), fold)
        @test Pioneer.nrows(test_view) == length(test_idx)
    end

    @testset "phase outputs are the correct arrays" begin
        @test Pioneer.get_phase_output(ws, Pioneer.TrainingPhase()) === ws.prob_train
        @test Pioneer.get_phase_output(ws, Pioneer.PredictionPhase()) === ws.prob_test
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Model storage
# ═══════════════════════════════════════════════════════════════════════════════

@testset "model storage" begin
    container, config = _make_workspace_psms(10; n_folds=2)
    ws = Pioneer.setup_scoring_workspace(container, config)
    fold = UInt8(0)

    @testset "init_fold_models TrainingPhase" begin
        models = Pioneer.init_fold_models(ws, Pioneer.TrainingPhase(), fold, 1)
        @test length(models) == 1
        @test models isa Vector{Any}
    end

    @testset "commit_fold! stores models" begin
        models = Any["fake_model"]
        Pioneer.commit_fold!(ws, Pioneer.TrainingPhase(), fold, models)
        @test Pioneer.get_fold_models(ws, fold) === models
    end

    @testset "init_fold_models PredictionPhase retrieves stored" begin
        stored = Pioneer.init_fold_models(ws, Pioneer.PredictionPhase(), fold, 1)
        @test stored[1] == "fake_model"
    end

    @testset "get_all_models" begin
        all_models = Pioneer.get_all_models(ws)
        @test haskey(all_models, UInt8(0))
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# store_final_predictions!
# ═══════════════════════════════════════════════════════════════════════════════

@testset "store_final_predictions!" begin
    container, config = _make_workspace_psms(10; n_folds=2)

    # Add trace_prob column (normally done by prepare_training_data!)
    Pioneer.set_column!(container, :trace_prob, zeros(Float32, 10))
    Pioneer.set_column!(container, :q_value, zeros(Float64, 10))

    ws = Pioneer.setup_scoring_workspace(container, config)
    fold = UInt8(0)
    test_idx = Pioneer.get_test_indices(ws, fold)

    # Simulate prediction: set trace_prob for test fold's rows
    probs = Float32.(0.1:0.1:1.0)
    for (i, idx) in enumerate(test_idx)
        container.data.trace_prob[idx] = probs[i]
    end

    Pioneer.store_final_predictions!(ws, fold)

    # prob_test should have the predictions at the test fold indices
    for (i, idx) in enumerate(test_idx)
        @test ws.prob_test[idx] == probs[i]
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# prepare_training_data!
# ═══════════════════════════════════════════════════════════════════════════════

@testset "prepare_training_data!" begin
    container, config = _make_workspace_psms(10; n_folds=2)

    Pioneer.prepare_training_data!(container, config)

    @testset "adds trace_prob column" begin
        @test Pioneer.has_column(container, :trace_prob)
        @test all(Pioneer.get_column(container, :trace_prob) .== 0.0f0)
    end

    @testset "adds q_value column" begin
        @test Pioneer.has_column(container, :q_value)
        @test all(Pioneer.get_column(container, :q_value) .== 0.0)
    end

    @testset "sorts by isotopes_captured, precursor_idx, ms_file_idx" begin
        # After sorting, precursor_idx should be sorted (since isotopes_captured
        # and ms_file_idx are constant in our test data)
        pidx = Pioneer.get_column(container, :precursor_idx)
        @test issorted(pidx)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Three-fold CV
# ═══════════════════════════════════════════════════════════════════════════════

@testset "three-fold CV" begin
    container, config = _make_workspace_psms(30; n_folds=3)
    ws = Pioneer.setup_scoring_workspace(container, config)

    folds = Pioneer.get_cv_folds(ws)
    @test length(folds) == 3

    # Each fold should have 10 rows
    for fold in folds
        @test length(Pioneer.get_test_indices(ws, fold)) == 10
        @test length(Pioneer.get_train_indices(ws, fold)) == 20
    end
end

end # top-level testset
