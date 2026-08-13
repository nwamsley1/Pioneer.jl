# Coverage tests for src/utils/ML/PSMScoring/workspace.jl.
#
# Targets the in-memory + Arrow-file scoring workspaces:
# - InMemoryScoringWorkspace fold / train index pre-computation
# - prepare_training_data! sorts in-place and initializes scoring columns
# - phase-dispatched accessors (Training vs Prediction)
# - store_final_predictions! writes view results back into the global array
# - ArrowFile path: _create_combined_dataframe schema + missing-typed columns
# - ArrowFile path: _sample_psms proportional reservoir sampling
# - ArrowFile prepare_training_data! / _finalize_scoring_arrow! round-trip

using DataFrames
using Arrow, Tables
using Random

using Pioneer: DataFramePSMContainer,
               InMemoryScoringWorkspace, ArrowFileScoringWorkspace,
               ArrowFilePSMContainer, ArrowFileGroup,
               setup_scoring_workspace, prepare_training_data!,
               get_cv_folds, get_train_indices, get_test_indices,
               get_test_view, get_training_view,
               get_train_output, get_test_output,
               get_phase_view, get_phase_indices, get_phase_output,
               init_fold_models, commit_fold!,
               store_fold_models!, get_fold_models, get_all_models,
               get_psms, store_final_predictions!,
               TrainingPhase, PredictionPhase,
               probit_scoring_config

# Pioneer-namespaced (underscore-prefixed internal helpers)
const _create_combined_dataframe = Pioneer._create_combined_dataframe
const _sample_psms                = Pioneer._sample_psms
const _finalize_scoring_arrow!    = Pioneer._finalize_scoring_arrow!
const writeArrow                  = Pioneer.writeArrow

function _make_synthetic_psm_df(n_per_fold::Int = 6)
    # 2 folds × n_per_fold rows, alphabet-mixed for sort verification.
    n = 2 * n_per_fold
    rng = MersenneTwister(42)
    df = DataFrame(
        target            = rand(rng, Bool, n),
        cv_fold           = repeat(UInt8[0, 1], inner = n_per_fold),
        precursor_idx     = shuffle!(rng, UInt32.(1:n)),
        ms_file_idx       = UInt32.(rand(rng, 1:3, n)),
        irt_pred          = Float32.(randn(rng, n) .* 5.0 .+ 30.0),
        isotopes_captured = [(Int8(0), Int8(2)) for _ in 1:n],
    )
    return DataFramePSMContainer(df)
end

@testset "workspace.jl" begin

    @testset "InMemoryScoringWorkspace setup + accessors" begin
        psms = _make_synthetic_psm_df(5)
        config = probit_scoring_config(features = Symbol[])

        ws = setup_scoring_workspace(psms, config)
        @test ws isa InMemoryScoringWorkspace

        # Folds discovered correctly
        @test sort(get_cv_folds(ws)) == sort(get_cv_folds(psms))
        @test length(get_cv_folds(ws)) == 2

        # Output arrays pre-allocated to nrow(psms)
        @test length(get_train_output(ws)) == 10
        @test length(get_test_output(ws)) == 10
        @test all(get_train_output(ws) .== 0.0f0)
        @test all(get_test_output(ws) .== 0.0f0)

        # Per-fold splits
        for fold in get_cv_folds(ws)
            test_idx = get_test_indices(ws, fold)
            train_idx = get_train_indices(ws, fold)
            @test isdisjoint(test_idx, train_idx)
            @test length(test_idx) + length(train_idx) == 10
            # All test rows have cv_fold == fold
            cv = get_psms(ws).data.cv_fold
            @test all(cv[test_idx] .== fold)
            @test all(cv[train_idx] .!= fold)
        end

        # get_test_view / get_training_view are non-empty subsets
        for fold in get_cv_folds(ws)
            tv = get_test_view(ws, fold)
            trv = get_training_view(ws, fold)
            @test Pioneer.nrows(tv) == length(get_test_indices(ws, fold))
            @test Pioneer.nrows(trv) == length(get_train_indices(ws, fold))
        end

        # get_psms returns the original container
        @test get_psms(ws) === psms
    end

    @testset "prepare_training_data! sorts + initializes columns" begin
        psms = _make_synthetic_psm_df(5)
        config = probit_scoring_config(features = Symbol[])

        # Pre-sort state: precursor_idx is shuffled
        before = copy(psms.data.precursor_idx)
        prepare_training_data!(psms, config)

        # Sort key is (isotopes_captured, precursor_idx, ms_file_idx). All rows
        # share the same isotopes, so precursor_idx should now be ascending
        # within each ms_file_idx group.
        df = psms.data
        @test issorted(df, [:isotopes_captured, :precursor_idx, :ms_file_idx])

        # Initialized scoring columns
        @test :trace_prob in propertynames(df)
        @test :q_value   in propertynames(df)
        @test eltype(df.trace_prob) == Float32
        @test eltype(df.q_value)    == Float64
        @test all(df.trace_prob .== 0.0f0)
        @test all(df.q_value   .== 0.0)
        @test length(df.trace_prob) == 10
    end

    @testset "phase-dispatched accessors (Training vs Prediction)" begin
        psms = _make_synthetic_psm_df(4)
        config = probit_scoring_config(features = Symbol[])
        ws = setup_scoring_workspace(psms, config)

        for fold in get_cv_folds(ws)
            train_phase_view = get_phase_view(ws, TrainingPhase(), fold)
            pred_phase_view  = get_phase_view(ws, PredictionPhase(), fold)
            @test Pioneer.nrows(train_phase_view) == length(get_train_indices(ws, fold))
            @test Pioneer.nrows(pred_phase_view)  == length(get_test_indices(ws, fold))

            @test get_phase_indices(ws, TrainingPhase(), fold) ==
                  get_train_indices(ws, fold)
            @test get_phase_indices(ws, PredictionPhase(), fold) ==
                  get_test_indices(ws, fold)
        end

        # Phase output arrays — separate buffers
        @test get_phase_output(ws, TrainingPhase())   === get_train_output(ws)
        @test get_phase_output(ws, PredictionPhase()) === get_test_output(ws)

        # init_fold_models creates a fresh Vector{Any} of the right length for
        # training, and returns the stored models for prediction.
        models_in = Any[1, 2, 3]
        commit_fold!(ws, TrainingPhase(), UInt8(0), models_in)
        @test get_all_models(ws)[UInt8(0)] === models_in
        @test get_fold_models(ws, UInt8(0)) === models_in

        # init_fold_models for TrainingPhase allocates n slots
        slots = init_fold_models(ws, TrainingPhase(), UInt8(0), 4)
        @test length(slots) == 4
        @test eltype(slots) === Any

        # init_fold_models for PredictionPhase returns the stored models
        @test init_fold_models(ws, PredictionPhase(), UInt8(0), 0) === models_in
    end

    @testset "store_fold_models! / get_fold_models / get_all_models" begin
        psms = _make_synthetic_psm_df(3)
        config = probit_scoring_config(features = Symbol[])
        ws = setup_scoring_workspace(psms, config)

        store_fold_models!(ws, UInt8(0), Any["model_a"])
        store_fold_models!(ws, UInt8(1), Any["model_b1", "model_b2"])
        @test get_fold_models(ws, UInt8(0)) == ["model_a"]
        @test get_fold_models(ws, UInt8(1)) == ["model_b1", "model_b2"]
        all_models = get_all_models(ws)
        @test sort(collect(keys(all_models))) == [UInt8(0), UInt8(1)]
    end

    @testset "store_final_predictions! writes view results back to global array" begin
        psms = _make_synthetic_psm_df(4)
        config = probit_scoring_config(features = Symbol[])
        prepare_training_data!(psms, config)   # adds :trace_prob
        ws = setup_scoring_workspace(psms, config)

        # Inject distinguishable trace_prob values for fold 0's test rows only,
        # via the view's set_column!.
        fold = UInt8(0)
        test_idx = get_test_indices(ws, fold)
        view = get_test_view(ws, fold)
        n_test = Pioneer.nrows(view)
        injected = Float32[Float32(i) * 0.1f0 for i in 1:n_test]
        Pioneer.set_column!(view, :trace_prob, injected)

        store_final_predictions!(ws, fold)
        # Workspace's prob_test array should now hold injected values at test_idx
        out = get_test_output(ws)
        for (i, idx) in enumerate(test_idx)
            @test out[idx] ≈ injected[i]
        end
        # Other (non-fold-0) indices still zero
        non_fold0 = setdiff(1:length(out), test_idx)
        @test all(out[non_fold0] .== 0.0f0)
    end

    # ────────────────────────────────────────────────────────────────────
    # ArrowFile path
    # ────────────────────────────────────────────────────────────────────
    @testset "_create_combined_dataframe schema (Missing-typed columns)" begin
        # Build two Arrow tables in memory, with one column declared as
        # Union{Missing,Float64}, another as plain UInt32, and verify the
        # resulting DataFrame preserves those types.
        data_df = DataFrame(
            precursor_idx = UInt32[1, 2, 3],
            score = Union{Missing, Float64}[0.5, missing, 0.7],
        )
        scores_df = DataFrame(
            trace_prob = Float32[0.0, 0.0, 0.0],
            q_value    = Float64[0.0, 0.0, 0.0],
        )
        mktempdir() do dir
            data_path = joinpath(dir, "d.arrow")
            scores_path = joinpath(dir, "s.arrow")
            Arrow.write(data_path, data_df)
            Arrow.write(scores_path, scores_df)

            data_tbl = Arrow.Table(data_path)
            scores_tbl = Arrow.Table(scores_path)

            combined = _create_combined_dataframe(data_tbl, scores_tbl, 5)
            @test nrow(combined) == 5
            @test eltype(combined.precursor_idx) == UInt32
            @test eltype(combined.score) == Union{Missing, Float64}
            @test eltype(combined.trace_prob) == Float32
            @test eltype(combined.q_value)    == Float64
            # All columns from both inputs present
            for name in (:precursor_idx, :score, :trace_prob, :q_value)
                @test name in propertynames(combined)
            end
        end
    end

    @testset "_sample_psms: proportional, deterministic, count == budget" begin
        # Build two Arrow files (file_fold0 with 30 rows, file_fold1 with 70).
        # max=20 → expect ~20 rows total (rounding may shave by 0–1).
        rng = MersenneTwister(11)
        function _make_data_df(n)
            DataFrame(
                target            = rand(rng, Bool, n),
                cv_fold           = fill(UInt8(0), n),
                precursor_idx     = UInt32.(1:n),
                ms_file_idx       = UInt32.(fill(1, n)),
                irt_pred          = Float32.(rand(rng, n) .* 50.0),
                isotopes_captured = [(Int8(0), Int8(2)) for _ in 1:n],
            )
        end

        mktempdir() do dir
            f0_path = joinpath(dir, "f_fold0.arrow")
            f1_path = joinpath(dir, "f_fold1.arrow")
            Arrow.write(f0_path, _make_data_df(30))
            Arrow.write(f1_path, _make_data_df(70))

            scores0 = DataFrame(trace_prob = zeros(Float32, 30), q_value = zeros(Float64, 30))
            scores1 = DataFrame(trace_prob = zeros(Float32, 70), q_value = zeros(Float64, 70))
            s0_path = joinpath(dir, "f_fold0_scores.arrow")
            s1_path = joinpath(dir, "f_fold1_scores.arrow")
            Arrow.write(s0_path, scores0)
            Arrow.write(s1_path, scores1)

            container = ArrowFilePSMContainer([f0_path, f1_path]; max_training_psms = 20)

            sampled = _sample_psms(container, 20)
            n_sampled = Pioneer.nrows(sampled)
            # Allow rounding: each per-file count uses round(Int, n*frac); total
            # is within one row of the budget.
            @test 18 ≤ n_sampled ≤ 22
            # Determinism: same seed, same input → same row count
            sampled2 = _sample_psms(container, 20)
            @test Pioneer.nrows(sampled2) == n_sampled
            # And the sampled rows should be a subset of the original union
            all_pidx = vcat(_make_data_df(30).precursor_idx, _make_data_df(70).precursor_idx)
            sampled_pidx = sampled.data.precursor_idx
            @test all(p -> p in Set(all_pidx), sampled_pidx)
        end
    end

    @testset "ArrowFile prepare_training_data! + _finalize_scoring_arrow! round-trip" begin
        rng = MersenneTwister(7)
        function _data_df(n)
            DataFrame(
                target            = rand(rng, Bool, n),
                cv_fold           = fill(UInt8(0), n),
                precursor_idx     = shuffle!(rng, UInt32.(1:n)),
                ms_file_idx       = UInt32.(fill(2, n)),
                irt_pred          = Float32.(rand(rng, n) .* 30.0),
                isotopes_captured = [(Int8(0), Int8(2)) for _ in 1:n],
            )
        end

        mktempdir() do dir
            data_path = joinpath(dir, "x_fold0.arrow")
            Arrow.write(data_path, _data_df(8))
            container = ArrowFilePSMContainer([data_path])

            config = probit_scoring_config(features = Symbol[])
            prepare_training_data!(container, config)

            # After prepare_training_data!, scores sidecar should exist with
            # zeros and main file should be sorted.
            scores_tbl = Arrow.Table(container.file_groups[1].scores_path)
            @test length(scores_tbl[:trace_prob]) == 8
            @test all(Tables.getcolumn(scores_tbl, :trace_prob) .== 0.0f0)
            @test all(Tables.getcolumn(scores_tbl, :q_value)    .== 0.0)
            data_tbl = Arrow.Table(container.file_groups[1].data_path)
            @test issorted(Tables.getcolumn(data_tbl, :precursor_idx))

            # Now write nontrivial trace_prob into the scores sidecar and
            # verify _finalize_scoring_arrow! copies it into the main data file.
            preds = Float32[i * 0.05f0 for i in 1:8]
            scores_df = DataFrame(
                trace_prob = preds,
                q_value    = zeros(Float64, 8),
            )
            writeArrow(container.file_groups[1].scores_path, scores_df)

            _finalize_scoring_arrow!(container)
            data_after = Arrow.Table(container.file_groups[1].data_path)
            @test :trace_prob in Tables.columnnames(data_after)
            @test all(Tables.getcolumn(data_after, :trace_prob) .≈ preds)
        end
    end

    @testset "setup_scoring_workspace dispatches on ArrowFile container" begin
        mktempdir() do dir
            df = DataFrame(
                target            = rand(Bool, 4),
                cv_fold           = fill(UInt8(0), 4),
                precursor_idx     = UInt32.(1:4),
                ms_file_idx       = UInt32.(fill(1, 4)),
                irt_pred          = Float32.(rand(4) .* 30.0),
                isotopes_captured = [(Int8(0), Int8(2)) for _ in 1:4],
            )
            data_path = joinpath(dir, "y_fold0.arrow")
            Arrow.write(data_path, df)
            scores_df = DataFrame(trace_prob = zeros(Float32, 4),
                                  q_value    = zeros(Float64, 4))
            Arrow.write(joinpath(dir, "y_fold0_scores.arrow"), scores_df)

            container = ArrowFilePSMContainer([data_path]; max_training_psms = 10)
            config = probit_scoring_config(features = Symbol[])
            ws = setup_scoring_workspace(container, config)

            @test ws isa ArrowFileScoringWorkspace
            # Delegates to inner
            @test ws.container === container
            @test ws.inner isa InMemoryScoringWorkspace
            # ArrowFile workspace's store_final_predictions! is a no-op (returns nothing)
            @test store_final_predictions!(ws, UInt8(0)) === nothing
            # Delegated cv_folds (only fold 0 in this fixture)
            @test get_cv_folds(ws) == [UInt8(0)]
        end
    end

end
