using Test, Pioneer, LightGBM, DataFrames, Random

# MainSearch fills LightGBM matrices straight from the PSM table's columns and predicts in fixed-size batches
# (fill_rows!, rows_matrix!, predict_rows!) instead of materializing every PSM. Every result here must equal
# (==) what the whole-table matrix (feature_matrix) gives.

function _batched_test_table(n; seed = 11)
    rng = MersenneTwister(seed)
    target = rand(rng, Bool, n)
    s = Float32.(target)
    DataFrame(
        target = target, cv_fold = UInt8.(rand(rng, 0:1, n)),
        f_f32 = s .+ randn(rng, Float32, n),
        f_f16 = Float16.(0.5f0 .* s .+ randn(rng, Float32, n)),
        f_u8 = UInt8.(rand(rng, 0:20, n)),
        f_bool = rand(rng, Bool, n),
        f_i64 = rand(rng, -5:5, n),
        f_miss = Union{Missing, Float32}[rand(rng) < 0.2 ? missing : x for x in (s .+ randn(rng, Float32, n))],
        f_f64 = Float64.(0.3 .* s .+ randn(rng, n)),
    )
end
const _BATCHED_FEATS = [:f_f32, :f_f16, :f_u8, :f_bool, :f_i64, :f_miss, :f_f64]
_cols(df) = AbstractVector[df[!, f] for f in _BATCHED_FEATS]

@testset "fill_rows! / rows_matrix! equal the whole-table matrix rows" begin
    df = _batched_test_table(5_000)
    full = Pioneer.feature_matrix(df, _BATCHED_FEATS)
    cols = _cols(df)
    for rows in (collect(1:5_000), shuffle(MersenneTwister(1), 1:5_000)[1:1_234], [7, 7, 3, 5_000, 1], Int[])
        @test Pioneer.rows_matrix(cols, rows) == full[rows, :]
        buf = Float32[]
        GC.@preserve buf begin
            @test Pioneer.rows_matrix!(buf, cols, rows) == full[rows, :]
        end
    end
    @test_throws BoundsError Pioneer.rows_matrix(cols, [1, 5_001])
    @test_throws BoundsError Pioneer.rows_matrix(cols, [0])
    @test_throws DimensionMismatch Pioneer.fill_rows!(Matrix{Float32}(undef, 2, 3), cols, [1, 2])
end

@testset "predict_rows! writes each batch's predictions to the right rows" begin
    df = _batched_test_table(20_000)
    full = Pioneer.feature_matrix(df, _BATCHED_FEATS)
    cols = _cols(df)
    model = Pioneer.build_lightgbm_classifier(num_iterations = 20, min_data_in_leaf = 20,
                                              min_gain_to_split = 0.0, num_threads = 2)
    LightGBM.fit!(model, full[1:10_000, :], Int.(df.target[1:10_000]); verbosity = -1)
    rows = shuffle(MersenneTwister(2), 1:20_000)[1:12_345]                 # unsorted, as fold indices can be
    expected = vec(LightGBM.predict(model, full[rows, :]))                 # one call on the whole-table rows
    for batch_rows in (1, 7, 4_096, 12_344, 12_345, 100_000)
        buf = Float32[]
        out = Pioneer.predict_rows!(Vector{Float64}(undef, length(rows)), X -> LightGBM.predict(model, X),
                                    buf, cols, rows; batch_rows = batch_rows)
        @test out == expected
    end
    # A row-identifying "model" catches any misplaced batch: score = the row's first feature.
    ident = Pioneer.predict_rows!(Vector{Float64}(undef, length(rows)), X -> Float64.(X[:, 1]),
                                  Float32[], cols, rows; batch_rows = 1_000)
    @test ident == Float64.(full[rows, 1])
    @test Pioneer.predict_rows!(Float64[], X -> error("not called"), Float32[], cols, Int[]) == Float64[]
    @test_throws DimensionMismatch Pioneer.predict_rows!(zeros(3), X -> zeros(size(X, 1)), Float32[], cols, [1, 2])
    @test_throws DimensionMismatch Pioneer.predict_rows!(zeros(2), X -> zeros(1), Float32[], cols, [1, 2])
end

@testset "MainSearch classifier scores equal whole-table predictions" begin
    df = _batched_test_table(60_000; seed = 5)
    Random.seed!(1844)
    scores, _, _, info = Pioneer.train_psm_classifier_with_fallback(
        df; features = _BATCHED_FEATS, buffers = Pioneer.LGBMMatrixBuffers())
    full = Pioneer.feature_matrix(df, info.available_features)
    for (fi, fold) in enumerate((0x00, 0x01))
        idx = findall(==(fold), df.cv_fold)
        p = info.predictor.fold_predictors[fi]
        @test p.kind == :lgbm
        @test scores[idx] == vec(LightGBM.predict(p.model, full[idx, :]))
    end
    @test Pioneer.predict_psm_classifier_scores(df, info.predictor; buffers = Pioneer.LGBMMatrixBuffers()) == scores
end
