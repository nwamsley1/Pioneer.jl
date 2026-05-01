# Unit tests for DataFramePSMContainer and DataFramePSMContainerView
# (src/utils/ML/PSMScoring/psm_container.jl).
#
# Tests the core PSM data abstraction used by the entire scoring pipeline.
#
# Run standalone: julia --project=. test/UnitTests/test_psm_container.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
    using DataFrames
end

# Helper: build a minimal valid PSM DataFrame with required columns
function _make_psm_df(n::Int; n_folds::Int=2)
    DataFrame(
        target           = rand(Bool, n),
        cv_fold          = UInt8.(mod.(1:n, n_folds)),
        precursor_idx    = UInt32.(1:n),
        ms_file_idx      = UInt32.(ones(Int, n)),
        irt_pred         = rand(Float32, n),
        isotopes_captured = [(Int8(0), Int8(1)) for _ in 1:n],
        score            = rand(Float32, n),
    )
end

@testset "PSMContainer" begin

# ═══════════════════════════════════════════════════════════════════════════════
# DataFramePSMContainer — construction and validation
# ═══════════════════════════════════════════════════════════════════════════════

@testset "DataFramePSMContainer construction" begin
    @testset "valid DataFrame" begin
        df = _make_psm_df(10)
        c = Pioneer.DataFramePSMContainer(df)
        @test Pioneer.nrows(c) == 10
    end

    @testset "missing required column errors" begin
        df = DataFrame(target = [true], score = [0.5f0])
        @test_throws ErrorException Pioneer.DataFramePSMContainer(df)
    end

    @testset "unsafe constructor skips validation" begin
        df = DataFrame(x = [1, 2, 3])
        c = Pioneer.DataFramePSMContainer(df, Val(:unsafe))
        @test Pioneer.nrows(c) == 3
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Column operations
# ═══════════════════════════════════════════════════════════════════════════════

@testset "column operations" begin
    df = _make_psm_df(5)
    c = Pioneer.DataFramePSMContainer(df)

    @testset "get_column" begin
        col = Pioneer.get_column(c, :score)
        @test length(col) == 5
        @test col === df[!, :score]  # same underlying array
    end

    @testset "set_column!" begin
        new_data = Float32[10, 20, 30, 40, 50]
        Pioneer.set_column!(c, :new_col, new_data)
        @test Pioneer.has_column(c, :new_col)
        @test Pioneer.get_column(c, :new_col) == new_data
    end

    @testset "has_column" begin
        @test Pioneer.has_column(c, :target) == true
        @test Pioneer.has_column(c, :nonexistent) == false
    end

    @testset "getproperty shorthand" begin
        @test c.target === Pioneer.get_column(c, :target)
        @test c.data === df
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# CV fold operations
# ═══════════════════════════════════════════════════════════════════════════════

@testset "CV fold operations" begin
    df = _make_psm_df(10; n_folds=2)
    c = Pioneer.DataFramePSMContainer(df)

    @testset "get_cv_folds" begin
        folds = Pioneer.get_cv_folds(c)
        @test sort(folds) == UInt8[0, 1]
    end

    @testset "get_fold_indices" begin
        idx0 = Pioneer.get_fold_indices(c, UInt8(0))
        idx1 = Pioneer.get_fold_indices(c, UInt8(1))
        # Every row accounted for
        @test sort(vcat(idx0, idx1)) == collect(1:10)
        # Fold values match
        @test all(df.cv_fold[idx0] .== UInt8(0))
        @test all(df.cv_fold[idx1] .== UInt8(1))
    end

    @testset "get_train_indices" begin
        train = Pioneer.get_train_indices(c, UInt8(0))
        test = Pioneer.get_fold_indices(c, UInt8(0))
        # Training = everything except the test fold
        @test sort(vcat(train, test)) == collect(1:10)
        @test isempty(intersect(train, test))
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# View operations
# ═══════════════════════════════════════════════════════════════════════════════

@testset "DataFramePSMContainerView" begin
    df = _make_psm_df(10; n_folds=2)
    c = Pioneer.DataFramePSMContainer(df)
    view = Pioneer.get_view(c, [1, 3, 5, 7, 9])

    @testset "nrows" begin
        @test Pioneer.nrows(view) == 5
    end

    @testset "get_column returns view of parent" begin
        col = Pioneer.get_column(view, :precursor_idx)
        @test length(col) == 5
        @test col[1] == df.precursor_idx[1]
        @test col[3] == df.precursor_idx[5]
    end

    @testset "set_column! writes through to parent" begin
        Pioneer.set_column!(c, :test_col, zeros(Float32, 10))
        view_col = Float32[10, 20, 30, 40, 50]
        Pioneer.set_column!(view, :test_col, view_col)
        # Parent DataFrame should be modified at the view indices
        @test df.test_col[1] == 10.0f0
        @test df.test_col[3] == 20.0f0
        @test df.test_col[5] == 30.0f0
        # Non-view indices unchanged
        @test df.test_col[2] == 0.0f0
    end

    @testset "has_column delegates to parent" begin
        @test Pioneer.has_column(view, :target) == true
        @test Pioneer.has_column(view, :nonexistent) == false
    end

    @testset "nested view composes indices" begin
        inner = Pioneer.get_view(view, [1, 3])
        @test Pioneer.nrows(inner) == 2
        col = Pioneer.get_column(inner, :precursor_idx)
        @test col[1] == df.precursor_idx[1]
        @test col[2] == df.precursor_idx[5]
    end

    @testset "copy_container creates independent copy" begin
        copied = Pioneer.copy_container(view)
        @test Pioneer.nrows(copied) == 5
        Pioneer.set_column!(copied, :score, ones(Float32, 5))
        # Original unchanged
        @test Pioneer.get_column(view, :score) != ones(Float32, 5)
    end

    @testset "to_dataframe returns view" begin
        result = Pioneer.to_dataframe(view)
        @test size(result, 1) == 5
    end

    @testset "get_cv_folds on view" begin
        folds = Pioneer.get_cv_folds(view)
        @test all(f -> f in UInt8[0, 1], folds)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Utility functions
# ═══════════════════════════════════════════════════════════════════════════════

@testset "select_subset" begin
    df = _make_psm_df(10)
    c = Pioneer.DataFramePSMContainer(df)
    mask = BitVector([i <= 5 for i in 1:10])
    subset = Pioneer.select_subset(c, mask)
    @test Pioneer.nrows(subset) == 5
    @test Pioneer.get_column(subset, :precursor_idx) == UInt32.(1:5)
end

@testset "get_feature_matrix" begin
    df = _make_psm_df(5)
    df[!, :feat1] = Float32[1, 2, 3, 4, 5]
    df[!, :feat2] = Union{Missing, Float32}[10, missing, 30, 40, 50]
    c = Pioneer.DataFramePSMContainer(df)

    X = Pioneer.get_feature_matrix(c, [:feat1, :feat2])
    @test size(X) == (5, 2)
    @test X[:, 1] == Float32[1, 2, 3, 4, 5]
    @test X[2, 2] == 0.0f0  # missing → 0
    @test X[1, 2] == 10.0f0
end

@testset "get_labels" begin
    df = _make_psm_df(5)
    df.target .= [true, false, true, false, true]
    c = Pioneer.DataFramePSMContainer(df)
    labels = Pioneer.get_labels(c)
    @test labels == [true, false, true, false, true]
    @test eltype(labels) == Bool
end

@testset "sort_container!" begin
    df = _make_psm_df(5)
    df.precursor_idx .= UInt32[5, 3, 1, 4, 2]
    c = Pioneer.DataFramePSMContainer(df)
    Pioneer.sort_container!(c, [:precursor_idx])
    @test Pioneer.get_column(c, :precursor_idx) == UInt32[1, 2, 3, 4, 5]
end

@testset "eachindex" begin
    df = _make_psm_df(7)
    c = Pioneer.DataFramePSMContainer(df)
    @test collect(eachindex(c)) == 1:7
end

end # top-level testset
