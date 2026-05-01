# Unit tests for ArrowFilePSMContainer and ArrowFileGroup.
#
# These provide file-backed PSM storage for out-of-memory scoring.
# ArrowFileGroup wraps a single (ms_file_idx, cv_fold) Arrow file pair;
# ArrowFilePSMContainer groups multiple file groups under the
# AbstractPSMContainer interface.
#
# Run standalone: julia --project=. test/UnitTests/test_arrow_psm_container.jl
# Run via suite:  julia --project=. test/runtests.jl

if !@isdefined(Pioneer)
    using Test
    using Pioneer
    using Arrow
end

# ═══════════════════════════════════════════════════════════════════════════════
# Helper: write a minimal Arrow file with the given columns
# ═══════════════════════════════════════════════════════════════════════════════

function _write_test_arrow(path::String, n::Int; extra_cols::Dict=Dict{Symbol,Vector}())
    cols = Dict{Symbol,Vector}(
        :score      => rand(Float32, n),
        :target     => rand(Bool, n),
        :scan_idx   => UInt32.(1:n),
    )
    merge!(cols, extra_cols)
    Arrow.write(path, (; (k => v for (k, v) in cols)...))
end

# ═══════════════════════════════════════════════════════════════════════════════
# _derive_scores_path
# ═══════════════════════════════════════════════════════════════════════════════

@testset "ArrowPSMContainer" begin

@testset "_derive_scores_path" begin
    derive = Pioneer._derive_scores_path

    # Standard fold-suffixed paths
    @test derive("/tmp/file_fold0.arrow") == "/tmp/file_fold0_scores.arrow"
    @test derive("/tmp/file_fold1.arrow") == "/tmp/file_fold1_scores.arrow"
    @test derive("/tmp/data_fold12.arrow") == "/tmp/data_fold12_scores.arrow"

    # Fallback: no _foldN suffix → insert _scores before .arrow
    @test derive("/tmp/plain_data.arrow") == "/tmp/plain_data_scores.arrow"
end

# ═══════════════════════════════════════════════════════════════════════════════
# get_fold_number
# ═══════════════════════════════════════════════════════════════════════════════

@testset "get_fold_number" begin
    mktempdir() do dir
        path0 = joinpath(dir, "psms_fold0.arrow")
        path1 = joinpath(dir, "psms_fold1.arrow")
        _write_test_arrow(path0, 5)
        _write_test_arrow(path1, 3)

        g0 = Pioneer.ArrowFileGroup(path0)
        g1 = Pioneer.ArrowFileGroup(path1)

        @test Pioneer.get_fold_number(g0) == UInt8(0)
        @test Pioneer.get_fold_number(g1) == UInt8(1)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# ArrowFileGroup constructor
# ═══════════════════════════════════════════════════════════════════════════════

@testset "ArrowFileGroup" begin
    mktempdir() do dir
        path = joinpath(dir, "test_fold0.arrow")
        _write_test_arrow(path, 42)

        g = Pioneer.ArrowFileGroup(path)
        @test g.data_path == path
        @test g.scores_path == joinpath(dir, "test_fold0_scores.arrow")
        @test g.n_rows == 42
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# ArrowFilePSMContainer — constructor and interface methods
# ═══════════════════════════════════════════════════════════════════════════════

@testset "ArrowFilePSMContainer" begin
    mktempdir() do dir
        # Create 4 files: 2 ms_files × 2 folds
        paths = String[]
        for ms in ["fileA", "fileB"], fold in [0, 1]
            p = joinpath(dir, "$(ms)_fold$(fold).arrow")
            n = fold == 0 ? 10 : 20
            _write_test_arrow(p, n)
            push!(paths, p)
        end

        container = Pioneer.ArrowFilePSMContainer(paths; max_training_psms=500)

        # ── nrows ──────────────────────────────────────────────────────────
        @testset "nrows" begin
            # 2 files × fold0(10) + 2 files × fold1(20) = 60
            @test Pioneer.nrows(container) == 60
        end

        # ── get_cv_folds ───────────────────────────────────────────────────
        @testset "get_cv_folds" begin
            folds = Pioneer.get_cv_folds(container)
            @test folds == UInt8[0, 1]
        end

        # ── get_file_groups_for_fold ───────────────────────────────────────
        @testset "get_file_groups_for_fold" begin
            groups0 = Pioneer.get_file_groups_for_fold(container, UInt8(0))
            groups1 = Pioneer.get_file_groups_for_fold(container, UInt8(1))
            @test length(groups0) == 2
            @test length(groups1) == 2
            @test all(g -> g.n_rows == 10, groups0)
            @test all(g -> g.n_rows == 20, groups1)
        end

        # ── has_column ─────────────────────────────────────────────────────
        @testset "has_column" begin
            @test Pioneer.has_column(container, :score) == true
            @test Pioneer.has_column(container, :target) == true
            @test Pioneer.has_column(container, :nonexistent_col) == false
        end

        # ── max_training_psms stored correctly ─────────────────────────────
        @testset "max_training_psms" begin
            @test container.max_training_psms == 500
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Error stubs — operations requiring in-memory data should error
# ═══════════════════════════════════════════════════════════════════════════════

@testset "error stubs" begin
    mktempdir() do dir
        path = joinpath(dir, "stub_fold0.arrow")
        _write_test_arrow(path, 5)
        container = Pioneer.ArrowFilePSMContainer([path])

        @test_throws ErrorException Pioneer.get_column(container, :score)
        @test_throws ErrorException Pioneer.set_column!(container, :score, Float32[1.0])
        @test_throws ErrorException Pioneer.to_dataframe(container)
        @test_throws ErrorException Pioneer.sort_container!(container, [:score])
        @test_throws ErrorException Pioneer.get_view(container, [1, 2])
        @test_throws ErrorException Pioneer.copy_container(container)
        @test_throws ErrorException Pioneer.get_fold_indices(container, UInt8(0))
        @test_throws ErrorException Pioneer.get_train_indices(container, UInt8(0))
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Edge cases
# ═══════════════════════════════════════════════════════════════════════════════

@testset "edge cases" begin
    # ── Empty container — sum over empty groups errors ────────────────────
    @testset "empty container" begin
        @test_throws ArgumentError Pioneer.ArrowFilePSMContainer(String[])
    end

    # ── Single file ────────────────────────────────────────────────────────
    @testset "single file" begin
        mktempdir() do dir
            path = joinpath(dir, "only_fold0.arrow")
            _write_test_arrow(path, 7)
            container = Pioneer.ArrowFilePSMContainer([path])
            @test Pioneer.nrows(container) == 7
            @test Pioneer.get_cv_folds(container) == UInt8[0]
            @test length(Pioneer.get_file_groups_for_fold(container, UInt8(0))) == 1
            @test isempty(Pioneer.get_file_groups_for_fold(container, UInt8(1)))
        end
    end

    # ── has_column with extra columns ──────────────────────────────────────
    @testset "has_column with extra columns" begin
        mktempdir() do dir
            path = joinpath(dir, "extra_fold0.arrow")
            _write_test_arrow(path, 3; extra_cols=Dict(:my_feature => Float32[1.0, 2.0, 3.0]))
            container = Pioneer.ArrowFilePSMContainer([path])
            @test Pioneer.has_column(container, :my_feature) == true
            @test Pioneer.has_column(container, :score) == true
            @test Pioneer.has_column(container, :missing_col) == false
        end
    end

    # ── get_fold_number parse error ────────────────────────────────────────
    @testset "get_fold_number bad path" begin
        # Manually construct a group with a path that has no fold suffix
        mktempdir() do dir
            path = joinpath(dir, "no_fold.arrow")
            _write_test_arrow(path, 2)
            g = Pioneer.ArrowFileGroup(path)
            @test_throws ErrorException Pioneer.get_fold_number(g)
        end
    end
end

end # top-level testset
