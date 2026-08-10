# Regression coverage for the portable application build.
#
# LoopVectorization decides its vector width from the build host before
# PackageCompiler creates target-specific clones. That made AVX-512 code from
# the build machine leak into the portable Windows sysimage. Keep this test at
# the start of a fresh test subprocess so the loaded-module check is meaningful.

const _PHASE1_FORBIDDEN_CPU_PACKAGES = Set([
    "HostCPUFeatures",
    "LinearSolve",
    "LoopVectorization",
    "RecursiveFactorization",
    "TriangularSolve",
    "VectorizationBase",
])

@testset "Portable CPU dependency boundary" begin
    loaded_packages = Set(pkgid.name for pkgid in keys(Base.loaded_modules))
    @test isempty(intersect(loaded_packages, _PHASE1_FORBIDDEN_CPU_PACKAGES))

    project_text = read(joinpath(@__DIR__, "..", "..", "Project.toml"), String)
    for package in ("HostCPUFeatures", "LinearSolve", "LoopVectorization")
        @test isnothing(match(Regex("(?m)^" * package * raw"\s*="), project_text))
    end

    source_root = joinpath(@__DIR__, "..", "..", "src")
    source_text = join(
        (read(joinpath(root, file), String)
         for (root, _, files) in walkdir(source_root)
         for file in files if endswith(file, ".jl")),
        '\n',
    )
    @test !occursin("@turbo", source_text)
end

@testset "Target-neutral loop replacements" begin
    n = 257
    feature = collect(range(-1.0, 1.0; length=n))
    X = DataFrame(intercept=ones(Float64, n), feature=feature)
    X_matrix = Matrix(X)
    W = collect(range(0.5, 1.5; length=n))
    Z = @. 0.25 + 0.75 * feature

    XWX = zeros(Float64, 2, 2)
    Pioneer.fillXWX!(XWX, X, W)
    @test XWX ≈ X_matrix' * (reshape(W, :, 1) .* X_matrix)

    Y = zeros(Float64, 2)
    Pioneer.fillY!(Y, X, W, Z)
    @test Y ≈ X_matrix' * (W .* Z)

    chunks = Iterators.partition(1:n, 31)
    β = Float64[0.25, -0.75]
    η = fill(99.0, n)
    Pioneer.fillη!(η, X, β, (-8.0, 8.0), chunks)
    expected_η = clamp.(X_matrix * β, -8.0, 8.0)
    @test η ≈ expected_η
    @test Pioneer.vecSum!(η, chunks) ≈ sum(expected_η)

    scores = fill(Float32(99), n)
    Pioneer.ModelPredictProbs!(scores, X, β, chunks)
    expected_scores = @. (1 + erf(expected_η / sqrt(2))) / 2
    @test scores ≈ Float32.(expected_scores)

    H = Pioneer.SparseArrayFused(UInt32(4))
    H.n_vals = 2
    H.m = 2
    H.n = 1
    H.rowval[1:2] .= UInt32[1, 2]
    H.nzval[1:2] .= Float32[1, 3]
    H.x[1:2] .= Float32[0.5, 1]
    H.colptr[1:2] .= UInt32[1, 3]

    residuals = fill(Float32(99), 4)
    Pioneer.initResiduals!(residuals, H, Float32[2])
    @test residuals == Float32[1.5, 5, 99, 99]

    residuals .= Float32(99)
    Pioneer._accumulate_residuals!(Float32[2], residuals, H)
    @test residuals == Float32[1.5, 5, 99, 99]
end
