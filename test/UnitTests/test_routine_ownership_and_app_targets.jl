using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const ROUTINE_SHIMS = Dict(
    joinpath(REPO_ROOT, "src", "Routines", "GenerateParams.jl") =>
        "packages\", \"PioneerParams\", \"src\", \"routines\", \"GenerateParams.jl",
    joinpath(REPO_ROOT, "src", "Routines", "BuildSpecLib.jl") =>
        "packages\", \"PioneerPredict\", \"src\", \"routines\", \"BuildSpecLib.jl",
    joinpath(REPO_ROOT, "src", "Routines", "SearchDIA.jl") =>
        "packages\", \"PioneerSearch\", \"src\", \"routines\", \"SearchDIA.jl"
)

const ROUTINE_IMPLEMENTATIONS = [
    joinpath(REPO_ROOT, "packages", "PioneerParams", "src", "routines", "GenerateParams.jl"),
    joinpath(REPO_ROOT, "packages", "PioneerPredict", "src", "routines", "BuildSpecLib.jl"),
    joinpath(REPO_ROOT, "packages", "PioneerSearch", "src", "routines", "SearchDIA.jl")
]

const WRAPPER_EXPECTATIONS = Dict(
    joinpath(REPO_ROOT, "src", "build", "CLI", "pioneer") => [
        "apps/params/bin",
        "apps/predict/bin",
        "apps/search/bin",
        "apps/convert/bin"
    ],
    joinpath(REPO_ROOT, "src", "build", "CLI", "pioneer.bat") => [
        "apps\\params\\bin",
        "apps\\predict\\bin",
        "apps\\search\\bin",
        "apps\\convert\\bin"
    ]
)

const WORKFLOW_EXPECTATIONS = Dict(
    joinpath(REPO_ROOT, ".github", "workflows", "build_app_linux.yml") => [
        "\"packages/PioneerParams\"",
        "\"packages/PioneerPredict\"",
        "\"packages/PioneerSearch\"",
        "\"packages/PioneerConvert\""
    ],
    joinpath(REPO_ROOT, ".github", "workflows", "build_app_windows.yml") => [
        "\"packages/PioneerParams\"",
        "\"packages/PioneerPredict\"",
        "\"packages/PioneerSearch\"",
        "\"packages/PioneerConvert\""
    ],
    joinpath(REPO_ROOT, ".github", "workflows", "build_app_macos.yml") => [
        "\"packages/PioneerParams\"",
        "\"packages/PioneerPredict\"",
        "\"packages/PioneerSearch\"",
        "\"packages/PioneerConvert\""
    ]
)

@testset "Routine implementations live in package routine folders" begin
    for routine_path in ROUTINE_IMPLEMENTATIONS
        @test isfile(routine_path)
    end

    for package_name in ("PioneerParams", "PioneerPredict", "PioneerSearch", "PioneerPlots")
        @test !isdir(joinpath(REPO_ROOT, "packages", package_name, "src", "owned"))
    end

    for (shim_path, routine_fragment) in ROUTINE_SHIMS
        @test isfile(shim_path)
        @test occursin(routine_fragment, read(shim_path, String))
    end
end

@testset "CLI wrapper dispatches to sub-app bundles" begin
    for (wrapper_path, fragments) in WRAPPER_EXPECTATIONS
        content = read(wrapper_path, String)
        for fragment in fragments
            @test occursin(fragment, content)
        end
    end
end

@testset "PackageCompiler workflows target split app packages" begin
    for (workflow_path, fragments) in WORKFLOW_EXPECTATIONS
        content = read(workflow_path, String)
        for fragment in fragments
            @test occursin(fragment, content)
        end
    end
end
