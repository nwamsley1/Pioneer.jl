using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const OWNED_ROUTINES = Dict(
    joinpath(REPO_ROOT, "src", "Routines", "GenerateParams.jl") =>
        "packages\", \"PioneerSearch\", \"src\", \"owned\", \"GenerateParams.jl",
    joinpath(REPO_ROOT, "src", "Routines", "BuildSpecLib.jl") =>
        "packages\", \"PioneerPredict\", \"src\", \"owned\", \"BuildSpecLib.jl",
    joinpath(REPO_ROOT, "src", "Routines", "SearchDIA.jl") =>
        "packages\", \"PioneerSearch\", \"src\", \"owned\", \"SearchDIA.jl"
)

const OWNED_IMPLEMENTATIONS = [
    joinpath(REPO_ROOT, "packages", "PioneerSearch", "src", "owned", "GenerateParams.jl"),
    joinpath(REPO_ROOT, "packages", "PioneerPredict", "src", "owned", "BuildSpecLib.jl"),
    joinpath(REPO_ROOT, "packages", "PioneerSearch", "src", "owned", "SearchDIA.jl")
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
        "\"apps/PioneerParamsApp\"",
        "\"apps/PioneerPredictApp\"",
        "\"apps/PioneerSearchApp\"",
        "\"apps/PioneerConvertApp\""
    ],
    joinpath(REPO_ROOT, ".github", "workflows", "build_app_windows.yml") => [
        "\"apps/PioneerParamsApp\"",
        "\"apps/PioneerPredictApp\"",
        "\"apps/PioneerSearchApp\"",
        "\"apps/PioneerConvertApp\""
    ],
    joinpath(REPO_ROOT, ".github", "workflows", "build_app_macos.yml") => [
        "\"apps/PioneerParamsApp\"",
        "\"apps/PioneerPredictApp\"",
        "\"apps/PioneerSearchApp\"",
        "\"apps/PioneerConvertApp\""
    ]
)

@testset "Routine ownership moved out of root" begin
    for owned_path in OWNED_IMPLEMENTATIONS
        @test isfile(owned_path)
    end

    for (shim_path, owned_fragment) in OWNED_ROUTINES
        @test isfile(shim_path)
        @test occursin(owned_fragment, read(shim_path, String))
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
