using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PREDICT_BOOTSTRAP = joinpath(REPO_ROOT, "packages", "PioneerPredict", "src", "bootstrap.jl")

@testset "PioneerPredict bootstrap loads mod types before Koina FASTA structs" begin
    bootstrap = read(PREDICT_BOOTSTRAP, String)

    mods_idx = findfirst("safe_include!(joinpath(root_path, \"structs\", \"mods.jl\"))", bootstrap)
    koina_idx = findfirst("safe_include_directory!(joinpath(REPO_ROOT, \"src\", \"structs\", \"KoinaStructs\"))", bootstrap)

    @test mods_idx !== nothing
    @test koina_idx !== nothing
    @test mods_idx < koina_idx
end
