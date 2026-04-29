using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

@testset "macOS codesigning filters bundled files" begin
    workflow = read(joinpath(REPO_ROOT, ".github", "workflows", "build_app_macos.yml"), String)

    @test occursin("is_macho_file()", workflow)
    @test occursin("file -b \"\$target\"", workflow)
    @test occursin("Skipping non-Mach-O file during codesign", workflow)
    @test occursin("find \"\$PKGROOT/usr/local/\$APP\" -type f -print0", workflow)
end
