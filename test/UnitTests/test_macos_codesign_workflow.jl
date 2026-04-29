using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

@testset "macOS codesigning normalizes bundled file permissions" begin
    workflow = read(joinpath(REPO_ROOT, ".github", "workflows", "build_app_macos.yml"), String)

    @test occursin("chmod -R u+rwX \"\$PKGROOT/usr/local/\$APP\"", workflow)
    @test occursin("find \"\$PKGROOT/usr/local/\$APP\" -type f -print0", workflow)
    @test occursin("codesign_with_retry \"\$f\" || exit 1", workflow)
    @test !occursin("is_macho_file()", workflow)
end
