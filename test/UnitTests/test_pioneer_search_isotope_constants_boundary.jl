using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const JULIA_BIN = joinpath(Sys.BINDIR, Base.julia_exename())

function run_search_package_check(script::String)
    cmd = Cmd([JULIA_BIN, "--startup-file=no", "--project=$(joinpath(REPO_ROOT, "packages", "PioneerSearch"))", "-e", script])
    run(cmd)
    return nothing
end

@testset "PioneerSearch exposes isotope mass constants in split package runtime" begin
    script = """
    using PioneerSearch
    @assert isdefined(PioneerSearch, :NEUTRON)
    @assert isdefined(PioneerSearch, :PROTON)
    @assert isdefined(PioneerSearch, :H2O)
    @assert PioneerSearch.NEUTRON > 0
    @assert PioneerSearch.PROTON > 0
    @assert PioneerSearch.H2O > 0
    """
    run_search_package_check(script)
    @test true
end
