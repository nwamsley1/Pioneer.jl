using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const JULIA_BIN = joinpath(Sys.BINDIR, Base.julia_exename())

function run_package_julia(project_path::String, script::String, args::Vector{String}=String[])
    cmd = Cmd(vcat([JULIA_BIN, "--startup-file=no", "--project=$project_path", "-e", script], args))
    run(cmd)
    return nothing
end

@testset "Split package library serialization stays interoperable" begin
    mktempdir() do tmpdir
        payload_path = joinpath(tmpdir, "detailed_fragments.jls")
        predict_project = joinpath(REPO_ROOT, "packages", "PioneerPredict")
        search_project = joinpath(REPO_ROOT, "packages", "PioneerSearch")

        write_script = """
        using PioneerPredict
        frags = PioneerPredict.DetailedFrag{Float32}[
            PioneerPredict.DetailedFrag(
                UInt32(7),
                Float32(500.25),
                Float16(1.5),
                UInt16(3),
                true,
                false,
                false,
                false,
                UInt8(1),
                UInt8(4),
                UInt8(2),
                UInt8(6),
                UInt8(0),
            ),
        ]
        PioneerPredict.serialize_library_detailed_frags(ARGS[1], frags)
        """

        read_script = """
        using PioneerSearch
        frags = PioneerSearch.load_detailed_frags(ARGS[1])
        @assert length(frags) == 1
        @assert frags[1] isa PioneerSearch.DetailedFrag{Float32}
        @assert PioneerSearch.getPID(frags[1]) == UInt32(7)
        @assert PioneerSearch.getMz(frags[1]) == Float32(500.25)
        @assert PioneerSearch.getIntensity(frags[1]) == Float16(1.5)
        """

        run_package_julia(predict_project, write_script, [payload_path])
        @test isfile(payload_path)
        run_package_julia(search_project, read_script, [payload_path])
        @test true
    end
end
