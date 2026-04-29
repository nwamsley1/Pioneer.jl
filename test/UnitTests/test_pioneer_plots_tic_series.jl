using Test
using Arrow
using DataFrames

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PIONEER_SEARCH_PROJECT = joinpath(REPO_ROOT, "packages", "PioneerSearch")

pushfirst!(LOAD_PATH, PIONEER_SEARCH_PROJECT)
try
    using PioneerPlots
finally
    popfirst!(LOAD_PATH)
end

@testset "PioneerPlots TIC reader handles Arrow record batches" begin
    mktempdir() do tmp_dir
        ms_table_path = joinpath(tmp_dir, "ms_table.arrow")
        open(Arrow.Writer, ms_table_path) do writer
            Arrow.write(writer, DataFrame(
                msOrder = UInt8[0x01, 0x02],
                retentionTime = Float32[1.0, 2.0],
                TIC = Float32[10.0, 20.0],
            ))
            Arrow.write(writer, DataFrame(
                msOrder = UInt8[0x01, 0x01],
                retentionTime = Float32[3.0, 4.0],
                TIC = Float32[30.0, 40.0],
            ))
        end

        retention_times, tics = PioneerPlots._read_tic_series(ms_table_path)
        @test retention_times == Float64[1.0, 3.0, 4.0]
        @test tics == Float64[10.0, 30.0, 40.0]
    end
end
