# A staged MBR integration input is a row selection of a scored table and its sidecars. Loading it
# must give exactly what the old staging copy gave: the selected rows, every column, in order.
using Test, Arrow, DataFrames
import Pioneer

@testset "staged MBR row selection" begin
    mktempdir() do dir
        n = 50
        source_path = joinpath(dir, "run1.arrow")
        Arrow.write(source_path, DataFrame(
            precursor_idx = UInt32.(101:100+n),
            rt = Float32.(1:n),
            note = Union{Missing, String}[isodd(i) ? "n$i" : missing for i in 1:n],
        ))
        source = Pioneer.PSMFileReference(source_path)
        Pioneer.add_columns_via_sidecar!(source, :qval => Float32.(1:n) ./ 100; tag = "scores")
        rows = [2, 3, 17, 40, 50]

        staged_path = joinpath(dir, "passing", "run1.arrow")
        mkpath(dirname(staged_path))
        Pioneer.write_staged_selection(staged_path, source, rows)
        @test !isfile(staged_path)

        expected = Pioneer.load_with_sidecars(source)[rows, :]
        @test isequal(Pioneer.load_staged_psms(staged_path), expected)
        @test isequal(Pioneer.load_staged_psms(staged_path, [:precursor_idx, :qval, :absent]),
                      expected[:, [:precursor_idx, :qval]])

        # Once integration writes the real table, it is read and the selection is gone.
        Pioneer.writeArrow(staged_path, expected[1:2, :])
        Pioneer.clear_staged_selection!(staged_path)
        @test !isfile(staged_path * Pioneer.MBR_SELECTION_SUFFIX)
        @test isequal(Pioneer.load_staged_psms(staged_path), expected[1:2, :])
        @test Pioneer.load_staged_psms(staged_path, [:rt]).rt == expected.rt[1:2]

        # An empty selection loads as an empty table with every column.
        empty_path = joinpath(dir, "passing", "run2.arrow")
        Pioneer.write_staged_selection(empty_path, source, Int[])
        empty = Pioneer.load_staged_psms(empty_path)
        @test nrow(empty) == 0 && names(empty) == names(expected)
    end
end
