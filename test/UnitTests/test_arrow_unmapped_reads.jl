using Test
using Pioneer
using Arrow
using DataFrames
using Tables

# Julia's rm refuses a file that is still mapped on Windows, so deleting a file right after a
# read checks that the read unmapped it. On other systems the delete always succeeds.
deleted_now(path) = (rm(path); !isfile(path))

function unmapped_reads_test_file(path)
    Arrow.write(path, (
        id = UInt32.(1:6),
        score = Union{Missing, Float32}[0.5, 0.25, missing, 1.0, 0.0, 0.75],
        species = ["HUMAN", "YEAST", "HUMAN", "ECOLI", "YEAST", "HUMAN"],
        peptides = [UInt32[1, 2], UInt32[], UInt32[3], UInt32[4, 5, 6], UInt32[7], UInt32[8, 9]],
        maybe_list = [Union{Missing, UInt32}[1, missing], missing, Union{Missing, UInt32}[2],
                      Union{Missing, UInt32}[], missing, Union{Missing, UInt32}[3, 4]],
    ); dictencode = true)
    return path
end

@testset "Arrow reads that unmap their file" begin
    mktempdir() do dir
        path = unmapped_reads_test_file(joinpath(dir, "table.arrow"))
        reference = DataFrame(Tables.columntable(Arrow.Table(read(path))))   # nothing mapped

        df = Pioneer.load_arrow_dataframe(path)
        @test names(df) == names(reference)
        for c in names(reference)
            @test isequal(collect(df[!, c]), collect(reference[!, c]))
        end
        # Same column types as a plain copying read, except list columns, whose cells are copied.
        @test typeof(df.species) == typeof(reference.species)
        @test typeof(df.score) == typeof(reference.score)
        @test df.peptides isa Vector{Vector{UInt32}}
        @test df.maybe_list isa Vector{Union{Missing, Vector{Union{Missing, UInt32}}}}
        @test deleted_now(path)
        @test df.peptides[4] == UInt32[4, 5, 6]          # still readable after the file is gone

        path = unmapped_reads_test_file(joinpath(dir, "cols.arrow"))
        part = Pioneer.load_arrow_dataframe(path; cols = [:peptides, :id, :not_a_column])
        @test names(part) == ["peptides", "id"]
        @test part.id == UInt32.(1:6)
        @test deleted_now(path)

        # A column that would still point into the file is an error, not a dangling read.
        nested = joinpath(dir, "nested.arrow")
        Arrow.write(nested, (x = [[UInt32[1], UInt32[2, 3]], [UInt32[4]]],))
        @test_throws ErrorException Pioneer.load_arrow_dataframe(nested)
        @test deleted_now(nested)

        a = unmapped_reads_test_file(joinpath(dir, "a.arrow"))
        b = unmapped_reads_test_file(joinpath(dir, "b.arrow"))
        total = Pioneer.with_arrow_tables() do open_table
            sum(open_table(a).id) + sum(open_table(b).id)
        end
        @test total == 42
        @test deleted_now(a) && deleted_now(b)

        s = unmapped_reads_test_file(joinpath(dir, "stream.arrow"))
        n = Pioneer.with_arrow_stream(s) do stream
            sum(length(batch.id) for batch in stream)
        end
        @test n == 6
        @test deleted_now(s)

        # File references read their schema and row count without leaving the file mapped.
        r = unmapped_reads_test_file(joinpath(dir, "ref.arrow"))
        ref = Pioneer.PSMFileReference(r)
        @test ref.row_count == 6
        @test deleted_now(r)
    end
end
