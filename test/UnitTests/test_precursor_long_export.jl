@testset "Precursor Arrow export across growing string dictionaries" begin
    cases = ((1_000, 2, true, false, Int16), (33_000, 128, false, false, Int32),
             (8, 2, false, true, Int8))
    for (nrows, first_chunk_rows, normalized, all_missing_mods, index_type) in cases
        @testset "$nrows distinct file names" begin
            mktempdir() do dir
                file_names = ["run_$i" for i in 1:nrows]
                rows = DataFrame(
                    precursor_idx = UInt32.(1:nrows),
                    ms_file_idx = UInt32.(1:nrows),
                    file_name = file_names,
                    species = [isodd(i) ? "human" : "mouse" for i in 1:nrows],
                    structural_mods = Union{Missing, String}[
                        all_missing_mods || i <= 2 || i == 500 ? missing : "mod_$i" for i in 1:nrows],
                    sequence = ["PEPTIDE_$i" for i in 1:nrows],
                    target = [i != 3 for i in 1:nrows],
                    peak_area = Float32[i % 5 == 0 ? 0 : i for i in 1:nrows],
                    peak_area_normalized = Float32[i % 5 == 0 ? 0 : 2i for i in 1:nrows],
                    irt_error = fill(Float16(1), nrows),
                    rt_fwhm = fill(Float16(0.1), nrows),
                    points_integrated = fill(UInt32(5), nrows),
                    charge = fill(UInt8(2), nrows),
                    missed_cleavage = fill(UInt8(0), nrows),
                    isotopic_mods = fill("", nrows),
                )
                refs = Pioneer.PSMFileReference[]
                for (chunk, indices) in enumerate((1:first_chunk_rows, first_chunk_rows+1:nrows))
                    path = joinpath(dir, "chunk_$chunk.arrow")
                    midpoint = first(indices) + length(indices) ÷ 2 - 1
                    # Merged MaxLFQ chunks contain multiple Arrow record batches.
                    open(Arrow.Writer, path; file=true) do writer
                        Arrow.write(writer, rows[first(indices):midpoint, :])
                        Arrow.write(writer, rows[midpoint+1:last(indices), :])
                    end
                    @test Arrow.Table(path).file_name isa SentinelArrays.ChainedVector
                    push!(refs, Pioneer.PSMFileReference(path))
                end

                policy = Pioneer.OutputSchemaPolicy(Dict("isotope_mod_groups" => []))
                output = joinpath(dir, "precursors_long.arrow")
                stats = Pioneer.write_precursor_long_arrow(output, refs, file_names, policy;
                    run_to_run_normalization = normalized)
                for batch in Arrow.Stream(output)
                    for (name, expected_type) in ((:file_name, index_type), (:species, Int8),
                            (:structural_mods, all_missing_mods ? Int8 : index_type))
                        column = Tables.getcolumn(batch, name)
                        @test column isa Arrow.DictEncoded
                        @test eltype(column.indices) == expected_type
                    end
                end
                actual = DataFrame(Arrow.Table(output))
                expected = select(rows, Not(:isotopic_mods))
                for name in (:peak_area, :peak_area_normalized)
                    expected[!, name] = Union{Missing, Float32}[v > 0 ? v : missing for v in rows[!, name]]
                end
                normalized || select!(expected, Not(:peak_area_normalized))
                @test names(actual) == names(expected)
                @test nrow(actual) == nrows
                for name in names(expected)
                    @test eltype(actual[!, name]) == eltype(expected[!, name])
                    @test isequal(actual[!, name], expected[!, name])
                end
                @test getfield.(stats, :file_name) == file_names
                @test sum(s.precursors_identified for s in stats) == nrows - 1
                @test sum(s.precursors_quantified for s in stats) == nrows - nrows ÷ 5 - 1
                @test stats[3].precursors_identified == 0
                @test stats[1].total_peak_area == 1.0
                @test stats[1].normalization_factors == Float32[2]

                if nrows == 1_000
                    repeated_refs = Pioneer.PSMFileReference[]
                    for (chunk, ref) in enumerate(refs)
                        path = joinpath(dir, "repeated_$chunk.arrow")
                        repeated = repeat(DataFrame(Arrow.Table(Pioneer.file_path(ref))); outer=10)
                        Arrow.write(path, repeated)
                        push!(repeated_refs, Pioneer.PSMFileReference(path))
                    end
                    encoded_path = joinpath(dir, "repeated_encoded.arrow")
                    Pioneer.write_precursor_long_arrow(encoded_path, repeated_refs, file_names, policy;
                        run_to_run_normalization = normalized)
                    plain_path = joinpath(dir, "repeated_plain.arrow")
                    plain_columns = map(collect, Tables.columntable(Arrow.Table(encoded_path)))
                    Arrow.write(plain_path, plain_columns)
                    @test filesize(encoded_path) < filesize(plain_path)
                end
            end
        end
    end
end
