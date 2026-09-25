using Test, DataFrames, Arrow, CSV

function protein_export_fixture(dir)
    library_path = joinpath(dir, "library.arrow")
    Arrow.write(library_path, DataFrame(accession=["A","B","D"],
        gene_name=["GA","","GD"], protein_name=["Alpha","Beta","Delta"]))
    proteins = Pioneer.LibraryProteins(Arrow.Table(library_path))
    rows = DataFrame(
        file_name=["r2","r1","r3","r1","r2","r3"],
        species=["human","human","yeast","human","human",""],
        protein=["A","A","B","C","C","D"], target=fill(true, 6), entrap_id=zeros(UInt8, 6),
        peptides=Vector{Union{Missing, UInt32}}[[1,2], [1], [1,2,3], [2,3], [2], []],
        n_precursors=UInt32[2,1,3,2,1,0], n_precursors_quantified=UInt32[2,1,3,2,0,0],
        n_modified_peptides=UInt32[2,1,3,2,1,0], n_peptides=UInt32[2,1,3,2,1,0],
        global_pg_score=Float32[10,10,5,20,20,40], pg_score=fill(8f0, 6),
        global_qval=fill(0.001f0, 6), qval=fill(0.002f0, 6), pg_pep=fill(0.003f0, 6),
        total_peak_area=Union{Missing,Float32}[20,10,30,40,missing,missing],
        abundance=Union{Missing,Float32}[20,10,30,40,missing,missing])
    path = joinpath(dir, "protein_groups_long.arrow")
    open(Arrow.Writer, path; file=true, ntasks=0) do writer
        for span in (1:1, 2:4, 5:6)
            Arrow.write(writer, rows[span, :])
        end
    end
    return path, proteins, rows
end

@testset "Bounded protein export" begin
    mktempdir() do dir
        path, proteins, rows = protein_export_fixture(dir)
        args = (["AAA","BBB","CCC"], Union{Missing,String}[missing,missing,missing], Union{Missing,String}[missing,missing,missing], fill(UInt8(2),3),
                ["r1","r2","r3","r4"], proteins)
        stats = [Pioneer.RunSummaryStats(name) for name in args[5]]
        Pioneer.add_protein_group_counts!(stats, path, args[5])
        @test getfield.(stats, :protein_groups_identified) == [2,2,2,0]
        @test getfield.(stats, :protein_groups_quantified) == [2,1,1,0]
        original = read(path)
        wide_path = Pioneer.writeProteinGroupsCSV(path, args...; memory_budget_bytes=4096, batch_size=1)
        wide = DataFrame(Arrow.Table(wide_path))
        @test wide.protein == ["B","C","A","D"]
        @test wide.gene_names == ["B","C","GA","GD"]
        @test wide.protein_names == ["Beta","C","Alpha","Delta"]
        @test isequal(wide.r1, Union{Missing,Float32}[missing,40,10,missing])
        @test isequal(wide.r2, Union{Missing,Float32}[missing,missing,20,missing])
        @test isequal(wide.r3, Union{Missing,Float32}[30,missing,missing,missing])
        @test all(ismissing, wide.r4)
        @test all(batch -> length(batch.protein) == 1, Arrow.Stream(wide_path))
        long = CSV.read(joinpath(dir, "protein_groups_long.tsv"), DataFrame; delim='\t')
        @test long.protein == ["B","C","C","A","A","D"]
        @test long.file_name == ["r3","r1","r2","r2","r1","r3"]
        @test long.peptides[1] == "_AAA_.2;_BBB_.2;_CCC_.2"
        @test ismissing(long.peptides[end])
        @test isequal(select(DataFrame(Arrow.Table(path)), names(rows)), rows)
        @test read(path) == original
        @test !any(startswith(".protein_export_"), readdir(dir))

        # Arrow-only export cannot dereference precursor annotations.
        Pioneer.writeProteinGroupsCSV(path, String[], Union{Missing,String}[],
            Union{Missing,String}[], UInt8[], args[5], proteins; write_csv=false, memory_budget_bytes=4096)
        @test isequal(DataFrame(Arrow.Table(wide_path)), wide)
        @test !isfile(joinpath(dir, "protein_groups_long.tsv"))
        @test !isfile(joinpath(dir, "protein_groups_wide.tsv"))

        # Failed exports preserve existing output and remove scratch files.
        previous = read(wide_path)
        bad = copy(rows[1:2, :]); bad.file_name .= "r1"
        Arrow.write(path, bad)
        @test_throws ArgumentError Pioneer.writeProteinGroupsCSV(path, args...)
        @test read(wide_path) == previous
        @test !any(startswith(".protein_export_"), readdir(dir))

        # Empty inputs still produce readable, consistently typed output.
        Arrow.write(path, rows[1:0, :])
        Pioneer.writeProteinGroupsCSV(path, args...; write_csv=false)
        @test nrow(DataFrame(Arrow.Table(wide_path))) == 0
        @test propertynames(Arrow.Table(wide_path)) == propertynames(wide)
    end
end

@testset "Protein export buffer sizing" begin
    df = DataFrame(peptides=[UInt32[1,2] for _ in 1:10])
    @test collect(Pioneer._protein_export_ranges(df, 550, 100)) == [1:2,3:4,5:6,7:8,9:10]
    @test collect(Pioneer._protein_export_ranges(df, 550, 1)) == [i:i for i in 1:10]
    @test collect(Pioneer._protein_export_ranges(df, 550, 100;
        sequence_lengths=Dict(UInt32(1)=>500, UInt32(2)=>500))) == [i:i for i in 1:10]
end

@testset "Protein export group identity" begin
    mktempdir() do dir
        path, proteins, rows = protein_export_fixture(dir)
        target = rows[1:2, :]
        decoy = copy(target)
        decoy.target .= false
        decoy.entrap_id .= 1
        decoy.global_pg_score .= 30
        Arrow.write(path, vcat(target, decoy))
        output = Pioneer.writeProteinGroupsCSV(path, String[], Union{Missing,String}[],
            Union{Missing,String}[], UInt8[], ["r1","r2"], proteins; write_csv=false)
        result = DataFrame(Arrow.Table(output))
        @test result.protein == ["A","A"]
        @test result.target == [false,true]
        @test result.entrap_id == [1,0]
        @test result.r1 == Float32[10,10]
        # Noncontiguous input groups are rejected rather than silently duplicated.
        Arrow.write(path, vcat(target[1:1,:], decoy, target[2:2,:]))
        @test_throws ArgumentError Pioneer.writeProteinGroupsCSV(path, String[],
            Union{Missing,String}[], Union{Missing,String}[], UInt8[], ["r1","r2"], proteins;
            write_csv=false)
        @test !any(startswith(".protein_export_"), readdir(dir))
    end
end
