using Test, Pioneer, DataFrames

@testset "DataFrame block store append after read" begin
    mktempdir() do dir
        @test_throws ArgumentError Pioneer.DataFrameBlockStore(joinpath(dir, "invalid"), 0)
        store = Pioneer.DataFrameBlockStore(joinpath(dir, "blocks.bin"), 1)
        try
            blocks = [DataFrame(x=fill(i, i)) for i in 1:4]
            for block in blocks[1:3]
                Pioneer._store_dataframe_block!(store, block)
            end
            @test Pioneer._dataframe_block(store, 1; cache=true) == blocks[1]
            @test Pioneer._dataframe_block(store, 1; cache=true) === store.cached_block
            Pioneer._store_dataframe_block!(store, blocks[4])
            flush(store)
            @test store.row_ends == [1, 3, 6, 10]
            @test isempty(store.blocks)
            @test readdir(dir) == ["blocks.bin"]
            for i in [4, 2, 3, 1]
                @test Pioneer._dataframe_block(store, i) == blocks[i]
            end
        finally
            close(store)
        end
    end
end

@testset "block byte estimate (cache budget) tracks summarysize" begin
    n = 1_000
    numeric = DataFrame(a = rand(Float32, n), b = rand(UInt32, n), c = allowmissing(rand(Float64, n)))
    # isbits / isbits-union columns are counted exactly (plus the union selector byte)
    @test Pioneer._approx_block_bytes(numeric) == 4n + 4n + (8 + 1)n
    # A protein-export-like block: strings, a missing-able string, and a per-row vector of peptide ids
    protein = DataFrame(
        protein = ["P$(i);Q$(i)" for i in 1:n],
        file_name = fill("20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E10H50Y40_30SPD_DIA_1", n),
        species = Union{Missing, String}[isodd(i) ? "HUMAN" : missing for i in 1:n],
        abundance = rand(Float32, n),
        peptides = [Union{Missing, UInt32}[isodd(j) ? UInt32(j) : missing for j in 1:20] for _ in 1:n],
    )
    est, exact = Pioneer._approx_block_bytes(protein), Base.summarysize(protein)
    @test 0.5exact <= est <= 2exact
    @test Pioneer._approx_block_bytes(DataFrame()) == 0
end
