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
