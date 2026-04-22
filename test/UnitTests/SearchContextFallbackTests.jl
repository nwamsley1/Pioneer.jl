@testset "Search context fallback models" begin
    struct DummySpecLib <: Pioneer.SpectralLibrary end
    struct DummyMsRef <: Pioneer.MassSpecDataReference end
    struct DummySearchData <: Pioneer.SearchDataStructures end

    ctx = Pioneer.SearchContext(DummySpecLib(), [DummySearchData()], DummyMsRef(), 1, 0, 1)

    @test Pioneer.getIrtError(ctx, 2) == typemax(Float32)
    @test Pioneer.getIrtErrors(ctx)[2] == typemax(Float32)
end
