using Test
using Pioneer

@testset "MainSearch pair competition removal" begin
    @test !isdefined(Pioneer, :apply_pair_competition!)
    @test !isdefined(Pioneer, :_build_paircomp_lookups)
    @test !isdefined(Pioneer, :_mark_paircomp_losers!)
end
