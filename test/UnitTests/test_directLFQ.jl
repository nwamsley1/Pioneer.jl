using Test, Statistics, JSON, Random
if !@isdefined(Pioneer)
    using Pioneer
end

@testset "directLFQ reference parity" begin
    cases = JSON.parsefile(joinpath(@__DIR__, "..", "fixtures", "directlfq", "reference.json"))
    for case in cases
        @testset "$(case["name"])" begin
            X = [isnothing(x) ? missing : Float64(x) for row in case["input"] for x in row]
            X = permutedims(reshape(X, length(case["input"][1]), length(case["input"])))
            result, selected, counts = Pioneer.solve_directlfq(X)
            @test selected == case["selected"]
            @test counts == case["counts"]
            @test ismissing.(result) == isnothing.(case["expected"])
            for i in eachindex(result)
                if !ismissing(result[i])
                    @test result[i] ≈ case["expected"][i] atol=3e-5 rtol=1e-6
                end
            end
            # Stable precursor identities make input row ordering irrelevant.
            order = randperm(MersenneTwister(44), size(X, 1))
            permuted, _, _ = Pioneer.solve_directlfq(X[order, :]; precursor_ids=order)
            @test isequal(permuted, result)
        end
    end
    @test all(ismissing, first(Pioneer.solve_directlfq(fill(missing, 3, 4))))
    @test isempty(first(Pioneer.solve_directlfq(zeros(3, 0))))
    @test_throws DimensionMismatch Pioneer.solve_directlfq(zeros(3, 4); precursor_ids=[1])
    # Numerically stable rescaling beyond the linear Float64 intensity range.
    @test first(Pioneer.solve_directlfq([1100.0 1101; 1101 1102])) ≈
          Float32[1101.5849625, 1102.5849625]
end

@testset "directLFQ selects traces before building matrices" begin
    rng = MersenneTwister(16)
    for nprec in (1, 11, 101, 180)
        areas = exp2.(20 .+ randn(rng, nprec, 21))
        areas[rand(rng, size(areas)...) .< 0.25] .= NaN
        peptides = repeat(UInt32.(1:nprec); outer=21)
        experiments = repeat(UInt32.(1:21); inner=nprec)
        expected = first(Pioneer.solve_directlfq(log2.(areas)))
        actual, counts = Pioneer.directlfq_from_observations(peptides, experiments,
            vec(areas), Dict(UInt32(i)=>i for i in 1:21))
        @test isequal(actual, expected)
        @test counts == last(Pioneer.solve_directlfq(log2.(areas)))
    end
    @test_throws ArgumentError Pioneer.directlfq_from_observations(UInt32[1,1],
        UInt32[1,1], [10.,11.], Dict(UInt32(1)=>1))
end

@testset "directLFQ invalid areas" begin
    result, counts = Pioneer.directlfq_from_observations(UInt32[1,1,1,1,1],
        UInt32[1,2,3,4,5], Union{Missing, Float64}[0,-1,missing,Inf,NaN],
        Dict(UInt32(i)=>i for i in 1:5))
    @test all(ismissing, result)
    @test all(iszero, counts)
end
