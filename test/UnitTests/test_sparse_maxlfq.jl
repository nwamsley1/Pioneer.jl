using Test, Random
using Pioneer: solve_sparse_maxlfq, solve_maxlfq, _sparse_lfq_pairs

@testset "Sparse MaxLFQ" begin
    rng = MersenneTwister(18)
    for nr in (0, 1, 2, 6, 20), np in (1, 8, 25), seed in 1:3
        X = Matrix{Union{Missing,Float64}}(20 .+ randn(rng, np, nr))
        X[rand(rng, np, nr) .< 0.3] .= missing
        priorities = Union{Missing,Float64}[rand(rng) for _ in 1:nr]
        expected, labels = solve_maxlfq(X, priorities)
        actual = solve_sparse_maxlfq(X, priorities; partners=max(0,nr-1), seed)
        @test actual.component_labels == labels
        @test ismissing.(actual.estimates) == ismissing.(expected)
        @test all(ismissing(a) || isapprox(a,b; atol=4e-6,rtol=0) for (a,b) in zip(actual.estimates,expected))
        for k in (0, 2)
            sparse = solve_sparse_maxlfq(X, priorities; partners=k, seed)
            @test sparse.component_labels == labels
            @test ismissing.(sparse.estimates) == ismissing.(expected)
            @test sparse.edge_count <= (k+1)*nr
        end
    end
    # A long chain, disconnected components, and priority-based tie breaking.
    X = Matrix{Union{Missing,Float64}}(missing, 9, 10)
    for p in 1:9
        X[p,p], X[p,p+1] = p, p+1
    end
    chain = solve_sparse_maxlfq(X, zeros(10); partners=0)
    @test chain.edge_count == 9
    @test diff(Float64.(chain.estimates)) ≈ ones(9) atol=2e-6
    @test_throws ErrorException solve_sparse_maxlfq(X, zeros(10); partners=0, maxiter=1)
    d = Union{Missing,Float64}[1 2 missing missing; missing missing 4 5]
    @test ismissing.(solve_sparse_maxlfq(d, [0,0,2,2]).estimates) == [true,true,false,false]
    full = fill(10.0, 3, 40)
    @test solve_sparse_maxlfq(full, zeros(40); partners=2).iterations == 0
    e1,_,_ = _sparse_lfq_pairs(full,2,7)
    e2,_,_ = _sparse_lfq_pairs(full .+ randn(rng,3,40),2,7)
    @test e1 == e2
    @test solve_sparse_maxlfq(full, zeros(40); partners=2,seed=3) ==
          solve_sparse_maxlfq(full, zeros(40); partners=2,seed=3)
    extreme = Union{Missing,Float64}[1100 1101; 1101 1102]
    @test solve_sparse_maxlfq(extreme, zeros(2)).estimates ≈
        first(solve_maxlfq(extreme, Union{Missing,Float64}[0,0]))
    @test all(ismissing, solve_sparse_maxlfq(fill(missing,0,5),zeros(5)).estimates)
    @test_throws ArgumentError solve_sparse_maxlfq(full, zeros(40); partners=-1)
    @test_throws DimensionMismatch solve_sparse_maxlfq(full, zeros(2))
    @test_throws ArgumentError solve_sparse_maxlfq([NaN 2.], zeros(2))
end
