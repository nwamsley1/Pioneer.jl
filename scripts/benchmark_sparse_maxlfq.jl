using Pioneer, Random, JSON
function sparse_scaling_benchmark()
    rng = MersenneTwister(63)
    warm = Matrix{Union{Missing,Float64}}(20 .+ randn(rng,8,100))
    Pioneer.solve_sparse_maxlfq(warm,zeros(100);partners=16)
    results=Any[]
    for n in (100,500,1000,3000,6000)
        X=Matrix{Union{Missing,Float64}}(20 .+ randn(rng,8,n))
        X[rand(rng,8,n).<0.2].=missing
        GC.gc()
        timed=@timed Pioneer.solve_sparse_maxlfq(X,zeros(n);partners=16,seed=17)
        r=timed.value
        push!(results,(;runs=n,precursors=8,seconds=timed.time,allocated_bytes=timed.bytes,
            edges=r.edge_count,iterations=r.iterations,residual_norm=r.residual_norm))
        println(last(results));flush(stdout)
    end
    isempty(ARGS) || write(only(ARGS),JSON.json(results,2))
end
sparse_scaling_benchmark()
