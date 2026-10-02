# How build_run_similarity scales with the number of runs, on synthetic observations.
# Each precursor ID has a per-run detection probability p (a mix of common, mid and rare IDs, roughly
# like a cohort); every run observes each ID independently with that probability.
# Usage: julia --project=. scripts/profiling/run_similarity_scaling.jl [n_ids] [N...]
using Pioneer, Random, Printf

function synthetic_observations(rng, n_ids, n_runs)
    p = similar(zeros(Float64), n_ids)
    for i in 1:n_ids
        u = rand(rng)
        p[i] = u < 0.4 ? 0.9 + 0.1rand(rng) : u < 0.7 ? 0.3 + 0.6rand(rng) : 0.3rand(rng)
    end
    obs = Dict{UInt32, Vector{UInt32}}()
    for run in UInt32(1):UInt32(n_runs)
        obs[run] = UInt32[UInt32(i) for i in 1:n_ids if rand(rng) < p[i]]
    end
    return obs
end

# The work the atlas loop does: per ID with document frequency d among N runs, min(d(d-1)/2, d(N-d)) Dict updates.
function predicted_updates(obs, n_runs)
    df = Dict{UInt32, Int}()
    for ids in values(obs), id in ids
        df[id] = get(df, id, 0) + 1
    end
    return sum(d -> d == n_runs ? 0 : min(d * (d - 1) ÷ 2, d * (n_runs - d)), values(df); init = 0)
end

n_ids = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 300_000
Ns = length(ARGS) >= 2 ? parse.(Int, ARGS[2:end]) : [25, 50, 100, 200, 400]
rng = MersenneTwister(7)
Pioneer.build_run_similarity(synthetic_observations(rng, 1000, 5))      # compile
@printf("%6s %12s %10s %12s %10s %14s\n", "runs", "postings", "seconds", "updates", "ns/update", "atlas_MB")
for N in Ns
    obs = synthetic_observations(rng, n_ids, N)
    postings = sum(length, values(obs))
    upd = predicted_updates(obs, N)
    GC.gc()
    t = @elapsed atlas = Pioneer.build_run_similarity(obs)
    @printf("%6d %12d %10.2f %12d %10.1f %14.1f\n", N, postings, t, upd, 1e9t / max(upd, 1),
            Base.summarysize(atlas) / 2^20)
    flush(stdout)
end
