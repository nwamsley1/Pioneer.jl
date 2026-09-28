# ZT multi-chunk main search: the k-way merge of sorted chunk files into precursor-complete files.

using Test, Arrow, DataFrames, Random
using Pioneer: zt_merge_by_precursor, _check_precursor_partitions

@testset "ZT merge into precursor-complete files" begin
    rng = MersenneTwister(7)
    dir = mktempdir()
    # 5 chunks (cycle ranges); precursors elute across neighbouring chunks, so most are split
    paths = String[]; all_rows = DataFrame()
    for c in 1:5
        pid = UInt32.(rand(rng, 1:200, 400)); scn = UInt32.(c * 1000 .+ rand(rng, 1:900, 400))
        df = unique(DataFrame(precursor_idx = pid, scan_idx = scn))
        sort!(df, [:precursor_idx, :scan_idx])
        df[!, :w] = rand(rng, Float32, nrow(df)); df[!, :t] = rand(rng, Bool, nrow(df))
        p = joinpath(dir, "chunk_$c.arrow"); Arrow.write(p, df); push!(paths, p)
        append!(all_rows, df)
    end
    sort!(all_rows, [:precursor_idx, :scan_idx])
    for target in (1, 97, 500, 10_000)
        out = zt_merge_by_precursor(paths, joinpath(dir, "merged_$target"); target_rows = target, batch_rows = 64)
        merged = reduce(vcat, (DataFrame(Arrow.Table(p)) for p in out))
        @test merged == all_rows                               # same rows, same order, same values
        owner = Dict{UInt32, Int}()
        split = false
        for (i, p) in enumerate(out), pid in Arrow.Table(p).precursor_idx
            get!(owner, pid, i) == i || (split = true)
        end
        @test !split                                           # each precursor in exactly one file
        target == 10_000 && @test length(out) == 1
        target == 1 && @test length(out) == length(unique(all_rows.precursor_idx))
    end
    # the invariant check rejects a precursor split across two files
    a = joinpath(dir, "a.arrow"); b = joinpath(dir, "b.arrow")
    Arrow.write(a, (precursor_idx = UInt32[1, 2], scan_idx = UInt32[1, 1]))
    Arrow.write(b, (precursor_idx = UInt32[2, 3], scan_idx = UInt32[2, 1]))
    @test_throws ErrorException _check_precursor_partitions([a, b], 4)
end
