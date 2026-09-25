# Count passing precursors (targets / decoys at q <= 0.01, by charge) in a run's temp_data/passing_psms.
# Usage: julia --project=<Pioneer checkout> count_passing.jl <run_dir> [<run_dir> ...]
using Arrow, DataFrames, Printf
for run_dir in ARGS
    d = joinpath(run_dir, "temp_data", "passing_psms")
    files = filter(f -> endswith(f, ".arrow") && !occursin("sidecar", f), readdir(d; join = true))
    for f in files
        t = DataFrame(Arrow.Table(f))
        ok = t[t.qval .<= 0.01, :]
        println(basename(run_dir), " / ", basename(f), ": rows ", nrow(t), ", q<=0.01 targets ", count(ok.target), " decoys ", count(.!ok.target))
        for z in sort(unique(ok.charge))
            s = ok[ok.charge .== z, :]
            println(@sprintf("    z=%d  targets %6d  decoys %4d", z, count(s.target), count(.!s.target)))
        end
        t01 = ok[ok.target .& (ok.qval .<= 0.001), :]
        println("    q<=0.001 targets ", nrow(t01))
    end
end
