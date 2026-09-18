# Run several centroiding conversions of one .d in ONE Julia session (JIT once).
# Usage: julia -t 12 --project=proto convert_many.jl <run.d> <out_dir> "S M K Q C [sum|avg] [minscans] [ms1k=K] [ms1q=Q] [split]" ...
#   S im sigma, M mz sigma, K stride, Q cull quantile, C wmean|gauss; sum = sum-scaled intensities; a bare integer = min scans;
#   ms1k / ms1q = MS1 stride / cull quantile (default: the MS2 values); split = per-level cull thresholds.
include(joinpath(@__DIR__, "tdf_centroid_to_arrow.jl"))

function main()
    dpath, out_dir = ARGS[1], ARGS[2]
    for spec in ARGS[3:end]
        toks = split(spec); S, M, K, Q, C = toks[1:5]
        sum_scale = false; min_scans = 1; ms1k = parse(Int, K); ms1q = parse(Float64, Q); split_cull = false
        for t in toks[6:end]
            if t == "sum";               sum_scale = true
            elseif t == "avg";           sum_scale = false
            elseif t == "split";         split_cull = true
            elseif startswith(t, "ms1k="); ms1k = parse(Int, t[6:end])
            elseif startswith(t, "ms1q="); ms1q = parse(Float64, t[6:end])
            else;                        min_scans = parse(Int, t)
            end
        end
        p = CParams(parse(Float64, S), parse(Float64, M), parse(Int, K), parse(Float64, Q), Symbol(C), 4, sum_scale, min_scans)
        t0 = time()
        println("\n=== convert $spec ===")
        convert_centroid(dpath, out_dir, p; ms1_stride = ms1k, ms1_cull_q = ms1q, split_cull = split_cull)
        println("=== done $spec in $(round(time() - t0, digits = 1)) s ===")
        GC.gc()
    end
end
main()
