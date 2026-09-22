# Does a precursor's mobility apex move between runs? Per precursor and run: the weight-weighted IM centroid,
# the argmax IM slice, and the best PSM's IM scan (the window's centre).
using Arrow, DataFrames, Statistics, Printf
const STRIDE = 8
const H = expanduser("~/BrukerTims/pride_hye")
runs = [("A", r) for r in 1:6]; append!(runs, [("B", r) for r in 1:6])
summ = Dict{String, Dict{UInt32, NTuple{3, Float64}}}()   # run => pid => (centroid, argmax, best-PSM IM)
for (C, R) in runs
    f = "h50_$(C)_REP$R"
    d = DataFrame(Arrow.Table(joinpath(H, "dump_h50_cond$(C)_wide", f * "_chrom_weights.arrow")))
    p = DataFrame(Arrow.Table(joinpath(H, "out_h50_cond$(C)_wide", "temp_data", "passing_psms", f * ".arrow")))
    p = p[(p.target) .& (p.qval .<= 0.01), :]
    bs = Dict{UInt32, Int}(r.precursor_idx => Int(r.scan_idx) for r in eachrow(p))
    scan_im = Dict{Int, Int}(); for r in eachrow(d); scan_im[Int(r.scan_idx)] = Int(r.im_scan); end
    m = Dict{UInt32, NTuple{3, Float64}}()
    for sub in groupby(d, :precursor_idx)
        w = Float64.(sub.weight); s = sum(w); s > 0 || continue
        ims = Float64.(sub.im_scan)
        pid = sub.precursor_idx[1]
        m[pid] = (sum(w .* ims) / s, ims[argmax(w)], Float64(get(scan_im, get(bs, pid, -1), NaN)))
    end
    summ[f] = m
    @printf("%s: %d precursors\n", f, length(m))
end
common = intersect([Set(keys(summ["h50_$(C)_REP$R"])) for (C, R) in runs]...)
@printf("\nprecursors in all 6 runs: %d\n", length(common))
cen_sd = Float64[]; arg_sd = Float64[]; psm_sd = Float64[]; cen_minus_psm = Float64[]; ab = Float64[]
for pid in common
    cs = [summ["h50_$(C)_REP$R"][pid][1] for (C, R) in runs]
    as = [summ["h50_$(C)_REP$R"][pid][2] for (C, R) in runs]
    ps = [summ["h50_$(C)_REP$R"][pid][3] for (C, R) in runs]
    push!(cen_sd, std(cs)); push!(arg_sd, std(as)); any(isnan, ps) || push!(psm_sd, std(ps))
    any(isnan, ps) || push!(cen_minus_psm, median(cs .- ps))
    push!(ab, mean(cs[1:3]) - mean(cs[4:6]))
end
@printf("across the 6 runs, per precursor (IM scans):\n")
@printf("  SD of the weighted IM centroid : median %.1f, p90 %.1f, p99 %.1f\n", median(cen_sd), quantile(cen_sd, .9), quantile(cen_sd, .99))
@printf("  SD of the argmax IM slice      : median %.1f, p90 %.1f\n", median(arg_sd), quantile(arg_sd, .9))
@printf("  SD of the best-PSM IM scan     : median %.1f, p90 %.1f\n", median(psm_sd), quantile(psm_sd, .9))
@printf("  centroid - best-PSM IM (offset): median %+.1f, IQR %.1f .. %.1f, |offset| > 16 scans in %.1f%%\n",
        median(cen_minus_psm), quantile(cen_minus_psm, .25), quantile(cen_minus_psm, .75), 100mean(abs.(cen_minus_psm) .> 16))
@printf("  condition A minus B centroid   : median %+.1f, p90 %.1f (a systematic shift would show here)\n", median(ab), quantile(abs.(ab), .9))
