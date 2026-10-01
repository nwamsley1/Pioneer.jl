# Is fragment m/z across adjacent Q1 bins iid jitter or drift? Track raw centroids along streaks.
using Arrow, Statistics, Printf
const t = Arrow.Table(ARGS[1]); const CYC = parse(Int, ARGS[2])
ms2 = [i for i in 1:length(t.msOrder) if t.msOrder[i] == 2 && t.cycle_idx[i] == CYC]
nb = length(ms2)
# seeds: peaks in bin j with intensity > 2000; follow ±6 bins within 15 ppm, pick nearest peak per bin
lag1 = Float64[]; span_ppm = Float64[]; step_ppm = Float64[]; slope_ppm_per_bin = Float64[]; nlen = Int[]
for j in 7:12:nb-7
    s = ms2[j]; mz = t.mz_array[s]; it = t.intensity_array[s]
    for (m0, x0) in zip(mz, it)
        (ismissing(m0) || ismissing(x0) || x0 < 2000) && continue
        devs = Float64[]
        for d in -6:6
            s2 = ms2[j+d]; best = NaN; bi = 0.0
            for (m, x) in zip(t.mz_array[s2], t.intensity_array[s2])
                (ismissing(m) || ismissing(x)) && continue
                abs(m - m0) / m0 * 1e6 < 15 && x > bi && (best = (m - m0) / m0 * 1e6; bi = x)
            end
            push!(devs, best)
        end
        ok = .!isnan.(devs); count(ok) >= 9 || continue
        v = devs[ok]; v .-= mean(v)
        n = length(v); push!(nlen, n)
        push!(lag1, sum(v[1:end-1] .* v[2:end]) / sum(v .^ 2))
        push!(span_ppm, maximum(v) - minimum(v))
        push!(step_ppm, median(abs.(diff(v))))
        xs = collect(1:n) .- (n+1)/2; push!(slope_ppm_per_bin, sum(xs .* v) / sum(xs .^ 2))
    end
end
q(v) = join([@sprintf("%.2f", quantile(v, p)) for p in (0.1, 0.5, 0.9)], " / ")
println("streaks: ", length(lag1))
println("lag-1 autocorrelation (p10/p50/p90): ", q(lag1), "   (iid jitter ~0; drift >0)")
println("bin-to-bin |step| ppm median:        ", q(step_ppm))
println("span (max-min) ppm over streak:      ", q(span_ppm), "   grid step 5 mDa = ", @sprintf("%.1f", 5e-3/600*1e6), " ppm at 600")
println("|linear slope| ppm/bin:              ", q(abs.(slope_ppm_per_bin)))
println("frac of streaks with span > 8 ppm (one 5 mDa row at 600): ", @sprintf("%.2f", mean(span_ppm .> 8)))
