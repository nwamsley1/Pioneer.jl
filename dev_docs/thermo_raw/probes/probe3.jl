# Find scan 1's packet start (data stream start) and read its centroid list vs the .arrow.
using Mmap, Arrow, Tables
const RAW = expanduser("~/BrukerTims/rawformat/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw")
const ARW = expanduser("~/BrukerTims/tuning_ab/data/OlsenExploris500ng/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.arrow")
b = Mmap.mmap(RAW); t = Arrow.Table(ARW)
rd(T, o) = reinterpret(T, b[o+1:o+sizeof(T)])[1]
pat = vcat(reinterpret(UInt8, [350.0f0]), reinterpret(UInt8, [1400.0f0]))
lo = 111456 - 1888 - 64
r = findnext(pat, b, lo + 1)
println("Float32 350,1400 pair at ", r === nothing ? "none" : first(r) - 1)
if r !== nothing
    h = first(r) - 1 - 32      # assume low/high are the last 8 bytes of a 40-byte header
    println("header candidate @$h words: ", [rd(UInt32, h + 4j) for j in 0:9])
    w = [rd(UInt32, h + 4j) for j in 0:7]
    println("  sum of words 1..6 *4 + 40 = ", 40 + 4 * sum(w[2:7]), " (packet size 1888)")
end
n = rd(UInt32, 111452)
pairs = [(rd(Float32, 111456 + 8k), rd(Float32, 111460 + 8k)) for k in 0:n-1]
mz = collect(skipmissing(t.mz_array[1])); it = collect(skipmissing(t.intensity_array[1]))
println("count field $n; arrow $(length(mz))")
println("raw:   ", pairs[1:min(end, 30)])
println("arrow: ", collect(zip(mz, it)))
println("after peak list @", 111456 + 8n, ": ", [rd(UInt32, 111456 + 8n + 4j) for j in 0:15])
