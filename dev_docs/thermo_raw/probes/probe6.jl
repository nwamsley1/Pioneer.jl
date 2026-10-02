# Decode per-scan scan events: find scans 2-4's isolation centres (Float64) after the index and dump the records.
using Mmap, Arrow, Tables
const RAW = expanduser("~/BrukerTims/rawformat/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw")
const ARW = expanduser("~/BrukerTims/tuning_ab/data/OlsenExploris500ng/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.arrow")
b = Mmap.mmap(RAW); t = Arrow.Table(ARW)
rd(T, o) = reinterpret(T, b[o+1:o+sizeof(T)])[1]
const EV0 = 1326252075
for i in 1:5
    println("scan $i: order=$(t.msOrder[i]) center=$(t.centerMz[i]) width=$(t.isolationWidthMz[i]) ce=$(t.collisionEnergyField[i]) ceEv=$(t.collisionEnergyEvField[i]) hdr=$(t.scanHeader[i])")
end
# the filter string says 367.8500@hcd27.00; search for 367.85 as Float64 near EV0
function near(x, lo, hi)
    pat = collect(reinterpret(UInt8, [x])); r = findnext(pat, b, lo + 1)
    r === nothing || first(r) > hi ? nothing : first(r) - 1
end
c2 = near(367.85, EV0, EV0 + 10_000_000)
println("367.85 Float64 first hit after EV0: ", c2, c2 === nothing ? "" : " (+$(c2 - EV0))")
c3 = c2 === nothing ? nothing : near(Float64(395.8), c2 + 1, c2 + 100_000)
println("395.8 Float64 (scan 3 centre?) hit: ", c3, c3 === nothing || c2 === nothing ? "" : " (delta $(c3 - c2))")
println("first 200 bytes at EV0:")
for o in EV0:16:EV0+200
    println(o - EV0, ": ", join([string(b[o+1+j], base = 16, pad = 2) for j in 0:15], " "))
end
if c2 !== nothing
    println("\naround scan-2 centre (Float64 grid relative to hit):")
    for j in -8:10
        o = c2 + 8j; println("  ", 8j, ": ", rd(Float64, o), "   u32 ", rd(UInt32, o), " ", rd(UInt32, o + 4))
    end
end
