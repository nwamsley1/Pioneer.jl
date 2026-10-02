# Dump UInt64 / Int32 fields around the two places that store the data-stream start.
using Mmap
const RAW = expanduser("~/BrukerTims/rawformat/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw")
b = Mmap.mmap(RAW)
rd(T, o) = reinterpret(T, b[o+1:o+sizeof(T)])[1]
for site in (2850, 1320260968)
    println("== around $site (UInt64 at 8-byte steps relative to the hit)")
    for j in -12:12
        o = site + 8j
        v = rd(UInt64, o)
        tag = v == 110544 ? "  <- data start" : (1_000_000 < v < length(b) ? "  (plausible file offset)" : "")
        println("  @", o, " (", j >= 0 ? "+" : "", 8j, "): ", v, tag)
    end
end
println("candidates for index start: IDX0-12 = ", 1320778831 - 12)
