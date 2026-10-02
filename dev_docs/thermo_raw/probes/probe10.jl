# Tail of the 272-byte scan-event record: count at +216 and what follows, for an MS1 and an MS2 scan.
using Mmap
const RAW = expanduser("~/BrukerTims/rawformat/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw")
b = Mmap.mmap(RAW)
@inline rd(T, o) = unsafe_load(Ptr{T}(pointer(b, o + 1)))
const EV0 = 1326252075
for s in (1, 2)
    r = EV0 + 4 + 272(s - 1)
    println("scan $s: +188 u32 ", rd(UInt32, r + 188), " ", rd(UInt32, r + 192), "; +216 count ", rd(UInt32, r + 216))
    println("  Float64 from +220: ", [rd(Float64, r + 220 + 8j) for j in 0:5])
    println("  UInt32 from +220: ", [rd(UInt32, r + 220 + 4j) for j in 0:12])
end
