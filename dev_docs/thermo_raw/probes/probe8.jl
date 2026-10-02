# Dump per-scan 272-byte scan-event records (scan 1 = MS1, scan 2 = MS2) as bytes and as Float64 at every 4-byte step.
using Mmap
const RAW = expanduser("~/BrukerTims/rawformat/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw")
b = Mmap.mmap(RAW)
@inline rd(T, o) = unsafe_load(Ptr{T}(pointer(b, o + 1)))
const EV0 = 1326252075
println("header UInt32: ", rd(UInt32, EV0))
for s in (1, 2, 3)
    r = EV0 + 4 + 272(s - 1)
    println("\n== scan $s record @$r")
    println("bytes 0..47: ", join([string(b[r+1+j], base = 16, pad = 2) for j in 0:47], " "))
    for o in 0:4:264
        x = rd(Float64, r + o); u = rd(UInt32, r + o)
        (isfinite(x) && 0.01 < abs(x) < 1e5) && println("  +$o Float64 = $x")
    end
    println("  nonzero UInt32 words: ", [(o, rd(UInt32, r + o)) for o in 0:4:268 if 0 < rd(UInt32, r + o) < 100000])
end
