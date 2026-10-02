# Byte-step scan of the section after the index for values near scan 2/3's isolation centres (Float64 and Float32).
using Mmap
const RAW = expanduser("~/BrukerTims/rawformat/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw")
b = Mmap.mmap(RAW)
@inline rd(T, o) = unsafe_load(Ptr{T}(pointer(b, o + 1)))
function scan(lo, hi, target; tol = 1e-3, limit = 6)
    h64 = Int[]; h32 = Int[]
    for o in lo:hi-8
        x = rd(Float64, o); (abs(x - target) < tol && length(h64) < limit) && push!(h64, o)
        y = rd(Float32, o); (abs(y - target) < tol && length(h32) < limit) && push!(h32, o)
        length(h64) >= limit && length(h32) >= limit && break
    end
    h64, h32
end
EV0 = 1326252075; EV1 = 1343169391
for (s, c) in ((2, 367.85), (3, 381.55), (4, 395.25))
    h64, h32 = scan(EV0, EV1 + 5_000_000, c)
    println("scan $s centre $c: Float64 @", h64, "  Float32 @", h32)
end
h64, _ = scan(EV0, EV1, 15.4; tol = 1e-4)
println("15.4 (CE eV of scan 2) Float64 @", h64)
