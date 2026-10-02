# Map RunHeader fields relative to its start (pointed to by the FileHeader at byte 2866).
using Mmap
const RAW = expanduser("~/BrukerTims/rawformat/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw")
b = Mmap.mmap(RAW)
@inline rd(T, o) = unsafe_load(Ptr{T}(pointer(b, o + 1)))
rh = Int(rd(UInt64, 2866))
println("FileHeader UInt64 @2850 (data addr) = $(rd(UInt64, 2850)); @2866 (run header addr) = $rh")
# first / last scan numbers (1, 62196) as Int32, and low/high time, within the first 8 KB of the run header
for o in 0:4:8191
    v = rd(Int32, rh + o)
    v == 62196 && println("Int32 62196 at RunHeader+$o; neighbours: ", [rd(Int32, rh + o + 4j) for j in -3:3])
end
for (name, v) in (("index", 1320778827), ("data", 110544), ("stream A", 1320272278), ("stream B", 1320759706),
                  ("scan events", 1326252075), ("trailer?", 1343169391), ("run header self", rh))
    for o in 0:8:8191
        rd(UInt64, rh + o) == v && println(rpad(name, 16), " UInt64 at RunHeader+$o")
    end
    for o in 0:4:8191
        rd(UInt32, rh + o) == v && v < typemax(UInt32) && println(rpad(name, 16), " UInt32 at RunHeader+$o")
    end
end
println("Float64 at RunHeader+ 0..40: ", [rd(Float64, rh + o) for o in 0:8:40])
