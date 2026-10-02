# Verify per-scan scan-event fields (MS order, isolation centre/width, collision energy) against the .arrow.
using Mmap, Arrow, Tables
const RAW = expanduser("~/BrukerTims/rawformat/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw")
const ARW = expanduser("~/BrukerTims/tuning_ab/data/OlsenExploris500ng/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.arrow")
b = Mmap.mmap(RAW); t = Arrow.Table(ARW)
@inline rd(T, o) = unsafe_load(Ptr{T}(pointer(b, o + 1)))
const EV0 = 1326252075; const N = length(t.msOrder); const REC = 272
bad = Dict(:order => 0, :center => 0, :width => 0, :ce => 0, :counts => 0)
for k in 0:N-1
    r = EV0 + 4 + REC * k; i = k + 1
    (rd(UInt32, r + 136) == 1 && rd(UInt32, r + 196) == 1) || (bad[:counts] += 1)
    order = b[r+7]
    order == t.msOrder[i] || (bad[:order] += 1)
    if order == 2
        Float32(rd(Float64, r + 140)) == t.centerMz[i] || (bad[:center] += 1)
        Float32(rd(Float64, r + 148)) == t.isolationWidthMz[i] || (bad[:width] += 1)
        Float32(rd(Float64, r + 156)) == t.collisionEnergyField[i] || (bad[:ce] += 1)
    else
        (ismissing(t.centerMz[i]) && ismissing(t.isolationWidthMz[i])) || (bad[:center] += 1)
    end
end
println("scans $N; scan-event mismatches: ", bad)
