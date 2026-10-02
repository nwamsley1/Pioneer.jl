# Decode every scan through the 88-byte index and compare centroids, RT, TIC, base peak, scan range with the .arrow.
using Mmap, Arrow, Tables
const RAW = expanduser("~/BrukerTims/rawformat/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw")
const ARW = expanduser("~/BrukerTims/tuning_ab/data/OlsenExploris500ng/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.arrow")
b = Mmap.mmap(RAW); t = Arrow.Table(ARW)
@inline rd(T, o) = ltoh(unsafe_load(Ptr{T}(pointer(b, o + 1))))
const IDX0 = 1320778831; const DATA0 = 110544; const N = length(t.msOrder)
bad = Dict(:index => 0, :rt => 0, :tic => 0, :bp => 0, :range => 0, :npk => 0, :mz => 0, :int => 0, :ptype => 0)
ptypes = Dict{Int32, Int}(); flagcount = Dict{UInt8, Int}()
for k in 0:N-1
    e = IDX0 + 88k; i = k + 1
    rd(Int32, e) == k || (bad[:index] += 1; continue)
    ptype = rd(Int32, e + 12); ptypes[ptype] = get(ptypes, ptype, 0) + 1
    ptype == t.packetType[i] || (bad[:ptype] += 1)
    Float32(rd(Float64, e + 20)) == t.retentionTime[i] || (bad[:rt] += 1)
    Float32(rd(Float64, e + 28)) == t.TIC[i] || (bad[:tic] += 1)
    (Float32(rd(Float64, e + 36)) == t.basePeakIntensity[i] && Float32(rd(Float64, e + 44)) == t.basePeakMz[i]) || (bad[:bp] += 1)
    (Float32(rd(Float64, e + 52)) == t.lowMz[i] && Float32(rd(Float64, e + 60)) == t.highMz[i]) || (bad[:range] += 1)
    p = DATA0 + Int(rd(UInt64, e + 68))
    prof = Int(rd(UInt32, p + 4)); pkw = Int(rd(UInt32, p + 8))
    q = p + 40 + 4prof
    n = pkw == 0 ? 0 : Int(rd(UInt32, q))
    d = q + 4pkw                                   # descriptor list: one UInt32 per peak
    keep = Int[]
    for j in 0:n-1
        f = UInt8((rd(UInt32, d + 4j) >> 16) & 0xff); flagcount[f] = get(flagcount, f, 0) + 1
        (f & 0x10) == 0 && push!(keep, j)
    end
    amz = collect(skipmissing(t.mz_array[i])); ait = collect(skipmissing(t.intensity_array[i]))
    length(keep) == length(amz) || (bad[:npk] += 1; continue)
    all(rd(Float32, q + 4 + 8keep[j]) == amz[j] for j in eachindex(keep)) || (bad[:mz] += 1)
    all(rd(Float32, q + 8 + 8keep[j]) == ait[j] for j in eachindex(keep)) || (bad[:int] += 1)
end
println("scans $N; mismatches: ", bad)
println("packet types: ", ptypes, "; descriptor flag bytes: ", sort(collect(flagcount)))
