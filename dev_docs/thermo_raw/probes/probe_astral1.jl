# Astral: packet types per MS order, and where a known MS2 scan's peaks sit inside its packet.
include(joinpath(@__DIR__, "ThermoRaw.jl")); using .ThermoRaw, Arrow, Tables
const D = expanduser("~/BrukerTims/rawformat/astral")
f = ThermoRaw.open_raw(joinpath(D, "20241211_bkc_25-0856_Goldfarb_Wamsley_Yeast_Alternating-v2_3min_Rep1.raw"))
t = Arrow.Table(joinpath(D, "20241211_bkc_25-0856_Goldfarb_Wamsley_Yeast_Alternating-v2_3min_Rep1.arrow"))
b = f.bytes; rd(T, o) = reinterpret(T, b[o+1:o+sizeof(T)])[1]
d = Dict{Tuple{Int32, UInt8}, Int}()
for s in f.scans; d[(s.packet_type, s.ms_order)] = get(d, (s.packet_type, s.ms_order), 0) + 1; end
println("(packet type, order) counts: ", d)
i = findfirst(k -> f.scans[k].ms_order == 2 && length(collect(skipmissing(t.mz_array[k]))) > 20, eachindex(f.scans))
s = f.scans[i]; mz = collect(skipmissing(t.mz_array[i])); it = collect(skipmissing(t.intensity_array[i]))
nxt = f.scans[i+1].packet_offset
println("scan $i: ptype $(s.packet_type), packet @$(s.packet_offset), next @$nxt (size $(nxt - s.packet_offset)); arrow npk $(length(mz)); hdr $(t.scanHeader[i])")
println("  first arrow peaks: ", collect(zip(mz[1:4], it[1:4])))
println("  header words: ", [rd(UInt32, s.packet_offset + 4j) for j in 0:15])
for (name, x) in (("mz1 F32", mz[1]), ("mz1 F64", Float64(mz[1])), ("int1 F32", it[1]), ("int1 F64", Float64(it[1])))
    pat = collect(reinterpret(UInt8, [x])); r = findnext(pat, b, s.packet_offset + 1)
    println("  $name at packet+", r === nothing || first(r) > nxt ? "none" : first(r) - 1 - s.packet_offset)
end
# nearest-value search for mz1 as Float64 within the packet (stored value may be higher precision)
best = (Inf, -1)
for o in s.packet_offset:nxt-8
    x = rd(Float64, o); abs(x - mz[1]) < abs(best[1] - mz[1]) && (best = (x, o - s.packet_offset))
end
println("  closest Float64 to mz1 in packet: $(best[1]) at +$(best[2])")
