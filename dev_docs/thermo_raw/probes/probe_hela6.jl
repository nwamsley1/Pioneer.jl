# HELA peak-count mismatches: which packet types / orders, and how our count compares to the .arrow.
include(joinpath(@__DIR__, "ThermoRaw.jl")); using .ThermoRaw, Arrow, Tables
f = ThermoRaw.open_raw(expanduser("~/Projects/pioneer-gui-test/rawtest/HELA_MOUSESIL_METHTEST_07312024_06.raw"))
t = Arrow.Table(expanduser("~/Projects/pioneer-gui-test/rawtest/arrow_out/HELA_MOUSESIL_METHTEST_07312024_06.arrow"))
cnt = Dict{Tuple{Int32, UInt8, Symbol}, Int}(); ex = Dict{Int32, Any}()
for (i, s) in enumerate(f.scans)
    mz, _ = ThermoRaw.centroids(f, i); a = length(collect(skipmissing(t.mz_array[i])))
    key = (s.packet_type, s.ms_order, length(mz) == a ? :match : (a == 0 ? :arrow_empty : :differ))
    cnt[key] = get(cnt, key, 0) + 1
    key[3] != :match && !haskey(ex, s.packet_type) && (ex[s.packet_type] = (i, length(mz), a, t.scanHeader[i]))
end
println(sort(collect(cnt))); println("examples: ", ex)
