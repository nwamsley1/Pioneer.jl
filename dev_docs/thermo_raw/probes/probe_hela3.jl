# HELA: locate the first scan-event records by known Float64 values (scan 1-3 mass ranges 300/1100) and dump them.
using Mmap, Arrow, Tables
const RAW = expanduser("~/Projects/pioneer-gui-test/rawtest/HELA_MOUSESIL_METHTEST_07312024_06.raw")
const ARW = expanduser("~/Projects/pioneer-gui-test/rawtest/arrow_out/HELA_MOUSESIL_METHTEST_07312024_06.arrow")
b = Mmap.mmap(RAW); t = Arrow.Table(ARW)
rd(T, o) = reinterpret(T, b[o+1:o+sizeof(T)])[1]
ev = 520189768
function hits(x, lo, hi; limit = 8)
    pat = collect(reinterpret(UInt8, [x])); out = Int[]; i = lo + 1
    while length(out) < limit
        r = findnext(pat, b, i); (r === nothing || first(r) > hi) && break
        push!(out, first(r) - 1); i = first(r) + 1
    end
    out
end
h300 = hits(300.0, ev, ev + 5000); h1100 = hits(1100.0, ev, ev + 5000)
println("300.0 at +", h300 .- ev, "\n1100.0 at +", h1100 .- ev)
println("record bytes from ev+4:")
for o in 0:16:543
    println("  +$o: ", join([string(b[ev+5+o+j], base = 16, pad = 2) for j in 0:15], " "))
end
for i in 1:3; println("arrow scan $i: ", t.scanHeader[i], " low/high ", t.lowMz[i], "/", t.highMz[i]); end
