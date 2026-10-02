# HELA (Orbitrap Eclipse/Lumos with ion-trap scans): run header, index entries, and scan-event record structure.
using Mmap, Arrow, Tables
const RAW = expanduser("~/Projects/pioneer-gui-test/rawtest/HELA_MOUSESIL_METHTEST_07312024_06.raw")
const ARW = expanduser("~/Projects/pioneer-gui-test/rawtest/arrow_out/HELA_MOUSESIL_METHTEST_07312024_06.arrow")
b = Mmap.mmap(RAW); t = Arrow.Table(ARW)
rd(T, o) = reinterpret(T, b[o+1:o+sizeof(T)])[1]          # bounds-checked
println("version ", rd(UInt32, 0x24), "; size ", length(b))
rh = Int(rd(UInt64, 2866)); println("run header @$rh; data addr (file hdr) ", rd(UInt64, 2850))
println("first/last scan ", rd(Int32, rh + 8), " ", rd(Int32, rh + 12), "; arrow scans ", length(t.msOrder))
addrs = [(o, rd(UInt64, rh + o)) for o in 7400:8:7480]; println("RunHeader UInt64s @7400..: ", addrs)
idx = Int(rd(UInt64, rh + 7408)); ev = Int(rd(UInt64, rh + 7448)); tr = Int(rd(UInt64, rh + 7456))
n = rd(Int32, rh + 12) - rd(Int32, rh + 8) + 1
println("index entry 0: idx=", rd(Int32, idx + 4), " ptype=", rd(Int32, idx + 16), " rt=", rd(Float64, idx + 24))
println("events section bytes ", tr - ev, "; per scan if uniform ", (tr - ev - 4) / n)
for s in 1:3
    println("arrow scan $s: order $(t.msOrder[s]) hdr $(t.scanHeader[s])")
end
r = ev + 4
println("first record bytes 0..159:")
for o in 0:16:159
    println("  +$o: ", join([string(b[r+1+o+j], base = 16, pad = 2) for j in 0:15], " "))
end
for o in 0:4:400
    x = rd(Float64, r + o); (isfinite(x) && 1.0 < abs(x) < 1e5) && println("  Float64 +$o = $x")
end
