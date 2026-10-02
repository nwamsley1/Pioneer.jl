# HELA: dump the first MS2 scan-event record (scan 937), reached by walking 936 records with P=132.
using Mmap, Arrow, Tables
const RAW = expanduser("~/Projects/pioneer-gui-test/rawtest/HELA_MOUSESIL_METHTEST_07312024_06.raw")
const ARW = expanduser("~/Projects/pioneer-gui-test/rawtest/arrow_out/HELA_MOUSESIL_METHTEST_07312024_06.arrow")
b = Mmap.mmap(RAW); t = Arrow.Table(ARW)
rd(T, o) = reinterpret(T, b[o+1:o+sizeof(T)])[1]
ev = 520189768; P = 132
r = ev + 4
for k in 1:936
    nreact = Int(rd(UInt32, r + P)); q = r + P + 4 + 52nreact + 4
    nrange = Int(rd(UInt32, q)); q += 4 + 16nrange; ncoef = Int(rd(UInt32, q))
    global r = q + 4 + 8ncoef + 12
end
println("scan 937 record @$r; arrow: ", t.scanHeader[937], " center ", t.centerMz[937], " width ", t.isolationWidthMz[937], " CE ", t.collisionEnergyField[937], " low/high ", t.lowMz[937], "/", t.highMz[937])
for o in 0:16:335
    println("  +$o: ", join([string(b[r+1+o+j], base = 16, pad = 2) for j in 0:15], " "))
end
for o in 0:4:330
    x = rd(Float64, r + o); (isfinite(x) && 0.5 < abs(x) < 1e5) && println("  Float64 +$o = $x")
end
