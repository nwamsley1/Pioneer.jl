# HELA MS controller: scan-event record sizes and fields for ion-trap MS1 (packet 18) vs Orbitrap MS2 (packet 21).
using Mmap, Arrow, Tables
const RAW = expanduser("~/Projects/pioneer-gui-test/rawtest/HELA_MOUSESIL_METHTEST_07312024_06.raw")
const ARW = expanduser("~/Projects/pioneer-gui-test/rawtest/arrow_out/HELA_MOUSESIL_METHTEST_07312024_06.arrow")
b = Mmap.mmap(RAW); t = Arrow.Table(ARW)
rd(T, o) = reinterpret(T, b[o+1:o+sizeof(T)])[1]
rh = 505935206
n = rd(Int32, rh + 12); ev = Int(rd(UInt64, rh + 7448)); tr = Int(rd(UInt64, rh + 7456)); idx = Int(rd(UInt64, rh + 7408))
println("scans $n; events @$ev, next section @$tr; bytes/scan ", (tr - ev - 4) / n)
# which scans are MS2, and packet types
i2 = findfirst(==(2), t.msOrder); println("first MS2 arrow scan $i2: ", t.scanHeader[i2], " center ", t.centerMz[i2], " ptype ", t.packetType[i2])
println("ptypes by order: ", Dict(o => unique(t.packetType[t.msOrder .== o]) for o in unique(t.msOrder)))
# walk records with the Exploris count rule and see whether it stays consistent
r = ev + 4
for k in 0:5
    order = b[r+7]; nreact = rd(UInt32, r + 136)
    q = r + 140 + 48 * Int(nreact) + 8; nr = rd(UInt32, q); q2 = q + 4 + 16 * Int(nr); nc = rd(UInt32, q2)
    println("rec $k @$r: order $order nreact $nreact nrange $nr ncoef $nc | arrow order $(t.msOrder[k+1])")
    k == 0 && println("  bytes 0..63: ", join([string(b[r+1+j], base = 16, pad = 2) for j in 0:63], " "))
    global r = q2 + 4 + 8 * Int(nc) + 12
end
