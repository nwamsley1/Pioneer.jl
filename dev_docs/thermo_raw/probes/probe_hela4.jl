# HELA: find scan 937 (first MS2) record by its isolation centre and test the record-length rule
# [preamble P][u32 nreact][R x nreact][u32 X][u32 nrange][16 x nrange][u32 ncoef][8 x ncoef][12].
using Mmap, Arrow, Tables
const RAW = expanduser("~/Projects/pioneer-gui-test/rawtest/HELA_MOUSESIL_METHTEST_07312024_06.raw")
const ARW = expanduser("~/Projects/pioneer-gui-test/rawtest/arrow_out/HELA_MOUSESIL_METHTEST_07312024_06.arrow")
b = Mmap.mmap(RAW); t = Arrow.Table(ARW)
rd(T, o) = reinterpret(T, b[o+1:o+sizeof(T)])[1]
ev = 520189768; n = length(t.msOrder)
# walk with candidate (P, R) and check: every MS2 record's centre == arrow centre, ranges plausible
function walk(P, R; upto = n)
    r = ev + 4; bad = 0; firstbad = 0
    for k in 1:upto
        nreact = rd(UInt32, r + P)
        nreact > 10 && return (k, :nreact, nreact)
        c = nreact > 0 ? rd(Float64, r + P + 4) : NaN
        if t.msOrder[k] == 2 && Float32(c) != t.centerMz[k]
            bad += 1; firstbad == 0 && (firstbad = k)
        end
        q = r + P + 4 + R * Int(nreact) + 4
        nrange = rd(UInt32, q); nrange > 10 && return (k, :nrange, nrange)
        lo = rd(Float64, q + 4); Float32(lo) == t.lowMz[k] || (bad += 1; firstbad == 0 && (firstbad = k))
        q += 4 + 16 * Int(nrange)
        ncoef = rd(UInt32, q); ncoef > 20 && return (k, :ncoef, ncoef)
        r = q + 4 + 8 * Int(ncoef) + 12
    end
    (upto, :done, (bad = bad, firstbad = firstbad, end_at = r))
end
for (P, R) in ((132, 52), (132, 48), (136, 48), (136, 52), (132, 56))
    println("P=$P R=$R => ", walk(P, R; upto = 5000))
end
println("events section end (next addr) = 553794148")
