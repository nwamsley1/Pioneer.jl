# Dump the candidate scan index (88-byte entries around the TIC hits) and the packet bytes before scan 1's peaks.
using Mmap
const RAW = expanduser("~/BrukerTims/rawformat/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw")
b = Mmap.mmap(RAW)
rd(T, o) = reinterpret(T, b[o+1:o+sizeof(T)])[1]
tic1 = 1320778859
# try every alignment of an 88-byte entry containing the TIC: print entries as mixed fields
for start in (tic1 - 88 + 1):(tic1)
    # look for an alignment where an Int32 field equals 0 (scan index 0) or 1 and a UInt64 looks like a file offset
end
base = tic1 - 40
println("== region around index (base $base), 3 entries")
for e in 0:3
    o = base + 88e
    print("entry $e @$o: ")
    for j in 0:4:84
        print(rd(UInt32, o + j), " ")
    end
    println()
    print("     as Float64 @+0,+8..: ")
    for j in 0:8:80
        x = rd(Float64, o + j); print(isfinite(x) && abs(x) < 1e12 ? round(x, sigdigits = 7) : "~", " ")
    end
    println()
end
println("\n== bytes before scan 1 peaks (111456)")
for o in 111456-160:16:111456-1
    print(o, ": ")
    for j in 0:4:12
        print(rd(UInt32, o + j), "(", round(rd(Float32, o + j), sigdigits = 6), ") ")
    end
    println()
end
