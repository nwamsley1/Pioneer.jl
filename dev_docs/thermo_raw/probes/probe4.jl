# Where does the file record the scan-index start (1320778831) and data-stream start (110544)?
using Mmap
const RAW = expanduser("~/BrukerTims/rawformat/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw")
b = Mmap.mmap(RAW)
rd(T, o) = reinterpret(T, b[o+1:o+sizeof(T)])[1]
function hits(x; limit = 10)
    pat = collect(reinterpret(UInt8, [x])); out = Int[]; i = 1
    while length(out) < limit
        r = findnext(pat, b, i); r === nothing && break
        push!(out, first(r) - 1); i = first(r) + 1
    end
    out
end
const IDX0 = 1320778831; const DATA0 = 110544; const N = 62196
idx_end = IDX0 + 88N
for (name, v) in (("index start", IDX0), ("data start", DATA0), ("index end", idx_end), ("n scans", N),
                  ("last scan number", N), ("data stream end?", 0))
    v == 0 && continue
    println(rpad(name, 18), " UInt64 hits: ", hits(UInt64(v)), "  UInt32 hits: ", hits(UInt32(v))[1:min(end, 6)])
end
println("bytes after index end ($idx_end): ", [rd(UInt32, idx_end + 4j) for j in 0:15])
println("file size $(length(b)); bytes left after index: $(length(b) - idx_end)")
