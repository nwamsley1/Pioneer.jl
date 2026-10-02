# Locate known values from the ground-truth .arrow inside the .raw: scan peaks (Float32/Float64), RTs, TICs.
using Arrow, Mmap, Tables
const RAW = expanduser("~/BrukerTims/rawformat/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw")
const ARW = expanduser("~/BrukerTims/tuning_ab/data/OlsenExploris500ng/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.arrow")
b = Mmap.mmap(RAW); t = Arrow.Table(ARW)
println("raw bytes $(length(b)); arrow scans $(length(t.msOrder)); MS1 $(count(==(1), t.msOrder)); columns $(Tables.columnnames(t))")
function findall_bytes(b, pat::Vector{UInt8}; limit = 5)
    hits = Int[]; n = length(pat); i = 1
    while length(hits) < limit
        r = findnext(pat, b, i); r === nothing && break
        push!(hits, first(r) - 1); i = first(r) + 1
    end
    hits
end
bytes(x) = collect(reinterpret(UInt8, [x]))
for i in (1, 2, findfirst(==(2), t.msOrder))
    mz = collect(skipmissing(t.mz_array[i])); it = collect(skipmissing(t.intensity_array[i]))
    println("\nscan $i order=$(t.msOrder[i]) rt=$(t.retentionTime[i]) npk=$(length(mz)) packet=$(t.packetType[i]) hdr=$(t.scanHeader[i])")
    println("  first peaks: ", collect(zip(mz[1:3], it[1:3])), " TIC=$(t.TIC[i]) low/high=$(t.lowMz[i])/$(t.highMz[i])")
    # m/z as Float32 sequence of 2 consecutive values, as Float64 single value, intensity Float32
    println("  mz[1] Float32 hits: ", findall_bytes(b, bytes(mz[1])))
    println("  mz[1:2] Float32 contiguous: ", findall_bytes(b, vcat(bytes(mz[1]), bytes(mz[2]))))
    println("  mz[1] as Float64 (exact from Float32): ", findall_bytes(b, bytes(Float64(mz[1]))))
    println("  int[1] Float32 hits: ", findall_bytes(b, bytes(it[1])))
    println("  (mz1,int1) Float32 interleaved: ", findall_bytes(b, vcat(bytes(mz[1]), bytes(it[1]))))
    println("  (int1,mz1)? : ", findall_bytes(b, vcat(bytes(it[1]), bytes(mz[1]))))
    println("  RT Float64: ", findall_bytes(b, bytes(Float64(t.retentionTime[i]))), "  RT Float32: ", findall_bytes(b, bytes(t.retentionTime[i])))
    println("  TIC Float64: ", findall_bytes(b, bytes(Float64(t.TIC[i]))), "  TIC Float32: ", findall_bytes(b, bytes(t.TIC[i])))
end
