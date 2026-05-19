# For missing precursors that reach second_pass_psms, distribution of their
# best q_value (the per-row q-value that gates passing_psms via threshold 0.01).

using Arrow, DataFrames, Printf, Statistics, Tables

const RUN06     = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/06_fragmatch_diag"
const BASE_PATH = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/04_develop_baseline/precursors_long.arrow"
const HUPO_PATH = "/Users/nathanwamsley/Data/HUPO2026/MtacThreeProteomeStandard_hupo2026_v0o6o6/precursors_long.arrow"
const LIB_PATH  = "/Users/nathanwamsley/Data/SPEC_LIBS/altimeter_3P_len7o40_ch2o3_mc1_OlsenAstral_mzsorted.poin/precursors_table.arrow"
const FILES = ["Sample-A_Rep1","Sample-A_Rep2","Sample-A_Rep3","Sample-B_Rep1","Sample-B_Rep2","Sample-B_Rep3"]

_norm(x) = ismissing(x) ? "" : String(x)
key_t(seq,ch,sm,im,dc) = (String(seq), UInt8(ch), _norm(sm), _norm(im), Bool(dc))
key3(seq,ch,sm,dc)     = (String(seq), UInt8(ch), _norm(sm), Bool(dc))

println("[load] library + maps")
lib = DataFrame(Arrow.Table(LIB_PATH))
lib_key2pid  = Dict{Tuple{String,UInt8,String,String,Bool}, UInt32}()
lib_key3_pid = Dict{Tuple{String,UInt8,String,Bool}, UInt32}()
@inbounds for i in 1:nrow(lib)
    lib_key2pid[key_t(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i], lib.isotopic_mods[i], lib.is_decoy[i])] = UInt32(i)
    lib_key3_pid[key3(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i], lib.is_decoy[i])] = UInt32(i)
end
hupo = DataFrame(Arrow.Table(HUPO_PATH)); hupo = hupo[(hupo.MBR_boosted_qval .< 0.01) .& hupo.target, :]
base = DataFrame(Arrow.Table(BASE_PATH)); base = base[(base.qval .< 0.01) .& base.target, :]

quant(x, p) = (xs = sort(x); xs[clamp(round(Int, p*length(xs)), 1, length(xs))])

println("\n─── Missing precursors' best q_value (spline qval from prec_prob) ───")
@printf "%-16s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n" "file" "miss_n" "<.01" "<.02" "<.05" "<.10" "p50" "p90"

for fn in FILES
    H = Set{UInt32}()
    for r in eachrow(hupo[hupo.file_name .== fn, :])
        pid = get(lib_key2pid, key_t(r.sequence, r.charge, r.structural_mods, r.isotopic_mods, r.is_decoy), UInt32(0))
        pid != 0 && push!(H, pid)
    end
    B = Set{UInt32}()
    for r in eachrow(base[base.file_name .== fn, :])
        pid = get(lib_key3_pid, key3(r.sequence, r.charge, r.structural_mods, r.is_decoy), UInt32(0))
        pid != 0 && push!(B, pid)
    end
    miss = setdiff(H, B)

    sp = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN06,"temp_data","second_pass_psms","$fn.arrow"))))
    qmin = Dict{UInt32, Float32}()
    for r in eachrow(sp)
        pid = UInt32(r.precursor_idx); q = Float32(r.q_value)
        if !haskey(qmin, pid) || q < qmin[pid]
            qmin[pid] = q
        end
    end
    qs = Float32[qmin[p] for p in miss if haskey(qmin, p)]
    isempty(qs) && continue
    @printf "%-16s  %8d  %8d  %8d  %8d  %8d  %8.4f  %8.4f\n" fn length(qs) count(<(0.01f0), qs) count(<(0.02f0), qs) count(<(0.05f0), qs) count(<(0.10f0), qs) quant(qs, 0.5) quant(qs, 0.9)
end
