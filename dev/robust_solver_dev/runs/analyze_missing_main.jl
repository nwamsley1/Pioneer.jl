# MainSearch (per-file LightGBM) score distribution for missing vs survived
# precursors, plus q-value at MainSearch fold level.

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

# Per-file MainSearch fold q-values: rebuild by ranking lgbm_prob with target/decoy.
# (MainSearch summarize_results uses per-fold per-file LightGBM probs; baseline
#  prescore filter applies a per-precursor aggregation.) Here just report
#  the per-file lgbm_prob distribution at the precursor level.

println("\n─── MainSearch per-file lgbm_prob: missing vs survived ───")
@printf "%-16s  %-8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n" "file" "set" "n" "p10" "p25" "p50" "p75" "p90" ">.99" ">.95"

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

    main = vcat(
        DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN06,"temp_data","main_search_psms","$(fn)_fold0.arrow")))),
        DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN06,"temp_data","main_search_psms","$(fn)_fold1.arrow"))));
        cols=:union)
    pid2lgbm = Dict{UInt32, Float32}()
    for r in eachrow(main)
        pid = UInt32(r.precursor_idx); p = Float32(r.lgbm_prob)
        if !haskey(pid2lgbm, pid) || p > pid2lgbm[pid]
            pid2lgbm[pid] = p
        end
    end

    function row(label, S)
        ps = Float32[pid2lgbm[p] for p in S if haskey(pid2lgbm, p)]
        isempty(ps) && return
        @printf "%-16s  %-8s  %8d  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n" fn label length(ps) quant(ps,0.1) quant(ps,0.25) quant(ps,0.5) quant(ps,0.75) quant(ps,0.9) (count(>(0.99f0), ps)/length(ps)) (count(>(0.95f0), ps)/length(ps))
    end
    row("missing", miss)
    row("survived", B)
end
