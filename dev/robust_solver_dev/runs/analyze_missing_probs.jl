# For HUPO N=6 q<0.01 targets baseline missed, profile their per-precursor
# scores in second_pass_psms vs the survivors of the prec_prob qval filter.

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
println("[load] HUPO + baseline")
hupo = DataFrame(Arrow.Table(HUPO_PATH))
hupo = hupo[(hupo.MBR_boosted_qval .< 0.01) .& hupo.target, :]
base = DataFrame(Arrow.Table(BASE_PATH))
base = base[(base.qval .< 0.01) .& base.target, :]

quant(x, p) = (xs = sort(x); xs[clamp(round(Int, p*length(xs)), 1, length(xs))])

println("\n─── Per-file score distribution: missing vs survived ───")
@printf "%-16s  %-8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n" "file" "set" "n" "trace_p50" "trace_p10" "prec_p50" "prec_p10" "gqv_p50" "gqv_p90"

for fn in FILES
    # missing set in this file
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

    # second_pass_psms: collapse to per-precursor (max trace_prob row)
    sp = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN06,"temp_data","second_pass_psms","$fn.arrow"))))
    # global_qval not in second_pass_psms but in passing_psms; use trace_prob, prec_prob, and per-row q_value
    sort!(sp, [:precursor_idx, :trace_prob], rev=[false, true])
    sp_best = combine(groupby(sp, :precursor_idx),
                      :trace_prob => first => :trace_prob,
                      :prec_prob => first => :prec_prob,
                      :q_value => first => :q_value)

    # passing_psms: load global_qval per precursor
    pp = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN06,"temp_data","passing_psms","$fn.arrow"))))
    gqv_dict = Dict{UInt32, Float32}()
    if hasproperty(pp, :global_qval)
        for r in eachrow(pp)
            pid = UInt32(r.precursor_idx)
            if !haskey(gqv_dict, pid) || Float32(r.global_qval) < gqv_dict[pid]
                gqv_dict[pid] = Float32(r.global_qval)
            end
        end
    end

    miss_in_sp = filter(r -> UInt32(r.precursor_idx) in miss, sp_best)
    surv_in_sp = filter(r -> UInt32(r.precursor_idx) in B,   sp_best)

    function row(label, df)
        n = nrow(df)
        n == 0 && return
        tps = df.trace_prob; pps = df.prec_prob
        # global_qval lookup — survived likely all <0.01; missing varies
        gqs = Float32[]
        for pid in df.precursor_idx
            haskey(gqv_dict, UInt32(pid)) && push!(gqs, gqv_dict[UInt32(pid)])
        end
        gqs_p50 = isempty(gqs) ? NaN32 : quant(gqs, 0.5)
        gqs_p90 = isempty(gqs) ? NaN32 : quant(gqs, 0.9)
        @printf "%-16s  %-8s  %8d  %8.3f  %8.3f  %8.3f  %8.3f  %8.4f  %8.4f\n" fn label n quant(tps, 0.5) quant(tps, 0.1) quant(pps, 0.5) quant(pps, 0.1) gqs_p50 gqs_p90
    end
    row("missing", miss_in_sp)
    row("survived", surv_in_sp)
end
