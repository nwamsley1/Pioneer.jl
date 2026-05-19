# Funnel: where do HUPO N=6 q<0.01 targets that baseline (04) MISSES at q<0.01
# die in the pipeline? Use 06_fragmatch_diag's temp_data for upstream stages
# (MainSearch is identical to baseline since both run on develop@same point).

using Arrow, DataFrames, Printf, Statistics

const RUN06     = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/06_fragmatch_diag"
const BASE_PATH = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/04_develop_baseline/precursors_long.arrow"
const HUPO_PATH = "/Users/nathanwamsley/Data/HUPO2026/MtacThreeProteomeStandard_hupo2026_v0o6o6/precursors_long.arrow"
const LIB_PATH  = "/Users/nathanwamsley/Data/SPEC_LIBS/altimeter_3P_len7o40_ch2o3_mc1_OlsenAstral_mzsorted.poin/precursors_table.arrow"
const FILES = ["Sample-A_Rep1","Sample-A_Rep2","Sample-A_Rep3","Sample-B_Rep1","Sample-B_Rep2","Sample-B_Rep3"]

_norm(x) = ismissing(x) ? "" : String(x)
key_t(seq,ch,sm,im,dc) = (String(seq), UInt8(ch), _norm(sm), _norm(im), Bool(dc))
key3(seq,ch,sm,dc) = (String(seq), UInt8(ch), _norm(sm), Bool(dc))

println("[load] library + maps")
lib = DataFrame(Arrow.Table(LIB_PATH))
lib_key2pid = Dict{Tuple{String,UInt8,String,String,Bool}, UInt32}()
sizehint!(lib_key2pid, nrow(lib))
@inbounds for i in 1:nrow(lib)
    lib_key2pid[key_t(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i],
                      lib.isotopic_mods[i], lib.is_decoy[i])] = UInt32(i)
end
# baseline output lacks isotopic_mods → 4-key map
lib_key3_to_pid = Dict{Tuple{String,UInt8,String,Bool}, UInt32}()
@inbounds for i in 1:nrow(lib)
    lib_key3_to_pid[key3(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i],
                         lib.is_decoy[i])] = UInt32(i)
end

# HUPO -> per-file set of OUR lib precursor_idx
println("[load] HUPO ref")
hupo = DataFrame(Arrow.Table(HUPO_PATH))
hupo = hupo[(hupo.MBR_boosted_qval .< 0.01) .& hupo.target, :]
hupo_by_file = Dict{String, Set{UInt32}}()
for fn in FILES
    sub = hupo[hupo.file_name .== fn, :]
    s = Set{UInt32}()
    for r in eachrow(sub)
        pid = get(lib_key2pid, key_t(r.sequence, r.charge, r.structural_mods, r.isotopic_mods, r.is_decoy), UInt32(0))
        pid != 0 && push!(s, pid)
    end
    hupo_by_file[fn] = s
end

# Baseline -> per-file set of OUR lib precursor_idx (target q<0.01)
println("[load] baseline 04")
base = DataFrame(Arrow.Table(BASE_PATH))
base = base[(base.qval .< 0.01) .& base.target, :]
base_by_file = Dict{String, Set{UInt32}}()
for fn in FILES
    sub = base[base.file_name .== fn, :]
    s = Set{UInt32}()
    for r in eachrow(sub)
        pid = get(lib_key3_to_pid, key3(r.sequence, r.charge, r.structural_mods, r.is_decoy), UInt32(0))
        pid != 0 && push!(s, pid)
    end
    base_by_file[fn] = s
end

# Global passing (q≤0.10)
passing = Set{UInt32}(DataFrame(Arrow.Table(joinpath(RUN06,"temp_data","passing_precursors.arrow"))).precursor_idx)

println("\n─── Funnel for HUPO_q<0.01 ∩ baseline-MISSED, per file ───")
@printf "%-16s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n" "file" "miss" "in_main" "med_lgbm" "in_pass10" "in_2nd" "2nd_qv01" "passing" "in_06any" "06qv01"

for fn in FILES
    H = hupo_by_file[fn]
    B = base_by_file[fn]
    miss = setdiff(H, B)
    # MainSearch PSMs (best per precursor across folds)
    main0 = DataFrame(Arrow.Table(joinpath(RUN06,"temp_data","main_search_psms","$(fn)_fold0.arrow")))
    main1 = DataFrame(Arrow.Table(joinpath(RUN06,"temp_data","main_search_psms","$(fn)_fold1.arrow")))
    main = vcat(main0, main1; cols=:union)
    pid_to_lgbm = Dict{UInt32, Float32}()
    for r in eachrow(main)
        pid = UInt32(r.precursor_idx); p = Float32(r.lgbm_prob)
        if !haskey(pid_to_lgbm, pid) || p > pid_to_lgbm[pid]
            pid_to_lgbm[pid] = p
        end
    end
    in_main = count(p -> haskey(pid_to_lgbm, p), miss)
    miss_lgbm = Float32[pid_to_lgbm[p] for p in miss if haskey(pid_to_lgbm, p)]
    med_lgbm = isempty(miss_lgbm) ? NaN32 : median(miss_lgbm)

    # passing precs at q≤0.10
    in_pass = count(p -> p in passing, miss)

    # second_pass_psms (06's path: single file per ms_file, has q_value)
    sp_pids = Set{UInt32}()
    sp_qmin = Dict{UInt32, Float32}()
    sp_path = joinpath(RUN06,"temp_data","second_pass_psms","$(fn).arrow")
    if isfile(sp_path)
        t = DataFrame(Arrow.Table(sp_path))
        for r in eachrow(t)
            pid = UInt32(r.precursor_idx); q = Float32(r.q_value)
            push!(sp_pids, pid)
            if !haskey(sp_qmin, pid) || q < sp_qmin[pid]
                sp_qmin[pid] = q
            end
        end
    end
    in_2nd = count(p -> p in sp_pids, miss)
    in_2nd_q01 = count(p -> haskey(sp_qmin, p) && sp_qmin[p] < 0.01, miss)

    # passing_psms: post prec_prob → qval<0.01 filter, fed to IntegrateChromatograms
    pp_pids = Set{UInt32}()
    pp_path = joinpath(RUN06,"temp_data","passing_psms","$(fn).arrow")
    if isfile(pp_path)
        t = DataFrame(Arrow.Table(pp_path))
        for p in t.precursor_idx; push!(pp_pids, UInt32(p)); end
    end
    in_passing = count(p -> p in pp_pids, miss)

    # 06 run's precursors_long: did it score q<0.01 in 06 (with secondpass)?
    pl06 = DataFrame(Arrow.Table(joinpath(RUN06,"precursors_long.arrow")))
    qcol = hasproperty(pl06, :MBR_boosted_qval) ? :MBR_boosted_qval : :qval
    pl06 = pl06[pl06.file_name .== fn, :]
    pid06 = Dict{UInt32, Float32}()
    for r in eachrow(pl06)
        if !ismissing(r[qcol])
            pid = get(lib_key3_to_pid, key3(r.sequence, r.charge, r.structural_mods, r.is_decoy), UInt32(0))
            pid != 0 && (pid06[pid] = Float32(r[qcol]))
        end
    end
    in_06_any = count(p -> haskey(pid06, p), miss)
    in_06_qv01 = count(p -> haskey(pid06, p) && pid06[p] < 0.01, miss)
    miss_qv = Float32[pid06[p] for p in miss if haskey(pid06, p)]
    med_qv = isempty(miss_qv) ? NaN32 : median(miss_qv)

    @printf "%-16s  %8d  %8d  %8.3f  %8d  %8d  %8d  %8d  %8d  %8d\n" fn length(miss) in_main med_lgbm in_pass in_2nd in_2nd_q01 in_passing in_06_any in_06_qv01
end
