# For each file: of the precursors additionally dropped by the per-file
# q≤0.10 filter (passed global, failed per-file), how many were ID'd by
# HUPO 2026 (q<0.01 target) in that same file?

using Arrow, DataFrames, Printf, Tables

const RUN07 = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/07_n_scans_feature"
const HUPO_PATH = "/Users/nathanwamsley/Data/HUPO2026/MtacThreeProteomeStandard_hupo2026_v0o6o6/precursors_long.arrow"
const LIB_PATH  = "/Users/nathanwamsley/Data/SPEC_LIBS/altimeter_3P_len7o40_ch2o3_mc1_OlsenAstral_mzsorted.poin/precursors_table.arrow"
const FILES = ["Sample-A_Rep1","Sample-A_Rep2","Sample-A_Rep3","Sample-B_Rep1","Sample-B_Rep2","Sample-B_Rep3"]

_norm(x) = ismissing(x) ? "" : String(x)
key_t(seq,ch,sm,im,dc) = (String(seq), UInt8(ch), _norm(sm), _norm(im), Bool(dc))

function qvalues_desc(scores::AbstractVector{<:Real}, targets::AbstractVector{<:Bool})
    perm = sortperm(scores; rev=true)
    n = length(scores); qv = Vector{Float32}(undef, n)
    t = 0; d = 0
    for i in perm; targets[i] ? (t += 1) : (d += 1); qv[i] = d / max(t, 1); end
    fdr = Inf32
    for i in reverse(perm); qv[i] = qv[i] > fdr ? fdr : (fdr = qv[i]); end
    return qv
end

# Library map
lib = DataFrame(Arrow.Table(LIB_PATH))
lib_keyt = Dict{Tuple{String,UInt8,String,String,Bool}, UInt32}()
@inbounds for i in 1:nrow(lib)
    lib_keyt[key_t(lib.sequence[i], lib.prec_charge[i], lib.structural_mods[i], lib.isotopic_mods[i], lib.is_decoy[i])] = UInt32(i)
end

# HUPO q<0.01 target per file → set of OUR pids
hupo = DataFrame(Arrow.Table(HUPO_PATH))
hupo = hupo[(hupo.MBR_boosted_qval .< 0.01) .& hupo.target, :]
hupo_by_file = Dict{String, Set{UInt32}}()
for fn in FILES
    sub = hupo[hupo.file_name .== fn, :]
    s = Set{UInt32}()
    for r in eachrow(sub)
        pid = get(lib_keyt, key_t(r.sequence, r.charge, r.structural_mods, r.isotopic_mods, r.is_decoy), UInt32(0))
        pid != 0 && push!(s, pid)
    end
    hupo_by_file[fn] = s
end

# Global passing per file (precursor_idx in second_pass_psms = post-global-filter)
println("─── Per-file additional drops vs HUPO IDs ───")
@printf "%-16s  %10s  %10s  %10s  %10s  %10s\n" "file" "addΔ_T" "addΔ_D" "addΔ_T∩HUPO" "%T_in_HUPO" "addΔ_T_notHUPO"

for fn in FILES
    f0 = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN07,"temp_data","main_search_psms","$(fn)_fold0.arrow"))))
    f1 = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN07,"temp_data","main_search_psms","$(fn)_fold1.arrow"))))
    main = vcat(f0, f1; cols=:union)

    # Global passing set for this file = unique precursor_idx in second_pass_psms
    sp = Arrow.Table(joinpath(RUN07,"temp_data","second_pass_psms","$(fn).arrow"))
    global_pass = Set{UInt32}(UInt32(p) for p in sp.precursor_idx)

    # Per-file q-value
    qv = qvalues_desc(Float32.(main.lgbm_prob), Bool.(main.target))
    perfile_pass_mask = qv .<= 0.10f0

    in_global  = [pid in global_pass for pid in main.precursor_idx]
    additional_drop_mask = in_global .& .!perfile_pass_mask  # passed global, failed per-file
    targets = Bool.(main.target)

    add_T_pids = Set{UInt32}(UInt32.(main.precursor_idx[additional_drop_mask .& targets]))
    add_D_pids = Set{UInt32}(UInt32.(main.precursor_idx[additional_drop_mask .& .!targets]))
    H = hupo_by_file[fn]
    add_T_in_hupo = length(intersect(add_T_pids, H))
    add_T_not_hupo = length(add_T_pids) - add_T_in_hupo
    pct = length(add_T_pids) == 0 ? 0.0 : 100*add_T_in_hupo/length(add_T_pids)
    @printf "%-16s  %10d  %10d  %10d  %9.2f%%  %10d\n" fn length(add_T_pids) length(add_D_pids) add_T_in_hupo pct add_T_not_hupo
end
