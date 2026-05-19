# Per-file 1% FDR target IDs from MainSearch best-PSM-per-precursor output
# (main_search_psms/<file>_fold{0,1}.arrow), BEFORE the global q-value filter.

using Arrow, DataFrames, Printf, Tables

const RUN06 = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/06_fragmatch_diag"
const FILES = ["Sample-A_Rep1","Sample-A_Rep2","Sample-A_Rep3","Sample-B_Rep1","Sample-B_Rep2","Sample-B_Rep3"]

# Standard target-decoy q-value (sorted desc by score)
function qvalues_desc(scores::AbstractVector{<:Real}, targets::AbstractVector{<:Bool})
    perm = sortperm(scores; rev=true)
    n = length(scores)
    qv = Vector{Float32}(undef, n)
    t = 0; d = 0
    for i in perm
        targets[i] ? (t += 1) : (d += 1)
        qv[i] = d / max(t, 1)
    end
    # monotone non-increasing as scores decrease
    fdr = Inf32
    for i in reverse(perm)
        qv[i] = qv[i] > fdr ? fdr : (fdr = qv[i])
    end
    return qv
end

println("─── Per-file IDs at 1% FDR from MainSearch best-PSM-per-precursor ───")
println("(scoring via lgbm_prob; concat fold0+fold1 = best PSM per precursor)")
@printf "%-16s  %10s  %10s  %10s  %10s  %10s  %10s\n" "file" "n_total" "n_targets" "n_decoys" "T_q<0.01" "T_q<0.05" "T_q<0.10"
for fn in FILES
    f0 = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN06,"temp_data","main_search_psms","$(fn)_fold0.arrow"))))
    f1 = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN06,"temp_data","main_search_psms","$(fn)_fold1.arrow"))))
    df = vcat(f0, f1; cols=:union)
    # MainSearch's main_search_psms is already best-scan-per-precursor per fold
    n_total = nrow(df)
    n_t = count(df.target); n_d = n_total - n_t
    qv = qvalues_desc(Float32.(df.lgbm_prob), Bool.(df.target))
    Tq01 = sum((qv .< 0.01) .& df.target)
    Tq05 = sum((qv .< 0.05) .& df.target)
    Tq10 = sum((qv .< 0.10) .& df.target)
    @printf "%-16s  %10d  %10d  %10d  %10d  %10d  %10d\n" fn n_total n_t n_d Tq01 Tq05 Tq10
end
