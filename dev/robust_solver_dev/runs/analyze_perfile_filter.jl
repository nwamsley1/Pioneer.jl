# For each file, compute per-file q-value from MainSearch's lgbm_prob and
# count additional targets/decoys removed by adding a per-file q≤0.10 filter
# ON TOP OF the existing global q≤0.10 prescore filter.
#
# Also reports counts in/out of intersection with global passing.

using Arrow, DataFrames, Printf, Tables

const RUN07 = "/Users/nathanwamsley/Data/RegressionTestsLite/MTAC_3P_Standard/sweep_l1robust/07_n_scans_feature"
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
    fdr = Inf32
    for i in reverse(perm)
        qv[i] = qv[i] > fdr ? fdr : (fdr = qv[i])
    end
    return qv
end

# Global passing precs (q≤0.10 prescore aggregation across all files):
# union of precursor_idx in second_pass_psms (which is the post-global-filter
# input fed to PrecursorScoringSearch).
passing = Set{UInt32}()
for fn in FILES
    t = Arrow.Table(joinpath(RUN07,"temp_data","second_pass_psms","$(fn).arrow"))
    for p in t.precursor_idx; push!(passing, UInt32(p)); end
end

println("─── Per-file q≤0.10 filter on top of global q≤0.10 ───")
@printf "%-16s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n" "file" "n_total" "global_T" "global_D" "perfile_T" "perfile_D" "both_T" "both_D" "addΔ_T" "addΔ_D"

for fn in FILES
    f0 = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN07,"temp_data","main_search_psms","$(fn)_fold0.arrow"))))
    f1 = DataFrame(Tables.columntable(Arrow.Table(joinpath(RUN07,"temp_data","main_search_psms","$(fn)_fold1.arrow"))))
    df = vcat(f0, f1; cols=:union)
    n = nrow(df)
    qv = qvalues_desc(Float32.(df.lgbm_prob), Bool.(df.target))
    in_global  = [pid in passing for pid in df.precursor_idx]
    in_perfile = qv .<= 0.10f0
    in_both    = in_global .& in_perfile
    targets    = Bool.(df.target)
    g_T = sum(in_global .& targets);  g_D = sum(in_global) - g_T
    p_T = sum(in_perfile .& targets); p_D = sum(in_perfile) - p_T
    b_T = sum(in_both .& targets);    b_D = sum(in_both) - b_T
    # additional dropped by adding the per-file filter on top of the global filter
    addΔ_T = g_T - b_T
    addΔ_D = g_D - b_D
    @printf "%-16s  %10d  %10d  %10d  %10d  %10d  %10d  %10d  %10d  %10d\n" fn n g_T g_D p_T p_D b_T b_D addΔ_T addΔ_D
end
