# Two-round experiment-wide scoring: cross-run consistency features.
# Design + CV protocol: see TWO_ROUND_SCORING.md.
#
# Enabled via ENV["TWO_ROUND"] == "1". After the round-1 Pass-1 OOM trainer has
# written each file's `.pass1_sidecar.arrow` (OOF s1 = trace_prob_prepass),
# `write_two_round_feature_columns!` derives two cross-run features from s1 and
# writes them as columns into each per-file fold Arrow, so a second
# `train_and_predict_pass1_oom!` pass can train on [ base ; twin_score ; delta_irt ].
#
# Leakage-safety rests on the existing precursor_idx-keyed CV folds (Invariant A):
# every PSM of a precursor shares one fold, so a PSM's twin/best-instance lookups
# land in the same held-out fold and their s1 is OOF-consistent. The round-2 trainer
# reuses the unchanged `cv_fold` column (Invariant B).
#
# Assumes one row per precursor_idx per file (current MainSearch invariant), so
# s1[file, precursor] is a single value — no trace->precursor aggregation needed.

const TWO_ROUND_ENABLED = get(ENV, "TWO_ROUND", "") == "1"
const TWO_ROUND_FEATURES = [:twin_score, :delta_irt]

# Most-similar run per file via cosine similarity of the per-file s1 profiles
# (profile[file] = vector of best-per-precursor s1 over the precursor vocabulary,
# absent = 0). Returns (twin, dicts) where twin[f] is the index of f's most-similar
# other file (0 if none) and dicts[f] maps precursor_idx -> s1 for lookups.
function _compute_twin_runs(file_pid::Vector{Vector{UInt32}},
                            file_s1::Vector{Vector{Float32}})
    nf = length(file_pid)
    dicts = Vector{Dict{UInt32,Float32}}(undef, nf)
    for f in 1:nf
        d = Dict{UInt32,Float32}(); pid = file_pid[f]; s1 = file_s1[f]
        @inbounds for i in eachindex(pid); d[pid[i]] = s1[i]; end
        dicts[f] = d
    end
    # Gram matrix of profile dot-products (G[f,f] = squared L2 norm).
    G = zeros(Float64, nf, nf)
    for f in 1:nf, g in f:nf
        small, big = length(dicts[f]) <= length(dicts[g]) ? (dicts[f], dicts[g]) : (dicts[g], dicts[f])
        acc = 0.0
        for (p, v) in small
            acc += Float64(v) * Float64(get(big, p, 0.0f0))
        end
        G[f, g] = acc; G[g, f] = acc
    end
    twin = zeros(Int, nf)
    for f in 1:nf
        norm_f = sqrt(G[f, f]); best = -Inf; bg = 0
        norm_f <= 0 && continue
        for g in 1:nf
            (g == f || G[g, g] <= 0) && continue
            c = G[f, g] / (norm_f * sqrt(G[g, g]))
            c > best && (best = c; bg = g)
        end
        twin[f] = bg
    end
    return twin, dicts
end

"""
    write_two_round_feature_columns!(file_paths)

Compute `twin_score` and `delta_irt` from the round-1 OOF scores in each file's
`.pass1_sidecar.arrow` and write them as columns into the per-file fold Arrow files,
in place. Must run after round-1 `train_and_predict_pass1_oom!` and before round-2.
"""
function write_two_round_feature_columns!(file_paths::Vector{String})
    nf = length(file_paths)
    file_pid = Vector{Vector{UInt32}}(undef, nf)
    file_s1  = Vector{Vector{Float32}}(undef, nf)
    file_irt = Vector{Vector{Float32}}(undef, nf)
    best_s1  = Dict{UInt32,Float32}()   # precursor -> highest s1 seen (for ref_irt)
    ref_irt  = Dict{UInt32,Float32}()   # precursor -> irt_obs at its highest-s1 instance

    # Pass A: read slim (precursor_idx, irt_obs, s1) for every file; track ref_irt.
    for f in 1:nf
        p = file_paths[f]
        tbl = Arrow.Table(p)
        side = Arrow.Table(p * PASS1_SIDECAR_SUFFIX)
        pid = collect(UInt32.(tbl.precursor_idx))
        irt = collect(Float32.(tbl.irt_obs))
        s1  = collect(Float32.(side.trace_prob_prepass))
        length(s1) == length(pid) ||
            error("write_two_round_feature_columns!: sidecar row mismatch at $p")
        file_pid[f] = pid; file_s1[f] = s1; file_irt[f] = irt
        @inbounds for i in eachindex(pid)
            if s1[i] > get(best_s1, pid[i], -1.0f0)
                best_s1[pid[i]] = s1[i]; ref_irt[pid[i]] = irt[i]
            end
        end
    end

    twin, dicts = _compute_twin_runs(file_pid, file_s1)

    # Pass B: per file, compute the two features and rewrite the Arrow with them.
    for f in 1:nf
        pid = file_pid[f]; irt = file_irt[f]; n = length(pid)
        tw = Vector{Float32}(undef, n); di = Vector{Float32}(undef, n)
        tdict = twin[f] > 0 ? dicts[twin[f]] : nothing
        @inbounds for i in 1:n
            tw[i] = tdict === nothing ? 0.0f0 : get(tdict, pid[i], 0.0f0)
            di[i] = abs(irt[i] - ref_irt[pid[i]])
        end
        main = DataFrame(Tables.columntable(Arrow.Table(file_paths[f])))
        nrow(main) == n ||
            error("write_two_round_feature_columns!: arrow row count changed at $(file_paths[f])")
        main[!, :twin_score] = tw
        main[!, :delta_irt]  = di
        writeArrow(file_paths[f], main)
    end
    n_twinned = count(>(0), twin)
    @user_info "two-round: wrote twin_score + delta_irt to $nf files ($n_twinned with a twin, cosine s1-profile)"
    return nothing
end
