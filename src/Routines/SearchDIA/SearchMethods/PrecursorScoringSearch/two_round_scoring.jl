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

# Read ENV at RUNTIME (not a top-level const): a `const = get(ENV,...)` is
# evaluated at precompile time and baked into the cached module, so it would
# ignore the env var on a cached load. Mirrors the GLOBAL_AGG runtime check.
# K for the k-nearest-neighbor cross-run feature (aggregation of top-K neighbors'
# same-precursor round-1 scores). To sweep K, change KNN_K. To switch aggregation,
# change KNN_AGG to :mean, :median, or :quantile (with KNN_Q). Column name auto-
# derives (e.g. :knn10_median, :knn10_q75).
const KNN_K = 5
const KNN_AGG = :mean      # :mean | :median | :quantile
const KNN_Q = 0.75f0       # only used when KNN_AGG === :quantile
# Zero out per-row round-1 s1 for rows whose LOCAL (binned target/decoy) PEP
# exceeds this threshold BEFORE computing cosine similarity + top-K lookup.
# 1.0 disables the filter. Only rows kept above this threshold contribute
# either to a run's cosine similarity profile or to the neighbor lookup mean.
const PEP_FILTER_THRESHOLD = 0.5
const PEP_FILTER_N_BINS = 100
const KNN_COL = KNN_AGG === :quantile ?
    Symbol("knn$(KNN_K)_q$(Int(round(KNN_Q * 100)))") :
    Symbol("knn$(KNN_K)_$(KNN_AGG)")
two_round_enabled() = get(ENV, "TWO_ROUND", "") == "1"
const TWO_ROUND_FEATURES = [KNN_COL, :delta_irt]

# Type-7 quantile of the first `n` elements of `buf` via in-place insertion sort.
# Allocation-free (caller owns `buf`). Insertion sort is fine for K ≤ 16. Passing
# q=0.5f0 gives the standard median.
@inline function _quantile_first_n!(buf::AbstractVector{Float32}, n::Int, q::Float32)
    @inbounds for i in 2:n
        v = buf[i]; j = i
        while j > 1 && buf[j-1] > v
            buf[j] = buf[j-1]; j -= 1
        end
        buf[j] = v
    end
    n == 1 && return @inbounds buf[1]
    h  = 1f0 + Float32(n - 1) * q          # 1-indexed fractional position
    lo = floor(Int, h)
    lo >= n && return @inbounds buf[n]
    lo <  1 && return @inbounds buf[1]
    frac = h - Float32(lo)
    return @inbounds buf[lo] + frac * (buf[lo+1] - buf[lo])
end

# Zero out per-row s1 across all files where the global local PEP > threshold.
# Local PEP is estimated by sorting all rows descending by s1, splitting into
# equal-count bins, and taking bin_pep = n_decoys / n_total_in_bin per bin. A
# row's PEP is its bin's PEP (rank-to-bin lookup). Simple, monotonically-noisy
# in principle but fine at bin counts ≥ 50 with the row counts we see here
# (millions of rows). Returns (n_filtered, n_kept).
function _filter_by_local_pep!(file_s1::Vector{Vector{Float32}},
                              file_tgt::Vector{Vector{Bool}},
                              pep_threshold::Float64,
                              n_bins::Int)
    # Concat rows into a single flat view for binning.
    nf = length(file_s1)
    offsets = Vector{Int}(undef, nf + 1); offsets[1] = 0
    for f in 1:nf
        offsets[f+1] = offsets[f] + length(file_s1[f])
    end
    ntot = offsets[end]
    ntot == 0 && return (0, 0)
    scores  = Vector{Float32}(undef, ntot)
    targets = Vector{Bool}(undef, ntot)
    for f in 1:nf, i in eachindex(file_s1[f])
        scores[offsets[f] + i]  = file_s1[f][i]
        targets[offsets[f] + i] = file_tgt[f][i]
    end
    perm = sortperm(scores; rev=true)  # rank 1 = highest score
    bin_size = max(1, cld(ntot, n_bins))
    n_actual_bins = cld(ntot, bin_size)
    bin_pep = Vector{Float64}(undef, n_actual_bins)
    for b in 1:n_actual_bins
        lo = 1 + (b-1) * bin_size
        hi = min(b * bin_size, ntot)
        n_target = 0
        @inbounds for i in lo:hi
            targets[perm[i]] && (n_target += 1)
        end
        bin_n = hi - lo + 1
        bin_pep[b] = 1.0 - n_target / bin_n
    end
    # rank[i] = 1..ntot; from perm we build the inverse
    rank = Vector{Int}(undef, ntot)
    @inbounds for r in 1:ntot
        rank[perm[r]] = r
    end
    n_filtered = 0
    @inbounds for i in 1:ntot
        b = min(cld(rank[i], bin_size), n_actual_bins)
        if bin_pep[b] > pep_threshold
            scores[i] = 0.0f0
            n_filtered += 1
        end
    end
    # Scatter filtered scores back into per-file vectors.
    for f in 1:nf, i in eachindex(file_s1[f])
        file_s1[f][i] = scores[offsets[f] + i]
    end
    return (n_filtered, ntot - n_filtered)
end

# Log a pass's LightGBM feature gains at @user_info (visible regardless of debug
# level), top-N sorted. Used to compare round-1 vs round-2 gains.
function log_pass_importance(pass, label::String; topn::Int = 30)
    pass.last_classifier === nothing && return nothing
    model = LightGBMModel(pass.last_classifier, pass.available_features, nothing)
    imp = importance(model)
    imp === nothing && return nothing
    sorted_imp = sort(imp, by = x -> -x[2])
    total = sum(x -> x[2], sorted_imp)
    lines = ["$label LGBM feature gains (top $(min(topn,length(sorted_imp))) of $(length(sorted_imp))):"]
    for (fname, gain) in first(sorted_imp, topn)
        pct = total > 0 ? round(100*gain/total, digits=1) : 0.0
        push!(lines, "    $(rpad(string(fname), 42)) $(rpad(round(Int, gain), 10)) ($(pct)%)")
    end
    @user_info join(lines, "\n")
    return nothing
end

# Top-K most-similar runs per file via cosine similarity of the per-file s1 profiles
# (profile[file] = vector of best-per-precursor s1 over the precursor vocabulary,
# absent = 0). Returns (topk, dicts) where topk[f] is a length-KNN_K Vector{Int} of
# neighbor file indices sorted descending by similarity (0 for slots with no valid
# neighbor) and dicts[f] maps precursor_idx -> s1 for lookups.
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
    topk = [zeros(Int, KNN_K) for _ in 1:nf]
    for f in 1:nf
        norm_f = sqrt(G[f, f])
        norm_f <= 0 && continue
        best_c = fill(-Inf, KNN_K)
        best_g = topk[f]  # zeros(Int, KNN_K), filled in place
        for g in 1:nf
            (g == f || G[g, g] <= 0) && continue
            c = G[f, g] / (norm_f * sqrt(G[g, g]))
            # Insert (c, g) into the sorted-descending top-K if c beats the worst kept.
            c <= best_c[KNN_K] && continue
            pos = KNN_K
            while pos > 1 && c > best_c[pos-1]
                best_c[pos] = best_c[pos-1]; best_g[pos] = best_g[pos-1]
                pos -= 1
            end
            best_c[pos] = c; best_g[pos] = g
        end
    end
    return topk, dicts
end

"""
    write_two_round_feature_columns!(file_paths)

Compute `KNN_COL` (mean of same-precursor round-1 scores in the top-`KNN_K` most-similar
other runs) and `delta_irt` from the `.pass1_sidecar.arrow` files, and write them as
columns into the per-file fold Arrow files, in place. Must run after round-1
`train_and_predict_pass1_oom!` and before round-2.

Fallback: if fewer than K other runs exist, the mean is taken over the neighbors that
do exist. If none exist (single-file experiment), the feature is 0.
"""
function write_two_round_feature_columns!(file_paths::Vector{String})
    nf = length(file_paths)
    file_pid = Vector{Vector{UInt32}}(undef, nf)
    file_s1  = Vector{Vector{Float32}}(undef, nf)
    file_irt = Vector{Vector{Float32}}(undef, nf)
    file_tgt = Vector{Vector{Bool}}(undef, nf)  # target flag per row, for PEP filter
    best_s1  = Dict{UInt32,Float32}()   # precursor -> highest s1 seen (for ref_irt)
    ref_irt  = Dict{UInt32,Float32}()   # precursor -> irt_obs at its highest-s1 instance

    # Pass A: read slim (precursor_idx, irt_obs, s1, target) for every file; track ref_irt.
    for f in 1:nf
        p = file_paths[f]
        tbl = Arrow.Table(p)
        side = Arrow.Table(p * PASS1_SIDECAR_SUFFIX)
        pid = collect(UInt32.(tbl.precursor_idx))
        irt = collect(Float32.(tbl.irt_obs))
        tgt = collect(Bool.(tbl.target))
        s1  = collect(Float32.(side.trace_prob_prepass))
        length(s1) == length(pid) ||
            error("write_two_round_feature_columns!: sidecar row mismatch at $p")
        file_pid[f] = pid; file_s1[f] = s1; file_irt[f] = irt; file_tgt[f] = tgt
        @inbounds for i in eachindex(pid)
            if s1[i] > get(best_s1, pid[i], -1.0f0)
                best_s1[pid[i]] = s1[i]; ref_irt[pid[i]] = irt[i]
            end
        end
    end

    # PEP filter: zero out per-row s1 where the global-binned local PEP exceeds
    # PEP_FILTER_THRESHOLD. This drops the row from BOTH the cosine similarity
    # profile (absent = 0) AND the neighbor lookup mean (get() returns 0).
    if PEP_FILTER_THRESHOLD < 1.0
        n_filt, n_kept = _filter_by_local_pep!(file_s1, file_tgt,
                                              PEP_FILTER_THRESHOLD, PEP_FILTER_N_BINS)
        pct = n_filt + n_kept > 0 ? round(100 * n_filt / (n_filt + n_kept), digits=1) : 0.0
        @user_info "two-round PEP filter: $(PEP_FILTER_THRESHOLD) threshold zeroed $(n_filt) of $(n_filt + n_kept) rows ($(pct)%); $(n_kept) rows kept"
    end

    topk, dicts = _compute_twin_runs(file_pid, file_s1)

    # Pass B: per file, compute the two features and rewrite the Arrow with them.
    scratch = Vector{Float32}(undef, KNN_K)   # reused per row for the median path
    for f in 1:nf
        pid = file_pid[f]; irt = file_irt[f]; n = length(pid)
        km = Vector{Float32}(undef, n); di = Vector{Float32}(undef, n)
        # Resolve up to KNN_K neighbor dicts (contiguous from front since sorted-descending;
        # a 0 slot means no more valid neighbors).
        ndicts = Dict{UInt32,Float32}[]
        for g in topk[f]
            g == 0 && break
            push!(ndicts, dicts[g])
        end
        nvalid = length(ndicts)
        if KNN_AGG === :mean
            inv_nvalid = nvalid > 0 ? 1.0f0 / Float32(nvalid) : 0.0f0
            @inbounds for i in 1:n
                s = 0.0f0
                for d in ndicts
                    s += get(d, pid[i], 0.0f0)
                end
                km[i] = s * inv_nvalid  # 0 if nvalid == 0
                di[i] = abs(irt[i] - ref_irt[pid[i]])
            end
        elseif KNN_AGG === :median
            @inbounds for i in 1:n
                for k in 1:nvalid
                    scratch[k] = get(ndicts[k], pid[i], 0.0f0)
                end
                km[i] = nvalid > 0 ? _quantile_first_n!(scratch, nvalid, 0.5f0) : 0.0f0
                di[i] = abs(irt[i] - ref_irt[pid[i]])
            end
        elseif KNN_AGG === :quantile
            @inbounds for i in 1:n
                for k in 1:nvalid
                    scratch[k] = get(ndicts[k], pid[i], 0.0f0)
                end
                km[i] = nvalid > 0 ? _quantile_first_n!(scratch, nvalid, KNN_Q) : 0.0f0
                di[i] = abs(irt[i] - ref_irt[pid[i]])
            end
        else
            error("KNN_AGG must be :mean, :median, or :quantile (got $KNN_AGG)")
        end
        main = DataFrame(Tables.columntable(Arrow.Table(file_paths[f])))
        nrow(main) == n ||
            error("write_two_round_feature_columns!: arrow row count changed at $(file_paths[f])")
        main[!, KNN_COL]    = km
        main[!, :delta_irt] = di
        writeArrow(file_paths[f], main)
    end
    n_full    = count(nb -> all(>(0), nb), topk)
    n_partial = count(nb -> nb[1] > 0 && !all(>(0), nb), topk)
    agg_desc = KNN_AGG === :quantile ? "$(KNN_Q)-quantile" : String(KNN_AGG)
    @user_info "two-round: wrote $(KNN_COL) ($(agg_desc) of top-$(KNN_K)) + delta_irt to $nf files ($n_full with all $KNN_K neighbors, $n_partial with fewer, cosine s1-profile)"
    return nothing
end
