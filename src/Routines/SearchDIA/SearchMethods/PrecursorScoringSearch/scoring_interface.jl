# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
Scoring interface for precursor probability aggregation and dictionary construction.
"""

#==========================================================
Per-File Precursor Aggregation
==========================================================#

"""
    _aggregate_trace_to_precursor_probs!(df)

Per-file Bayesian aggregation of trace-level → precursor-level probabilities.
Groups by (precursor_idx, ms_file_idx). Since ms_file_idx is constant within
a single file, this is effectively grouping by precursor_idx alone.
"""
function _aggregate_trace_to_precursor_probs!(df::DataFrame)
    prob_agg = p -> begin
        trace_prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
        clamp(trace_prob, eps(Float32), 1.0f0 - eps(Float32))
    end
    transform!(groupby(df, [:precursor_idx, :ms_file_idx]),
               :trace_prob => prob_agg => :prec_prob)
    return df
end

"""
    aggregate_per_file!(refs)

Per-file precursor probability aggregation (no MBR filtering).
"""
function aggregate_per_file!(refs::Vector{PSMFileReference})
    for ref in refs
        df = load_with_sidecars(ref)
        _aggregate_trace_to_precursor_probs!(df)
        # Write only :prec_prob as a row-aligned sidecar instead of rewriting
        # the entire main file. Downstream readers locate :prec_prob via the
        # PSMFileReference's sidecar registry.
        side_path = file_path(ref) * ".prec_prob.sidecar.arrow"
        writeArrow(side_path, DataFrame(prec_prob = df.prec_prob))
        register_sidecar!(ref, side_path, [:prec_prob])
    end
    return nothing
end

#==========================================================
Additional Interface Functions
==========================================================#

"""
    logodds(probs::AbstractVector{T}, top_n::Int) where {T<:AbstractFloat}

Combine probabilities using a log-odds average.
The final value is converted back to a probability via the logistic function.
"""
function logodds(probs::AbstractVector{T}, top_n::Int) where {T<:AbstractFloat}
    isempty(probs) && return 0.0f0
    n = min(length(probs), top_n)
    # Sort descending and select the top n probabilities
    sorted = sort(probs; rev=true)
    selected = sorted[1:n]
    eps = 1f-6
    # Convert to log-odds, clip to avoid Inf or negative contribution
    logodds = log.(clamp.(selected, 0.1f0, 1 - eps) ./ (1 .- clamp.(selected, 0.1f0, 1 - eps)))
    avg = sum(logodds) / n
    return 1.0f0 / (1 + exp(-avg))
end

#==========================================================
Dictionary + Sidecar Helper Functions for OOM Scoring Pipeline
==========================================================#

"""
    build_precursor_global_prob_dicts(refs, sqrt_n_runs, n_precursors)
    → (global_prob_dict, target_dict)

Stream per-file to build precursor_idx → global_prob dictionaries.
Reads only (precursor_idx, prec_prob, target) via mmap.
"""
function build_precursor_global_prob_dicts(
    refs::Vector{PSMFileReference},
    sqrt_n_runs::Int,
    n_precursors::Int
)
    # Pre-allocate accumulation dictionaries with known upper bound
    prob_acc = Dict{UInt32, Vector{Float32}}()
    sizehint!(prob_acc, n_precursors)
    target_dict = Dict{UInt32, Bool}()
    sizehint!(target_dict, n_precursors)

    for ref in refs
        # prec_prob lives in a sidecar after aggregate_per_file!; pull it
        # from wherever it's registered.
        cols_df = materialize_columns(ref, [:precursor_idx, :prec_prob, :target])
        n = nrow(cols_df)
        n == 0 && continue
        prec_ids = cols_df.precursor_idx
        prec_probs = cols_df.prec_prob
        targets = cols_df.target

        @inbounds for i in 1:n
            pid = prec_ids[i]
            if !haskey(prob_acc, pid)
                prob_acc[pid] = Float32[]
                target_dict[pid] = targets[i]
            end
            push!(prob_acc[pid], prec_probs[i])
        end
    end

    # Compute logodds per precursor
    global_prob_dict = Dict{UInt32, Float32}()
    sizehint!(global_prob_dict, length(prob_acc))
    for (pid, probs) in prob_acc
        global_prob_dict[pid] = logodds(probs, sqrt_n_runs)
    end

    return global_prob_dict, target_dict
end

# Learned global model (ENV["GLOBAL_MODEL"]=1): replace the fixed top-√N log-odds with a
# precursor-level LightGBM on round-1 OOF-score aggregates, 2-fold precursor-keyed CV.
global_model_enabled() = get(ENV, "GLOBAL_MODEL", "") == "1"
const _GLOBAL_MODEL_NFEAT = 10   # n, top1, top2, top3, min, mean, std, frac>.99, range, logodds
const _GLOBAL_MODEL_TOPK = 8     # running top-K per precursor (covers top-3 + logodds for ≤64 runs)

@inline function _logit_clamp(x::Float32)
    c = clamp(x, 0.1f0, 1.0f0 - 1.0f-6); return log(c / (1.0f0 - c))
end

# 2-fold entity-keyed CV for the learned global models: train on one fold's rows, score the other
# (clean OOF since each entity — precursor or protein — is in exactly one fold). Falls back to the
# log-odds feature (column `logodds_col`) for thin/one-class folds or entities with no fold (-1).
function _fit_global_model_2fold!(X::Matrix{Float32}, y::Vector{Bool}, foldv::Vector{Int8},
                                  scores::Vector{Float32}, logodds_col::Int)
    fill!(scores, 0.0f0)
    idx0 = findall(==(Int8(0)), foldv); idx1 = findall(==(Int8(1)), foldv)
    other = findall(f -> f != Int8(0) && f != Int8(1), foldv)   # unassigned → logodds fallback
    @inbounds for j in other; scores[j] = X[j, logodds_col]; end
    function _fs!(train_idx, test_idx)
        isempty(test_idx) && return
        if length(train_idx) < 100 || length(unique(@view y[train_idx])) < 2
            @inbounds for j in test_idx; scores[j] = X[j, logodds_col]; end
            return
        end
        cls = _fit_pass1_booster(X[train_idx, :], y[train_idx], SCORING_LGBM_HP)
        if cls isa NamedTuple && cls.kind === :constant
            @inbounds for j in test_idx; scores[j] = Float32(cls.value); end
        else
            raw = LightGBM.predict(cls, X[test_idx, :]; num_threads = Threads.nthreads())
            pr = ndims(raw) == 2 ? dropdims(raw; dims = 2) : raw
            @inbounds for (t, j) in enumerate(test_idx); scores[j] = Float32(pr[t]); end
        end
    end
    _fs!(idx1, idx0)   # train fold-1 → score fold-0
    _fs!(idx0, idx1)   # train fold-0 → score fold-1
    return scores
end

# Same 10 aggregate features as the precursor global model, computed from a per-entity score VECTOR
# (used for proteins, which are few enough that a vector + sort is cheap). Cols must match:
# [n, top1, top2, top3, min, mean, std, frac>.99, range, logodds].
function _global_model_feats_vec!(out::Vector{Float32}, v::Vector{Float32}, sqrt_n_runs::Int)
    n = length(v); K = _GLOBAL_MODEL_TOPK
    mx = v[1]; mn = v[1]; s = 0.0f0; s2 = 0.0f0; c99 = 0
    @inbounds for x in v
        x > mx && (mx = x); x < mn && (mn = x)
        s += x; s2 += x * x; x > 0.99f0 && (c99 += 1)
    end
    sv = sort(v; rev = true)
    klog = min(sqrt_n_runs, K, n); lo = 0.0f0
    @inbounds for j in 1:klog; lo += _logit_clamp(sv[j]); end
    lo /= klog
    nf = Float32(n)
    out[1]  = nf
    out[2]  = sv[1]
    out[3]  = n >= 2 ? sv[2] : sv[1]
    out[4]  = n >= 3 ? sv[3] : out[3]
    out[5]  = mn
    out[6]  = s / nf
    out[7]  = n > 1 ? sqrt(max(0.0f0, (s2 - s * s / nf) / (nf - 1.0f0))) : 0.0f0
    out[8]  = Float32(c99) / nf
    out[9]  = mx - mn
    out[10] = lo
    return out
end

"""
    build_precursor_global_model_dict(refs, sqrt_n_runs, n_precursors) → (global_prob_dict, target_dict)

Learned replacement for the fixed top-√N log-odds global score. STREAMING + memory-efficient: one
pass over the refs into fixed-size per-precursor accumulators (indexed by precursor id, bounded by
the library — independent of #runs/#PSMs; no per-precursor vector, no sort). Features per precursor
from the round-2 per-run precursor score (`:prec_prob`, same input the fixed log-odds aggregates):
n_runs, top-1/2/3, min, mean, std, frac>0.99, range, and the top-√N log-odds itself (as a feature, so
it can't do worse). Trains a precursor-level LightGBM with 2-fold precursor-keyed CV (train one fold,
score the other — clean OOF since each precursor is in one fold); falls back to the log-odds feature
on thin/one-class folds.
"""
function build_precursor_global_model_dict(
    refs::Vector{PSMFileReference}, sqrt_n_runs::Int, n_precursors::Int,
)
    P = n_precursors
    K = _GLOBAL_MODEL_TOPK
    cnt  = zeros(Int32, P + 1); c99 = zeros(Int32, P + 1)
    ssum = zeros(Float32, P + 1); ssq = zeros(Float32, P + 1)
    mn   = fill(2.0f0, P + 1)
    topM = zeros(Float32, K, P + 1)          # per-precursor running top-K (descending)
    fold = fill(Int8(-1), P + 1); tgt = falses(P + 1)

    # ---- streaming pass ----
    for ref in refs
        cols = materialize_columns(ref, [:precursor_idx, :prec_prob, :cv_fold, :target])
        n = nrow(cols); n == 0 && continue
        pid = cols.precursor_idx; s1 = cols.prec_prob; cf = cols.cv_fold; tg = cols.target
        @inbounds for i in 1:n
            p = Int(pid[i]) + 1; x = Float32(s1[i])
            fold[p] == Int8(-1) && (fold[p] = Int8(cf[i]); tgt[p] = Bool(tg[i]))
            cnt[p] += Int32(1); ssum[p] += x; ssq[p] += x * x
            x > 0.99f0 && (c99[p] += Int32(1))
            x < mn[p] && (mn[p] = x)
            if x > topM[K, p]                # top-K insert (early exit for the common low score)
                j = K
                while j > 1 && x > topM[j - 1, p]
                    topM[j, p] = topM[j - 1, p]; j -= 1
                end
                topM[j, p] = x
            end
        end
    end

    # ---- build feature matrix over observed precursors ----
    obs = findall(>(Int32(0)), cnt)          # array indices (p = pid + 1)
    nprec = length(obs)
    X = Matrix{Float32}(undef, nprec, _GLOBAL_MODEL_NFEAT)
    y = Vector{Bool}(undef, nprec); foldv = Vector{Int8}(undef, nprec)
    pids = Vector{UInt32}(undef, nprec)
    @inbounds for (i, p) in enumerate(obs)
        pids[i] = UInt32(p - 1)
        nf = Float32(cnt[p]); ni = Int(cnt[p])
        mean_v = ssum[p] / nf
        var_v  = ni > 1 ? max(0.0f0, (ssq[p] - ssum[p]^2 / nf) / (nf - 1.0f0)) : 0.0f0
        klog = min(sqrt_n_runs, K, ni)       # top-√N log-odds (mean of top logits), capped at K
        lo = 0.0f0
        for j in 1:klog; lo += _logit_clamp(topM[j, p]); end
        lo /= klog
        X[i, 1]  = nf
        X[i, 2]  = topM[1, p]
        X[i, 3]  = topM[2, p]
        X[i, 4]  = topM[3, p]
        X[i, 5]  = mn[p]
        X[i, 6]  = mean_v
        X[i, 7]  = sqrt(var_v)
        X[i, 8]  = Float32(c99[p]) / nf
        X[i, 9]  = topM[1, p] - mn[p]
        X[i, 10] = lo
        y[i] = tgt[p]; foldv[i] = fold[p]
    end

    # ---- 2-fold precursor-keyed CV: train one fold, score the other ----
    scores = Vector{Float32}(undef, nprec)
    _fit_global_model_2fold!(X, y, foldv, scores, 10)   # logodds is feature column 10

    global_prob_dict = Dict{UInt32, Float32}(); sizehint!(global_prob_dict, nprec)
    target_dict = Dict{UInt32, Bool}(); sizehint!(target_dict, nprec)
    @inbounds for i in 1:nprec
        global_prob_dict[pids[i]] = scores[i]; target_dict[pids[i]] = y[i]
    end
    return global_prob_dict, target_dict
end

"""
    build_global_qval_dict_from_scores(score_dict, target_dict, fdr_scale) → Dict{UInt32, Float32}

Compute global q-values from a score dictionary without any file I/O.
"""
function build_global_qval_dict_from_scores(
    score_dict::Dict{UInt32, Float32},
    target_dict::Dict{UInt32, Bool},
    fdr_scale::Float32
)
    n = length(score_dict)
    pids = collect(keys(score_dict))
    scores = Float32[score_dict[pid] for pid in pids]
    targets = Bool[get(target_dict, pid, false) for pid in pids]

    # Sort descending by score
    perm = sortperm(scores; rev=true)
    permute!(pids, perm)
    permute!(scores, perm)
    permute!(targets, perm)

    # Compute q-values
    qvals = Vector{Float32}(undef, n)
    get_qvalues!(scores, targets, qvals; fdr_scale_factor=fdr_scale)

    # Build dictionary
    qval_dict = Dict{UInt32, Float32}()
    sizehint!(qval_dict, n)
    for i in 1:n
        qval_dict[pids[i]] = qvals[i]
    end
    return qval_dict
end

"""
    build_global_pep_dict_from_scores(score_dict, target_dict, fdr_scale) → Dict{UInt32, Float32}

Compute global posterior error probabilities (local FDR) from a score dictionary
without any file I/O. Parallel to `build_global_qval_dict_from_scores` but uses
`get_PEP!` (PAVA-fit) instead of cumulative q-values.
"""
function build_global_pep_dict_from_scores(
    score_dict::Dict{UInt32, Float32},
    target_dict::Dict{UInt32, Bool},
    fdr_scale::Float32
)
    n = length(score_dict)
    pids = collect(keys(score_dict))
    scores = Float32[score_dict[pid] for pid in pids]
    targets = Bool[get(target_dict, pid, false) for pid in pids]

    peps = Vector{Float32}(undef, n)
    get_PEP!(scores, targets, peps; doSort=true, fdr_scale_factor=fdr_scale)

    pep_dict = Dict{UInt32, Float32}()
    sizehint!(pep_dict, n)
    for i in 1:n
        pep_dict[pids[i]] = peps[i]
    end
    return pep_dict
end


"""
    write_score_sidecars(refs, columns; temp_prefix) → Vector{PSMFileReference}

Extract only the named columns from each file into a temporary Arrow sidecar file.
"""
function write_score_sidecars(
    refs::Vector{<:FileReference},
    columns::Vector{Symbol};
    temp_prefix::String = "sidecar"
)
    sidecar_refs = PSMFileReference[]
    for ref in refs
        # Use materialize_columns so columns are pulled from main OR any
        # registered sidecar (e.g. :prec_prob now lives in a sidecar after
        # aggregate_per_file!).
        df = ref isa PSMFileReference ? materialize_columns(ref, columns) :
             DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))[!, columns]
        nrow(df) == 0 && continue
        temp_path = tempname() * "_$(temp_prefix).arrow"
        writeArrow(temp_path, df)
        push!(sidecar_refs, PSMFileReference(temp_path))
    end
    return sidecar_refs
end

"""
    build_qvalue_spline_from_refs(refs, score_col, merged_path; ...) → Union{Nothing, NamedTuple}

Encapsulates the full sidecar lifecycle: write → sort → merge → cleanup → spline computation.
"""
function build_qvalue_spline_from_refs(
    refs::Vector{<:FileReference},
    score_col::Symbol,
    merged_path::String;
    batch_size::Int = 10_000_000,
    compute_pep::Bool = false,
    min_pep_points_per_bin::Int = 100,
    fdr_scale_factor::Float32 = 1.0f0,
    temp_prefix::String = "sidecar"
)
    sidecar_refs = write_score_sidecars(refs, [score_col, :target]; temp_prefix=temp_prefix)
    isempty(sidecar_refs) && return nothing

    try
        sort_file_by_keys!(sidecar_refs, score_col, :target; reverse=[true, true])
        stream_sorted_merge(sidecar_refs, merged_path, score_col, :target;
                           batch_size=batch_size, reverse=[true, true])
    finally
        GC.gc(false)
        for ref in sidecar_refs
            safeRm(file_path(ref), nothing; force=true)
        end
    end

    qval_spline = get_qvalue_spline(merged_path, score_col, false;
        min_pep_points_per_bin=min_pep_points_per_bin,
        fdr_scale_factor=fdr_scale_factor)

    pep_interp = if compute_pep
        get_pep_interpolation(merged_path, score_col;
            fdr_scale_factor=fdr_scale_factor)
    else
        nothing
    end

    return (; qval_spline, pep_interp)
end
