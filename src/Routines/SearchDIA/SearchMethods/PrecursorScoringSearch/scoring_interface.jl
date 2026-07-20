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
const _GLOBAL_MODEL_NFEAT = 10

# Per-precursor aggregate features of the round-1 OOF scores `v` (across the runs it was seen in).
function _global_model_feats!(out::Vector{Float32}, v::Vector{Float32}, sqrt_n_runs::Int)
    n = length(v)
    mx = v[1]; mn = v[1]; s = 0.0f0; s2 = 0.0f0; c99 = 0
    @inbounds for x in v
        x > mx && (mx = x); x < mn && (mn = x)
        s += x; s2 += x * x; x > 0.99f0 && (c99 += 1)
    end
    mean_v = s / n
    var_v  = n > 1 ? max(0.0f0, (s2 - s * s / n) / (n - 1)) : 0.0f0
    sv = sort(v; rev = true)
    med = isodd(n) ? sv[(n + 1) ÷ 2] : 0.5f0 * (sv[n ÷ 2] + sv[n ÷ 2 + 1])
    out[1]  = Float32(n)
    out[2]  = mx
    out[3]  = mean_v
    out[4]  = mn
    out[5]  = sqrt(var_v)
    out[6]  = med
    out[7]  = logodds(v, sqrt_n_runs)      # the current fixed global score, as a feature (can't do worse)
    out[8]  = Float32(c99) / n
    out[9]  = mx - mn
    out[10] = n >= 2 ? sv[2] : sv[1]
    return out
end

"""
    build_precursor_global_model_dict(refs, sqrt_n_runs, n_precursors) → (global_prob_dict, target_dict)

Learned replacement for the fixed top-√N log-odds global score. Aggregates each precursor's ROUND-1
OOF scores (`:s1_round1`) across runs into features, then trains a precursor-level LightGBM with
2-fold precursor-keyed CV (train one fold, score the other — clean OOF because each precursor is in
exactly one fold). Falls back to the log-odds feature when a fold has too little/one-class data.
Requires the `:s1_round1` column (written by the two-round group-feature step).
"""
function build_precursor_global_model_dict(
    refs::Vector{PSMFileReference}, sqrt_n_runs::Int, n_precursors::Int,
)
    s1_acc = Dict{UInt32, Vector{Float32}}(); sizehint!(s1_acc, n_precursors)
    fold_of = Dict{UInt32, UInt8}(); sizehint!(fold_of, n_precursors)
    target_dict = Dict{UInt32, Bool}(); sizehint!(target_dict, n_precursors)
    for ref in refs
        cols = materialize_columns(ref, [:precursor_idx, :s1_round1, :cv_fold, :target])
        n = nrow(cols); n == 0 && continue
        pid = cols.precursor_idx; s1 = cols.s1_round1; cf = cols.cv_fold; tg = cols.target
        @inbounds for i in 1:n
            p = pid[i]
            if !haskey(s1_acc, p)
                s1_acc[p] = Float32[]; fold_of[p] = UInt8(cf[i]); target_dict[p] = tg[i]
            end
            push!(s1_acc[p], Float32(s1[i]))
        end
    end
    pids = collect(keys(s1_acc)); nprec = length(pids)
    X = Matrix{Float32}(undef, nprec, _GLOBAL_MODEL_NFEAT)
    y = Vector{Bool}(undef, nprec); fold = Vector{UInt8}(undef, nprec)
    feat = Vector{Float32}(undef, _GLOBAL_MODEL_NFEAT)
    @inbounds for i in 1:nprec
        _global_model_feats!(feat, s1_acc[pids[i]], sqrt_n_runs)
        for j in 1:_GLOBAL_MODEL_NFEAT; X[i, j] = feat[j]; end
        y[i] = target_dict[pids[i]]; fold[i] = fold_of[pids[i]]
    end
    scores = Vector{Float32}(undef, nprec)
    idx0 = findall(==(UInt8(0)), fold); idx1 = findall(==(UInt8(1)), fold)
    function _fit_score!(train_idx, test_idx)
        isempty(test_idx) && return
        if length(train_idx) < 100 || length(unique(@view y[train_idx])) < 2
            @inbounds for j in test_idx; scores[j] = X[j, 7]; end     # fallback: logodds feature
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
    _fit_score!(idx1, idx0)   # train fold-1 → score fold-0
    _fit_score!(idx0, idx1)   # train fold-0 → score fold-1
    global_prob_dict = Dict{UInt32, Float32}(); sizehint!(global_prob_dict, nprec)
    @inbounds for i in 1:nprec; global_prob_dict[pids[i]] = scores[i]; end
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
