# Two-round experiment-wide scoring: cross-run consistency features.
# Design + CV protocol: see TWO_ROUND_SCORING.md.
#
# Runs when `match_between_runs` is true (the second round). After the round-1 Pass-1 OOM trainer has
# written each file's `.pass1_sidecar.arrow` (OOF s1 = trace_prob_prepass), `write_two_round_feature_columns!`
# clusters the runs per CV fold and derives the GROUP cross-run features (global/cluster
# max, top3_logodds, n_present, n_confident) + delta_irt from s1, writing them as columns into each
# per-file fold Arrow. The round-2 `train_and_predict_pass1_oom!` pass then builds symmetric shadow
# rows INSIDE its training sample (each grafted with the original row's cross-run features) and trains
# on [ ADVANCED_FEATURE_SET ; GROUP_COLS ; delta_irt ] with a 1:1 target/decoy marginal on the grafts.
# Shadows never touch disk, so the per-file Arrows keep their schema and row count throughout.
#
# Leakage-safety rests on the existing precursor_idx-keyed CV folds (Invariant A): every PSM of a
# precursor shares one fold, so a PSM's cross-run lookups land in the same held-out fold and their s1
# is OOF-consistent; run clustering is done independently per fold. The round-2 trainer reuses the
# unchanged `cv_fold` column (Invariant B).
#
# Assumes one row per precursor_idx per file (current MainSearch invariant), so
# s1[file, precursor] is a single value — no trace->precursor aggregation needed.

# Whether the second (cross-run) round runs is gated by `match_between_runs` in
# `score_precursor_isotope_traces` — there is no separate two-round toggle.

# Symmetric shadow-decoys: for every sampled real row the trainer builds a paired shadow of the
# opposite class carrying the source's other features but the ORIGINAL row's grafted cross-run
# feature(s). See TWO_ROUND_SCORING.md.
# Positive logit: 0 for s1<=0.5, else log(s1/(1-s1)) clamped away from the asymptote.
@inline function _pos_logit(x::Float32)
    x <= 0.5f0 && return 0.0f0
    c = min(x, 1.0f0 - 1.0f-6)
    return log(c / (1.0f0 - c))
end

# GROUP cross-run features: scalable global + hard-cluster consistency, clustered PER CV FOLD.
# Eight grafted features — global_{max,top3_logodds,n_present,n_confident} over ALL other runs +
# cluster_{max,top3_logodds,n_present,n_confident} over the run's fold-specific cluster. Built in
# one O(N) per-fold pass via top-4 (score,run) records with leave-one-out (see
# _write_group_features!). Cluster size follows _group_size(R): OFF for R<9 (cluster cols mirror
# global), else min(9, floor(R/2)).
const GROUP_COLS = [:global_max, :global_top3_logodds, :cluster_max, :cluster_top3_logodds,
                    :global_n_present, :global_n_confident, :cluster_n_present, :cluster_n_confident]

# Round-2 feature columns appended to ADVANCED_FEATURE_SET: the GROUP cross-run features + delta_irt.
two_round_features() = vcat(GROUP_COLS, [:delta_irt])

# Cross-run feature columns grafted onto each shadow twin so every one keeps a 1:1 target/decoy
# marginal — forcing round-2 LGBM to use them jointly with the base features, never as a standalone
# signal. delta_irt is NOT grafted: each shadow keeps its source row's own iRT self-consistency.
_mbr_graft_cols() = Tuple(c for c in GROUP_COLS)

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

# ================= Cross-run GROUP features: per-fold clustering + global/cluster =================
const GROUP_EMBED_DIM   = 16      # SVD/eigen embedding dimension for run clustering
const GROUP_PRESENCE_S1 = 0.9f0   # round-1 OOF score threshold for "present" (clustering input)
const GROUP_CONF_S1     = 0.99f0  # round-1 OOF score threshold for the n_confident count (q<0.01 proxy)

# Group-size schedule: 0 = cluster OFF (global only); else min(9, floor(R/2)). Grows as R/2 (2
# groups) until it caps at 9 at R=18, then group size stays 9 and the group COUNT grows (~R/9).
@inline _group_size(R::Int) = R < 9 ? 0 : min(9, R ÷ 2)

# Top-4 (score, run) record for one precursor, kept descending — enough for leave-one-out of top-3.
struct TopK4
    v::NTuple{4,Float32}
    r::NTuple{4,UInt32}
end
const EMPTY_TOPK4 = TopK4((0f0, 0f0, 0f0, 0f0), (UInt32(0), UInt32(0), UInt32(0), UInt32(0)))
@inline function _tk_insert(t::TopK4, s::Float32, run::UInt32)
    v = t.v; r = t.r
    s <= v[4] && return t
    if     s > v[1]; return TopK4((s, v[1], v[2], v[3]), (run, r[1], r[2], r[3]))
    elseif s > v[2]; return TopK4((v[1], s, v[2], v[3]), (r[1], run, r[2], r[3]))
    elseif s > v[3]; return TopK4((v[1], v[2], s, v[3]), (r[1], r[2], run, r[3]))
    else             return TopK4((v[1], v[2], v[3], s), (r[1], r[2], r[3], run))
    end
end
@inline _tk_loo_max(t::TopK4, own::UInt32) = @inbounds (t.r[1] != own) ? t.v[1] : t.v[2]
@inline function _tk_loo_top3(t::TopK4, own::UInt32)   # Σ positive logits of top-3 excluding own run
    acc = 0f0; taken = 0
    @inbounds for i in 1:4
        taken == 3 && break
        t.v[i] <= 0f0 && break
        t.r[i] == own && continue
        acc += _pos_logit(t.v[i]); taken += 1
    end
    return acc
end

# Randomized truncated SVD of the L2-row-normalized presence matrix M (runs × precursors), returning
# the top-`dim` left-singular run embedding scaled by singular values. O(nnz·ℓ) time and O(nnz + P·ℓ)
# memory — NO R×R Gram (which is O(R²·density) to build and O(R³) to eigendecompose). Works directly
# from the sparse presence sets (run→precursors), never materializing a dense matrix. Drops singleton
# precursors (<2 runs; no pairwise signal) and ubiquitous ones (>90% of runs; no discriminative
# signal). Halko et al. range-finder with power iterations; qr/svd act only on small (R×ℓ / ℓ×P) tiles.
function _rsvd_run_embed(run_present::Vector{Vector{UInt32}}, R::Int;
                         dim::Int = GROUP_EMBED_DIM, oversample::Int = 8, power_iters::Int = 2)
    # postings: precursor -> runs; frequency-filter to informative precursors; compress ids.
    postings = Dict{UInt32,Vector{Int}}()
    for r in 1:R, p in run_present[r]
        push!(get!(() -> Int[], postings, p), r)
    end
    maxfreq = 0.9 * R
    pres = [Int[] for _ in 1:R]            # run -> compressed kept-precursor ids
    P = 0
    for (_, runs) in postings
        L = length(runs)
        (L < 2 || L > maxfreq) && continue
        P += 1
        for r in runs; push!(pres[r], P); end
    end
    (P == 0 || R < 2) && return zeros(Float64, R, 1)
    w = Vector{Float64}(undef, R)          # L2 row-normalization weight
    @inbounds for r in 1:R; n = length(pres[r]); w[r] = n > 0 ? 1.0 / sqrt(n) : 0.0; end

    ℓ = min(dim + oversample, R, P)
    rng = MersenneTwister(0x51ED270F)
    # Y = Mw * Ω, Ω ~ N(0,1) of size P×ℓ (materialized once).
    Ω = randn(rng, P, ℓ)
    Y = zeros(Float64, R, ℓ)
    @inbounds for r in 1:R
        wr = w[r]; wr == 0.0 && continue
        for c in pres[r], j in 1:ℓ; Y[r, j] += Ω[c, j]; end
        for j in 1:ℓ; Y[r, j] *= wr; end
    end
    Z = Matrix{Float64}(undef, P, ℓ)       # reused scratch for power iterations / B
    for _ in 1:power_iters
        Q = Matrix(qr!(Y).Q)               # R×ℓ orthonormal
        fill!(Z, 0.0)                      # Z = Mwᵀ Q  (P×ℓ)
        @inbounds for r in 1:R
            wr = w[r]; wr == 0.0 && continue
            for c in pres[r], j in 1:ℓ; Z[c, j] += wr * Q[r, j]; end
        end
        fill!(Y, 0.0)                      # Y = Mw Z  (R×ℓ)
        @inbounds for r in 1:R
            wr = w[r]; wr == 0.0 && continue
            for c in pres[r], j in 1:ℓ; Y[r, j] += Z[c, j]; end
            for j in 1:ℓ; Y[r, j] *= wr; end
        end
    end
    Q = Matrix(qr!(Y).Q)                   # R×ℓ
    B = zeros(Float64, ℓ, P)               # B = Qᵀ Mw  (ℓ×P)
    @inbounds for r in 1:R
        wr = w[r]; wr == 0.0 && continue
        for c in pres[r], j in 1:ℓ; B[j, c] += wr * Q[r, j]; end
    end
    F = svd!(B)                            # B = U_B Σ V_Bᵀ ; U_B is ℓ×ℓ
    U = Q * F.U                            # R×ℓ left singular vectors of Mw
    dd = min(dim, ℓ)
    E = Matrix{Float64}(undef, R, dd)
    @inbounds for j in 1:dd, r in 1:R; E[r, j] = U[r, j] * F.S[j]; end
    return E
end

# Small Lloyd's k-means (R runs, low dim, few clusters) with fixed-seed restarts.
function _kmeans_small(E::Matrix{Float64}, ncl::Int; iters::Int = 30, restarts::Int = 5)
    R, d = size(E); ncl = clamp(ncl, 1, R)
    ncl == 1 && return ones(Int, R)
    rng = MersenneTwister(0x9E3779B9)
    best_lab = ones(Int, R); best_cost = Inf
    cent = Matrix{Float64}(undef, ncl, d); lab = Vector{Int}(undef, R); cnt = Vector{Int}(undef, ncl)
    for _ in 1:restarts
        init = randperm(rng, R)[1:ncl]
        @inbounds for c in 1:ncl, j in 1:d; cent[c, j] = E[init[c], j]; end
        cost = Inf
        for _ in 1:iters
            newcost = 0.0
            @inbounds for r in 1:R
                bd = Inf; bc = 1
                for c in 1:ncl
                    s = 0.0
                    for j in 1:d; δ = E[r, j] - cent[c, j]; s += δ * δ; end
                    s < bd && (bd = s; bc = c)
                end
                lab[r] = bc; newcost += bd
            end
            fill!(cent, 0.0); fill!(cnt, 0)
            @inbounds for r in 1:R
                c = lab[r]; cnt[c] += 1
                for j in 1:d; cent[c, j] += E[r, j]; end
            end
            @inbounds for c in 1:ncl, j in 1:d
                cnt[c] > 0 && (cent[c, j] /= cnt[c])
            end
            abs(cost - newcost) < 1e-9 && (cost = newcost; break)
            cost = newcost
        end
        cost < best_cost && (best_cost = cost; copyto!(best_lab, lab))
    end
    return best_lab
end

# Split any cluster larger than k into two (recursively), so all groups end ≤ k.
function _split_oversized!(lab::Vector{Int}, E::Matrix{Float64}, k::Int)
    changed = true
    while changed
        changed = false
        for c in unique(lab)
            idx = findall(==(c), lab)
            if length(idx) > k
                sub = _kmeans_small(E[idx, :], 2; restarts = 2)
                nxt = maximum(lab) + 1
                @inbounds for t in eachindex(idx); sub[t] == 2 && (lab[idx[t]] = nxt); end
                changed = true
            end
        end
    end
    return lab
end

# run -> cluster labels (1-based, contiguous). k<=0 -> single cluster (cluster feature == global).
function _cluster_runs(run_present::Vector{Vector{UInt32}}, R::Int, k::Int)
    k <= 0 && return ones(Int, R)
    E = _rsvd_run_embed(run_present, R)     # randomized SVD on sparse presence (no R×R Gram)
    lab = _kmeans_small(E, max(1, round(Int, R / k)))
    _split_oversized!(lab, E, k)
    u = unique(lab); remap = Dict(u[i] => i for i in eachindex(u))
    return [remap[l] for l in lab]
end

# The two-round path receives FOLD-SPLIT files (`*_fold{0,1}.arrow`): one Arrow per (run, CV fold),
# each single-fold. So the run count for a fold = the number of that fold's files, and clustering must
# see only that fold's files. `_file_fold` reads a file's fold from its (constant) cv_fold column.
@inline function _file_fold(fold_col)::Int
    isempty(fold_col) ? -1 : Int(fold_col[1])
end

"""
    _write_group_features!(file_paths, max_pid)

GROUP_SCORING path. Files are fold-split (one per run×fold). Per CV fold, over ONLY that fold's files
(= the runs): cluster the runs from their own OOF presence, then build per-precursor top-4 (score,run)
records — `grec` over all the fold's runs, `crec` per cluster (reset via `touched`) — and write
leave-one-out `global_{max,top3_logodds}`, `cluster_{max,top3_logodds}`, and `delta_irt`.
Run size follows `_group_size(R)`; R<9 → cluster OFF (cluster cols mirror global).

STREAMING: peak memory is O(`max_pid`) accumulators + the fold's sparse presence + one cluster's
(≤`k` ≤ 9) files — independent of the number of runs. Files are re-read per pass (fold + emit)
instead of the whole experiment being held in RAM, so this scales to hundreds of runs. `max_pid` is
the library precursor count (an upper bound on `precursor_idx`); accumulators are sized once from it.
Bit-identical to the prior all-in-memory version (same read order, same insertion/emit order).
"""
function _write_group_features!(file_paths::Vector{String}, max_pid::Int)
    nfiles = length(file_paths)

    # Cheap metadata pass: each file's CV fold (constant column, one value read).
    fold_of = Vector{Int}(undef, nfiles)
    for f in 1:nfiles
        fold_of[f] = _file_fold(Arrow.Table(file_paths[f]).cv_fold)
    end

    # Read one file's (precursor_idx, s1, irt_obs). s1 = round-1 OOF score from the Pass-1 sidecar.
    read_file = function(f::Int)
        tbl  = Arrow.Table(file_paths[f])
        side = Arrow.Table(file_paths[f] * PASS1_SIDECAR_SUFFIX)
        pid  = collect(UInt32.(tbl.precursor_idx))
        s1   = collect(Float32.(side.trace_prob_prepass))
        length(s1) == length(pid) ||
            error("_write_group_features!: sidecar row mismatch at $(file_paths[f])")
        irt  = collect(Float32.(tbl.irt_obs))
        return pid, s1, irt
    end

    grec    = fill(EMPTY_TOPK4, max_pid + 1)        # global records (reused per fold)
    crec    = fill(EMPTY_TOPK4, max_pid + 1)        # cluster records (reused per cluster via touched)
    gcnt_p  = zeros(Int32, max_pid + 1)             # global count: runs present (reused per fold)
    gcnt_c  = zeros(Int32, max_pid + 1)             # global count: runs with s1>=GROUP_CONF_S1
    ccnt_p  = zeros(Int32, max_pid + 1)             # cluster count: present (reset per cluster via touched)
    ccnt_c  = zeros(Int32, max_pid + 1)             # cluster count: confident
    best_s1 = fill(-1.0f0, max_pid + 1)             # for ref_irt (highest-s1 instance per precursor)
    ref_irt = zeros(Float32, max_pid + 1)
    touched = UInt32[]
    kseen = -1

    for fold in 0:1
        files = [f for f in 1:nfiles if fold_of[f] == fold]   # this fold's files == its runs
        R = length(files); R == 0 && continue
        k = _group_size(R); kseen = k

        # Pass 1 (stream one file at a time): global records + counts + ref_irt, and the fold's
        # per-run presence for clustering.
        fill!(grec, EMPTY_TOPK4); fill!(best_s1, -1.0f0); fill!(gcnt_p, 0); fill!(gcnt_c, 0)
        run_present = Vector{Vector{UInt32}}(undef, R)
        for ri in 1:R
            pid, s1, irt = read_file(files[ri])
            pres = UInt32[]
            @inbounds for i in eachindex(pid)
                p = pid[i]; s = s1[i]
                grec[p + 1] = _tk_insert(grec[p + 1], s, UInt32(ri))
                gcnt_p[p + 1] += Int32(1)
                s >= GROUP_CONF_S1 && (gcnt_c[p + 1] += Int32(1))
                if s > best_s1[p + 1]; best_s1[p + 1] = s; ref_irt[p + 1] = irt[i]; end
                s >= GROUP_PRESENCE_S1 && push!(pres, p)
            end
            run_present[ri] = unique!(sort!(pres))
        end

        labels = _cluster_runs(run_present, R, k)   # length R; run id = index into `files`
        ncl = maximum(labels)
        runs_in = [Int[] for _ in 1:ncl]
        for ri in 1:R; push!(runs_in[labels[ri]], ri); end
        run_present = nothing   # release sparse presence

        # Pass 2: per cluster (≤k runs), hold just that cluster's files, build cluster records,
        # then emit LOO features for its runs' rows and write each file's GROUP columns in place.
        for c in 1:ncl
            cl_ris = runs_in[c]
            buf = Dict{Int, Tuple{Vector{UInt32}, Vector{Float32}, Vector{Float32}}}()
            empty!(touched)
            for ri in cl_ris
                pid, s1, irt = read_file(files[ri])
                buf[ri] = (pid, s1, irt)
                @inbounds for i in eachindex(pid)
                    p = pid[i]; s = s1[i]
                    crec[p + 1] === EMPTY_TOPK4 && push!(touched, p)
                    crec[p + 1] = _tk_insert(crec[p + 1], s, UInt32(ri))
                    ccnt_p[p + 1] += Int32(1)
                    s >= GROUP_CONF_S1 && (ccnt_c[p + 1] += Int32(1))
                end
            end
            for ri in cl_ris
                f = files[ri]; ru = UInt32(ri)
                pid, s1, irt = buf[ri]
                n = length(pid)
                gmax = Vector{Float32}(undef, n); gt3 = Vector{Float32}(undef, n)
                cmax = Vector{Float32}(undef, n); ct3 = Vector{Float32}(undef, n)
                gnp  = Vector{Float32}(undef, n); gnc = Vector{Float32}(undef, n)
                cnp  = Vector{Float32}(undef, n); cnc = Vector{Float32}(undef, n)
                di   = Vector{Float32}(undef, n)
                @inbounds for i in 1:n
                    p = pid[i]; s = s1[i]
                    ownc = s >= GROUP_CONF_S1 ? Int32(1) : Int32(0)
                    gmax[i] = _tk_loo_max(grec[p + 1], ru)
                    gt3[i]  = _tk_loo_top3(grec[p + 1], ru)
                    gnp[i]  = Float32(gcnt_p[p + 1] - Int32(1))       # exclude own row
                    gnc[i]  = Float32(gcnt_c[p + 1] - ownc)
                    if k <= 0
                        cmax[i] = gmax[i]; ct3[i] = gt3[i]
                        cnp[i]  = gnp[i];  cnc[i] = gnc[i]
                    else
                        cmax[i] = _tk_loo_max(crec[p + 1], ru)
                        ct3[i]  = _tk_loo_top3(crec[p + 1], ru)
                        cnp[i]  = Float32(ccnt_p[p + 1] - Int32(1))
                        cnc[i]  = Float32(ccnt_c[p + 1] - ownc)
                    end
                    di[i] = abs(irt[i] - ref_irt[p + 1])
                end
                # Write the 9 GROUP columns to a row-aligned sidecar instead of reading and
                # rewriting the whole ~80-column fold file. The round-2 read sites resolve
                # them via `_feature_columns` (pass1_oom.jl), and Step 1b's PSMFileReference
                # auto-discovers `.group.sidecar.arrow` so they still reach the merged output.
                ref = PSMFileReference(file_paths[f])
                ref.row_count == n ||
                    error("_write_group_features!: arrow row count changed at $(file_paths[f])")
                add_columns_via_sidecar!(ref,
                    :global_max           => gmax, :global_top3_logodds  => gt3,
                    :cluster_max          => cmax, :cluster_top3_logodds => ct3,
                    :global_n_present     => gnp,  :global_n_confident   => gnc,
                    :cluster_n_present    => cnp,  :cluster_n_confident  => cnc,
                    :delta_irt            => di;
                    tag = "group")
            end
            for p in touched; crec[p + 1] = EMPTY_TOPK4; ccnt_p[p + 1] = Int32(0); ccnt_c[p + 1] = Int32(0); end
        end
    end
    nruns = count(!=(-1), fold_of) ÷ 2
    @user_info "group-scoring (streaming): wrote global+cluster (max, top3_logodds, n_present, n_confident) + delta_irt for ~$nruns runs/fold ($nfiles fold-files, k=$kseen)"
    return nothing
end

"""
    write_two_round_feature_columns!(file_paths, max_pid)

Round-2 cross-run features. After round-1 `train_and_predict_pass1_oom!` has written each file's
`.pass1_sidecar.arrow` (OOF s1 = `trace_prob_prepass`), cluster the runs per CV fold and write the
GROUP cross-run features + `delta_irt` to each fold file's row-aligned
`.group.sidecar.arrow` (NOT back into the fold Arrow itself) — so the round-2 pass can train on
[ADVANCED_FEATURE_SET ; GROUP_COLS ; delta_irt]. `max_pid` is the library
precursor count (sizes the streaming accumulators; see `_write_group_features!`).
"""
write_two_round_feature_columns!(file_paths::Vector{String}, max_pid::Int) =
    _write_group_features!(file_paths, max_pid)

# ---------- Symmetric shadow-decoy pairing ----------
#
# For every sampled real row the trainer builds one shadow of the opposite class:
#   real target T  → shadow_decoy  = nearest-iRT decoy's columns except the grafted cols = T's
#   real decoy D   → shadow_target = nearest-iRT target's columns except the grafted cols = D's
#
# The training sample doubles. Shadows inherit their SOURCE's cv_fold (files are fold-split, so
# source and original always share one), keeping precursor-keyed CV folds valid.
#
# The shadows destroy the marginal signal on each grafted feature (at every value taken from a
# target there is 1 target-labeled and 1 decoy-labeled row). Round-2 LGBM is therefore forced to
# use those features jointly with the base features rather than trusting them standalone.
# `_nearest_sorted_idx` and `_mbr_graft_cols` below are the pieces the sampler consumes.

# Vector of indices into `decoy_irt` (sorted asc) whose value is nearest to `q`.
# Binary search on the sorted array; ties break to the lower index. `perm` maps
# sorted-position → original row index.
@inline function _nearest_sorted_idx(sorted_vals::AbstractVector{Float32}, q::Float32)
    n = length(sorted_vals)
    n == 0 && return 0
    # searchsortedfirst returns the first index whose value >= q, or n+1 if none.
    hi = searchsortedfirst(sorted_vals, q)
    if hi > n
        return n
    elseif hi == 1
        return 1
    end
    lo = hi - 1
    return (q - sorted_vals[lo]) <= (sorted_vals[hi] - q) ? lo : hi
end

"""
    log_shadow_fdr_diagnostics(scores, is_target, is_shadow; q_threshold = 0.01)

Log (A) real-only FDR count and (B) shadow-inclusive FDR count of REAL targets passing at
the threshold. Inputs are the round-2 training sample's OOF scores and its target /
shadow flags — shadows never reach disk, so this is computed in memory.
"""
function log_shadow_fdr_diagnostics(
    scores::AbstractVector{<:Real},
    is_target::AbstractVector{Bool},
    is_shadow::AbstractVector{Bool};
    q_threshold::Float64 = 0.01,
)
    isempty(scores) && (@user_info "shadow FDR diagnostics: no rows"; return nothing)

    perm = sortperm(scores; rev=true)  # rank 1 = highest score
    ntot = length(scores)
    # A: real-only FDR — sweep high→low, at each real row compute (real_D / real_T); the
    # q-value at that row is the min FDR from here up, but for "max targets passing ≤ q"
    # the answer is the largest real_T seen at any point where FDR <= q_threshold.
    countA_real_target_pass = let rT = 0, rD = 0, best = 0
        for r in 1:ntot
            i = perm[r]
            is_shadow[i] && continue
            if is_target[i]; rT += 1; else; rD += 1; end
            fdr = rT > 0 ? rD / rT : 1.0
            if fdr <= q_threshold; best = max(best, rT); end
        end
        best
    end
    # B: shadow-inclusive FDR (all rows count). Count REAL targets above the score where FDR = q_threshold.
    countB_real_target_pass = let allT = 0, allD = 0, real_T = 0, best = 0
        for r in 1:ntot
            i = perm[r]
            if is_target[i]; allT += 1; else; allD += 1; end
            !is_shadow[i] && is_target[i] && (real_T += 1)
            fdr = allT > 0 ? allD / allT : 1.0
            if fdr <= q_threshold; best = max(best, real_T); end
        end
        best
    end
    n_real_targets = count(i -> is_target[i] && !is_shadow[i], 1:ntot)
    n_real_decoys  = count(i -> !is_target[i] && !is_shadow[i], 1:ntot)
    n_shadow_targets = count(i -> is_target[i] && is_shadow[i], 1:ntot)
    n_shadow_decoys  = count(i -> !is_target[i] && is_shadow[i], 1:ntot)
    @user_info "shadow FDR diagnostics @ q ≤ $(q_threshold):"
    @user_info "  real:   $(n_real_targets) targets, $(n_real_decoys) decoys"
    @user_info "  shadow: $(n_shadow_targets) shadow-targets, $(n_shadow_decoys) shadow-decoys"
    @user_info "  A (real-only FDR):        real targets passing = $(countA_real_target_pass)"
    @user_info "  B (shadow-included FDR):  real targets passing = $(countB_real_target_pass)"
    return nothing
end
