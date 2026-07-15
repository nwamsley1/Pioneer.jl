# Two-round experiment-wide scoring — CLUSTER-CONSENSUS variant.
# Design + CV protocol: see TWO_ROUND_SCORING.md; clustering: CLUSTER_CONSENSUS_CLUSTERING.md.
#
# Like the twin variant, but the cross-run feature is `cluster_best` = the BEST (max)
# round-1 score of the precursor among the OTHER runs in its CLUSTER (leave-one-out),
# instead of the score in the single most-similar run.
#
# Clustering (fully unsupervised): per-run stringent presence (q<0.001 on the round-1
# score) -> cosine Gram of the binary presence vectors -> eigen-embedding -> KMeans +
# silhouette, with a K=1 floor (best silhouette < 0.2 -> pool to one cluster). Singleton
# cluster -> the run is its own reference (cluster_best = self). Leakage-safe via the
# precursor_idx-keyed CV folds (all instances of a precursor share a fold, so the LOO
# lookups are OOF-consistent). Assumes one row per precursor_idx per file.

using Statistics, LinearAlgebra, Random

# Runtime ENV (a precompile-time const would ignore the env var on a cached load).
two_round_enabled() = get(ENV, "TWO_ROUND", "") == "1"
const TWO_ROUND_FEATURES = [:cluster_best, :delta_irt]
const CLUSTER_PRESENCE_QVAL = 0.001f0     # "present" = round-1 score passes this per-run FDR
const CLUSTER_SILHOUETTE_FLOOR = 0.2      # best silhouette below this -> K=1 (pool)

function log_pass_importance(pass, label::String; topn::Int = 30)
    pass.last_classifier === nothing && return nothing
    model = LightGBMModel(pass.last_classifier, pass.available_features, nothing)
    imp = importance(model)
    imp === nothing && return nothing
    sorted_imp = sort(imp, by = x -> -x[2]); total = sum(x -> x[2], sorted_imp)
    lines = ["$label LGBM feature gains (top $(min(topn,length(sorted_imp))) of $(length(sorted_imp))):"]
    for (fname, gain) in first(sorted_imp, topn)
        pct = total > 0 ? round(100*gain/total, digits=1) : 0.0
        push!(lines, "    $(rpad(string(fname), 42)) $(rpad(round(Int, gain), 10)) ($(pct)%)")
    end
    @user_info join(lines, "\n"); return nothing
end

# Per-run target/decoy q-values from scores (higher score = better).
function _td_qvalues(s1::Vector{Float32}, tgt::Vector{Bool})
    o = sortperm(s1; rev = true); n = length(s1)
    q = Vector{Float32}(undef, n); ct = 0; cd = 0; raw = Vector{Float32}(undef, n)
    @inbounds for i in 1:n
        j = o[i]; tgt[j] ? (ct += 1) : (cd += 1)
        raw[i] = ct > 0 ? Float32(cd) / Float32(ct) : 1.0f0
    end
    m = Inf32
    @inbounds for i in n:-1:1; m = min(m, raw[i]); q[o[i]] = m; end
    return q
end

# Lloyd's k-means on the embedding (few restarts); returns labels.
function _kmeans(X::Matrix{Float64}, K::Int; restarts = 5, iters = 50)
    n, d = size(X); best_lab = ones(Int, n); best_inertia = Inf
    rng = MersenneTwister(1234)
    for _ in 1:restarts
        cent = X[randperm(rng, n)[1:K], :]; lab = ones(Int, n)
        for _ in 1:iters
            changed = false
            @inbounds for i in 1:n
                bd = Inf; bl = 1
                for c in 1:K
                    dd = 0.0; for j in 1:d; e = X[i,j] - cent[c,j]; dd += e*e; end
                    dd < bd && (bd = dd; bl = c)
                end
                lab[i] != bl && (lab[i] = bl; changed = true)
            end
            fill!(cent, 0.0); cnt = zeros(Int, K)
            @inbounds for i in 1:n; cnt[lab[i]] += 1; for j in 1:d; cent[lab[i],j] += X[i,j]; end; end
            for c in 1:K, j in 1:d; cnt[c] > 0 && (cent[c,j] /= cnt[c]); end
            !changed && break
        end
        inertia = 0.0
        @inbounds for i in 1:n, j in 1:d; e = X[i,j] - cent[lab[i],j]; inertia += e*e; end
        inertia < best_inertia && (best_inertia = inertia; best_lab = copy(lab))
    end
    return best_lab
end

# Mean silhouette (euclidean) over points.
function _silhouette(X::Matrix{Float64}, lab::Vector{Int})
    n, d = size(X); length(unique(lab)) < 2 && return -1.0
    D = zeros(n, n)
    @inbounds for i in 1:n, j in i+1:n
        s = 0.0; for k in 1:d; e = X[i,k] - X[j,k]; s += e*e; end
        s = sqrt(s); D[i,j] = s; D[j,i] = s
    end
    tot = 0.0; valid = 0
    for i in 1:n
        own = lab[i]; a_sum = 0.0; a_n = 0; others = Dict{Int,Tuple{Float64,Int}}()
        for j in 1:n
            i == j && continue
            if lab[j] == own; a_sum += D[i,j]; a_n += 1
            else; t = get(others, lab[j], (0.0,0)); others[lab[j]] = (t[1]+D[i,j], t[2]+1); end
        end
        # singleton-cluster point: silhouette defined as 0 (matches sklearn) — counts
        # toward the mean, so degenerate high-K (many singletons) is penalized, not skipped.
        if a_n == 0
            valid += 1; continue
        end
        a = a_sum / a_n; b = Inf
        for (_, t) in others; t[2] > 0 && (b = min(b, t[1]/t[2])); end
        b == Inf && (valid += 1; continue)
        tot += (b - a) / max(a, b); valid += 1
    end
    return valid > 0 ? tot / valid : -1.0
end

# Cluster the runs from per-run presence sets. Returns cluster_id (1-based) per run.
function _cluster_runs(presence::Vector{Vector{UInt32}})
    nr = length(presence)
    nr <= 2 && return ones(Int, nr)
    counts = Float64[length(p) for p in presence]; norms = sqrt.(counts)
    pid_runs = Dict{UInt32, Vector{Int}}()
    for r in 1:nr, p in presence[r]; push!(get!(pid_runs, p, Int[]), r); end
    G = zeros(Float64, nr, nr)
    for (_, rs) in pid_runs
        for a in 1:length(rs), b in a:length(rs); G[rs[a], rs[b]] += 1.0; end
    end
    for r in 1:nr, s in r:nr
        den = norms[r]*norms[s]; v = den > 0 ? G[r,s]/den : 0.0; G[r,s] = v; G[s,r] = v
    end
    ev = eigen(Symmetric(G)); ord = sortperm(ev.values; rev = true)
    k = min(nr-1, 20); tix = ord[1:k]
    emb = ev.vectors[:, tix] .* sqrt.(max.(ev.values[tix], 0.0))'
    best_sil = -Inf; best_lab = ones(Int, nr)
    for K in 2:min(nr-1, 40)
        lab = _kmeans(emb, K); sil = _silhouette(emb, lab)
        sil > best_sil && (best_sil = sil; best_lab = lab)
    end
    @user_info "cluster-consensus: best silhouette $(round(best_sil, digits=3)) (floor $CLUSTER_SILHOUETTE_FLOOR)"
    best_sil < CLUSTER_SILHOUETTE_FLOOR && return ones(Int, nr)
    return best_lab
end

"""
    write_two_round_feature_columns!(file_paths)

Compute `cluster_best` (LOO max round-1 score within the run's cluster) and `delta_irt`
from the round-1 OOF sidecars and write them as columns into the per-file fold Arrows.
"""
function write_two_round_feature_columns!(file_paths::Vector{String})
    nf = length(file_paths)
    file_pid = Vector{Vector{UInt32}}(undef, nf)
    file_s1  = Vector{Vector{Float32}}(undef, nf)
    file_irt = Vector{Vector{Float32}}(undef, nf)
    file_tgt = Vector{Vector{Bool}}(undef, nf)
    file_run = Vector{Int}(undef, nf)   # ms_file_idx (biological run) per fold-file
    best_s1 = Dict{UInt32,Float32}(); ref_irt = Dict{UInt32,Float32}()

    for f in 1:nf
        p = file_paths[f]; tbl = Arrow.Table(p); side = Arrow.Table(p * PASS1_SIDECAR_SUFFIX)
        pid = collect(UInt32.(tbl.precursor_idx)); irt = collect(Float32.(tbl.irt_obs))
        s1 = collect(Float32.(side.trace_prob_prepass)); tgt = collect(Bool.(tbl.target))
        length(s1) == length(pid) || error("cluster-consensus: sidecar row mismatch at $p")
        file_pid[f] = pid; file_s1[f] = s1; file_irt[f] = irt; file_tgt[f] = tgt
        # CV fold files are per-(run,fold); group them back to the biological run via ms_file_idx.
        file_run[f] = length(pid) > 0 ? Int(tbl.ms_file_idx[1]) : -f
        @inbounds for i in eachindex(pid)
            if s1[i] > get(best_s1, pid[i], -1.0f0); best_s1[pid[i]] = s1[i]; ref_irt[pid[i]] = irt[i]; end
        end
    end

    # --- aggregate fold-files to RUNS ---
    runs = sort(unique(file_run)); nrun = length(runs); rix = Dict(r => i for (i, r) in enumerate(runs))
    run_ffs = [Int[] for _ in 1:nrun]; for f in 1:nf; push!(run_ffs[rix[file_run[f]]], f); end
    run_presence = Vector{Vector{UInt32}}(undef, nrun)   # union of confident pids over the run's folds
    run_dict = Vector{Dict{UInt32,Float32}}(undef, nrun) # pid -> best s1 over the run's folds
    for u in 1:nrun
        d = Dict{UInt32,Float32}(); pres = UInt32[]
        for f in run_ffs[u]
            pid = file_pid[f]; s1 = file_s1[f]; tgt = file_tgt[f]; q = _td_qvalues(s1, tgt)
            @inbounds for i in eachindex(pid)
                d[pid[i]] = max(get(d, pid[i], -1.0f0), s1[i])
                (tgt[i] && q[i] < CLUSTER_PRESENCE_QVAL) && push!(pres, pid[i])
            end
        end
        run_presence[u] = pres; run_dict[u] = d
    end

    cluster_id = _cluster_runs(run_presence); nclust = maximum(cluster_id)  # per RUN
    cluster_runs = [Int[] for _ in 1:nclust]; for u in 1:nrun; push!(cluster_runs[cluster_id[u]], u); end
    # per cluster, per precursor: top1 (val, run-unit) and top2 val (for LOO max over OTHER runs)
    top1v = [Dict{UInt32,Float32}() for _ in 1:nclust]
    top1r = [Dict{UInt32,Int}()     for _ in 1:nclust]
    top2v = [Dict{UInt32,Float32}() for _ in 1:nclust]
    for c in 1:nclust, u in cluster_runs[c], (p, v) in run_dict[u]
        if v > get(top1v[c], p, -1.0f0)
            haskey(top1v[c], p) && (top2v[c][p] = top1v[c][p])
            top1v[c][p] = v; top1r[c][p] = u
        elseif v > get(top2v[c], p, -1.0f0)
            top2v[c][p] = v
        end
    end

    for f in 1:nf
        pid = file_pid[f]; irt = file_irt[f]; s1 = file_s1[f]; n = length(pid)
        u = rix[file_run[f]]; c = cluster_id[u]; singleton = length(cluster_runs[c]) <= 1
        cb = Vector{Float32}(undef, n); di = Vector{Float32}(undef, n)
        @inbounds for i in 1:n
            p = pid[i]
            if singleton
                cb[i] = s1[i]
            else
                loo = get(top1r[c], p, 0) == u ? get(top2v[c], p, -1.0f0) : get(top1v[c], p, -1.0f0)
                cb[i] = loo < 0 ? 0.0f0 : loo
            end
            di[i] = abs(irt[i] - ref_irt[p])
        end
        main = DataFrame(Tables.columntable(Arrow.Table(file_paths[f])))
        nrow(main) == n || error("cluster-consensus: arrow row count changed at $(file_paths[f])")
        main[!, :cluster_best] = cb
        main[!, :delta_irt] = di
        writeArrow(file_paths[f], main)
    end
    @user_info "cluster-consensus: $nrun runs ($nf fold-files) -> $nclust cluster(s); wrote cluster_best + delta_irt"
    return nothing
end
