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
const KNN_SCOPE = :all_runs   # :nearest_k | :all_runs — with :all_runs, KNN_K is ignored.
const KNN_AGG = :max          # :mean | :median | :quantile | :max
const KNN_Q = 0.75f0          # only used when KNN_AGG === :quantile
# PEP filter on round-1 s1 before the cross-run feature computation. 1.0 = off.
const PEP_FILTER_THRESHOLD = 1.0
const PEP_FILTER_N_BINS = 100
# Symmetric shadow-decoy injection: after mbr_score is written, for every real
# row inject a paired shadow of the opposite class carrying the source's other
# features but the ORIGINAL row's mbr_score. See TWO_ROUND_SCORING.md.
const SHADOW_DECOY_MODE = :symmetric   # :none | :symmetric
const KNN_COL = if KNN_SCOPE === :all_runs
    Symbol("global_$(KNN_AGG)")   # :global_mean, :global_max, ...
elseif KNN_AGG === :quantile
    Symbol("knn$(KNN_K)_q$(Int(round(KNN_Q * 100)))")
else
    Symbol("knn$(KNN_K)_$(KNN_AGG)")
end
two_round_enabled() = get(ENV, "TWO_ROUND", "") == "1"
# Cross-run shadow-spectrum agreement: donor<->acceptor top-8 smoothed-fragment Hellinger.
# For each (precursor, run), the donor is the neighbor run with the highest s1 (the run supplying
# the KNN_COL value). A true transfer (same real peptide) -> low Hellinger (high agreement); a
# false transfer (precursor absent in this run) -> the acceptor is noise -> high Hellinger. Gives
# round-2 the signal KNN_COL lacks, so it can borrow breadth without leaking. Toggle with
# shadow_hel_enabled().
# Runtime toggle (default off) so base shadow+global_max and +shadow_hel can be A/B'd on the
# same compiled build via ENV["SHADOW_HEL"]=1, without a recompile.
shadow_hel_enabled() = get(ENV, "SHADOW_HEL", "") == "1"
const SHADOW_HEL_COL = Symbol("$(KNN_COL)_shadow_hel")

# GATED_HEL mode (ENV["GATED_HEL"]=1): instead of feeding global_max and shadow_hel as separate
# additive features, feed the single GATED score  global_max_gated = global_max * (1 - shadow_hel).
# Spectral disagreement (high Hellinger) can then ONLY withhold a transfer, never add it — a guard,
# not a booster. See shadow_hel_guard.pdf.
gated_hel_enabled() = get(ENV, "GATED_HEL", "") == "1"
const GATED_COL = Symbol("$(KNN_COL)_gated")

# MULTI_FEATURE mode (ENV["MULTI_FEATURE"]=1): compute several cross-run summaries at once and let
# round-2 select. Safe BECAUSE of the shadow-decoy regularization — each feature's target/decoy
# marginal is forced to 1:1, so the LGBM can't over-trust any single one; it uses whichever add
# signal jointly with the base features. All of these are grafted onto the shadows.
#   global_max/global_mean : max/mean s1 over ALL other runs
#   cluster_max/cluster_mean: max/mean s1 over the top-KNN_K cosine-similar (condition) runs
#   shadow_hel             : donor<->acceptor top-8 fragment Hellinger (donor = global_max run)
multi_feature_enabled() = get(ENV, "MULTI_FEATURE", "") == "1"
const MULTI_FEATURE_COLS = [:global_max, :global_mean, :cluster_max, :cluster_mean, :shadow_hel]

# GATE_GLOBAL mode (with MULTI_FEATURE=1, ENV["GATE_GLOBAL"]=1): apply the shadow_hel gate to ONLY
# the leakage-prone GLOBAL channel — feed global_{max,mean}_gated = global_{max,mean}*(1-shadow_hel)
# — while leaving the condition-scoped cluster_{max,mean} RAW (leak-resistant by scope). Combines the
# gate's transfer PRECISION on the leaky channel with the cluster features' breadth. The gate folds
# agreement in, so no separate shadow_hel term. Graft is selective (the gated global cols only).
gate_global_enabled() = get(ENV, "GATE_GLOBAL", "") == "1"
const MULTI_GATED_COLS = [:global_max_gated, :global_mean_gated, :cluster_max, :cluster_mean]

two_round_features() =
    gated_hel_enabled()     ? [GATED_COL, :delta_irt] :
    (multi_feature_enabled() && gate_global_enabled()) ? vcat(MULTI_GATED_COLS, [:delta_irt]) :
    multi_feature_enabled() ? vcat(MULTI_FEATURE_COLS, [:delta_irt]) :
    shadow_hel_enabled()    ? [KNN_COL, :delta_irt, SHADOW_HEL_COL] :
                              [KNN_COL, :delta_irt]

# Cross-run mbr feature columns to graft onto shadow twins (so each stays 1:1-regularized).
# In MULTI_FEATURE mode, GRAFT_SCOPE selects which get regularized:
#   "all"    (default) — graft all cross-run features (uniform regularization).
#   "global" — graft only the LEAKAGE-PRONE global_* (+shadow_hel) features; leave the
#              condition-scoped cluster_* ungrafted (they're leak-resistant by scope, so the
#              model can lean on them fully). Regularize by leakage risk, not uniformly.
_mbr_graft_cols() =
    gated_hel_enabled()     ? (GATED_COL,) :
    (multi_feature_enabled() && gate_global_enabled()) ? (:global_max_gated, :global_mean_gated) :
    multi_feature_enabled() ?
        (get(ENV, "GRAFT_SCOPE", "all") == "global" ? (:global_max, :global_mean, :shadow_hel) :
                                                       Tuple(c for c in MULTI_FEATURE_COLS)) :
    shadow_hel_enabled()    ? (KNN_COL, SHADOW_HEL_COL) :
                              (KNN_COL,)

# Hellinger distance in [0,1] between two 8-fragment intensity vectors (each L1-normalized to a
# probability, then sqrt). 1.0 (max distance) when either vector has no positive intensity.
@inline function _hellinger8(a::NTuple{8,Float32}, b::NTuple{8,Float32})
    sa = 0.0f0; sb = 0.0f0
    @inbounds for k in 1:8; sa += max(a[k], 0.0f0); sb += max(b[k], 0.0f0); end
    (sa <= 0.0f0 || sb <= 0.0f0) && return 1.0f0
    acc = 0.0f0
    @inbounds for k in 1:8
        p = sqrt(max(a[k], 0.0f0) / sa); q = sqrt(max(b[k], 0.0f0) / sb)
        acc += (p - q)^2
    end
    return sqrt(0.5f0 * acc)
end

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

# Neighbor selection per file. Returns (topk, dicts) where topk[f] is a Vector{Int}
# of neighbor file indices (0 for slots with no valid neighbor) and dicts[f] maps
# precursor_idx -> s1 for lookups. Dispatch:
#   KNN_SCOPE === :nearest_k: top-KNN_K by cosine similarity of s1 profiles (per-file
#     profile = best-per-precursor s1 over the precursor vocabulary, absent = 0).
#   KNN_SCOPE === :all_runs:  every non-self file, unranked. KNN_K is ignored.
function _compute_twin_runs(file_pid::Vector{Vector{UInt32}},
                            file_s1::Vector{Vector{Float32}})
    nf = length(file_pid)
    dicts = Vector{Dict{UInt32,Float32}}(undef, nf)
    for f in 1:nf
        d = Dict{UInt32,Float32}(); pid = file_pid[f]; s1 = file_s1[f]
        @inbounds for i in eachindex(pid); d[pid[i]] = s1[i]; end
        dicts[f] = d
    end
    if KNN_SCOPE === :all_runs
        topk = [Int[g for g in 1:nf if g != f] for f in 1:nf]
        return topk, dicts
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
# Per-file top-K cosine-similar runs (condition neighbors) on the s1 profiles. Mirrors the
# Gram-matrix path in _compute_twin_runs; used by MULTI_FEATURE for the cluster_* features.
function _cosine_topk(dicts::Vector{Dict{UInt32,Float32}}, K::Int)
    nf = length(dicts)
    G = zeros(Float64, nf, nf)
    for f in 1:nf, g in f:nf
        small, big = length(dicts[f]) <= length(dicts[g]) ? (dicts[f], dicts[g]) : (dicts[g], dicts[f])
        acc = 0.0
        for (p, v) in small; acc += Float64(v) * Float64(get(big, p, 0.0f0)); end
        G[f, g] = acc; G[g, f] = acc
    end
    out = [Int[] for _ in 1:nf]
    for f in 1:nf
        nrm = sqrt(G[f, f]); nrm <= 0 && continue
        cand = Tuple{Int,Float64}[]
        for g in 1:nf
            (g == f || G[g, g] <= 0) && continue
            push!(cand, (g, G[f, g] / (nrm * sqrt(G[g, g]))))
        end
        sort!(cand; by = x -> -x[2])
        out[f] = [c[1] for c in first(cand, min(K, length(cand)))]
    end
    return out
end

function write_two_round_feature_columns!(file_paths::Vector{String})
    nf = length(file_paths)
    file_pid = Vector{Vector{UInt32}}(undef, nf)
    file_s1  = Vector{Vector{Float32}}(undef, nf)
    file_irt = Vector{Vector{Float32}}(undef, nf)
    file_tgt = Vector{Vector{Bool}}(undef, nf)  # target flag per row, for PEP filter
    file_frag = Vector{Matrix{Float32}}(undef, nf)                    # n x 8 shadow spectrum
    frag_dicts = Vector{Dict{UInt32,NTuple{8,Float32}}}(undef, nf)    # per run: pid -> frag8
    best_s1  = Dict{UInt32,Float32}()   # precursor -> highest s1 seen (for ref_irt)
    ref_irt  = Dict{UInt32,Float32}()   # precursor -> irt_obs at its highest-s1 instance
    _FRAG_COLS = (:frag1_smoothed_intensity, :frag2_smoothed_intensity, :frag3_smoothed_intensity,
                  :frag4_smoothed_intensity, :frag5_smoothed_intensity, :frag6_smoothed_intensity,
                  :frag7_smoothed_intensity, :frag8_smoothed_intensity)

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
        # Shadow spectrum: top-8 smoothed fragment intensities (guarded; zeros if columns absent).
        fr = zeros(Float32, length(pid), 8); fd = Dict{UInt32,NTuple{8,Float32}}()
        if (shadow_hel_enabled() || gated_hel_enabled() || multi_feature_enabled()) && all(c -> hasproperty(tbl, c), _FRAG_COLS)
            @inbounds for (k, c) in enumerate(_FRAG_COLS)
                col = getproperty(tbl, c)
                for i in eachindex(pid); fr[i, k] = Float32(col[i]); end
            end
            @inbounds for i in eachindex(pid); fd[pid[i]] = ntuple(k -> fr[i, k], 8); end
        end
        file_frag[f] = fr; frag_dicts[f] = fd
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

    # MULTI_FEATURE: compute several cross-run summaries at once (all shadow-regularized downstream).
    if multi_feature_enabled()
        clus_topk = _cosine_topk(dicts, KNN_K)   # per-file top-K cosine-similar (condition) runs
        for f in 1:nf
            pid = file_pid[f]; irt = file_irt[f]; n = length(pid); fr = file_frag[f]
            gmax = Vector{Float32}(undef, n); gmean = Vector{Float32}(undef, n)
            cmax = Vector{Float32}(undef, n); cmean = Vector{Float32}(undef, n)
            sh   = Vector{Float32}(undef, n); di    = Vector{Float32}(undef, n)
            all_runs = [g for g in 1:nf if g != f]; clus = clus_topk[f]
            @inbounds for i in 1:n
                p = pid[i]
                mx = 0.0f0; sm = 0.0f0; cnt = 0; dmx = -1.0f0; donor = 0
                for g in all_runs
                    v = get(dicts[g], p, -1.0f0)
                    v < 0.0f0 && continue
                    v > mx && (mx = v); sm += v; cnt += 1
                    v > dmx && (dmx = v; donor = g)
                end
                gmax[i] = mx; gmean[i] = cnt > 0 ? sm / cnt : 0.0f0
                cm = 0.0f0; cs = 0.0f0; ccnt = 0
                for g in clus
                    v = get(dicts[g], p, -1.0f0)
                    v < 0.0f0 && continue
                    v > cm && (cm = v); cs += v; ccnt += 1
                end
                cmax[i] = cm; cmean[i] = ccnt > 0 ? cs / ccnt : 0.0f0
                sh[i] = donor == 0 ? 1.0f0 :
                    _hellinger8(ntuple(k -> fr[i, k], 8), get(frag_dicts[donor], p, ntuple(_ -> 0.0f0, 8)))
                di[i] = abs(irt[i] - ref_irt[p])
            end
            main = DataFrame(Tables.columntable(Arrow.Table(file_paths[f])))
            nrow(main) == n || error("write_two_round_feature_columns!: arrow row count changed at $(file_paths[f])")
            main[!, :global_max] = gmax; main[!, :global_mean] = gmean
            main[!, :cluster_max] = cmax; main[!, :cluster_mean] = cmean
            main[!, :shadow_hel] = sh;   main[!, :delta_irt] = di
            # GATE_GLOBAL: gate ONLY the global channel by donor agreement; cluster_* stay raw.
            if gate_global_enabled()
                main[!, :global_max_gated]  = gmax  .* (1.0f0 .- sh)
                main[!, :global_mean_gated] = gmean .* (1.0f0 .- sh)
            end
            writeArrow(file_paths[f], main)
        end
        @user_info "two-round MULTI_FEATURE: wrote global_{max,mean}$(gate_global_enabled() ? "(+_gated)" : "") + cluster_{max,mean} (top-$KNN_K cosine) + shadow_hel + delta_irt to $nf files"
        return
    end

    # Pass B: per file, compute the two features and rewrite the Arrow with them.
    max_neighbors = KNN_SCOPE === :all_runs ? max(nf - 1, 1) : KNN_K
    scratch = Vector{Float32}(undef, max_neighbors)   # reused per row for the median/quantile path
    for f in 1:nf
        pid = file_pid[f]; irt = file_irt[f]; n = length(pid)
        km = Vector{Float32}(undef, n); di = Vector{Float32}(undef, n)
        sh = (shadow_hel_enabled() || gated_hel_enabled()) ? Vector{Float32}(undef, n) : Float32[]
        # Resolve up to KNN_K neighbor dicts (contiguous from front since sorted-descending;
        # a 0 slot means no more valid neighbors). Track the run index of each for the donor lookup.
        ndicts = Dict{UInt32,Float32}[]; nruns = Int[]
        for g in topk[f]
            g == 0 && break
            push!(ndicts, dicts[g]); push!(nruns, g)
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
        elseif KNN_AGG === :max
            @inbounds for i in 1:n
                m = 0.0f0
                for d in ndicts
                    v = get(d, pid[i], 0.0f0)
                    v > m && (m = v)
                end
                km[i] = m   # 0 if nvalid == 0
                di[i] = abs(irt[i] - ref_irt[pid[i]])
            end
        else
            error("KNN_AGG must be :mean, :median, :quantile, or :max (got $KNN_AGG)")
        end
        # Cross-run shadow-spectrum agreement: donor = neighbor run with the highest s1 for this
        # precursor (the run supplying KNN_COL); Hellinger between this row's shadow spectrum and
        # the donor's. No donor -> 1.0 (max distance / no agreement).
        if shadow_hel_enabled() || gated_hel_enabled()
            fr = file_frag[f]
            @inbounds for i in 1:n
                p = pid[i]; bv = -1.0f0; bo = 0
                for k in 1:nvalid
                    v = get(ndicts[k], p, -1.0f0); v > bv && (bv = v; bo = k)
                end
                sh[i] = bo == 0 ? 1.0f0 :
                    _hellinger8(ntuple(k -> fr[i, k], 8),
                                get(frag_dicts[nruns[bo]], p, ntuple(_ -> 0.0f0, 8)))
            end
        end
        main = DataFrame(Tables.columntable(Arrow.Table(file_paths[f])))
        nrow(main) == n ||
            error("write_two_round_feature_columns!: arrow row count changed at $(file_paths[f])")
        main[!, KNN_COL]    = km
        main[!, :delta_irt] = di
        shadow_hel_enabled() && (main[!, SHADOW_HEL_COL] = sh)
        # GATED_HEL guard: global_max_gated = global_max * (1 - shadow_hel). Spectral disagreement
        # (high Hellinger) can only withhold the transfer, never add it. See shadow_hel_guard.pdf.
        gated_hel_enabled() && (main[!, GATED_COL] = km .* (1.0f0 .- sh))
        writeArrow(file_paths[f], main)
    end
    agg_desc = KNN_AGG === :quantile ? "$(KNN_Q)-quantile" : String(KNN_AGG)
    if KNN_SCOPE === :all_runs
        @user_info "two-round: wrote $(KNN_COL) ($(agg_desc) over all $(nf-1) other runs) + delta_irt to $nf files"
    else
        n_full    = count(nb -> all(>(0), nb), topk)
        n_partial = count(nb -> nb[1] > 0 && !all(>(0), nb), topk)
        @user_info "two-round: wrote $(KNN_COL) ($(agg_desc) of top-$(KNN_K)) + delta_irt to $nf files ($n_full with all $KNN_K neighbors, $n_partial with fewer, cosine s1-profile)"
    end
    return nothing
end

# ---------- Symmetric shadow-decoy injection ----------
#
# For every real row in a fold-file, inject one shadow of the opposite class:
#   real target T  → shadow_decoy  = nearest-iRT decoy's columns except mbr_score = T.mbr_score
#   real decoy D   → shadow_target = nearest-iRT target's columns except mbr_score = D.mbr_score
# `is_shadow` flags the shadow (false on originals).
#
# Row count exactly doubles when both classes are present in a file. Shadows
# inherit their SOURCE's precursor_idx / cv_fold, so precursor-keyed CV folds
# stay valid. LGBM trains transparently on the enlarged set.
#
# The shadows destroy the marginal-on-mbr_score signal (at every mbr_score
# value from a target, there is 1 target-labeled and 1 decoy-labeled row).
# Round-2 LGBM is therefore forced to use mbr_score jointly with the base
# features rather than trusting it as a standalone signal.

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

# Sidecar suffix + column: mirrors the round-1 sidecar so we can also
# extend/truncate the sidecar rows to stay aligned with the fold-file.
const SHADOW_SIDECAR_SCORE_COL = :trace_prob_prepass

"""
    inject_shadow_decoys!(file_paths::Vector{String}) -> Int

For each fold-file, symmetrically inject shadow rows (see module comment). Also
extends `.pass1_sidecar.arrow` so its row count stays aligned; shadow sidecar
scores are copied from the SOURCE row (they'll be overwritten when round-2
predicts). Returns total shadows injected across all files.
"""
function inject_shadow_decoys!(file_paths::Vector{String})
    SHADOW_DECOY_MODE === :none && return 0
    total = 0
    n_skipped_target_side = 0
    n_skipped_decoy_side  = 0
    for p in file_paths
        main = DataFrame(Tables.columntable(Arrow.Table(p)))
        # If the sidecar exists (round-1 already wrote it), we mirror row count.
        side_path = p * PASS1_SIDECAR_SUFFIX
        side_present = isfile(side_path)
        side = side_present ? DataFrame(Tables.columntable(Arrow.Table(side_path))) : nothing
        side_present && nrow(side) != nrow(main) &&
            error("inject_shadow_decoys!: sidecar row count mismatch at $p ($(nrow(side)) vs $(nrow(main)))")

        n = nrow(main)
        n == 0 && continue
        tgt = main.target::AbstractVector{Bool}
        irt = Float32.(main.irt_obs)
        target_idx = findall(tgt)
        decoy_idx  = findall(.!tgt)

        if isempty(target_idx) || isempty(decoy_idx)
            # No opposite class present — can't build symmetric shadows.
            n_skipped_target_side += length(target_idx)
            n_skipped_decoy_side  += length(decoy_idx)
            continue
        end

        # For each real target, source_row = nearest-iRT decoy → shadow_decoy(source=decoy).
        # For each real decoy,  source_row = nearest-iRT target → shadow_target(source=target).
        decoy_perm  = sortperm(view(irt, decoy_idx))
        target_perm = sortperm(view(irt, target_idx))
        decoy_irt_sorted  = irt[decoy_idx[decoy_perm]]
        target_irt_sorted = irt[target_idx[target_perm]]

        shadow_src_row = Vector{Int}(undef, n)   # source row index in main
        shadow_mbr_row = Vector{Int}(undef, n)   # row whose mbr_score to graft
        @inbounds for r in 1:n
            if tgt[r]
                pos = _nearest_sorted_idx(decoy_irt_sorted, irt[r])
                pos == 0 && (n_skipped_target_side += 1; shadow_src_row[r] = 0; continue)
                shadow_src_row[r] = decoy_idx[decoy_perm[pos]]
                shadow_mbr_row[r] = r
            else
                pos = _nearest_sorted_idx(target_irt_sorted, irt[r])
                pos == 0 && (n_skipped_decoy_side += 1; shadow_src_row[r] = 0; continue)
                shadow_src_row[r] = target_idx[target_perm[pos]]
                shadow_mbr_row[r] = r
            end
        end
        valid_rows = findall(!=(0), shadow_src_row)
        isempty(valid_rows) && continue

        # Build the shadow subset in bulk by selecting rows [shadow_src_row[i] for i in valid_rows]
        # from main, then overwrite mbr_score.
        src_positions = shadow_src_row[valid_rows]
        mbr_positions = shadow_mbr_row[valid_rows]
        shadows = main[src_positions, :]

        # Overwrite the cross-run transfer feature column(s) with the ORIGINAL row's value(s),
        # so each feature's target/decoy marginal is 1:1 (forces joint use with base features).
        # Both KNN_COL (the score) and SHADOW_HEL_COL (the donor<->acceptor agreement) are grafted.
        for mbr_col in _mbr_graft_cols()
            hasproperty(main, mbr_col) && (shadows[!, mbr_col] = main[mbr_positions, mbr_col])
        end

        # Set is_shadow flag (add column to `main` too if absent — schema union).
        if !hasproperty(main, :is_shadow)
            main[!, :is_shadow] = falses(nrow(main))
        end
        shadows[!, :is_shadow] = trues(nrow(shadows))

        # Concat.
        main_out = vcat(main, shadows)
        writeArrow(p, main_out)

        # Mirror sidecar: append sidecar rows from src_positions.
        if side_present
            side_shadow = side[src_positions, :]
            side_out = vcat(side, side_shadow)
            writeArrow(side_path, side_out)
        end

        total += nrow(shadows)
    end
    @user_info "shadow-decoy injection: injected $total shadows across $(length(file_paths)) files (skipped $(n_skipped_target_side) target-side + $(n_skipped_decoy_side) decoy-side for missing opposite class)"
    return total
end

"""
    remove_shadow_decoys!(file_paths::Vector{String}) -> Int

Filter out `is_shadow == true` rows from each fold-file and its sidecar. Called
after round-2 scoring completes and after any shadow-inclusive diagnostics are
computed. Returns total shadows removed.
"""
function remove_shadow_decoys!(file_paths::Vector{String})
    total = 0
    for p in file_paths
        main = DataFrame(Tables.columntable(Arrow.Table(p)))
        !hasproperty(main, :is_shadow) && continue
        side_path = p * PASS1_SIDECAR_SUFFIX
        side_present = isfile(side_path)
        side = side_present ? DataFrame(Tables.columntable(Arrow.Table(side_path))) : nothing
        side_present && nrow(side) != nrow(main) &&
            error("remove_shadow_decoys!: sidecar row count mismatch at $p")

        real_mask = .!Bool.(main.is_shadow)
        n_removed = sum(.!real_mask)
        n_removed == 0 && continue
        main_kept = main[real_mask, :]
        # Drop the flag column too so the schema matches the pre-injection layout.
        select!(main_kept, Not(:is_shadow))
        writeArrow(p, main_kept)
        if side_present
            writeArrow(side_path, side[real_mask, :])
        end
        total += n_removed
    end
    @user_info "shadow-decoy removal: removed $total shadow rows from $(length(file_paths)) files"
    return total
end

"""
    log_shadow_fdr_diagnostics(file_paths::Vector{String}; q_threshold::Float64 = 0.01)

BEFORE removing shadows: log (A) real-only FDR count and (B) shadow-inclusive
FDR count of REAL targets passing at the threshold. Assumes the round-2 sidecar
(`trace_prob_prepass`) is the current OOF score.
"""
function log_shadow_fdr_diagnostics(file_paths::Vector{String}; q_threshold::Float64 = 0.01)
    scores = Float32[]; is_target = Bool[]; is_shadow = Bool[]
    for p in file_paths
        main = Arrow.Table(p)
        side = Arrow.Table(p * PASS1_SIDECAR_SUFFIX)
        n = length(main.precursor_idx)
        length(side.trace_prob_prepass) == n ||
            error("log_shadow_fdr_diagnostics: sidecar/main row count mismatch at $p")
        append!(scores, Float32.(side.trace_prob_prepass))
        append!(is_target, Bool.(main.target))
        shadow_col = hasproperty(main, :is_shadow) ? Bool.(main.is_shadow) : falses(n)
        append!(is_shadow, shadow_col)
    end
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
