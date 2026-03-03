#=
FRAGCORR Top-N Fragment Comparison
===================================
Compares the effect of using top 5 vs 6 vs 7 fragments for FRAGCORR scoring.
Uses a single search run (with 7 intensity columns) and subsets to 5, 6, or 7
columns for each analysis variant.

For each fragment count:
  1. Compute unsupervised FRAGCORR scores (eigengap, lambda1_frac, etc.)
  2. Join with PSM data → train LightGBM baseline vs augmented
  3. Compute global probs → simulate filter → compare precursor recovery
  4. DIA-NN validation

Usage:
    # After running Pioneer with the 7-fragment code:
    #   SearchDIA("/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/OlsenEclipse_fragcorr_7frag/config.json")
    # Then:
    include("scripts/fragcorr_topn_comparison.jl")
=#

using Arrow, DataFrames, Tables
using Statistics, LinearAlgebra
using StatsBase: mad
using Plots
using LightGBM
using CSV

# ============================================================
# Configuration
# ============================================================

const TN_DATA_DIR = "/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/OlsenEclipse_fragcorr_7frag/temp_data"
const TN_PSM_DIR  = joinpath(TN_DATA_DIR, "first_pass_psms")
const TN_CORR_DIR = joinpath(TN_DATA_DIR, "first_pass_corr")

const TN_LIBRARY_PATH = "/Users/nathanwamsley/Data/SPEC_LIBS/3P_rank_based.poin"
const TN_DIANN_PR_MATRIX = joinpath("/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse",
    "DIA-NN_results", "OlsenExplorisThreeProteome500ng-11-24-2025-report.pr_matrix.tsv")
const TN_PIONEER_FILES = ["E45H50Y5_2", "E45H50Y5_3", "E45H50Y5_4",
                           "E5H50Y45_1", "E5H50Y45_2", "E5H50Y45_4"]

# Fragment count variants to compare
const FRAG_COUNTS = [5, 6, 7]

# ============================================================
# Scoring functions (adapted from fragcorr_unsupervised_scores.jl)
# ============================================================

"""Extract n_scans × n_frags intensity matrix using the specified columns."""
function tn_get_intensity_matrix(sub::SubDataFrame, intensity_cols::Vector{Symbol})
    n = nrow(sub)
    p = length(intensity_cols)
    X = Matrix{Float64}(undef, n, p)
    @inbounds for (j, col) in enumerate(intensity_cols)
        copyto!(view(X, :, j), sub[!, col])
    end
    return X
end

"""Clamp negatives to zero."""
tn_clamp_nn!(X) = (X .= max.(X, 0.0); X)

"""Consensus apex index: argmax of mean across max-normalized columns."""
function tn_consensus_apex_idx(X::Matrix{Float64})
    n, p = size(X)
    best_i, best_v = 1, -Inf
    @inbounds for i in 1:n
        v = 0.0
        for j in 1:p
            mx = 0.0
            for ii in 1:n; mx = max(mx, X[ii, j]); end
            v += mx > 0.0 ? X[i, j] / mx : 0.0
        end
        if v > best_v; best_v = v; best_i = i; end
    end
    return best_i
end

"""Pairwise correlation eigenvalue scores."""
function tn_score_corr_eigen(X::Matrix{Float64}; eps=1e-12)
    n, p = size(X)
    C = Matrix{Float64}(undef, p, p)
    stds = Vector{Float64}(undef, p)
    means = Vector{Float64}(undef, p)
    @inbounds for j in 1:p
        s = 0.0; ss = 0.0
        for i in 1:n; s += X[i, j]; end
        means[j] = s / n
        for i in 1:n; ss += (X[i, j] - means[j])^2; end
        stds[j] = sqrt(ss / max(n - 1, 1))
    end
    @inbounds for j1 in 1:p
        C[j1, j1] = 1.0
        for j2 in (j1+1):p
            if stds[j1] > eps && stds[j2] > eps
                cov = 0.0
                for i in 1:n
                    cov += (X[i, j1] - means[j1]) * (X[i, j2] - means[j2])
                end
                c = cov / (max(n - 1, 1) * stds[j1] * stds[j2])
                C[j1, j2] = c; C[j2, j1] = c
            else
                C[j1, j2] = 0.0; C[j2, j1] = 0.0
            end
        end
    end
    evals = eigvals(Symmetric(C))
    sort!(evals, rev=true)
    tr = sum(evals)
    lambda1_frac = tr > eps ? evals[1] / tr : 0.0
    eigengap = p >= 2 ? evals[1] - evals[2] : 0.0
    corrs = Vector{Float64}(undef, p * (p - 1) ÷ 2)
    k = 0
    @inbounds for j1 in 1:p, j2 in (j1+1):p
        k += 1; corrs[k] = C[j1, j2]
    end
    return (lambda1_frac=lambda1_frac, eigengap=eigengap, median_corr=median(corrs))
end

"""Rank-1 explained variance (SVD)."""
function tn_score_ev1(X::Matrix{Float64}; eps=1e-12)
    all_zero = true
    @inbounds for v in X; v != 0.0 && (all_zero = false; break); end
    all_zero && return 0.0
    S = svdvals(X)
    total = sum(abs2, S)
    return total > eps ? S[1]^2 / total : 0.0
end

"""3-point weighted smoothing."""
function tn_smooth3!(dst::AbstractVector, src::AbstractVector)
    n = length(src)
    n <= 1 && (copyto!(dst, src); return dst)
    first_val = (2.0/3.0) * src[1] + (1.0/3.0) * src[2]
    last_val  = (2.0/3.0) * src[n] + (1.0/3.0) * src[n-1]
    @inbounds for i in 2:n-1
        dst[i] = 0.5 * src[i] + 0.25 * (src[i-1] + src[i+1])
    end
    dst[1] = first_val; dst[n] = last_val
    return dst
end

"""Find best fragment by pairwise correlation sum."""
function tn_best_fragment_idx(X::Matrix{Float64}; eps=1e-12)
    n, p = size(X)
    means = Vector{Float64}(undef, p)
    stds  = Vector{Float64}(undef, p)
    @inbounds for j in 1:p
        s = 0.0
        for i in 1:n; s += X[i, j]; end
        means[j] = s / n
        ss = 0.0
        for i in 1:n; ss += (X[i, j] - means[j])^2; end
        stds[j] = sqrt(ss / max(n - 1, 1))
    end
    C = zeros(p, p)
    @inbounds for j1 in 1:p
        C[j1, j1] = 1.0
        for j2 in (j1+1):p
            if stds[j1] > eps && stds[j2] > eps
                cov_val = 0.0
                for i in 1:n
                    cov_val += (X[i, j1] - means[j1]) * (X[i, j2] - means[j2])
                end
                c = cov_val / (max(n - 1, 1) * stds[j1] * stds[j2])
                C[j1, j2] = c; C[j2, j1] = c
            end
        end
    end
    best_j = 1; best_sum = -Inf
    @inbounds for j in 1:p
        s = 0.0
        for j2 in 1:p; j2 == j && continue; s += C[j, j2]; end
        if s > best_sum; best_sum = s; best_j = j; end
    end
    return best_j
end

"""Best-fragment correlation score."""
function tn_score_best_frag_corr(X::Matrix{Float64}; eps=1e-12)
    n, p = size(X)
    n < 3 && return (corr_sum=0.0, corr_min=0.0, corr_mean=0.0)
    best_j = tn_best_fragment_idx(X)
    elution = Vector{Float64}(undef, n)
    tn_smooth3!(elution, @view X[:, best_j])
    el_std = std(elution)
    corrs = Vector{Float64}(undef, p)
    @inbounds for j in 1:p
        col = @view X[:, j]
        s = std(col)
        corrs[j] = (s > eps && el_std > eps) ? cor(col, elution) : 0.0
    end
    return (corr_sum=sum(corrs), corr_min=minimum(corrs), corr_mean=mean(corrs))
end

"""Compute FRAGCORR scores for one precursor group using specified intensity columns."""
function tn_compute_scores(sub::SubDataFrame, intensity_cols::Vector{Symbol}; min_scans=3)
    X = tn_get_intensity_matrix(sub, intensity_cols)
    n = size(X, 1)
    n < min_scans && return nothing
    Xnn = tn_clamp_nn!(copy(X))
    ev1 = tn_score_ev1(Xnn)
    corr = tn_score_corr_eigen(Xnn)
    bfc  = tn_score_best_frag_corr(Xnn)
    return (
        precursor_idx    = first(sub.precursor_idx),
        target           = first(sub.target),
        n_scans          = n,
        ev1_raw          = ev1,
        lambda1_frac_raw = corr.lambda1_frac,
        eigengap_raw     = corr.eigengap,
        median_corr_raw  = corr.median_corr,
        best_frag_corr_mean = bfc.corr_mean,
    )
end

"""Score all precursor groups for a given fragment count."""
function tn_score_all(gdf, intensity_cols::Vector{Symbol}; min_scans=3)
    results = Vector{NamedTuple}()
    sizehint!(results, length(gdf))
    for sub in gdf
        r = tn_compute_scores(sub, intensity_cols; min_scans)
        r !== nothing && push!(results, r)
    end
    return DataFrame(results)
end

# ============================================================
# Q-value utilities (from fragcorr_unsupervised_scores.jl)
# ============================================================

function tn_compute_qvalues(scores::AbstractVector{<:Real},
                            is_target::AbstractVector{Bool};
                            higher_is_better::Bool=true)
    n = length(scores)
    order = sortperm(scores, rev=higher_is_better)
    cum_targets = 0; cum_decoys = 0
    fdr = Vector{Float64}(undef, n)
    for idx in order
        if is_target[idx]; cum_targets += 1; else; cum_decoys += 1; end
        fdr[idx] = cum_targets > 0 ? cum_decoys / cum_targets : 1.0
    end
    qvals = Vector{Float64}(undef, n)
    min_fdr = 1.0
    for idx in reverse(order)
        min_fdr = min(min_fdr, fdr[idx])
        qvals[idx] = min_fdr
    end
    return qvals
end

function tn_targets_at_qvalue(qvals, is_target, threshold)
    c = 0
    @inbounds for i in eachindex(qvals)
        c += (is_target[i] && qvals[i] <= threshold)
    end
    return c
end

# ============================================================
# Log-odds combine (from getBestPrecursorsAccrossRuns.jl)
# ============================================================

function tn_logodds_combine(probs::Vector{Float64}, top_n::Int)::Float64
    isempty(probs) && return 0.0
    n = min(length(probs), top_n)
    sorted = sort(probs; rev=true)
    selected = @view sorted[1:n]
    eps = 1e-6
    lo = log.(clamp.(selected, 0.1, 1 - eps) ./ (1 .- clamp.(selected, 0.1, 1 - eps)))
    return 1.0 / (1 + exp(-sum(lo) / n))
end

# ============================================================
# CV LightGBM (from fragcorr_model_comparison.jl)
# ============================================================

tn_cv_fold(precursor_idx::AbstractVector) = Int.(precursor_idx .% 3 .+ 1)

function tn_run_cv_lightgbm(df::DataFrame, features::Vector{Symbol}; n_folds::Int=3)
    folds = tn_cv_fold(df.precursor_idx)
    preds = fill(NaN, nrow(df))
    for k in 1:n_folds
        train_mask = folds .!= k
        test_mask  = folds .== k
        X_train = Matrix{Float32}(df[train_mask, features])
        y_train = Int.(df[train_mask, :target])
        X_test  = Matrix{Float32}(df[test_mask, features])
        estimator = LGBMClassification(
            objective = "binary", metric = ["binary_logloss"],
            num_iterations = 200, learning_rate = 0.1,
            max_depth = 4, num_leaves = 15, num_class = 1,
            feature_fraction = 0.8, min_data_in_leaf = 200,
            min_gain_to_split = 0.5, num_threads = Threads.nthreads(),
            verbosity = -1, seed = 1776, deterministic = true, force_row_wise = true,
        )
        fit!(estimator, X_train, y_train; verbosity=-1)
        raw = LightGBM.predict(estimator, X_test)
        ŷ = ndims(raw) == 2 ? dropdims(raw; dims=2) : raw
        preds[test_mask] .= Float64.(ŷ)
    end
    return preds
end

# ============================================================
# DIA-NN key matching (from fragcorr_diann_recovery.jl)
# ============================================================

const TN_ModTuple = Tuple{Int, Int}
const TN_CanonicalKey = Tuple{String, Vector{TN_ModTuple}, Int}

function tn_parse_diann_modified_seq(mod_seq::AbstractString, charge::Int)::TN_CanonicalKey
    stripped = IOBuffer(); mods = TN_ModTuple[]; pos = 0; i = 1
    while i <= ncodeunits(mod_seq)
        c = mod_seq[i]
        if c == '('
            close_idx = findnext(')', mod_seq, i)
            mod_str = mod_seq[i+1:close_idx-1]
            colon_idx = findlast(':', mod_str)
            push!(mods, (pos, parse(Int, mod_str[colon_idx+1:end])))
            i = close_idx + 1
        else
            pos += 1; write(stripped, c); i += 1
        end
    end
    sort!(mods)
    return (String(take!(stripped)), mods, charge)
end

function tn_parse_pioneer_structural_mods(mods_str::AbstractString)::Vector{TN_ModTuple}
    isempty(mods_str) && return TN_ModTuple[]
    mods = TN_ModTuple[]
    for m in eachmatch(r"\((\d+),\w,(\w+:\d+)\)", mods_str)
        pos = parse(Int, m.captures[1])
        mod_id_str = m.captures[2]
        colon_idx = findlast(':', mod_id_str)
        push!(mods, (pos, parse(Int, mod_id_str[colon_idx+1:end])))
    end
    sort!(mods)
    return mods
end

function tn_pioneer_canonical_key(seq, structural_mods, charge)::TN_CanonicalKey
    mods = tn_parse_pioneer_structural_mods(ismissing(structural_mods) ? "" : structural_mods)
    return (String(seq), mods, Int(charge))
end

# ============================================================
# Global prob computation (from fragcorr_global_prob_sim.jl)
# ============================================================

function tn_compute_global_probs(df::DataFrame, scores::Vector{Float64}, n_files::Int)
    top_n = max(1, floor(Int, sqrt(n_files)))
    prec_file_best = Dict{UInt32, Dict{UInt32, Float64}}()
    prec_target = Dict{UInt32, Bool}()
    prec_diann  = Dict{UInt32, Bool}()
    for i in 1:nrow(df)
        pidx = UInt32(df.precursor_idx[i])
        fidx = UInt32(df.ms_file_idx[i])
        s    = scores[i]
        if !haskey(prec_file_best, pidx)
            prec_file_best[pidx] = Dict{UInt32, Float64}(fidx => s)
            prec_target[pidx] = df.target[i]
            prec_diann[pidx]  = hasproperty(df, :is_diann) ? df.is_diann[i] : false
        else
            fb = prec_file_best[pidx]
            fb[fidx] = max(get(fb, fidx, -Inf), s)
            if hasproperty(df, :is_diann) && df.is_diann[i]
                prec_diann[pidx] = true
            end
        end
    end
    n_prec = length(prec_file_best)
    out = DataFrame(
        precursor_idx = Vector{UInt32}(undef, n_prec),
        target        = Vector{Bool}(undef, n_prec),
        is_diann      = Vector{Bool}(undef, n_prec),
        global_prob   = Vector{Float64}(undef, n_prec),
    )
    for (i, (pidx, file_dict)) in enumerate(prec_file_best)
        out.precursor_idx[i] = pidx
        out.target[i]        = prec_target[pidx]
        out.is_diann[i]      = get(prec_diann, pidx, false)
        out.global_prob[i]   = tn_logodds_combine(collect(values(file_dict)), top_n)
    end
    return out
end

# ============================================================
# Main: Load data, run comparison
# ============================================================

println("\n╔══════════════════════════════════════════════════════════╗")
println("║  FRAGCORR Top-N Fragment Comparison (5 vs 6 vs 7)      ║")
println("╚══════════════════════════════════════════════════════════╝")

# ── Load PSM + correlation data ──
println("\nLoading data...")
psm_files = filter(f -> endswith(f, ".arrow"), readdir(TN_PSM_DIR))
corr_files_set = Set(readdir(TN_CORR_DIR))
common_files = filter(f -> f in corr_files_set, psm_files)
println("  Found $(length(common_files)) file pairs")

all_psm_dfs  = DataFrame[]
all_corr_dfs = Dict{Int, DataFrame}()  # n_frags → concatenated corr df

for f in common_files
    println("\n── $f ──")
    psm_df = DataFrame(Tables.columntable(Arrow.Table(joinpath(TN_PSM_DIR, f))))
    push!(all_psm_dfs, psm_df)

    # Load correlation data (has intensity_1 through intensity_7)
    corr_raw = DataFrame(Tables.columntable(Arrow.Table(joinpath(TN_CORR_DIR, f))))
    sort!(corr_raw, [:precursor_idx, :scan_idx])

    # Score with each fragment count
    for nf in FRAG_COUNTS
        cols = [Symbol("intensity_$i") for i in 1:nf]
        gdf = groupby(corr_raw, :precursor_idx)
        scored = tn_score_all(gdf, cols)
        if !haskey(all_corr_dfs, nf)
            all_corr_dfs[nf] = scored
        else
            append!(all_corr_dfs[nf], scored)
        end
    end
    println("  PSM rows: $(nrow(psm_df))")
end

psm_combined = vcat(all_psm_dfs...; cols=:union)
println("\n════════════════════════════════════════")
println("Total PSMs: $(nrow(psm_combined))")
println("════════════════════════════════════════")

# ── DIA-NN annotation ──
println("\nBuilding DIA-NN annotation...")
lib_df = DataFrame(Arrow.Table(joinpath(TN_LIBRARY_PATH, "precursors_table.arrow")))

mc_pidx_set = Set{UInt32}(psm_combined[psm_combined.target, :precursor_idx])
pidx_to_key = Dict{UInt32, TN_CanonicalKey}()
for i in 1:nrow(lib_df)
    pidx = UInt32(i)
    pidx ∈ mc_pidx_set || continue
    lib_df.is_decoy[i] && continue
    pidx_to_key[pidx] = tn_pioneer_canonical_key(
        lib_df.sequence[i], lib_df.structural_mods[i], lib_df.prec_charge[i])
end

diann_df = CSV.read(TN_DIANN_PR_MATRIX, DataFrame; delim='\t')
metadata_cols = Set(["Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                     "First.Protein.Description", "Proteotypic",
                     "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge", "Precursor.Id"])
run_cols = [c for c in names(diann_df) if c ∉ metadata_cols]

function tn_diann_col_to_pioneer(col)
    m = match(r"(E\d+H\d+Y\d+)_\d+SPD_DIA_(\d+)", col)
    m === nothing && return nothing
    return "$(m.captures[1])_$(m.captures[2])"
end

pioneer_run_cols = filter(run_cols) do col
    pf = tn_diann_col_to_pioneer(col)
    pf !== nothing && pf ∈ Set(TN_PIONEER_FILES)
end

diann_fair_keys = Set{TN_CanonicalKey}()
for i in 1:nrow(diann_df)
    key = tn_parse_diann_modified_seq(diann_df[i, "Modified.Sequence"],
                                       Int(diann_df[i, "Precursor.Charge"]))
    for col in pioneer_run_cols
        val = diann_df[i, col]
        if !ismissing(val) && val isa Number && val > 0
            push!(diann_fair_keys, key); break
        end
    end
end
println("  DIA-NN fair keys: $(length(diann_fair_keys))")

# Annotate psm_combined with is_diann
is_diann = falses(nrow(psm_combined))
for i in 1:nrow(psm_combined)
    psm_combined.target[i] || continue
    key = get(pidx_to_key, UInt32(psm_combined.precursor_idx[i]), nothing)
    key !== nothing && key ∈ diann_fair_keys && (is_diann[i] = true)
end
psm_combined.is_diann = is_diann

# ── Feature definitions ──
const TN_BASELINE_FEATURES = [:scribe, :spectral_contrast, :city_block, :entropy_score, :y_count, :err_norm]
const TN_FRAGCORR_FEATURES = [:eigengap_raw, :lambda1_frac_raw, :ev1_raw, :median_corr_raw, :best_frag_corr_mean]

n_files = length(common_files)
top_n_logodds = max(1, floor(Int, sqrt(n_files)))

# ============================================================
# Run comparison for each fragment count
# ============================================================

# Storage for results across fragment counts
all_results = Dict{Int, Dict{String, NamedTuple}}()  # n_frags → condition → results
all_psm_qvals = Dict{Int, Dict{String, Vector{Float64}}}()  # PSM-level q-values

for nf in FRAG_COUNTS
    println("\n\n" * "╔" * "═"^56 * "╗")
    println("║  Fragment count: $nf" * " "^(56 - 20 - ndigits(nf)) * "║")
    println("╚" * "═"^56 * "╝")

    # Join PSM data with FRAGCORR scores for this fragment count
    fragcorr_df = all_corr_dfs[nf]
    mc_df = innerjoin(psm_combined, fragcorr_df;
                       on=[:precursor_idx, :target], makeunique=true)
    println("  Joined: $(nrow(psm_combined)) PSMs × $(nrow(fragcorr_df)) FRAGCORR → $(nrow(mc_df))")

    # Replace NaN with 0.0
    for col in vcat(TN_BASELINE_FEATURES, TN_FRAGCORR_FEATURES)
        if hasproperty(mc_df, col) && eltype(mc_df[!, col]) <: AbstractFloat
            mc_df[!, col] = replace(mc_df[!, col], NaN => 0.0)
        end
    end

    # Propagate is_diann from psm_combined
    if !hasproperty(mc_df, :is_diann)
        mc_df.is_diann = falses(nrow(mc_df))
    end

    augmented_features = vcat(TN_BASELINE_FEATURES, TN_FRAGCORR_FEATURES)
    is_target = mc_df.target

    # ── Train models ──
    conditions = Dict{String, Vector{Float64}}()

    println("  Training LightGBM baseline...")
    conditions["LightGBM baseline"] = tn_run_cv_lightgbm(mc_df, TN_BASELINE_FEATURES)

    println("  Training LightGBM augmented...")
    conditions["LightGBM augmented"] = tn_run_cv_lightgbm(mc_df, augmented_features)

    conditions["Pipeline prob (ref)"] = Float64.(mc_df.prob)

    # ── PSM-level q-values ──
    psm_qv = Dict{String, Vector{Float64}}()
    for (name, scores) in conditions
        psm_qv[name] = tn_compute_qvalues(scores, is_target; higher_is_better=true)
    end
    all_psm_qvals[nf] = psm_qv

    # ── Global prob simulation ──
    gp_res = Dict{String, NamedTuple}()
    for (name, scores) in conditions
        gp_df = tn_compute_global_probs(mc_df, scores, n_files)
        qvals = tn_compute_qvalues(gp_df.global_prob, gp_df.target; higher_is_better=true)

        thresholds = [0.01, 0.05, 0.10, 0.15]
        qval_counts = Dict{Float64, Int}()
        qval_diann  = Dict{Float64, Int}()
        for t in thresholds
            passing = qvals .<= t
            qval_counts[t] = count(passing .& gp_df.target)
            qval_diann[t]  = count(passing .& gp_df.target .& gp_df.is_diann)
        end

        gp_res[name] = (gp_df=gp_df, qvals=qvals, qval_counts=qval_counts, qval_diann=qval_diann)
    end
    all_results[nf] = gp_res
end

# ============================================================
# Summary tables
# ============================================================

thresholds = [0.01, 0.05, 0.10, 0.15]
thresh_strs = [string(round(t * 100, digits=0)) * "%" for t in thresholds]
cond_names = ["LightGBM baseline", "LightGBM augmented", "Pipeline prob (ref)"]

# Table 1: PSM-level targets at q-value thresholds
println("\n\n" * "═"^90)
println("TABLE 1: PSM-level target counts at q-value thresholds")
println("═"^90)
println(rpad("Frags", 8), rpad("Condition", 28), join(rpad.("q≤" .* thresh_strs, 12)))
println("─"^90)

for nf in FRAG_COUNTS
    for name in cond_names
        qv = all_psm_qvals[nf][name]
        is_tgt = nothing  # we need the is_target from this run
        # Re-derive counts from stored q-values
        # We stored q-values, so count targets directly
        # But we don't have is_target per nf variant... let's use a simpler approach
        # Actually the PSM-level is_target is the same for all, just the q-values differ
    end
end

# Simpler approach: just report global-level results (the key comparison)
println("\n\n" * "═"^100)
println("TABLE: Global-level unique target precursors at q-value thresholds")
println("═"^100)
println(rpad("Frags", 8), rpad("Condition", 28), join(rpad.("q≤" .* thresh_strs, 14)))
println("─"^100)

for nf in FRAG_COUNTS
    for name in cond_names
        res = all_results[nf][name]
        counts = [res.qval_counts[t] for t in thresholds]
        println(rpad(string(nf), 8), rpad(name, 28), join(rpad.(string.(counts), 14)))
    end
    println("─"^100)
end

# DIA-NN validation
println("\n\n" * "═"^100)
println("TABLE: DIA-NN validated targets at global q ≤ 1%")
println("═"^100)
println(rpad("Frags", 8), rpad("Condition", 28), rpad("Total", 12), rpad("DIA-NN", 12), "DIA-NN%")
println("─"^100)

for nf in FRAG_COUNTS
    for name in cond_names
        res = all_results[nf][name]
        total = res.qval_counts[0.01]
        diann = res.qval_diann[0.01]
        pct = total > 0 ? round(100.0 * diann / total, digits=1) : 0.0
        println(rpad(string(nf), 8), rpad(name, 28), rpad(string(total), 12),
                rpad(string(diann), 12), "$(pct)%")
    end
    println("─"^100)
end

# Venn: baseline vs augmented at q ≤ 1% for each fragment count
println("\n\n" * "═"^80)
println("TABLE: Venn Decomposition (LightGBM baseline vs augmented, global q ≤ 1%)")
println("═"^80)
println(rpad("Frags", 8), rpad("Category", 24), rpad("Count", 10), rpad("DIA-NN", 10), "DIA-NN%")
println("─"^80)

for nf in FRAG_COUNTS
    base_res = all_results[nf]["LightGBM baseline"]
    aug_res  = all_results[nf]["LightGBM augmented"]

    pass_base = Set{UInt32}()
    pass_aug  = Set{UInt32}()
    diann_base = Set{UInt32}()
    diann_aug  = Set{UInt32}()

    gp_b = base_res.gp_df; qv_b = base_res.qvals
    gp_a = aug_res.gp_df;  qv_a = aug_res.qvals

    for i in 1:nrow(gp_b)
        gp_b.target[i] || continue
        if qv_b[i] <= 0.01
            push!(pass_base, gp_b.precursor_idx[i])
            gp_b.is_diann[i] && push!(diann_base, gp_b.precursor_idx[i])
        end
    end
    for i in 1:nrow(gp_a)
        gp_a.target[i] || continue
        if qv_a[i] <= 0.01
            push!(pass_aug, gp_a.precursor_idx[i])
            gp_a.is_diann[i] && push!(diann_aug, gp_a.precursor_idx[i])
        end
    end

    both     = intersect(pass_base, pass_aug)
    aug_only = setdiff(pass_aug, pass_base)
    base_only = setdiff(pass_base, pass_aug)

    for (label, s, ds) in [("Both", both, intersect(both, diann_base ∪ diann_aug)),
                            ("Aug only (FRAGCORR)", aug_only, intersect(aug_only, diann_aug)),
                            ("Base only", base_only, intersect(base_only, diann_base))]
        n = length(s); nd = length(ds)
        pct = n > 0 ? round(100.0 * nd / n, digits=1) : 0.0
        println(rpad(string(nf), 8), rpad(label, 24), rpad(string(n), 10),
                rpad(string(nd), 10), "$(pct)%")
    end
    println("─"^80)
end

# ============================================================
# Plots
# ============================================================

println("\nGenerating plots...")

n_points = 200
qv_range = range(0.0, 0.15, length=n_points)

frag_colors = Dict(5 => :blue, 6 => :green, 7 => :red)
frag_styles = Dict("LightGBM baseline" => :dash, "LightGBM augmented" => :solid,
                    "Pipeline prob (ref)" => :dot)

# Plot 1: Global-level targets for LightGBM augmented across fragment counts
p_nfrag = plot(xlabel="q-value (global)", ylabel="Unique target precursors",
               title="FRAGCORR: top-N fragments × LightGBM augmented (global probs)",
               legend=:topleft, size=(900, 550))

for nf in FRAG_COUNTS
    for name in ["LightGBM baseline", "LightGBM augmented"]
        res = all_results[nf][name]
        gp_df = res.gp_df; qvals = res.qvals
        counts = [count((qvals .<= t) .& gp_df.target) for t in qv_range]
        plot!(p_nfrag, collect(qv_range), counts;
              label="$(nf) frags — $name", lw=2.5,
              color=frag_colors[nf], ls=frag_styles[name])
    end
end
display(p_nfrag)
savefig(p_nfrag, "p_topn_global_targets.png")

# Plot 2: DIA-NN validated for LightGBM augmented across fragment counts
p_nfrag_diann = plot(xlabel="q-value (global)", ylabel="DIA-NN validated unique targets",
                     title="FRAGCORR: top-N fragments × LightGBM augmented (DIA-NN validated)",
                     legend=:topleft, size=(900, 550))

for nf in FRAG_COUNTS
    for name in ["LightGBM baseline", "LightGBM augmented"]
        res = all_results[nf][name]
        gp_df = res.gp_df; qvals = res.qvals
        counts = [count((qvals .<= t) .& gp_df.target .& gp_df.is_diann) for t in qv_range]
        plot!(p_nfrag_diann, collect(qv_range), counts;
              label="$(nf) frags — $name", lw=2.5,
              color=frag_colors[nf], ls=frag_styles[name])
    end
end
display(p_nfrag_diann)
savefig(p_nfrag_diann, "p_topn_global_diann.png")

# Plot 3: Augmented gain (aug - baseline) at each fragment count
p_gain = plot(xlabel="q-value (global)", ylabel="Augmented - Baseline targets",
              title="FRAGCORR gain from augmentation by fragment count",
              legend=:topleft, size=(900, 550))

for nf in FRAG_COUNTS
    base_res = all_results[nf]["LightGBM baseline"]
    aug_res  = all_results[nf]["LightGBM augmented"]
    gp_b = base_res.gp_df; qv_b = base_res.qvals
    gp_a = aug_res.gp_df;  qv_a = aug_res.qvals

    gains = Int[]
    for t in qv_range
        nb = count((qv_b .<= t) .& gp_b.target)
        na = count((qv_a .<= t) .& gp_a.target)
        push!(gains, na - nb)
    end
    plot!(p_gain, collect(qv_range), gains;
          label="$(nf) frags", lw=2.5, color=frag_colors[nf])
end
hline!(p_gain, [0]; color=:gray, ls=:dot, label=nothing)
display(p_gain)
savefig(p_gain, "p_topn_augmented_gain.png")

println("\nPlots saved:")
println("  p_topn_global_targets.png")
println("  p_topn_global_diann.png")
println("  p_topn_augmented_gain.png")

println("\nPlot objects available: p_nfrag, p_nfrag_diann, p_gain")
println("Results dict: all_results[n_frags][condition_name]")
println("\nDone! Top-N fragment comparison complete.")
