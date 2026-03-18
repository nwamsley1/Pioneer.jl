#=
FRAGCORR Chromatogram Viewer
=============================
Interactive plotting of fragment chromatograms from the focused FRAGCORR pass.
Each plot shows 6 fragment traces annotated with target/decoy status, sequence,
scoring metrics, DIA-NN validation + RT, and Pioneer global scoring status.

Supports B2 precursor browsing: DIA-NN-detected precursors that had first-pass
PSMs but were filtered at the first-pass hybrid PEP/q-value stage.

Usage:
    include("scripts/fragcorr_plot_chromatograms.jl")

    # B2 precursors (DIA-NN YES, first-pass filtered)
    b2 = b2_precursors(cp_data)
    plot_chrom_pair(cp_data, b2, 1)    # Browse by first-pass PEP (best first)
    b2_grid(cp_data; n=12)             # Grid view

    # Plot a single precursor by index into the grouped data
    plot_chrom(cp_data, 50)

    # Browse the top-scoring targets
    plot_top_targets(cp_data; n=12)

    # DIA-NN identified but Pioneer global scoring REJECTED
    missed = diann_not_pioneer(cp_data)
    plot_chrom_pair(cp_data, missed, 1)
=#

using Arrow, DataFrames, Tables, CSV, Parquet2
using Statistics, LinearAlgebra
using Plots

# ============================================================
# Constants
# ============================================================

const CP_RESULTS_DIR = "/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/OlsenEclipse_diagnostic"
const CP_DATA_DIR = joinpath(CP_RESULTS_DIR, "temp_data")
const CP_CORR_DIR = joinpath(CP_DATA_DIR, "first_pass_corr")
const CP_PSM_DIR  = joinpath(CP_DATA_DIR, "first_pass_psms")
const CP_LIBRARY_PATH = "/Users/nathanwamsley/Data/SPEC_LIBS/3P_rank_based.poin"
const CP_DIANN_REPORT = joinpath("/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse",
    "DIA-NN_results", "OlsenExplorisThreeProteome500ng-11-24-2025-report.parquet")
const CP_PIONEER_GLOBAL = joinpath(CP_RESULTS_DIR, "precursors_long.arrow")
const CP_PIONEER_FILES = ["E45H50Y5_2", "E45H50Y5_3", "E45H50Y5_4",
                          "E5H50Y45_1", "E5H50Y45_2", "E5H50Y45_4"]
const CP_INTENSITY_COLS = [:intensity_1, :intensity_2, :intensity_3,
                           :intensity_4, :intensity_5, :intensity_6]
const CP_FRAG_COLORS = [1, 2, 3, 4, 5, 6]

# ============================================================
# DIA-NN / Pioneer Key Helpers
# ============================================================

const _CP_CanonicalKey = Tuple{String, Vector{Tuple{Int,Int}}, Int}

function _cp_parse_diann_seq(mod_seq::AbstractString, charge::Int)::_CP_CanonicalKey
    stripped = IOBuffer()
    mods = Tuple{Int,Int}[]
    pos = 0
    i = 1
    while i <= ncodeunits(mod_seq)
        c = mod_seq[i]
        if c == '('
            close_idx = findnext(')', mod_seq, i)
            mod_str = mod_seq[i+1:close_idx-1]
            colon_idx = findlast(':', mod_str)
            unimod_id = parse(Int, mod_str[colon_idx+1:end])
            push!(mods, (pos, unimod_id))
            i = close_idx + 1
        else
            pos += 1
            write(stripped, c)
            i += 1
        end
    end
    sort!(mods)
    return (String(take!(stripped)), mods, charge)
end

function _cp_pioneer_key(seq::AbstractString, structural_mods::AbstractString, charge::Integer)::_CP_CanonicalKey
    mods = Tuple{Int,Int}[]
    for m in eachmatch(r"\((\d+),\w,(\w+:\d+)\)", structural_mods)
        pos = parse(Int, m.captures[1])
        mod_str = m.captures[2]
        colon_idx = findlast(':', mod_str)
        unimod_id = parse(Int, mod_str[colon_idx+1:end])
        push!(mods, (pos, unimod_id))
    end
    sort!(mods)
    return (String(seq), mods, Int(charge))
end

"""Map a DIA-NN Run string to a Pioneer file stem, or nothing."""
function _cp_diann_run_to_pioneer(run::AbstractString)::Union{String, Nothing}
    m = match(r"(E\d+H\d+Y\d+)_\d+SPD_DIA_(\d+)", run)
    m === nothing && return nothing
    return "$(m.captures[1])_$(m.captures[2])"
end

# ============================================================
# Global PEP Reconstruction (mirrors Pioneer pipeline)
# ============================================================

"""Log-odds average of top-N probabilities → probability. Mirrors Pioneer's _logodds_combine."""
function _cp_logodds_combine(probs::Vector{Float32}, top_n::Int)::Float32
    isempty(probs) && return 0.0f0
    n = min(length(probs), top_n)
    sorted = sort(probs; rev=true)
    selected = @view sorted[1:n]
    eps = 1f-6
    lo = log.(clamp.(selected, 0.1f0, 1 - eps) ./ (1 .- clamp.(selected, 0.1f0, 1 - eps)))
    avg = sum(lo) / n
    return 1.0f0 / (1 + exp(-avg))
end

"""Weighted PAVA for isotonic regression. Mirrors Pioneer's _weighted_pava."""
function _cp_weighted_pava(y::Vector{Float64}, w::Vector{Float64})
    n = length(y)
    v    = Vector{Float64}(undef, n)
    wt   = Vector{Float64}(undef, n)
    lenv = Vector{Int}(undef, n)
    m = 0
    for i in 1:n
        m += 1
        @inbounds begin
            v[m]    = y[i]
            wt[m]   = w[i]
            lenv[m] = 1
        end
        while m > 1 && v[m-1] > v[m]
            @inbounds begin
                new_w   = wt[m-1] + wt[m]
                new_v   = (v[m-1]*wt[m-1] + v[m]*wt[m]) / new_w
                new_len = lenv[m-1] + lenv[m]
                m -= 1
                v[m]    = new_v
                wt[m]   = new_w
                lenv[m] = new_len
            end
        end
    end
    result = Vector{Float64}(undef, n)
    idx = 1
    for j in 1:m
        @inbounds for _ in 1:lenv[j]
            result[idx] = v[j]
            idx += 1
        end
    end
    return result
end

"""
    _cp_reconstruct_global_pep(psm_df) -> Dict{UInt32, Float64}

Reconstruct the global PEP for each precursor_idx from first-pass PSM probs.
Mirrors Pioneer's get_best_precursors_accross_runs logic:
  1. Best prob per precursor per file
  2. Log-odds combine top-sqrt(n_files)
  3. Isotonic regression (PAVA) → calibrated PEP
"""
function _cp_reconstruct_global_pep(psm_df::DataFrame)
    # Step 1: collect best prob per precursor per file
    # source_file acts as the file identifier
    files = unique(psm_df.source_file)
    n_files = length(files)
    top_n = max(1, floor(Int, sqrt(n_files)))

    prec_probs_by_file = Dict{UInt32, Vector{Float32}}()
    prec_is_target = Dict{UInt32, Bool}()

    for f in files
        file_mask = psm_df.source_file .== f
        file_df = psm_df[file_mask, :]
        # Best prob per precursor in this file
        file_best = Dict{UInt32, Float32}()
        for i in 1:nrow(file_df)
            pidx = UInt32(file_df.precursor_idx[i])
            p = Float32(file_df.prob[i])
            if !haskey(file_best, pidx) || p > file_best[pidx]
                file_best[pidx] = p
            end
            if !haskey(prec_is_target, pidx)
                prec_is_target[pidx] = file_df.target[i]
            end
        end
        for (pidx, p) in file_best
            if haskey(prec_probs_by_file, pidx)
                push!(prec_probs_by_file[pidx], p)
            else
                prec_probs_by_file[pidx] = Float32[p]
            end
        end
    end

    # Step 2: log-odds combine
    n_unique = length(prec_probs_by_file)
    pidxs = Vector{UInt32}(undef, n_unique)
    global_probs = Vector{Float32}(undef, n_unique)
    targets = Vector{Bool}(undef, n_unique)

    for (i, (pidx, probs)) in enumerate(prec_probs_by_file)
        pidxs[i] = pidx
        global_probs[i] = _cp_logodds_combine(probs, top_n)
        targets[i] = get(prec_is_target, pidx, false)
    end

    # Step 3: isotonic regression → PEP
    order = sortperm(global_probs; rev=true)
    N = n_unique

    labels  = Vector{Float64}(undef, N + 1)
    weights = Vector{Float64}(undef, N + 1)
    labels[1]  = 0.0  # pseudo observation
    weights[1] = 1.0
    @inbounds for j in 1:N
        idx = order[j]
        labels[j+1]  = targets[idx] ? 0.0 : 1.0
        weights[j+1] = targets[idx] ? 1.0 : 1.0
    end

    fitted = _cp_weighted_pava(labels, weights)
    fitted = fitted[2:end]  # remove pseudo row
    pep_vals = clamp.(fitted ./ (1.0 .- fitted), 0.0, 1.0)

    # Build lookup
    result = Dict{UInt32, Float64}()
    @inbounds for (j, idx) in enumerate(order)
        result[pidxs[idx]] = pep_vals[j]
    end
    return result
end

# ============================================================
# Data Loading
# ============================================================

function cp_load_all_data()
    println("Loading chromatogram data...")

    # --- Load raw chromatograms ---
    corr_files = filter(f -> endswith(f, ".arrow"), readdir(CP_CORR_DIR))
    psm_files  = Set(filter(f -> endswith(f, ".arrow"), readdir(CP_PSM_DIR)))
    corr_files = filter(f -> f in psm_files, corr_files)

    all_corr = DataFrame[]
    all_psm  = DataFrame[]
    for f in corr_files
        println("  Loading $f ...")
        corr = DataFrame(Tables.columntable(Arrow.Table(joinpath(CP_CORR_DIR, f))))
        psm  = DataFrame(Tables.columntable(Arrow.Table(joinpath(CP_PSM_DIR, f))))
        corr.source_file .= f
        psm.source_file  .= f
        push!(all_corr, corr)
        push!(all_psm, psm)
    end

    corr_df = vcat(all_corr...; cols=:union)
    psm_df  = vcat(all_psm...; cols=:union)
    sort!(corr_df, [:precursor_idx, :scan_idx])

    println("  Chromatogram rows: $(nrow(corr_df))")
    println("  PSM rows: $(nrow(psm_df))")

    # --- Load library for sequence/charge ---
    println("  Loading precursor library...")
    lib_df = DataFrame(Arrow.Table(joinpath(CP_LIBRARY_PATH, "precursors_table.arrow")))

    pidx_set = Set{UInt32}(corr_df.precursor_idx)
    pidx_meta = Dict{UInt32, NamedTuple{(:sequence, :charge, :mods, :is_decoy, :species),
                                         Tuple{String, Int, String, Bool, String}}}()
    for i in 1:nrow(lib_df)
        pidx = UInt32(i)
        pidx in pidx_set || continue
        smods = ismissing(lib_df.structural_mods[i]) ? "" : lib_df.structural_mods[i]
        sp = ismissing(lib_df.proteome_identifiers[i]) ? "?" : String(lib_df.proteome_identifiers[i])
        pidx_meta[pidx] = (sequence=lib_df.sequence[i],
                           charge=Int(lib_df.prec_charge[i]),
                           mods=smods,
                           is_decoy=lib_df.is_decoy[i],
                           species=sp)
    end
    println("  Library entries matched: $(length(pidx_meta))")

    # --- Build canonical_key -> precursor_idx reverse lookup (targets only) ---
    key_to_pidxs = Dict{_CP_CanonicalKey, Vector{UInt32}}()
    for (pidx, meta) in pidx_meta
        meta.is_decoy && continue
        k = _cp_pioneer_key(meta.sequence, meta.mods, meta.charge)
        push!(get!(Vector{UInt32}, key_to_pidxs, k), pidx)
    end

    # --- Load DIA-NN report parquet for per-run RTs ---
    println("  Loading DIA-NN report parquet (for retention times)...")
    diann_report = DataFrame(Parquet2.Dataset(CP_DIANN_REPORT); copycols=false)
    println("  DIA-NN report rows: $(nrow(diann_report))")

    # Build lookup: (precursor_idx, pioneer_file) -> DIA-NN RT
    # Also track fair-key set (detected in Pioneer runs)
    pioneer_file_set = Set(CP_PIONEER_FILES)
    diann_fair_keys = Set{_CP_CanonicalKey}()
    # (canonical_key, pioneer_file) -> RT
    diann_key_file_rt = Dict{Tuple{_CP_CanonicalKey, String}, Float64}()

    for i in 1:nrow(diann_report)
        pf = _cp_diann_run_to_pioneer(diann_report[i, "Run"])
        pf === nothing && continue
        pf in pioneer_file_set || continue

        qty = diann_report[i, "Precursor.Quantity"]
        (ismissing(qty) || !(qty isa Number) || qty <= 0) && continue

        key = _cp_parse_diann_seq(diann_report[i, "Modified.Sequence"],
                                   Int(diann_report[i, "Precursor.Charge"]))
        push!(diann_fair_keys, key)
        diann_key_file_rt[(key, pf)] = Float64(diann_report[i, "RT"])
    end
    println("  DIA-NN fair keys: $(length(diann_fair_keys))")
    println("  DIA-NN (key, file) -> RT entries: $(length(diann_key_file_rt))")

    # Build (precursor_idx, source_file) -> DIA-NN RT lookup
    pidx_file_diann_rt = Dict{Tuple{UInt32, String}, Float64}()
    for ((key, pf), rt) in diann_key_file_rt
        pidxs = get(key_to_pidxs, key, UInt32[])
        for pidx in pidxs
            src = pf * ".arrow"
            pidx_file_diann_rt[(pidx, src)] = rt
        end
    end
    println("  (precursor_idx, file) -> DIA-NN RT entries: $(length(pidx_file_diann_rt))")

    # Build precursor_idx -> is_diann (across any file)
    pidx_is_diann = Dict{UInt32, Bool}()
    for (pidx, meta) in pidx_meta
        meta.is_decoy && continue
        k = _cp_pioneer_key(meta.sequence, meta.mods, meta.charge)
        pidx_is_diann[pidx] = k in diann_fair_keys
    end

    # --- Load first-pass global scores (before hybrid filter) ---
    fp_pep_lookup = Dict{UInt32, Float32}()
    fp_prob_lookup = Dict{UInt32, Float32}()
    fp_pidx_set = Set{UInt32}()
    fp_scores_path = joinpath(CP_DATA_DIR, "first_pass_global_scores", "global_scores.arrow")
    if isfile(fp_scores_path)
        println("  Loading first-pass global scores...")
        fp_global_scores = DataFrame(Tables.columntable(Arrow.Table(fp_scores_path)))
        for i in 1:nrow(fp_global_scores)
            pidx = UInt32(fp_global_scores.precursor_idx[i])
            push!(fp_pidx_set, pidx)
            fp_pep_lookup[pidx] = fp_global_scores.global_pep[i]
            fp_prob_lookup[pidx] = fp_global_scores.global_prob[i]
        end
        println("    First-pass precursors: $(length(fp_pidx_set))")
    else
        println("  WARNING: first_pass_global_scores not found — B2 classification unavailable")
    end

    # --- Classify B2 precursors by reconstructing the hybrid filter ---
    # A precursor is truly B2 if it was *filtered out* by the first-pass hybrid
    # PEP/q-value filter. We reconstruct that filter from the dump data rather
    # than checking `pre_filter_global_scores`, because a precursor can pass the
    # first-pass filter yet still be absent from the second-pass scoring dump
    # (different parameters/tolerances → no PSM produced).
    b2_pidxs = Set{UInt32}()
    if !isempty(fp_pidx_set)
        # Reconstruct hybrid filter from first-pass dump
        fp_n = nrow(fp_global_scores)
        fp_order = sortperm(fp_global_scores.global_pep)

        # Q-value floor: last position where cumulative D/T ≤ 0.15
        cum_t, cum_d, n_qvalue_based = 0, 0, 0
        for (rank, i) in enumerate(fp_order)
            if fp_global_scores.target[i]
                cum_t += 1
            else
                cum_d += 1
            end
            if cum_t > 0 && cum_d / cum_t <= 0.15
                n_qvalue_based = rank
            end
        end

        n_passing_threshold = count(fp_global_scores.global_pep .<= 0.5f0)
        n_min_floor = min(50_000, fp_n)
        n_to_keep = max(n_passing_threshold, n_qvalue_based, n_min_floor)

        # Precursors beyond n_to_keep are truly filtered at first pass
        fp_filtered_set = Set{UInt32}()
        if n_to_keep < fp_n
            for rank in (n_to_keep + 1):fp_n
                push!(fp_filtered_set, UInt32(fp_global_scores.precursor_idx[fp_order[rank]]))
            end
        end

        # B2 = DIA-NN detected AND truly filtered at first pass
        for pidx in fp_filtered_set
            get(pidx_is_diann, pidx, false) || continue
            push!(b2_pidxs, pidx)
        end
        println("  Hybrid filter reconstruction: keeping $n_to_keep / $fp_n")
        println("  B2 precursors (DIA-NN YES, truly first-pass filtered): $(length(b2_pidxs))")
    end

    # --- Load Pioneer global scoring output ---
    pidx_file_global_qval = Dict{Tuple{UInt32, String}, Float64}()
    if isfile(CP_PIONEER_GLOBAL)
        println("  Loading Pioneer global scoring (precursors_long.arrow)...")
        global_df = DataFrame(Tables.columntable(Arrow.Table(CP_PIONEER_GLOBAL)))
        println("  Pioneer global rows: $(nrow(global_df))")

        # Build (precursor_idx, source_file) -> global_qval lookup
        for i in 1:nrow(global_df)
            global_df.target[i] || continue
            pidx = UInt32(global_df.precursor_idx[i])
            fn = global_df.file_name[i]
            src = fn * ".arrow"
            qv = global_df.MBR_boosted_global_qval[i]
            if !ismissing(qv)
                pidx_file_global_qval[(pidx, src)] = Float64(qv)
            end
        end
        println("  (precursor_idx, file) -> global_qval entries: $(length(pidx_file_global_qval))")
    else
        println("  WARNING: precursors_long.arrow not found — global scoring status unavailable")
    end

    # --- Group chromatograms ---
    gdf = groupby(corr_df, [:precursor_idx, :source_file])

    # --- Build best-PSM lookup ---
    sort!(psm_df, :prob; rev=true)
    psm_best = combine(groupby(psm_df, [:precursor_idx, :source_file]), first)

    # --- Reconstruct global PEP from first-pass probs ---
    # Mirrors Pioneer's getBestPrecursorsAccrossRuns.jl logic:
    #   1. Best prob per precursor per file
    #   2. Log-odds combine top-sqrt(n_files)
    #   3. Isotonic regression (PAVA) → PEP
    println("  Reconstructing global PEP from first-pass probs...")
    recon_global_pep = _cp_reconstruct_global_pep(psm_df)
    println("  Reconstructed PEP for $(length(recon_global_pep)) precursors")

    println("\nReady! $(length(gdf)) chromatogram groups loaded.")
    return (corr_gdf=gdf, psm_best=psm_best, pidx_meta=pidx_meta,
            pidx_is_diann=pidx_is_diann, corr_df=corr_df,
            pidx_file_diann_rt=pidx_file_diann_rt,
            pidx_file_global_qval=pidx_file_global_qval,
            recon_global_pep=recon_global_pep,
            fp_pep_lookup=fp_pep_lookup, fp_prob_lookup=fp_prob_lookup,
            b2_pidxs=b2_pidxs)
end

# ============================================================
# Eigenvalue / Correlation Computation
# ============================================================

"""
    _cp_extract_intensity_matrix(chrom) -> Matrix{Float64}

Extract the n×6 intensity matrix from a chromatogram group, clamped to ≥ 0.
"""
function _cp_extract_intensity_matrix(chrom)
    n = nrow(chrom)
    p = length(CP_INTENSITY_COLS)
    X = Matrix{Float64}(undef, n, p)
    @inbounds for (j, col) in enumerate(CP_INTENSITY_COLS)
        copyto!(view(X, :, j), chrom[!, col])
    end
    X .= max.(X, 0.0)
    return X
end

"""
    _cp_pearson_corrmat(X; eps=1e-12) -> Matrix{Float64}

Compute the p×p Pearson correlation matrix from an n×p data matrix.
"""
function _cp_pearson_corrmat(X::Matrix{Float64}; eps=1e-12)
    n, p = size(X)
    means = vec(mean(X; dims=1))
    stds  = Vector{Float64}(undef, p)
    @inbounds for j in 1:p
        ss = 0.0
        for i in 1:n
            ss += (X[i, j] - means[j])^2
        end
        stds[j] = sqrt(ss / max(n - 1, 1))
    end

    C = Matrix{Float64}(undef, p, p)
    @inbounds for j1 in 1:p
        C[j1, j1] = 1.0
        for j2 in (j1+1):p
            if stds[j1] > eps && stds[j2] > eps
                cov_val = 0.0
                for i in 1:n
                    cov_val += (X[i, j1] - means[j1]) * (X[i, j2] - means[j2])
                end
                c = cov_val / (max(n - 1, 1) * stds[j1] * stds[j2])
                C[j1, j2] = c
                C[j2, j1] = c
            else
                C[j1, j2] = 0.0
                C[j2, j1] = 0.0
            end
        end
    end
    return C
end

"""
    _cp_rank_matrix(X) -> Matrix{Float64}

Replace each column of X with average ranks (handling ties).
"""
function _cp_rank_matrix(X::Matrix{Float64})
    n, p = size(X)
    R = Matrix{Float64}(undef, n, p)
    @inbounds for j in 1:p
        col = view(X, :, j)
        ord = sortperm(col)
        i = 1
        while i <= n
            # Find run of tied values
            val = col[ord[i]]
            k = i
            while k <= n && col[ord[k]] == val
                k += 1
            end
            avg_rank = (i + k - 1) / 2.0
            for t in i:(k-1)
                R[ord[t], j] = avg_rank
            end
            i = k
        end
    end
    return R
end

"""
    _cp_eigen_from_corrmat(C; eps=1e-12) -> NamedTuple

Compute eigenvalues and derived scores from a correlation matrix.
"""
function _cp_eigen_from_corrmat(C::Matrix{Float64}; eps=1e-12)
    p = size(C, 1)
    evals = eigvals(Symmetric(C))
    sort!(evals; rev=true)
    tr = sum(evals)

    lambda1_frac = tr > eps ? evals[1] / tr : 0.0
    eigengap     = p >= 2 ? evals[1] - evals[2] : 0.0

    corrs = Float64[]
    @inbounds for j1 in 1:p, j2 in (j1+1):p
        push!(corrs, C[j1, j2])
    end
    med_corr = isempty(corrs) ? 0.0 : median(corrs)

    return (eigenvalues=evals, corr_matrix=C,
            lambda1_frac=lambda1_frac, eigengap=eigengap, median_corr=med_corr)
end

"""
    compute_corr_eigenvalues(chrom)

Build the 6×6 Pearson correlation matrix from fragment traces and return
eigenvalues (descending), correlation matrix, and derived scores.
"""
function compute_corr_eigenvalues(chrom; eps=1e-12)
    X = _cp_extract_intensity_matrix(chrom)
    C = _cp_pearson_corrmat(X; eps)
    return _cp_eigen_from_corrmat(C; eps)
end

"""
    compute_spearman_eigenvalues(chrom)

Build the 6×6 Spearman (rank) correlation matrix from fragment traces and
return eigenvalues (descending), correlation matrix, and derived scores.
"""
function compute_spearman_eigenvalues(chrom; eps=1e-12)
    X = _cp_extract_intensity_matrix(chrom)
    R = _cp_rank_matrix(X)
    C = _cp_pearson_corrmat(R; eps)
    return _cp_eigen_from_corrmat(C; eps)
end

# ============================================================
# Plotting
# ============================================================

"""
    plot_chrom(data, idx; show_points=true)

Plot the chromatogram for the `idx`-th group in `data.corr_gdf`.

Left panel: 6 fragment traces. If DIA-NN identified this precursor in this
file, a vertical dashed green line marks DIA-NN's retention time.

Right panel: eigenvalue bar chart of the 6x6 correlation matrix with
lambda1_frac, eigengap, and median_corr.

Title: sequence, charge, target/decoy, prob, scribe, q-value, DIA-NN status,
and Pioneer global scoring status.
"""
function plot_chrom(data, idx::Int; show_points::Bool=true)
    gdf = data.corr_gdf
    1 <= idx <= length(gdf) || error("idx must be in 1:$(length(gdf))")

    chrom = gdf[idx]
    pidx  = UInt32(first(chrom.precursor_idx))
    src   = first(chrom.source_file)
    is_target = first(chrom.target)

    # Library metadata
    meta = get(data.pidx_meta, pidx, nothing)
    seq_str = meta !== nothing ? meta.sequence : "?"
    charge  = meta !== nothing ? meta.charge : 0
    species = meta !== nothing ? meta.species : "?"

    seq_display = length(seq_str) > 25 ? seq_str[1:22] * "..." : seq_str

    # PSM scores
    psm_row = data.psm_best[
        (data.psm_best.precursor_idx .== pidx) .& (data.psm_best.source_file .== src), :]
    if nrow(psm_row) > 0
        prob_val   = round(Float64(psm_row.prob[1]); digits=4)
        scribe_val = round(Float64(psm_row.scribe[1]); digits=3)
        qval       = round(Float64(psm_row.q_value[1]); digits=4)
    else
        prob_val = scribe_val = qval = NaN
    end

    # DIA-NN status + RT
    diann_detected = get(data.pidx_is_diann, pidx, false)
    diann_rt = get(data.pidx_file_diann_rt, (pidx, src), nothing)
    diann_str = if diann_rt !== nothing
        "YES (RT=$(round(diann_rt; digits=2)))"
    elseif diann_detected
        "YES (other file)"
    else
        "no"
    end

    # Pioneer global scoring status
    global_qv = get(data.pidx_file_global_qval, (pidx, src), nothing)
    if global_qv !== nothing
        global_str = global_qv <= 0.01 ? "PASS (q=$(round(global_qv; digits=4)))" : "FAIL (q=$(round(global_qv; digits=4)))"
    else
        global_str = is_target ? "FILTERED" : "n/a"
    end

    # First-pass global PEP from diagnostic dump (actual value used for filtering)
    fp_pep = get(data.fp_pep_lookup, pidx, nothing)
    fp_pep_str = fp_pep !== nothing ? "$(round(Float64(fp_pep); digits=4))" : "n/a"

    # Reconstructed global PEP from first-pass log-odds + PAVA
    recon_pep = get(data.recon_global_pep, pidx, nothing)
    pep_str = recon_pep !== nothing ? "$(round(recon_pep; digits=4))" : "n/a"

    # B2 status
    is_b2 = pidx in data.b2_pidxs
    b2_str = is_b2 ? " [B2]" : ""

    # Scan count
    n_scans = nrow(chrom)

    # Compute eigenvalues (Pearson and Spearman)
    pearson_result  = compute_corr_eigenvalues(chrom)
    spearman_result = compute_spearman_eigenvalues(chrom)
    p_evals = pearson_result.eigenvalues
    s_evals = spearman_result.eigenvalues

    # --- Left panel: chromatogram traces ---
    td_str = is_target ? "TARGET" : "DECOY"
    title_chrom = "$seq_display $(charge)+ | $species | $td_str | $(n_scans) scans$b2_str\n" *
                  "prob=$prob_val  scribe=$scribe_val  q=$qval\n" *
                  "DIA-NN=$diann_str  Global=$global_str  fpPEP=$fp_pep_str"

    p_chrom = plot(; title=title_chrom, titlefontsize=8,
                   xlabel="RT (min)", ylabel="Intensity",
                   legend=:topright, legendfontsize=7)

    rt = chrom[!, :rt]
    for (j, col) in enumerate(CP_INTENSITY_COLS)
        intensities = chrom[!, col]
        lab = "frag $j"
        if show_points
            scatter!(p_chrom, rt, intensities; label=nothing,
                     color=CP_FRAG_COLORS[j], markersize=3, markerstrokewidth=0, alpha=0.7)
        end
        plot!(p_chrom, rt, intensities; label=lab, color=CP_FRAG_COLORS[j], lw=1.5)
    end

    # Draw DIA-NN RT vertical line if available
    if diann_rt !== nothing
        vline!(p_chrom, [diann_rt]; color=:green, lw=2, ls=:dash,
               label="DIA-NN RT")
    end

    # --- Right panel: Pearson vs Spearman eigenvalue bar chart ---
    n_ev = length(p_evals)
    title_eigen = "Eigenvalues: Pearson vs Spearman\n" *
                  "P: λ1f=$(round(pearson_result.lambda1_frac; digits=3)) " *
                  "gap=$(round(pearson_result.eigengap; digits=3)) " *
                  "med=$(round(pearson_result.median_corr; digits=3))\n" *
                  "S: λ1f=$(round(spearman_result.lambda1_frac; digits=3)) " *
                  "gap=$(round(spearman_result.eigengap; digits=3)) " *
                  "med=$(round(spearman_result.median_corr; digits=3))"

    max_eval = max(maximum(p_evals), maximum(s_evals))
    min_eval = min(minimum(p_evals), minimum(s_evals))

    # Grouped bars: Pearson left, Spearman right at each component position
    x_pearson  = (1:n_ev) .- 0.17
    x_spearman = (1:n_ev) .+ 0.17

    p_eigen = bar(x_pearson, p_evals;
                  title=title_eigen, titlefontsize=8,
                  xlabel="Component", ylabel="Eigenvalue",
                  color=:steelblue, label="Pearson", alpha=0.8,
                  ylims=(min(0.0, min_eval - 0.1), max_eval * 1.2),
                  bar_width=0.3, xticks=1:n_ev, legendfontsize=7)
    bar!(p_eigen, x_spearman, s_evals;
         color=:darkorange, label="Spearman", alpha=0.8, bar_width=0.3)

    # Annotate top two eigenvalues
    annotate!(p_eigen, x_pearson[1], p_evals[1] + max_eval*0.04,
              text("$(round(p_evals[1]; digits=2))", 7, :center, :steelblue))
    annotate!(p_eigen, x_spearman[1], s_evals[1] + max_eval*0.04,
              text("$(round(s_evals[1]; digits=2))", 7, :center, :darkorange))
    if n_ev >= 2
        annotate!(p_eigen, x_pearson[2], p_evals[2] + max_eval*0.04,
                  text("$(round(p_evals[2]; digits=2))", 7, :center, :steelblue))
        annotate!(p_eigen, x_spearman[2], s_evals[2] + max_eval*0.04,
                  text("$(round(s_evals[2]; digits=2))", 7, :center, :darkorange))
    end

    # Combine
    p = plot(p_chrom, p_eigen; layout=grid(1, 2, widths=[0.6, 0.4]),
             size=(1100, 450), margin=4Plots.mm)

    return p
end

# ============================================================
# Lookup helper: find group index for a (precursor_idx, source_file)
# ============================================================

function _cp_find_group(data, pidx::UInt32, src::String)
    for (gi, g) in enumerate(data.corr_gdf)
        if UInt32(first(g.precursor_idx)) == pidx && first(g.source_file) == src
            return gi
        end
    end
    return nothing
end

# ============================================================
# Convenience Plotting Functions
# ============================================================

"""
    plot_chrom_by_pidx(data, precursor_idx; file=nothing)

Plot chromatogram for a specific precursor_idx. If `file` is given, use that
source file; otherwise pick the one with the highest prob.
"""
function plot_chrom_by_pidx(data, precursor_idx::Integer; file::Union{Nothing,String}=nothing)
    pidx = UInt32(precursor_idx)
    if file === nothing
        rows = data.psm_best[data.psm_best.precursor_idx .== pidx, :]
        nrow(rows) == 0 && error("precursor_idx $pidx not found in PSM data")
        best_row = argmax(rows.prob)
        file = rows.source_file[best_row]
    end
    gi = _cp_find_group(data, pidx, file)
    gi === nothing && error("No chromatogram found for precursor_idx=$pidx, file=$file")
    return plot_chrom(data, gi)
end

"""Build a grid plot from a list of (pidx, source_file) pairs."""
function _cp_plot_grid(data, pairs::Vector{Tuple{UInt32, String}};
                       ncols::Int=3, show_points::Bool=false)
    plots = Plots.Plot[]
    for (pidx, src) in pairs
        gi = _cp_find_group(data, pidx, src)
        gi === nothing && continue
        push!(plots, plot_chrom(data, gi; show_points=show_points))
    end
    isempty(plots) && error("No matching chromatograms found")
    nrows = ceil(Int, length(plots) / ncols)
    return plot(plots...; layout=(nrows, ncols), size=(ncols * 550, nrows * 380))
end

"""
    plot_top_targets(data; n=12, ncols=3)

Plot the `n` highest-prob target chromatograms in a grid.
"""
function plot_top_targets(data; n::Int=12, ncols::Int=3)
    psm = data.psm_best
    targets = sort(psm[psm.target .== true, :], :prob; rev=true)
    n = min(n, nrow(targets))
    pairs = [(UInt32(targets.precursor_idx[i]), targets.source_file[i]) for i in 1:n]
    return _cp_plot_grid(data, pairs; ncols)
end

"""
    plot_high_scoring_decoys(data; n=12, ncols=3)

Plot the `n` highest-prob decoy chromatograms in a grid.
"""
function plot_high_scoring_decoys(data; n::Int=12, ncols::Int=3)
    psm = data.psm_best
    decoys = sort(psm[psm.target .== false, :], :prob; rev=true)
    n = min(n, nrow(decoys))
    pairs = [(UInt32(decoys.precursor_idx[i]), decoys.source_file[i]) for i in 1:n]
    return _cp_plot_grid(data, pairs; ncols)
end

"""
    plot_random_pairs(data; n=6, ncols=2)

Plot `n` target + decoy pairs side by side, sampled across the prob range
and matched by similar prob score.
"""
function plot_random_pairs(data; n::Int=6, ncols::Int=2)
    psm = data.psm_best
    targets = sort(psm[psm.target .== true, :], :prob; rev=true)
    decoys  = sort(psm[psm.target .== false, :], :prob; rev=true)

    n_t = min(n, nrow(targets))
    step = max(1, nrow(targets) ÷ n_t)
    tidxs = [1 + (i-1)*step for i in 1:n_t]
    tidxs = tidxs[tidxs .<= nrow(targets)]

    pairs = Tuple{UInt32, String}[]
    for ti in tidxs
        push!(pairs, (UInt32(targets.precursor_idx[ti]), targets.source_file[ti]))
        diffs = abs.(decoys.prob .- targets.prob[ti])
        di = argmin(diffs)
        push!(pairs, (UInt32(decoys.precursor_idx[di]), decoys.source_file[di]))
    end
    return _cp_plot_grid(data, pairs; ncols)
end

# ============================================================
# DIA-NN vs Pioneer Global Scoring Comparisons
# ============================================================

"""
    _cp_build_candidate_list(data, mode; max_global_qval=0.01)

Build a sorted list of (precursor_idx, source_file) candidates.
`mode` is `:diann_not_pioneer` or `:diann_and_pioneer`.
Returns a Vector of (pidx, src) sorted by first-pass prob (highest first).
"""
function _cp_build_candidate_list(data, mode::Symbol;
                                   max_global_qval::Float64=0.01)
    psm = data.psm_best
    candidates = Tuple{UInt32, String, Float64}[]

    for i in 1:nrow(psm)
        psm.target[i] || continue
        pidx = UInt32(psm.precursor_idx[i])
        src  = psm.source_file[i]

        # Must be DIA-NN detected in this file
        haskey(data.pidx_file_diann_rt, (pidx, src)) || continue

        gqv = get(data.pidx_file_global_qval, (pidx, src), nothing)
        pioneer_pass = gqv !== nothing && gqv <= max_global_qval

        if mode === :diann_not_pioneer && !pioneer_pass
            push!(candidates, (pidx, src, Float64(psm.prob[i])))
        elseif mode === :diann_and_pioneer && pioneer_pass
            push!(candidates, (pidx, src, Float64(psm.prob[i])))
        end
    end

    sort!(candidates; by=x -> x[3], rev=true)
    return [(c[1], c[2]) for c in candidates]
end

"""
    diann_not_pioneer(data; max_global_qval=0.01) -> Vector

Build and return a sorted list of (precursor_idx, source_file) for precursors
DIA-NN identified but Pioneer's global scoring REJECTED.
Sorted by first-pass prob descending.

Then browse with:
    list = diann_not_pioneer(cp_data)
    plot_chrom_pair(cp_data, list, 1)   # 1st entry
    plot_chrom_pair(cp_data, list, 2)   # 2nd entry, etc.
"""
function diann_not_pioneer(data; max_global_qval::Float64=0.01)
    list = _cp_build_candidate_list(data, :diann_not_pioneer; max_global_qval)
    println("DIA-NN found / Pioneer rejected: $(length(list)) precursors (sorted by prob)")
    return list
end

"""
    diann_and_pioneer(data; max_global_qval=0.01) -> Vector

Build and return a sorted list of (precursor_idx, source_file) for precursors
both DIA-NN identified and Pioneer's global scoring ACCEPTED.
Sorted by first-pass prob descending.

Then browse with:
    list = diann_and_pioneer(cp_data)
    plot_chrom_pair(cp_data, list, 1)
"""
function diann_and_pioneer(data; max_global_qval::Float64=0.01)
    list = _cp_build_candidate_list(data, :diann_and_pioneer; max_global_qval)
    println("DIA-NN found / Pioneer accepted: $(length(list)) precursors (sorted by prob)")
    return list
end

"""
    plot_chrom_pair(data, list, idx)

Plot a single full-window chromatogram from a precomputed candidate list.
`list` is a Vector of (precursor_idx, source_file) from `diann_not_pioneer()`
or `diann_and_pioneer()`. `idx` selects which entry (1-based).
"""
function plot_chrom_pair(data, list::Vector{Tuple{UInt32, String}}, idx::Int)
    1 <= idx <= length(list) || error("idx must be in 1:$(length(list))")
    pidx, src = list[idx]
    gi = _cp_find_group(data, pidx, src)
    gi === nothing && error("No chromatogram found for precursor_idx=$pidx, file=$src")
    println("[$idx / $(length(list))]  precursor_idx=$pidx  file=$(src)")
    return plot_chrom(data, gi)
end

# ============================================================
# Score Distribution Plots
# ============================================================

"""
    plot_score_distributions(data)

Plot scribe score and Pearson/Spearman eigenvalue distributions for targets,
decoys, and DIA-NN precursors. Returns a 4-row × 2-column panel figure:
  Row 1: Scribe score | (empty/scribe duplicate)
  Row 2: Pearson λ1 fraction | Spearman λ1 fraction
  Row 3: Pearson λ1           | Spearman λ1
  Row 4: Pearson λ2           | Spearman λ2
"""
function plot_score_distributions(data)
    psm = data.psm_best

    # --- Scribe score ---
    target_mask = psm.target .== true
    diann_mask = Bool[get(data.pidx_is_diann, UInt32(p), false) for p in psm.precursor_idx]
    diann_target_mask = target_mask .& diann_mask

    target_scribe = Float64.(psm[target_mask, :scribe])
    decoy_scribe = Float64.(psm[.!target_mask, :scribe])
    diann_scribe = Float64.(psm[diann_target_mask, :scribe])

    p_scribe = histogram(decoy_scribe; bins=100, normalize=:probability,
                         label="Decoys (n=$(length(decoy_scribe)))",
                         alpha=0.4, color=:red, linewidth=0)
    histogram!(p_scribe, target_scribe; bins=100, normalize=:probability,
               label="Targets (n=$(length(target_scribe)))",
               alpha=0.4, color=:steelblue, linewidth=0)
    histogram!(p_scribe, diann_scribe; bins=100, normalize=:probability,
               label="DIA-NN (n=$(length(diann_scribe)))",
               alpha=0.4, color=:green, linewidth=0)
    title!(p_scribe, "Scribe Score (first-pass PSMs)")
    xlabel!(p_scribe, "Scribe")
    ylabel!(p_scribe, "Fraction")

    # --- Compute Pearson and Spearman eigenvalues for all chromatogram groups ---
    println("  Computing Pearson + Spearman eigenvalues for all chromatogram groups...")
    gdf = data.corr_gdf
    n_groups = length(gdf)

    # Storage: [pearson/spearman] × [target/decoy/diann] × [lam_frac/lam1/lam2]
    vecs = Dict{Symbol, Vector{Float64}}()
    for prefix in (:p, :s), cat in (:target, :decoy, :diann), metric in (:lam_frac, :lam1, :lam2)
        k = Symbol(prefix, :_, cat, :_, metric)
        vecs[k] = Float64[]
        sizehint!(vecs[k], n_groups)
    end

    for gi in 1:n_groups
        chrom = gdf[gi]
        nrow(chrom) < 3 && continue
        is_target = first(chrom.target)
        pidx = UInt32(first(chrom.precursor_idx))
        is_diann = is_target && get(data.pidx_is_diann, pidx, false)

        pr = compute_corr_eigenvalues(chrom)
        sr = compute_spearman_eigenvalues(chrom)

        for (prefix, res) in ((:p, pr), (:s, sr))
            evals = res.eigenvalues
            if is_target
                push!(vecs[Symbol(prefix, :_target_lam_frac)], res.lambda1_frac)
                push!(vecs[Symbol(prefix, :_target_lam1)], evals[1])
                length(evals) >= 2 && push!(vecs[Symbol(prefix, :_target_lam2)], evals[2])
                if is_diann
                    push!(vecs[Symbol(prefix, :_diann_lam_frac)], res.lambda1_frac)
                    push!(vecs[Symbol(prefix, :_diann_lam1)], evals[1])
                    length(evals) >= 2 && push!(vecs[Symbol(prefix, :_diann_lam2)], evals[2])
                end
            else
                push!(vecs[Symbol(prefix, :_decoy_lam_frac)], res.lambda1_frac)
                push!(vecs[Symbol(prefix, :_decoy_lam1)], evals[1])
                length(evals) >= 2 && push!(vecs[Symbol(prefix, :_decoy_lam2)], evals[2])
            end
        end
    end
    println("    Computed: $(length(vecs[:p_target_lam_frac])) target, " *
            "$(length(vecs[:p_decoy_lam_frac])) decoy, $(length(vecs[:p_diann_lam_frac])) DIA-NN groups")

    # --- Helper to build a 3-layer histogram panel ---
    function _hist3(decoy_v, target_v, diann_v, title_str, xlabel_str)
        ph = histogram(decoy_v; bins=100, normalize=:probability,
                       label="Decoys (n=$(length(decoy_v)))",
                       alpha=0.4, color=:red, linewidth=0)
        histogram!(ph, target_v; bins=100, normalize=:probability,
                   label="Targets (n=$(length(target_v)))",
                   alpha=0.4, color=:steelblue, linewidth=0)
        histogram!(ph, diann_v; bins=100, normalize=:probability,
                   label="DIA-NN (n=$(length(diann_v)))",
                   alpha=0.4, color=:green, linewidth=0)
        title!(ph, title_str)
        xlabel!(ph, xlabel_str)
        ylabel!(ph, "Fraction")
        return ph
    end

    # Pearson panels
    p_p_frac = _hist3(vecs[:p_decoy_lam_frac], vecs[:p_target_lam_frac], vecs[:p_diann_lam_frac],
                      "Pearson: λ1 Fraction", "λ1 / Trace")
    p_p_lam1 = _hist3(vecs[:p_decoy_lam1], vecs[:p_target_lam1], vecs[:p_diann_lam1],
                      "Pearson: λ1", "λ1")
    p_p_lam2 = _hist3(vecs[:p_decoy_lam2], vecs[:p_target_lam2], vecs[:p_diann_lam2],
                      "Pearson: λ2", "λ2")

    # Spearman panels
    p_s_frac = _hist3(vecs[:s_decoy_lam_frac], vecs[:s_target_lam_frac], vecs[:s_diann_lam_frac],
                      "Spearman: λ1 Fraction", "λ1 / Trace")
    p_s_lam1 = _hist3(vecs[:s_decoy_lam1], vecs[:s_target_lam1], vecs[:s_diann_lam1],
                      "Spearman: λ1", "λ1")
    p_s_lam2 = _hist3(vecs[:s_decoy_lam2], vecs[:s_target_lam2], vecs[:s_diann_lam2],
                      "Spearman: λ2", "λ2")

    # Layout: 4 rows × 2 cols
    # Row 1: Scribe (spans conceptually, but we duplicate or use empty)
    # Rows 2-4: Pearson (left) | Spearman (right)
    p_empty = plot(; framestyle=:none, legend=false, title="")

    p = plot(p_scribe, p_empty,
             p_p_frac, p_s_frac,
             p_p_lam1, p_s_lam1,
             p_p_lam2, p_s_lam2;
             layout=(4, 2), size=(1200, 1500), margin=5Plots.mm)
    return p
end

"""
    plot_score_density(data; corr_type=:both, nbins=80)

2D density heatmaps of λ1 fraction vs scribe score, with separate panels for
targets, decoys, and DIA-NN precursors. `corr_type` can be `:pearson`,
`:spearman`, or `:both` (default: 2 rows).
"""
function plot_score_density(data; corr_type::Symbol=:both, nbins::Int=80)
    gdf = data.corr_gdf
    psm = data.psm_best
    n_groups = length(gdf)

    # Build fast (precursor_idx, source_file) -> scribe lookup
    scribe_lookup = Dict{Tuple{UInt32, String}, Float64}()
    for i in 1:nrow(psm)
        scribe_lookup[(UInt32(psm.precursor_idx[i]), psm.source_file[i])] = Float64(psm.scribe[i])
    end

    # Collect (scribe, pearson_lambda2, spearman_lambda2) per group
    println("  Computing eigenvalues and matching scribe scores...")
    target_scribe = Float64[]; target_p_lam2 = Float64[]; target_s_lam2 = Float64[]
    decoy_scribe  = Float64[]; decoy_p_lam2  = Float64[]; decoy_s_lam2  = Float64[]
    diann_scribe  = Float64[]; diann_p_lam2  = Float64[]; diann_s_lam2  = Float64[]
    for v in (target_scribe, target_p_lam2, target_s_lam2,
              decoy_scribe, decoy_p_lam2, decoy_s_lam2,
              diann_scribe, diann_p_lam2, diann_s_lam2)
        sizehint!(v, n_groups)
    end

    for gi in 1:n_groups
        chrom = gdf[gi]
        nrow(chrom) < 3 && continue
        pidx = UInt32(first(chrom.precursor_idx))
        src  = first(chrom.source_file)
        sc = get(scribe_lookup, (pidx, src), nothing)
        sc === nothing && continue

        is_target = first(chrom.target)
        pr = compute_corr_eigenvalues(chrom)
        sr = compute_spearman_eigenvalues(chrom)

        p_lam2 = length(pr.eigenvalues) >= 2 ? max(pr.eigenvalues[2], 0.0) : 0.0
        s_lam2 = length(sr.eigenvalues) >= 2 ? max(sr.eigenvalues[2], 0.0) : 0.0

        if is_target
            push!(target_scribe, sc)
            push!(target_p_lam2, p_lam2)
            push!(target_s_lam2, s_lam2)
            if get(data.pidx_is_diann, pidx, false)
                push!(diann_scribe, sc)
                push!(diann_p_lam2, p_lam2)
                push!(diann_s_lam2, s_lam2)
            end
        else
            push!(decoy_scribe, sc)
            push!(decoy_p_lam2, p_lam2)
            push!(decoy_s_lam2, s_lam2)
        end
    end
    println("    Matched: $(length(target_scribe)) target, $(length(decoy_scribe)) decoy, $(length(diann_scribe)) DIA-NN")

    # --- Build panels using manual 2D binning with log color scale ---
    all_scribe = vcat(target_scribe, decoy_scribe)
    scribe_lo = isempty(all_scribe) ? 0.0 : floor(minimum(all_scribe))
    scribe_hi = isempty(all_scribe) ? 15.0 : ceil(maximum(all_scribe))

    all_lam2 = vcat(target_p_lam2, decoy_p_lam2, target_s_lam2, decoy_s_lam2)
    lam2_hi = isempty(all_lam2) ? 1.0 : quantile(all_lam2, 0.99)  # clip top 1% outliers
    lam2_hi = max(lam2_hi, 0.01)  # avoid degenerate range

    function _density_panel(x, y, title_str; nbins=nbins)
        isempty(x) && return plot(; title=title_str, framestyle=:none)
        # Manual 2D histogram
        counts = zeros(Float64, nbins, nbins)
        dx = (scribe_hi - scribe_lo) / nbins
        dy = lam2_hi / nbins
        for k in eachindex(x)
            xi = clamp(floor(Int, (x[k] - scribe_lo) / dx) + 1, 1, nbins)
            yi = clamp(floor(Int, y[k] / dy) + 1, 1, nbins)
            counts[yi, xi] += 1.0
        end
        # Log scale: log10(count + 1) so empty bins → 0
        log_counts = log10.(counts .+ 1.0)
        x_centers = collect(range(scribe_lo + dx/2, scribe_hi - dx/2; length=nbins))
        y_centers = collect(range(dy/2, lam2_hi - dy/2; length=nbins))
        heatmap(x_centers, y_centers, log_counts;
                color=:viridis, colorbar=true, colorbar_title="log10(n+1)",
                title=title_str, titlefontsize=9,
                xlabel="Scribe", ylabel="λ2 (raw)",
                xlims=(scribe_lo, scribe_hi), ylims=(0, lam2_hi))
    end

    panels = Plots.Plot[]

    if corr_type in (:pearson, :both)
        push!(panels, _density_panel(target_scribe, target_p_lam2,
              "Pearson λ2: Targets (n=$(length(target_scribe)))"))
        push!(panels, _density_panel(decoy_scribe, decoy_p_lam2,
              "Pearson λ2: Decoys (n=$(length(decoy_scribe)))"))
        push!(panels, _density_panel(diann_scribe, diann_p_lam2,
              "Pearson λ2: DIA-NN (n=$(length(diann_scribe)))"))
    end

    if corr_type in (:spearman, :both)
        push!(panels, _density_panel(target_scribe, target_s_lam2,
              "Spearman λ2: Targets (n=$(length(target_scribe)))"))
        push!(panels, _density_panel(decoy_scribe, decoy_s_lam2,
              "Spearman λ2: Decoys (n=$(length(decoy_scribe)))"))
        push!(panels, _density_panel(diann_scribe, diann_s_lam2,
              "Spearman λ2: DIA-NN (n=$(length(diann_scribe)))"))
    end

    nrows = corr_type == :both ? 2 : 1
    p = plot(panels...; layout=(nrows, 3), size=(1500, nrows * 450), margin=5Plots.mm)
    return p
end

# ============================================================
# B2 Precursor Browsing (DIA-NN YES, first-pass filtered)
# ============================================================

"""
    b2_precursors(data) -> Vector{Tuple{UInt32, String}}

Build a sorted list of B2 precursors (DIA-NN detected, had first-pass PSM,
filtered at first-pass hybrid PEP/q-value stage). Sorted by first-pass
global PEP ascending (best near-misses first).

Browse with:
    b2 = b2_precursors(cp_data)
    plot_chrom_pair(cp_data, b2, 1)
    plot_chrom_pair(cp_data, b2, 2)
"""
function b2_precursors(data)
    isempty(data.b2_pidxs) && error("No B2 precursors — check that first_pass_global_scores dump exists")

    psm = data.psm_best
    # Collect (pidx, src, fp_pep) for B2 precursors with chromatogram data
    candidates = Tuple{UInt32, String, Float32}[]
    for i in 1:nrow(psm)
        pidx = UInt32(psm.precursor_idx[i])
        pidx in data.b2_pidxs || continue
        src = psm.source_file[i]
        fp_pep = get(data.fp_pep_lookup, pidx, 1.0f0)
        push!(candidates, (pidx, src, fp_pep))
    end

    # Sort by first-pass PEP ascending (best near-misses first)
    sort!(candidates; by=x -> x[3])

    # Keep best file per precursor (lowest PEP already sorted first)
    seen = Set{UInt32}()
    result = Tuple{UInt32, String}[]
    for (pidx, src, pep) in candidates
        pidx in seen && continue
        push!(seen, pidx)
        push!(result, (pidx, src))
    end

    println("B2 precursors with chromatograms: $(length(result)) (sorted by first-pass PEP)")
    if !isempty(result)
        best_pep = get(data.fp_pep_lookup, result[1][1], NaN32)
        worst_pep = get(data.fp_pep_lookup, result[end][1], NaN32)
        println("  PEP range: $(round(Float64(best_pep); digits=4)) to $(round(Float64(worst_pep); digits=4))")
    end
    return result
end

"""
    b2_grid(data; n=12, ncols=3)

Plot a grid of B2 precursor chromatograms, sorted by first-pass PEP (best first).
"""
function b2_grid(data; n::Int=12, ncols::Int=3)
    b2 = b2_precursors(data)
    n = min(n, length(b2))
    return _cp_plot_grid(data, b2[1:n]; ncols, show_points=false)
end

# ============================================================
# Main: load data on include
# ============================================================

println("\n+========================================================+")
println("|  FRAGCORR Chromatogram Viewer                          |")
println("+========================================================+")

cp_data = cp_load_all_data()

println("\nUsage:")
println("  plot_chrom(cp_data, 50)                       # Plot group #50")
println("  plot_chrom_by_pidx(cp_data, 12345)            # Plot by precursor_idx")
println("  plot_top_targets(cp_data; n=12)                # Top-prob targets grid")
println("  plot_high_scoring_decoys(cp_data; n=12)        # Top-prob decoys grid")
println("  plot_random_pairs(cp_data; n=6)                # Target/decoy pairs")
println()
println("  # B2 precursors (DIA-NN YES, first-pass filtered):")
println("  b2 = b2_precursors(cp_data)                   # Build sorted list")
println("  plot_chrom_pair(cp_data, b2, 1)               # Browse one at a time")
println("  b2_grid(cp_data; n=12)                        # Grid view")
println()
println("  # DIA-NN vs Pioneer global scoring:")
println("  missed = diann_not_pioneer(cp_data)            # Build list: DIA-NN YES, Pioneer NO")
println("  found  = diann_and_pioneer(cp_data)            # Build list: DIA-NN YES, Pioneer YES")
println("  plot_chrom_pair(cp_data, missed, 1)            # Browse one at a time")
println("  plot_chrom_pair(cp_data, found, 1)")
