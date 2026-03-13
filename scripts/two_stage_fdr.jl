#!/usr/bin/env julia
#=
Two-Stage FDR: Experiment-Wide Q-Values on Globally Passing Precursors

Stage 1: Replicate current pipeline (PEP-calibrated log-odds → global q-values)
Stage 2: Filter per-file rows to Stage 1 passing precursors, concatenate,
         compute experiment-wide q-values, compare per-file counts.

Usage:
  julia --threads 4 --project=. scripts/two_stage_fdr.jl <results_dir>
  julia --threads 4 --project=. scripts/two_stage_fdr.jl /path/to/exp_depth3_leaves10
=#

using Pioneer
using Arrow, DataFrames, Tables
using Printf

const GLOBAL_Q_THRESHOLD = 0.01f0
const STAGE2_Q_THRESHOLD = 0.01f0

# ── Log-odds combination (reimplemented from prescore_aggregation.jl) ──────
function logodds_combine(probs::Vector{Float32}, top_n::Int, floor::Float32)::Float32
    isempty(probs) && return 0.0f0
    eps = 1f-6
    n = min(length(probs), top_n)
    if n == 1
        p = length(probs) == 1 ? probs[1] : maximum(probs)
        return clamp(p, floor, 1 - eps)
    end
    partialsort!(probs, 1:n; rev=true)
    lo_sum = 0.0f0
    @inbounds for i in 1:n
        c = clamp(probs[i], floor, 1 - eps)
        lo_sum += log(c / (1 - c))
    end
    avg = lo_sum / n
    return 1.0f0 / (1 + exp(-avg))
end

# ── Main ───────────────────────────────────────────────────────────────────
function main(results_dir::String)
    prescore_dir = joinpath(results_dir, "temp_data", "prescore_scores")
    if !isdir(prescore_dir)
        error("Prescore directory not found: $prescore_dir")
    end

    arrow_files = filter(f -> endswith(f, ".arrow"), readdir(prescore_dir))
    if isempty(arrow_files)
        error("No .arrow files in $prescore_dir")
    end
    sort!(arrow_files)

    n_files = length(arrow_files)
    sqrt_n = max(1, floor(Int, sqrt(n_files)))
    println("Found $n_files files, top_n = $sqrt_n")
    for f in arrow_files
        println("  $f")
    end

    # ── Load all per-file data ─────────────────────────────────────────────
    file_dfs = DataFrame[]
    for fname in arrow_files
        tbl = Arrow.Table(joinpath(prescore_dir, fname))
        df = DataFrame(
            precursor_idx = UInt32.(tbl[:precursor_idx]),
            lgbm_prob     = Float32.(tbl[:lgbm_prob]),
            target        = Bool.(tbl[:target]),
        )
        df.file_name .= replace(fname, ".arrow" => "")
        push!(file_dfs, df)
    end

    # ── Stage 1: PEP calibration + global aggregation ─────────────────────
    println("\n═══ STAGE 1: Global FDR (replicating pipeline) ═══")

    # PEP-calibrate each file
    for df in file_dfs
        n_t = count(df.target)
        n_d = count(!, df.target)
        if n_t >= 50 && n_d >= 20
            peps = Vector{Float32}(undef, nrow(df))
            Pioneer.get_PEP!(df.lgbm_prob, df.target, peps; doSort=true)
            df.calibrated_prob = 1.0f0 .- peps
        else
            df.calibrated_prob = copy(df.lgbm_prob)
        end
    end

    # Collect calibrated probs per precursor across files
    prec_probs = Dict{UInt32, Vector{Float32}}()
    prec_target = Dict{UInt32, Bool}()
    for df in file_dfs
        for row in eachrow(df)
            pid = row.precursor_idx
            p = row.calibrated_prob
            if haskey(prec_probs, pid)
                push!(prec_probs[pid], p)
            else
                prec_probs[pid] = Float32[p]
                prec_target[pid] = row.target
            end
        end
    end

    # Combine via log-odds of top-√n
    n_unique = length(prec_probs)
    global_pids = Vector{UInt32}(undef, n_unique)
    global_probs = Vector{Float32}(undef, n_unique)
    global_targets = Vector{Bool}(undef, n_unique)

    for (i, (pid, probs)) in enumerate(prec_probs)
        global_pids[i] = pid
        global_probs[i] = logodds_combine(copy(probs), sqrt_n, 0.01f0)
        global_targets[i] = prec_target[pid]
    end

    # Global q-values
    global_qvals = Vector{Float32}(undef, n_unique)
    Pioneer.get_qvalues!(global_probs, global_targets, global_qvals; doSort=true)

    # Count targets at various thresholds
    for qt in [0.01f0, 0.05f0, 0.10f0]
        n = count(i -> global_qvals[i] <= qt && global_targets[i], eachindex(global_qvals))
        @printf("  Global targets at q ≤ %.2f: %d\n", qt, n)
    end

    # Build passing set at GLOBAL_Q_THRESHOLD
    passing_precs = Set{UInt32}()
    for i in eachindex(global_pids)
        if global_qvals[i] <= GLOBAL_Q_THRESHOLD
            push!(passing_precs, global_pids[i])
        end
    end
    n_passing_targets = count(pid -> prec_target[pid], passing_precs)
    println("  Passing precursors at q ≤ $GLOBAL_Q_THRESHOLD: $(length(passing_precs)) ($(n_passing_targets) targets)")

    # ── Stage 2: Experiment-wide FDR on passing precursors ────────────────
    println("\n═══ STAGE 2: Experiment-wide FDR on globally passing precursors ═══")

    # Filter per-file data to passing precursors, concatenate
    filtered_dfs = DataFrame[]
    for df in file_dfs
        mask = [pid in passing_precs for pid in df.precursor_idx]
        push!(filtered_dfs, df[mask, :])
    end
    concat = vcat(filtered_dfs...)
    println("  Concatenated rows (all files): $(nrow(concat))")
    println("  Targets: $(count(concat.target)), Decoys: $(count(!, concat.target))")

    # Try both raw lgbm_prob and calibrated_prob as scoring columns
    for (score_name, score_col) in [("lgbm_prob", concat.lgbm_prob),
                                     ("calibrated_prob", concat.calibrated_prob)]
        println("\n  ── Stage 2 scores: $score_name ──")

        qvals = Vector{Float32}(undef, nrow(concat))
        Pioneer.get_qvalues!(score_col, concat.target, qvals; doSort=true)

        # Per-file counts at Stage 2 threshold
        concat_with_q = DataFrame(
            file_name = concat.file_name,
            precursor_idx = concat.precursor_idx,
            target = concat.target,
            q_value = qvals,
        )

        println(@sprintf("  %-20s %8s %8s %8s", "File", "Total", "Pass", "Targets"))
        println("  " * "─"^48)

        total_pass = 0
        total_targets = 0
        for gdf in groupby(concat_with_q, :file_name)
            fname = first(gdf.file_name)
            passing_mask = gdf.q_value .<= STAGE2_Q_THRESHOLD
            n_pass = count(passing_mask)
            n_targets_pass = count(passing_mask .& gdf.target)
            total_pass += n_pass
            total_targets += n_targets_pass
            println(@sprintf("  %-20s %8d %8d %8d", fname, nrow(gdf), n_pass, n_targets_pass))
        end
        println("  " * "─"^48)
        println(@sprintf("  %-20s %8s %8d %8d", "TOTAL", "", total_pass, total_targets))
    end

    # ── Comparison with current pipeline ──────────────────────────────────
    println("\n═══ COMPARISON SUMMARY ═══")
    println("  Stage 1 global threshold: q ≤ $GLOBAL_Q_THRESHOLD")
    println("  Stage 2 experiment-wide threshold: q ≤ $STAGE2_Q_THRESHOLD")
    println("  Stage 1 passing target precursors (unique): $n_passing_targets")

    # Current pipeline: per-file counts of globally passing targets
    println("\n  Current pipeline (global pass → per-file count of observations):")
    firstpass_per_file = Dict{String, Int}()
    for df in file_dfs
        fname = first(df.file_name)
        mask = [pid in passing_precs for pid in df.precursor_idx]
        n_obs = count(mask .& df.target)
        firstpass_per_file[fname] = n_obs
        println(@sprintf("    %-20s %d target observations", fname, n_obs))
    end

    # ── Path B: Full pipeline (ScoringSearch LightGBM) ─────────────────────
    passing_psms_dir = joinpath(results_dir, "temp_data", "passing_psms")
    if isdir(passing_psms_dir)
        println("\n═══ PATH B: Full Pipeline (ScoringSearch LightGBM, 20 features) ═══")

        psm_files = filter(f -> endswith(f, ".arrow"), readdir(passing_psms_dir))
        sort!(psm_files)

        scoring_per_file = Dict{String, Int}()
        scoring_unique_targets = Set{UInt32}()

        for fname in psm_files
            tbl = Arrow.Table(joinpath(passing_psms_dir, fname))
            pids = UInt32.(tbl[:precursor_idx])
            targets = Bool.(tbl[:target])
            qvals = Float32.(tbl[:qval])

            file_label = replace(fname, ".arrow" => "")
            n_targets = count(i -> qvals[i] <= 0.01f0 && targets[i], eachindex(qvals))
            scoring_per_file[file_label] = n_targets

            for i in eachindex(pids)
                if qvals[i] <= 0.01f0 && targets[i]
                    push!(scoring_unique_targets, pids[i])
                end
            end
        end

        println("  Unique target precursors at qval ≤ 0.01: $(length(scoring_unique_targets))")

        # ── Side-by-side comparison ────────────────────────────────────────
        println("\n═══ SIDE-BY-SIDE: FirstPass-Only vs Full Pipeline ═══")
        println(@sprintf("  %-20s %12s %12s %8s", "File", "FirstPass", "FullPipeline", "Delta"))
        println("  " * "─"^56)

        all_files = sort(collect(union(keys(firstpass_per_file), keys(scoring_per_file))))
        total_fp = 0
        total_sp = 0
        for fname in all_files
            fp = get(firstpass_per_file, fname, 0)
            sp = get(scoring_per_file, fname, 0)
            delta = sp - fp
            total_fp += fp
            total_sp += sp
            println(@sprintf("  %-20s %12d %12d %+8d", fname, fp, sp, delta))
        end
        println("  " * "─"^56)
        println(@sprintf("  %-20s %12d %12d %+8d", "TOTAL", total_fp, total_sp, total_sp - total_fp))
        println("\n  FirstPass-only unique targets (global q ≤ 0.01): $n_passing_targets")
        println("  Full pipeline unique targets (per-file qval ≤ 0.01): $(length(scoring_unique_targets))")
        println(@sprintf("  Unique target delta: %+d", length(scoring_unique_targets) - n_passing_targets))
    else
        println("\n  [Skipping Path B: $passing_psms_dir not found]")
    end
end

# ── Entry point ───────────────────────────────────────────────────────────
if isempty(ARGS)
    error("Usage: julia --project=. scripts/two_stage_fdr.jl <results_dir>")
end
main(ARGS[1])
