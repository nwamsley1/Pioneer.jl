#=
Borderline Precursor Analysis (4-Category)
==========================================
Classifies DIA-NN-passing precursors into 4 categories based on WHERE Pioneer
loses them in the pipeline:

  A)  Pioneer passes final scoring (q ≤ 0.01) — agreement with DIA-NN
  B1) Has first-pass PSM, passes first-pass filter, but fails final scoring q-value
      — "recoverable at scoring stage"
  B2) Has first-pass PSM, but filtered at first-pass global PEP stage
      — "recoverable at first-pass filter"
  C)  Never matched in first pass at all — not in library or spectral match failed

Requires two diagnostic dumps:
  - temp_data/first_pass_global_scores/global_scores.arrow  (first-pass, before hybrid filter)
  - temp_data/pre_filter_global_scores/global_scores.arrow   (final scoring, before q=0.01 filter)

Usage:
    include("scripts/analyze_borderline_precursors.jl")

    # Summary statistics
    borderline_summary(bl_data)

    # Browse B1 precursors (passed first-pass, failed final scoring)
    plot_borderline(bl_data, 1)
    plot_near_miss(bl_data, 1)    # q ∈ (0.01, 0.05] — almost passed
    plot_hopeless(bl_data, 1)     # q > 0.25 — far from passing

    # Grid view
    plot_borderline_grid(bl_data; n=6)
=#

using Arrow, DataFrames, Tables, Parquet2
using Statistics
using Plots

# ============================================================
# Configuration — update RESULTS_DIR after running the search
# ============================================================

const BL_RESULTS_DIR = "/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/OlsenEclipse_diagnostic"
const BL_TEMP_DIR = joinpath(BL_RESULTS_DIR, "temp_data")
const BL_LIBRARY_PATH = "/Users/nathanwamsley/Data/SPEC_LIBS/3P_rank_based.poin"
const BL_DIANN_REPORT = "/Users/nathanwamsley/Data/MS_DATA/OlsenEclipse/DIA-NN_results/OlsenExplorisThreeProteome500ng-11-24-2025-report.parquet"
const BL_Q_THRESHOLD = 0.01

# Pioneer files (stem names matching second_pass_psms/*.arrow)
const BL_PIONEER_FILES = ["E45H50Y5_2", "E45H50Y5_3", "E45H50Y5_4",
                          "E5H50Y45_1", "E5H50Y45_2", "E5H50Y45_4"]

# ============================================================
# Canonical key helpers (same as fragcorr script)
# ============================================================

const _BL_CanonicalKey = Tuple{String, Vector{Tuple{Int,Int}}, Int}

function _bl_parse_diann_seq(mod_seq::AbstractString, charge::Int)::_BL_CanonicalKey
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

function _bl_pioneer_key(seq::AbstractString, structural_mods::AbstractString, charge::Integer)::_BL_CanonicalKey
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

function _bl_diann_run_to_pioneer(run::AbstractString)::Union{String, Nothing}
    m = match(r"(E\d+H\d+Y\d+)_\d+SPD_DIA_(\d+)", run)
    m === nothing && return nothing
    return "$(m.captures[1])_$(m.captures[2])"
end

# ============================================================
# Data Loading
# ============================================================

function bl_load_data()
    println("Loading borderline precursor analysis data...")

    # --- 1a. First-pass global scores (before hybrid PEP filter) ---
    fp_scores_path = joinpath(BL_TEMP_DIR, "first_pass_global_scores", "global_scores.arrow")
    has_fp_scores = isfile(fp_scores_path)
    fp_scores = nothing
    fp_lookup = Dict{UInt32, NamedTuple{(:global_prob, :global_pep, :global_qval, :target),
                                         Tuple{Float32, Float32, Float32, Bool}}}()
    if has_fp_scores
        println("  Loading first-pass global scores...")
        fp_scores = DataFrame(Tables.columntable(Arrow.Table(fp_scores_path)))
        println("    $(nrow(fp_scores)) precursors total ($(count(fp_scores.target)) targets)")
        for i in 1:nrow(fp_scores)
            fp_lookup[UInt32(fp_scores.precursor_idx[i])] = (
                global_prob = fp_scores.global_prob[i],
                global_pep = fp_scores.global_pep[i],
                global_qval = fp_scores.global_qval[i],
                target = fp_scores.target[i]
            )
        end
    else
        println("  WARNING: First-pass global scores not found at $fp_scores_path")
        println("           Run SearchDIA with updated code to generate this dump.")
        println("           Falling back to 3-category classification (A/B/C).")
    end

    # --- 1b. Final scoring global scores (before q=0.01 filter) ---
    global_scores_path = joinpath(BL_TEMP_DIR, "pre_filter_global_scores", "global_scores.arrow")
    println("  Loading final scoring global scores...")
    global_scores = DataFrame(Tables.columntable(Arrow.Table(global_scores_path)))
    println("    $(nrow(global_scores)) precursors total")
    println("    $(count(global_scores.target)) targets")
    n_pass = count(r -> r.target && r.global_qval <= BL_Q_THRESHOLD, eachrow(global_scores))
    n_fail = count(r -> r.target && r.global_qval > BL_Q_THRESHOLD, eachrow(global_scores))
    println("    Targets passing global q <= $(BL_Q_THRESHOLD): $n_pass")
    println("    Targets failing global q > $(BL_Q_THRESHOLD): $n_fail")

    # Build fast lookup: precursor_idx -> (global_prob, global_qval, target)
    global_lookup = Dict{UInt32, NamedTuple{(:global_prob, :global_qval, :target),
                                             Tuple{Float32, Float32, Bool}}}()
    for i in 1:nrow(global_scores)
        global_lookup[UInt32(global_scores.precursor_idx[i])] = (
            global_prob = global_scores.global_prob[i],
            global_qval = global_scores.global_qval[i],
            target = global_scores.target[i]
        )
    end

    # --- 2. Load library for sequence/charge/mods ---
    println("  Loading precursor library...")
    lib_df = DataFrame(Arrow.Table(joinpath(BL_LIBRARY_PATH, "precursors_table.arrow")))
    pidx_meta = Dict{UInt32, NamedTuple{(:sequence, :charge, :mods, :is_decoy),
                                         Tuple{String, Int, String, Bool}}}()
    for i in 1:nrow(lib_df)
        pidx = UInt32(i)
        smods = ismissing(lib_df.structural_mods[i]) ? "" : lib_df.structural_mods[i]
        pidx_meta[pidx] = (sequence=lib_df.sequence[i],
                           charge=Int(lib_df.prec_charge[i]),
                           mods=smods,
                           is_decoy=lib_df.is_decoy[i])
    end

    # Build canonical_key -> precursor_idx reverse lookup (targets only)
    key_to_pidxs = Dict{_BL_CanonicalKey, Vector{UInt32}}()
    for (pidx, meta) in pidx_meta
        meta.is_decoy && continue
        k = _bl_pioneer_key(meta.sequence, meta.mods, meta.charge)
        push!(get!(Vector{UInt32}, key_to_pidxs, k), pidx)
    end

    # --- 3. Load DIA-NN report ---
    println("  Loading DIA-NN report...")
    pioneer_file_set = Set(BL_PIONEER_FILES)
    diann_report = DataFrame(Parquet2.Dataset(BL_DIANN_REPORT); copycols=false)

    # Build DIA-NN fair keys and per-file RT lookup
    diann_fair_keys = Set{_BL_CanonicalKey}()
    diann_key_file_rt = Dict{Tuple{_BL_CanonicalKey, String}, Float64}()

    for i in 1:nrow(diann_report)
        pf = _bl_diann_run_to_pioneer(diann_report[i, "Run"])
        pf === nothing && continue
        pf in pioneer_file_set || continue
        qty = diann_report[i, "Precursor.Quantity"]
        (ismissing(qty) || !(qty isa Number) || qty <= 0) && continue
        key = _bl_parse_diann_seq(diann_report[i, "Modified.Sequence"],
                                   Int(diann_report[i, "Precursor.Charge"]))
        push!(diann_fair_keys, key)
        diann_key_file_rt[(key, pf)] = Float64(diann_report[i, "RT"])
    end
    println("    DIA-NN fair keys (quantified in >=1 Pioneer file): $(length(diann_fair_keys))")

    # Map precursor_idx -> is_diann_detected
    pidx_is_diann = Dict{UInt32, Bool}()
    for (pidx, meta) in pidx_meta
        meta.is_decoy && continue
        k = _bl_pioneer_key(meta.sequence, meta.mods, meta.charge)
        pidx_is_diann[pidx] = k in diann_fair_keys
    end

    # Map (precursor_idx, file_stem) -> DIA-NN RT
    pidx_file_diann_rt = Dict{Tuple{UInt32, String}, Float64}()
    for ((key, pf), rt) in diann_key_file_rt
        for pidx in get(key_to_pidxs, key, UInt32[])
            pidx_file_diann_rt[(pidx, pf)] = rt
        end
    end

    # --- 4. Load second_pass_psms (per-file PSM scores before filtering) ---
    println("  Loading second-pass PSMs...")
    second_pass_dir = joinpath(BL_TEMP_DIR, "second_pass_psms")
    all_psms = DataFrame[]
    for f in readdir(second_pass_dir)
        endswith(f, ".arrow") || continue
        stem = replace(f, ".arrow" => "")
        stem in pioneer_file_set || continue
        df = DataFrame(Tables.columntable(Arrow.Table(joinpath(second_pass_dir, f))))
        df[!, :file_stem] .= stem
        push!(all_psms, df)
    end
    psm_df = vcat(all_psms...; cols=:union)
    println("    $(nrow(psm_df)) PSMs across $(length(all_psms)) files")

    # Best PSM per precursor per file
    sort!(psm_df, :prec_prob; rev=true)
    psm_best = combine(groupby(psm_df, [:precursor_idx, :file_stem]), first)

    # --- 5. Classify precursors (4 categories) ---
    println("\n  Classifying precursors...")

    # Categories for DIA-NN-detected target precursors:
    # A)  Pioneer passes final scoring (global_qval <= threshold)
    # B1) Has first-pass PSM, passes first-pass filter, but fails final scoring q-value
    # B2) Has first-pass PSM, but filtered at first-pass global PEP stage (never enters second pass)
    # C)  Never matched in first pass at all

    diann_pidxs = Set{UInt32}()
    for (pidx, is_d) in pidx_is_diann
        is_d && push!(diann_pidxs, pidx)
    end

    # --- Reconstruct the first-pass hybrid filter decision ---
    # Instead of checking `pre_filter_global_scores` (which misses precursors
    # that passed the first-pass filter but got no second-pass PSM), we
    # reconstruct the hybrid filter from the first-pass dump itself.
    fp_precursor_set = Set(keys(fp_lookup))          # seen in first pass

    # Determine which precursors were truly filtered at first pass
    fp_filtered_set = Set{UInt32}()
    if has_fp_scores
        fp_n = nrow(fp_scores)
        fp_order = sortperm(fp_scores.global_pep)

        # Q-value floor: last position where cumulative D/T ≤ 0.15
        cum_t, cum_d, n_qvalue_based = 0, 0, 0
        for (rank, i) in enumerate(fp_order)
            if fp_scores.target[i]
                cum_t += 1
            else
                cum_d += 1
            end
            if cum_t > 0 && cum_d / cum_t <= 0.15
                n_qvalue_based = rank
            end
        end

        n_passing_threshold = count(fp_scores.global_pep .<= 0.5f0)
        n_min_floor = min(50_000, fp_n)
        n_to_keep = max(n_passing_threshold, n_qvalue_based, n_min_floor)

        if n_to_keep < fp_n
            for rank in (n_to_keep + 1):fp_n
                push!(fp_filtered_set, UInt32(fp_scores.precursor_idx[fp_order[rank]]))
            end
        end
        println("    Hybrid filter reconstruction: keeping $n_to_keep / $fp_n")
        println("    Truly filtered at first pass: $(length(fp_filtered_set))")
    end

    cat_a      = UInt32[]  # Pioneer passes final scoring (q ≤ threshold)
    cat_b1     = UInt32[]  # In final scoring, fails q-value
    cat_b2     = UInt32[]  # Truly filtered at first-pass hybrid filter
    cat_b_lost = UInt32[]  # Passed first-pass filter but lost at second pass (no second-pass PSM)
    cat_c      = UInt32[]  # Never matched in first pass

    for pidx in diann_pidxs
        in_fp = pidx in fp_precursor_set
        in_final = haskey(global_lookup, pidx)

        if in_final && global_lookup[pidx].global_qval <= BL_Q_THRESHOLD
            push!(cat_a, pidx)
        elseif in_final
            # In final scoring but fails q-value — passed first-pass filter
            push!(cat_b1, pidx)
        elseif has_fp_scores && in_fp && pidx in fp_filtered_set
            # In first-pass dump AND truly filtered by hybrid filter
            push!(cat_b2, pidx)
        elseif has_fp_scores && in_fp
            # In first-pass dump, passed hybrid filter, but NOT in final scoring
            # — lost because second pass (different params/tolerances) didn't produce a PSM
            push!(cat_b_lost, pidx)
        else
            # Not seen at all (or first-pass dump unavailable and not in final)
            push!(cat_c, pidx)
        end
    end

    println("    DIA-NN detected precursors: $(length(diann_pidxs))")
    println("    A)      Pioneer PASSES final scoring (q <= $(BL_Q_THRESHOLD)): $(length(cat_a))")
    println("    B1)     Passes first-pass filter, FAILS final scoring q: $(length(cat_b1))")
    if has_fp_scores
        println("    B2)     Truly filtered at first-pass hybrid PEP: $(length(cat_b2))")
        println("    B_lost) Passed first-pass filter, lost at second pass: $(length(cat_b_lost))")
    end
    println("    C)      Never matched in first pass: $(length(cat_c))")
    println()
    println("    Sanity check: A + B1 + B2 + B_lost + C = $(length(cat_a) + length(cat_b1) + length(cat_b2) + length(cat_b_lost) + length(cat_c))")

    # Show first-pass PEP distribution for B2 precursors
    if has_fp_scores && !isempty(cat_b2)
        b2_peps = Float32[fp_lookup[pidx].global_pep for pidx in cat_b2]
        b2_probs = Float32[fp_lookup[pidx].global_prob for pidx in cat_b2]
        println("\n    B2 first-pass global PEP distribution:")
        for thresh in (0.1f0, 0.25f0, 0.5f0, 0.75f0, 0.9f0)
            n = count(<=(thresh), b2_peps)
            println("      PEP <= $thresh: $n / $(length(cat_b2))")
        end
        println("    B2 first-pass global prob: median=$(round(median(b2_probs); digits=4)), " *
                "mean=$(round(mean(b2_probs); digits=4))")
    end

    # Show first-pass PEP distribution for B_lost precursors
    if has_fp_scores && !isempty(cat_b_lost)
        bl_peps = Float32[fp_lookup[pidx].global_pep for pidx in cat_b_lost]
        bl_probs = Float32[fp_lookup[pidx].global_prob for pidx in cat_b_lost]
        println("\n    B_lost first-pass global PEP distribution:")
        for thresh in (0.1f0, 0.25f0, 0.5f0, 0.75f0, 0.9f0)
            n = count(<=(thresh), bl_peps)
            println("      PEP <= $thresh: $n / $(length(cat_b_lost))")
        end
        println("    B_lost first-pass global prob: median=$(round(median(bl_probs); digits=4)), " *
                "mean=$(round(mean(bl_probs); digits=4))")
    end

    # --- 6. Build borderline detail table (B1 precursors — in final scoring, fail q) ---
    borderline_set = Set(cat_b1)
    borderline_rows = filter(r -> UInt32(r.precursor_idx) in borderline_set, psm_best)

    # Add global scores and DIA-NN info
    borderline_rows[!, :global_prob] = Float32[
        get(global_lookup, UInt32(p), (global_prob=0f0, global_qval=1f0, target=false)).global_prob
        for p in borderline_rows.precursor_idx]
    borderline_rows[!, :global_qval] = Float32[
        get(global_lookup, UInt32(p), (global_prob=0f0, global_qval=1f0, target=false)).global_qval
        for p in borderline_rows.precursor_idx]
    borderline_rows[!, :diann_rt] = Union{Float64, Missing}[
        get(pidx_file_diann_rt, (UInt32(r.precursor_idx), r.file_stem), missing)
        for r in eachrow(borderline_rows)]
    borderline_rows[!, :sequence] = String[
        get(pidx_meta, UInt32(p), (sequence="?", charge=0, mods="", is_decoy=false)).sequence
        for p in borderline_rows.precursor_idx]
    borderline_rows[!, :charge] = Int[
        get(pidx_meta, UInt32(p), (sequence="?", charge=0, mods="", is_decoy=false)).charge
        for p in borderline_rows.precursor_idx]

    # Sort by global_qval ascending (near-misses first)
    sort!(borderline_rows, :global_qval)

    println("\n  B1 borderline PSM entries (precursor x file): $(nrow(borderline_rows))")

    # --- Sub-categories of B1 borderline ---
    near_miss = filter(r -> r.global_qval <= 0.05, borderline_rows)
    moderate  = filter(r -> 0.05 < r.global_qval <= 0.25, borderline_rows)
    hopeless  = filter(r -> r.global_qval > 0.25, borderline_rows)

    println("    Near-miss (q in ($(BL_Q_THRESHOLD), 0.05]): $(nrow(near_miss)) entries ($(length(unique(near_miss.precursor_idx))) unique)")
    println("    Moderate  (q in (0.05, 0.25]):  $(nrow(moderate)) entries ($(length(unique(moderate.precursor_idx))) unique)")
    println("    Hopeless  (q > 0.25):           $(nrow(hopeless)) entries ($(length(unique(hopeless.precursor_idx))) unique)")

    println("\nReady!")
    return (global_lookup=global_lookup, fp_lookup=fp_lookup,
            pidx_meta=pidx_meta,
            pidx_is_diann=pidx_is_diann, pidx_file_diann_rt=pidx_file_diann_rt,
            psm_best=psm_best, borderline=borderline_rows,
            near_miss=near_miss, moderate=moderate, hopeless=hopeless,
            cat_a=cat_a, cat_b1=cat_b1, cat_b2=cat_b2, cat_b_lost=cat_b_lost, cat_c=cat_c,
            has_fp_scores=has_fp_scores,
            fp_scores=fp_scores, global_scores=global_scores)
end

# ============================================================
# Summary Statistics
# ============================================================

function borderline_summary(data)
    bl = data.borderline

    println("\n", "="^80)
    println("BORDERLINE PRECURSOR SUMMARY (5-Category)")
    println("="^80)

    # Overall category counts
    println("\nCategory breakdown (DIA-NN-detected precursors):")
    println("  A)      Pioneer passes final scoring:              $(length(data.cat_a))")
    println("  B1)     Passes first-pass filter, fails final q:   $(length(data.cat_b1))")
    if data.has_fp_scores
        println("  B2)     Truly filtered at first-pass hybrid PEP:   $(length(data.cat_b2))")
        println("  B_lost) Passed first-pass filter, lost at 2nd pass: $(length(data.cat_b_lost))")
    end
    println("  C)      Never matched in first pass:               $(length(data.cat_c))")
    total = length(data.cat_a) + length(data.cat_b1) + length(data.cat_b2) + length(data.cat_b_lost) + length(data.cat_c)
    println("  Total: $total")

    # B2 details (if available)
    if data.has_fp_scores && !isempty(data.cat_b2)
        println("\nB2 first-pass global PEP distribution (truly filtered at first pass):")
        b2_peps = Float32[data.fp_lookup[pidx].global_pep for pidx in data.cat_b2]
        for thresh in (0.1f0, 0.25f0, 0.5f0, 0.75f0, 0.9f0, 0.95f0)
            n = count(<=(thresh), b2_peps)
            pct = round(100.0 * n / length(b2_peps); digits=1)
            println("  PEP <= $thresh: $(rpad(n, 8)) ($pct%)")
        end
    end

    # B_lost details (if available)
    if data.has_fp_scores && !isempty(data.cat_b_lost)
        println("\nB_lost first-pass global PEP distribution (passed filter, lost at 2nd pass):")
        bl_peps = Float32[data.fp_lookup[pidx].global_pep for pidx in data.cat_b_lost]
        for thresh in (0.1f0, 0.25f0, 0.5f0, 0.75f0, 0.9f0, 0.95f0)
            n = count(<=(thresh), bl_peps)
            pct = round(100.0 * n / length(bl_peps); digits=1)
            println("  PEP <= $thresh: $(rpad(n, 8)) ($pct%)")
        end
    end

    # B1 borderline details
    println("\nB1 global q-value distribution:")
    edges = [0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 1.0]
    for i in 1:length(edges)-1
        lo, hi = edges[i], edges[i+1]
        n = count(r -> lo < r.global_qval <= hi, eachrow(bl))
        n_unique = length(unique(filter(r -> lo < r.global_qval <= hi, bl).precursor_idx))
        println("  q in ($(lpad(lo, 4)), $(lpad(hi, 4))]: $(rpad(n, 8)) entries  $(rpad(n_unique, 6)) unique precursors")
    end

    println("\nB1 global probability distribution:")
    prob_edges = [0.0, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.0]
    for i in 1:length(prob_edges)-1
        lo, hi = prob_edges[i], prob_edges[i+1]
        n = count(r -> lo <= r.global_prob < hi, eachrow(bl))
        println("  prob in [$(lpad(lo, 4)), $(lpad(hi, 4))): $n entries")
    end

    println("\nB1 per-file prec_prob distribution:")
    pp_edges = [0.0, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.0]
    for i in 1:length(pp_edges)-1
        lo, hi = pp_edges[i], pp_edges[i+1]
        n = count(r -> lo <= r.prec_prob < hi, eachrow(bl))
        println("  prec_prob in [$(lpad(lo, 4)), $(lpad(hi, 4))): $n entries")
    end

    println("\nFiles per B1 borderline precursor:")
    pidx_counts = combine(groupby(bl, :precursor_idx), nrow => :n_files)
    for nf in sort(unique(pidx_counts.n_files))
        n = count(==(nf), pidx_counts.n_files)
        println("  $nf file(s): $n precursors")
    end
end

# ============================================================
# Plotting
# ============================================================

"""
    _bl_plot_one(data, row)

Plot a borderline precursor: left panel shows per-file prec_prob vs DIA-NN RT context,
right panel shows the global scoring position.
"""
function _bl_plot_one(data, row)
    pidx = UInt32(row.precursor_idx)
    meta = get(data.pidx_meta, pidx, nothing)
    seq = meta !== nothing ? meta.sequence : "?"
    charge = meta !== nothing ? meta.charge : 0
    seq_display = length(seq) > 30 ? seq[1:27] * "..." : seq

    gs = data.global_lookup[pidx]

    # Get all PSM entries for this precursor across files
    pidx_psms = filter(r -> UInt32(r.precursor_idx) == pidx, data.psm_best)

    # --- Left panel: per-file scores ---
    title_left = "$seq_display $(charge)+\n" *
                 "global_prob=$(round(gs.global_prob; digits=4))  " *
                 "global_qval=$(round(gs.global_qval; digits=4))"

    p_left = plot(; title=title_left, titlefontsize=9,
                  xlabel="File", ylabel="Score", legend=:bottomright, legendfontsize=7)

    files = sort(unique(pidx_psms.file_stem))
    x_pos = 1:length(files)
    prec_probs = Float64[pidx_psms[pidx_psms.file_stem .== f, :prec_prob][1] for f in files]

    bar!(p_left, x_pos, prec_probs; label="prec_prob", color=:steelblue,
         bar_width=0.6, alpha=0.8)

    # Mark which files have DIA-NN detection
    for (xi, f) in enumerate(files)
        diann_rt = get(data.pidx_file_diann_rt, (pidx, f), nothing)
        if diann_rt !== nothing
            scatter!(p_left, [xi], [prec_probs[xi]]; label=(xi == 1 ? "DIA-NN detected" : nothing),
                     color=:green, markersize=8, markershape=:diamond, markerstrokewidth=2)
        end
    end

    hline!(p_left, [0.5]; color=:red, ls=:dash, label="prob=0.5", lw=1)
    xticks!(p_left, x_pos, files; rotation=15, tickfontsize=6)

    # --- Right panel: global score context ---
    # Show where this precursor sits in the global_prob distribution
    target_scores = data.global_scores
    tgt = filter(r -> r.target, target_scores)

    title_right = "Global Score Distribution\n(this precursor marked in red)"
    p_right = histogram(tgt.global_prob; bins=50, label="all targets",
                        color=:lightblue, alpha=0.7, normalize=:probability,
                        title=title_right, titlefontsize=9,
                        xlabel="global_prob", ylabel="fraction")
    vline!(p_right, [gs.global_prob]; color=:red, lw=3, label="this precursor")
    vline!(p_right, [0.5]; color=:orange, lw=1, ls=:dash, label="prob=0.5")

    p = plot(p_left, p_right; layout=(1, 2), size=(1000, 400), margin=5Plots.mm)
    return p
end

"""
    plot_borderline(data, idx)

Plot the `idx`-th borderline precursor (sorted by global_qval, near-misses first).
"""
function plot_borderline(data, idx::Int)
    bl = data.borderline
    1 <= idx <= nrow(bl) || error("idx must be in 1:$(nrow(bl))")
    row = bl[idx, :]
    meta = get(data.pidx_meta, UInt32(row.precursor_idx), nothing)
    seq = meta !== nothing ? meta.sequence : "?"
    println("[$idx / $(nrow(bl))]  precursor_idx=$(row.precursor_idx)  " *
            "$(seq) $(row.charge)+  global_qval=$(round(row.global_qval; digits=4))")
    return _bl_plot_one(data, row)
end

"""
    plot_near_miss(data, idx)

Plot the `idx`-th near-miss precursor (global_qval ∈ (0.01, 0.05]).
"""
function plot_near_miss(data, idx::Int)
    nm = data.near_miss
    1 <= idx <= nrow(nm) || error("idx must be in 1:$(nrow(nm))")
    row = nm[idx, :]
    meta = get(data.pidx_meta, UInt32(row.precursor_idx), nothing)
    seq = meta !== nothing ? meta.sequence : "?"
    println("[near-miss $idx / $(nrow(nm))]  precursor_idx=$(row.precursor_idx)  " *
            "$(seq) $(row.charge)+  global_qval=$(round(row.global_qval; digits=4))")
    return _bl_plot_one(data, row)
end

"""
    plot_hopeless(data, idx)

Plot the `idx`-th hopeless precursor (global_qval > 0.25).
"""
function plot_hopeless(data, idx::Int)
    hp = data.hopeless
    1 <= idx <= nrow(hp) || error("idx must be in 1:$(nrow(hp))")
    row = hp[idx, :]
    meta = get(data.pidx_meta, UInt32(row.precursor_idx), nothing)
    seq = meta !== nothing ? meta.sequence : "?"
    println("[hopeless $idx / $(nrow(hp))]  precursor_idx=$(row.precursor_idx)  " *
            "$(seq) $(row.charge)+  global_qval=$(round(row.global_qval; digits=4))")
    return _bl_plot_one(data, row)
end

"""
    plot_borderline_grid(data; n=12, category=:all)

Plot a grid of borderline precursors. category: :all, :near_miss, :moderate, :hopeless
"""
function plot_borderline_grid(data; n::Int=12, category::Symbol=:all)
    df = if category == :near_miss
        data.near_miss
    elseif category == :moderate
        data.moderate
    elseif category == :hopeless
        data.hopeless
    else
        data.borderline
    end

    n = min(n, nrow(df))
    plots = [_bl_plot_one(data, df[i, :]) for i in 1:n]
    nrows = ceil(Int, n / 2)
    return plot(plots...; layout=(nrows, 2), size=(1050, nrows * 380))
end

# ============================================================
# Main
# ============================================================

println("\n+========================================================+")
println("|  Borderline Precursor Analysis (5-Category)           |")
println("+========================================================+")

bl_data = bl_load_data()

println("\nUsage:")
println("  borderline_summary(bl_data)                  # Full 5-category breakdown")
println("  plot_borderline(bl_data, 1)                  # Browse B1 (by global_qval)")
println("  plot_near_miss(bl_data, 1)                   # B1 near-misses (q in (0.01, 0.05])")
println("  plot_hopeless(bl_data, 1)                    # B1 hopeless (q > 0.25)")
println("  plot_borderline_grid(bl_data; n=6)           # Grid view")
println("  plot_borderline_grid(bl_data; n=6, category=:near_miss)")
println()
println("Categories: A=$(length(bl_data.cat_a))  B1=$(length(bl_data.cat_b1))  " *
        "B2=$(length(bl_data.cat_b2))  B_lost=$(length(bl_data.cat_b_lost))  C=$(length(bl_data.cat_c))")
