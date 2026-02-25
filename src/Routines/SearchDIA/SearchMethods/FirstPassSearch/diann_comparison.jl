"""
Compare Pioneer first-pass global precursor probabilities against DIA-NN results.

When a `benchmark_results` path (DIA-NN .parquet) is specified in the
`first_search` parameter section, this module reads the DIA-NN precursors
at 1% global Q-value, maps them to Pioneer precursor indices, and reports:
  - How many DIA-NN precursors are found vs missing in Pioneer
  - The distribution of Pioneer global probabilities for DIA-NN precursors
  - A cumulative curve: global prob cutoff → DIA-NN precursors retained
  - Separate PDFs for LightGBM-rescored and probit-based global scores
"""

# ─────────────────────────────────────────────────────────────────────────────
# Pioneer → DIA-NN precursor ID conversion
# ─────────────────────────────────────────────────────────────────────────────

"""
    _pioneer_to_precursor_id(sequence, structural_mods, charge) -> String

Convert Pioneer library representation (sequence, structural_mods, charge)
to the DIA-NN `Precursor.Id` format:  `SEQM(UniMod:35)ENCE2`.

Pioneer stores mods as `"(pos,AA,Unimod:N)"` tuples; DIA-NN inlines
`"(UniMod:N)"` after the residue at that position and appends the charge.
"""
function _pioneer_to_precursor_id(
    sequence::AbstractString,
    structural_mods,
    charge::UInt8
)
    if ismissing(structural_mods) || structural_mods == ""
        return string(sequence, Int(charge))
    end

    # Parse "(pos,AA,Unimod:N)" groups
    mods = Tuple{Int, String}[]
    for m in eachmatch(r"\((\d+),[A-Z],([Uu]nimod:\d+)\)", structural_mods)
        pos = parse(Int, m.captures[1])
        unimod_tag = replace(m.captures[2], "Unimod:" => "UniMod:")
        push!(mods, (pos, "($(unimod_tag))"))
    end

    # Insert from right to left so positions stay valid
    sort!(mods; by=first, rev=true)
    result = string(sequence)
    for (pos, tag) in mods
        result = string(result[1:pos], tag, result[pos+1:end])
    end

    return string(result, Int(charge))
end

# ─────────────────────────────────────────────────────────────────────────────
# Build Pioneer global probability map
# ─────────────────────────────────────────────────────────────────────────────

"""
    _build_pioneer_global_probs(all_psms_paths, valid_indices, prec_is_decoy, fdr_scale_factor)
        -> (global_prob::Dict{UInt32,Float32}, global_pep::Dict{UInt32,Float32})

Read the first-pass Arrow files and compute a log-odds–combined global
probability and global PEP for every unique precursor (targets + decoys).
"""
function _build_pioneer_global_probs(
    all_psms_paths::Vector{String},
    valid_indices::Vector{Int},
    prec_is_decoy::AbstractVector{Bool},
    fdr_scale_factor::Float32
)
    prec_probs_by_run = Dictionary{UInt32, Vector{Float32}}()
    n_valid_files = 0

    for ms_file_idx in valid_indices
        path = all_psms_paths[ms_file_idx]
        (isempty(path) || !isfile(path)) && continue

        tbl = Arrow.Table(path)
        prec_col = tbl[:precursor_idx]
        prob_col = tbl[:prob]
        isempty(prec_col) && continue
        n_valid_files += 1

        # Best prob per precursor in this file
        file_best = Dictionary{UInt32, Float32}()
        for i in eachindex(prec_col)
            pid = prec_col[i]
            p   = Float32(prob_col[i])
            if haskey(file_best, pid)
                file_best[pid] = max(file_best[pid], p)
            else
                insert!(file_best, pid, p)
            end
        end
        for (pid, p) in pairs(file_best)
            if haskey(prec_probs_by_run, pid)
                push!(prec_probs_by_run[pid], p)
            else
                insert!(prec_probs_by_run, pid, Float32[p])
            end
        end
    end

    n_valid_files == 0 && return (Dict{UInt32,Float32}(), Dict{UInt32,Float32}())

    sqrt_n = max(1, floor(Int, sqrt(n_valid_files)))
    n_unique = length(prec_probs_by_run)

    pids   = Vector{UInt32}(undef, n_unique)
    gprobs = Vector{Float32}(undef, n_unique)
    targets = Vector{Bool}(undef, n_unique)
    for (i, (pid, probs)) in enumerate(pairs(prec_probs_by_run))
        pids[i]    = pid
        gprobs[i]  = _logodds_combine(probs, sqrt_n)
        targets[i] = !prec_is_decoy[pid]
    end

    # Global PEP
    gpeps = Vector{Float32}(undef, n_unique)
    get_PEP!(gprobs, targets, gpeps; doSort=true, fdr_scale_factor=fdr_scale_factor)

    prob_dict = Dict{UInt32,Float32}()
    pep_dict  = Dict{UInt32,Float32}()
    for i in 1:n_unique
        prob_dict[pids[i]] = gprobs[i]
        pep_dict[pids[i]]  = gpeps[i]
    end
    return (prob_dict, pep_dict)
end

# ─────────────────────────────────────────────────────────────────────────────
# Plot + log helper for one set of global scores
# ─────────────────────────────────────────────────────────────────────────────

"""
    _diann_score_analysis(diann_set, id_to_pidx, global_prob, global_pep,
                          global_pep_threshold, label, qc_plot_folder,
                          n_diann) -> nothing

Match DIA-NN precursors against one set of Pioneer global scores, log
diagnostics, and generate a 4-page comparison PDF. `label` is used for
log lines and file names (e.g. "LightGBM" or "probit").
"""
function _diann_score_analysis(
    diann_set::Set{String},
    id_to_pidx::Dict{String, UInt32},
    global_prob::Dict{UInt32, Float32},
    global_pep::Dict{UInt32, Float32},
    global_pep_threshold::Float32,
    label::String,
    qc_plot_folder::String,
    n_diann::Int
)
    # ── Match DIA-NN precursors ──
    ids_not_in_library = String[]
    ids_no_psm         = String[]
    matched_ids   = String[]
    matched_probs = Float32[]
    matched_peps  = Float32[]

    for prec_id in diann_set
        pidx = get(id_to_pidx, prec_id, nothing)
        if pidx === nothing
            push!(ids_not_in_library, prec_id)
            continue
        end
        prob = get(global_prob, pidx, nothing)
        if prob === nothing
            push!(ids_no_psm, prec_id)
            continue
        end
        push!(matched_ids, prec_id)
        push!(matched_probs, prob)
        push!(matched_peps, get(global_pep, pidx, 1.0f0))
    end

    not_in_library = length(ids_not_in_library)
    not_scored     = length(ids_no_psm)
    n_matched      = length(matched_probs)

    @user_info "  [$label] DIA-NN precursor matching:"
    @user_info "    Matched (scored):   $n_matched / $n_diann ($(round(100*n_matched/n_diann; digits=1))%)"
    @user_info "    Not in library:     $not_in_library / $n_diann ($(round(100*not_in_library/n_diann; digits=1))%)"
    @user_info "    In library, no PSM: $not_scored / $n_diann ($(round(100*not_scored/n_diann; digits=1))%)"

    # PEP distribution
    pep_thresholds = [0.01f0, 0.05f0, 0.10f0, 0.25f0, 0.50f0, 0.75f0, 0.90f0, 0.99f0, 1.0f0]
    parts = String[]
    for t in pep_thresholds
        n = count(<=(t), matched_peps)
        push!(parts, "≤$t: $n")
    end
    @user_info "  [$label] DIA-NN by global PEP: " * join(parts, " | ")

    # Prob distribution
    prob_thresholds = [0.50f0, 0.60f0, 0.70f0, 0.80f0, 0.90f0, 0.95f0, 0.99f0]
    prob_parts = String[]
    for t in prob_thresholds
        n = count(>=(t), matched_probs)
        push!(prob_parts, "≥$t: $n")
    end
    @user_info "  [$label] DIA-NN by global prob: " * join(prob_parts, " | ")

    # Threshold report
    n_pass_threshold = count(<=(global_pep_threshold), matched_peps)
    n_fail_threshold = n_matched - n_pass_threshold
    @user_info "  [$label] Global PEP threshold = $global_pep_threshold: " *
        "$n_pass_threshold / $n_matched pass " *
        "($(round(100*n_pass_threshold/n_matched; digits=1))%), " *
        "$n_fail_threshold filtered out"

    # ── Write TSV files ──
    tag = lowercase(label)

    not_in_lib_path = joinpath(qc_plot_folder, "diann_$(tag)_not_in_library.tsv")
    open(not_in_lib_path, "w") do io
        println(io, "Precursor.Id")
        for pid in sort(ids_not_in_library)
            println(io, pid)
        end
    end

    no_psm_path = joinpath(qc_plot_folder, "diann_$(tag)_in_library_no_psm.tsv")
    open(no_psm_path, "w") do io
        println(io, "Precursor.Id")
        for pid in sort(ids_no_psm)
            println(io, pid)
        end
    end

    matched_path = joinpath(qc_plot_folder, "diann_$(tag)_matched_scores.tsv")
    perm = sortperm(matched_peps)
    open(matched_path, "w") do io
        println(io, "Precursor.Id\tglobal_prob\tglobal_PEP\tpass_threshold")
        for i in perm
            pass = matched_peps[i] <= global_pep_threshold ? "true" : "false"
            println(io, matched_ids[i], "\t",
                    round(matched_probs[i]; digits=6), "\t",
                    round(matched_peps[i]; digits=6), "\t", pass)
        end
    end
    @user_info "  [$label] Wrote TSV files to $qc_plot_folder"

    # ── Generate plots ──
    if isempty(matched_probs)
        @user_warn "  [$label] No matched precursors — skipping plots"
        return
    end

    # Plot 1: Cumulative prob curve
    cutoffs = range(0.0f0, 1.0f0; length=200)
    retained = [count(>=(c), matched_probs) for c in cutoffs]

    p1 = Plots.plot(
        Float64.(cutoffs), retained;
        xlabel = "Pioneer global prob cutoff (≥)",
        ylabel = "DIA-NN precursors retained",
        title  = "[$label] DIA-NN vs Pioneer global prob cutoff",
        legend = false, lw = 2, color = :steelblue, size = (700, 450)
    )
    Plots.hline!(p1, [n_diann]; linestyle=:dash, color=:gray,
                  label="Total DIA-NN ($n_diann)")
    Plots.hline!(p1, [n_matched]; linestyle=:dot, color=:orange,
                  label="Matched ($n_matched)")
    n_gap = n_diann - n_matched
    Plots.annotate!(p1, 0.05, n_diann - n_gap/2,
                     Plots.text("$not_in_library not in lib + $not_scored no PSM",
                                8, :left, :gray))

    # Plot 2: Prob histogram
    p2 = Plots.histogram(
        Float64.(matched_probs);
        bins = 50, xlabel = "Pioneer global prob", ylabel = "Count",
        title = "[$label] Global prob for DIA-NN precursors",
        legend = false, color = :steelblue, alpha = 0.7, size = (700, 450)
    )

    # Plot 3: Cumulative PEP curve
    pep_cutoffs = range(0.0f0, 1.0f0; length=200)
    pep_retained = [count(<=(c), matched_peps) for c in pep_cutoffs]

    p3 = Plots.plot(
        Float64.(pep_cutoffs), pep_retained;
        xlabel = "Pioneer global PEP cutoff (≤)",
        ylabel = "DIA-NN precursors retained",
        title  = "[$label] DIA-NN vs Pioneer global PEP cutoff",
        legend = false, lw = 2, color = :firebrick, size = (700, 450)
    )
    Plots.hline!(p3, [n_diann]; linestyle=:dash, color=:gray)
    Plots.hline!(p3, [n_matched]; linestyle=:dot, color=:orange)
    Plots.vline!(p3, [Float64(global_pep_threshold)]; linestyle=:dash, color=:red, lw=2)
    Plots.annotate!(p3, Float64(global_pep_threshold) + 0.02,
                     n_pass_threshold * 0.5,
                     Plots.text("PEP≤$global_pep_threshold\n$n_pass_threshold pass",
                                8, :left, :red))

    # Plot 4: PEP histogram
    p4 = Plots.histogram(
        Float64.(matched_peps);
        bins = 50, xlabel = "Pioneer global PEP", ylabel = "Count",
        title = "[$label] Global PEP for DIA-NN precursors",
        legend = false, color = :firebrick, alpha = 0.7, size = (700, 450)
    )
    Plots.vline!(p4, [Float64(global_pep_threshold)]; linestyle=:dash, color=:red, lw=2)

    # Save
    diann_plots = Plots.Plot[p1, p2, p3, p4]
    output_path = joinpath(qc_plot_folder, "diann_comparison_$(tag).pdf")
    try
        isfile(output_path) && safeRm(output_path, nothing)
    catch e
        @user_warn "Could not clear existing file: $e"
    end
    save_multipage_pdf(diann_plots, output_path)
    @user_info "  [$label] Saved comparison plots to $output_path"
end

# ─────────────────────────────────────────────────────────────────────────────
# Write annotated first-pass PSMs as Arrow table
# ─────────────────────────────────────────────────────────────────────────────

"""
    _write_annotated_first_pass_psms(all_psms_paths, valid_indices,
        diann_pidx_set, global_prob, global_pep, global_pep_threshold,
        label, qc_plot_folder) -> String

Read all first-pass PSM Arrow files, combine into one DataFrame, and annotate
each row with DIA-NN membership, global scores, and threshold pass/fail.
Returns the path to the written Arrow file.
"""
function _write_annotated_first_pass_psms(
    all_psms_paths::Vector{String},
    valid_indices::Vector{Int},
    diann_pidx_set::Set{UInt32},
    global_prob::Dict{UInt32, Float32},
    global_pep::Dict{UInt32, Float32},
    global_pep_threshold::Float32,
    label::String,
    qc_plot_folder::String
)
    # Read and combine all first-pass PSM Arrow files
    dfs = DataFrame[]
    for ms_file_idx in valid_indices
        path = all_psms_paths[ms_file_idx]
        (isempty(path) || !isfile(path)) && continue
        tbl = Arrow.Table(path)
        isempty(Tables.columns(tbl)) && continue
        push!(dfs, DataFrame(tbl))
    end

    if isempty(dfs)
        @user_warn "  [$label] No first-pass PSM files found — skipping Arrow table"
        return ""
    end

    combined = vcat(dfs...; cols=:union)
    @user_info "  [$label] Combined $(length(dfs)) files → $(nrow(combined)) total PSMs"

    # Annotate: is_diann, global_prob, global_pep, pass_threshold
    n = nrow(combined)
    prec_col = combined[!, :precursor_idx]

    is_diann_col    = Vector{Bool}(undef, n)
    gprob_col       = Vector{Union{Missing, Float32}}(undef, n)
    gpep_col        = Vector{Union{Missing, Float32}}(undef, n)
    pass_thresh_col = Vector{Union{Missing, Bool}}(undef, n)

    for i in 1:n
        pidx = prec_col[i]
        is_diann_col[i] = pidx in diann_pidx_set

        gp = get(global_prob, pidx, nothing)
        if gp !== nothing
            gprob_col[i] = gp
            pe = get(global_pep, pidx, 1.0f0)
            gpep_col[i] = pe
            pass_thresh_col[i] = pe <= global_pep_threshold
        else
            gprob_col[i] = missing
            gpep_col[i] = missing
            pass_thresh_col[i] = missing
        end
    end

    combined[!, :is_diann]        = is_diann_col
    combined[!, :global_prob]     = gprob_col
    combined[!, :global_pep]      = gpep_col
    combined[!, :pass_threshold]  = pass_thresh_col

    tag = lowercase(label)
    out_path = joinpath(qc_plot_folder, "diann_annotated_first_pass_psms_$(tag).arrow")
    Arrow.write(out_path, combined)

    n_diann_psms = count(is_diann_col)
    n_diann_pass = count(i -> is_diann_col[i] && !ismissing(pass_thresh_col[i]) && pass_thresh_col[i], 1:n)
    @user_info "  [$label] Annotated Arrow: $n PSMs total, $n_diann_psms DIA-NN PSMs, $n_diann_pass DIA-NN passing threshold"
    @user_info "  [$label] Wrote annotated PSMs to $out_path"

    return out_path
end

# ─────────────────────────────────────────────────────────────────────────────
# Main comparison function
# ─────────────────────────────────────────────────────────────────────────────

"""
    compare_first_pass_with_diann(search_context, fdr_scale_factor, diann_path,
                                  qc_plot_folder, global_pep_threshold,
                                  probit_global_probs)

Read a DIA-NN `.parquet` report, match precursors to Pioneer's first-pass
global probabilities (LightGBM-rescored), log diagnostics, save comparison
plots. If `probit_global_probs` is provided, also generates a separate PDF
using the pre-LightGBM probit-based global scores for comparison.
"""
function compare_first_pass_with_diann(
    search_context::SearchContext,
    fdr_scale_factor::Float32,
    diann_path::String,
    qc_plot_folder::String,
    global_pep_threshold::Float32,
    probit_global_probs  # (prob_dict, pep_dict) or nothing
)
    if !isfile(diann_path)
        @user_warn "Benchmark comparison: DIA-NN file not found at $diann_path"
        return
    end

    @user_info "Benchmark comparison: loading DIA-NN results from $diann_path"

    # ── 1. Read DIA-NN precursors at 1% global Q ──
    diann_df = DataFrame(Parquet2.Dataset(diann_path))
    diann_pass = diann_df[diann_df[:, "Global.Q.Value"] .<= 0.01f0, :]
    diann_prec_ids = unique(diann_pass[:, "Precursor.Id"])
    diann_set = Set{String}(diann_prec_ids)
    n_diann = length(diann_set)
    @user_info "  DIA-NN: $n_diann unique precursors at 1% global Q-value"

    # ── 2. Build Pioneer precursor_id → precursor_idx lookup ──
    precursors = getPrecursors(getSpecLib(search_context))
    seqs       = getSequence(precursors)
    struct_mods = getStructuralMods(precursors)
    charges    = getCharge(precursors)
    is_decoy   = getIsDecoy(precursors)
    entrap_ids = precursors.data[:entrapment_group_id]

    id_to_pidx = Dict{String, UInt32}()
    for pidx in UInt32(1):UInt32(length(seqs))
        is_decoy[pidx] && continue
        entrap_ids[pidx] != 0x00 && continue
        prec_id = _pioneer_to_precursor_id(seqs[pidx], struct_mods[pidx], charges[pidx])
        if !haskey(id_to_pidx, prec_id)
            id_to_pidx[prec_id] = pidx
        end
    end
    @user_info "  Pioneer library: $(length(id_to_pidx)) unique target precursor IDs"

    mkpath(qc_plot_folder)

    # ── 3. LightGBM-rescored analysis (current Arrow files) ──
    all_psms_paths = getFirstPassPsms(getMSData(search_context))
    valid_indices  = get_valid_file_indices(search_context)
    prec_is_decoy  = getIsDecoy(precursors)

    (lgbm_prob, lgbm_pep) = _build_pioneer_global_probs(
        all_psms_paths, valid_indices, prec_is_decoy, fdr_scale_factor)

    n_lgbm_scored = count(pid -> !prec_is_decoy[pid], keys(lgbm_prob))
    @user_info "  Pioneer (LightGBM): $n_lgbm_scored target precursors with global scores"

    _diann_score_analysis(diann_set, id_to_pidx, lgbm_prob, lgbm_pep,
                          global_pep_threshold, "LightGBM", qc_plot_folder, n_diann)

    # Build set of precursor indices that correspond to DIA-NN identifications
    diann_pidx_set = Set{UInt32}()
    for prec_id in diann_set
        pidx = get(id_to_pidx, prec_id, nothing)
        pidx !== nothing && push!(diann_pidx_set, pidx)
    end

    # Write annotated Arrow table with LightGBM global scores
    _write_annotated_first_pass_psms(
        all_psms_paths, valid_indices, diann_pidx_set,
        lgbm_prob, lgbm_pep, global_pep_threshold,
        "LightGBM", qc_plot_folder)

    # ── 4. Probit-based analysis (pre-LightGBM snapshot) ──
    if probit_global_probs !== nothing
        (probit_prob, probit_pep) = probit_global_probs
        n_probit_scored = count(pid -> !prec_is_decoy[pid], keys(probit_prob))
        @user_info "  Pioneer (probit): $n_probit_scored target precursors with global scores"

        _diann_score_analysis(diann_set, id_to_pidx, probit_prob, probit_pep,
                              global_pep_threshold, "probit", qc_plot_folder, n_diann)

        # Write annotated Arrow table with probit global scores
        _write_annotated_first_pass_psms(
            all_psms_paths, valid_indices, diann_pidx_set,
            probit_prob, probit_pep, global_pep_threshold,
            "probit", qc_plot_folder)
    end
end
