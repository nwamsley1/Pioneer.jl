# Scanning-quad (ZT) main search with more than one chunk: the meta-PSMs live on disk, merged into
# precursor-complete files, and MainSearch scores them file by file.
#
# The result is the same as scoring the whole table in memory. The in-memory table is the chunks in
# cycle order, stably sorted by precursor, i.e. ordered by (precursor_idx, scan_idx); the merged
# files hold exactly that order. Every step is per row or per precursor, apart from LightGBM's
# training sample, which is drawn by the same global row positions, and the iRT correction fit,
# which uses one refinement row per precursor.

"""Where a ZT file's on-disk meta-PSMs live during its main search."""
zt_psm_dir(search_context::SearchContext, ms_file_idx::Integer) =
    joinpath(getDataOutDir(search_context), "temp_data", "zt_meta_psms", "file_$(ms_file_idx)")

"""
Working bytes per meta-PSM while a merged file is scored: the featured row (~311 B for 97 columns)
plus its LightGBM feature matrix and the fold gather in prediction. Sizes the merged files from
`max_psm_memory_mb`.
"""
const ZT_PARTITION_BYTES_PER_ROW = 900

"""Rows per merged meta-PSM file: the PSM memory budget over `ZT_PARTITION_BYTES_PER_ROW`."""
zt_partition_rows(search_context::SearchContext, ms_file_idx::Integer) =
    max(250_000, floor(Int, zt_psm_memory_mb(search_context, ms_file_idx) * 1e6 / ZT_PARTITION_BYTES_PER_ROW))

"""Below this many meta-PSMs a multi-chunk file is scored in memory: it is small, and LightGBM's
small-data model competition needs the whole table."""
const ZT_PARTITION_MIN_ROWS = 400_000

@inline _pid_scan_key(pc, sc, r) = (UInt64(pc[r]) << 32) | UInt64(sc[r])

"""Throw unless `df` is sorted by (precursor_idx, scan_idx): the merge assumes it."""
function _assert_precursor_scan_sorted(df::DataFrame)
    pc = df[!, :precursor_idx]::Vector{UInt32}; sc = df[!, :scan_idx]::Vector{UInt32}
    @inbounds for r in 2:length(pc)
        _pid_scan_key(pc, sc, r - 1) <= _pid_scan_key(pc, sc, r) ||
            error("ZT meta-PSM chunk is not sorted by (precursor_idx, scan_idx) at row $r")
    end
    return nothing
end

"""
    precursor_part_ranges(pcs, target_rows) -> Vector{Tuple{UInt32, UInt32}}

The precursor-id ranges [lo, hi] of the parts a merge of files with precursor columns `pcs` (each
sorted) is cut into: walking all rows in precursor order, a part is closed before the next
precursor once it holds `target_rows` rows, so it ends only where the precursor changes.
"""
function precursor_part_ranges(pcs::Vector{Vector{UInt32}}, target_rows::Int)
    allp = reduce(vcat, pcs; init = UInt32[])
    sort!(allp)
    ranges = Tuple{UInt32, UInt32}[]
    acc = 0; lo = UInt32(0); prev = UInt32(0); open = false
    i = 1; n = length(allp)
    @inbounds while i <= n
        pid = allp[i]; j = i
        while j < n && allp[j + 1] == pid
            j += 1
        end
        if open && acc >= target_rows          # close the part before this precursor
            push!(ranges, (lo, prev)); acc = 0; open = false
        end
        open || (lo = pid; open = true)
        acc += j - i + 1; prev = pid; i = j + 1
    end
    open && push!(ranges, (lo, prev))
    return ranges
end

# function barrier: the slices of one column (every element of `cols` the same concrete Arrow column
# type), concatenated in file order, then put in `perm` order
function _slice_gather(cols::Vector{C}, slices::Vector{UnitRange{Int}}, perm::AbstractVector{<:Integer},
                       n::Int) where {C}
    T = eltype(first(cols))
    cat = Vector{T}(undef, n); off = 0
    @inbounds for (f, rng) in enumerate(slices)
        col = cols[f]
        for r in rng
            off += 1; cat[off] = col[r]
        end
    end
    out = Vector{T}(undef, n)
    @inbounds for i in 1:n
        out[i] = cat[perm[i]]
    end
    return out
end

"""
    zt_merge_by_precursor(on_part, paths; target_rows) -> Int

Merge Arrow files, each sorted by (precursor_idx, scan_idx), and call `on_part(df)` with
consecutive parts of about `target_rows` rows. Each part is sorted the same way and holds every row
of each of its precursors: a part ends only where the precursor changes (`precursor_part_ranges`).
Returns the number of parts.

A part covering precursors [lo, hi] is one contiguous slice of every file (two binary searches on
its sorted precursor column), so each part is built from slices: concatenated, then sorted by the
(precursor_idx, scan_idx) key, which is unique per row; the columns are gathered in parallel. The
parts and their rows are exactly those of a streaming k-way merge, without a heap operation per row.

Checks every part before handing it on: sorted, and starting after the previous part's last
precursor, so each precursor is in exactly one part; and that the parts hold every input row.
"""
function zt_merge_by_precursor(on_part, paths::Vector{String}; target_rows::Int)
    tabs = [Arrow.Table(p) for p in paths]
    pcs = [Vector{UInt32}(t.precursor_idx) for t in tabs]
    lens = [length(p) for p in pcs]
    names_ = collect(Tables.columnnames(first(tabs)))
    colsets = [[Tables.getcolumn(t, nm) for t in tabs] for nm in names_]
    scan_cols = [Tables.getcolumn(t, :scan_idx) for t in tabs]
    ranges = precursor_part_ranges(pcs, target_rows)
    workspace = Int32SortPermWorkspace()
    n_parts = 0; n_out = 0; last_pid = UInt32(0)
    for (lo, hi) in ranges
        slices = [searchsortedfirst(pcs[f], lo):searchsortedlast(pcs[f], hi) for f in eachindex(pcs)]
        n = sum(length, slices)
        keys = Vector{UInt64}(undef, n); off = 0
        @inbounds for (f, rng) in enumerate(slices)
            pc = pcs[f]; sc = scan_cols[f]
            for r in rng
                off += 1; keys[off] = (UInt64(pc[r]) << 32) | UInt64(sc[r])
            end
        end
        perm = parallel_sortperm_int32!(workspace, keys)
        cols = Vector{AbstractVector}(undef, length(names_))
        Threads.@threads for j in eachindex(names_)
            cols[j] = _slice_gather(colsets[j], slices, perm, n)
        end
        df = DataFrame(cols, names_; copycols = false)
        _check_precursor_part(df, last_pid, n_parts + 1)
        last_pid = last(df[!, :precursor_idx]); n_parts += 1; n_out += n
        on_part(df)
    end
    n_out == sum(lens) || error("ZT merge produced $n_out rows, expected $(sum(lens))")
    return n_parts
end

"""Throw unless part `i` is sorted by (precursor_idx, scan_idx) and starts after `last_pid`, the
previous part's last precursor."""
function _check_precursor_part(df::DataFrame, last_pid::UInt32, i::Int)
    pc = df[!, :precursor_idx]::Vector{UInt32}; sc = df[!, :scan_idx]::Vector{UInt32}
    i > 1 && first(pc) <= last_pid &&
        error("ZT merge: precursor $(first(pc)) is split across parts $(i - 1) and $i")
    @inbounds for r in 2:length(pc)
        _pid_scan_key(pc, sc, r - 1) <= _pid_scan_key(pc, sc, r) ||
            error("ZT merge: part $i is not sorted at row $r")
    end
    return nothing
end

"""A merged or featured meta-PSM file as an in-memory DataFrame with ordinary vector columns."""
_zt_load(path::String) = DataFrame(Arrow.Table(path); copycols = true)

"""
    zt_mainsearch_best_partitioned!(parts, results, params, search_context, ms_file_idx, spectra,
                                    center_mzs, isolation_widths, bitvec_rank_table)

First half of MainSearch's per-file scoring (`_mainsearch_best_in_memory!`) over precursor-complete
meta-PSM files, one file in memory at a time:

- pass A: prescore and chromatogram features; write the featured file
- LightGBM: the same 250k-per-fold training sample, gathered by global row position
- pass B: out-of-fold scores; the iRT refinement rows (one per precursor); fit the correction
- pass C: apply the correction, re-predict, best meta-PSM per precursor

Returns the same stage as `_mainsearch_best_in_memory!`; its `trace_input` is pass D, which
gathers only the rows of the precursors the trace features see.
"""
function zt_mainsearch_best_partitioned!(parts::Vector{String}, results::MainSearchResults,
                                         params::MainSearchParameters,
                                         search_context::SearchContext, ms_file_idx::Int64,
                                         spectra::MassSpecData, center_mzs, isolation_widths,
                                         bitvec_rank_table)
    dir = dirname(first(parts))                       # the chunk files, sorted by (precursor, scan)
    cleanup = () -> (_zt_rm_dir(dir); setZTPsmPartitions!(search_context, ms_file_idx, nothing))
    n_meta = sum(p -> length(Arrow.Table(p).precursor_idx), parts)
    if n_meta < ZT_PARTITION_MIN_ROWS
        # the in-memory table: chunks in cycle order, stably sorted by precursor (as process_file!)
        psms = n_meta == 0 ? DataFrame() : reduce(vcat, (_zt_load(p) for p in parts))
        nrow(psms) > 0 && permute_psms_by_precursor_idx!(psms, results.sortperm_workspace)
        results.psms[] = psms
        cleanup()
        return _mainsearch_best_in_memory!(results, params, search_context, ms_file_idx, spectra,
                                           center_mzs, isolation_widths, bitvec_rank_table)
    end
    buffers = results.lgbm_buffers

    # ---- merge + pass A: each precursor-complete part gets its prescore and chromatogram
    # features in memory and is written once, featured ----
    featured = String[]; n_rows = Int[]; t_prepare = 0.0; t_ms1 = 0.0
    t_merge = @elapsed zt_merge_by_precursor(parts; target_rows = zt_partition_rows(search_context, ms_file_idx)) do df
        t_prepare += @elapsed prepare_psm_features!(df, params, search_context, ms_file_idx, spectra)
        t_ms1 += @elapsed add_chromatogram_features!(df, spectra; bitvec_rank_table = bitvec_rank_table)
        f = joinpath(dir, "part_$(lpad(length(featured) + 1, 4, '0')).featured.arrow")
        Arrow.write(f, df); push!(featured, f); push!(n_rows, nrow(df))
    end
    foreach(safeRm, parts)                            # Windows: chunks may still be memory-mapped
    n_total = sum(n_rows)
    @user_info "ZT main search: merged $(length(parts)) chunk files into $(length(featured)) " *
               "precursor-complete parts in $(round(t_merge; digits = 1))s " *
               "(features $(round(t_prepare + t_ms1; digits = 1))s)"

    # ---- LightGBM: as _train_psm_classifier_with_fallback, without the whole-table matrix ----
    t_lgbm_start = time()
    lazy = DataFrame(Arrow.Table(featured); copycols = false)       # memory-mapped, read-only
    available_features = filter(f -> hasproperty(lazy, f),
                                model_features(PRESCORE_FEATURES, has_ion_mobility(search_context)))
    if :num_enzymatic_termini in available_features
        et = lazy[!, :num_enzymatic_termini]; fv = isempty(et) ? nothing : first(et)
        (isempty(et) || all(v -> isequal(v, fv), et)) &&
            deleteat!(available_features, findfirst(==(:num_enzymatic_termini), available_features))
    end
    all_targets = Vector{Bool}(lazy[!, :target])
    cv_fold = Vector{UInt8}(lazy[!, :cv_fold])
    idx0 = findall(cv_fold .== 0); idx1 = findall(cv_fold .== 1)
    fold_pairs = [(idx1, idx0), (idx0, idx1)]                         # (fit rows, scored rows)
    _sample_pos(n_avail) = n_avail > MAIN_LGBM_MAX_TRAIN ? randperm(n_avail)[1:MAIN_LGBM_MAX_TRAIN] :
                                                           collect(1:n_avail)
    sub_positions = [_sample_pos(length(fit_idx)) for (fit_idx, _) in fold_pairs]
    min_fit = minimum(length(fit_idx) for (fit_idx, _) in fold_pairs)
    min_fit > 50_000 || error("ZT partitioned scoring needs > 50,000 meta-PSMs per fold, got $min_fit")
    @debug_l1 "  LightGBM CV (partitioned, $(length(parts)) files): fold0=$(length(idx0)) fold1=$(length(idx1)) PSMs; " *
              "train=$(length.(sub_positions))  MAX_TRAIN=$MAIN_LGBM_MAX_TRAIN"
    fold_predictors = Vector{Any}(undef, 2)
    t_fit = 0.0
    for (fi, (fit_idx, _)) in enumerate(fold_pairs)
        sub_pos = fit_idx[sub_positions[fi]]
        y_lbl = _prepare_labels(all_targets[sub_pos])
        if isempty(y_lbl) || length(unique(y_lbl)) == 1
            v = isempty(y_lbl) || y_lbl[1] == 0 ? 0.0 : 1.0
            fold_predictors[fi] = (kind = :constant, value = v, model = nothing, beta = nothing)
            continue
        end
        cls = build_lightgbm_classifier(; MAINSEARCH_LGBM_HP...)
        b_train = buffers.train
        GC.@preserve b_train begin
            X_tr = feature_matrix!(b_train, lazy[sub_pos, available_features], available_features)
            t_fit += @elapsed LightGBM.fit!(cls, X_tr, y_lbl; verbosity = -1)
        end
        fold_predictors[fi] = (kind = :lgbm, value = NaN, model = cls, beta = nothing)
    end
    predictor = (available_features = available_features, fold_predictors = fold_predictors)
    lazy = nothing

    # ---- pass B: out-of-fold scores and the iRT refinement rows ----
    scores_oof = Vector{Float32}(undef, n_total); refinement_parts = DataFrame[]; off = 0
    for (i, f) in enumerate(featured)
        df = _zt_load(f)
        df[!, :lgbm_score] = Float32.(predict_psm_classifier_scores(df, predictor; buffers = buffers))
        scores_oof[(off + 1):(off + n_rows[i])] = df[!, :lgbm_score]
        push!(refinement_parts, _select_irt_refinement_psms(df))
        off += n_rows[i]
    end
    refinement_psms = reduce(vcat, refinement_parts)
    t_train_cv = time() - t_lgbm_start
    precursors = getPrecursors(getSpecLib(search_context))
    strategy = MainSearchIrtRefinement(precursors;
                                       q_value_threshold = PRESCORE_QVALUE_THRESHOLD,
                                       min_precursors = MAIN_IRT_REFINEMENT_MIN_PRECURSORS)
    models, training_target_precursors, _, _ = fit_mainsearch_irt_refinement(refinement_psms, strategy)
    refined = !isempty(models)
    t_lgbm_end = time()

    # ---- pass C: refined iRT, re-predicted scores, best meta-PSM per precursor ----
    t_c = time()
    all_scores = Vector{Float32}(undef, n_total); best_parts = DataFrame[]; off = 0
    for (i, f) in enumerate(featured)
        df = _zt_load(f)
        rows = (off + 1):(off + n_rows[i])
        if refined
            apply_mainsearch_irt_refinement_model!(df, strategy, models)
            df[!, :lgbm_score] = Float32.(predict_psm_classifier_scores(df, predictor; buffers = buffers))
        else
            df[!, :lgbm_score] = scores_oof[rows]
        end
        all_scores[rows] = df[!, :lgbm_score]
        push!(best_parts, select_best_per_precursor(df; center_mzs = center_mzs,
                                                    isolation_widths = isolation_widths))
        off += n_rows[i]
    end
    best_psms = reduce(vcat, best_parts)
    best_psms[!, :lgbm_prob] = copy(best_psms[!, :lgbm_score])
    @debug_l1 "  iRT refinement (file_idx=$ms_file_idx, partitioned): " *
              (refined ? "$(length(training_target_precursors)) training precursors" :
                         "skipped ($(length(training_target_precursors)) high-confidence target precursors)") *
              "; pass C $(round(time() - t_c, digits = 2))s"

    trace_input = (best, mask, peps) ->
        _zt_trace_rows(featured, n_rows, best, mask, peps, refined, strategy, models, all_scores)
    return (
        best_psms = best_psms,
        all_scores = all_scores,
        all_targets = all_targets,
        trace_input = trace_input,
        cleanup = cleanup,
        n_total_psms = n_total,
        timings = (prepare = t_prepare, competition = 0.0, apex = 0.0, ms1 = t_ms1,
                   lgbm_start = t_lgbm_start, lgbm_end = t_lgbm_end,
                   lgbm = (train_cv = t_train_cv, best = t_lgbm_end - t_lgbm_start - t_train_cv)),
    )
end

"""
Pass D: the rows the trace features read, i.e. every meta-PSM of the precursors in `best`, in their
post-refinement state, with the matching entries of the per-row pass mask and PEPs.
"""
function _zt_trace_rows(featured, n_rows, best::DataFrame, mask::AbstractVector{Bool}, peps,
                        refined, strategy, models, all_scores)
    keep = Set(best[!, :precursor_idx]::Vector{UInt32})
    tabs = DataFrame[]; sel = Int[]; off = 0
    for (i, f) in enumerate(featured)
        pid = Arrow.Table(f).precursor_idx
        rows = findall(in(keep), pid)
        if !isempty(rows)
            df = _zt_load(f)[rows, :]
            refined && apply_mainsearch_irt_refinement_model!(df, strategy, models)
            df[!, :lgbm_score] = all_scores[off .+ rows]
            push!(tabs, df); append!(sel, off .+ rows)
        end
        off += n_rows[i]
    end
    return (isempty(tabs) ? DataFrame() : reduce(vcat, tabs)), mask[sel], peps[sel]
end

"""Remove the ZT chunk directory; on Windows a still-mapped Arrow file can block this, so a failure
only warns (the directory is under temp_data)."""
function _zt_rm_dir(dir::String)
    try
        rm(dir; recursive = true, force = true)
    catch e
        GC.gc(true)
        try rm(dir; recursive = true, force = true) catch
            @user_warn "Could not remove ZT chunk directory $dir: $(sprint(showerror, e))"
        end
    end
end
