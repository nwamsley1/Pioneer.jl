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

# function barrier: every element of `cols` is the same concrete Arrow column type
function _zt_gather!(dst::Vector, cols::Vector{C}, bf::Vector{Int32}, br::Vector{Int32}, n::Int) where {C}
    @inbounds for i in 1:n; dst[i] = cols[bf[i]][br[i]]; end
    return dst
end

"""
    zt_merge_by_precursor(paths, outdir; target_rows, batch_rows = 1_000_000) -> Vector{String}

k-way merge of Arrow files, each sorted by (precursor_idx, scan_idx), into files of about
`target_rows` rows that are sorted the same way and hold every row of each of their precursors: a
file ends only where the precursor changes. Keys are packed into a UInt64 and merged through a
(key, file) heap; rows are copied a batch at a time, column by column.

Checks the result: every file sorted, and each file's precursors all after the previous file's,
so each precursor is in exactly one file.
"""
function zt_merge_by_precursor(paths::Vector{String}, outdir::String;
                               target_rows::Int, batch_rows::Int = 1_000_000)
    rm(outdir; recursive = true, force = true); mkpath(outdir)
    tabs = [Arrow.Table(p) for p in paths]
    pcs = [t.precursor_idx for t in tabs]; scs = [t.scan_idx for t in tabs]
    lens = [length(p) for p in pcs]
    names_ = collect(Tables.columnnames(first(tabs)))
    colsets = [[Tables.getcolumn(t, nm) for t in tabs] for nm in names_]
    bufs = [Vector{eltype(first(cs))}(undef, batch_rows) for cs in colsets]
    heap = BinaryMinHeap{Tuple{UInt64, Int32}}()
    pos = ones(Int, length(tabs))
    for f in eachindex(tabs)
        lens[f] > 0 && push!(heap, (_pid_scan_key(pcs[f], scs[f], 1), Int32(f)))
    end
    bf = Vector{Int32}(undef, batch_rows); br = Vector{Int32}(undef, batch_rows)
    out_paths = String[]; writer = nothing; rows_in_file = 0; prev_pid = typemax(UInt32); n = 0
    function flush_batch!()
        n == 0 && return
        for j in eachindex(names_); _zt_gather!(bufs[j], colsets[j], bf, br, n); end
        # copies: Arrow's writer reads the columns after this returns, and the buffers are reused
        Arrow.write(writer, DataFrame([nm => bufs[j][1:n] for (j, nm) in enumerate(names_)]; copycols = false))
        n = 0
    end
    function new_file!()
        writer === nothing || close(writer)
        p = joinpath(outdir, "part_$(lpad(length(out_paths) + 1, 4, '0')).arrow"); push!(out_paths, p)
        writer = open(Arrow.Writer, p); rows_in_file = 0
    end
    new_file!()
    while !isempty(heap)
        key, f = pop!(heap)
        pid = UInt32(key >> 32)
        if rows_in_file >= target_rows && pid != prev_pid          # only between precursors
            flush_batch!(); new_file!()
        end
        r = pos[f]; n += 1; bf[n] = f; br[n] = Int32(r)
        rows_in_file += 1; prev_pid = pid; pos[f] = r + 1
        r + 1 <= lens[f] && push!(heap, (_pid_scan_key(pcs[f], scs[f], r + 1), f))
        n == batch_rows && flush_batch!()
    end
    flush_batch!(); close(writer)
    _check_precursor_partitions(out_paths, sum(lens))
    return out_paths
end

"""Throw unless the merged files are each sorted, hold `n_expected` rows in total, and have
disjoint, increasing precursor ranges."""
function _check_precursor_partitions(paths::Vector{String}, n_expected::Int)
    last_pid = UInt32(0); n = 0
    for (i, p) in enumerate(paths)
        t = Arrow.Table(p); pc = t.precursor_idx; sc = t.scan_idx
        isempty(pc) && continue
        n += length(pc)
        i > 1 && first(pc) <= last_pid &&
            error("ZT merge: precursor $(first(pc)) is split across files $(i - 1) and $i")
        @inbounds for r in 2:length(pc)
            _pid_scan_key(pc, sc, r - 1) <= _pid_scan_key(pc, sc, r) ||
                error("ZT merge: $(basename(p)) is not sorted at row $r")
        end
        last_pid = last(pc)
    end
    n == n_expected || error("ZT merge wrote $n rows, expected $n_expected")
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
    dir = dirname(dirname(first(parts)))
    cleanup = () -> (rm(dir; recursive = true, force = true); setZTPsmPartitions!(search_context, ms_file_idx, nothing))
    n_meta = sum(p -> length(Arrow.Table(p).precursor_idx), parts)
    if n_meta < ZT_PARTITION_MIN_ROWS
        results.psms[] = n_meta == 0 ? DataFrame() : reduce(vcat, (_zt_load(p) for p in parts))
        cleanup()
        return _mainsearch_best_in_memory!(results, params, search_context, ms_file_idx, spectra,
                                           center_mzs, isolation_widths, bitvec_rank_table)
    end
    buffers = results.lgbm_buffers

    # ---- pass A: prescore + chromatogram features, per precursor-complete file ----
    featured = [replace(p, r"\.arrow$" => ".featured.arrow") for p in parts]
    n_rows = zeros(Int, length(parts)); t_prepare = 0.0; t_ms1 = 0.0
    for (i, p) in enumerate(parts)
        df = _zt_load(p)
        t_prepare += @elapsed prepare_psm_features!(df, params, search_context, ms_file_idx, spectra)
        t_ms1 += @elapsed add_chromatogram_features!(df, spectra; bitvec_rank_table = bitvec_rank_table)
        n_rows[i] = nrow(df)
        Arrow.write(featured[i], df); rm(p)
    end
    n_total = sum(n_rows)

    # ---- LightGBM: as _train_psm_classifier_with_fallback, without the whole-table matrix ----
    t_lgbm_start = time()
    lazy = DataFrame(Arrow.Table(featured); copycols = false)       # memory-mapped, read-only
    available_features = filter(f -> hasproperty(lazy, f), collect(PRESCORE_FEATURES))
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
