# MBR streaming pipeline. Replaces an older in-memory chain that held a single
# big best_psms DataFrame across MBR feature compute. Normal rows are narrowed
# to candidate-only FTR rows before MBR donor features are computed; only the
# final recovery decisions are expanded back to per-file sidecars.
#
# Two sidecar types per normal main file `<f>.arrow`:
#   <f>.pass1_sidecar.arrow      Pass-1 LightGBM scores (trace_prob_prepass +
#                                trace_prob_infold). Written first so best_psms
#                                can be freed before the heavy MBR stage.
#   <f>.recovery_sidecar.arrow   FTR controller output (mbr_recovered, ftr_qval,
#                                ftr_pep, MBR_transfer_candidate).
#
# Rescue files skip dense Pass-1/MBR/recovery sidecars. They are streamed into
# candidate-only FTR rows once the normal-row donor dict and donor score floor
# are known, and only recovered rescue rows are written back as
# `<f>.arrow.mbr_recovered.arrow` for the final fold merge.
#
# Two main-file sweeps drive MBR feature compute:
#   Sweep 1 (build_mbr_donor_dict_streaming_with_pass1)
#     Streams each (main + pass1 sidecar) pair to build the donor dict keyed
#     by precursor_idx. Each pid keeps the two best donor rows from distinct
#     files, sorted by score, so `_donor_for_pid` can choose the best cross-file
#     donor without storage growing with file count.
#
#   Sweep 2 (load_normal_mbr_candidate_slim_dataframe)
#     Computes pre-MBR row q-values from Pass-1 scores, keeps only rows that
#     failed that first q-value gate and have both cross-file donors, then
#     computes MBR_* donor features for those sparse candidates.
#
# Alignment invariant: row N of every sidecar ↔ row N of the main file
# (positional join). `:precursor_idx` and `:scan_idx` are written into each
# sidecar redundantly and asserted equal at every read/merge, so any sort or
# filter that breaks alignment fails loud rather than silently corrupting data.
#
# merge_mbr_sidecars_into_main! at the end folds the normal sidecars into
# a single mbr_outputs sidecar and deletes the intermediates. Sparse rescue
# recovered files are already full rows and are consumed directly by the final
# fold merge.

# A single cross-run donor entry: one row's score + the per-row donor
# features we need to compute MBR features for a recipient row in
# another file. Holds the bare minimum so the donor dict stays small.
struct _MBRDonorEntry
    trace_prob::Float32
    precursor_idx::UInt32
    weight::Float32
    log2_intensity_explained::Float32
    irt_residual::Float32  # irt_pred − irt_obs of the donor row
    irt_obs::Float32       # raw observed iRT of the donor row
    log_by_ratio::Float32  # log(b_int+1) − log(y_int+1) of the donor row
    n_scans::Float32
    smoothed_frag_sqrt::NTuple{8, Float32}
    ms_file_idx::UInt32
    is_decoy::Bool
end

const MBR_MAX_DONOR_FILES_PER_PRECURSOR = 2
const MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT = ntuple(_ -> 0.0f0, 8)
const MBR_COUNTERFACTUAL_LOCAL_IRT_WINDOW = 1.0f0

struct _MBRFragmentAnnotationKeys
    keys::Vector{UInt16}
end

@inline function _mbr_smoothed_fragment_intensity(frag_cols, row::Integer, rank::Integer)
    value = Float32(frag_cols[rank][row])
    return isfinite(value) ? max(value, 0.0f0) : 0.0f0
end

@inline function _mbr_smoothed_spectrum_sqrt_tuple(frag_cols, row::Integer)
    f1 = _mbr_smoothed_fragment_intensity(frag_cols, row, 1)
    f2 = _mbr_smoothed_fragment_intensity(frag_cols, row, 2)
    f3 = _mbr_smoothed_fragment_intensity(frag_cols, row, 3)
    f4 = _mbr_smoothed_fragment_intensity(frag_cols, row, 4)
    f5 = _mbr_smoothed_fragment_intensity(frag_cols, row, 5)
    f6 = _mbr_smoothed_fragment_intensity(frag_cols, row, 6)
    f7 = _mbr_smoothed_fragment_intensity(frag_cols, row, 7)
    f8 = _mbr_smoothed_fragment_intensity(frag_cols, row, 8)
    frag_sum = f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8
    frag_sum > 0.0f0 || return MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT

    inv_sum = 1.0f0 / frag_sum
    return (
        sqrt(f1 * inv_sum),
        sqrt(f2 * inv_sum),
        sqrt(f3 * inv_sum),
        sqrt(f4 * inv_sum),
        sqrt(f5 * inv_sum),
        sqrt(f6 * inv_sum),
        sqrt(f7 * inv_sum),
        sqrt(f8 * inv_sum),
    )
end

@inline function _mbr_smoothed_spectrum_sqrt_is_valid(frag_sqrt::NTuple{8, Float32})
    @inbounds for rank in 1:8
        frag_sqrt[rank] > 0.0f0 && return true
    end
    return false
end

@inline function _mbr_fragment_annotation_key(frag)
    ion_type = isY(frag) ? UInt16(1) : isB(frag) ? UInt16(2) : isP(frag) ? UInt16(3) : UInt16(0)
    ion_type == UInt16(0) && return UInt16(0)
    return ion_type |
           (UInt16(getFragCharge(frag)) << 2) |
           (UInt16(getIonPosition(frag)) << 4)
end

function build_mbr_fragment_annotation_keys(frag_lookup::LibraryFragmentLookup)
    n_precursors = length(frag_lookup.prec_frag_ranges) - 1
    keys = zeros(UInt16, 8 * n_precursors)
    fragments = getFragments(frag_lookup)
    @inbounds for pid in 1:n_precursors
        offset = 8 * (pid - 1)
        for frag_idx in getPrecFragRange(frag_lookup, pid)
            frag = fragments[Int(frag_idx)]
            isIso(frag) && continue
            rank = Int(getRank(frag))
            1 <= rank <= 8 || continue
            keys[offset + rank] = _mbr_fragment_annotation_key(frag)
        end
    end
    return _MBRFragmentAnnotationKeys(keys)
end

@inline function _mbr_fragment_annotation_key(
    fragment_keys::_MBRFragmentAnnotationKeys,
    pid::UInt32,
    rank::Integer,
)
    return fragment_keys.keys[8 * (Int(pid) - 1) + rank]
end

@inline function _mbr_smoothed_spectrum_hellinger_from_sqrt(
    recipient_sqrt::NTuple{8, Float32},
    donor_sqrt::NTuple{8, Float32},
    fragment_keys::_MBRFragmentAnnotationKeys,
    recipient_pid::UInt32,
    donor_pid::UInt32,
)
    (_mbr_smoothed_spectrum_sqrt_is_valid(recipient_sqrt) &&
     _mbr_smoothed_spectrum_sqrt_is_valid(donor_sqrt)) || return 1.0f0

    dist2 = 0.0f0
    donor_matched = UInt16(0)
    @inbounds for recipient_rank in 1:8
        recipient_key = _mbr_fragment_annotation_key(fragment_keys, recipient_pid, recipient_rank)
        recipient_key == UInt16(0) && continue
        donor_rank = 0
        for candidate_rank in 1:8
            if _mbr_fragment_annotation_key(fragment_keys, donor_pid, candidate_rank) == recipient_key
                donor_rank = candidate_rank
                break
            end
        end
        if donor_rank == 0
            dist2 += recipient_sqrt[recipient_rank] * recipient_sqrt[recipient_rank]
        else
            delta = recipient_sqrt[recipient_rank] - donor_sqrt[donor_rank]
            dist2 += delta * delta
            donor_matched |= UInt16(1 << (donor_rank - 1))
        end
    end
    @inbounds for donor_rank in 1:8
        _mbr_fragment_annotation_key(fragment_keys, donor_pid, donor_rank) == UInt16(0) && continue
        (donor_matched & UInt16(1 << (donor_rank - 1))) != UInt16(0) && continue
        dist2 += donor_sqrt[donor_rank] * donor_sqrt[donor_rank]
    end
    return sqrt(max(0.0f0, 0.5f0 * dist2))
end

# Columns the donor dict reads from each (main + pass1 sidecar) pair.
# Listed once so reader and consumer stay in sync. `is_decoy` is hardcoded
# false on donor entries (the field is unused downstream), so we don't read
# it from the file.
const _MBR_DONOR_COLS = (:precursor_idx, :trace_prob_prepass, :weight,
    :log2_intensity_explained, :irt_pred, :irt_obs, :log_by_ratio_m0,
    :n_scans, :ms_file_idx, SMOOTHED_FRAGMENT_INTENSITY_COLUMNS...)

# Columns the per-file MBR sidecar emits. precursor_idx + scan_idx are
# redundant with the main file (same positions) but kept as alignment
# check-keys (asserted equal at downstream join time).
const _MBR_SIDECAR_OUT_COLS = (:precursor_idx, :scan_idx,
    :MBR_max_pair_prob_true,        :MBR_max_pair_prob_false,
    :MBR_log2_weight_ratio_true,    :MBR_log2_weight_ratio_false,
    :MBR_log2_explained_ratio_true, :MBR_log2_explained_ratio_false,
    :MBR_abs_n_scans_diff_true,     :MBR_abs_n_scans_diff_false,
    :MBR_log2_n_scans_ratio_true,   :MBR_log2_n_scans_ratio_false,
    :MBR_best_irt_diff_true,        :MBR_best_irt_diff_false,
    :MBR_is_missing_true,           :MBR_is_missing_false,
    :MBR_log_by_diff_true,          :MBR_log_by_diff_false,
    :MBR_smoothed_frag_hellinger_true,
    :MBR_smoothed_frag_hellinger_false,
)

# Suffix conventions for the three sidecar types used by the streaming MBR
# pipeline. All sidecars carry (:precursor_idx, :scan_idx) as alignment keys.
const PASS1_SIDECAR_SUFFIX    = ".pass1_sidecar.arrow"
const RECOVERY_SIDECAR_SUFFIX = ".recovery_sidecar.arrow"

@inline function _mbr_main_pep_confidence(pep)
    p = Float64(pep)
    isfinite(p) || return 0.0f0
    return Float32(round(clamp(1.0 - p, 0.0, 1.0), digits = 6))
end

# Slice [offset+1 .. offset+n] out of the four Pass-1 columns and write the
# sidecar at `side_path`. Concrete `AbstractVector{T}` signatures force Julia
# to specialise on the eltype (UInt32 / Float32 — see caller for the
# rationale) and act as a schema guard. View-based, no per-column copies.
#
# `main_pid_first` / `main_scn_first` are the first-row pid/scan of the main
# file we're slicing for; we assert they match best_psms[offset+1, :] so
# producer/consumer drift fails loud.
function _emit_pass1_sidecar!(
    side_path::String,
    pid_v::AbstractVector{UInt32},
    scan_v::AbstractVector{UInt32},
    score_v::AbstractVector{Float32},
    infold_v::Union{AbstractVector{Float32}, Nothing},
    main_pid::AbstractVector{UInt32},
    main_scn::AbstractVector{UInt32},
    offset::Int,
    n::Int,
    file_label::String,
)
    # Alignment guard at the slice boundary (first row's pid+scan vs main file).
    (pid_v[offset + 1] == first(main_pid) &&
     scan_v[offset + 1] == first(main_scn)) || error(
        "write_pass1_score_sidecars!: row $(offset + 1) of best_psms " *
        "misaligned with $file_label"
    )
    rng = (offset + 1):(offset + n)
    side_df = DataFrame(
        precursor_idx      = view(pid_v,   rng),
        scan_idx           = view(scan_v,  rng),
        trace_prob_prepass = view(score_v, rng),
        # No infold column when match_between_runs is false (rare; keeps the
        # sidecar schema constant so downstream readers don't branch).
        trace_prob_infold  = infold_v === nothing ?
                              fill(NaN32, n) :
                              view(infold_v, rng),
    )
    writeArrow(side_path, side_df)
    return nothing
end

function _empty_mbr_rescue_candidate_slim_columns()
    base_cols = (
        precursor_idx = UInt32[],
        scan_idx = UInt32[],
        ms_file_idx = UInt32[],
        cv_fold = UInt8[],
        target = Bool[],
        decoy = Bool[],
        mbr_rescue_candidate = Bool[],
        trace_prob_prepass = Float32[],
        main_search_prob = Float32[],
        trace_prob_infold = Float32[],
        _mbr_normal_file_idx = UInt32[],
        _mbr_normal_row_idx = UInt32[],
        _mbr_rescue_file_idx = UInt32[],
        _mbr_rescue_row_idx = UInt32[],
        MBR_max_pair_prob_true = Float32[],
        MBR_max_pair_prob_false = Float32[],
        MBR_log2_weight_ratio_true = Float32[],
        MBR_log2_weight_ratio_false = Float32[],
        MBR_log2_explained_ratio_true = Float32[],
        MBR_log2_explained_ratio_false = Float32[],
        MBR_abs_n_scans_diff_true = Float32[],
        MBR_abs_n_scans_diff_false = Float32[],
        MBR_log2_n_scans_ratio_true = Float32[],
        MBR_log2_n_scans_ratio_false = Float32[],
        MBR_best_irt_diff_true = Float32[],
        MBR_best_irt_diff_false = Float32[],
        MBR_is_missing_true = Bool[],
        MBR_is_missing_false = Bool[],
        MBR_log_by_diff_true = Float32[],
        MBR_log_by_diff_false = Float32[],
        MBR_smoothed_frag_hellinger_true = Float32[],
        MBR_smoothed_frag_hellinger_false = Float32[],
    )
    pass1_cols = NamedTuple{Tuple(ADVANCED_FEATURE_SET)}(
        ntuple(_ -> Float32[], length(ADVANCED_FEATURE_SET))
    )
    return merge(base_cols, pass1_cols)
end

function _mbr_candidate_columns_to_dataframe(cols)
    return DataFrame(cols; copycols = false)
end

function _empty_mbr_rescue_candidate_slim_dataframe()
    return _mbr_candidate_columns_to_dataframe(_empty_mbr_rescue_candidate_slim_columns())
end

function _mbr_pass1_feature_columns(main)
    return NamedTuple{Tuple(ADVANCED_FEATURE_SET)}(
        ntuple(i -> begin
            feature = ADVANCED_FEATURE_SET[i]
            hasproperty(main, feature) ? getproperty(main, feature) : nothing
        end, length(ADVANCED_FEATURE_SET))
    )
end

function _append_mbr_candidate_pass1_features!(out, pass1_feature_cols, row_idx::Int)
    for feature in ADVANCED_FEATURE_SET
        col = getproperty(pass1_feature_cols, feature)
        value = col === nothing ? 0.0f0 : Float32(col[row_idx])
        push!(getproperty(out, feature), value)
    end
    return nothing
end

function _append_mbr_candidate_row!(
    out,
    normal_file_idx::UInt32,
    normal_row_idx::UInt32,
    rescue_file_idx::UInt32,
    rescue_row_idx::UInt32,
    is_rescue::Bool,
    pid::UInt32,
    scan_idx::UInt32,
    ms_file_idx::UInt32,
    cv_fold::UInt8,
    target::Bool,
    prepass_score::Float32,
    main_search_prob::Float32,
    cross_run_prob::Float32,
    weight::Float32,
    log2_intensity_explained::Float32,
    irt_residual::Float32,
    log_by_ratio::Float32,
    n_scans::Float32,
    recipient_sqrt::NTuple{8, Float32},
    donor_t::_MBRDonorEntry,
    donor_f::_MBRDonorEntry,
    fragment_keys::_MBRFragmentAnnotationKeys,
)
    log2_weight_ratio_t = -1.0f0
    log2_weight_ratio_f = -1.0f0
    if weight > 0.0f0
        donor_t.weight > 0.0f0 && (log2_weight_ratio_t = log2(weight / donor_t.weight))
        donor_f.weight > 0.0f0 && (log2_weight_ratio_f = log2(weight / donor_f.weight))
    end

    hellinger_t = _mbr_smoothed_spectrum_hellinger_from_sqrt(
        recipient_sqrt,
        donor_t.smoothed_frag_sqrt,
        fragment_keys,
        pid,
        donor_t.precursor_idx,
    )
    hellinger_f = _mbr_smoothed_spectrum_hellinger_from_sqrt(
        recipient_sqrt,
        donor_f.smoothed_frag_sqrt,
        fragment_keys,
        pid,
        donor_f.precursor_idx,
    )

    push!(out.precursor_idx, pid)
    push!(out.scan_idx, scan_idx)
    push!(out.ms_file_idx, ms_file_idx)
    push!(out.cv_fold, cv_fold)
    push!(out.target, target)
    push!(out.decoy, !target)
    push!(out.mbr_rescue_candidate, is_rescue)
    push!(out.trace_prob_prepass, prepass_score)
    push!(out.main_search_prob, main_search_prob)
    push!(out.trace_prob_infold, cross_run_prob)
    push!(out._mbr_normal_file_idx, normal_file_idx)
    push!(out._mbr_normal_row_idx, normal_row_idx)
    push!(out._mbr_rescue_file_idx, rescue_file_idx)
    push!(out._mbr_rescue_row_idx, rescue_row_idx)
    push!(out.MBR_max_pair_prob_true, donor_t.trace_prob)
    push!(out.MBR_max_pair_prob_false, donor_f.trace_prob)
    push!(out.MBR_log2_weight_ratio_true, log2_weight_ratio_t)
    push!(out.MBR_log2_weight_ratio_false, log2_weight_ratio_f)
    push!(out.MBR_log2_explained_ratio_true, log2_intensity_explained - donor_t.log2_intensity_explained)
    push!(out.MBR_log2_explained_ratio_false, log2_intensity_explained - donor_f.log2_intensity_explained)
    push!(out.MBR_abs_n_scans_diff_true, abs(n_scans - donor_t.n_scans))
    push!(out.MBR_abs_n_scans_diff_false, abs(n_scans - donor_f.n_scans))
    push!(out.MBR_log2_n_scans_ratio_true, log2((n_scans + 1.0f0) / (donor_t.n_scans + 1.0f0)))
    push!(out.MBR_log2_n_scans_ratio_false, log2((n_scans + 1.0f0) / (donor_f.n_scans + 1.0f0)))
    push!(out.MBR_best_irt_diff_true, abs(irt_residual - donor_t.irt_residual))
    push!(out.MBR_best_irt_diff_false, abs(irt_residual - donor_f.irt_residual))
    push!(out.MBR_is_missing_true, false)
    push!(out.MBR_is_missing_false, false)
    push!(out.MBR_log_by_diff_true, log_by_ratio - donor_t.log_by_ratio)
    push!(out.MBR_log_by_diff_false, log_by_ratio - donor_f.log_by_ratio)
    push!(out.MBR_smoothed_frag_hellinger_true, hellinger_t)
    push!(out.MBR_smoothed_frag_hellinger_false, hellinger_f)
    return nothing
end

# Distribute Pass-1 LightGBM scores (trace_prob_prepass, trace_prob_infold)
# back to per-file Pass-1 sidecars.
#
# Invariant we rely on (no groupby needed): best_psms was built by
# load_psms_for_lightgbm via `Arrow.Table(file_paths)` over the readdir-
# sorted contents of second_pass_folder, then concatenated by
# `Tables.columntable + DataFrame`. That preserves row order, so rows
# 1..n1 belong to file 1, n1+1..n1+n2 belong to file 2, etc. (alphabetical
# readdir order). We sort `file_paths` the same way to match that layout,
# walk linearly with a cumulative offset, and slice via `view` — no
# groupby, no per-column copies, no Type.(vector) allocations.
#
# Columns we read on best_psms are already in their target eltypes
# (precursor_idx/scan_idx UInt32 from Arrow; trace_prob_prepass/_infold
# Float32 from score_psms.jl), so no narrowing conversion is needed.
# An alignment guard at every file boundary (first-row pid+scan_idx
# against the main file) fails loud if the producer/consumer ever drift
# apart (e.g. someone passes file_paths in a different order or filters
# rows mid-pipeline).
function write_pass1_score_sidecars!(best_psms::DataFrame, file_paths::Vector{String})
    # Match load_psms_for_lightgbm's row order: readdir is alphabetical,
    # so `sort(file_paths)` gives the same order without depending on the
    # caller's convention.
    sorted_paths = sort(file_paths)

    # Pull DataFrame columns once and hand them to a typed helper. The helper's
    # `AbstractVector{T}` signature triggers compiler specialisation on the
    # concrete eltype and acts as a schema guard — bails immediately if a
    # column type ever drifts. (Arrow.Primitive{UInt32, ...} and Vector{UInt32}
    # both satisfy AbstractVector{UInt32}, so no copy is forced.)
    pid_v      = best_psms.precursor_idx
    scan_v     = best_psms.scan_idx
    score_v    = best_psms.trace_prob_prepass
    has_infold = hasproperty(best_psms, :trace_prob_infold)
    infold_v   = has_infold ? best_psms.trace_prob_infold : nothing

    n_total   = nrow(best_psms)
    offset    = 0
    n_written = 0
    for fpath in sorted_paths
        main = Arrow.Table(fpath)
        n = length(main.precursor_idx)
        n == 0 && continue
        offset + n <= n_total || error(
            "write_pass1_score_sidecars!: best_psms has $n_total rows but " *
            "cumulative file row count exceeds it at $(basename(fpath)) " *
            "(offset=$offset, this file n=$n). Likely best_psms / file_paths drift."
        )
        _emit_pass1_sidecar!(fpath * PASS1_SIDECAR_SUFFIX,
                             pid_v, scan_v, score_v, infold_v,
                             main.precursor_idx, main.scan_idx,
                             offset, n, basename(fpath))
        offset    += n
        n_written += 1
    end
    offset == n_total || error(
        "write_pass1_score_sidecars!: distributed $offset rows but best_psms has $n_total " *
        "(some file_paths missing from second_pass_folder?)"
    )
    @debug_l1 "  Wrote $n_written Pass-1 score sidecars (positional, no groupby)"
    return n_written
end

function _insert_sorted_donor_entry!(
    entries::Vector{_MBRDonorEntry},
    entry::_MBRDonorEntry,
)
    insert_idx = length(entries) + 1
    @inbounds for idx in eachindex(entries)
        if entry.trace_prob > entries[idx].trace_prob
            insert_idx = idx
            break
        end
    end
    insert!(entries, insert_idx, entry)
    return nothing
end

# Inner per-file row loop for build_mbr_donor_dict_streaming_with_pass1.
# Extracted so Julia specializes on the concrete Arrow column types and the
# Union{Nothing, ...} branch on logby_c becomes a one-time call-site
# decision instead of per-row dispatch.
@inline function _accumulate_donor_entries!(
    all_entries::Dict{UInt32, Vector{_MBRDonorEntry}},
    pid_c::AbstractVector{UInt32},
    side_pid::AbstractVector{UInt32},
    main_scn::AbstractVector{UInt32},
    side_scn::AbstractVector{UInt32},
    score_c::AbstractVector{Float32},
    w_c::AbstractVector{Float32},
    l2ie_c::AbstractVector{Float16},
    irtp_c::AbstractVector{Float32},
    irto_c::AbstractVector{Float32},
    logby_c::Union{Nothing, AbstractVector{Float16}},
    nscans_c::AbstractVector,
    smoothed_frag_cols,
    fidx_c::AbstractVector{UInt32},
    side_path::String,
)
    n = length(pid_c)
    has_logby = logby_c !== nothing
    @inbounds for i in 1:n
        # Sanity-check alignment (cheap: one cmp per row)
        (pid_c[i] == side_pid[i] && main_scn[i] == side_scn[i]) ||
            error("Pass-1 sidecar misaligned at row $i of $side_path")
        pid = pid_c[i]
        e = _MBRDonorEntry(
            score_c[i], pid, w_c[i], Float32(l2ie_c[i]),
            irtp_c[i] - irto_c[i],
            irto_c[i],
            has_logby ? Float32(logby_c[i]) : 0f0,
            Float32(nscans_c[i]),
            _mbr_smoothed_spectrum_sqrt_tuple(smoothed_frag_cols, i),
            fidx_c[i], false,
        )
        entries = get!(() -> _MBRDonorEntry[], all_entries, pid)
        same_file_idx = 0
        for idx in eachindex(entries)
            if entries[idx].ms_file_idx == e.ms_file_idx
                same_file_idx = idx
                break
            end
        end
        if same_file_idx != 0
            e.trace_prob > entries[same_file_idx].trace_prob || continue
            deleteat!(entries, same_file_idx)
        end
        _insert_sorted_donor_entry!(entries, e)
        length(entries) > MBR_MAX_DONOR_FILES_PER_PRECURSOR && pop!(entries)
    end
    return nothing
end

# Streaming version: like build_mbr_donor_dict_streaming, but reads
# `trace_prob_prepass` from the Pass-1 sidecar (rather than expecting it
# in the main file). All other donor columns come from the main file.
# Alignment between main and sidecar is asserted via (:precursor_idx, :scan_idx).
function build_mbr_donor_dict_streaming_with_pass1(file_paths::Vector{String})
    all_entries = Dict{UInt32, Vector{_MBRDonorEntry}}()
    for main_path in file_paths
        side_path = main_path * PASS1_SIDECAR_SUFFIX
        isfile(side_path) || error("Missing Pass-1 sidecar at $side_path")
        main = Arrow.Table(main_path)
        side = Arrow.Table(side_path)
        n = length(main.precursor_idx)
        n == length(side.precursor_idx) ||
            error("Pass-1 sidecar row count mismatch at $side_path")
        _accumulate_donor_entries!(
            all_entries,
            main.precursor_idx, side.precursor_idx,
            main.scan_idx, side.scan_idx,
            side.trace_prob_prepass,
            main.weight,
            main.log2_intensity_explained,
            main.irt_pred, main.irt_obs,
            hasproperty(main, :log_by_ratio_m0) ? main.log_by_ratio_m0 : nothing,
            main.n_scans,
            ntuple(rank -> getproperty(main, SMOOTHED_FRAGMENT_INTENSITY_COLUMNS[rank]), 8),
            main.ms_file_idx,
            side_path,
        )
    end
    return all_entries
end

# Find the donor entry for `pid` from a file OTHER than `my_file`. Pulled
# out of the inner loop so Julia specializes on the donor_dict value type.
@inline function _donor_for_pid(donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
                                 pid::UInt32, my_file::UInt32)
    entries = get(donor_dict, pid, nothing)
    entries === nothing && return nothing
    @inbounds for e in entries
        e.ms_file_idx != my_file && return e
    end
    return nothing
end

@inline function _nearest_cross_file_donor(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    pool::_IrtPool,
    target_irt::Float32,
    my_file::UInt32,
    my_pid::UInt32,
)
    n = length(pool.irts)
    n == 0 && return nothing

    right = searchsortedfirst(pool.irts, target_irt)
    left = right - 1
    @inbounds while left >= 1 || right <= n
        use_left = if right > n
            true
        elseif left < 1
            false
        else
            abs(target_irt - pool.irts[left]) <= abs(pool.irts[right] - target_irt)
        end
        pool_idx = use_left ? left : right
        use_left ? (left -= 1) : (right += 1)
        pool.pids[pool_idx] == my_pid && continue
        donor = _donor_for_pid(donor_dict, pool.pids[pool_idx], my_file)
        donor !== nothing && return donor
    end
    return nothing
end

struct _MBRDonorRangeIndex
    irts::Vector{Float32}
    pid_to_index::Dict{UInt32, Int}
    leaf_count::Int
    best1::Vector{Union{Nothing, _MBRDonorEntry}}
    best2::Vector{Union{Nothing, _MBRDonorEntry}}
end

struct _MBRHardCounterfactualIndexes
    pools::Dict{Tuple{Int, Int}, _MBRDonorRangeIndex}
end

@inline function _mbr_donor_is_better(a::_MBRDonorEntry, b)
    b === nothing && return true
    return a.trace_prob > b.trace_prob ||
           (a.trace_prob == b.trace_prob && a.precursor_idx < b.precursor_idx)
end

@inline _mbr_donor_is_better(::Nothing, _) = false

@inline function _mbr_insert_top2(
    best1::Union{Nothing, _MBRDonorEntry},
    best2::Union{Nothing, _MBRDonorEntry},
    donor::_MBRDonorEntry,
)
    if best1 !== nothing && donor.ms_file_idx == best1.ms_file_idx
        _mbr_donor_is_better(donor, best1) && (best1 = donor)
    elseif best2 !== nothing && donor.ms_file_idx == best2.ms_file_idx
        _mbr_donor_is_better(donor, best2) && (best2 = donor)
    elseif _mbr_donor_is_better(donor, best1)
        best2 = best1
        best1 = donor
    elseif _mbr_donor_is_better(donor, best2)
        best2 = donor
    end
    if best1 !== nothing && best2 !== nothing && _mbr_donor_is_better(best2, best1)
        best1, best2 = best2, best1
    end
    return best1, best2
end

@inline function _mbr_merge_top2!(
    best1::Vector{Union{Nothing, _MBRDonorEntry}},
    best2::Vector{Union{Nothing, _MBRDonorEntry}},
    node::Int,
    donor,
)
    donor === nothing && return nothing
    best1[node], best2[node] = _mbr_insert_top2(best1[node], best2[node], donor)
    return nothing
end

function _build_mbr_donor_range_index(
    pool::_IrtPool,
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
)
    n = length(pool.pids)
    leaf_count = 1
    while leaf_count < max(n, 1)
        leaf_count <<= 1
    end
    best1 = Vector{Union{Nothing, _MBRDonorEntry}}(undef, 2 * leaf_count)
    best2 = Vector{Union{Nothing, _MBRDonorEntry}}(undef, 2 * leaf_count)
    fill!(best1, nothing)
    fill!(best2, nothing)
    pid_to_index = Dict{UInt32, Int}()

    @inbounds for pool_idx in 1:n
        node = leaf_count + pool_idx - 1
        pid_to_index[pool.pids[pool_idx]] = pool_idx
        entries = get(donor_dict, pool.pids[pool_idx], nothing)
        entries === nothing && continue
        for donor in entries
            _mbr_merge_top2!(best1, best2, node, donor)
        end
    end

    @inbounds for node in (leaf_count - 1):-1:1
        left = node << 1
        right = left + 1
        _mbr_merge_top2!(best1, best2, node, best1[left])
        _mbr_merge_top2!(best1, best2, node, best2[left])
        _mbr_merge_top2!(best1, best2, node, best1[right])
        _mbr_merge_top2!(best1, best2, node, best2[right])
    end

    return _MBRDonorRangeIndex(pool.irts, pid_to_index, leaf_count, best1, best2)
end

function build_mbr_hard_counterfactual_indexes(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
)
    pools = Dict{Tuple{Int, Int}, _MBRDonorRangeIndex}(
        key => _build_mbr_donor_range_index(pool, donor_dict)
        for (key, pool) in partner_pools.pools
    )
    return _MBRHardCounterfactualIndexes(pools)
end

@inline function _best_cross_file_from_range_node(
    index::_MBRDonorRangeIndex,
    node::Int,
    my_file::UInt32,
)
    donor = index.best1[node]
    donor !== nothing && donor.ms_file_idx != my_file && return donor
    donor = index.best2[node]
    donor !== nothing && donor.ms_file_idx != my_file && return donor
    return nothing
end

@inline function _best_scoring_cross_file_donor_in_range(
    index::_MBRDonorRangeIndex,
    first_idx::Int,
    last_idx::Int,
    my_file::UInt32,
)
    first_idx <= last_idx || return nothing

    left = index.leaf_count + first_idx - 1
    right = index.leaf_count + last_idx - 1
    best_donor = nothing
    @inbounds while left <= right
        if isodd(left)
            donor = _best_cross_file_from_range_node(index, left, my_file)
            donor !== nothing && _mbr_donor_is_better(donor, best_donor) && (best_donor = donor)
            left += 1
        end
        if iseven(right)
            donor = _best_cross_file_from_range_node(index, right, my_file)
            donor !== nothing && _mbr_donor_is_better(donor, best_donor) && (best_donor = donor)
            right -= 1
        end
        left >>= 1
        right >>= 1
    end
    return best_donor
end

@inline function _best_scoring_cross_file_donor_in_irt_window(
    index::_MBRDonorRangeIndex,
    target_irt::Float32,
    my_file::UInt32,
    my_pid::UInt32,
    irt_window::Float32,
)
    n = length(index.irts)
    n == 0 && return nothing

    first_idx = searchsortedfirst(index.irts, target_irt - irt_window)
    last_idx = searchsortedlast(index.irts, target_irt + irt_window)
    first_idx <= last_idx || return nothing

    self_idx = get(index.pid_to_index, my_pid, 0)
    if first_idx <= self_idx <= last_idx
        left_donor = _best_scoring_cross_file_donor_in_range(index, first_idx, self_idx - 1, my_file)
        right_donor = _best_scoring_cross_file_donor_in_range(index, self_idx + 1, last_idx, my_file)
        return _mbr_donor_is_better(right_donor, left_donor) ? right_donor : left_donor
    end
    return _best_scoring_cross_file_donor_in_range(index, first_idx, last_idx, my_file)
end

@inline function _best_scoring_cross_file_donor_in_irt_window(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    pool::_IrtPool,
    target_irt::Float32,
    my_file::UInt32,
    my_pid::UInt32,
    irt_window::Float32,
)
    n = length(pool.irts)
    n == 0 && return nothing

    first_idx = searchsortedfirst(pool.irts, target_irt - irt_window)
    last_idx = searchsortedlast(pool.irts, target_irt + irt_window)
    first_idx <= last_idx || return nothing

    best_donor = nothing
    best_score = -Inf32
    best_irt_delta = Inf32
    @inbounds for pool_idx in first_idx:last_idx
        donor_pid = pool.pids[pool_idx]
        donor_pid == my_pid && continue
        donor = _donor_for_pid(donor_dict, donor_pid, my_file)
        donor === nothing && continue
        score = donor.trace_prob
        irt_delta = abs(pool.irts[pool_idx] - target_irt)
        if best_donor === nothing ||
           score > best_score ||
           (score == best_score && irt_delta < best_irt_delta) ||
           (score == best_score && irt_delta == best_irt_delta &&
            donor.precursor_idx < best_donor.precursor_idx)
            best_donor = donor
            best_score = score
            best_irt_delta = irt_delta
        end
    end
    return best_donor
end

@inline function _hard_counterfactual_donor_from_pool(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    pool::_IrtPool,
    index::Union{Nothing, _MBRDonorRangeIndex},
    target_irt::Float32,
    my_file::UInt32,
    my_pid::UInt32,
)
    donor = index === nothing ?
            _best_scoring_cross_file_donor_in_irt_window(
                donor_dict,
                pool,
                target_irt,
                my_file,
                my_pid,
                MBR_COUNTERFACTUAL_LOCAL_IRT_WINDOW,
            ) :
            _best_scoring_cross_file_donor_in_irt_window(
                index,
                target_irt,
                my_file,
                my_pid,
                MBR_COUNTERFACTUAL_LOCAL_IRT_WINDOW,
            )
    donor !== nothing && return donor
    return _nearest_cross_file_donor(donor_dict, pool, target_irt, my_file, my_pid)
end

@inline function _resolve_false_donor_for_pid(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
    hard_indexes::Union{Nothing, _MBRHardCounterfactualIndexes},
    my_pid::UInt32,
    my_file::UInt32,
)
    pid_idx = Int(my_pid)
    my_irt = partner_pools.irt_by_pid[pid_idx]
    my_fold = Int(partner_pools.fold_by_pid[pid_idx])
    my_mz = Int(partner_pools.mz_bin_by_pid[pid_idx])
    pool_key = (my_fold, my_mz)

    donor = _hard_counterfactual_donor_from_pool(
        donor_dict,
        get(partner_pools.pools, pool_key, _empty_pool()),
        hard_indexes === nothing ? nothing : get(hard_indexes.pools, pool_key, nothing),
        my_irt,
        my_file,
        my_pid,
    )
    donor !== nothing && return donor
    donor = _nearest_cross_file_donor(donor_dict, partner_pools.fold_pool[my_fold], my_irt, my_file, my_pid)
    donor !== nothing && return donor
    return _nearest_cross_file_donor(donor_dict, partner_pools.global_pool, my_irt, my_file, my_pid)
end

@inline function _false_donor_for_pid(
    false_donor_cache::Dict{UInt32, Union{Nothing, _MBRDonorEntry}},
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
    hard_indexes::Union{Nothing, _MBRHardCounterfactualIndexes},
    my_pid::UInt32,
    my_file::UInt32,
)
    if haskey(false_donor_cache, my_pid)
        return false_donor_cache[my_pid]
    end
    donor = _resolve_false_donor_for_pid(donor_dict, partner_pools, hard_indexes, my_pid, my_file)
    false_donor_cache[my_pid] = donor
    return donor
end

function _normal_mbr_prepass_qvals_and_threshold(
    file_paths::Vector{String};
    donor_q_threshold::Float32 = MBR_DONOR_Q_THRESHOLD,
)
    scores = Float32[]
    targets = Bool[]
    row_counts = Int[]

    for main_path in file_paths
        pass1_path = main_path * PASS1_SIDECAR_SUFFIX
        isfile(pass1_path) || error("Missing Pass-1 sidecar at $pass1_path")
        main = Arrow.Table(main_path)
        pass1 = Arrow.Table(pass1_path)
        n = length(main.precursor_idx)
        n == length(pass1.precursor_idx) ||
            error("Pass-1 sidecar row count mismatch at $pass1_path")
        push!(row_counts, n)
        @inbounds for i in 1:n
            (main.precursor_idx[i] == pass1.precursor_idx[i] &&
             main.scan_idx[i] == pass1.scan_idx[i]) ||
                error("Pass-1 sidecar misaligned at row $i of $pass1_path")
            push!(scores, Float32(pass1.trace_prob_prepass[i]))
            push!(targets, Bool(main.target[i]))
        end
    end

    qvals = fill(1.0f0, length(scores))
    if !isempty(scores)
        get_qvalues!(scores, targets, qvals)
    end
    donor_target_pass = (qvals .<= donor_q_threshold) .& targets
    prob_thresh = any(donor_target_pass) ? minimum(scores[donor_target_pass]) : Float32(Inf)
    return (
        scores = scores,
        qvals = qvals,
        row_counts = row_counts,
        prob_thresh = prob_thresh,
    )
end

function load_normal_mbr_candidate_slim_dataframe(
    file_paths::Vector{String},
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
    fragment_keys::_MBRFragmentAnnotationKeys;
    q_thresh::Float32 = 0.01f0,
    donor_q_threshold::Float32 = MBR_DONOR_Q_THRESHOLD,
    hard_indexes::Union{Nothing, _MBRHardCounterfactualIndexes} = nothing,
)
    prepass = _normal_mbr_prepass_qvals_and_threshold(
        file_paths;
        donor_q_threshold = donor_q_threshold,
    )
    out = _empty_mbr_rescue_candidate_slim_columns()
    offset = 0

    for (file_idx, main_path) in pairs(file_paths)
        pass1_path = main_path * PASS1_SIDECAR_SUFFIX
        main = Arrow.Table(main_path)
        pass1 = Arrow.Table(pass1_path)
        n = prepass.row_counts[file_idx]
        n == 0 && continue

        pid_v = main.precursor_idx
        scan_v = main.scan_idx
        file_v = main.ms_file_idx
        fold_v = main.cv_fold
        target_v = main.target
        weight_v = main.weight
        l2ie_v = main.log2_intensity_explained
        irtp_v = main.irt_pred
        irto_v = main.irt_obs
        logby_v = hasproperty(main, :log_by_ratio_m0) ? main.log_by_ratio_m0 : nothing
        has_logby = logby_v !== nothing
        nscans_v = main.n_scans
        smoothed_frag_cols = ntuple(rank -> getproperty(main, SMOOTHED_FRAGMENT_INTENSITY_COLUMNS[rank]), 8)
        pass1_feature_cols = _mbr_pass1_feature_columns(main)
        false_donor_cache = Dict{UInt32, Union{Nothing, _MBRDonorEntry}}()

        @inbounds for i in 1:n
            global_row = offset + i
            prepass.qvals[global_row] > q_thresh || continue
            pid = UInt32(pid_v[i])
            ms_file_idx = UInt32(file_v[i])
            donor_t = _donor_for_pid(donor_dict, pid, ms_file_idx)
            donor_t === nothing && continue
            donor_t.trace_prob >= prepass.prob_thresh || continue
            donor_f = _false_donor_for_pid(false_donor_cache, donor_dict, partner_pools, hard_indexes, pid, ms_file_idx)
            donor_f === nothing && continue

            recipient_sqrt = _mbr_smoothed_spectrum_sqrt_tuple(smoothed_frag_cols, i)
            _append_mbr_candidate_row!(
                out,
                UInt32(file_idx),
                UInt32(i),
                UInt32(0),
                UInt32(0),
                false,
                pid,
                UInt32(scan_v[i]),
                ms_file_idx,
                UInt8(fold_v[i]),
                Bool(target_v[i]),
                Float32(pass1.trace_prob_prepass[i]),
                _mbr_main_pep_confidence(main.main_pep[i]),
                Float32(pass1.trace_prob_infold[i]),
                Float32(weight_v[i]),
                Float32(l2ie_v[i]),
                Float32(irtp_v[i]) - Float32(irto_v[i]),
                has_logby ? Float32(logby_v[i]) : 0.0f0,
                Float32(nscans_v[i]),
                recipient_sqrt,
                donor_t,
                donor_f,
                fragment_keys,
            )
            _append_mbr_candidate_pass1_features!(out, pass1_feature_cols, i)
        end
        offset += n
    end

    return (
        candidates = _mbr_candidate_columns_to_dataframe(out),
        prob_thresh = prepass.prob_thresh,
        n_rows = length(prepass.scores),
    )
end

function load_mbr_rescue_candidate_slim_dataframe(
    file_paths::Vector{String},
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
    fragment_keys::_MBRFragmentAnnotationKeys,
    prob_thresh::Float32;
    hard_indexes::Union{Nothing, _MBRHardCounterfactualIndexes} = nothing,
)
    out = _empty_mbr_rescue_candidate_slim_columns()
    for (file_idx, main_path) in pairs(file_paths)
        isfile(main_path) || continue
        main = Arrow.Table(main_path)
        hasproperty(main, :main_pep) ||
            error("MBR rescue file $main_path is missing required :main_pep column")
        n = length(main.precursor_idx)
        n == 0 && continue

        pid_v = main.precursor_idx
        scan_v = main.scan_idx
        file_v = main.ms_file_idx
        fold_v = main.cv_fold
        target_v = main.target
        pep_v = main.main_pep
        weight_v = main.weight
        l2ie_v = main.log2_intensity_explained
        irtp_v = main.irt_pred
        irto_v = main.irt_obs
        logby_v = hasproperty(main, :log_by_ratio_m0) ? main.log_by_ratio_m0 : nothing
        has_logby = logby_v !== nothing
        nscans_v = main.n_scans
        smoothed_frag_cols = ntuple(rank -> getproperty(main, SMOOTHED_FRAGMENT_INTENSITY_COLUMNS[rank]), 8)
        pass1_feature_cols = _mbr_pass1_feature_columns(main)
        false_donor_cache = Dict{UInt32, Union{Nothing, _MBRDonorEntry}}()

        @inbounds for i in 1:n
            pid = UInt32(pid_v[i])
            ms_file_idx = UInt32(file_v[i])
            donor_t = _donor_for_pid(donor_dict, pid, ms_file_idx)
            donor_t === nothing && continue
            donor_t.trace_prob >= prob_thresh || continue
            donor_f = _false_donor_for_pid(false_donor_cache, donor_dict, partner_pools, hard_indexes, pid, ms_file_idx)
            donor_f === nothing && continue

            recipient_sqrt = _mbr_smoothed_spectrum_sqrt_tuple(smoothed_frag_cols, i)
            main_search_prob = _mbr_main_pep_confidence(pep_v[i])
            _append_mbr_candidate_row!(
                out,
                UInt32(0),
                UInt32(0),
                UInt32(file_idx),
                UInt32(i),
                true,
                pid,
                UInt32(scan_v[i]),
                ms_file_idx,
                UInt8(fold_v[i]),
                Bool(target_v[i]),
                main_search_prob,
                main_search_prob,
                0.0f0,
                Float32(weight_v[i]),
                Float32(l2ie_v[i]),
                Float32(irtp_v[i]) - Float32(irto_v[i]),
                has_logby ? Float32(logby_v[i]) : 0.0f0,
                Float32(nscans_v[i]),
                recipient_sqrt,
                donor_t,
                donor_f,
                fragment_keys,
            )
            _append_mbr_candidate_pass1_features!(out, pass1_feature_cols, i)
        end
    end
    return _mbr_candidate_columns_to_dataframe(out)
end

function write_sparse_normal_recovery_sidecars!(slim_df::DataFrame, file_paths::Vector{String})
    rows_by_file = Dict{UInt32, Vector{Int}}()
    @inbounds for row in 1:nrow(slim_df)
        file_idx = UInt32(slim_df._mbr_normal_file_idx[row])
        file_idx == UInt32(0) && continue
        1 <= Int(file_idx) <= length(file_paths) ||
            error("write_sparse_normal_recovery_sidecars!: invalid normal file index $file_idx")
        push!(get!(() -> Int[], rows_by_file, file_idx), row)
    end

    n_files_written = 0
    n_rows_written = 0
    for file_idx in eachindex(file_paths)
        main_path = file_paths[file_idx]
        main = Arrow.Table(main_path)
        n = length(main.precursor_idx)
        n == 0 && continue

        mbr_recovered = falses(n)
        transfer_candidate = falses(n)
        target_decoy_prob = fill(NaN32, n)
        ftr_qval = fill(NaN32, n)
        ftr_pep = fill(NaN32, n)

        for row in get(rows_by_file, UInt32(file_idx), Int[])
            source_row = Int(slim_df._mbr_normal_row_idx[row])
            1 <= source_row <= n ||
                error("write_sparse_normal_recovery_sidecars!: invalid row $source_row for $main_path")
            (UInt32(main.precursor_idx[source_row]) == UInt32(slim_df.precursor_idx[row]) &&
             UInt32(main.scan_idx[source_row]) == UInt32(slim_df.scan_idx[row])) ||
                error("write_sparse_normal_recovery_sidecars!: slim_df misaligned at row $row for $main_path")
            mbr_recovered[source_row] = Bool(slim_df.mbr_recovered[row])
            transfer_candidate[source_row] = Bool(slim_df.MBR_transfer_candidate[row])
            target_decoy_prob[source_row] = Float32(slim_df.mbr_target_decoy_prob[row])
            ftr_qval[source_row] = Float32(slim_df.ftr_qval_true[row])
            ftr_pep[source_row] = Float32(slim_df.ftr_pep_true[row])
        end

        side_path = main_path * RECOVERY_SIDECAR_SUFFIX
        writeArrow(side_path, DataFrame(
            precursor_idx = collect(UInt32.(main.precursor_idx)),
            scan_idx = collect(UInt32.(main.scan_idx)),
            mbr_recovered = mbr_recovered,
            MBR_transfer_candidate = transfer_candidate,
            mbr_target_decoy_prob = target_decoy_prob,
            ftr_qval_true = ftr_qval,
            ftr_pep_true = ftr_pep,
        ))
        n_files_written += 1
        n_rows_written += n
    end
    return (n_files = n_files_written, n_rows = n_rows_written)
end

function write_mbr_rescue_recovered_files!(slim_df::DataFrame, file_paths::Vector{String})
    rows_by_file = Dict{UInt32, Vector{Int}}()
    n_rows = nrow(slim_df)
    if n_rows == 0
        for main_path in file_paths
            safeRm(mbr_rescue_recovered_path(main_path), nothing; force = true)
        end
        return (n_files = 0, n_rows = 0)
    end

    @inbounds for row in 1:n_rows
        Bool(slim_df.mbr_recovered[row]) || continue
        file_idx = UInt32(slim_df._mbr_rescue_file_idx[row])
        1 <= Int(file_idx) <= length(file_paths) ||
            error("write_mbr_rescue_recovered_files!: invalid rescue file index $file_idx")
        push!(get!(() -> Int[], rows_by_file, file_idx), row)
    end

    n_files_written = 0
    n_recovered_rows = 0
    for file_idx in eachindex(file_paths)
        main_path = file_paths[file_idx]
        rows = get(rows_by_file, UInt32(file_idx), Int[])
        if isempty(rows)
            safeRm(mbr_rescue_recovered_path(main_path), nothing; force = true)
            continue
        end
        main = DataFrame(Arrow.Table(main_path))
        source_rows = Int[slim_df._mbr_rescue_row_idx[row] for row in rows]
        recovered = main[source_rows, :]
        recovered[!, :trace_prob_prepass] = Float32[slim_df.trace_prob_prepass[row] for row in rows]
        recovered[!, :trace_prob_infold] = Float32[slim_df.trace_prob_infold[row] for row in rows]
        recovered[!, :trace_prob] = copy(recovered[!, :trace_prob_prepass])
        recovered[!, :mbr_recovered] = Bool[slim_df.mbr_recovered[row] for row in rows]
        recovered[!, :MBR_transfer_candidate] = Bool[slim_df.MBR_transfer_candidate[row] for row in rows]
        recovered[!, :mbr_target_decoy_prob] = Float32[slim_df.mbr_target_decoy_prob[row] for row in rows]
        recovered[!, :ftr_qval_true] = Float32[slim_df.ftr_qval_true[row] for row in rows]
        recovered[!, :ftr_pep_true] = Float32[slim_df.ftr_pep_true[row] for row in rows]

        writeArrow(mbr_rescue_recovered_path(main_path), recovered)
        n_files_written += 1
        n_recovered_rows += nrow(recovered)
    end
    return (n_files = n_files_written, n_rows = n_recovered_rows)
end

# Per-file merge step: read main + Pass-1 sidecar + recovery sidecar and write
# the downstream MBR output sidecar. MBR feature sidecars are intentionally not
# required here: FTR consumes sparse candidate rows directly, and the MBR
# feature columns are not needed after FTR scoring.
function merge_mbr_sidecars_into_main!(file_paths::Vector{String}; cleanup::Bool = true)
    n_merged = 0
    for path in file_paths
        pass1_path = path * PASS1_SIDECAR_SUFFIX
        rec_path   = path * RECOVERY_SIDECAR_SUFFIX
        all(isfile, (pass1_path, rec_path)) || continue

        main  = Arrow.Table(path)
        pass1 = Arrow.Table(pass1_path)
        rec   = Arrow.Table(rec_path)
        n = length(main.precursor_idx)
        (length(pass1.precursor_idx) == n &&
         length(rec.precursor_idx)   == n) ||
            error("Sidecar row-count mismatch at $path")
        @inbounds for i in 1:n
            (main.precursor_idx[i] == pass1.precursor_idx[i] == rec.precursor_idx[i] &&
             main.scan_idx[i]      == pass1.scan_idx[i]      == rec.scan_idx[i]) ||
                error("Sidecar misaligned at row $i of $path")
        end

        # Write the FTR/Pass-1 outputs as a single narrow sidecar tagged
        # "mbr_outputs" rather than rewriting the main file. Downstream
        # PSMFileReference construction auto-discovers and joins it.
        # Note: :decoy is already populated by MainSearch (features.jl:108)
        # — including it here would trigger the auto-discovery collision
        # check and reject the entire sidecar.
        side_df = DataFrame()
        side_df[!, :trace_prob_prepass]   = collect(Float32.(pass1.trace_prob_prepass))
        side_df[!, :trace_prob_infold]    = collect(Float32.(pass1.trace_prob_infold))
        side_df[!, :trace_prob]           = side_df[!, :trace_prob_prepass]
        side_df[!, :mbr_recovered]        = collect(Bool.(rec.mbr_recovered))
        side_df[!, :MBR_transfer_candidate] = collect(Bool.(rec.MBR_transfer_candidate))
        side_df[!, :mbr_target_decoy_prob] = collect(Float32.(rec.mbr_target_decoy_prob))
        side_df[!, :ftr_qval_true]        = collect(Float32.(rec.ftr_qval_true))
        side_df[!, :ftr_pep_true]         = collect(Float32.(rec.ftr_pep_true))

        new_sidecar_path = path * ".mbr_outputs.sidecar.arrow"
        writeArrow(new_sidecar_path, side_df)

        if cleanup
            # Release the Arrow mmap handles before deleting: on Windows a raw
            # rm() of a memory-mapped file throws EACCES until the mapping is
            # released. Mirror the safeRm pattern used elsewhere in this module.
            main = nothing; pass1 = nothing; rec = nothing
            GC.gc(false)
            safeRm(pass1_path, nothing; force=true)
            safeRm(rec_path, nothing; force=true)
            safeRm(path * ".mbr_sidecar.arrow", nothing; force=true)
        end
        n_merged += 1
    end
    @debug_l1 "  Wrote $n_merged consolidated mbr_outputs sidecars"
    return n_merged
end
