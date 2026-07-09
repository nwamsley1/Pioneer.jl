# MBR streaming pipeline. Replaces an older in-memory chain that held a single
# big best_psms DataFrame across MBR feature compute. We instead emit three
# per-file sidecar Arrow files and merge them back into the main file at the
# end (see score_psms.jl::score_precursor_isotope_traces for the orchestrator).
#
# Three sidecar types per main file `<f>.arrow`:
#   <f>.pass1_sidecar.arrow      Pass-1 LightGBM scores (trace_prob_prepass +
#                                trace_prob_infold). Written first so best_psms
#                                can be freed before the heavy MBR stage.
#   <f>.mbr_sidecar.arrow        Per-row MBR_* features (the _true/_false donor
#                                comparisons fed to the FTR LightGBM).
#   <f>.recovery_sidecar.arrow   FTR controller output (mbr_recovered, ftr_qval,
#                                ftr_pep, MBR_transfer_candidate).
#
# Two main-file sweeps drive MBR feature compute:
#   Sweep 1 (build_mbr_donor_dict_streaming_with_pass1)
#     Streams each (main + pass1 sidecar) pair to build the donor dict keyed
#     by precursor_idx. Each pid keeps the two best donor rows from distinct
#     files, sorted by score, so `_donor_for_pid` can choose the best cross-file
#     donor without storage growing with file count.
#
#   Sweep 2 (compute_mbr_features_per_file_to_sidecar_with_pass1!)
#     Loads ONE file's main + pass1 sidecar at a time, computes the MBR_*
#     features using `donor_dict`, writes the per-file MBR sidecar.
#     Main file is never rewritten during sweeps.
#
# Alignment invariant: row N of every sidecar ↔ row N of the main file
# (positional join). `:precursor_idx` and `:scan_idx` are written into each
# sidecar redundantly and asserted equal at every read/merge, so any sort or
# filter that breaks alignment fails loud rather than silently corrupting data.
#
# merge_mbr_sidecars_into_main! at the end folds all three sidecars back into
# the main file and deletes them, restoring the per-file Arrow as the single
# source of truth for downstream stages.

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
    rt_obs::Float32        # literal scan RT (minutes) of the donor row
    n_scans::Float32
    smoothed_frag_sqrt::NTuple{8, Float32}
    library_hellinger::Float32
    ms_file_idx::UInt32
    is_decoy::Bool
end

const MBR_SIDECAR_SUFFIX = ".mbr_sidecar.arrow"
const MBR_MAX_DONOR_FILES_PER_PRECURSOR = 2
const MBR_MAX_WORST_DONOR_FILES_PER_PRECURSOR = 2
const MBR_DEFAULT_N_COUNTERFACTUALS = 2
const MBR_MAX_COUNTERFACTUALS = 8
const MBR_LOD_WEIGHT_QUANTILE = Float32(0.05)
const MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT = ntuple(_ -> 0.0f0, 8)
const MBR_ANY_DONOR_FILE = typemax(UInt32)
const _MBRFalseDonorTuple = NTuple{MBR_MAX_COUNTERFACTUALS, Union{Nothing, _MBRDonorEntry}}
const _MBRFalseDonorCacheKey = Tuple{UInt32, UInt32, Float32, UInt32}

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
    _fragment_keys::_MBRFragmentAnnotationKeys,
    _recipient_pid::UInt32,
    _donor_pid::UInt32,
)
    (_mbr_smoothed_spectrum_sqrt_is_valid(recipient_sqrt) &&
     _mbr_smoothed_spectrum_sqrt_is_valid(donor_sqrt)) || return 1.0f0

    dist2 = 0.0f0
    @inbounds for rank in 1:8
        delta = recipient_sqrt[rank] - donor_sqrt[rank]
        dist2 += delta * delta
    end
    return sqrt(max(0.0f0, 0.5f0 * dist2))
end

# Columns the donor dict reads from each (main + pass1 sidecar) pair.
# Listed once so reader and consumer stay in sync. `is_decoy` is hardcoded
# false on donor entries (the field is unused downstream), so we don't read
# it from the file.
const _MBR_DONOR_COLS = (:precursor_idx, :trace_prob_prepass, :weight,
    :log2_intensity_explained, :irt_pred, :irt_obs, :log_by_ratio_m0, :rt,
    :n_scans, :ms_file_idx, :smoothed_2d_shadow_hellinger,
    SMOOTHED_FRAGMENT_INTENSITY_COLUMNS...)

# Columns the per-file MBR sidecar emits. precursor_idx + scan_idx are
# redundant with the main file (same positions) but kept as alignment
# check-keys (asserted equal at downstream join time).
@inline _mbr_counterfactual_suffix(counterfactual_idx::Int) =
    counterfactual_idx == 1 ? "_false" : "_false$(counterfactual_idx)"

@inline _mbr_counterfactual_column(base::AbstractString, counterfactual_idx::Int) =
    Symbol(base * _mbr_counterfactual_suffix(counterfactual_idx))

function _mbr_true_false_sidecar_columns(base::AbstractString)
    cols = Symbol[Symbol(base * "_true")]
    for counterfactual_idx in 1:MBR_MAX_COUNTERFACTUALS
        push!(cols, _mbr_counterfactual_column(base, counterfactual_idx))
    end
    return cols
end

function _mbr_sidecar_output_columns()
    cols = Symbol[:precursor_idx, :scan_idx]
    append!(cols, _mbr_true_false_sidecar_columns("MBR_best_pair_prob"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_worst_pair_prob"))
    push!(cols, :MBR_log2_weight_lod_ratio)
    append!(cols, _mbr_true_false_sidecar_columns("MBR_best_log2_weight_ratio"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_worst_log2_weight_ratio"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_best_log2_explained_ratio"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_worst_log2_explained_ratio"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_best_abs_n_scans_diff"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_worst_abs_n_scans_diff"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_best_log2_n_scans_ratio"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_worst_log2_n_scans_ratio"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_best_irt_diff"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_worst_irt_diff"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_best_observed_irt_diff"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_worst_observed_irt_diff"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_single_donor"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_best_is_missing"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_best_rt_diff"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_best_log_by_diff"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_best_smoothed_frag_hellinger"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_worst_smoothed_frag_hellinger"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_best_donor_library_hellinger"))
    append!(cols, _mbr_true_false_sidecar_columns("MBR_worst_donor_library_hellinger"))
    return Tuple(cols)
end

const _MBR_SIDECAR_OUT_COLS = _mbr_sidecar_output_columns()

# Suffix conventions for the three sidecar types used by the streaming MBR
# pipeline. All sidecars carry (:precursor_idx, :scan_idx) as alignment keys.
const PASS1_SIDECAR_SUFFIX    = ".pass1_sidecar.arrow"
const RECOVERY_SIDECAR_SUFFIX = ".recovery_sidecar.arrow"

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

@inline function _mbr_lower_weight_donor_first(a::_MBRDonorEntry, b::_MBRDonorEntry)
    a.weight < b.weight && return true
    a.weight > b.weight && return false
    a.trace_prob < b.trace_prob && return true
    a.trace_prob > b.trace_prob && return false
    return a.ms_file_idx < b.ms_file_idx
end

function _mbr_lod_log2_weight!(samples::Vector{Float32})
    isempty(samples) && return NaN32
    sort!(samples)
    idx = clamp(ceil(Int, Float64(MBR_LOD_WEIGHT_QUANTILE) * length(samples)), 1, length(samples))
    return samples[idx]
end

@inline function _mbr_log2_weight_lod_ratio(
    weight::Float32,
    ms_file_idx::UInt32,
    lod_log2_weight_by_file::Dict{UInt32, Float32},
    lod_log2_weight_global::Float32,
)
    weight > 0.0f0 || return -1.0f0
    lod_log2_weight = get(lod_log2_weight_by_file, ms_file_idx, lod_log2_weight_global)
    isfinite(lod_log2_weight) || return -1.0f0
    return log2(weight) - lod_log2_weight
end

@inline function _mbr_log2_weight_ratio(weight::Float32, donor::Union{Nothing, _MBRDonorEntry})
    (donor !== nothing && weight > 0.0f0 && donor.weight > 0.0f0) || return -1.0f0
    return log2(weight / donor.weight)
end

@inline function _mbr_log2_explained_ratio(
    log2_intensity_explained::Float32,
    donor::Union{Nothing, _MBRDonorEntry},
)
    donor === nothing && return -1.0f0
    return log2_intensity_explained - donor.log2_intensity_explained
end

@inline function _mbr_abs_n_scans_diff(
    n_scans::Float32,
    donor::Union{Nothing, _MBRDonorEntry},
)
    donor === nothing && return -1.0f0
    return abs(n_scans - donor.n_scans)
end

@inline function _mbr_log2_n_scans_ratio(
    n_scans::Float32,
    donor::Union{Nothing, _MBRDonorEntry},
)
    donor === nothing && return -1.0f0
    return log2((n_scans + 1.0f0) / (donor.n_scans + 1.0f0))
end

@inline function _mbr_smoothed_spectrum_hellinger_from_donor(
    recipient_sqrt::NTuple{8, Float32},
    donor::Union{Nothing, _MBRDonorEntry},
    fragment_keys::_MBRFragmentAnnotationKeys,
    recipient_pid::UInt32,
)
    donor === nothing && return 1.0f0
    return _mbr_smoothed_spectrum_hellinger_from_sqrt(
        recipient_sqrt,
        donor.smoothed_frag_sqrt,
        fragment_keys,
        recipient_pid,
        donor.precursor_idx,
    )
end

@inline function _mbr_donor_library_hellinger(donor::Union{Nothing, _MBRDonorEntry})
    donor === nothing && return 1.0f0
    return donor.library_hellinger
end

@inline function _mbr_same_donor(a::_MBRDonorEntry, b::_MBRDonorEntry)
    return a.precursor_idx == b.precursor_idx && a.ms_file_idx == b.ms_file_idx
end

# Inner per-file row loop for build_mbr_donor_dict_streaming_with_pass1.
# Extracted so Julia specializes on the concrete Arrow column types and the
# Union{Nothing, ...} branch on logby_c / rt_c becomes a one-time call-site
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
    rt_c::Union{Nothing, AbstractVector{Float32}},
    nscans_c::AbstractVector,
    library_hellinger_c::AbstractVector{Float32},
    smoothed_frag_cols,
    fidx_c::AbstractVector{UInt32},
    side_path::String,
    ;
    passing_score_floor::Float32 = Float32(Inf),
)
    n = length(pid_c)
    has_logby = logby_c !== nothing
    has_rt    = rt_c    !== nothing
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
            has_rt    ? rt_c[i]             : 0f0,
            Float32(nscans_c[i]),
            _mbr_smoothed_spectrum_sqrt_tuple(smoothed_frag_cols, i),
            Float32(library_hellinger_c[i]),
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
        length(entries) > MBR_MAX_DONOR_FILES_PER_PRECURSOR &&
            _prune_donor_entries!(entries, passing_score_floor)
    end
    return nothing
end

# Streaming version: like build_mbr_donor_dict_streaming, but reads
# `trace_prob_prepass` from the Pass-1 sidecar (rather than expecting it
# in the main file). All other donor columns come from the main file.
# Alignment between main and sidecar is asserted via (:precursor_idx, :scan_idx).
function build_mbr_donor_dict_streaming_with_pass1(
    file_paths::Vector{String};
    passing_score_floor::Float32 = Float32(Inf),
)
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
            hasproperty(main, :rt) ? main.rt : nothing,
            main.n_scans,
            main.smoothed_2d_shadow_hellinger,
            ntuple(rank -> getproperty(main, SMOOTHED_FRAGMENT_INTENSITY_COLUMNS[rank]), 8),
            main.ms_file_idx,
            side_path,
            passing_score_floor = passing_score_floor,
        )
    end
    return all_entries
end

@inline function _min_weight_alternate_donor_for_pid(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    pid::UInt32,
    my_file::UInt32,
    best_donor::_MBRDonorEntry,
    passing_score_floor::Float32,
)
    entries = get(donor_dict, pid, nothing)
    entries === nothing && return nothing
    min_weight = nothing
    @inbounds for e in entries
        e.ms_file_idx == my_file && continue
        e.trace_prob >= passing_score_floor || continue
        _mbr_same_donor(e, best_donor) && continue
        if min_weight === nothing || _mbr_lower_weight_donor_first(e, min_weight)
            min_weight = e
        end
    end
    return min_weight
end

function _prune_donor_entries!(
    entries::Vector{_MBRDonorEntry},
    passing_score_floor::Float32,
)
    sort!(entries, by = entry -> entry.trace_prob, rev = true)
    keep = falses(length(entries))
    @inbounds for idx in 1:min(length(entries), MBR_MAX_DONOR_FILES_PER_PRECURSOR)
        keep[idx] = true
    end

    kept_low_weight = 0
    while kept_low_weight < MBR_MAX_WORST_DONOR_FILES_PER_PRECURSOR
        low_weight_idx = 0
        @inbounds for idx in eachindex(entries)
            keep[idx] && continue
            entries[idx].trace_prob >= passing_score_floor || continue
            if low_weight_idx == 0 ||
               _mbr_lower_weight_donor_first(entries[idx], entries[low_weight_idx])
                low_weight_idx = idx
            end
        end
        low_weight_idx == 0 && break
        keep[low_weight_idx] = true
        kept_low_weight += 1
    end

    write_idx = 1
    @inbounds for read_idx in eachindex(entries)
        if keep[read_idx]
            entries[write_idx] = entries[read_idx]
            write_idx += 1
        end
    end
    resize!(entries, write_idx - 1)
    return nothing
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

@inline function _donor_for_pid_in_file(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    pid::UInt32,
    donor_file::UInt32,
)
    entries = get(donor_dict, pid, nothing)
    entries === nothing && return nothing
    @inbounds for e in entries
        e.ms_file_idx == donor_file && return e
    end
    return nothing
end

@inline function _counterfactual_receiver_eligible(
    ::Nothing,
    receiver_file::UInt32,
    pid::UInt32,
)
    return true
end

@inline function _counterfactual_receiver_eligible(
    eligibility_by_file::_CounterfactualEligibilityByFile,
    receiver_file::UInt32,
    pid::UInt32,
)
    pid_int = Int(pid)
    pid_int in eligibility_by_file.dataset_passed || return false
    run_passed = get(eligibility_by_file.run_passed_by_file, receiver_file, nothing)
    run_passed === nothing && return true
    return !(pid_int in run_passed)
end

@inline function _empty_false_donors()::_MBRFalseDonorTuple
    return ntuple(_ -> nothing, MBR_MAX_COUNTERFACTUALS)
end

@inline function _mbr_false_donor_seen(
    donor::_MBRDonorEntry,
    donors::Vector{Union{Nothing, _MBRDonorEntry}},
    n_found::Int,
)
    @inbounds for i in 1:n_found
        existing = donors[i]
        existing !== nothing && donor.precursor_idx == existing.precursor_idx && return true
    end
    return false
end

@inline function _nearest_irt_cross_file_donors(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    pool::_IrtPool,
    target_irt::Float32,
    my_file::UInt32,
    my_pid::UInt32,
    eligibility_by_file::Union{Nothing, _CounterfactualEligibilityByFile},
)
    n = length(pool.irts)
    n == 0 && return _empty_false_donors()

    right = searchsortedfirst(pool.irts, target_irt)
    left = right - 1
    donors = Union{Nothing, _MBRDonorEntry}[nothing for _ in 1:MBR_MAX_COUNTERFACTUALS]
    n_found = 0
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
        pool_pid = pool.pids[pool_idx]
        pool_pid == my_pid && continue
        _counterfactual_receiver_eligible(eligibility_by_file, my_file, pool_pid) || continue
        donor = _donor_for_pid(donor_dict, pool_pid, my_file)
        donor === nothing && continue
        _mbr_false_donor_seen(donor, donors, n_found) && continue
        n_found += 1
        donors[n_found] = donor
        n_found == MBR_MAX_COUNTERFACTUALS &&
            return ntuple(idx -> donors[idx], MBR_MAX_COUNTERFACTUALS)
    end
    return ntuple(idx -> donors[idx], MBR_MAX_COUNTERFACTUALS)
end

@inline function _nearest_irt_same_file_donors(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    pool::_IrtPool,
    target_irt::Float32,
    my_file::UInt32,
    my_pid::UInt32,
    donor_file::UInt32,
    eligibility_by_file::Union{Nothing, _CounterfactualEligibilityByFile},
)
    n = length(pool.irts)
    n == 0 && return _empty_false_donors()

    right = searchsortedfirst(pool.irts, target_irt)
    left = right - 1
    donors = Union{Nothing, _MBRDonorEntry}[nothing for _ in 1:MBR_MAX_COUNTERFACTUALS]
    n_found = 0
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
        pool_pid = pool.pids[pool_idx]
        pool_pid == my_pid && continue
        _counterfactual_receiver_eligible(eligibility_by_file, my_file, pool_pid) || continue
        donor = _donor_for_pid_in_file(donor_dict, pool_pid, donor_file)
        donor === nothing && continue
        _mbr_false_donor_seen(donor, donors, n_found) && continue
        n_found += 1
        donors[n_found] = donor
        n_found == MBR_MAX_COUNTERFACTUALS &&
            return ntuple(idx -> donors[idx], MBR_MAX_COUNTERFACTUALS)
    end
    return ntuple(idx -> donors[idx], MBR_MAX_COUNTERFACTUALS)
end

@inline function _resolve_false_donors_for_pid(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
    my_pid::UInt32,
    my_file::UInt32,
    eligibility_by_file::Union{Nothing, _CounterfactualEligibilityByFile},
    target_irt::Float32 = NaN32,
)
    pid_idx = Int(my_pid)
    my_irt = isfinite(target_irt) ? target_irt : partner_pools.irt_by_pid[pid_idx]
    my_charge = Int(partner_pools.charge_by_pid[pid_idx])
    my_length = Int(partner_pools.length_by_pid[pid_idx])

    return _nearest_irt_cross_file_donors(
        donor_dict,
        get(partner_pools.charge_length_pool, (my_charge, my_length), _empty_irt_pool()),
        my_irt,
        my_file,
        my_pid,
        eligibility_by_file,
    )
end

@inline function _resolve_false_donors_for_pid(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
    my_pid::UInt32,
    my_file::UInt32,
    eligibility_by_file::Union{Nothing, _CounterfactualEligibilityByFile},
    donor_file::UInt32,
    target_irt::Float32 = NaN32,
)
    pid_idx = Int(my_pid)
    my_irt = isfinite(target_irt) ? target_irt : partner_pools.irt_by_pid[pid_idx]
    my_charge = Int(partner_pools.charge_by_pid[pid_idx])
    my_length = Int(partner_pools.length_by_pid[pid_idx])

    return _nearest_irt_same_file_donors(
        donor_dict,
        get(partner_pools.file_charge_length_pool, (donor_file, my_charge, my_length), _empty_irt_pool()),
        my_irt,
        my_file,
        my_pid,
        donor_file,
        eligibility_by_file,
    )
end

@inline function _resolve_false_donor_for_pid(
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
    my_pid::UInt32,
    my_file::UInt32,
    eligibility_by_file::Union{Nothing, _CounterfactualEligibilityByFile},
    target_irt::Float32 = NaN32,
)
    return _resolve_false_donors_for_pid(
        donor_dict,
        partner_pools,
        my_pid,
        my_file,
        eligibility_by_file,
        target_irt,
    )[1]
end

@inline function _false_donors_for_pid(
    false_donor_cache::Dict{_MBRFalseDonorCacheKey, _MBRFalseDonorTuple},
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
    my_pid::UInt32,
    my_file::UInt32,
    eligibility_by_file::Union{Nothing, _CounterfactualEligibilityByFile},
    target_irt::Float32 = NaN32,
)
    cache_irt = isfinite(target_irt) ? target_irt : partner_pools.irt_by_pid[Int(my_pid)]
    cache_key = (my_file, my_pid, cache_irt, MBR_ANY_DONOR_FILE)
    if haskey(false_donor_cache, cache_key)
        return false_donor_cache[cache_key]
    end
    donors = _resolve_false_donors_for_pid(
        donor_dict,
        partner_pools,
        my_pid,
        my_file,
        eligibility_by_file,
        target_irt,
    )
    false_donor_cache[cache_key] = donors
    return donors
end

@inline function _false_donors_for_pid(
    false_donor_cache::Dict{_MBRFalseDonorCacheKey, _MBRFalseDonorTuple},
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
    my_pid::UInt32,
    my_file::UInt32,
    eligibility_by_file::Union{Nothing, _CounterfactualEligibilityByFile},
    donor_file::UInt32,
    target_irt::Float32 = NaN32,
)
    cache_irt = isfinite(target_irt) ? target_irt : partner_pools.irt_by_pid[Int(my_pid)]
    cache_key = (my_file, my_pid, cache_irt, donor_file)
    if haskey(false_donor_cache, cache_key)
        return false_donor_cache[cache_key]
    end
    donors = _resolve_false_donors_for_pid(
        donor_dict,
        partner_pools,
        my_pid,
        my_file,
        eligibility_by_file,
        donor_file,
        target_irt,
    )
    false_donor_cache[cache_key] = donors
    return donors
end

@inline function _false_donor_for_pid(
    false_donor_cache::Dict{_MBRFalseDonorCacheKey, _MBRFalseDonorTuple},
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
    my_pid::UInt32,
    my_file::UInt32,
    eligibility_by_file::Union{Nothing, _CounterfactualEligibilityByFile},
    target_irt::Float32 = NaN32,
)
    return _false_donors_for_pid(
        false_donor_cache,
        donor_dict,
        partner_pools,
        my_pid,
        my_file,
        eligibility_by_file,
        target_irt,
    )[1]
end

@inline function _false_donor_for_pid(
    false_donor_cache::Dict{UInt32, Union{Nothing, _MBRDonorEntry}},
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
    my_pid::UInt32,
    my_file::UInt32,
)
    if haskey(false_donor_cache, my_pid)
        return false_donor_cache[my_pid]
    end
    donor = _resolve_false_donor_for_pid(
        donor_dict,
        partner_pools,
        my_pid,
        my_file,
        nothing,
    )
    false_donor_cache[my_pid] = donor
    return donor
end

function _mbr_prepass_donor_summary(
    file_paths::Vector{String};
    donor_q_threshold::Float32 = 0.01f0,
)
    scores = Float32[]
    targets = Bool[]
    weights = Float32[]
    file_indices = UInt32[]
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
            push!(weights, Float32(main.weight[i]))
            push!(file_indices, UInt32(main.ms_file_idx[i]))
        end
    end

    qvals = fill(1.0f0, length(scores))
    if !isempty(scores)
        get_qvalues!(scores, targets, qvals)
    end
    donor_target_pass = (qvals .<= donor_q_threshold) .& targets
    prob_thresh = any(donor_target_pass) ? minimum(scores[donor_target_pass]) : Float32(Inf)

    lod_samples_by_file = Dict{UInt32, Vector{Float32}}()
    lod_samples_global = Float32[]
    @inbounds for i in eachindex(scores)
        donor_target_pass[i] || continue
        weights[i] > 0.0f0 || continue
        log2_weight = log2(weights[i])
        push!(get!(() -> Float32[], lod_samples_by_file, file_indices[i]), log2_weight)
        push!(lod_samples_global, log2_weight)
    end
    lod_log2_weight_by_file = Dict{UInt32, Float32}()
    for (file_idx, samples) in lod_samples_by_file
        lod_log2_weight_by_file[file_idx] = _mbr_lod_log2_weight!(samples)
    end

    return (
        scores = scores,
        qvals = qvals,
        row_counts = row_counts,
        prob_thresh = prob_thresh,
        lod_log2_weight_by_file = lod_log2_weight_by_file,
        lod_log2_weight_global = _mbr_lod_log2_weight!(lod_samples_global),
    )
end

# Inner per-row MBR-feature compute. Extracted from
# compute_mbr_features_per_file_to_sidecar_with_pass1! so Julia specializes
# on the concrete column types (Arrow.Primitive{T} or Vector{T}) and the
# Union{Nothing, ...} branch on rt_v becomes a one-time call-site
# decision instead of per-row dispatch.
@inline function _compute_mbr_inner!(
    out_best_pair_t::Vector{Float32}, out_best_pair_f::Vector{Float32},
    out_best_pair_f2::Vector{Float32},
    out_best_pair_f3::Vector{Float32}, out_best_pair_f4::Vector{Float32},
    out_worst_pair_t::Vector{Float32}, out_worst_pair_f::Vector{Float32},
    out_worst_pair_f2::Vector{Float32},
    out_worst_pair_f3::Vector{Float32}, out_worst_pair_f4::Vector{Float32},
    out_lw_lod::Vector{Float32},
    out_lw_t::Vector{Float32},  out_lw_f::Vector{Float32},
    out_lw_f2::Vector{Float32},
    out_lw_f3::Vector{Float32}, out_lw_f4::Vector{Float32},
    out_lw_worst_t::Vector{Float32}, out_lw_worst_f::Vector{Float32},
    out_lw_worst_f2::Vector{Float32},
    out_lw_worst_f3::Vector{Float32}, out_lw_worst_f4::Vector{Float32},
    out_le_t::Vector{Float32},  out_le_f::Vector{Float32},
    out_le_f2::Vector{Float32},
    out_le_f3::Vector{Float32}, out_le_f4::Vector{Float32},
    out_le_worst_t::Vector{Float32}, out_le_worst_f::Vector{Float32},
    out_le_worst_f2::Vector{Float32},
    out_le_worst_f3::Vector{Float32}, out_le_worst_f4::Vector{Float32},
    out_nscans_t::Vector{Float32}, out_nscans_f::Vector{Float32},
    out_nscans_f2::Vector{Float32},
    out_nscans_f3::Vector{Float32}, out_nscans_f4::Vector{Float32},
    out_nscans_worst_t::Vector{Float32}, out_nscans_worst_f::Vector{Float32},
    out_nscans_worst_f2::Vector{Float32},
    out_nscans_worst_f3::Vector{Float32}, out_nscans_worst_f4::Vector{Float32},
    out_l2_nscans_t::Vector{Float32}, out_l2_nscans_f::Vector{Float32},
    out_l2_nscans_f2::Vector{Float32},
    out_l2_nscans_f3::Vector{Float32}, out_l2_nscans_f4::Vector{Float32},
    out_l2_nscans_worst_t::Vector{Float32}, out_l2_nscans_worst_f::Vector{Float32},
    out_l2_nscans_worst_f2::Vector{Float32},
    out_l2_nscans_worst_f3::Vector{Float32}, out_l2_nscans_worst_f4::Vector{Float32},
    out_ir_t::Vector{Float32},  out_ir_f::Vector{Float32},
    out_ir_f2::Vector{Float32},
    out_ir_f3::Vector{Float32}, out_ir_f4::Vector{Float32},
    out_ir_worst_t::Vector{Float32}, out_ir_worst_f::Vector{Float32},
    out_ir_worst_f2::Vector{Float32},
    out_ir_worst_f3::Vector{Float32}, out_ir_worst_f4::Vector{Float32},
    out_obs_ir_t::Vector{Float32}, out_obs_ir_f::Vector{Float32},
    out_obs_ir_f2::Vector{Float32},
    out_obs_ir_f3::Vector{Float32}, out_obs_ir_f4::Vector{Float32},
    out_obs_ir_worst_t::Vector{Float32}, out_obs_ir_worst_f::Vector{Float32},
    out_obs_ir_worst_f2::Vector{Float32},
    out_obs_ir_worst_f3::Vector{Float32}, out_obs_ir_worst_f4::Vector{Float32},
    out_single_donor_t::Vector{Float32}, out_single_donor_f::Vector{Float32},
    out_single_donor_f2::Vector{Float32},
    out_single_donor_f3::Vector{Float32}, out_single_donor_f4::Vector{Float32},
    out_miss_t::BitVector,      out_miss_f::BitVector,
    out_miss_f2::BitVector,
    out_miss_f3::BitVector, out_miss_f4::BitVector,
    out_rt_t::Vector{Float32},   out_rt_f::Vector{Float32},
    out_rt_f2::Vector{Float32},
    out_rt_f3::Vector{Float32}, out_rt_f4::Vector{Float32},
    out_log_by_t::Vector{Float32}, out_log_by_f::Vector{Float32},
    out_log_by_f2::Vector{Float32},
    out_log_by_f3::Vector{Float32}, out_log_by_f4::Vector{Float32},
    out_smoothed_hellinger_t::Vector{Float32},
    out_smoothed_hellinger_f::Vector{Float32},
    out_smoothed_hellinger_f2::Vector{Float32},
    out_smoothed_hellinger_f3::Vector{Float32},
    out_smoothed_hellinger_f4::Vector{Float32},
    out_smoothed_hellinger_worst_t::Vector{Float32},
    out_smoothed_hellinger_worst_f::Vector{Float32},
    out_smoothed_hellinger_worst_f2::Vector{Float32},
    out_smoothed_hellinger_worst_f3::Vector{Float32},
    out_smoothed_hellinger_worst_f4::Vector{Float32},
    out_donor_library_hellinger_t::Vector{Float32},
    out_donor_library_hellinger_f::Vector{Float32},
    out_donor_library_hellinger_f2::Vector{Float32},
    out_donor_library_hellinger_f3::Vector{Float32},
    out_donor_library_hellinger_f4::Vector{Float32},
    out_donor_library_hellinger_worst_t::Vector{Float32},
    out_donor_library_hellinger_worst_f::Vector{Float32},
    out_donor_library_hellinger_worst_f2::Vector{Float32},
    out_donor_library_hellinger_worst_f3::Vector{Float32},
    out_donor_library_hellinger_worst_f4::Vector{Float32},
    pid_v::AbstractVector{UInt32},
    file_v::AbstractVector{UInt32},
    weight_v::AbstractVector{Float32},
    l2ie_v::AbstractVector{Float16},
    irtp_v::AbstractVector{Float32},
    irto_v::AbstractVector{Float32},
    logby_v::Union{Nothing, AbstractVector},
    rt_v::Union{Nothing, AbstractVector{Float32}},
    nscans_v::AbstractVector,
    smoothed_frag_cols,
    fragment_keys::_MBRFragmentAnnotationKeys,
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
    false_donor_cache::Dict{_MBRFalseDonorCacheKey, _MBRFalseDonorTuple},
    counterfactual_eligibility_by_file::Union{Nothing, _CounterfactualEligibilityByFile},
    lod_log2_weight_by_file::Dict{UInt32, Float32},
    lod_log2_weight_global::Float32,
    passing_score_floor::Float32,
)
    n = length(pid_v)
    has_logby = logby_v !== nothing
    has_rt    = rt_v    !== nothing
    @inbounds for i in 1:n
        my_file = file_v[i]
        my_pid  = pid_v[i]
        my_weight = weight_v[i]
        my_l2ie = Float32(l2ie_v[i])
        my_n_scans = Float32(nscans_v[i])
        out_lw_lod[i] = _mbr_log2_weight_lod_ratio(
            my_weight,
            my_file,
            lod_log2_weight_by_file,
            lod_log2_weight_global,
        )
        recipient_sqrt = MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT
        has_recipient_spectrum = false

        donor_t = _donor_for_pid(donor_dict, my_pid, my_file)
        if donor_t !== nothing
            out_best_pair_t[i] = donor_t.trace_prob
            out_lw_t[i] = _mbr_log2_weight_ratio(my_weight, donor_t)
            out_le_t[i] = _mbr_log2_explained_ratio(my_l2ie, donor_t)
            out_nscans_t[i] = _mbr_abs_n_scans_diff(my_n_scans, donor_t)
            out_l2_nscans_t[i] = _mbr_log2_n_scans_ratio(my_n_scans, donor_t)
            out_ir_t[i] = abs((irtp_v[i] - irto_v[i]) - donor_t.irt_residual)
            out_obs_ir_t[i] = abs(irto_v[i] - donor_t.irt_obs)
            if has_rt
                out_rt_t[i] = abs(rt_v[i] - donor_t.rt_obs)
            end
            if has_logby
                out_log_by_t[i] = Float32(logby_v[i]) - donor_t.log_by_ratio
            end
            recipient_sqrt = _mbr_smoothed_spectrum_sqrt_tuple(smoothed_frag_cols, i)
            has_recipient_spectrum = true
            out_smoothed_hellinger_t[i] = _mbr_smoothed_spectrum_hellinger_from_sqrt(
                recipient_sqrt,
                donor_t.smoothed_frag_sqrt,
                fragment_keys,
                my_pid,
                donor_t.precursor_idx,
            )
            out_donor_library_hellinger_t[i] = donor_t.library_hellinger
            donor_worst_t = _min_weight_alternate_donor_for_pid(
                donor_dict,
                my_pid,
                my_file,
                donor_t,
                passing_score_floor,
            )
            if donor_worst_t !== nothing
                out_worst_pair_t[i] = donor_worst_t.trace_prob
                out_lw_worst_t[i] = _mbr_log2_weight_ratio(my_weight, donor_worst_t)
                out_le_worst_t[i] = _mbr_log2_explained_ratio(my_l2ie, donor_worst_t)
                out_nscans_worst_t[i] = _mbr_abs_n_scans_diff(my_n_scans, donor_worst_t)
                out_l2_nscans_worst_t[i] = _mbr_log2_n_scans_ratio(my_n_scans, donor_worst_t)
                out_ir_worst_t[i] = abs((irtp_v[i] - irto_v[i]) - donor_worst_t.irt_residual)
                out_obs_ir_worst_t[i] = abs(irto_v[i] - donor_worst_t.irt_obs)
                out_smoothed_hellinger_worst_t[i] = _mbr_smoothed_spectrum_hellinger_from_donor(
                    recipient_sqrt,
                    donor_worst_t,
                    fragment_keys,
                    my_pid,
                )
                out_donor_library_hellinger_worst_t[i] =
                    _mbr_donor_library_hellinger(donor_worst_t)
            else
                out_single_donor_t[i] = 1.0f0
            end
            out_miss_t[i] = false
        end
        donor_fs = donor_t === nothing ? _empty_false_donors() : _false_donors_for_pid(
            false_donor_cache,
            donor_dict,
            partner_pools,
            my_pid,
            my_file,
            counterfactual_eligibility_by_file,
            donor_t.ms_file_idx,
            Float32(irtp_v[i]),
        )
        donor_f = donor_fs[1]
        if donor_f !== nothing
            out_best_pair_f[i] = donor_f.trace_prob
            out_lw_f[i] = _mbr_log2_weight_ratio(my_weight, donor_f)
            out_le_f[i] = _mbr_log2_explained_ratio(my_l2ie, donor_f)
            out_nscans_f[i] = _mbr_abs_n_scans_diff(my_n_scans, donor_f)
            out_l2_nscans_f[i] = _mbr_log2_n_scans_ratio(my_n_scans, donor_f)
            out_ir_f[i] = abs((irtp_v[i] - irto_v[i]) - donor_f.irt_residual)
            out_obs_ir_f[i] = abs(irto_v[i] - donor_f.irt_obs)
            if has_rt
                out_rt_f[i] = abs(rt_v[i] - donor_f.rt_obs)
            end
            if has_logby
                out_log_by_f[i] = Float32(logby_v[i]) - donor_f.log_by_ratio
            end
            if !has_recipient_spectrum
                recipient_sqrt = _mbr_smoothed_spectrum_sqrt_tuple(smoothed_frag_cols, i)
            end
            out_smoothed_hellinger_f[i] = _mbr_smoothed_spectrum_hellinger_from_sqrt(
                recipient_sqrt,
                donor_f.smoothed_frag_sqrt,
                fragment_keys,
                my_pid,
                donor_f.precursor_idx,
            )
            out_donor_library_hellinger_f[i] = donor_f.library_hellinger
            donor_worst_f = _min_weight_alternate_donor_for_pid(
                donor_dict,
                donor_f.precursor_idx,
                my_file,
                donor_f,
                passing_score_floor,
            )
            if donor_worst_f !== nothing
                out_worst_pair_f[i] = donor_worst_f.trace_prob
                out_lw_worst_f[i] = _mbr_log2_weight_ratio(my_weight, donor_worst_f)
                out_le_worst_f[i] = _mbr_log2_explained_ratio(my_l2ie, donor_worst_f)
                out_nscans_worst_f[i] = _mbr_abs_n_scans_diff(my_n_scans, donor_worst_f)
                out_l2_nscans_worst_f[i] = _mbr_log2_n_scans_ratio(my_n_scans, donor_worst_f)
                out_ir_worst_f[i] = abs((irtp_v[i] - irto_v[i]) - donor_worst_f.irt_residual)
                out_obs_ir_worst_f[i] = abs(irto_v[i] - donor_worst_f.irt_obs)
                out_smoothed_hellinger_worst_f[i] = _mbr_smoothed_spectrum_hellinger_from_donor(
                    recipient_sqrt,
                    donor_worst_f,
                    fragment_keys,
                    my_pid,
                )
                out_donor_library_hellinger_worst_f[i] =
                    _mbr_donor_library_hellinger(donor_worst_f)
            else
                out_single_donor_f[i] = 1.0f0
            end
            out_miss_f[i] = false
        end
        donor_f2 = donor_fs[2]
        if donor_f2 !== nothing
            out_best_pair_f2[i] = donor_f2.trace_prob
            out_lw_f2[i] = _mbr_log2_weight_ratio(my_weight, donor_f2)
            out_le_f2[i] = _mbr_log2_explained_ratio(my_l2ie, donor_f2)
            out_nscans_f2[i] = _mbr_abs_n_scans_diff(my_n_scans, donor_f2)
            out_l2_nscans_f2[i] = _mbr_log2_n_scans_ratio(my_n_scans, donor_f2)
            out_ir_f2[i] = abs((irtp_v[i] - irto_v[i]) - donor_f2.irt_residual)
            out_obs_ir_f2[i] = abs(irto_v[i] - donor_f2.irt_obs)
            if has_rt
                out_rt_f2[i] = abs(rt_v[i] - donor_f2.rt_obs)
            end
            if has_logby
                out_log_by_f2[i] = Float32(logby_v[i]) - donor_f2.log_by_ratio
            end
            if !has_recipient_spectrum
                recipient_sqrt = _mbr_smoothed_spectrum_sqrt_tuple(smoothed_frag_cols, i)
                has_recipient_spectrum = true
            end
            out_smoothed_hellinger_f2[i] = _mbr_smoothed_spectrum_hellinger_from_sqrt(
                recipient_sqrt,
                donor_f2.smoothed_frag_sqrt,
                fragment_keys,
                my_pid,
                donor_f2.precursor_idx,
            )
            out_donor_library_hellinger_f2[i] = donor_f2.library_hellinger
            donor_worst_f2 = _min_weight_alternate_donor_for_pid(
                donor_dict,
                donor_f2.precursor_idx,
                my_file,
                donor_f2,
                passing_score_floor,
            )
            if donor_worst_f2 !== nothing
                out_worst_pair_f2[i] = donor_worst_f2.trace_prob
                out_lw_worst_f2[i] = _mbr_log2_weight_ratio(my_weight, donor_worst_f2)
                out_le_worst_f2[i] = _mbr_log2_explained_ratio(my_l2ie, donor_worst_f2)
                out_nscans_worst_f2[i] = _mbr_abs_n_scans_diff(my_n_scans, donor_worst_f2)
                out_l2_nscans_worst_f2[i] = _mbr_log2_n_scans_ratio(my_n_scans, donor_worst_f2)
                out_ir_worst_f2[i] = abs((irtp_v[i] - irto_v[i]) - donor_worst_f2.irt_residual)
                out_obs_ir_worst_f2[i] = abs(irto_v[i] - donor_worst_f2.irt_obs)
                out_smoothed_hellinger_worst_f2[i] = _mbr_smoothed_spectrum_hellinger_from_donor(
                    recipient_sqrt,
                    donor_worst_f2,
                    fragment_keys,
                    my_pid,
                )
                out_donor_library_hellinger_worst_f2[i] =
                    _mbr_donor_library_hellinger(donor_worst_f2)
            else
                out_single_donor_f2[i] = 1.0f0
            end
            out_miss_f2[i] = false
        end
        donor_f3 = donor_fs[3]
        if donor_f3 !== nothing
            out_best_pair_f3[i] = donor_f3.trace_prob
            out_lw_f3[i] = _mbr_log2_weight_ratio(my_weight, donor_f3)
            out_le_f3[i] = _mbr_log2_explained_ratio(my_l2ie, donor_f3)
            out_nscans_f3[i] = _mbr_abs_n_scans_diff(my_n_scans, donor_f3)
            out_l2_nscans_f3[i] = _mbr_log2_n_scans_ratio(my_n_scans, donor_f3)
            out_ir_f3[i] = abs((irtp_v[i] - irto_v[i]) - donor_f3.irt_residual)
            out_obs_ir_f3[i] = abs(irto_v[i] - donor_f3.irt_obs)
            if has_rt
                out_rt_f3[i] = abs(rt_v[i] - donor_f3.rt_obs)
            end
            if has_logby
                out_log_by_f3[i] = Float32(logby_v[i]) - donor_f3.log_by_ratio
            end
            if !has_recipient_spectrum
                recipient_sqrt = _mbr_smoothed_spectrum_sqrt_tuple(smoothed_frag_cols, i)
                has_recipient_spectrum = true
            end
            out_smoothed_hellinger_f3[i] = _mbr_smoothed_spectrum_hellinger_from_sqrt(
                recipient_sqrt,
                donor_f3.smoothed_frag_sqrt,
                fragment_keys,
                my_pid,
                donor_f3.precursor_idx,
            )
            out_donor_library_hellinger_f3[i] = donor_f3.library_hellinger
            donor_worst_f3 = _min_weight_alternate_donor_for_pid(
                donor_dict,
                donor_f3.precursor_idx,
                my_file,
                donor_f3,
                passing_score_floor,
            )
            if donor_worst_f3 !== nothing
                out_worst_pair_f3[i] = donor_worst_f3.trace_prob
                out_lw_worst_f3[i] = _mbr_log2_weight_ratio(my_weight, donor_worst_f3)
                out_le_worst_f3[i] = _mbr_log2_explained_ratio(my_l2ie, donor_worst_f3)
                out_nscans_worst_f3[i] = _mbr_abs_n_scans_diff(my_n_scans, donor_worst_f3)
                out_l2_nscans_worst_f3[i] = _mbr_log2_n_scans_ratio(my_n_scans, donor_worst_f3)
                out_ir_worst_f3[i] = abs((irtp_v[i] - irto_v[i]) - donor_worst_f3.irt_residual)
                out_obs_ir_worst_f3[i] = abs(irto_v[i] - donor_worst_f3.irt_obs)
                out_smoothed_hellinger_worst_f3[i] = _mbr_smoothed_spectrum_hellinger_from_donor(
                    recipient_sqrt,
                    donor_worst_f3,
                    fragment_keys,
                    my_pid,
                )
                out_donor_library_hellinger_worst_f3[i] =
                    _mbr_donor_library_hellinger(donor_worst_f3)
            else
                out_single_donor_f3[i] = 1.0f0
            end
            out_miss_f3[i] = false
        end
        donor_f4 = donor_fs[4]
        if donor_f4 !== nothing
            out_best_pair_f4[i] = donor_f4.trace_prob
            out_lw_f4[i] = _mbr_log2_weight_ratio(my_weight, donor_f4)
            out_le_f4[i] = _mbr_log2_explained_ratio(my_l2ie, donor_f4)
            out_nscans_f4[i] = _mbr_abs_n_scans_diff(my_n_scans, donor_f4)
            out_l2_nscans_f4[i] = _mbr_log2_n_scans_ratio(my_n_scans, donor_f4)
            out_ir_f4[i] = abs((irtp_v[i] - irto_v[i]) - donor_f4.irt_residual)
            out_obs_ir_f4[i] = abs(irto_v[i] - donor_f4.irt_obs)
            if has_rt
                out_rt_f4[i] = abs(rt_v[i] - donor_f4.rt_obs)
            end
            if has_logby
                out_log_by_f4[i] = Float32(logby_v[i]) - donor_f4.log_by_ratio
            end
            if !has_recipient_spectrum
                recipient_sqrt = _mbr_smoothed_spectrum_sqrt_tuple(smoothed_frag_cols, i)
                has_recipient_spectrum = true
            end
            out_smoothed_hellinger_f4[i] = _mbr_smoothed_spectrum_hellinger_from_sqrt(
                recipient_sqrt,
                donor_f4.smoothed_frag_sqrt,
                fragment_keys,
                my_pid,
                donor_f4.precursor_idx,
            )
            out_donor_library_hellinger_f4[i] = donor_f4.library_hellinger
            donor_worst_f4 = _min_weight_alternate_donor_for_pid(
                donor_dict,
                donor_f4.precursor_idx,
                my_file,
                donor_f4,
                passing_score_floor,
            )
            if donor_worst_f4 !== nothing
                out_worst_pair_f4[i] = donor_worst_f4.trace_prob
                out_lw_worst_f4[i] = _mbr_log2_weight_ratio(my_weight, donor_worst_f4)
                out_le_worst_f4[i] = _mbr_log2_explained_ratio(my_l2ie, donor_worst_f4)
                out_nscans_worst_f4[i] = _mbr_abs_n_scans_diff(my_n_scans, donor_worst_f4)
                out_l2_nscans_worst_f4[i] = _mbr_log2_n_scans_ratio(my_n_scans, donor_worst_f4)
                out_ir_worst_f4[i] = abs((irtp_v[i] - irto_v[i]) - donor_worst_f4.irt_residual)
                out_obs_ir_worst_f4[i] = abs(irto_v[i] - donor_worst_f4.irt_obs)
                out_smoothed_hellinger_worst_f4[i] = _mbr_smoothed_spectrum_hellinger_from_donor(
                    recipient_sqrt,
                    donor_worst_f4,
                    fragment_keys,
                    my_pid,
                )
                out_donor_library_hellinger_worst_f4[i] =
                    _mbr_donor_library_hellinger(donor_worst_f4)
            else
                out_single_donor_f4[i] = 1.0f0
            end
            out_miss_f4[i] = false
        end
    end
    return nothing
end

@inline function _compute_mbr_inner_v2!(
    out_best_pair_t::Vector{Float32}, out_best_pair_f::Vector{Vector{Float32}},
    out_worst_pair_t::Vector{Float32}, out_worst_pair_f::Vector{Vector{Float32}},
    out_lw_lod::Vector{Float32},
    out_lw_t::Vector{Float32}, out_lw_f::Vector{Vector{Float32}},
    out_lw_worst_t::Vector{Float32}, out_lw_worst_f::Vector{Vector{Float32}},
    out_le_t::Vector{Float32}, out_le_f::Vector{Vector{Float32}},
    out_le_worst_t::Vector{Float32}, out_le_worst_f::Vector{Vector{Float32}},
    out_nscans_t::Vector{Float32}, out_nscans_f::Vector{Vector{Float32}},
    out_nscans_worst_t::Vector{Float32}, out_nscans_worst_f::Vector{Vector{Float32}},
    out_l2_nscans_t::Vector{Float32}, out_l2_nscans_f::Vector{Vector{Float32}},
    out_l2_nscans_worst_t::Vector{Float32}, out_l2_nscans_worst_f::Vector{Vector{Float32}},
    out_ir_t::Vector{Float32}, out_ir_f::Vector{Vector{Float32}},
    out_ir_worst_t::Vector{Float32}, out_ir_worst_f::Vector{Vector{Float32}},
    out_obs_ir_t::Vector{Float32}, out_obs_ir_f::Vector{Vector{Float32}},
    out_obs_ir_worst_t::Vector{Float32}, out_obs_ir_worst_f::Vector{Vector{Float32}},
    out_single_donor_t::Vector{Float32}, out_single_donor_f::Vector{Vector{Float32}},
    out_miss_t::BitVector, out_miss_f::Vector{BitVector},
    out_rt_t::Vector{Float32}, out_rt_f::Vector{Vector{Float32}},
    out_log_by_t::Vector{Float32}, out_log_by_f::Vector{Vector{Float32}},
    out_smoothed_hellinger_t::Vector{Float32},
    out_smoothed_hellinger_f::Vector{Vector{Float32}},
    out_smoothed_hellinger_worst_t::Vector{Float32},
    out_smoothed_hellinger_worst_f::Vector{Vector{Float32}},
    out_donor_library_hellinger_t::Vector{Float32},
    out_donor_library_hellinger_f::Vector{Vector{Float32}},
    out_donor_library_hellinger_worst_t::Vector{Float32},
    out_donor_library_hellinger_worst_f::Vector{Vector{Float32}},
    pid_v::AbstractVector{UInt32},
    file_v::AbstractVector{UInt32},
    weight_v::AbstractVector{Float32},
    l2ie_v::AbstractVector{Float16},
    irtp_v::AbstractVector{Float32},
    irto_v::AbstractVector{Float32},
    logby_v::Union{Nothing, AbstractVector},
    rt_v::Union{Nothing, AbstractVector{Float32}},
    nscans_v::AbstractVector,
    smoothed_frag_cols,
    fragment_keys::_MBRFragmentAnnotationKeys,
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_pools::_CounterfactualPartnerPools,
    false_donor_cache::Dict{_MBRFalseDonorCacheKey, _MBRFalseDonorTuple},
    counterfactual_eligibility_by_file::Union{Nothing, _CounterfactualEligibilityByFile},
    lod_log2_weight_by_file::Dict{UInt32, Float32},
    lod_log2_weight_global::Float32,
    passing_score_floor::Float32,
)
    n = length(pid_v)
    has_logby = logby_v !== nothing
    has_rt = rt_v !== nothing
    @inbounds for i in 1:n
        my_file = file_v[i]
        my_pid = pid_v[i]
        my_weight = weight_v[i]
        my_l2ie = Float32(l2ie_v[i])
        my_n_scans = Float32(nscans_v[i])
        out_lw_lod[i] = _mbr_log2_weight_lod_ratio(
            my_weight,
            my_file,
            lod_log2_weight_by_file,
            lod_log2_weight_global,
        )
        recipient_sqrt = MBR_SMOOTHED_SPECTRUM_EMPTY_SQRT
        has_recipient_spectrum = false

        donor_t = _donor_for_pid(donor_dict, my_pid, my_file)
        if donor_t !== nothing
            out_best_pair_t[i] = donor_t.trace_prob
            out_lw_t[i] = _mbr_log2_weight_ratio(my_weight, donor_t)
            out_le_t[i] = _mbr_log2_explained_ratio(my_l2ie, donor_t)
            out_nscans_t[i] = _mbr_abs_n_scans_diff(my_n_scans, donor_t)
            out_l2_nscans_t[i] = _mbr_log2_n_scans_ratio(my_n_scans, donor_t)
            out_ir_t[i] = abs((irtp_v[i] - irto_v[i]) - donor_t.irt_residual)
            out_obs_ir_t[i] = abs(irto_v[i] - donor_t.irt_obs)
            if has_rt
                out_rt_t[i] = abs(rt_v[i] - donor_t.rt_obs)
            end
            if has_logby
                out_log_by_t[i] = Float32(logby_v[i]) - donor_t.log_by_ratio
            end
            recipient_sqrt = _mbr_smoothed_spectrum_sqrt_tuple(smoothed_frag_cols, i)
            has_recipient_spectrum = true
            out_smoothed_hellinger_t[i] = _mbr_smoothed_spectrum_hellinger_from_sqrt(
                recipient_sqrt,
                donor_t.smoothed_frag_sqrt,
                fragment_keys,
                my_pid,
                donor_t.precursor_idx,
            )
            out_donor_library_hellinger_t[i] = donor_t.library_hellinger
            donor_worst_t = _min_weight_alternate_donor_for_pid(
                donor_dict,
                my_pid,
                my_file,
                donor_t,
                passing_score_floor,
            )
            if donor_worst_t !== nothing
                out_worst_pair_t[i] = donor_worst_t.trace_prob
                out_lw_worst_t[i] = _mbr_log2_weight_ratio(my_weight, donor_worst_t)
                out_le_worst_t[i] = _mbr_log2_explained_ratio(my_l2ie, donor_worst_t)
                out_nscans_worst_t[i] = _mbr_abs_n_scans_diff(my_n_scans, donor_worst_t)
                out_l2_nscans_worst_t[i] = _mbr_log2_n_scans_ratio(my_n_scans, donor_worst_t)
                out_ir_worst_t[i] = abs((irtp_v[i] - irto_v[i]) - donor_worst_t.irt_residual)
                out_obs_ir_worst_t[i] = abs(irto_v[i] - donor_worst_t.irt_obs)
                out_smoothed_hellinger_worst_t[i] = _mbr_smoothed_spectrum_hellinger_from_donor(
                    recipient_sqrt,
                    donor_worst_t,
                    fragment_keys,
                    my_pid,
                )
                out_donor_library_hellinger_worst_t[i] =
                    _mbr_donor_library_hellinger(donor_worst_t)
            else
                out_single_donor_t[i] = 1.0f0
            end
            out_miss_t[i] = false
        end

        donor_fs = donor_t === nothing ? _empty_false_donors() : _false_donors_for_pid(
            false_donor_cache,
            donor_dict,
            partner_pools,
            my_pid,
            my_file,
            counterfactual_eligibility_by_file,
            donor_t.ms_file_idx,
            Float32(irtp_v[i]),
        )
        for counterfactual_idx in 1:MBR_MAX_COUNTERFACTUALS
            donor_f = donor_fs[counterfactual_idx]
            donor_f === nothing && continue

            out_best_pair_f[counterfactual_idx][i] = donor_f.trace_prob
            out_lw_f[counterfactual_idx][i] = _mbr_log2_weight_ratio(my_weight, donor_f)
            out_le_f[counterfactual_idx][i] = _mbr_log2_explained_ratio(my_l2ie, donor_f)
            out_nscans_f[counterfactual_idx][i] = _mbr_abs_n_scans_diff(my_n_scans, donor_f)
            out_l2_nscans_f[counterfactual_idx][i] = _mbr_log2_n_scans_ratio(my_n_scans, donor_f)
            out_ir_f[counterfactual_idx][i] = abs((irtp_v[i] - irto_v[i]) - donor_f.irt_residual)
            out_obs_ir_f[counterfactual_idx][i] = abs(irto_v[i] - donor_f.irt_obs)
            if has_rt
                out_rt_f[counterfactual_idx][i] = abs(rt_v[i] - donor_f.rt_obs)
            end
            if has_logby
                out_log_by_f[counterfactual_idx][i] = Float32(logby_v[i]) - donor_f.log_by_ratio
            end
            if !has_recipient_spectrum
                recipient_sqrt = _mbr_smoothed_spectrum_sqrt_tuple(smoothed_frag_cols, i)
                has_recipient_spectrum = true
            end
            out_smoothed_hellinger_f[counterfactual_idx][i] = _mbr_smoothed_spectrum_hellinger_from_sqrt(
                recipient_sqrt,
                donor_f.smoothed_frag_sqrt,
                fragment_keys,
                my_pid,
                donor_f.precursor_idx,
            )
            out_donor_library_hellinger_f[counterfactual_idx][i] = donor_f.library_hellinger
            donor_worst_f = _min_weight_alternate_donor_for_pid(
                donor_dict,
                donor_f.precursor_idx,
                my_file,
                donor_f,
                passing_score_floor,
            )
            if donor_worst_f !== nothing
                out_worst_pair_f[counterfactual_idx][i] = donor_worst_f.trace_prob
                out_lw_worst_f[counterfactual_idx][i] = _mbr_log2_weight_ratio(my_weight, donor_worst_f)
                out_le_worst_f[counterfactual_idx][i] = _mbr_log2_explained_ratio(my_l2ie, donor_worst_f)
                out_nscans_worst_f[counterfactual_idx][i] = _mbr_abs_n_scans_diff(my_n_scans, donor_worst_f)
                out_l2_nscans_worst_f[counterfactual_idx][i] = _mbr_log2_n_scans_ratio(my_n_scans, donor_worst_f)
                out_ir_worst_f[counterfactual_idx][i] = abs((irtp_v[i] - irto_v[i]) - donor_worst_f.irt_residual)
                out_obs_ir_worst_f[counterfactual_idx][i] = abs(irto_v[i] - donor_worst_f.irt_obs)
                out_smoothed_hellinger_worst_f[counterfactual_idx][i] = _mbr_smoothed_spectrum_hellinger_from_donor(
                    recipient_sqrt,
                    donor_worst_f,
                    fragment_keys,
                    my_pid,
                )
                out_donor_library_hellinger_worst_f[counterfactual_idx][i] =
                    _mbr_donor_library_hellinger(donor_worst_f)
            else
                out_single_donor_f[counterfactual_idx][i] = 1.0f0
            end
            out_miss_f[counterfactual_idx][i] = false
        end
    end
    return nothing
end

@inline _mbr_float_counterfactual_vectors(n::Int; sentinel::Float32 = -1.0f0) =
    [fill(sentinel, n) for _ in 1:MBR_MAX_COUNTERFACTUALS]

@inline _mbr_bit_counterfactual_vectors(n::Int) =
    [trues(n) for _ in 1:MBR_MAX_COUNTERFACTUALS]

function _mbr_add_counterfactual_columns!(
    df::DataFrame,
    base::AbstractString,
    values::Union{Vector{Vector{Float32}}, Vector{BitVector}},
)
    @inbounds for counterfactual_idx in 1:MBR_MAX_COUNTERFACTUALS
        df[!, _mbr_counterfactual_column(base, counterfactual_idx)] = values[counterfactual_idx]
    end
    return nothing
end

# Per-file MBR feature compute + sidecar write, reading trace_prob_prepass
# from the Pass-1 sidecar. Mirrors compute_mbr_features_per_file_to_sidecar!
# but for the streaming-with-Pass-1 flow.
function compute_mbr_features_per_file_to_sidecar_with_pass1!(
        main_path::AbstractString,
        donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
        partner_pools::_CounterfactualPartnerPools,
        fragment_keys::_MBRFragmentAnnotationKeys;
        counterfactual_eligibility_by_file::Union{Nothing, _CounterfactualEligibilityByFile} = nothing,
        passing_score_floor::Float32 = Float32(Inf),
        lod_log2_weight_by_file::Dict{UInt32, Float32} = Dict{UInt32, Float32}(),
        lod_log2_weight_global::Float32 = NaN32)
    pass1_path = main_path * PASS1_SIDECAR_SUFFIX
    isfile(pass1_path) || error("Missing Pass-1 sidecar at $pass1_path")
    main = Arrow.Table(main_path)
    pass1 = Arrow.Table(pass1_path)
    n = length(main.precursor_idx)
    n == length(pass1.precursor_idx) ||
        error("Pass-1 sidecar row count mismatch at $pass1_path")

    pid_v   = main.precursor_idx
    scan_v  = main.scan_idx
    file_v  = main.ms_file_idx
    weight_v= main.weight
    l2ie_v  = main.log2_intensity_explained
    irtp_v  = main.irt_pred
    irto_v  = main.irt_obs
    logby_v = hasproperty(main, :log_by_ratio_m0) ? main.log_by_ratio_m0 : nothing
    rt_v    = hasproperty(main, :rt) ? main.rt : nothing
    nscans_v = main.n_scans
    smoothed_frag_cols = ntuple(rank -> getproperty(main, SMOOTHED_FRAGMENT_INTENSITY_COLUMNS[rank]), 8)

    @inbounds for i in 1:n
        (pid_v[i] == pass1.precursor_idx[i] && scan_v[i] == pass1.scan_idx[i]) ||
            error("Pass-1 sidecar misaligned at row $i of $pass1_path")
    end

    out_best_pair_t = fill(-1f0, n); out_best_pair_f = _mbr_float_counterfactual_vectors(n)
    out_worst_pair_t = fill(-1f0, n); out_worst_pair_f = _mbr_float_counterfactual_vectors(n)
    out_lw_lod = fill(-1f0, n)
    out_lw_t  = fill(-1f0, n); out_lw_f  = _mbr_float_counterfactual_vectors(n)
    out_lw_worst_t = fill(-1f0, n); out_lw_worst_f = _mbr_float_counterfactual_vectors(n)
    out_le_t  = fill(-1f0, n); out_le_f  = _mbr_float_counterfactual_vectors(n)
    out_le_worst_t = fill(-1f0, n); out_le_worst_f = _mbr_float_counterfactual_vectors(n)
    out_nscans_t = fill(-1f0, n); out_nscans_f = _mbr_float_counterfactual_vectors(n)
    out_nscans_worst_t = fill(-1f0, n); out_nscans_worst_f = _mbr_float_counterfactual_vectors(n)
    out_l2_nscans_t = fill(-1f0, n); out_l2_nscans_f = _mbr_float_counterfactual_vectors(n)
    out_l2_nscans_worst_t = fill(-1f0, n); out_l2_nscans_worst_f = _mbr_float_counterfactual_vectors(n)
    out_ir_t  = fill(-1f0, n); out_ir_f  = _mbr_float_counterfactual_vectors(n)
    out_ir_worst_t = fill(-1f0, n); out_ir_worst_f = _mbr_float_counterfactual_vectors(n)
    out_obs_ir_t = fill(-1f0, n); out_obs_ir_f = _mbr_float_counterfactual_vectors(n)
    out_obs_ir_worst_t = fill(-1f0, n); out_obs_ir_worst_f = _mbr_float_counterfactual_vectors(n)
    out_single_donor_t = fill(0f0, n); out_single_donor_f = _mbr_float_counterfactual_vectors(n, sentinel = 0.0f0)
    out_miss_t = trues(n); out_miss_f = _mbr_bit_counterfactual_vectors(n)
    out_rt_t = fill(-1f0, n); out_rt_f = _mbr_float_counterfactual_vectors(n)
    out_log_by_t = fill(-1f0, n); out_log_by_f = _mbr_float_counterfactual_vectors(n)
    out_smoothed_h_t = fill(-1f0, n); out_smoothed_h_f = _mbr_float_counterfactual_vectors(n)
    out_smoothed_h_worst_t = fill(1f0, n); out_smoothed_h_worst_f = _mbr_float_counterfactual_vectors(n, sentinel = 1.0f0)
    out_donor_library_h_t = fill(1f0, n); out_donor_library_h_f = _mbr_float_counterfactual_vectors(n, sentinel = 1.0f0)
    out_donor_library_h_worst_t = fill(1f0, n)
    out_donor_library_h_worst_f = _mbr_float_counterfactual_vectors(n, sentinel = 1.0f0)
    false_donor_cache = Dict{_MBRFalseDonorCacheKey, _MBRFalseDonorTuple}()

    _compute_mbr_inner_v2!(
        out_best_pair_t, out_best_pair_f,
        out_worst_pair_t, out_worst_pair_f,
        out_lw_lod,
        out_lw_t, out_lw_f,
        out_lw_worst_t, out_lw_worst_f,
        out_le_t, out_le_f,
        out_le_worst_t, out_le_worst_f,
        out_nscans_t, out_nscans_f,
        out_nscans_worst_t, out_nscans_worst_f,
        out_l2_nscans_t, out_l2_nscans_f,
        out_l2_nscans_worst_t, out_l2_nscans_worst_f,
        out_ir_t, out_ir_f,
        out_ir_worst_t, out_ir_worst_f,
        out_obs_ir_t, out_obs_ir_f,
        out_obs_ir_worst_t, out_obs_ir_worst_f,
        out_single_donor_t, out_single_donor_f,
        out_miss_t, out_miss_f,
        out_rt_t, out_rt_f,
        out_log_by_t, out_log_by_f,
        out_smoothed_h_t, out_smoothed_h_f,
        out_smoothed_h_worst_t, out_smoothed_h_worst_f,
        out_donor_library_h_t, out_donor_library_h_f,
        out_donor_library_h_worst_t, out_donor_library_h_worst_f,
        pid_v, file_v, weight_v, l2ie_v, irtp_v, irto_v,
        logby_v,
        rt_v, nscans_v, smoothed_frag_cols, fragment_keys,
        donor_dict, partner_pools, false_donor_cache,
        counterfactual_eligibility_by_file,
        lod_log2_weight_by_file, lod_log2_weight_global, passing_score_floor,
    )

    side_path = main_path * MBR_SIDECAR_SUFFIX
    side_df = DataFrame(
        precursor_idx = collect(UInt32.(pid_v)),
        scan_idx = collect(UInt32.(scan_v)),
    )
    side_df[!, :MBR_best_pair_prob_true] = out_best_pair_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_best_pair_prob", out_best_pair_f)
    side_df[!, :MBR_worst_pair_prob_true] = out_worst_pair_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_worst_pair_prob", out_worst_pair_f)
    side_df[!, :MBR_log2_weight_lod_ratio] = out_lw_lod
    side_df[!, :MBR_best_log2_weight_ratio_true] = out_lw_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_best_log2_weight_ratio", out_lw_f)
    side_df[!, :MBR_worst_log2_weight_ratio_true] = out_lw_worst_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_worst_log2_weight_ratio", out_lw_worst_f)
    side_df[!, :MBR_best_log2_explained_ratio_true] = out_le_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_best_log2_explained_ratio", out_le_f)
    side_df[!, :MBR_worst_log2_explained_ratio_true] = out_le_worst_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_worst_log2_explained_ratio", out_le_worst_f)
    side_df[!, :MBR_best_abs_n_scans_diff_true] = out_nscans_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_best_abs_n_scans_diff", out_nscans_f)
    side_df[!, :MBR_worst_abs_n_scans_diff_true] = out_nscans_worst_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_worst_abs_n_scans_diff", out_nscans_worst_f)
    side_df[!, :MBR_best_log2_n_scans_ratio_true] = out_l2_nscans_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_best_log2_n_scans_ratio", out_l2_nscans_f)
    side_df[!, :MBR_worst_log2_n_scans_ratio_true] = out_l2_nscans_worst_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_worst_log2_n_scans_ratio", out_l2_nscans_worst_f)
    side_df[!, :MBR_best_irt_diff_true] = out_ir_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_best_irt_diff", out_ir_f)
    side_df[!, :MBR_worst_irt_diff_true] = out_ir_worst_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_worst_irt_diff", out_ir_worst_f)
    side_df[!, :MBR_best_observed_irt_diff_true] = out_obs_ir_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_best_observed_irt_diff", out_obs_ir_f)
    side_df[!, :MBR_worst_observed_irt_diff_true] = out_obs_ir_worst_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_worst_observed_irt_diff", out_obs_ir_worst_f)
    side_df[!, :MBR_single_donor_true] = out_single_donor_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_single_donor", out_single_donor_f)
    side_df[!, :MBR_best_is_missing_true] = out_miss_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_best_is_missing", out_miss_f)
    side_df[!, :MBR_best_rt_diff_true] = out_rt_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_best_rt_diff", out_rt_f)
    side_df[!, :MBR_best_log_by_diff_true] = out_log_by_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_best_log_by_diff", out_log_by_f)
    side_df[!, :MBR_best_smoothed_frag_hellinger_true] = out_smoothed_h_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_best_smoothed_frag_hellinger", out_smoothed_h_f)
    side_df[!, :MBR_worst_smoothed_frag_hellinger_true] = out_smoothed_h_worst_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_worst_smoothed_frag_hellinger", out_smoothed_h_worst_f)
    side_df[!, :MBR_best_donor_library_hellinger_true] = out_donor_library_h_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_best_donor_library_hellinger", out_donor_library_h_f)
    side_df[!, :MBR_worst_donor_library_hellinger_true] = out_donor_library_h_worst_t
    _mbr_add_counterfactual_columns!(side_df, "MBR_worst_donor_library_hellinger", out_donor_library_h_worst_f)
    writeArrow(side_path, side_df)
    return side_path
end

# Slim DataFrame loader for the FTR controller: pulls only the columns
# apply_mbr_filter_paired! needs, from main + Pass-1 sidecar + MBR sidecar,
# in main-file row order across all files. Substantially smaller than
# best_psms (≈20 cols vs ≈80).
function load_ftr_slim_dataframe(file_paths::Vector{String})
    parts = DataFrame[]
    for path in file_paths
        pass1_path = path * PASS1_SIDECAR_SUFFIX
        mbr_path   = path * MBR_SIDECAR_SUFFIX
        isfile(pass1_path) || error("Missing Pass-1 sidecar at $pass1_path")
        isfile(mbr_path)   || error("Missing MBR sidecar at $mbr_path")
        main  = Arrow.Table(path)
        pass1 = Arrow.Table(pass1_path)
        mbr   = Arrow.Table(mbr_path)
        n = length(main.precursor_idx)
        n == length(pass1.precursor_idx) || error("Pass-1 sidecar size mismatch at $pass1_path")
        n == length(mbr.precursor_idx)   || error("MBR sidecar size mismatch at $mbr_path")
        @inbounds for i in 1:n
            (main.precursor_idx[i] == pass1.precursor_idx[i] == mbr.precursor_idx[i] &&
             main.scan_idx[i]      == pass1.scan_idx[i]      == mbr.scan_idx[i]) ||
                error("Sidecar misaligned at row $i for $path")
        end
        d = DataFrame(
            precursor_idx       = collect(UInt32.(main.precursor_idx)),
            scan_idx            = collect(UInt32.(main.scan_idx)),
            ms_file_idx         = collect(UInt32.(main.ms_file_idx)),
            cv_fold             = collect(UInt8.(main.cv_fold)),
            target              = collect(Bool.(main.target)),
            decoy               = .!collect(Bool.(main.target)),
            trace_prob_prepass  = collect(Float32.(pass1.trace_prob_prepass)),
            trace_prob_infold   = collect(Float32.(pass1.trace_prob_infold)),
        )
        if hasproperty(main, :qval)
            d[!, :qval] = collect(Float32.(main.qval))
        end
        if hasproperty(main, :global_qval)
            d[!, :global_qval] = collect(Float32.(main.global_qval))
        end
        for c in MBR_CROSS_RUN_FTR_FEATURES
            hasproperty(main, c) || continue
            d[!, c] = collect(Tables.getcolumn(main, c))
        end
        # Pull all MBR_* columns from the MBR sidecar.
        for c in _MBR_SIDECAR_OUT_COLS
            c === :precursor_idx && continue
            c === :scan_idx      && continue
            d[!, c] = collect(Tables.getcolumn(mbr, c))
        end
        push!(parts, d)
    end
    return vcat(parts...)
end

# After apply_mbr_filter_paired! adds the four recovery columns to the slim
# DataFrame, distribute them back to per-file recovery sidecars in the SAME
# row order as the main files.
function write_recovery_sidecars(slim_df::DataFrame, file_paths::Vector{String})
    # Build (ms_file_idx, cv_fold) → main path map.
    is_fold_split = any(p -> occursin("_fold", p), file_paths)
    key_to_path = if is_fold_split
        d = Dict{Tuple{UInt32, UInt8}, String}()
        for fpath in file_paths
            fold_match = match(r"_fold(\d)\.arrow$", fpath)
            fold_match === nothing && continue
            fold_num = parse(UInt8, fold_match.captures[1])
            orig = Arrow.Table(fpath)
            length(orig.ms_file_idx) == 0 && continue
            d[(UInt32(first(orig.ms_file_idx)), fold_num)] = fpath
        end
        d
    else
        d = Dict{UInt32, String}()
        for fpath in file_paths
            orig = Arrow.Table(fpath)
            length(orig.ms_file_idx) == 0 && continue
            d[UInt32(first(orig.ms_file_idx))] = fpath
        end
        d
    end

    n_written = 0
    if is_fold_split && hasproperty(slim_df, :cv_fold)
        for (key, g) in pairs(groupby(slim_df, [:ms_file_idx, :cv_fold]))
            lookup_key = (UInt32(key[:ms_file_idx]), UInt8(key[:cv_fold]))
            haskey(key_to_path, lookup_key) || continue
            main_path = key_to_path[lookup_key]
            side_path = main_path * RECOVERY_SIDECAR_SUFFIX
            d = DataFrame(
                precursor_idx          = collect(UInt32.(g[!, :precursor_idx])),
                scan_idx               = collect(UInt32.(g[!, :scan_idx])),
                mbr_recovered          = collect(Bool.(g[!, :mbr_recovered])),
                MBR_transfer_candidate = collect(Bool.(g[!, :MBR_transfer_candidate])),
                mbr_target_decoy_prob  = collect(Float32.(g[!, :mbr_target_decoy_prob])),
                ftr_qval_true          = collect(Float32.(g[!, :ftr_qval_true])),
                ftr_pep_true           = collect(Float32.(g[!, :ftr_pep_true])),
            )
            writeArrow(side_path, d)
            n_written += 1
        end
    else
        for (key, g) in pairs(groupby(slim_df, :ms_file_idx))
            lookup_key = UInt32(key[:ms_file_idx])
            haskey(key_to_path, lookup_key) || continue
            main_path = key_to_path[lookup_key]
            side_path = main_path * RECOVERY_SIDECAR_SUFFIX
            d = DataFrame(
                precursor_idx          = collect(UInt32.(g[!, :precursor_idx])),
                scan_idx               = collect(UInt32.(g[!, :scan_idx])),
                mbr_recovered          = collect(Bool.(g[!, :mbr_recovered])),
                MBR_transfer_candidate = collect(Bool.(g[!, :MBR_transfer_candidate])),
                mbr_target_decoy_prob  = collect(Float32.(g[!, :mbr_target_decoy_prob])),
                ftr_qval_true          = collect(Float32.(g[!, :ftr_qval_true])),
                ftr_pep_true           = collect(Float32.(g[!, :ftr_pep_true])),
            )
            writeArrow(side_path, d)
            n_written += 1
        end
    end
    @debug_l1 "  Wrote $n_written recovery sidecars"
    return n_written
end

function _write_pass1_sidecars_from_main!(file_paths::Vector{String})
    n_written = 0
    for path in file_paths
        main = Arrow.Table(path)
        for col in (:precursor_idx, :scan_idx, :trace_prob_prepass, :trace_prob_infold)
            hasproperty(main, col) ||
                error("Post-qvalue MBR requires column $col in filtered file $path")
        end
        side_df = DataFrame(
            precursor_idx = collect(UInt32.(main.precursor_idx)),
            scan_idx = collect(UInt32.(main.scan_idx)),
            trace_prob_prepass = collect(Float32.(main.trace_prob_prepass)),
            trace_prob_infold = collect(Float32.(main.trace_prob_infold)),
        )
        writeArrow(path * PASS1_SIDECAR_SUFFIX, side_df)
        n_written += 1
    end
    @debug_l1 "  Wrote $n_written temporary post-qvalue Pass-1 sidecars"
    return n_written
end

function merge_mbr_recovery_sidecars_into_main!(
    file_paths::Vector{String};
    cleanup::Bool = true,
)
    n_merged = 0
    for path in file_paths
        rec_path = path * RECOVERY_SIDECAR_SUFFIX
        isfile(rec_path) || continue

        main = DataFrame(Tables.columntable(Arrow.Table(path)))
        rec = Arrow.Table(rec_path)
        n = nrow(main)
        length(rec.precursor_idx) == n ||
            error("Recovery sidecar row-count mismatch at $rec_path")
        @inbounds for i in 1:n
            (main.precursor_idx[i] == rec.precursor_idx[i] &&
             main.scan_idx[i]      == rec.scan_idx[i]) ||
                error("Recovery sidecar misaligned at row $i of $path")
        end

        rec_recovered = collect(Bool.(rec.mbr_recovered))
        old_recovered = hasproperty(main, :mbr_recovered) ?
                        collect(Bool.(main[!, :mbr_recovered])) :
                        falses(n)
        old_candidate = hasproperty(main, :MBR_transfer_candidate) ?
                        collect(Bool.(main[!, :MBR_transfer_candidate])) :
                        falses(n)
        old_prob = hasproperty(main, :mbr_target_decoy_prob) ?
                   collect(Float32.(main[!, :mbr_target_decoy_prob])) :
                   fill(NaN32, n)
        old_qval = hasproperty(main, :ftr_qval_true) ?
                   collect(Float32.(main[!, :ftr_qval_true])) :
                   fill(NaN32, n)
        old_pep = hasproperty(main, :ftr_pep_true) ?
                  collect(Float32.(main[!, :ftr_pep_true])) :
                  fill(NaN32, n)

        rec_candidate = collect(Bool.(rec.MBR_transfer_candidate))
        rec_prob = collect(Float32.(rec.mbr_target_decoy_prob))
        rec_qval = collect(Float32.(rec.ftr_qval_true))
        rec_pep = collect(Float32.(rec.ftr_pep_true))

        main[!, :mbr_recovered] = old_recovered .| rec_recovered
        main[!, :MBR_transfer_candidate] = old_candidate .| rec_candidate
        main[!, :mbr_target_decoy_prob] = ifelse.(rec_recovered, rec_prob, old_prob)
        main[!, :ftr_qval_true] = ifelse.(rec_recovered, rec_qval, old_qval)
        main[!, :ftr_pep_true] = ifelse.(rec_recovered, rec_pep, old_pep)

        rec = nothing
        GC.gc(false)
        writeArrow(path, main)

        if cleanup
            for suffix in (PASS1_SIDECAR_SUFFIX, MBR_SIDECAR_SUFFIX, RECOVERY_SIDECAR_SUFFIX)
                side_path = path * suffix
                isfile(side_path) && safeRm(side_path, nothing; force=true)
            end
        end
        n_merged += 1
    end
    @debug_l1 "  Merged $n_merged post-qvalue recovery sidecars into filtered files"
    return n_merged
end

function _cleanup_mbr_temporary_sidecars!(file_paths::Vector{String})
    for path in file_paths
        for suffix in (PASS1_SIDECAR_SUFFIX, MBR_SIDECAR_SUFFIX, RECOVERY_SIDECAR_SUFFIX)
            side_path = path * suffix
            isfile(side_path) && safeRm(side_path, nothing; force=true)
        end
    end
    return nothing
end

function run_mbr_after_qvalue_filter!(
    candidate_refs::Vector{PSMFileReference},
    donor_refs::Vector{PSMFileReference},
    precursors::LibraryPrecursors,
    fragment_lookup::LibraryFragmentLookup;
    q_value_threshold::Float32 = 0.01f0,
    donor_q_threshold::Float32 = MBR_DONOR_Q_THRESHOLD,
)
    candidate_paths = String[file_path(ref) for ref in candidate_refs if exists(ref)]
    donor_paths = String[file_path(ref) for ref in donor_refs if exists(ref)]
    if isempty(candidate_paths) || isempty(donor_paths)
        return (
            n_files = 0,
            n_donor_files = length(donor_paths),
            n_pass1_sidecars = 0,
            n_recovery_sidecars = 0,
            n_merged = 0,
            n_candidates = 0,
            n_recovered = 0,
        )
    end

    @debug_l1 "MBR Batch F: running after initial run/global q-value filters..."
    sidecar_paths = unique(vcat(candidate_paths, donor_paths))
    n_pass1_sidecars = _write_pass1_sidecars_from_main!(sidecar_paths)

    cf_partner_pools = build_counterfactual_partner_pools(
        sidecar_paths,
        precursors;
        receiver_file_paths = candidate_paths,
    )
    fragment_keys = build_mbr_fragment_annotation_keys(fragment_lookup)

    @debug_l1 "MBR Batch F: computing post-qvalue donor score floor..."
    mbr_prepass = _mbr_prepass_donor_summary(
        donor_paths;
        donor_q_threshold = donor_q_threshold,
    )
    @debug_l1 "  donor prob floor: $(round(mbr_prepass.prob_thresh, digits=4)); " *
              "LOD files: $(length(mbr_prepass.lod_log2_weight_by_file))"
    donor_prob_thresh = mbr_prepass.prob_thresh

    @debug_l1 "MBR Batch F: building post-qvalue donor dict via sweep-1..."
    donor_dict = build_mbr_donor_dict_streaming_with_pass1(
        donor_paths;
        passing_score_floor = donor_prob_thresh,
    )
    @debug_l1 "  donor dict pids: $(length(donor_dict))"

    @debug_l1 "MBR Batch F: writing post-qvalue per-file MBR sidecars..."
    parallel_foreach!(length(candidate_paths)) do chunk
        for f_idx in chunk
            compute_mbr_features_per_file_to_sidecar_with_pass1!(
                candidate_paths[f_idx],
                donor_dict,
                cf_partner_pools,
                fragment_keys,
                passing_score_floor = donor_prob_thresh,
                lod_log2_weight_by_file = mbr_prepass.lod_log2_weight_by_file,
                lod_log2_weight_global = mbr_prepass.lod_log2_weight_global,
            )
        end
    end

    donor_dict = Dict{UInt32, Vector{_MBRDonorEntry}}()
    mbr_prepass = nothing
    cf_partner_pools = nothing
    fragment_keys = nothing
    GC.gc()

    @debug_l1 "MBR Batch F: loading post-qvalue slim FTR DataFrame..."
    psms = load_ftr_slim_dataframe(candidate_paths)
    @debug_l1 "  slim FTR rows: $(nrow(psms))"

    psms[!, :trace_prob] = psms[!, :trace_prob_prepass]
    ftr_summary = apply_mbr_filter_paired!(
        psms;
        alpha = _mbr_recovery_alpha_from_env(),
        q_thresh = q_value_threshold,
        prob_thresh_override = donor_prob_thresh,
    )

    @debug_l1 "MBR Batch F: writing post-qvalue recovery sidecars..."
    n_recovery_sidecars = write_recovery_sidecars(psms, candidate_paths)
    psms = DataFrame()
    GC.gc()

    @debug_l1 "MBR Batch F: merging post-qvalue recovery columns..."
    keep_sidecars = _mbr_keep_temporary_sidecars_from_env()
    n_merged = merge_mbr_recovery_sidecars_into_main!(
        candidate_paths;
        cleanup = !keep_sidecars,
    )
    if keep_sidecars
        @debug_l1 "MBR Batch F: keeping temporary sidecars because PIONEER_MBR_KEEP_SIDECARS is enabled"
    else
        _cleanup_mbr_temporary_sidecars!(sidecar_paths)
    end

    return (
        n_files = length(candidate_paths),
        n_donor_files = length(donor_paths),
        n_pass1_sidecars = n_pass1_sidecars,
        n_recovery_sidecars = n_recovery_sidecars,
        n_merged = n_merged,
        n_candidates = ftr_summary.n_candidates,
        n_recovered = ftr_summary.n_recovered,
    )
end

function run_mbr_after_qvalue_filter!(
    refs::Vector{PSMFileReference},
    precursors::LibraryPrecursors,
    fragment_lookup::LibraryFragmentLookup;
    q_value_threshold::Float32 = 0.01f0,
    donor_q_threshold::Float32 = MBR_DONOR_Q_THRESHOLD,
)
    return run_mbr_after_qvalue_filter!(
        refs,
        refs,
        precursors,
        fragment_lookup;
        q_value_threshold = q_value_threshold,
        donor_q_threshold = donor_q_threshold,
    )
end

# Per-file merge step: read main + Pass-1 sidecar + MBR sidecar + recovery
# sidecar, write back the augmented main file (matching what
# write_scored_psms_to_files! produces in the in-memory path). The :decoy
# column is added defensively (some downstream code reads it). Cleans up
# all three sidecars on success.
function merge_mbr_sidecars_into_main!(file_paths::Vector{String}; cleanup::Bool = true)
    n_merged = 0
    for path in file_paths
        pass1_path = path * PASS1_SIDECAR_SUFFIX
        mbr_path   = path * MBR_SIDECAR_SUFFIX
        rec_path   = path * RECOVERY_SIDECAR_SUFFIX
        all(isfile, (pass1_path, mbr_path, rec_path)) || continue

        main  = Arrow.Table(path)
        pass1 = Arrow.Table(pass1_path)
        mbr   = Arrow.Table(mbr_path)
        rec   = Arrow.Table(rec_path)
        n = length(main.precursor_idx)
        (length(pass1.precursor_idx) == n &&
         length(mbr.precursor_idx)   == n &&
         length(rec.precursor_idx)   == n) ||
            error("Sidecar row-count mismatch at $path")
        @inbounds for i in 1:n
            (main.precursor_idx[i] == pass1.precursor_idx[i] == mbr.precursor_idx[i] == rec.precursor_idx[i] &&
             main.scan_idx[i]      == pass1.scan_idx[i]      == mbr.scan_idx[i]      == rec.scan_idx[i]) ||
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
            main = nothing; pass1 = nothing; mbr = nothing; rec = nothing
            GC.gc(false)
            safeRm(pass1_path, nothing; force=true)
            safeRm(mbr_path, nothing; force=true)
            safeRm(rec_path, nothing; force=true)
        end
        n_merged += 1
    end
    @debug_l1 "  Wrote $n_merged consolidated mbr_outputs sidecars"
    return n_merged
end
