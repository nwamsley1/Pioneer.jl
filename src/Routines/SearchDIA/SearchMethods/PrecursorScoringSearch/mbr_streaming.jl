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
#     Streams each (main + pass1 sidecar) pair to build the top-2 donor dict
#     keyed by precursor_idx. Top-2 (not top-1) covers the same-file fallback
#     in `_donor_for` — top-1 might be the recipient's own file.
#     Memory: ~50k pids × 2 entries × ~40 B ≈ a few MB.
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
# another file. Holds the bare minimum so the donor dict stays small
# (each entry is ~40 bytes; ~50k pids × 2 ≈ a few MB).
struct _MBRDonorEntry
    trace_prob::Float32
    weight::Float32
    log2_intensity_explained::Float32
    irt_residual::Float32  # irt_pred − irt_obs of the donor row
    irt_obs::Float32       # raw observed iRT of the donor row
    log_by_ratio::Float32  # log(b_int+1) − log(y_int+1) of the donor row
    rt_obs::Float32        # literal scan RT (minutes) of the donor row
    ms_file_idx::UInt32
    is_decoy::Bool
    # Donor smoothed-apex fragment vector (top-8, library predicted-intensity
    # rank order — same ranking as the acceptor, so comparable element-wise).
    # Used by the env-gated MBR_frag_smoothed_* features.
    sm::NTuple{8,Float32}
    # Donor FITTED (deconvolved, interference-removed) fragment vector, same
    # rank order. Used by the env-gated MBR_fitfrag_* features.
    ft::NTuple{8,Float32}
end

const _MBR_ZERO8 = ntuple(_ -> 0f0, 8)

# Read row i of 8 smoothed-intensity columns into a stack NTuple (constant
# tuple indices keep this type-stable).
@inline _mbr_sm8(sc::NTuple{8,AbstractVector{Float32}}, i::Int) =
    (Float32(sc[1][i]), Float32(sc[2][i]), Float32(sc[3][i]), Float32(sc[4][i]),
     Float32(sc[5][i]), Float32(sc[6][i]), Float32(sc[7][i]), Float32(sc[8][i]))

# Fragment-vector comparison metrics between an acceptor smoothed 8-vector `p`
# and a donor smoothed 8-vector `q` (library-rank aligned). Returns a 5-tuple;
# -1f0 is the "no comparison" sentinel (matches the other MBR features).
#   cosine     : spectral angle (scale-invariant — magnitude-weighted shape)
#   sum_lr     : log2 ratio of total fragment signal (Σp / Σq)
#   top1_lr    : log2 ratio of the rank-1 fragment (p1 / q1)
#   disp       : std of per-fragment log2 ratios over co-detected ranks — low =
#                all fragments agree on one scale (real); high = inconsistent.
#                Equal-weights every fragment (unlike cosine), so it sees the
#                minor fragments cosine ignores. Sentinel if < 2 co-detected.
#   codetect   : number of ranks where both vectors are > 0 (0..8)
@inline function _mbr_frag_metrics(p::NTuple{8,Float32}, q::NTuple{8,Float32})
    eps = 1f-6
    dot = 0f0; np2 = 0f0; nq2 = 0f0; sp = 0f0; sq = 0f0
    n_cod = 0; m = 0f0; m2 = 0f0
    @inbounds for k in 1:8
        pk = p[k]; qk = q[k]
        dot += pk*qk; np2 += pk^2; nq2 += qk^2; sp += pk; sq += qk
        if pk > 0f0 && qk > 0f0
            r = log2(pk/qk); n_cod += 1; m += r; m2 += r*r
        end
    end
    d = sqrt(np2*nq2)
    cosine  = d  > 0f0 ? Float32(dot/d) : -1f0
    sum_lr  = (sp > 0f0 && sq > 0f0) ? log2((sp+eps)/(sq+eps)) : -1f0
    top1_lr = (p[1] > 0f0 && q[1] > 0f0) ? log2((p[1]+eps)/(q[1]+eps)) : -1f0
    disp    = n_cod >= 2 ? sqrt(max(m2/n_cod - (m/n_cod)^2, 0f0)) : -1f0
    return (cosine, sum_lr, top1_lr, disp, Float32(n_cod))
end

const MBR_SIDECAR_SUFFIX = ".mbr_sidecar.arrow"

# Columns the donor dict reads from each (main + pass1 sidecar) pair.
# Listed once so reader and consumer stay in sync. `is_decoy` is hardcoded
# false on donor entries (the field is unused downstream), so we don't read
# it from the file.
const _MBR_DONOR_COLS = (:precursor_idx, :trace_prob_prepass, :weight,
    :log2_intensity_explained, :irt_pred, :irt_obs, :log_by_ratio_m0, :rt,
    :ms_file_idx,
    # smoothed-apex fragment vector for MBR_frag_smoothed_cosine (env-gated)
    :frag1_smoothed_intensity, :frag2_smoothed_intensity, :frag3_smoothed_intensity,
    :frag4_smoothed_intensity, :frag5_smoothed_intensity, :frag6_smoothed_intensity,
    :frag7_smoothed_intensity, :frag8_smoothed_intensity,
    # fitted (deconvolved, interference-removed) fragment vector for MBR_fitfrag_* (env-gated)
    :fitted_frag1_int, :fitted_frag2_int, :fitted_frag3_int, :fitted_frag4_int,
    :fitted_frag5_int, :fitted_frag6_int, :fitted_frag7_int, :fitted_frag8_int)

# Columns the per-file MBR sidecar emits. precursor_idx + scan_idx are
# redundant with the main file (same positions) but kept as alignment
# check-keys (asserted equal at downstream join time).
const _MBR_SIDECAR_OUT_COLS = (:precursor_idx, :scan_idx,
    :MBR_max_pair_prob_true,        :MBR_max_pair_prob_false,
    :MBR_log2_weight_ratio_true,    :MBR_log2_weight_ratio_false,
    :MBR_log2_explained_ratio_true, :MBR_log2_explained_ratio_false,
    :MBR_best_irt_diff_true,        :MBR_best_irt_diff_false,
    :MBR_is_missing_true,           :MBR_is_missing_false,
    :MBR_log_by_diff_true,          :MBR_log_by_diff_false,
    :MBR_best_rt_diff_true,         :MBR_best_rt_diff_false,
    :MBR_frag_smoothed_cosine_true, :MBR_frag_smoothed_cosine_false,
    :MBR_frag_sum_log_ratio_true,   :MBR_frag_sum_log_ratio_false,
    :MBR_frag_top1_log_ratio_true,  :MBR_frag_top1_log_ratio_false,
    :MBR_frag_ratio_disp_true,      :MBR_frag_ratio_disp_false,
    :MBR_frag_codetect_true,        :MBR_frag_codetect_false,
    :MBR_donor_weight_true,         :MBR_donor_weight_false,
    :MBR_fitfrag_cosine_true,       :MBR_fitfrag_cosine_false,
    :MBR_fitfrag_sum_log_ratio_true,:MBR_fitfrag_sum_log_ratio_false,
    :MBR_fitfrag_top1_log_ratio_true,:MBR_fitfrag_top1_log_ratio_false,
    :MBR_fitfrag_ratio_disp_true,   :MBR_fitfrag_ratio_disp_false,
    :MBR_fitfrag_codetect_true,     :MBR_fitfrag_codetect_false,
)

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
    fidx_c::AbstractVector{UInt32},
    sm_cols::Union{Nothing, NTuple{8, AbstractVector{Float32}}},
    ft_cols::Union{Nothing, NTuple{8, AbstractVector{Float32}}},
    side_path::String,
)
    n = length(pid_c)
    K = 2
    has_logby = logby_c !== nothing
    has_rt    = rt_c    !== nothing
    @inbounds for i in 1:n
        # Sanity-check alignment (cheap: one cmp per row)
        (pid_c[i] == side_pid[i] && main_scn[i] == side_scn[i]) ||
            error("Pass-1 sidecar misaligned at row $i of $side_path")
        pid = pid_c[i]
        sm = sm_cols === nothing ? _MBR_ZERO8 : _mbr_sm8(sm_cols, i)
        ft = ft_cols === nothing ? _MBR_ZERO8 : _mbr_sm8(ft_cols, i)
        e = _MBRDonorEntry(
            score_c[i], w_c[i], Float32(l2ie_c[i]),
            irtp_c[i] - irto_c[i],
            irto_c[i],
            has_logby ? Float32(logby_c[i]) : 0f0,
            has_rt    ? rt_c[i]             : 0f0,
            fidx_c[i], false, sm, ft,
        )
        entries = get!(() -> _MBRDonorEntry[], all_entries, pid)
        if length(entries) < K
            push!(entries, e)
            if length(entries) > 1 && entries[end-1].trace_prob < e.trace_prob
                entries[end-1], entries[end] = entries[end], entries[end-1]
            end
        elseif e.trace_prob > entries[end].trace_prob
            entries[end] = e
            if entries[1].trace_prob < e.trace_prob
                entries[1], entries[end] = entries[end], entries[1]
            end
        end
    end
    return nothing
end

# Streaming version: like build_mbr_donor_dict_streaming, but reads
# `trace_prob_prepass` from the Pass-1 sidecar (rather than expecting it
# in the main file). All other donor columns come from the main file.
# Alignment between main and sidecar is asserted via (:precursor_idx, :scan_idx).
# The 8 smoothed-apex fragment columns as a tuple, or nothing if any is absent.
@inline function _mbr_smoothed_cols(tbl)
    all(c -> hasproperty(tbl, c), SMOOTHED_FRAGMENT_INTENSITY_COLUMNS) || return nothing
    return ntuple(k -> getproperty(tbl, SMOOTHED_FRAGMENT_INTENSITY_COLUMNS[k])::AbstractVector{Float32}, 8)
end

# The 8 fitted (deconvolved) fragment columns as a tuple, or nothing if any is
# absent (e.g. an older fold file produced before they were persisted).
@inline function _mbr_fitted_cols(tbl)
    all(c -> hasproperty(tbl, c), FITTED_FRAGMENT_INTENSITY_COLUMNS) || return nothing
    return ntuple(k -> getproperty(tbl, FITTED_FRAGMENT_INTENSITY_COLUMNS[k])::AbstractVector{Float32}, 8)
end

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
            hasproperty(main, :rt) ? main.rt : nothing,
            main.ms_file_idx,
            _mbr_smoothed_cols(main),
            _mbr_fitted_cols(main),
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

# Inner per-row MBR-feature compute. Extracted from
# compute_mbr_features_per_file_to_sidecar_with_pass1! so Julia specializes
# on the concrete column types (Arrow.Primitive{T} or Vector{T}) and the
# Union{Nothing, ...} branch on logby_v / rt_v becomes a one-time call-site
# decision instead of per-row dispatch.
@inline function _compute_mbr_inner!(
    out_max_t::Vector{Float32}, out_max_f::Vector{Float32},
    out_lw_t::Vector{Float32},  out_lw_f::Vector{Float32},
    out_le_t::Vector{Float32},  out_le_f::Vector{Float32},
    out_ir_t::Vector{Float32},  out_ir_f::Vector{Float32},
    out_miss_t::BitVector,      out_miss_f::BitVector,
    out_log_by_t::Vector{Float32}, out_log_by_f::Vector{Float32},
    out_rt_t::Vector{Float32},   out_rt_f::Vector{Float32},
    pid_v::AbstractVector{UInt32},
    file_v::AbstractVector{UInt32},
    weight_v::AbstractVector{Float32},
    l2ie_v::AbstractVector{Float16},
    irtp_v::AbstractVector{Float32},
    irto_v::AbstractVector{Float32},
    logby_v::Union{Nothing, AbstractVector{Float16}},
    rt_v::Union{Nothing, AbstractVector{Float32}},
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_col::Vector{UInt32},
)
    n = length(pid_v)
    plen = length(partner_col)
    has_logby = logby_v !== nothing
    has_rt    = rt_v    !== nothing
    @inbounds for i in 1:n
        my_file = file_v[i]
        my_pid  = pid_v[i]
        my_partner = my_pid <= UInt32(plen) ? partner_col[Int(my_pid)] : UInt32(0)

        donor_t = _donor_for_pid(donor_dict, my_pid, my_file)
        if donor_t !== nothing
            out_max_t[i] = donor_t.trace_prob
            wi = weight_v[i]
            if donor_t.weight > 0 && wi > 0
                out_lw_t[i] = log2(wi / donor_t.weight)
            end
            out_le_t[i] = Float32(l2ie_v[i]) - donor_t.log2_intensity_explained
            out_ir_t[i] = abs((irtp_v[i] - irto_v[i]) - donor_t.irt_residual)
            if has_logby
                out_log_by_t[i] = Float32(logby_v[i]) - donor_t.log_by_ratio
            end
            if has_rt
                out_rt_t[i] = abs(rt_v[i] - donor_t.rt_obs)
            end
            out_miss_t[i] = false
        end
        if my_partner != UInt32(0)
            donor_f = _donor_for_pid(donor_dict, my_partner, my_file)
            if donor_f !== nothing
                out_max_f[i] = donor_f.trace_prob
                wi = weight_v[i]
                if donor_f.weight > 0 && wi > 0
                    out_lw_f[i] = log2(wi / donor_f.weight)
                end
                out_le_f[i] = Float32(l2ie_v[i]) - donor_f.log2_intensity_explained
                out_ir_f[i] = abs((irtp_v[i] - irto_v[i]) - donor_f.irt_residual)
                if has_logby
                    out_log_by_f[i] = Float32(logby_v[i]) - donor_f.log_by_ratio
                end
                if has_rt
                    out_rt_f[i] = abs(rt_v[i] - donor_f.rt_obs)
                end
                out_miss_f[i] = false
            end
        end
    end
    return nothing
end

# Donor↔acceptor smoothed-apex fragment-vector metrics (env-gated MBR features).
# Kept in its own pass so the established `_compute_mbr_inner!` kernel is
# untouched; when no fragvec feature is gated on these columns simply go unused.
# Fills paired _true (real donor) / _false (counterfactual donor) arrays for the
# five metrics from `_mbr_frag_metrics`.
@inline function _compute_mbr_frag_metrics!(
    cos_t, slr_t, t1_t, disp_t, cod_t, dw_t,
    cos_f, slr_f, t1_f, disp_f, cod_f, dw_f,
    pid_v::AbstractVector{UInt32}, file_v::AbstractVector{UInt32},
    sm_self::NTuple{8, AbstractVector{Float32}},
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_col::Vector{UInt32},
)
    n = length(pid_v)
    plen = length(partner_col)
    @inbounds for i in 1:n
        my_file = file_v[i]
        my_pid  = pid_v[i]
        v_self = _mbr_sm8(sm_self, i)
        self_sig = (v_self[1]+v_self[2]+v_self[3]+v_self[4]+
                    v_self[5]+v_self[6]+v_self[7]+v_self[8]) > 0f0
        donor_t = _donor_for_pid(donor_dict, my_pid, my_file)
        if donor_t !== nothing
            # donor weight: log2-scaled (heavy-tailed); needs no acceptor signal
            dw_t[i] = donor_t.weight > 0f0 ? log2(donor_t.weight) : -1f0
            if self_sig
                (cos_t[i], slr_t[i], t1_t[i], disp_t[i], cod_t[i]) =
                    _mbr_frag_metrics(v_self, donor_t.sm)
            end
        end
        my_partner = my_pid <= UInt32(plen) ? partner_col[Int(my_pid)] : UInt32(0)
        if my_partner != UInt32(0)
            donor_f = _donor_for_pid(donor_dict, my_partner, my_file)
            if donor_f !== nothing
                dw_f[i] = donor_f.weight > 0f0 ? log2(donor_f.weight) : -1f0
                if self_sig
                    (cos_f[i], slr_f[i], t1_f[i], disp_f[i], cod_f[i]) =
                        _mbr_frag_metrics(v_self, donor_f.sm)
                end
            end
        end
    end
    return nothing
end

# Same as _compute_mbr_frag_metrics! but on the FITTED (deconvolved,
# interference-removed) vectors — acceptor fitted vs donor.ft. No donor_weight
# here (already emitted by the smoothed pass). Separate pass keeps the working
# smoothed kernel untouched.
@inline function _compute_mbr_fitted_metrics!(
    cos_t, slr_t, t1_t, disp_t, cod_t,
    cos_f, slr_f, t1_f, disp_f, cod_f,
    pid_v::AbstractVector{UInt32}, file_v::AbstractVector{UInt32},
    ft_self::NTuple{8, AbstractVector{Float32}},
    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
    partner_col::Vector{UInt32},
)
    n = length(pid_v)
    plen = length(partner_col)
    @inbounds for i in 1:n
        v = _mbr_sm8(ft_self, i)
        (v[1]+v[2]+v[3]+v[4]+v[5]+v[6]+v[7]+v[8]) > 0f0 || continue
        my_file = file_v[i]
        my_pid  = pid_v[i]
        donor_t = _donor_for_pid(donor_dict, my_pid, my_file)
        if donor_t !== nothing
            (cos_t[i], slr_t[i], t1_t[i], disp_t[i], cod_t[i]) =
                _mbr_frag_metrics(v, donor_t.ft)
        end
        my_partner = my_pid <= UInt32(plen) ? partner_col[Int(my_pid)] : UInt32(0)
        if my_partner != UInt32(0)
            donor_f = _donor_for_pid(donor_dict, my_partner, my_file)
            if donor_f !== nothing
                (cos_f[i], slr_f[i], t1_f[i], disp_f[i], cod_f[i]) =
                    _mbr_frag_metrics(v, donor_f.ft)
            end
        end
    end
    return nothing
end

# Per-file MBR feature compute + sidecar write, reading trace_prob_prepass
# from the Pass-1 sidecar. Mirrors compute_mbr_features_per_file_to_sidecar!
# but for the streaming-with-Pass-1 flow.
function compute_mbr_features_per_file_to_sidecar_with_pass1!(
        main_path::AbstractString,
        donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
        partner_col::Vector{UInt32})
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

    @inbounds for i in 1:n
        (pid_v[i] == pass1.precursor_idx[i] && scan_v[i] == pass1.scan_idx[i]) ||
            error("Pass-1 sidecar misaligned at row $i of $pass1_path")
    end

    out_max_t = fill(-1f0, n); out_max_f = fill(-1f0, n)
    out_lw_t  = fill(-1f0, n); out_lw_f  = fill(-1f0, n)
    out_le_t  = fill(-1f0, n); out_le_f  = fill(-1f0, n)
    out_ir_t  = fill(-1f0, n); out_ir_f  = fill(-1f0, n)
    out_miss_t = trues(n);     out_miss_f = trues(n)
    out_log_by_t = fill(-1f0, n); out_log_by_f = fill(-1f0, n)
    out_rt_t = fill(-1f0, n); out_rt_f = fill(-1f0, n)
    out_cos_t  = fill(-1f0, n); out_cos_f  = fill(-1f0, n)
    out_slr_t  = fill(-1f0, n); out_slr_f  = fill(-1f0, n)
    out_t1_t   = fill(-1f0, n); out_t1_f   = fill(-1f0, n)
    out_disp_t = fill(-1f0, n); out_disp_f = fill(-1f0, n)
    out_cod_t  = fill(-1f0, n); out_cod_f  = fill(-1f0, n)
    out_dw_t   = fill(-1f0, n); out_dw_f   = fill(-1f0, n)
    out_fcos_t  = fill(-1f0, n); out_fcos_f  = fill(-1f0, n)
    out_fslr_t  = fill(-1f0, n); out_fslr_f  = fill(-1f0, n)
    out_ft1_t   = fill(-1f0, n); out_ft1_f   = fill(-1f0, n)
    out_fdisp_t = fill(-1f0, n); out_fdisp_f = fill(-1f0, n)
    out_fcod_t  = fill(-1f0, n); out_fcod_f  = fill(-1f0, n)

    _compute_mbr_inner!(
        out_max_t, out_max_f, out_lw_t, out_lw_f,
        out_le_t, out_le_f, out_ir_t, out_ir_f,
        out_miss_t, out_miss_f,
        out_log_by_t, out_log_by_f, out_rt_t, out_rt_f,
        pid_v, file_v, weight_v, l2ie_v, irtp_v, irto_v,
        logby_v, rt_v,
        donor_dict, partner_col,
    )

    sm_self = _mbr_smoothed_cols(main)
    sm_self !== nothing && _compute_mbr_frag_metrics!(
        out_cos_t, out_slr_t, out_t1_t, out_disp_t, out_cod_t, out_dw_t,
        out_cos_f, out_slr_f, out_t1_f, out_disp_f, out_cod_f, out_dw_f,
        pid_v, file_v, sm_self, donor_dict, partner_col,
    )

    ft_self = _mbr_fitted_cols(main)
    ft_self !== nothing && _compute_mbr_fitted_metrics!(
        out_fcos_t, out_fslr_t, out_ft1_t, out_fdisp_t, out_fcod_t,
        out_fcos_f, out_fslr_f, out_ft1_f, out_fdisp_f, out_fcod_f,
        pid_v, file_v, ft_self, donor_dict, partner_col,
    )

    side_path = main_path * MBR_SIDECAR_SUFFIX
    side_df = DataFrame(
        precursor_idx                  = collect(UInt32.(pid_v)),
        scan_idx                       = collect(UInt32.(scan_v)),
        MBR_max_pair_prob_true         = out_max_t,
        MBR_max_pair_prob_false        = out_max_f,
        MBR_log2_weight_ratio_true     = out_lw_t,
        MBR_log2_weight_ratio_false    = out_lw_f,
        MBR_log2_explained_ratio_true  = out_le_t,
        MBR_log2_explained_ratio_false = out_le_f,
        MBR_best_irt_diff_true         = out_ir_t,
        MBR_best_irt_diff_false        = out_ir_f,
        MBR_is_missing_true            = out_miss_t,
        MBR_is_missing_false           = out_miss_f,
        MBR_log_by_diff_true           = out_log_by_t,
        MBR_log_by_diff_false          = out_log_by_f,
        MBR_best_rt_diff_true          = out_rt_t,
        MBR_best_rt_diff_false         = out_rt_f,
        MBR_frag_smoothed_cosine_true  = out_cos_t,
        MBR_frag_smoothed_cosine_false = out_cos_f,
        MBR_frag_sum_log_ratio_true    = out_slr_t,
        MBR_frag_sum_log_ratio_false   = out_slr_f,
        MBR_frag_top1_log_ratio_true   = out_t1_t,
        MBR_frag_top1_log_ratio_false  = out_t1_f,
        MBR_frag_ratio_disp_true       = out_disp_t,
        MBR_frag_ratio_disp_false      = out_disp_f,
        MBR_frag_codetect_true         = out_cod_t,
        MBR_frag_codetect_false        = out_cod_f,
        MBR_donor_weight_true          = out_dw_t,
        MBR_donor_weight_false         = out_dw_f,
        MBR_fitfrag_cosine_true        = out_fcos_t,
        MBR_fitfrag_cosine_false       = out_fcos_f,
        MBR_fitfrag_sum_log_ratio_true = out_fslr_t,
        MBR_fitfrag_sum_log_ratio_false= out_fslr_f,
        MBR_fitfrag_top1_log_ratio_true = out_ft1_t,
        MBR_fitfrag_top1_log_ratio_false= out_ft1_f,
        MBR_fitfrag_ratio_disp_true    = out_fdisp_t,
        MBR_fitfrag_ratio_disp_false   = out_fdisp_f,
        MBR_fitfrag_codetect_true      = out_fcod_t,
        MBR_fitfrag_codetect_false     = out_fcod_f,
    )
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

