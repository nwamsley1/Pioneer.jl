# MBR streaming pipeline — Phase 1: per-file MBR features as sidecars.
#
# Replaces the in-memory `best_psms → compute_mbr_features_dual! →
# write_scored_psms_to_files!` rewrite chain with two streaming sweeps:
#
#   Sweep 1: build_mbr_donor_dict_streaming(refs)
#     - Streams per-file `Arrow.Table` mmap, projects 14 columns, builds a
#       Dict{UInt32, Vector{_MBRDonorEntry}} with the top-2 donors per
#       precursor_idx. Top-2 (not top-1) covers the same-file fallback in
#       _donor_for (top-1 might be from the row's own file).
#     - Memory: ~50k pids × 2 × 60 B ≈ ~6 MB on Olsen 23-file scale.
#
#   Sweep 2: compute_mbr_features_per_file_to_sidecar!(ref, donor_dict, ...)
#     - Loads ONE file's PSMs at a time (~50k rows), computes 14 MBR_* columns
#       using `donor_dict`, writes a per-file sidecar Arrow with
#       (:precursor_idx, :scan_idx, MBR_*) in main-file row order.
#     - Main file is never rewritten.
#
# The alignment invariant for sidecars: row N of the sidecar ↔ row N of the
# main file (positional join). `:precursor_idx` and `:scan_idx` are written
# into the sidecar redundantly so a downstream consumer can ASSERT alignment
# at read time and fail loud if anything reordered between writes.
#
# Activated via env var PIONEER_MBR_STREAMING=1. When unset, the existing
# in-memory `compute_mbr_features_dual!` path is used.

const MBR_SIDECAR_SUFFIX = ".mbr_sidecar.arrow"

# Columns required to materialise a `_MBRDonorEntry` from a row of the
# per-file Arrow table.  Listed once so reader and consumer stay in sync.
# Note: the dual-MBR variant hardcodes `is_decoy=false` on donor entries
# (the field is unused downstream), so we don't read :is_decoy from the file.
const _MBR_DONOR_COLS = (:precursor_idx, :trace_prob_prepass, :weight,
    :log2_intensity_explained, :irt_pred, :irt_obs, :log_by_ratio_m0, :rt,
    :frag1_int, :frag2_int, :frag3_int, :frag4_int, :frag5_int, :frag6_int,
    :ms_file_idx)

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
    :MBR_frag_pattern_cosine_true,  :MBR_frag_pattern_cosine_false,
    :MBR_frag_pattern_scribe_true,  :MBR_frag_pattern_scribe_false,
    :MBR_best_rt_diff_true,         :MBR_best_rt_diff_false,
    :MBR_frag_rank_corr_true,       :MBR_frag_rank_corr_false,
)

# Build the top-2 donor entries per precursor_idx by streaming `refs`. Reads
# only the donor columns from each file; never materialises the full
# DataFrame.  Returns a Dict{UInt32, Vector{_MBRDonorEntry}}.
function build_mbr_donor_dict_streaming(refs::Vector{<:Any})
    K = 2  # top-2 covers the same-file fallback inside `_donor_for`
    all_entries = Dict{UInt32, Vector{_MBRDonorEntry}}()
    for ref in refs
        path = ref isa AbstractString ? ref : file_path(ref)
        tbl = Arrow.Table(path)
        n = length(tbl.precursor_idx)
        n == 0 && continue
        # Project the columns we need (mmap views — no copy).
        pid_c    = tbl.precursor_idx
        score_c  = tbl.trace_prob_prepass
        w_c      = tbl.weight
        l2ie_c   = tbl.log2_intensity_explained
        irtp_c   = tbl.irt_pred
        irto_c   = tbl.irt_obs
        logby_c  = hasproperty(tbl, :log_by_ratio_m0) ? tbl.log_by_ratio_m0 : nothing
        rt_c     = hasproperty(tbl, :rt) ? tbl.rt : nothing
        f1_c     = tbl.frag1_int
        f2_c     = tbl.frag2_int
        f3_c     = tbl.frag3_int
        f4_c     = tbl.frag4_int
        f5_c     = tbl.frag5_int
        f6_c     = tbl.frag6_int
        fidx_c   = tbl.ms_file_idx
        @inbounds for i in 1:n
            pid = UInt32(pid_c[i])
            e = _MBRDonorEntry(
                Float32(score_c[i]), Float32(w_c[i]), Float32(l2ie_c[i]),
                Float32(irtp_c[i] - irto_c[i]),                # irt_residual
                Float32(irto_c[i]),
                logby_c === nothing ? 0f0 : Float32(logby_c[i]),
                rt_c    === nothing ? 0f0 : Float32(rt_c[i]),
                Float32(f1_c[i]), Float32(f2_c[i]), Float32(f3_c[i]),
                Float32(f4_c[i]), Float32(f5_c[i]), Float32(f6_c[i]),
                UInt32(fidx_c[i]), false,    # is_decoy unused in dual variant
            )
            entries = get!(() -> _MBRDonorEntry[], all_entries, pid)
            # Streaming top-K insertion: keep up to K entries sorted by
            # trace_prob desc.  Cheaper than push+sort+trim per row when K is
            # small.
            if length(entries) < K
                push!(entries, e)
                if length(entries) > 1 && entries[end-1].trace_prob < e.trace_prob
                    # 2-element bubble: ensure descending
                    entries[end-1], entries[end] = entries[end], entries[end-1]
                end
            elseif e.trace_prob > entries[end].trace_prob
                # Replace the current 2nd-best
                entries[end] = e
                if entries[1].trace_prob < e.trace_prob
                    entries[1], entries[end] = entries[end], entries[1]
                end
            end
        end
    end
    return all_entries
end

# Compute the 14 MBR_* feature columns for ONE file by reading the file's
# main Arrow row-by-row and looking up donors in `donor_dict`. Writes the
# per-file sidecar (`<main>.mbr_sidecar.arrow`) and returns the path.
#
# `partner_col` is the spectral-library lookup
# (UInt32(0) → no counterfactual). Must be indexed by precursor_idx.
function compute_mbr_features_per_file_to_sidecar!(ref::Any,
                                                    donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
                                                    partner_col::AbstractVector)
    main_path = ref isa AbstractString ? ref : file_path(ref)
    tbl  = Arrow.Table(main_path)
    n    = length(tbl.precursor_idx)
    pid_v   = tbl.precursor_idx
    scan_v  = tbl.scan_idx
    file_v  = tbl.ms_file_idx
    score_v = tbl.trace_prob_prepass
    weight_v= tbl.weight
    l2ie_v  = tbl.log2_intensity_explained
    irtp_v  = tbl.irt_pred
    irto_v  = tbl.irt_obs
    logby_v = hasproperty(tbl, :log_by_ratio_m0) ? tbl.log_by_ratio_m0 : nothing
    rt_v    = hasproperty(tbl, :rt) ? tbl.rt : nothing
    f1_v    = tbl.frag1_int; f2_v = tbl.frag2_int; f3_v = tbl.frag3_int
    f4_v    = tbl.frag4_int; f5_v = tbl.frag5_int; f6_v = tbl.frag6_int

    out_max_t = fill(-1f0, n); out_max_f = fill(-1f0, n)
    out_lw_t  = fill(-1f0, n); out_lw_f  = fill(-1f0, n)
    out_le_t  = fill(-1f0, n); out_le_f  = fill(-1f0, n)
    out_ir_t  = fill(-1f0, n); out_ir_f  = fill(-1f0, n)
    out_miss_t = trues(n);     out_miss_f = trues(n)
    out_log_by_t = fill(-1f0, n); out_log_by_f = fill(-1f0, n)
    out_frag_cos_t = fill(-1f0, n); out_frag_cos_f = fill(-1f0, n)
    out_frag_scr_t = fill(-1f0, n); out_frag_scr_f = fill(-1f0, n)
    out_rt_t = fill(-1f0, n); out_rt_f = fill(-1f0, n)
    out_frag_kt_t = fill(-2f0, n); out_frag_kt_f = fill(-2f0, n)

    has_logby = logby_v !== nothing
    has_rt    = rt_v    !== nothing

    # Inline fragment-pattern helpers (kept here for streaming-path use).
    @inline _norm_l1(a, b, c, d, e, f) = begin
        s = a + b + c + d + e + f
        s > 0f0 ? (a/s, b/s, c/s, d/s, e/s, f/s) : (0f0, 0f0, 0f0, 0f0, 0f0, 0f0)
    end
    @inline _cosine6(p::NTuple{6,Float32}, q::NTuple{6,Float32}) = begin
        dot = p[1]*q[1]+p[2]*q[2]+p[3]*q[3]+p[4]*q[4]+p[5]*q[5]+p[6]*q[6]
        np2 = p[1]^2+p[2]^2+p[3]^2+p[4]^2+p[5]^2+p[6]^2
        nq2 = q[1]^2+q[2]^2+q[3]^2+q[4]^2+q[5]^2+q[6]^2
        d = sqrt(np2*nq2); d > 0f0 ? Float32(dot/d) : 0f0
    end
    @inline _scribe6(p::NTuple{6,Float32}, q::NTuple{6,Float32}) = begin
        ad = abs(p[1]-q[1])+abs(p[2]-q[2])+abs(p[3]-q[3])+abs(p[4]-q[4])+abs(p[5]-q[5])+abs(p[6]-q[6])
        Float32(-log2(ad/2f0 + 1f-6))
    end
    @inline _kendall6(p::NTuple{6,Float32}, q::NTuple{6,Float32}) = begin
        nc = 0; nd = 0
        @inbounds for i in 1:5, j in (i+1):6
            di = p[i] - p[j]; dj = q[i] - q[j]
            si = di > 0f0 ? 1 : (di < 0f0 ? -1 : 0)
            sj = dj > 0f0 ? 1 : (dj < 0f0 ? -1 : 0)
            (si != 0 && sj != 0) && (si == sj ? (nc += 1) : (nd += 1))
        end
        nt = nc + nd; nt > 0 ? Float32((nc - nd) / nt) : -2f0
    end
    @inline function _donor_for(pid::UInt32, my_file::UInt32)
        entries = get(donor_dict, pid, nothing)
        entries === nothing && return nothing
        @inbounds for e in entries
            e.ms_file_idx != my_file && return e
        end
        return nothing
    end

    @inbounds for i in 1:n
        my_file = UInt32(file_v[i])
        my_pid  = UInt32(pid_v[i])
        my_pidx = Int(my_pid)
        my_partner = my_pidx <= length(partner_col) ?
                      (ismissing(partner_col[my_pidx]) ? UInt32(0) : UInt32(partner_col[my_pidx])) :
                      UInt32(0)

        v_self_norm = _norm_l1(Float32(f1_v[i]), Float32(f2_v[i]), Float32(f3_v[i]),
                               Float32(f4_v[i]), Float32(f5_v[i]), Float32(f6_v[i]))
        v_self_has_signal = (v_self_norm[1]+v_self_norm[2]+v_self_norm[3]+
                             v_self_norm[4]+v_self_norm[5]+v_self_norm[6]) > 0f0

        # True donor (this precursor in another file)
        donor_t = _donor_for(my_pid, my_file)
        if donor_t !== nothing
            out_max_t[i] = donor_t.trace_prob
            if donor_t.weight > 0 && weight_v[i] > 0
                out_lw_t[i] = log2(Float32(weight_v[i]) / donor_t.weight)
            end
            out_le_t[i] = Float32(l2ie_v[i]) - donor_t.log2_intensity_explained
            out_ir_t[i] = abs(Float32(irtp_v[i] - irto_v[i]) - donor_t.irt_residual)
            if has_logby
                out_log_by_t[i] = Float32(logby_v[i]) - donor_t.log_by_ratio
            end
            if has_rt
                out_rt_t[i] = abs(Float32(rt_v[i]) - donor_t.rt_obs)
            end
            if v_self_has_signal
                v_donor_norm = _norm_l1(donor_t.frag1_int, donor_t.frag2_int, donor_t.frag3_int,
                                         donor_t.frag4_int, donor_t.frag5_int, donor_t.frag6_int)
                v_donor_has_signal = (v_donor_norm[1]+v_donor_norm[2]+v_donor_norm[3]+
                                      v_donor_norm[4]+v_donor_norm[5]+v_donor_norm[6]) > 0f0
                if v_donor_has_signal
                    out_frag_cos_t[i] = _cosine6(v_self_norm, v_donor_norm)
                    out_frag_scr_t[i] = _scribe6(v_self_norm, v_donor_norm)
                    out_frag_kt_t[i]  = _kendall6(v_self_norm, v_donor_norm)
                end
            end
            out_miss_t[i] = false
        end

        # False donor (counterfactual partner precursor in any file ≠ mine)
        if my_partner != UInt32(0)
            donor_f = _donor_for(my_partner, my_file)
            if donor_f !== nothing
                out_max_f[i] = donor_f.trace_prob
                if donor_f.weight > 0 && weight_v[i] > 0
                    out_lw_f[i] = log2(Float32(weight_v[i]) / donor_f.weight)
                end
                out_le_f[i] = Float32(l2ie_v[i]) - donor_f.log2_intensity_explained
                out_ir_f[i] = abs(Float32(irtp_v[i] - irto_v[i]) - donor_f.irt_residual)
                if has_logby
                    out_log_by_f[i] = Float32(logby_v[i]) - donor_f.log_by_ratio
                end
                if has_rt
                    out_rt_f[i] = abs(Float32(rt_v[i]) - donor_f.rt_obs)
                end
                if v_self_has_signal
                    v_donor_norm = _norm_l1(donor_f.frag1_int, donor_f.frag2_int, donor_f.frag3_int,
                                             donor_f.frag4_int, donor_f.frag5_int, donor_f.frag6_int)
                    v_donor_has_signal = (v_donor_norm[1]+v_donor_norm[2]+v_donor_norm[3]+
                                          v_donor_norm[4]+v_donor_norm[5]+v_donor_norm[6]) > 0f0
                    if v_donor_has_signal
                        out_frag_cos_f[i] = _cosine6(v_self_norm, v_donor_norm)
                        out_frag_scr_f[i] = _scribe6(v_self_norm, v_donor_norm)
                        out_frag_kt_f[i]  = _kendall6(v_self_norm, v_donor_norm)
                    end
                end
                out_miss_f[i] = false
            end
        end
    end

    sidecar_path = main_path * MBR_SIDECAR_SUFFIX
    sidecar_df = DataFrame(
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
        MBR_frag_pattern_cosine_true   = out_frag_cos_t,
        MBR_frag_pattern_cosine_false  = out_frag_cos_f,
        MBR_frag_pattern_scribe_true   = out_frag_scr_t,
        MBR_frag_pattern_scribe_false  = out_frag_scr_f,
        MBR_best_rt_diff_true          = out_rt_t,
        MBR_best_rt_diff_false         = out_rt_f,
        MBR_frag_rank_corr_true        = out_frag_kt_t,
        MBR_frag_rank_corr_false       = out_frag_kt_f,
    )
    writeArrow(sidecar_path, sidecar_df)
    return sidecar_path
end

# Load a per-file MBR sidecar and verify alignment with the main file via
# (:precursor_idx, :scan_idx). Returns a DataFrame with the MBR_* columns
# only (the alignment keys are dropped after verification).
function load_mbr_sidecar_aligned(main_path::AbstractString,
                                   sidecar_path::AbstractString=main_path*MBR_SIDECAR_SUFFIX)
    main = Arrow.Table(main_path)
    side = Arrow.Table(sidecar_path)
    n_main = length(main.precursor_idx)
    n_side = length(side.precursor_idx)
    n_main == n_side ||
        error("MBR sidecar row-count mismatch: main=$n_main side=$n_side at $sidecar_path")
    main_pid = main.precursor_idx; side_pid = side.precursor_idx
    main_scn = main.scan_idx;      side_scn = side.scan_idx
    @inbounds for i in 1:n_main
        if main_pid[i] != side_pid[i] || main_scn[i] != side_scn[i]
            error("MBR sidecar misaligned at row $i (precursor_idx/scan_idx mismatch) at $sidecar_path")
        end
    end
    df = DataFrame()
    for c in _MBR_SIDECAR_OUT_COLS
        c === :precursor_idx && continue
        c === :scan_idx      && continue
        df[!, c] = collect(Tables.getcolumn(side, c))
    end
    return df
end

# Suffix conventions for the three sidecar types used by the streaming MBR
# pipeline. All sidecars carry (:precursor_idx, :scan_idx) as alignment keys.
const PASS1_SIDECAR_SUFFIX    = ".pass1_sidecar.arrow"
const RECOVERY_SIDECAR_SUFFIX = ".recovery_sidecar.arrow"

# Distribute Pass-1 LightGBM scores (trace_prob_prepass, trace_prob_infold)
# back to per-(file, fold) Pass-1 sidecars. best_psms is grouped by
# (ms_file_idx, cv_fold) — within-group row order matches each main file's
# row order because load_psms_for_lightgbm preserves Arrow.Table row order.
function write_pass1_score_sidecars!(best_psms::DataFrame, file_paths::Vector{String})
    # Build (ms_file_idx, cv_fold) → main path map by peeking at each file.
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
    if is_fold_split && hasproperty(best_psms, :cv_fold)
        for (key, gpsms) in pairs(groupby(best_psms, [:ms_file_idx, :cv_fold]))
            lookup_key = (UInt32(key[:ms_file_idx]), UInt8(key[:cv_fold]))
            haskey(key_to_path, lookup_key) || continue
            main_path = key_to_path[lookup_key]
            side_path = main_path * PASS1_SIDECAR_SUFFIX
            side_df = DataFrame(
                precursor_idx       = collect(UInt32.(gpsms[!, :precursor_idx])),
                scan_idx            = collect(UInt32.(gpsms[!, :scan_idx])),
                trace_prob_prepass  = Float32.(gpsms[!, :trace_prob_prepass]),
                trace_prob_infold   = hasproperty(gpsms, :trace_prob_infold) ?
                                       Float32.(gpsms[!, :trace_prob_infold]) :
                                       fill(NaN32, nrow(gpsms)),
            )
            writeArrow(side_path, side_df)
            n_written += 1
        end
    else
        for (key, gpsms) in pairs(groupby(best_psms, :ms_file_idx))
            lookup_key = UInt32(key[:ms_file_idx])
            haskey(key_to_path, lookup_key) || continue
            main_path = key_to_path[lookup_key]
            side_path = main_path * PASS1_SIDECAR_SUFFIX
            side_df = DataFrame(
                precursor_idx       = collect(UInt32.(gpsms[!, :precursor_idx])),
                scan_idx            = collect(UInt32.(gpsms[!, :scan_idx])),
                trace_prob_prepass  = Float32.(gpsms[!, :trace_prob_prepass]),
                trace_prob_infold   = hasproperty(gpsms, :trace_prob_infold) ?
                                       Float32.(gpsms[!, :trace_prob_infold]) :
                                       fill(NaN32, nrow(gpsms)),
            )
            writeArrow(side_path, side_df)
            n_written += 1
        end
    end
    @user_info "  Wrote $n_written Pass-1 score sidecars"
    return n_written
end

# Streaming version: like build_mbr_donor_dict_streaming, but reads
# `trace_prob_prepass` from the Pass-1 sidecar (rather than expecting it
# in the main file). All other donor columns come from the main file.
# Alignment between main and sidecar is asserted via (:precursor_idx, :scan_idx).
function build_mbr_donor_dict_streaming_with_pass1(refs::Vector{<:Any})
    K = 2
    all_entries = Dict{UInt32, Vector{_MBRDonorEntry}}()
    for ref in refs
        main_path = ref isa AbstractString ? ref : file_path(ref)
        side_path = main_path * PASS1_SIDECAR_SUFFIX
        isfile(side_path) || error("Missing Pass-1 sidecar at $side_path")
        main = Arrow.Table(main_path)
        side = Arrow.Table(side_path)
        n = length(main.precursor_idx)
        n == length(side.precursor_idx) ||
            error("Pass-1 sidecar row count mismatch at $side_path")
        pid_c   = main.precursor_idx
        side_pid = side.precursor_idx
        side_scn = side.scan_idx
        main_scn = main.scan_idx
        score_c = side.trace_prob_prepass    # from sidecar
        w_c     = main.weight
        l2ie_c  = main.log2_intensity_explained
        irtp_c  = main.irt_pred
        irto_c  = main.irt_obs
        logby_c = hasproperty(main, :log_by_ratio_m0) ? main.log_by_ratio_m0 : nothing
        rt_c    = hasproperty(main, :rt)             ? main.rt             : nothing
        f1_c    = main.frag1_int; f2_c = main.frag2_int; f3_c = main.frag3_int
        f4_c    = main.frag4_int; f5_c = main.frag5_int; f6_c = main.frag6_int
        fidx_c  = main.ms_file_idx
        @inbounds for i in 1:n
            # Sanity-check alignment (cheap: one cmp per row)
            (pid_c[i] == side_pid[i] && main_scn[i] == side_scn[i]) ||
                error("Pass-1 sidecar misaligned at row $i of $side_path")
            pid = UInt32(pid_c[i])
            e = _MBRDonorEntry(
                Float32(score_c[i]), Float32(w_c[i]), Float32(l2ie_c[i]),
                Float32(irtp_c[i] - irto_c[i]),
                Float32(irto_c[i]),
                logby_c === nothing ? 0f0 : Float32(logby_c[i]),
                rt_c    === nothing ? 0f0 : Float32(rt_c[i]),
                Float32(f1_c[i]), Float32(f2_c[i]), Float32(f3_c[i]),
                Float32(f4_c[i]), Float32(f5_c[i]), Float32(f6_c[i]),
                UInt32(fidx_c[i]), false,
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
    end
    return all_entries
end

# Per-file MBR feature compute + sidecar write, reading trace_prob_prepass
# from the Pass-1 sidecar. Mirrors compute_mbr_features_per_file_to_sidecar!
# but for the streaming-with-Pass-1 flow.
function compute_mbr_features_per_file_to_sidecar_with_pass1!(ref::Any,
        donor_dict::Dict{UInt32, Vector{_MBRDonorEntry}},
        partner_col::AbstractVector)
    main_path = ref isa AbstractString ? ref : file_path(ref)
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
    f1_v    = main.frag1_int; f2_v = main.frag2_int; f3_v = main.frag3_int
    f4_v    = main.frag4_int; f5_v = main.frag5_int; f6_v = main.frag6_int

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
    out_frag_cos_t = fill(-1f0, n); out_frag_cos_f = fill(-1f0, n)
    out_frag_scr_t = fill(-1f0, n); out_frag_scr_f = fill(-1f0, n)
    out_rt_t = fill(-1f0, n); out_rt_f = fill(-1f0, n)
    out_frag_kt_t = fill(-2f0, n); out_frag_kt_f = fill(-2f0, n)
    has_logby = logby_v !== nothing
    has_rt    = rt_v    !== nothing

    @inline _norm_l1(a, b, c, d, e, f) = begin
        s = a + b + c + d + e + f
        s > 0f0 ? (a/s, b/s, c/s, d/s, e/s, f/s) : (0f0, 0f0, 0f0, 0f0, 0f0, 0f0)
    end
    @inline _cosine6(p::NTuple{6,Float32}, q::NTuple{6,Float32}) = begin
        dot = p[1]*q[1]+p[2]*q[2]+p[3]*q[3]+p[4]*q[4]+p[5]*q[5]+p[6]*q[6]
        np2 = p[1]^2+p[2]^2+p[3]^2+p[4]^2+p[5]^2+p[6]^2
        nq2 = q[1]^2+q[2]^2+q[3]^2+q[4]^2+q[5]^2+q[6]^2
        d = sqrt(np2*nq2); d > 0f0 ? Float32(dot/d) : 0f0
    end
    @inline _scribe6(p::NTuple{6,Float32}, q::NTuple{6,Float32}) = begin
        ad = abs(p[1]-q[1])+abs(p[2]-q[2])+abs(p[3]-q[3])+abs(p[4]-q[4])+abs(p[5]-q[5])+abs(p[6]-q[6])
        Float32(-log2(ad/2f0 + 1f-6))
    end
    @inline _kendall6(p::NTuple{6,Float32}, q::NTuple{6,Float32}) = begin
        nc = 0; nd = 0
        @inbounds for i in 1:5, j in (i+1):6
            di = p[i] - p[j]; dj = q[i] - q[j]
            si = di > 0f0 ? 1 : (di < 0f0 ? -1 : 0)
            sj = dj > 0f0 ? 1 : (dj < 0f0 ? -1 : 0)
            (si != 0 && sj != 0) && (si == sj ? (nc += 1) : (nd += 1))
        end
        nt = nc + nd; nt > 0 ? Float32((nc - nd) / nt) : -2f0
    end
    @inline function _donor_for(pid::UInt32, my_file::UInt32)
        entries = get(donor_dict, pid, nothing)
        entries === nothing && return nothing
        @inbounds for e in entries
            e.ms_file_idx != my_file && return e
        end
        return nothing
    end

    @inbounds for i in 1:n
        my_file = UInt32(file_v[i])
        my_pid  = UInt32(pid_v[i])
        my_partner = my_pid <= UInt32(length(partner_col)) ?
                      (ismissing(partner_col[Int(my_pid)]) ? UInt32(0) : UInt32(partner_col[Int(my_pid)])) :
                      UInt32(0)

        v_self_norm = _norm_l1(Float32(f1_v[i]), Float32(f2_v[i]), Float32(f3_v[i]),
                               Float32(f4_v[i]), Float32(f5_v[i]), Float32(f6_v[i]))
        v_self_has_signal = (v_self_norm[1]+v_self_norm[2]+v_self_norm[3]+
                             v_self_norm[4]+v_self_norm[5]+v_self_norm[6]) > 0f0

        donor_t = _donor_for(my_pid, my_file)
        if donor_t !== nothing
            out_max_t[i] = donor_t.trace_prob
            if donor_t.weight > 0 && weight_v[i] > 0
                out_lw_t[i] = log2(Float32(weight_v[i]) / donor_t.weight)
            end
            out_le_t[i] = Float32(l2ie_v[i]) - donor_t.log2_intensity_explained
            out_ir_t[i] = abs(Float32(irtp_v[i] - irto_v[i]) - donor_t.irt_residual)
            has_logby && (out_log_by_t[i] = Float32(logby_v[i]) - donor_t.log_by_ratio)
            has_rt    && (out_rt_t[i]     = abs(Float32(rt_v[i]) - donor_t.rt_obs))
            if v_self_has_signal
                v_donor_norm = _norm_l1(donor_t.frag1_int, donor_t.frag2_int, donor_t.frag3_int,
                                         donor_t.frag4_int, donor_t.frag5_int, donor_t.frag6_int)
                v_donor_has_signal = (v_donor_norm[1]+v_donor_norm[2]+v_donor_norm[3]+
                                      v_donor_norm[4]+v_donor_norm[5]+v_donor_norm[6]) > 0f0
                if v_donor_has_signal
                    out_frag_cos_t[i] = _cosine6(v_self_norm, v_donor_norm)
                    out_frag_scr_t[i] = _scribe6(v_self_norm, v_donor_norm)
                    out_frag_kt_t[i]  = _kendall6(v_self_norm, v_donor_norm)
                end
            end
            out_miss_t[i] = false
        end
        if my_partner != UInt32(0)
            donor_f = _donor_for(my_partner, my_file)
            if donor_f !== nothing
                out_max_f[i] = donor_f.trace_prob
                if donor_f.weight > 0 && weight_v[i] > 0
                    out_lw_f[i] = log2(Float32(weight_v[i]) / donor_f.weight)
                end
                out_le_f[i] = Float32(l2ie_v[i]) - donor_f.log2_intensity_explained
                out_ir_f[i] = abs(Float32(irtp_v[i] - irto_v[i]) - donor_f.irt_residual)
                has_logby && (out_log_by_f[i] = Float32(logby_v[i]) - donor_f.log_by_ratio)
                has_rt    && (out_rt_f[i]     = abs(Float32(rt_v[i]) - donor_f.rt_obs))
                if v_self_has_signal
                    v_donor_norm = _norm_l1(donor_f.frag1_int, donor_f.frag2_int, donor_f.frag3_int,
                                             donor_f.frag4_int, donor_f.frag5_int, donor_f.frag6_int)
                    v_donor_has_signal = (v_donor_norm[1]+v_donor_norm[2]+v_donor_norm[3]+
                                          v_donor_norm[4]+v_donor_norm[5]+v_donor_norm[6]) > 0f0
                    if v_donor_has_signal
                        out_frag_cos_f[i] = _cosine6(v_self_norm, v_donor_norm)
                        out_frag_scr_f[i] = _scribe6(v_self_norm, v_donor_norm)
                        out_frag_kt_f[i]  = _kendall6(v_self_norm, v_donor_norm)
                    end
                end
                out_miss_f[i] = false
            end
        end
    end

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
        MBR_frag_pattern_cosine_true   = out_frag_cos_t,
        MBR_frag_pattern_cosine_false  = out_frag_cos_f,
        MBR_frag_pattern_scribe_true   = out_frag_scr_t,
        MBR_frag_pattern_scribe_false  = out_frag_scr_f,
        MBR_best_rt_diff_true          = out_rt_t,
        MBR_best_rt_diff_false         = out_rt_f,
        MBR_frag_rank_corr_true        = out_frag_kt_t,
        MBR_frag_rank_corr_false       = out_frag_kt_f,
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
    @user_info "  Wrote $n_written recovery sidecars"
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

        df = DataFrame(main)
        # decoy mirror of target (added in in-memory path for downstream consumers)
        df[!, :decoy] = .!Bool.(df.target)
        df[!, :trace_prob_prepass] = collect(Float32.(pass1.trace_prob_prepass))
        df[!, :trace_prob_infold]  = collect(Float32.(pass1.trace_prob_infold))
        # trace_prob = NON-MBR pass-1 score (same as in-memory path line 138).
        df[!, :trace_prob] = df[!, :trace_prob_prepass]
        for c in _MBR_SIDECAR_OUT_COLS
            c === :precursor_idx && continue
            c === :scan_idx      && continue
            df[!, c] = collect(Tables.getcolumn(mbr, c))
        end
        df[!, :mbr_recovered]          = collect(Bool.(rec.mbr_recovered))
        df[!, :MBR_transfer_candidate] = collect(Bool.(rec.MBR_transfer_candidate))
        df[!, :ftr_qval_true]          = collect(Float32.(rec.ftr_qval_true))
        df[!, :ftr_pep_true]           = collect(Float32.(rec.ftr_pep_true))

        writeArrow(path, df)
        if cleanup
            rm(pass1_path); rm(mbr_path); rm(rec_path)
        end
        n_merged += 1
    end
    @user_info "  Merged $n_merged files (Pass-1 + MBR + recovery sidecars)"
    return n_merged
end

