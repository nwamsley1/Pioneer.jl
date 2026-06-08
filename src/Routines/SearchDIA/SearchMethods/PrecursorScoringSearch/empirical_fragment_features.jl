# Empirical observed-fragment shape features for cross-run scoring.
#
# MainSearch writes one best PSM per precursor/run. For each precursor, keep
# the top 2 per-run PSMs by MainSearch score and score each row against the
# best one that is not itself. References use deconvolved observed ("shadow")
# top-8 fragment trace intensities.

const EMPIRICAL_FRAGMENT_REF_K = 2
const EMPIRICAL_FRAGMENT_COLUMNS = (
    :shadow_frag1_int, :shadow_frag2_int, :shadow_frag3_int, :shadow_frag4_int,
    :shadow_frag5_int, :shadow_frag6_int, :shadow_frag7_int, :shadow_frag8_int,
)
const EMPIRICAL_FRAGMENT_SENTINEL_HELLINGER = 1.0f0
const EMPIRICAL_FRAGMENT_SENTINEL_REF_PEP = 1.0f0
const EMPIRICAL_FRAGMENT_EMPTY_VEC = ntuple(_ -> 0.0f0, 8)
const EMPIRICAL_FRAGMENT_EMPTY_ROWS = ntuple(_ -> UInt64(0), EMPIRICAL_FRAGMENT_REF_K)
const EMPIRICAL_FRAGMENT_EMPTY_SCORES = ntuple(_ -> typemin(Float32), EMPIRICAL_FRAGMENT_REF_K)
const EMPIRICAL_FRAGMENT_EMPTY_PEPS = ntuple(_ -> EMPIRICAL_FRAGMENT_SENTINEL_REF_PEP, EMPIRICAL_FRAGMENT_REF_K)
const EMPIRICAL_FRAGMENT_EMPTY_FRAGS = ntuple(_ -> EMPIRICAL_FRAGMENT_EMPTY_VEC, EMPIRICAL_FRAGMENT_REF_K)

struct EmpiricalFragmentHuberRefreshStats
    reference_rows::Int64
    replaced_rows::Int64
    fallback_rows::Int64
    unique_scans::Int64
    huber_rows::Int64
    runtime_sec::Float64
end

struct EmpiricalFragmentReferenceRow
    row_id::UInt64
    ms_file_idx::UInt32
    scan_idx::UInt32
    precursor_idx::UInt32
end

mutable struct EmpiricalFragmentTopKRef
    rows::NTuple{EMPIRICAL_FRAGMENT_REF_K, UInt64}
    scores::NTuple{EMPIRICAL_FRAGMENT_REF_K, Float32}
    peps::NTuple{EMPIRICAL_FRAGMENT_REF_K, Float32}
    sqrt_frags::NTuple{EMPIRICAL_FRAGMENT_REF_K, NTuple{8, Float32}}
end

function _empirical_fragment_empty_ref()
    return EmpiricalFragmentTopKRef(
        EMPIRICAL_FRAGMENT_EMPTY_ROWS,
        EMPIRICAL_FRAGMENT_EMPTY_SCORES,
        EMPIRICAL_FRAGMENT_EMPTY_PEPS,
        EMPIRICAL_FRAGMENT_EMPTY_FRAGS,
    )
end

function _empirical_fragment_empty_huber_stats()
    return EmpiricalFragmentHuberRefreshStats(0, 0, 0, 0, 0, 0.0)
end

@inline function _empirical_fragment_safe_float(x, default::Float32)
    ismissing(x) && return default
    y = Float32(x)
    return isfinite(y) ? y : default
end

@inline function _empirical_fragment_tuple(frag_cols, row::Integer)
    return ntuple(i -> max(_empirical_fragment_safe_float(frag_cols[i][row], 0.0f0), 0.0f0), 8)
end

@inline function _empirical_fragment_sqrt_tuple(frag::NTuple{8, Float32})
    frag_sum = 0.0f0
    @inbounds for i in 1:8
        frag_sum += frag[i]
    end
    frag_sum > 0.0f0 || return EMPIRICAL_FRAGMENT_EMPTY_VEC

    inv_sum = 1.0f0 / frag_sum
    return ntuple(i -> sqrt(frag[i] * inv_sum), 8)
end

@inline function _empirical_fragment_sqrt_is_valid(frag_sqrt::NTuple{8, Float32})
    sumsq = 0.0f0
    @inbounds for i in 1:8
        sumsq += frag_sqrt[i] * frag_sqrt[i]
    end
    return sumsq > 0.0f0
end

@inline function _empirical_fragment_hellinger_from_sqrt(
    obs_sqrt::NTuple{8, Float32},
    ref_sqrt::NTuple{8, Float32},
)
    dist2 = 0.0f0
    @inbounds for i in 1:8
        d = obs_sqrt[i] - ref_sqrt[i]
        dist2 += d * d
    end
    return sqrt(max(0.0f0, 0.5f0 * dist2))
end

function _empirical_fragment_apply_reference_overrides!(
    refs::Dict{UInt32, EmpiricalFragmentTopKRef},
    overrides::Dict{UInt64, NTuple{8, Float32}};
    unique_scans::Integer = 0,
    huber_rows::Integer = 0,
    runtime_sec::Real = 0.0,
)
    reference_rows = 0
    replaced_rows = 0

    for ref in values(refs)
        old_sqrt_frags = ref.sqrt_frags
        new_sqrt_frags = Vector{NTuple{8, Float32}}(undef, EMPIRICAL_FRAGMENT_REF_K)
        @inbounds for j in 1:EMPIRICAL_FRAGMENT_REF_K
            row_id = ref.rows[j]
            old_sqrt = old_sqrt_frags[j]
            if row_id == UInt64(0) || !_empirical_fragment_sqrt_is_valid(old_sqrt)
                new_sqrt_frags[j] = old_sqrt
                continue
            end

            reference_rows += 1
            override = get(overrides, row_id, nothing)
            if override !== nothing && _empirical_fragment_sqrt_is_valid(override)
                new_sqrt_frags[j] = override
                replaced_rows += 1
            else
                new_sqrt_frags[j] = old_sqrt
            end
        end
        ref.sqrt_frags = ntuple(j -> new_sqrt_frags[j], EMPIRICAL_FRAGMENT_REF_K)
    end

    return EmpiricalFragmentHuberRefreshStats(
        reference_rows,
        replaced_rows,
        reference_rows - replaced_rows,
        Int64(unique_scans),
        Int64(huber_rows),
        Float64(runtime_sec),
    )
end

@inline function _empirical_fragment_hellinger_distance(
    obs::NTuple{8, Float32},
    ref::NTuple{8, Float32},
)
    obs_sqrt = _empirical_fragment_sqrt_tuple(obs)
    ref_sqrt = _empirical_fragment_sqrt_tuple(ref)
    (_empirical_fragment_sqrt_is_valid(obs_sqrt) &&
     _empirical_fragment_sqrt_is_valid(ref_sqrt)) ||
        return EMPIRICAL_FRAGMENT_SENTINEL_HELLINGER
    return _empirical_fragment_hellinger_from_sqrt(obs_sqrt, ref_sqrt)
end

@inline function _empirical_fragment_update_topk!(
    refs::Dict{UInt32, EmpiricalFragmentTopKRef},
    pid::UInt32,
    row_id::UInt64,
    score::Float32,
    pep::Float32,
    sqrt_frag::NTuple{8, Float32},
)
    _empirical_fragment_sqrt_is_valid(sqrt_frag) || return nothing
    ref = get!(refs, pid) do
        _empirical_fragment_empty_ref()
    end

    insert_pos = 0
    @inbounds for j in 1:EMPIRICAL_FRAGMENT_REF_K
        if score > ref.scores[j]
            insert_pos = j
            break
        end
    end
    insert_pos == 0 && return nothing

    old_rows = ref.rows
    old_scores = ref.scores
    old_peps = ref.peps
    old_sqrt_frags = ref.sqrt_frags
    ref.rows = ntuple(
        j -> j < insert_pos ? old_rows[j] :
             j == insert_pos ? row_id :
             old_rows[j - 1],
        EMPIRICAL_FRAGMENT_REF_K,
    )
    ref.scores = ntuple(
        j -> j < insert_pos ? old_scores[j] :
             j == insert_pos ? score :
             old_scores[j - 1],
        EMPIRICAL_FRAGMENT_REF_K,
    )
    ref.peps = ntuple(
        j -> j < insert_pos ? old_peps[j] :
             j == insert_pos ? pep :
             old_peps[j - 1],
        EMPIRICAL_FRAGMENT_REF_K,
    )
    ref.sqrt_frags = ntuple(
        j -> j < insert_pos ? old_sqrt_frags[j] :
             j == insert_pos ? sqrt_frag :
             old_sqrt_frags[j - 1],
        EMPIRICAL_FRAGMENT_REF_K,
    )
    return nothing
end

function _empirical_fragment_topk_refs(
    precursor_idx::AbstractVector,
    row_ids::AbstractVector{UInt64},
    scores::AbstractVector,
    peps::AbstractVector,
    frag_cols,
)
    n = length(precursor_idx)
    length(row_ids) == n || throw(ArgumentError("row_ids length must match precursor_idx"))
    length(scores) == n || throw(ArgumentError("scores length must match precursor_idx"))
    length(peps) == n || throw(ArgumentError("peps length must match precursor_idx"))
    all(c -> length(c) == n, frag_cols) ||
        throw(ArgumentError("fragment columns must match precursor_idx length"))

    refs = Dict{UInt32, EmpiricalFragmentTopKRef}()
    sizehint!(refs, min(n, 1_000_000))
    @inbounds for i in 1:n
        pid = UInt32(precursor_idx[i])
        score = _empirical_fragment_safe_float(scores[i], typemin(Float32))
        pep = clamp(_empirical_fragment_safe_float(peps[i], EMPIRICAL_FRAGMENT_SENTINEL_REF_PEP), 0.0f0, 1.0f0)
        sqrt_frag = _empirical_fragment_sqrt_tuple(_empirical_fragment_tuple(frag_cols, i))
        _empirical_fragment_update_topk!(refs, pid, row_ids[i], score, pep, sqrt_frag)
    end
    return refs
end

@inline function _empirical_fragment_reference(ref::EmpiricalFragmentTopKRef, row_id::UInt64)
    @inbounds for j in 1:EMPIRICAL_FRAGMENT_REF_K
        rid = ref.rows[j]
        (rid == UInt64(0) || rid == row_id) && continue
        sqrt_frag = ref.sqrt_frags[j]
        _empirical_fragment_sqrt_is_valid(sqrt_frag) || continue
        return (
            found = true,
            sqrt_frag = sqrt_frag,
            pep = ref.peps[j],
        )
    end

    return (
        found = false,
        sqrt_frag = EMPIRICAL_FRAGMENT_EMPTY_VEC,
        pep = EMPIRICAL_FRAGMENT_SENTINEL_REF_PEP,
    )
end

function _empirical_fragment_reference_features(
    precursor_idx::AbstractVector,
    row_ids::AbstractVector{UInt64},
    frag_cols,
    refs::Dict{UInt32, EmpiricalFragmentTopKRef},
)
    n = length(precursor_idx)
    length(row_ids) == n || throw(ArgumentError("row_ids length must match precursor_idx"))
    all(c -> length(c) == n, frag_cols) ||
        throw(ArgumentError("fragment columns must match precursor_idx length"))

    hellinger = fill(EMPIRICAL_FRAGMENT_SENTINEL_HELLINGER, n)
    ref_pep = fill(EMPIRICAL_FRAGMENT_SENTINEL_REF_PEP, n)
    @inbounds for i in 1:n
        ref = get(refs, UInt32(precursor_idx[i]), nothing)
        ref === nothing && continue
        picked = _empirical_fragment_reference(ref, row_ids[i])
        picked.found || continue
        obs_sqrt = _empirical_fragment_sqrt_tuple(_empirical_fragment_tuple(frag_cols, i))
        _empirical_fragment_sqrt_is_valid(obs_sqrt) || continue
        hellinger[i] = _empirical_fragment_hellinger_from_sqrt(obs_sqrt, picked.sqrt_frag)
        ref_pep[i] = picked.pep
    end
    return hellinger, ref_pep
end

function _empirical_fragment_score_column(tbl)
    hasproperty(tbl, :lgbm_prob) && return tbl.lgbm_prob
    hasproperty(tbl, :lgbm_score) && return tbl.lgbm_score
    return nothing
end

function _empirical_fragment_has_required_columns(tbl)
    return hasproperty(tbl, :precursor_idx) &&
           hasproperty(tbl, :main_pep) &&
           _empirical_fragment_score_column(tbl) !== nothing &&
           all(c -> hasproperty(tbl, c), EMPIRICAL_FRAGMENT_COLUMNS)
end

function _empirical_fragment_build_refs(file_paths::Vector{String})
    refs = Dict{UInt32, EmpiricalFragmentTopKRef}()
    row_starts = Vector{UInt64}(undef, length(file_paths))
    row_id = UInt64(1)
    n_rows = 0

    for (file_idx, fpath) in enumerate(file_paths)
        tbl = Arrow.Table(fpath)
        n = length(tbl.precursor_idx)
        row_starts[file_idx] = row_id
        row_id += UInt64(n)
        n_rows += n
        n == 0 && continue
        _empirical_fragment_has_required_columns(tbl) || continue

        score_col = _empirical_fragment_score_column(tbl)
        frag_cols = ntuple(i -> getproperty(tbl, EMPIRICAL_FRAGMENT_COLUMNS[i]), 8)
        @inbounds for i in 1:n
            pid = UInt32(tbl.precursor_idx[i])
            rid = row_starts[file_idx] + UInt64(i - 1)
            score = _empirical_fragment_safe_float(score_col[i], typemin(Float32))
            pep = clamp(
                _empirical_fragment_safe_float(tbl.main_pep[i], EMPIRICAL_FRAGMENT_SENTINEL_REF_PEP),
                0.0f0,
                1.0f0,
            )
            sqrt_frag = _empirical_fragment_sqrt_tuple(_empirical_fragment_tuple(frag_cols, i))
            _empirical_fragment_update_topk!(refs, pid, rid, score, pep, sqrt_frag)
        end
    end

    return refs, row_starts, n_rows
end

function _empirical_fragment_reference_row_ids(refs::Dict{UInt32, EmpiricalFragmentTopKRef})
    row_ids = Set{UInt64}()
    for ref in values(refs)
        @inbounds for j in 1:EMPIRICAL_FRAGMENT_REF_K
            row_id = ref.rows[j]
            row_id == UInt64(0) && continue
            _empirical_fragment_sqrt_is_valid(ref.sqrt_frags[j]) || continue
            push!(row_ids, row_id)
        end
    end
    return row_ids
end

function _empirical_fragment_collect_reference_rows(
    file_paths::Vector{String},
    row_starts::Vector{UInt64},
    refs::Dict{UInt32, EmpiricalFragmentTopKRef},
)
    selected_row_ids = _empirical_fragment_reference_row_ids(refs)
    isempty(selected_row_ids) && return EmpiricalFragmentReferenceRow[]

    rows = EmpiricalFragmentReferenceRow[]
    sizehint!(rows, length(selected_row_ids))
    for (file_idx, fpath) in enumerate(file_paths)
        isfile(fpath) || continue
        tbl = Arrow.Table(fpath)
        n = length(tbl.precursor_idx)
        n == 0 && continue
        if !(hasproperty(tbl, :ms_file_idx) &&
             hasproperty(tbl, :scan_idx) &&
             hasproperty(tbl, :precursor_idx))
            continue
        end

        row_start = row_starts[file_idx]
        @inbounds for i in 1:n
            row_id = row_start + UInt64(i - 1)
            row_id in selected_row_ids || continue
            push!(
                rows,
                EmpiricalFragmentReferenceRow(
                    row_id,
                    UInt32(tbl.ms_file_idx[i]),
                    UInt32(tbl.scan_idx[i]),
                    UInt32(tbl.precursor_idx[i]),
                ),
            )
        end
    end
    return rows
end

@inline function _empirical_fragment_ref_key(ms_file_idx, scan_idx, precursor_idx)
    return (UInt32(ms_file_idx), UInt32(scan_idx), UInt32(precursor_idx))
end

function _empirical_fragment_huber_reference_overrides(
    search_context::SearchContext,
    main_search_params,
    reference_rows::Vector{EmpiricalFragmentReferenceRow},
)
    overrides = Dict{UInt64, NTuple{8, Float32}}()
    isempty(reference_rows) && return (
        overrides = overrides,
        unique_scans = 0,
        huber_rows = 0,
        runtime_sec = 0.0,
    )

    file_to_scans = Dict{Int64, Set{Int64}}()
    key_to_row = Dict{Tuple{UInt32, UInt32, UInt32}, UInt64}()
    for row in reference_rows
        scans = get!(file_to_scans, Int64(row.ms_file_idx)) do
            Set{Int64}()
        end
        push!(scans, Int64(row.scan_idx))
        key = _empirical_fragment_ref_key(row.ms_file_idx, row.scan_idx, row.precursor_idx)
        haskey(key_to_row, key) || (key_to_row[key] = row.row_id)
    end

    ms_ref = getMSData(search_context)
    huber_solver = default_chromatogram_integration_huber_solver()
    huber_rows = 0
    unique_scans = 0
    runtime_sec = @elapsed begin
        for ms_file_idx in sort!(collect(keys(file_to_scans)))
            scans = sort!(collect(file_to_scans[ms_file_idx]))
            isempty(scans) && continue
            unique_scans += length(scans)
            spectra = getMSData(ms_ref, ms_file_idx)
            huber_psms = library_search(
                spectra,
                search_context,
                main_search_params,
                ms_file_idx;
                scan_indices = scans,
                deconvolution_solver_override = huber_solver,
            )
            huber_rows += nrow(huber_psms)
            nrow(huber_psms) == 0 && continue
            all(c -> hasproperty(huber_psms, c), EMPIRICAL_FRAGMENT_COLUMNS) || continue
            frag_cols = ntuple(i -> huber_psms[!, EMPIRICAL_FRAGMENT_COLUMNS[i]], 8)

            @inbounds for i in 1:nrow(huber_psms)
                key = _empirical_fragment_ref_key(
                    huber_psms.ms_file_idx[i],
                    huber_psms.scan_idx[i],
                    huber_psms.precursor_idx[i],
                )
                row_id = get(key_to_row, key, UInt64(0))
                row_id == UInt64(0) && continue
                haskey(overrides, row_id) && continue
                sqrt_frag = _empirical_fragment_sqrt_tuple(_empirical_fragment_tuple(frag_cols, i))
                _empirical_fragment_sqrt_is_valid(sqrt_frag) || continue
                overrides[row_id] = sqrt_frag
            end
        end
    end

    return (
        overrides = overrides,
        unique_scans = unique_scans,
        huber_rows = huber_rows,
        runtime_sec = runtime_sec,
    )
end

function _empirical_fragment_refresh_refs_with_huber!(
    refs::Dict{UInt32, EmpiricalFragmentTopKRef},
    file_paths::Vector{String},
    row_starts::Vector{UInt64},
    search_context::SearchContext,
    main_search_params,
)
    reference_rows = _empirical_fragment_collect_reference_rows(file_paths, row_starts, refs)
    refresh = _empirical_fragment_huber_reference_overrides(
        search_context,
        main_search_params,
        reference_rows,
    )
    return _empirical_fragment_apply_reference_overrides!(
        refs,
        refresh.overrides;
        unique_scans = refresh.unique_scans,
        huber_rows = refresh.huber_rows,
        runtime_sec = refresh.runtime_sec,
    )
end

function add_empirical_fragment_features_to_fold_file!(
    fpath::String,
    refs::Dict{UInt32, EmpiricalFragmentTopKRef},
    row_start::UInt64,
)
    tbl = DataFrame(Tables.columntable(Arrow.Table(fpath)); copycols = true)
    n = nrow(tbl)
    hellinger = fill(EMPIRICAL_FRAGMENT_SENTINEL_HELLINGER, n)
    ref_pep = fill(EMPIRICAL_FRAGMENT_SENTINEL_REF_PEP, n)

    if n > 0 && _empirical_fragment_has_required_columns(tbl)
        frag_cols = ntuple(i -> tbl[!, EMPIRICAL_FRAGMENT_COLUMNS[i]], 8)
        @inbounds for i in 1:n
            ref = get(refs, UInt32(tbl.precursor_idx[i]), nothing)
            ref === nothing && continue
            picked = _empirical_fragment_reference(ref, row_start + UInt64(i - 1))
            picked.found || continue
            obs_sqrt = _empirical_fragment_sqrt_tuple(_empirical_fragment_tuple(frag_cols, i))
            _empirical_fragment_sqrt_is_valid(obs_sqrt) || continue
            hellinger[i] = _empirical_fragment_hellinger_from_sqrt(obs_sqrt, picked.sqrt_frag)
            ref_pep[i] = picked.pep
        end
    end

    tbl[!, :empirical_frag_best_hellinger] = hellinger
    tbl[!, :empirical_frag_ref_pep] = ref_pep
    writeArrow(fpath, tbl)
    return n
end

function add_empirical_fragment_features_to_fold_files!(
    file_paths::Vector{String};
    search_context = nothing,
    main_search_params = nothing,
)
    isempty(file_paths) && return nothing

    n_files = 0
    n_rows = 0
    n_refs = 0
    huber_stats = _empirical_fragment_empty_huber_stats()
    t = @elapsed begin
        refs, row_starts, n_rows_total = _empirical_fragment_build_refs(file_paths)
        n_refs = length(refs)
        n_rows = n_rows_total
        if search_context !== nothing && main_search_params !== nothing
            huber_stats = _empirical_fragment_refresh_refs_with_huber!(
                refs,
                file_paths,
                row_starts,
                search_context,
                main_search_params,
            )
        end
        for (file_idx, fpath) in enumerate(file_paths)
            isfile(fpath) || continue
            add_empirical_fragment_features_to_fold_file!(fpath, refs, row_starts[file_idx])
            n_files += 1
        end
    end

    huber_msg = if search_context !== nothing && main_search_params !== nothing
        "; huber_refs=$(huber_stats.replaced_rows)/$(huber_stats.reference_rows) " *
        "fallback=$(huber_stats.fallback_rows) scans=$(huber_stats.unique_scans) " *
        "huber_rows=$(huber_stats.huber_rows) huber_time=$(round(huber_stats.runtime_sec, digits = 2))s"
    else
        ""
    end
    @debug_l1 "Empirical fragment cross-run features: added $(length(EMPIRICAL_FRAGMENT_FEATURES)) features " *
              "to $n_files fold files ($n_rows rows; refs=$n_refs; topK=$(EMPIRICAL_FRAGMENT_REF_K), LOO$huber_msg) " *
              "in $(round(t, digits = 2))s"
    return nothing
end

function add_empirical_fragment_features_to_extra_fold_files!(
    reference_paths::Vector{String},
    extra_paths::Vector{String},
)
    isempty(reference_paths) && return nothing
    isempty(extra_paths) && return nothing

    n_files = 0
    n_rows = 0
    n_refs = 0
    t = @elapsed begin
        refs, _, reference_rows = _empirical_fragment_build_refs(reference_paths)
        n_refs = length(refs)
        row_start = UInt64(reference_rows) + UInt64(1)
        for fpath in extra_paths
            isfile(fpath) || continue
            n = add_empirical_fragment_features_to_fold_file!(fpath, refs, row_start)
            row_start += UInt64(n)
            n_rows += n
            n_files += 1
        end
    end

    @debug_l1 "Empirical fragment cross-run features: added $(length(EMPIRICAL_FRAGMENT_FEATURES)) features " *
              "to $n_files rescue fold files ($n_rows rows; refs=$n_refs; topK=$(EMPIRICAL_FRAGMENT_REF_K), LOO) " *
              "in $(round(t, digits = 2))s"
    return nothing
end
