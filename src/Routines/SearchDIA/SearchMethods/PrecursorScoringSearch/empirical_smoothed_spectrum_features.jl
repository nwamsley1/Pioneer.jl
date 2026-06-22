# Cross-run empirical spectrum features for selected precursor traces.
#
# MainSearch writes one best PSM per precursor/run into fold files. For each
# precursor, keep the top 2 per-run smoothed fragment spectra by MainSearch
# probability and score each row against the best one that is not itself.

const EMPIRICAL_SMOOTHED_SPECTRUM_REF_K = 2
const EMPIRICAL_SMOOTHED_SPECTRUM_SENTINEL_HELLINGER = 1.0f0
const EMPIRICAL_SMOOTHED_SPECTRUM_EMPTY_VEC = ntuple(_ -> 0.0f0, 8)
const EMPIRICAL_SMOOTHED_SPECTRUM_EMPTY_ROWS =
    ntuple(_ -> UInt64(0), EMPIRICAL_SMOOTHED_SPECTRUM_REF_K)
const EMPIRICAL_SMOOTHED_SPECTRUM_EMPTY_SCORES =
    ntuple(_ -> typemin(Float32), EMPIRICAL_SMOOTHED_SPECTRUM_REF_K)
const EMPIRICAL_SMOOTHED_SPECTRUM_EMPTY_FRAGS =
    ntuple(_ -> EMPIRICAL_SMOOTHED_SPECTRUM_EMPTY_VEC, EMPIRICAL_SMOOTHED_SPECTRUM_REF_K)

mutable struct EmpiricalSmoothedSpectrumTopKRef
    rows::NTuple{EMPIRICAL_SMOOTHED_SPECTRUM_REF_K, UInt64}
    scores::NTuple{EMPIRICAL_SMOOTHED_SPECTRUM_REF_K, Float32}
    sqrt_frags::NTuple{EMPIRICAL_SMOOTHED_SPECTRUM_REF_K, NTuple{8, Float32}}
end

function _empirical_smoothed_spectrum_empty_ref()
    return EmpiricalSmoothedSpectrumTopKRef(
        EMPIRICAL_SMOOTHED_SPECTRUM_EMPTY_ROWS,
        EMPIRICAL_SMOOTHED_SPECTRUM_EMPTY_SCORES,
        EMPIRICAL_SMOOTHED_SPECTRUM_EMPTY_FRAGS,
    )
end

@inline function _empirical_smoothed_spectrum_tuple(frag_cols, row::Integer)
    return ntuple(i -> begin
        value = Float32(frag_cols[i][row])
        isfinite(value) ? max(value, 0.0f0) : 0.0f0
    end, 8)
end

@inline function _empirical_smoothed_spectrum_sqrt_tuple(frag::NTuple{8, Float32})
    frag_sum = 0.0f0
    @inbounds for i in 1:8
        frag_sum += frag[i]
    end
    frag_sum > 0.0f0 || return EMPIRICAL_SMOOTHED_SPECTRUM_EMPTY_VEC

    inv_sum = 1.0f0 / frag_sum
    return ntuple(i -> sqrt(frag[i] * inv_sum), 8)
end

@inline function _empirical_smoothed_spectrum_sqrt_is_valid(frag_sqrt::NTuple{8, Float32})
    @inbounds for i in 1:8
        frag_sqrt[i] > 0.0f0 && return true
    end
    return false
end

@inline function _empirical_smoothed_spectrum_hellinger_from_sqrt(
    obs_sqrt::NTuple{8, Float32},
    ref_sqrt::NTuple{8, Float32},
)
    dist2 = 0.0f0
    @inbounds for i in 1:8
        delta = obs_sqrt[i] - ref_sqrt[i]
        dist2 += delta * delta
    end
    return sqrt(max(0.0f0, 0.5f0 * dist2))
end

@inline function _empirical_smoothed_spectrum_update_topk!(
    refs::Dict{UInt32, EmpiricalSmoothedSpectrumTopKRef},
    pid::UInt32,
    row_id::UInt64,
    score::Float32,
    sqrt_frag::NTuple{8, Float32},
)
    _empirical_smoothed_spectrum_sqrt_is_valid(sqrt_frag) || return nothing
    ref = get!(refs, pid) do
        _empirical_smoothed_spectrum_empty_ref()
    end

    insert_pos = 0
    @inbounds for j in 1:EMPIRICAL_SMOOTHED_SPECTRUM_REF_K
        if score > ref.scores[j]
            insert_pos = j
            break
        end
    end
    insert_pos == 0 && return nothing

    old_rows = ref.rows
    old_scores = ref.scores
    old_sqrt_frags = ref.sqrt_frags
    ref.rows = ntuple(
        j -> j < insert_pos ? old_rows[j] :
             j == insert_pos ? row_id :
             old_rows[j - 1],
        EMPIRICAL_SMOOTHED_SPECTRUM_REF_K,
    )
    ref.scores = ntuple(
        j -> j < insert_pos ? old_scores[j] :
             j == insert_pos ? score :
             old_scores[j - 1],
        EMPIRICAL_SMOOTHED_SPECTRUM_REF_K,
    )
    ref.sqrt_frags = ntuple(
        j -> j < insert_pos ? old_sqrt_frags[j] :
             j == insert_pos ? sqrt_frag :
             old_sqrt_frags[j - 1],
        EMPIRICAL_SMOOTHED_SPECTRUM_REF_K,
    )
    return nothing
end

function _empirical_smoothed_spectrum_topk_refs(
    precursor_idx::AbstractVector,
    row_ids::AbstractVector{UInt64},
    scores::AbstractVector,
    frag_cols,
)
    refs = Dict{UInt32, EmpiricalSmoothedSpectrumTopKRef}()
    sizehint!(refs, min(length(precursor_idx), 1_000_000))
    @inbounds for i in eachindex(precursor_idx)
        sqrt_frag = _empirical_smoothed_spectrum_sqrt_tuple(
            _empirical_smoothed_spectrum_tuple(frag_cols, i),
        )
        _empirical_smoothed_spectrum_update_topk!(
            refs,
            UInt32(precursor_idx[i]),
            row_ids[i],
            Float32(scores[i]),
            sqrt_frag,
        )
    end
    return refs
end

@inline function _empirical_smoothed_spectrum_reference(
    ref::EmpiricalSmoothedSpectrumTopKRef,
    row_id::UInt64,
)
    @inbounds for j in 1:EMPIRICAL_SMOOTHED_SPECTRUM_REF_K
        rid = ref.rows[j]
        (rid == UInt64(0) || rid == row_id) && continue
        sqrt_frag = ref.sqrt_frags[j]
        _empirical_smoothed_spectrum_sqrt_is_valid(sqrt_frag) || continue
        return (found = true, sqrt_frag = sqrt_frag)
    end

    return (found = false, sqrt_frag = EMPIRICAL_SMOOTHED_SPECTRUM_EMPTY_VEC)
end

function _empirical_smoothed_spectrum_reference_features(
    precursor_idx::AbstractVector,
    row_ids::AbstractVector{UInt64},
    frag_cols,
    refs::Dict{UInt32, EmpiricalSmoothedSpectrumTopKRef},
)
    hellinger = fill(EMPIRICAL_SMOOTHED_SPECTRUM_SENTINEL_HELLINGER, length(precursor_idx))
    @inbounds for i in eachindex(precursor_idx)
        ref = get(refs, UInt32(precursor_idx[i]), nothing)
        ref === nothing && continue
        picked = _empirical_smoothed_spectrum_reference(ref, row_ids[i])
        picked.found || continue
        obs_sqrt = _empirical_smoothed_spectrum_sqrt_tuple(
            _empirical_smoothed_spectrum_tuple(frag_cols, i),
        )
        _empirical_smoothed_spectrum_sqrt_is_valid(obs_sqrt) || continue
        hellinger[i] = _empirical_smoothed_spectrum_hellinger_from_sqrt(
            obs_sqrt,
            picked.sqrt_frag,
        )
    end
    return hellinger
end

function _empirical_smoothed_spectrum_build_refs(file_paths::Vector{String})
    refs = Dict{UInt32, EmpiricalSmoothedSpectrumTopKRef}()
    row_starts = Vector{UInt64}(undef, length(file_paths))
    row_id = UInt64(1)
    n_rows = 0

    for (file_idx, fpath) in enumerate(file_paths)
        tbl = Arrow.Table(fpath)
        n = length(tbl.precursor_idx)
        row_starts[file_idx] = row_id
        row_id += UInt64(n)
        n_rows += n

        frag_cols = ntuple(i -> getproperty(tbl, SMOOTHED_FRAGMENT_INTENSITY_COLUMNS[i]), 8)
        @inbounds for i in 1:n
            sqrt_frag = _empirical_smoothed_spectrum_sqrt_tuple(
                _empirical_smoothed_spectrum_tuple(frag_cols, i),
            )
            _empirical_smoothed_spectrum_update_topk!(
                refs,
                UInt32(tbl.precursor_idx[i]),
                row_starts[file_idx] + UInt64(i - 1),
                Float32(tbl.lgbm_prob[i]),
                sqrt_frag,
            )
        end
    end

    return refs, row_starts, n_rows
end

function add_empirical_smoothed_spectrum_features_to_fold_file!(
    fpath::String,
    refs::Dict{UInt32, EmpiricalSmoothedSpectrumTopKRef},
    row_start::UInt64,
)
    tbl = DataFrame(Tables.columntable(Arrow.Table(fpath)); copycols = true)
    n = nrow(tbl)
    hellinger = fill(EMPIRICAL_SMOOTHED_SPECTRUM_SENTINEL_HELLINGER, n)
    frag_cols = ntuple(i -> tbl[!, SMOOTHED_FRAGMENT_INTENSITY_COLUMNS[i]], 8)

    @inbounds for i in 1:n
        ref = get(refs, UInt32(tbl.precursor_idx[i]), nothing)
        ref === nothing && continue
        picked = _empirical_smoothed_spectrum_reference(ref, row_start + UInt64(i - 1))
        picked.found || continue
        obs_sqrt = _empirical_smoothed_spectrum_sqrt_tuple(
            _empirical_smoothed_spectrum_tuple(frag_cols, i),
        )
        _empirical_smoothed_spectrum_sqrt_is_valid(obs_sqrt) || continue
        hellinger[i] = _empirical_smoothed_spectrum_hellinger_from_sqrt(
            obs_sqrt,
            picked.sqrt_frag,
        )
    end

    tbl[!, :empirical_smoothed_frag_hellinger] = hellinger
    writeArrow(fpath, tbl)
    return n
end

function add_empirical_smoothed_spectrum_features_to_fold_files!(file_paths::Vector{String})
    isempty(file_paths) && return nothing

    n_files = 0
    n_rows = 0
    n_refs = 0
    t = @elapsed begin
        refs, row_starts, n_rows_total = _empirical_smoothed_spectrum_build_refs(file_paths)
        n_refs = length(refs)
        n_rows = n_rows_total
        for (file_idx, fpath) in enumerate(file_paths)
            add_empirical_smoothed_spectrum_features_to_fold_file!(
                fpath,
                refs,
                row_starts[file_idx],
            )
            n_files += 1
        end
    end

    @debug_l1 "Empirical smoothed spectrum features: added " *
              "$(length(EMPIRICAL_SMOOTHED_SPECTRUM_FEATURES)) features " *
              "to $n_files fold files ($n_rows rows; refs=$n_refs; " *
              "topK=$(EMPIRICAL_SMOOTHED_SPECTRUM_REF_K), LOO) " *
              "in $(round(t, digits = 2))s"
    return nothing
end
