# Empirical observed-fragment shape features for cross-run scoring.
#
# MainSearch writes one best PSM per precursor/run. For each precursor, keep
# the top K per-run PSMs by MainSearch score and score each row against a
# leave-one-out consensus of those empirical spectra. The consensus uses
# deconvolved observed ("shadow") top-8 fragment trace intensities, weighted
# by per-run confidence (1 - PEP).

const EMPIRICAL_FRAGMENT_REF_K = 5
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

@inline function _empirical_fragment_consensus_reference(ref::EmpiricalFragmentTopKRef, row_id::UInt64)
    h1 = 0.0f0
    h2 = 0.0f0
    h3 = 0.0f0
    h4 = 0.0f0
    h5 = 0.0f0
    h6 = 0.0f0
    h7 = 0.0f0
    h8 = 0.0f0
    weight_sum = 0.0f0
    pep_weighted_sum = 0.0f0

    @inbounds for j in 1:EMPIRICAL_FRAGMENT_REF_K
        rid = ref.rows[j]
        (rid == UInt64(0) || rid == row_id) && continue
        weight = max(1.0f0 - ref.peps[j], 0.0f0)
        weight > 0.0f0 || continue
        sqrt_frag = ref.sqrt_frags[j]
        _empirical_fragment_sqrt_is_valid(sqrt_frag) || continue

        h1 += weight * sqrt_frag[1]
        h2 += weight * sqrt_frag[2]
        h3 += weight * sqrt_frag[3]
        h4 += weight * sqrt_frag[4]
        h5 += weight * sqrt_frag[5]
        h6 += weight * sqrt_frag[6]
        h7 += weight * sqrt_frag[7]
        h8 += weight * sqrt_frag[8]
        weight_sum += weight
        pep_weighted_sum += weight * ref.peps[j]
    end

    norm = h1*h1 + h2*h2 + h3*h3 + h4*h4 + h5*h5 + h6*h6 + h7*h7 + h8*h8
    if weight_sum <= 0.0f0 || norm <= 0.0f0
        return (
            found = false,
            sqrt_frag = EMPIRICAL_FRAGMENT_EMPTY_VEC,
            pep = EMPIRICAL_FRAGMENT_SENTINEL_REF_PEP,
        )
    end

    inv_norm = 1.0f0 / sqrt(norm)
    return (
        found = true,
        sqrt_frag = (
            h1 * inv_norm, h2 * inv_norm, h3 * inv_norm, h4 * inv_norm,
            h5 * inv_norm, h6 * inv_norm, h7 * inv_norm, h8 * inv_norm,
        ),
        pep = pep_weighted_sum / weight_sum,
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
        picked = _empirical_fragment_consensus_reference(ref, row_ids[i])
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
            picked = _empirical_fragment_consensus_reference(ref, row_start + UInt64(i - 1))
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

function add_empirical_fragment_features_to_fold_files!(file_paths::Vector{String})
    isempty(file_paths) && return nothing

    n_files = 0
    n_rows = 0
    n_refs = 0
    t = @elapsed begin
        refs, row_starts, n_rows_total = _empirical_fragment_build_refs(file_paths)
        n_refs = length(refs)
        n_rows = n_rows_total
        for (file_idx, fpath) in enumerate(file_paths)
            isfile(fpath) || continue
            add_empirical_fragment_features_to_fold_file!(fpath, refs, row_starts[file_idx])
            n_files += 1
        end
    end

    @debug_l1 "Empirical fragment cross-run features: added $(length(EMPIRICAL_FRAGMENT_FEATURES)) features " *
              "to $n_files fold files ($n_rows rows; refs=$n_refs; topK=$(EMPIRICAL_FRAGMENT_REF_K), LOO, PEP-weighted) " *
              "in $(round(t, digits = 2))s"
    return nothing
end
