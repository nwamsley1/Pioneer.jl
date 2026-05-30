# Empirical observed-fragment shape features for cross-run scoring.
#
# MainSearch writes one best PSM per precursor/run. For each precursor, keep
# the top two per-run PSMs by MainSearch score so each row can be scored
# against the best empirical reference that is not itself.

const EMPIRICAL_FRAGMENT_COLUMNS = (
    :frag1_int, :frag2_int, :frag3_int, :frag4_int,
    :frag5_int, :frag6_int, :frag7_int, :frag8_int,
)
const EMPIRICAL_FRAGMENT_SENTINEL_HELLINGER = 1.0f0
const EMPIRICAL_FRAGMENT_SENTINEL_REF_PEP = 1.0f0
const EMPIRICAL_FRAGMENT_EMPTY_VEC = ntuple(_ -> 0.0f0, 8)

mutable struct EmpiricalFragmentTop2Ref
    row1::UInt64
    score1::Float32
    pep1::Float32
    frag1::NTuple{8, Float32}
    row2::UInt64
    score2::Float32
    pep2::Float32
    frag2::NTuple{8, Float32}
end

function _empirical_fragment_empty_ref()
    return EmpiricalFragmentTop2Ref(
        UInt64(0),
        typemin(Float32),
        EMPIRICAL_FRAGMENT_SENTINEL_REF_PEP,
        EMPIRICAL_FRAGMENT_EMPTY_VEC,
        UInt64(0),
        typemin(Float32),
        EMPIRICAL_FRAGMENT_SENTINEL_REF_PEP,
        EMPIRICAL_FRAGMENT_EMPTY_VEC,
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

@inline function _empirical_fragment_hellinger_distance(
    obs::NTuple{8, Float32},
    ref::NTuple{8, Float32},
)
    obs_sum = 0.0f0
    ref_sum = 0.0f0
    @inbounds for i in 1:8
        obs_sum += obs[i]
        ref_sum += ref[i]
    end
    (obs_sum > 0.0f0 && ref_sum > 0.0f0) || return EMPIRICAL_FRAGMENT_SENTINEL_HELLINGER

    dist2 = 0.0f0
    @inbounds for i in 1:8
        po = obs[i] / obs_sum
        pr = ref[i] / ref_sum
        d = sqrt(po) - sqrt(pr)
        dist2 += d * d
    end
    return sqrt(max(0.0f0, 0.5f0 * dist2))
end

@inline function _empirical_fragment_update_top2!(
    refs::Dict{UInt32, EmpiricalFragmentTop2Ref},
    pid::UInt32,
    row_id::UInt64,
    score::Float32,
    pep::Float32,
    frag::NTuple{8, Float32},
)
    ref = get!(refs, pid) do
        _empirical_fragment_empty_ref()
    end

    if score > ref.score1
        ref.row2 = ref.row1
        ref.score2 = ref.score1
        ref.pep2 = ref.pep1
        ref.frag2 = ref.frag1
        ref.row1 = row_id
        ref.score1 = score
        ref.pep1 = pep
        ref.frag1 = frag
    elseif score > ref.score2
        ref.row2 = row_id
        ref.score2 = score
        ref.pep2 = pep
        ref.frag2 = frag
    end
    return nothing
end

function _empirical_fragment_top2_refs(
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

    refs = Dict{UInt32, EmpiricalFragmentTop2Ref}()
    sizehint!(refs, min(n, 1_000_000))
    @inbounds for i in 1:n
        pid = UInt32(precursor_idx[i])
        score = _empirical_fragment_safe_float(scores[i], typemin(Float32))
        pep = clamp(_empirical_fragment_safe_float(peps[i], EMPIRICAL_FRAGMENT_SENTINEL_REF_PEP), 0.0f0, 1.0f0)
        frag = _empirical_fragment_tuple(frag_cols, i)
        _empirical_fragment_update_top2!(refs, pid, row_ids[i], score, pep, frag)
    end
    return refs
end

@inline function _empirical_fragment_pick_reference(ref::EmpiricalFragmentTop2Ref, row_id::UInt64)
    if ref.row1 != UInt64(0) && ref.row1 != row_id
        return (found = true, frag = ref.frag1, pep = ref.pep1)
    elseif ref.row2 != UInt64(0) && ref.row2 != row_id
        return (found = true, frag = ref.frag2, pep = ref.pep2)
    else
        return (
            found = false,
            frag = EMPIRICAL_FRAGMENT_EMPTY_VEC,
            pep = EMPIRICAL_FRAGMENT_SENTINEL_REF_PEP,
        )
    end
end

function _empirical_fragment_reference_features(
    precursor_idx::AbstractVector,
    row_ids::AbstractVector{UInt64},
    frag_cols,
    refs::Dict{UInt32, EmpiricalFragmentTop2Ref},
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
        picked = _empirical_fragment_pick_reference(ref, row_ids[i])
        picked.found || continue
        obs = _empirical_fragment_tuple(frag_cols, i)
        hellinger[i] = _empirical_fragment_hellinger_distance(obs, picked.frag)
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
    refs = Dict{UInt32, EmpiricalFragmentTop2Ref}()
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
            frag = _empirical_fragment_tuple(frag_cols, i)
            _empirical_fragment_update_top2!(refs, pid, rid, score, pep, frag)
        end
    end

    return refs, row_starts, n_rows
end

function add_empirical_fragment_features_to_fold_file!(
    fpath::String,
    refs::Dict{UInt32, EmpiricalFragmentTop2Ref},
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
            picked = _empirical_fragment_pick_reference(ref, row_start + UInt64(i - 1))
            picked.found || continue
            obs = _empirical_fragment_tuple(frag_cols, i)
            hellinger[i] = _empirical_fragment_hellinger_distance(obs, picked.frag)
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
              "to $n_files fold files ($n_rows rows; refs=$n_refs) in $(round(t, digits = 2))s"
    return nothing
end
