# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

#==========================================================
Experimental protein-group max feature (PG_MAX) for Pass-1 scoring.

After the first semi-supervised iteration, every precursor in a run gets the best
iteration-1 out-of-fold score among the OTHER peptides that share its library
accession string (the FASTA digest's "A;B;C", not an inferred protein group).
Targets, decoys and entrapment precursors are pooled, so a decoy on a present
protein gets the same protein evidence a false target there would. Iterations
2..N train with the column; the final prediction uses the kept iteration's columns.

CV folds are assigned by the same accession string, so a group never spans folds:
a fold's PG_MAX values come only from the model trained on the other fold.

Selected with ENV["PIONEER_PG_MAX"] (benchmark arms; default off):
  off  no feature (identical to no PG_MAX code)
  v1   PG_MAX from iteration-1 scores, own sequence left out
  v0   as v1 but the precursor's own sequence is included
  nc   as v1 on accession strings permuted across groups (negative control)
  v2   v1 plus the group size (distinct sequences with that string in the library)
  d1   as v1 but entrapment precursors form their own groups (diagnostic only)
  vr   as v1 but recomputed after every iteration from that iteration's scores
==========================================================#

const PG_MAX_MODES = (:off, :v1, :v0, :nc, :v2, :d1, :vr)
const PG_MAX_FEATURE = :pg_max
const PG_SIZE_FEATURE = :pg_size
# VR stops once fewer than this share of pool rows move PG_MAX by more than PG_MAX_VR_DELTA.
const PG_MAX_VR_DELTA = 0.05f0
const PG_MAX_VR_MIN_CHANGED = 0.01

function pg_max_mode()
    mode = Symbol(lowercase(strip(get(ENV, "PIONEER_PG_MAX", "off"))))
    mode in PG_MAX_MODES || throw(ArgumentError(
        "PIONEER_PG_MAX must be one of $(join(PG_MAX_MODES, ", ")); got $mode"))
    return mode
end

"""
Per-precursor grouping for PG_MAX, indexed by precursor_idx.

`group` is the accession-string id (or its arm-specific variant), `peptide` the
unmodified-sequence id used to leave out a precursor's own peptide (all charges and
modified forms), and `group_size` the distinct sequences per group (v2 only).
"""
struct ProteinGroupKeys
    mode::Symbol
    group::Vector{UInt32}
    peptide::Vector{UInt32}
    group_size::Vector{Float32}
end

pg_exclude_own(keys::ProteinGroupKeys) = keys.mode != :v0
pg_refresh(keys::ProteinGroupKeys) = keys.mode == :vr
pg_feature_names(keys::ProteinGroupKeys) =
    keys.mode == :v2 ? [PG_MAX_FEATURE, PG_SIZE_FEATURE] : [PG_MAX_FEATURE]

function build_protein_group_keys(precursors::LibraryPrecursors, mode::Symbol)
    ids = getPrecursorTextIds(precursors)
    group = copy(ids.accession_set_id)
    peptide = ids.sequence_id
    if mode == :d1
        entrap = getEntrapmentGroupId(precursors)
        n_sets = UInt32(maximum(group; init = UInt32(0)))
        @inbounds for i in eachindex(group)
            e = entrap[i]
            e == 0 || (group[i] += UInt32(e) * n_sets)
        end
    elseif mode == :nc
        group = _permute_groups_within_fold(group, getBasePepId(precursors), precursors.pid_to_cv_fold)
    end
    group_size = mode == :v2 ? _distinct_sequences_per_group(group, peptide) : Float32[]
    @user_info "PG_MAX mode=$mode: $(length(unique(group))) groups over $(length(group)) precursors"
    return ProteinGroupKeys(mode, group, peptide, group_size)
end

# Negative control: shuffle which group each base peptide belongs to, within its CV fold, so group
# sizes are kept, a peptide's target/decoy/entrapment versions stay together, and groups stay in one fold.
function _permute_groups_within_fold(group::Vector{UInt32}, base_pep_id::AbstractVector,
                                     pid_to_cv_fold::Vector{UInt8})
    unit_of = Dict{UInt32, Int}()
    unit_group = UInt32[]
    unit_fold = UInt8[]
    pid_unit = Vector{Int}(undef, length(group))
    @inbounds for i in eachindex(group)
        u = get!(unit_of, UInt32(base_pep_id[i])) do
            push!(unit_group, group[i])
            push!(unit_fold, pid_to_cv_fold[i])
            length(unit_group)
        end
        pid_unit[i] = u
    end
    rng = MersenneTwister(1776)
    for f in unique(unit_fold)
        idx = findall(==(f), unit_fold)
        unit_group[idx] = shuffle(rng, unit_group[idx])
    end
    return UInt32[unit_group[pid_unit[i]] for i in eachindex(group)]
end

function _distinct_sequences_per_group(group::Vector{UInt32}, peptide::Vector{UInt32})
    pairs = unique!(sort!(UInt64[(UInt64(group[i]) << 32) | peptide[i] for i in eachindex(group)]))
    counts = Dict{UInt32, Int}()
    for key in pairs
        g = UInt32(key >> 32)
        counts[g] = get(counts, g, 0) + 1
    end
    return Float32[counts[group[i]] for i in eachindex(group)]
end

"""
    protein_group_columns(keys, precursor_idx, is_target, scores) -> Matrix{Float32}

One file's PG_MAX column (and the group-size column for v2), one row per input row.
A precursor's score is its best row; a peptide's is its best precursor. PG_MAX is the
best peptide score in the precursor's group, leaving out its own peptide (except v0);
with no such peptide it is the file's median decoy precursor score.
"""
function protein_group_columns(keys::ProteinGroupKeys, precursor_idx::AbstractVector,
                               is_target::AbstractVector, scores::AbstractVector{Float32})
    n = length(precursor_idx)
    ncol = length(pg_feature_names(keys))
    out = Matrix{Float32}(undef, n, ncol)
    n == 0 && return out

    best_prec = Dict{UInt32, Float32}()
    decoy_prec = Dict{UInt32, Float32}()
    @inbounds for i in 1:n
        pid = UInt32(precursor_idx[i])
        s = scores[i]
        best_prec[pid] = max(get(best_prec, pid, -Inf32), s)
        is_target[i] || (decoy_prec[pid] = best_prec[pid])
    end
    floor = isempty(decoy_prec) ? minimum(scores) : Float32(median(collect(values(decoy_prec))))

    # Best score per (group, peptide), then each group's top two over distinct peptides.
    best_pep = Dict{UInt64, Float32}()
    for (pid, s) in best_prec
        key = (UInt64(keys.group[pid]) << 32) | keys.peptide[pid]
        best_pep[key] = max(get(best_pep, key, -Inf32), s)
    end
    top = Dict{UInt32, Tuple{Float32, UInt32, Float32}}()   # (best, its peptide, second best)
    for (key, s) in best_pep
        g = UInt32(key >> 32)
        pep = UInt32(key & 0xffffffff)
        v1, k1, v2 = get(top, g, (-Inf32, UInt32(0), -Inf32))
        top[g] = s > v1 ? (s, pep, v1) : (v1, k1, max(v2, s))
    end

    exclude_own = pg_exclude_own(keys)
    @inbounds for i in 1:n
        pid = UInt32(precursor_idx[i])
        v1, k1, v2 = top[keys.group[pid]]
        v = exclude_own && k1 == keys.peptide[pid] ? v2 : v1
        out[i, 1] = isfinite(v) ? v : floor
        ncol == 2 && (out[i, 2] = keys.group_size[pid])
    end
    return out
end

"""
    _pg_max_columns(file_paths, features, cls_trained_on, extra_by_file, keys, plan0, plan1, k0, k1)

Score every file's rows out-of-fold with `cls_trained_on` (whose models were trained on
`features` plus `extra_by_file`, if given), compute each file's PG columns, and gather them
for the training pool using the sampler's (source row, pool row) plans.
"""
function _pg_max_columns(file_paths::Vector{String}, features::Vector{Symbol},
                         cls_trained_on::NamedTuple,
                         extra_by_file::Union{Nothing, Vector{Matrix{Float32}}},
                         keys::ProteinGroupKeys, plan0, plan1, k0::Int, k1::Int)
    ncol = length(pg_feature_names(keys))
    by_file = Vector{Matrix{Float32}}(undef, length(file_paths))
    pool0 = Matrix{Float32}(undef, k0, ncol)
    pool1 = Matrix{Float32}(undef, k1, ncol)
    for (file_idx, path) in enumerate(file_paths)
        cols = with_arrow_table(path) do tbl
            extra = extra_by_file === nothing ? nothing : extra_by_file[file_idx]
            scores = _pass1_file_oof_scores(tbl, features, cls_trained_on, extra)
            protein_group_columns(keys, tbl.precursor_idx, tbl.target, scores)
        end
        by_file[file_idx] = cols
        for (src_i, dst_i) in plan0[file_idx]
            @views pool0[dst_i, :] .= cols[src_i, :]
        end
        for (src_i, dst_i) in plan1[file_idx]
            @views pool1[dst_i, :] .= cols[src_i, :]
        end
    end
    return (pool0 = pool0, pool1 = pool1, by_file = by_file)
end
