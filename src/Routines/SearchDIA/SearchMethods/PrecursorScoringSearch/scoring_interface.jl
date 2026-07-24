# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
Scoring interface for precursor probability aggregation and dictionary construction.
"""

#==========================================================
Per-File Precursor Aggregation
==========================================================#

"""
    _aggregate_trace_to_precursor_probs!(df)

Per-file Bayesian aggregation of trace-level → precursor-level probabilities.
Groups by (precursor_idx, ms_file_idx). Since ms_file_idx is constant within
a single file, this is effectively grouping by precursor_idx alone.
"""
function _aggregate_trace_to_precursor_probs!(df::DataFrame)
    prob_agg = p -> begin
        trace_prob = 1.0f0 - eps(Float32) - exp(sum(log1p.(-p)))
        clamp(trace_prob, eps(Float32), 1.0f0 - eps(Float32))
    end
    transform!(groupby(df, [:precursor_idx, :ms_file_idx]),
               :trace_prob => prob_agg => :prec_prob)
    return df
end

# `_aggregate_trace_to_precursor_probs!` is now called directly in Step 1b
# (PrecursorScoringSearch.jl) on the merged DataFrame before it is written, so
# `:prec_prob` ships as a main column and no separate per-file aggregation pass is needed.
# The former `aggregate_per_file!` (which re-read every merged file to add one sidecar
# column) has been removed.

#==========================================================
Dictionary + Sidecar Helper Functions for OOM Scoring Pipeline
==========================================================#

"""
    build_global_qval_dict_from_scores(score_dict, target_dict, fdr_scale) → Dict{UInt32, Float32}

Compute global q-values from a score dictionary without any file I/O.
"""
function build_global_qval_dict_from_scores(
    score_dict::Dict{UInt32, Float32},
    target_dict::Dict{UInt32, Bool},
    fdr_scale::Float32
)
    n = length(score_dict)
    pids = collect(keys(score_dict))
    scores = Float32[score_dict[pid] for pid in pids]
    targets = Bool[target_dict[pid] for pid in pids]

    # Sort descending by score
    perm = sortperm(scores; rev=true)
    permute!(pids, perm)
    permute!(scores, perm)
    permute!(targets, perm)

    # Compute q-values
    qvals = Vector{Float32}(undef, n)
    get_qvalues!(scores, targets, qvals; fdr_scale_factor=fdr_scale)

    # Build dictionary
    qval_dict = Dict{UInt32, Float32}()
    sizehint!(qval_dict, n)
    for i in 1:n
        qval_dict[pids[i]] = qvals[i]
    end
    return qval_dict
end

"""
    build_global_pep_dict_from_scores(score_dict, target_dict, fdr_scale) → Dict{UInt32, Float32}

Compute global posterior error probabilities (local FDR) from a score dictionary
without any file I/O. Parallel to `build_global_qval_dict_from_scores` but uses
`get_PEP!` (PAVA-fit) instead of cumulative q-values.
"""
function build_global_pep_dict_from_scores(
    score_dict::Dict{UInt32, Float32},
    target_dict::Dict{UInt32, Bool},
    fdr_scale::Float32
)
    n = length(score_dict)
    pids = collect(keys(score_dict))
    scores = Float32[score_dict[pid] for pid in pids]
    targets = Bool[target_dict[pid] for pid in pids]

    peps = Vector{Float32}(undef, n)
    get_PEP!(scores, targets, peps; doSort=true, fdr_scale_factor=fdr_scale)

    pep_dict = Dict{UInt32, Float32}()
    sizehint!(pep_dict, n)
    for i in 1:n
        pep_dict[pids[i]] = peps[i]
    end
    return pep_dict
end


"""
    write_score_sidecars(refs, columns; temp_prefix) → Vector{PSMFileReference}

Extract only the named columns from each file into a temporary Arrow sidecar file.
"""
function write_score_sidecars(
    refs::Vector{<:FileReference},
    columns::Vector{Symbol};
    temp_prefix::String = "sidecar"
)
    sidecar_refs = PSMFileReference[]
    for ref in refs
        # Use materialize_columns so columns are pulled from main OR any
        # registered sidecar (e.g. :prec_prob now lives in a sidecar after
        # aggregate_per_file!).
        df = ref isa PSMFileReference ? materialize_columns(ref, columns) :
             DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))[!, columns]
        nrow(df) == 0 && continue
        temp_path = tempname() * "_$(temp_prefix).arrow"
        writeArrow(temp_path, df)
        push!(sidecar_refs, PSMFileReference(temp_path))
    end
    return sidecar_refs
end

"""
    build_qvalue_spline_from_refs(refs, score_col, merged_path; ...) → Union{Nothing, NamedTuple}

Encapsulates the full sidecar lifecycle: write → sort → merge → cleanup → spline computation.
"""
function build_qvalue_spline_from_refs(
    refs::Vector{<:FileReference},
    score_col::Symbol,
    merged_path::String;
    batch_size::Int = 10_000_000,
    compute_pep::Bool = false,
    min_pep_points_per_bin::Int = 100,
    fdr_scale_factor::Float32 = 1.0f0,
    temp_prefix::String = "sidecar"
)
    sidecar_refs = @score_phase "    A3.1 write-score-sidecars" write_score_sidecars(
        refs, [score_col, :target]; temp_prefix=temp_prefix)
    isempty(sidecar_refs) && return nothing

    try
        @score_phase "    A3.2 sort-sidecars" sort_file_by_keys!(
            sidecar_refs, score_col, :target; reverse=[true, true])
        @score_phase "    A3.3 stream-merge" stream_sorted_merge(
            sidecar_refs, merged_path, score_col, :target;
            batch_size=batch_size, reverse=[true, true])
    finally
        GC.gc(false)
        for ref in sidecar_refs
            safeRm(file_path(ref), nothing; force=true)
        end
    end

    # Read the merged (score-descending) file ONCE and feed both the q-spline and the PEP interp,
    # instead of each re-opening it (T1: A3 double-read). `use_unique=false` here, so no per-column
    # dedup is needed — the two columns below are all both builders consume.
    scores, targets = @score_phase "    A3.4 read-merged-once" begin
        merged_df = DataFrame(Arrow.Table(merged_path))
        (Vector{Float32}(merged_df[!, score_col]), Vector{Bool}(merged_df[!, :target]))
    end
    qval_spline = @score_phase "    A3.5 qvalue-spline" _qvalue_spline_from_vectors(
        targets, scores;
        min_pep_points_per_bin=min_pep_points_per_bin, fdr_scale_factor=fdr_scale_factor)

    pep_interp = if compute_pep
        @score_phase "    A3.6 pep-interp" _pep_interp_from_vectors(
            scores, targets; fdr_scale_factor=fdr_scale_factor)
    else
        nothing
    end

    return (; qval_spline, pep_interp)
end
