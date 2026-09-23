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


"""
    _annotate_precursor_scores_via_sidecar!(refs, global_prob_dict, global_qval_dict,
                                            global_pep_dict, qval_spline, pep_interp) -> refs

Attach the five per-row precursor score annotations as a row-aligned sidecar instead of rewriting
every file. All five are additive, so no full materialisation is required here; the next pass that
removes rows (the initial q-value filter) consolidates them via load_with_sidecars.

Reads only :precursor_idx and :prec_prob -- :prec_prob itself lives in a sidecar after
aggregate_per_file!, which materialize_columns resolves.
"""
function _annotate_precursor_scores_via_sidecar!(
    refs::Vector{PSMFileReference},
    global_prob_dict::Dict{UInt32, Float32},
    global_qval_dict::Dict{UInt32, Float32},
    global_pep_dict::Dict{UInt32, Float32},
    qval_spline,
    pep_interp,
)
    started = last_progress = time()
    rows_processed = 0
    @debug_l1 "Precursor score annotation starting: files=$(length(refs))"
    for (file_idx, ref) in enumerate(refs)
        exists(ref) || continue
        cols = materialize_columns(ref, Symbol[:precursor_idx, :prec_prob])
        pids = cols[!, :precursor_idx]
        probs = cols[!, :prec_prob]
        n = length(pids)

        # add_dict_column produces Vector{Union{V,Missing}} for absent keys; match that exactly.
        global_prob = Vector{Union{Float32, Missing}}(undef, n)
        global_qval = Vector{Union{Float32, Missing}}(undef, n)
        global_pep = Vector{Union{Float32, Missing}}(undef, n)
        qvals = Vector{Float32}(undef, n)
        peps = Vector{Float32}(undef, n)
        @inbounds for row in 1:n
            pid = UInt32(pids[row])
            global_prob[row] = get(global_prob_dict, pid, missing)
            global_qval[row] = get(global_qval_dict, pid, missing)
            global_pep[row] = get(global_pep_dict, pid, missing)
            score = Float32(probs[row])
            qvals[row] = Float32(qval_spline(score))
            peps[row] = Float32(pep_interp(score))
        end

        add_columns_via_sidecar!(
            ref,
            :global_prob => global_prob,
            :global_qval => global_qval,
            :global_pep => global_pep,
            :qval => qvals,
            :pep => peps;
            tag = "precursor_scores",
        )
        rows_processed += n
        if time() - last_progress >= 60
            @debug_l1 "Precursor score annotation: files=$file_idx/$(length(refs)) rows=$rows_processed elapsed=$(round(time() - started, digits=2))s"
            last_progress = time()
        end
    end
    @debug_l1 "Precursor score annotation complete: rows=$rows_processed elapsed=$(round(time() - started, digits=2))s"
    return refs
end

"""
    aggregate_per_file!(refs)

Per-file precursor probability aggregation (no MBR filtering).
"""
function aggregate_per_file!(refs::Vector{PSMFileReference})
    started = last_progress = time()
    rows_processed = 0
    @debug_l1 "Precursor probability aggregation starting: files=$(length(refs))"
    for (file_idx, ref) in enumerate(refs)
        df = load_with_sidecars(ref)
        _aggregate_trace_to_precursor_probs!(df)
        # Write only :prec_prob as a row-aligned sidecar instead of rewriting
        # the entire main file. Downstream readers locate :prec_prob via the
        # PSMFileReference's sidecar registry.
        side_path = file_path(ref) * ".prec_prob.sidecar.arrow"
        writeArrow(side_path, DataFrame(prec_prob = df.prec_prob))
        register_sidecar!(ref, side_path, [:prec_prob])
        rows_processed += nrow(df)
        if time() - last_progress >= 60
            @debug_l1 "Precursor probability aggregation: files=$file_idx/$(length(refs)) rows=$rows_processed elapsed=$(round(time() - started, digits=2))s"
            last_progress = time()
        end
    end
    @debug_l1 "Precursor probability aggregation complete: files=$(length(refs)) rows=$rows_processed elapsed=$(round(time() - started, digits=2))s"
    return nothing
end

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
    _score_floor_for_qvalue(qval_spline, q_value_threshold)

Find the lowest run-level score whose pooled experiment-wide q-value passes
the requested threshold.
"""
function _score_floor_for_qvalue(
    qval_spline,
    q_value_threshold::Float32,
)::Float32
    low = eps(Float32)
    high = 1.0f0 - eps(Float32)
    Float32(qval_spline(low)) <= q_value_threshold && return low
    Float32(qval_spline(high)) > q_value_threshold && return high

    for _ in 1:32
        midpoint = (low + high) / 2.0f0
        if Float32(qval_spline(midpoint)) <= q_value_threshold
            high = midpoint
        else
            low = midpoint
        end
    end
    return high
end

"""
    build_qvalue_spline_from_refs(refs, score_col, merged_path; ...) → Union{Nothing, NamedTuple}

Build grouped q-value/PEP mappings using bounded sorting and disk-backed calibration.
`merged_path` selects the scratch directory; batch size is retained for caller compatibility.
"""
function build_qvalue_spline_from_refs(
    refs::Vector{<:FileReference},
    score_col::Symbol,
    merged_path::String;
    batch_size::Int = 10_000_000,
    compute_pep::Bool = false,
    fdr_scale_factor::Float32 = 1.0f0,
    temp_prefix::String = "sidecar",
    memory_budget_bytes::Int = SCORE_WORKSPACE_BYTES,
)
    started = time()
    context = "Score calibration ($temp_prefix, $score_col)"
    rows = 0
    @debug_l1 "$context grouping starting: files=$(length(refs))"
    result = build_score_calibration(; compute_pep, fdr_scale_factor, memory_budget_bytes,
        temp_parent=dirname(merged_path)) do emit
        for ref in refs
            if ref isa PSMFileReference
                table = materialize_columns(ref, [score_col, :target])
                _emit_score_arrays(emit, table[!, score_col], table[!, :target])
                rows += nrow(table)
            else
                for table in Arrow.Stream(file_path(ref))
                    scores, targets = Tables.getcolumn(table, score_col), Tables.getcolumn(table, :target)
                    _emit_score_arrays(emit, scores, targets)
                    rows += length(scores)
                end
            end
        end
    end
    @debug_l1 "$context complete: files=$(length(refs)) rows=$rows elapsed=$(round(time()-started, digits=2))s"
    return result
end
