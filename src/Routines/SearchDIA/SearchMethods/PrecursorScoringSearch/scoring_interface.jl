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
end

"""
    aggregate_per_file!(refs)

Per-file precursor probability aggregation (no MBR filtering).
"""
function aggregate_per_file!(refs::Vector{PSMFileReference})
    for ref in refs
        df = DataFrame(Tables.columntable(Arrow.Table(file_path(ref))))
        _aggregate_trace_to_precursor_probs!(df)
        write_arrow_file(ref, df)
    end
    return nothing
end

#==========================================================
Additional Interface Functions
==========================================================#

"""
    logodds(probs::AbstractVector{T}, top_n::Int) where {T<:AbstractFloat}

Combine probabilities using a log-odds average.
The final value is converted back to a probability via the logistic function.
"""
function logodds(probs::AbstractVector{T}, top_n::Int) where {T<:AbstractFloat}
    isempty(probs) && return 0.0f0
    n = min(length(probs), top_n)
    # Sort descending and select the top n probabilities
    sorted = sort(probs; rev=true)
    selected = sorted[1:n]
    eps = 1f-6
    # Convert to log-odds, clip to avoid Inf or negative contribution
    logodds = log.(clamp.(selected, 0.1f0, 1 - eps) ./ (1 .- clamp.(selected, 0.1f0, 1 - eps)))
    avg = sum(logodds) / n
    return 1.0f0 / (1 + exp(-avg))
end

#==========================================================
Dictionary + Sidecar Helper Functions for OOM Scoring Pipeline
==========================================================#

"""
    build_precursor_global_prob_dicts(refs, sqrt_n_runs, n_precursors)
    → (global_prob_dict, target_dict)

Stream per-file to build precursor_idx → global_prob dictionaries.
Reads only (precursor_idx, prec_prob, target) via mmap.
"""
function build_precursor_global_prob_dicts(
    refs::Vector{PSMFileReference},
    sqrt_n_runs::Int,
    n_precursors::Int
)
    # Pre-allocate accumulation dictionaries with known upper bound
    prob_acc = Dict{UInt32, Vector{Float32}}()
    sizehint!(prob_acc, n_precursors)
    target_dict = Dict{UInt32, Bool}()
    sizehint!(target_dict, n_precursors)

    n_swapped = 0
    n_total_pushed = 0
    for ref in refs
        tbl = Arrow.Table(file_path(ref))
        n = length(tbl.precursor_idx)
        n == 0 && continue
        prec_ids = tbl.precursor_idx
        prec_probs = tbl.prec_prob
        targets = tbl.target

        # Phase 8h: hybrid global aggregation. For rows with :mbr_recovered=true,
        # contribute the MBR-boosted :trace_prob_mbr instead of the non-MBR
        # :prec_prob to the global logodds. Lets B1-style precursors clear
        # global_qval ≤ 0.01 when FTR rescued them in some files. Both T←T
        # targets and D←D / D←T decoys recovered by FTR get the swap →
        # target/decoy symmetry preserved at the global level. Falls back to
        # plain :prec_prob if either column is absent (e.g., MBR off).
        has_boost  = hasproperty(tbl, :trace_prob_mbr) && hasproperty(tbl, :mbr_recovered)
        boosted    = has_boost ? tbl.trace_prob_mbr : nothing
        recovered  = has_boost ? tbl.mbr_recovered  : nothing

        @inbounds for i in 1:n
            pid = prec_ids[i]
            if !haskey(prob_acc, pid)
                prob_acc[pid] = Float32[]
                target_dict[pid] = targets[i]
            end
            score = if has_boost && recovered[i]
                n_swapped += 1
                Float32(boosted[i])
            else
                prec_probs[i]
            end
            n_total_pushed += 1
            push!(prob_acc[pid], score)
        end
    end

    # Compute logodds per precursor
    global_prob_dict = Dict{UInt32, Float32}()
    sizehint!(global_prob_dict, length(prob_acc))
    for (pid, probs) in prob_acc
        global_prob_dict[pid] = logodds(probs, sqrt_n_runs)
    end

    pct = n_total_pushed == 0 ? 0.0 : round(100 * n_swapped / n_total_pushed, digits=2)
    @user_info "Phase 8h global aggregator: $n_swapped of $n_total_pushed contributions ($pct%) used :trace_prob_mbr (MBR-recovered)"

    return global_prob_dict, target_dict
end

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
    targets = Bool[get(target_dict, pid, false) for pid in pids]

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
        tbl = Arrow.Table(file_path(ref))
        n = length(Tables.getcolumn(tbl, first(columns)))
        n == 0 && continue

        col_data = NamedTuple{Tuple(columns)}(Tuple(collect(Tables.getcolumn(tbl, c)) for c in columns))
        temp_path = tempname() * "_$(temp_prefix).arrow"
        writeArrow(temp_path, DataFrame(; col_data...))
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
    sidecar_refs = write_score_sidecars(refs, [score_col, :target]; temp_prefix=temp_prefix)
    isempty(sidecar_refs) && return nothing

    try
        sort_file_by_keys!(sidecar_refs, score_col, :target; reverse=[true, true])
        stream_sorted_merge(sidecar_refs, merged_path, score_col, :target;
                           batch_size=batch_size, reverse=[true, true])
    finally
        GC.gc(false)
        for ref in sidecar_refs
            safeRm(file_path(ref), nothing; force=true)
        end
    end

    qval_spline = get_qvalue_spline(merged_path, score_col, false;
        min_pep_points_per_bin=min_pep_points_per_bin,
        fdr_scale_factor=fdr_scale_factor)

    pep_interp = if compute_pep
        get_pep_interpolation(merged_path, score_col;
            fdr_scale_factor=fdr_scale_factor)
    else
        nothing
    end

    return (; qval_spline, pep_interp)
end
