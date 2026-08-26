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
Minimum PSMs per RT bin. A bin median is only meaningful when taken over enough
values; the previous fixed-bin-count scheme let this fall to 1, at which point
`median_rts` was just the sorted per-PSM retention times and the fitted
correction was noise. See dev_docs/quant_spline_guard.
"""
const QUANT_MIN_BIN_OCCUPANCY = 100

"""
Minimum fraction of runs in which a precursor must be quantified to serve as a
normalization anchor. With two or more runs, at least two observations are
always required; a one-run analysis uses its own precursors and therefore
reduces to the identity correction.
"""
const QUANT_MIN_ANCHOR_RUN_FRACTION = 0.5

"""
    _min_bins_for_spline(spline_n_knots) -> Int

Fewest RT bins that may be handed to `UniformSpline`.

One more than the coefficient count, so the least-squares system is
over-determined rather than exactly determined. The design matrix is full rank as
of the `_n_spline_coeffs` fix, so a square system is no longer a failure -- but
fitting `n` coefficients through exactly `n` points interpolates the bin medians
instead of smoothing them, which is not what this estimator is for.
"""
_min_bins_for_spline(spline_n_knots::Int) = _n_spline_coeffs(spline_n_knots) + 1

"""
    _occupancy_bins(n, min_occupancy, max_bins) -> Vector{UnitRange{Int}}

Contiguous bin ranges over `n` RT-sorted rows, each holding at least
`min_occupancy` rows, and at most `max_bins` bins.

Rather than forming `max_bins` bins and merging the undersized ones pairwise,
this computes how many bins the data can actually fill (`n ÷ min_occupancy`) and
lays them out evenly -- the same outcome with no merge loop to get wrong. Bin
sizes differ by at most one and every row is covered, unlike the previous
`bin_stop = bin_idx * bin_size` scheme, which silently dropped the final
`n mod max_bins` rows.

Returns a single range when `n` supports only one bin, and an empty vector when
`n < min_occupancy`; the caller treats both as "cannot fit a spline here".
"""
function _occupancy_bins(n::Int, min_occupancy::Int, max_bins::Int)
    (n <= 0 || min_occupancy <= 0 || max_bins <= 0) && return UnitRange{Int}[]
    n >= min_occupancy || return UnitRange{Int}[]
    k = min(max_bins, fld(n, min_occupancy))
    k <= 1 && return [1:n]
    base, rem = divrem(n, k)
    bins = Vector{UnitRange{Int}}(undef, k)
    start = 1
    @inbounds for i in 1:k
        len = base + (i <= rem ? 1 : 0)
        bins[i] = start:(start + len - 1)
        start += len
    end
    return bins
end

@inline function _required_anchor_runs(n_runs::Int, min_run_fraction::Real)
    n_runs <= 0 && return 0
    0.0 < min_run_fraction <= 1.0 || throw(ArgumentError(
        "min_run_fraction must be in (0, 1], got $min_run_fraction"
    ))
    n_runs == 1 && return 1
    return min(n_runs, max(2, ceil(Int, n_runs * min_run_fraction)))
end

@inline function _is_normalization_anchor_row(
    row::Int,
    target,
    mbr_recovered,
)
    if target !== nothing
        value = target[row]
        (ismissing(value) || !Bool(value)) && return false
    end
    if mbr_recovered !== nothing
        value = mbr_recovered[row]
        !ismissing(value) && Bool(value) && return false
    end
    return true
end

function _file_precursor_log2_quant(
    psms::DataFrame,
    quant_col_name::Symbol,
    precursor_col_name::Symbol,
)
    hasproperty(psms, precursor_col_name) || throw(ArgumentError(
        "Matched-precursor quant normalization requires column " *
        "$(precursor_col_name)."
    ))
    hasproperty(psms, quant_col_name) || throw(ArgumentError(
        "Quant normalization requires column $(quant_col_name)."
    ))

    precursor_idx = psms[!, precursor_col_name]
    quant = psms[!, quant_col_name]
    target = hasproperty(psms, :target) ? psms[!, :target] : nothing
    mbr_recovered = hasproperty(psms, :mbr_recovered) ?
        psms[!, :mbr_recovered] : nothing
    values = Dict{UInt32, Float32}()
    sizehint!(values, length(precursor_idx))

    @inbounds for i in eachindex(precursor_idx, quant)
        _is_normalization_anchor_row(i, target, mbr_recovered) || continue
        q_raw = quant[i]
        ismissing(q_raw) && continue
        q = Float64(q_raw)
        (isfinite(q) && q > 0.0) || continue
        pid = UInt32(precursor_idx[i])
        logq = Float32(log2(q))

        # Passing-PSM tables normally contain one row per precursor. If a
        # duplicate is present, keep one run-level value so that run cannot
        # receive extra weight in the cross-run precursor consensus.
        values[pid] = max(get(values, pid, -Inf32), logq)
    end
    return values
end

struct QuantRunAnchor
    log2_quant::Float32
    irt::Float32
end

function _file_precursor_quant_anchors(
    psms::DataFrame,
    quant_col_name::Symbol,
    precursor_col_name::Symbol,
)
    hasproperty(psms, :irt_obs) || throw(ArgumentError(
        "Pairwise quant normalization requires column irt_obs."
    ))
    hasproperty(psms, precursor_col_name) || throw(ArgumentError(
        "Pairwise quant normalization requires column $(precursor_col_name)."
    ))
    hasproperty(psms, quant_col_name) || throw(ArgumentError(
        "Quant normalization requires column $(quant_col_name)."
    ))

    precursor_idx = psms[!, precursor_col_name]
    quant = psms[!, quant_col_name]
    irt_obs = psms[!, :irt_obs]
    target = hasproperty(psms, :target) ? psms[!, :target] : nothing
    mbr_recovered = hasproperty(psms, :mbr_recovered) ?
        psms[!, :mbr_recovered] : nothing
    anchors = Dict{UInt32, QuantRunAnchor}()
    sizehint!(anchors, length(precursor_idx))

    @inbounds for row in eachindex(precursor_idx, quant, irt_obs)
        _is_normalization_anchor_row(row, target, mbr_recovered) || continue
        quant_raw = quant[row]
        irt_raw = irt_obs[row]
        (ismissing(quant_raw) || ismissing(irt_raw)) && continue
        abundance = Float64(quant_raw)
        irt = Float64(irt_raw)
        (isfinite(abundance) && abundance > 0.0 && isfinite(irt)) || continue

        pid = UInt32(precursor_idx[row])
        anchor = QuantRunAnchor(Float32(log2(abundance)), Float32(irt))
        previous = get(anchors, pid, nothing)
        if previous === nothing || anchor.log2_quant > previous.log2_quant
            anchors[pid] = anchor
        end
    end
    return anchors
end

"""
    getPrecursorQuantConsensus(psms_paths, quant_col_name;
        precursor_col_name=:precursor_idx,
        min_run_fraction=QUANT_MIN_ANCHOR_RUN_FRACTION)

Build a cross-run reference abundance for matched precursors. Each precursor's
reference is the median log2 abundance across runs in which it was quantified.
Only precursors observed in at least `min_run_fraction` of runs (and at least two
runs when multiple runs are present) are retained as normalization anchors.

The run-level de-duplication is deliberate: a precursor contributes at most one
value per run to the consensus.
"""
function getPrecursorQuantConsensus(
    psms_paths::Vector{String},
    quant_col_name::Symbol;
    precursor_col_name::Symbol = :precursor_idx,
    min_run_fraction::Real = QUANT_MIN_ANCHOR_RUN_FRACTION,
)
    required_runs = _required_anchor_runs(length(psms_paths), min_run_fraction)
    required_runs == 0 && return Dict{UInt32, Float32}()

    values_by_precursor = Dict{UInt32, Vector{Float32}}()
    for fpath in psms_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        file_values = _file_precursor_log2_quant(
            psms,
            quant_col_name,
            precursor_col_name,
        )
        for (pid, logq) in file_values
            push!(get!(() -> Float32[], values_by_precursor, pid), logq)
        end
    end

    consensus = Dict{UInt32, Float32}()
    sizehint!(consensus, length(values_by_precursor))
    for (pid, values) in values_by_precursor
        length(values) >= required_runs || continue
        consensus[pid] = Float32(median(values))
    end
    return consensus
end

"""
    getQuantSplines(psms_paths, quant_col_name; N=100, spline_n_knots=7,
                    min_bin_occupancy=QUANT_MIN_BIN_OCCUPANCY,
                    min_anchor_run_fraction=QUANT_MIN_ANCHOR_RUN_FRACTION)

Fit matched-precursor RT-dependent quantification splines for each MS file.
First builds a cross-run median log2 abundance for precursors observed in a
configurable fraction of runs. For each file, it then computes

`log2(file abundance) - precursor consensus`

for those matched anchors, bins the residuals by observed iRT into bins of at
least `min_bin_occupancy` anchors (at most `N` bins), and fits a `UniformSpline`
mapping iRT to the median matched-precursor residual.

Files whose matched-anchor count cannot support
`_min_bins_for_spline(spline_n_knots)` such bins get no entry in the returned
dictionary: they are left uncorrected and excluded from the cross-file median.
`applyNormalization!` handles the absent key.

Returns `(splines_dict, (min_rt, max_rt))` where `splines_dict` maps file path to
fitted spline, and the RT range spans only the files that got one.
"""
function getQuantSplines(psms_paths::Vector{String},
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7,
    min_bin_occupancy::Int = QUANT_MIN_BIN_OCCUPANCY,
    min_anchor_run_fraction::Real = QUANT_MIN_ANCHOR_RUN_FRACTION,
    precursor_col_name::Symbol = :precursor_idx)
    precursor_consensus = getPrecursorQuantConsensus(
        psms_paths,
        quant_col_name;
        precursor_col_name = precursor_col_name,
        min_run_fraction = min_anchor_run_fraction,
    )
    quant_splines = Dictionary{String, UniformSpline}()
    min_rt, max_rt = typemax(Float32), typemin(Float32)
    min_bins = _min_bins_for_spline(spline_n_knots)
    for fpath in psms_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        hasproperty(psms, :irt_obs) || throw(ArgumentError(
            "Matched-precursor quant normalization requires column irt_obs."
        ))
        file_values = _file_precursor_log2_quant(
            psms,
            quant_col_name,
            precursor_col_name,
        )

        anchor_rts = Float64[]
        anchor_residuals = Float64[]
        sizehint!(anchor_rts, min(length(file_values), length(precursor_consensus)))
        sizehint!(anchor_residuals, min(length(file_values), length(precursor_consensus)))
        precursor_idx = psms[!, precursor_col_name]
        irt_obs = psms[!, :irt_obs]
        target = hasproperty(psms, :target) ? psms[!, :target] : nothing
        mbr_recovered = hasproperty(psms, :mbr_recovered) ?
            psms[!, :mbr_recovered] : nothing
        seen = Set{UInt32}()
        sizehint!(seen, length(file_values))
        @inbounds for i in eachindex(precursor_idx, irt_obs)
            _is_normalization_anchor_row(i, target, mbr_recovered) || continue
            pid = UInt32(precursor_idx[i])
            pid in seen && continue
            reference_logq = get(precursor_consensus, pid, nothing)
            reference_logq === nothing && continue
            file_logq = get(file_values, pid, nothing)
            file_logq === nothing && continue
            rt_raw = irt_obs[i]
            ismissing(rt_raw) && continue
            rt = Float64(rt_raw)
            isfinite(rt) || continue
            push!(seen, pid)
            push!(anchor_rts, rt)
            push!(anchor_residuals, Float64(file_logq - reference_logq))
        end

        order = sortperm(anchor_rts)
        anchor_rts = anchor_rts[order]
        anchor_residuals = anchor_residuals[order]
        nanchors = length(anchor_rts)

        bins = _occupancy_bins(nanchors, min_bin_occupancy, N)
        if length(bins) < min_bins
            @user_warn "Skipping quant normalization for $(basename(fpath)): " *
                       "$nanchors matched precursor anchors support only " *
                       "$(length(bins)) RT bin(s) at >= $min_bin_occupancy " *
                       "anchors each, and a " *
                       "$(_n_spline_coeffs(spline_n_knots))-coefficient spline needs at least " *
                       "$min_bins. Abundances are left uncorrected and excluded " *
                       "from the cross-file median."
            continue
        end

        # Only files that get a spline define the correction grid. A skipped file
        # receives no correction, so widening the grid to its RT range would only
        # extrapolate the surviving splines further than their data supports.
        if first(anchor_rts) < min_rt
            min_rt = first(anchor_rts)
        end
        if last(anchor_rts) > max_rt
            max_rt = last(anchor_rts)
        end

        median_quant = Vector{Float64}(undef, length(bins))
        median_rts = Vector{Float64}(undef, length(bins))
        for (i, b) in enumerate(bins)
            median_rts[i] = median(@view(anchor_rts[b]))
            median_quant[i] = median(@view(anchor_residuals[b]))
        end

        splinefit = UniformSpline(median_quant, median_rts, 3, spline_n_knots)
        insert!(quant_splines, fpath, splinefit)
    end
    return quant_splines, (min_rt, max_rt)
end

struct PairwiseQuantSpline
    left_position::Int
    right_position::Int
    similarity::Float32
    nanchors::Int
    min_rt::Float64
    max_rt::Float64
    spline::UniformSpline
end

function _fit_pairwise_quant_spline(
    left_position::Int,
    right_position::Int,
    similarity::Float32,
    left_anchors::Dict{UInt32, QuantRunAnchor},
    right_anchors::Dict{UInt32, QuantRunAnchor};
    N::Int,
    spline_n_knots::Int,
    min_bin_occupancy::Int,
)
    anchor_rts = Float64[]
    anchor_residuals = Float64[]
    anchor_pair_log2_abundances = Float64[]
    npossible = min(length(left_anchors), length(right_anchors))
    sizehint!(anchor_rts, npossible)
    sizehint!(anchor_residuals, npossible)
    sizehint!(anchor_pair_log2_abundances, npossible)

    iterated_anchors, lookup_anchors = length(left_anchors) <= length(right_anchors) ?
        (left_anchors, right_anchors) : (right_anchors, left_anchors)
    @inbounds for pid in keys(iterated_anchors)
        haskey(lookup_anchors, pid) || continue
        left_anchor = left_anchors[pid]
        right_anchor = right_anchors[pid]
        push!(anchor_rts, (Float64(left_anchor.irt) + Float64(right_anchor.irt)) / 2.0)
        push!(
            anchor_residuals,
            Float64(left_anchor.log2_quant - right_anchor.log2_quant),
        )
        push!(
            anchor_pair_log2_abundances,
            (Float64(left_anchor.log2_quant) + Float64(right_anchor.log2_quant)) / 2.0,
        )
    end

    order = sortperm(anchor_rts)
    anchor_rts = anchor_rts[order]
    anchor_residuals = anchor_residuals[order]
    anchor_pair_log2_abundances = anchor_pair_log2_abundances[order]
    nanchors = length(anchor_rts)
    bins = _occupancy_bins(nanchors, min_bin_occupancy, N)
    length(bins) >= _min_bins_for_spline(spline_n_knots) || return nothing

    median_residuals = Vector{Float64}(undef, length(bins))
    median_rts = Vector{Float64}(undef, length(bins))
    for (bin_idx, bin) in enumerate(bins)
        median_rts[bin_idx] = median(@view(anchor_rts[bin]))

        # Rank abundance within this matched pair and RT neighborhood. The
        # average log2 abundance is the log2 geometric mean of the two raw
        # intensities, so it measures combined signal without directly
        # weighting on the left-minus-right ratio. Give the middle of that
        # distribution the most influence: low-abundance ratios can be biased
        # by detection limits, while the highest-abundance ratios can be
        # compressed by saturation or other high-signal effects.
        pair_abundances = @view(anchor_pair_log2_abundances[bin])
        sorted_pair_abundances = sort!(collect(pair_abundances))
        bin_weights = Vector{Float64}(undef, length(bin))
        for local_idx in eachindex(pair_abundances)
            lower_rank = searchsortedfirst(
                sorted_pair_abundances,
                pair_abundances[local_idx],
            )
            upper_rank = searchsortedlast(
                sorted_pair_abundances,
                pair_abundances[local_idx],
            )
            percentile = (
                (lower_rank + upper_rank) / 2.0 - 0.5
            ) / length(bin)
            bin_weights[local_idx] = 4.0 * percentile * (1.0 - percentile)
        end
        median_residuals[bin_idx] = median(
            @view(anchor_residuals[bin]),
            weights(bin_weights),
        )
    end

    return PairwiseQuantSpline(
        left_position,
        right_position,
        similarity,
        nanchors,
        first(anchor_rts),
        last(anchor_rts),
        UniformSpline(median_residuals, median_rts, 3, spline_n_knots),
    )
end

@inline function _quant_tree_root!(parent::Vector{Int}, node::Int)
    root = node
    while parent[root] != root
        root = parent[root]
    end
    while parent[node] != node
        next_node = parent[node]
        parent[node] = root
        node = next_node
    end
    return root
end

"""
    getPairwiseQuantTree(psms_paths, quant_col_name, run_ids,
                         run_similarity_atlas; ...)

Build a maximum-similarity spanning forest for pairwise quant normalization.
Candidate run pairs are ordered by the maximum directional containment in the
existing run-similarity atlas. An edge is accepted only when exact, non-MBR,
target precursor matches can support the requested RT spline. Within each RT
bin, matched precursor ratios receive middle-peaked weights based on the
percentile rank of their pairwise mean log2 abundance (equivalently, their
raw-intensity geometric mean). This downweights both detection-limited and
potentially high-signal-compressed measurements without directly weighting on
the ratio being estimated. Kruskal's algorithm therefore returns a maximum
spanning tree when the supported graph is connected, or a maximum spanning
forest otherwise.

Only candidate edges that could connect two current components are fitted, so
the expensive exact-match work normally stops after `n_runs - 1` successful
pairwise fits.
"""
function getPairwiseQuantTree(
    psms_paths::Vector{String},
    quant_col_name::Symbol,
    run_ids::Vector{UInt32},
    run_similarity_atlas;
    N::Int = 100,
    spline_n_knots::Int = 7,
    min_bin_occupancy::Int = QUANT_MIN_BIN_OCCUPANCY,
    precursor_col_name::Symbol = :precursor_idx,
)
    length(psms_paths) == length(run_ids) || throw(DimensionMismatch(
        "psms_paths and run_ids must have the same length"
    ))
    length(unique(run_ids)) == length(run_ids) || throw(ArgumentError(
        "run_ids must be unique"
    ))
    run_similarity_atlas === nothing && throw(ArgumentError(
        "pairwise quant normalization requires a run-similarity atlas"
    ))

    run_anchors = Vector{Dict{UInt32, QuantRunAnchor}}(undef, length(psms_paths))
    for (position, fpath) in enumerate(psms_paths)
        psms = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        run_anchors[position] = _file_precursor_quant_anchors(
            psms,
            quant_col_name,
            precursor_col_name,
        )
    end

    candidates = Tuple{Float32, Int, Int}[]
    n_runs = length(psms_paths)
    n_runs <= 1 && return PairwiseQuantSpline[]
    sizehint!(candidates, n_runs * max(n_runs - 1, 0) ÷ 2)
    for left_position in 1:(n_runs - 1)
        left_run = run_ids[left_position]
        for right_position in (left_position + 1):n_runs
            right_run = run_ids[right_position]
            similarity = max(
                run_similarity(run_similarity_atlas, left_run, right_run),
                run_similarity(run_similarity_atlas, right_run, left_run),
            )
            push!(candidates, (similarity, left_position, right_position))
        end
    end
    sort!(
        candidates;
        by = candidate -> (
            -Float64(candidate[1]),
            run_ids[candidate[2]],
            run_ids[candidate[3]],
        ),
    )

    parent = collect(1:n_runs)
    component_size = ones(Int, n_runs)
    tree = PairwiseQuantSpline[]
    sizehint!(tree, max(n_runs - 1, 0))
    for (similarity, left_position, right_position) in candidates
        left_root = _quant_tree_root!(parent, left_position)
        right_root = _quant_tree_root!(parent, right_position)
        left_root == right_root && continue

        pairwise_spline = _fit_pairwise_quant_spline(
            left_position,
            right_position,
            similarity,
            run_anchors[left_position],
            run_anchors[right_position];
            N = N,
            spline_n_knots = spline_n_knots,
            min_bin_occupancy = min_bin_occupancy,
        )
        pairwise_spline === nothing && continue
        push!(tree, pairwise_spline)

        if component_size[left_root] < component_size[right_root]
            left_root, right_root = right_root, left_root
        end
        parent[right_root] = left_root
        component_size[left_root] += component_size[right_root]
        length(tree) == n_runs - 1 && break
    end
    return tree
end

function _quant_tree_traversal(
    n_runs::Int,
    tree::Vector{PairwiseQuantSpline},
)
    adjacency = [Tuple{Int, Int}[] for _ in 1:n_runs]
    for (edge_idx, edge) in enumerate(tree)
        push!(adjacency[edge.left_position], (edge.right_position, edge_idx))
        push!(adjacency[edge.right_position], (edge.left_position, edge_idx))
    end

    parent = zeros(Int, n_runs)
    parent_edge = zeros(Int, n_runs)
    components = Vector{Vector{Int}}()
    orders = Vector{Vector{Int}}()
    for root in 1:n_runs
        parent[root] == 0 || continue
        parent[root] = root
        component = Int[root]
        order = Int[root]
        next_idx = 1
        while next_idx <= length(order)
            node = order[next_idx]
            next_idx += 1
            for (neighbor, edge_idx) in adjacency[node]
                parent[neighbor] == 0 || continue
                parent[neighbor] = node
                parent_edge[neighbor] = edge_idx
                push!(component, neighbor)
                push!(order, neighbor)
            end
        end
        push!(components, component)
        push!(orders, order)
    end
    return parent, parent_edge, components, orders
end

"""
    getPairwiseQuantCorrections(psms_paths, tree; N=100)

Propagate pairwise log2 differences through a maximum spanning forest and
return per-file correction functions. For every RT grid point, each connected
component is centered independently at its median; this fixes the additive
degree of freedom without inventing a relationship between disconnected run
groups.

Returns `(corrections, components)` where components contain path positions.
"""
function getPairwiseQuantCorrections(
    psms_paths::Vector{String},
    tree::Vector{PairwiseQuantSpline};
    N::Int = 100,
)
    n_runs = length(psms_paths)
    n_runs == 0 && return Dictionary{String, Any}(), Vector{Vector{Int}}()
    parent, parent_edge, components, orders = _quant_tree_traversal(n_runs, tree)
    isempty(tree) && return Dictionary{String, Any}(), components
    N >= 2 || throw(ArgumentError("N must be at least 2, got $N"))

    min_rt = minimum(edge.min_rt for edge in tree)
    max_rt = maximum(edge.max_rt for edge in tree)
    rt_grid = collect(LinRange(min_rt, max_rt, N))
    offsets = zeros(Float64, n_runs, N)
    for (rt_idx, rt) in enumerate(rt_grid)
        for (component, order) in zip(components, orders)
            root = first(order)
            offsets[root, rt_idx] = 0.0
            for node in @view(order[2:end])
                parent_node = parent[node]
                edge = tree[parent_edge[node]]
                edge_rt = clamp(rt, edge.min_rt, edge.max_rt)
                difference = Float64(edge.spline(edge_rt))
                if parent_node == edge.left_position
                    offsets[node, rt_idx] = offsets[parent_node, rt_idx] - difference
                else
                    offsets[node, rt_idx] = offsets[parent_node, rt_idx] + difference
                end
            end
            center = median(@view(offsets[component, rt_idx]))
            for node in component
                offsets[node, rt_idx] -= center
            end
        end
    end

    corrections = Dictionary{String, Any}()
    for (position, fpath) in enumerate(psms_paths)
        insert!(corrections, fpath, linear_interpolation(
            rt_grid,
            copy(@view(offsets[position, :]));
            extrapolation_bc = Interpolations.Flat(),
        ))
    end
    return corrections, components
end

"""
    getQuantCorrections(quant_splines, rt_range; N=100)

Compute per-file RT-dependent correction offsets. Evaluates each file's spline
on a grid, computes the cross-file median at each RT point, and returns a
dictionary of interpolated offset functions (file_spline - global_median).

Returns an empty dictionary when no file got a spline: there is no cross-file
median to correct toward, and `rt_range` is still the empty-range sentinel
`(Inf, -Inf)`. `applyNormalization!` reads that as "no correction for this file"
and copies the abundances through unchanged (correction factor 1.0).
"""
function getQuantCorrections(
    quant_splines::Dictionary{String, UniformSpline},
    rt_range::Tuple{AbstractFloat, AbstractFloat};
    N = 100)
    isempty(quant_splines) && return Dictionary{String, Any}()
    median_quant = zeros(Float32, (length(quant_splines), N))
    rt_grid = collect(LinRange(first(rt_range), last(rt_range), N))
    for (i, rt) in enumerate(rt_grid)
        j = 1
        for (key, spline) in pairs(quant_splines)
            median_quant[j, i] = spline(rt_grid[i])
            j += 1
        end
    end
    median_quant = reshape(median(median_quant, dims = 1), (N,))
    median_interp = linear_interpolation(rt_grid, median_quant)
    corrections = Dictionary{String, Any}()
    for (key, spline) in pairs(quant_splines)
        offset = zeros(Float32, N)
        for (i, rt) in enumerate(rt_grid)
            offset[i] = spline(rt) - median_interp(rt)
        end
        # Censoring can leave no matched anchors at one or both RT boundaries.
        # Rows outside the anchor-supported range still need a correction; hold
        # the nearest supported value constant instead of extrapolating a trend
        # that was never observed (or throwing on an out-of-range lookup).
        insert!(corrections, key, linear_interpolation(
            rt_grid,
            offset;
            extrapolation_bc = Interpolations.Flat(),
        ))
    end
    return corrections
end

"""
    applyNormalization!(psms_paths, quant_col, corrections)

Apply RT-dependent normalization corrections to Arrow files in place.
For each PSM, subtracts the file-specific log2-offset at its observed iRT,
writing a new `{quant_col}_normalized` column.
"""
function applyNormalization!(
    psms_paths::Vector{String},
    quant_col::Symbol,
    corrections::Dictionary{String, Any}
)
    for fpath in psms_paths
        psms = DataFrame(Tables.columntable(Arrow.Table(fpath)))
        norm_quant_col = Symbol(string(quant_col) * "_normalized")
        # getQuantSplines skips files that cannot support a spline, so this key may
        # be absent. Emit the column unchanged rather than omitting it: downstream
        # readers assume every file carries the same schema.
        correction_spline = get(corrections, fpath, nothing)
        if correction_spline === nothing
            psms[!, norm_quant_col] = Float64.(psms[!, quant_col])
            writeArrow(fpath, psms)
            continue
        end
        n = size(psms, 1)
        norm_vals = Vector{Float64}(undef, n)
        for i in range(1, n)
            hc = correction_spline(psms[i, :irt_obs])
            norm_vals[i] = 2^(log2(max(psms[i, quant_col], 0.0)) - hc)
        end
        psms[!, norm_quant_col] = norm_vals
        writeArrow(fpath, psms)
    end
end

"""
    normalizeQuant(psms_paths, quant_col_name; N=100, spline_n_knots=7,
                   run_ids=nothing, run_similarity_atlas=nothing)

End-to-end RT-dependent quantification normalization. When a run-similarity
atlas is supplied, fits exact-match pairwise splines only along a supported
maximum spanning tree, propagates their differences to run corrections, and
median-centers the tree at each RT. MBR-recovered rows never serve as anchors.

Without an atlas, retains the experiment-wide matched-precursor estimator as a
fallback.

This corrects systematic RT-dependent intensity differences between MS runs
without forcing the marginal observed intensity distributions to match.
"""
function normalizeQuant(
    psms_paths::Vector{String},
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7,
    min_anchor_run_fraction::Real = QUANT_MIN_ANCHOR_RUN_FRACTION,
    run_ids::Union{Nothing, Vector{UInt32}} = nothing,
    run_similarity_atlas = nothing)

    if run_similarity_atlas !== nothing
        resolved_run_ids = run_ids === nothing ?
            UInt32.(eachindex(psms_paths)) : run_ids
        pairwise_tree = getPairwiseQuantTree(
            psms_paths,
            quant_col_name,
            resolved_run_ids,
            run_similarity_atlas;
            N = N,
            spline_n_knots = spline_n_knots,
        )
        if !isempty(pairwise_tree) || length(psms_paths) <= 1
            quant_corrections, components = getPairwiseQuantCorrections(
                psms_paths,
                pairwise_tree;
                N = N,
            )
            if length(components) > 1
                @user_warn "Pairwise quant normalization formed " *
                           "$(length(components)) disconnected run groups. " *
                           "Each group is normalized and centered independently."
            end
            applyNormalization!(psms_paths, quant_col_name, quant_corrections)
            return nothing
        end

        @user_warn "No run pair had enough exact non-MBR precursor matches " *
                   "to support an RT spline. Falling back to the " *
                   "experiment-wide matched-precursor normalizer."
    end

    quant_splines_dict, rt_range = getQuantSplines(
        psms_paths,
        quant_col_name;
        N = N,
        spline_n_knots = spline_n_knots,
        min_anchor_run_fraction = min_anchor_run_fraction,
    )

    quant_corrections_dict = getQuantCorrections(
        quant_splines_dict, rt_range, N = N)

    applyNormalization!(psms_paths, quant_col_name, quant_corrections_dict)
    return nothing
end

function normalizeQuant(
    second_quant_folder::String,
    quant_col_name::Symbol;
    N::Int = 100,
    spline_n_knots::Int = 7,
    min_anchor_run_fraction::Real = QUANT_MIN_ANCHOR_RUN_FRACTION,
    run_ids::Union{Nothing, Vector{UInt32}} = nothing,
    run_similarity_atlas = nothing)

    psms_paths = [fpath for fpath in readdir(second_quant_folder, join=true) if endswith(fpath, ".arrow")]
    normalizeQuant(
        psms_paths,
        quant_col_name;
        N = N,
        spline_n_knots = spline_n_knots,
        min_anchor_run_fraction = min_anchor_run_fraction,
        run_ids = run_ids,
        run_similarity_atlas = run_similarity_atlas,
    )
end
