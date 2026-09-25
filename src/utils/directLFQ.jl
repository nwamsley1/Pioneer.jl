# Julia adaptation of MannLabs/directlfq (Apache-2.0), revision
# c1b5b650b61557c00fafc1f436eeeff6552b25de. See licenses/directlfq-LICENSE.

const DIRECTLFQ_REFERENCE_REVISION = "c1b5b650b61557c00fafc1f436eeeff6552b25de"
const DIRECTLFQ_MAX_PRECURSORS = 100
const DIRECTLFQ_REFERENCE_PRECURSORS = 10

function directlfq_distance!(buffer, X, i, j)
    empty!(buffer)
    for r in axes(X, 2)
        a, b = X[i, r], X[j, r]
        isfinite(a) && isfinite(b) && push!(buffer, a - b)
    end
    isempty(buffer) && return (Inf, Inf)
    return (median!(buffer), var(buffer; corrected=false))
end

function directlfq_connected_mask(distances)
    n = size(distances, 1)
    degrees = zeros(Int, n)
    for i in 1:n, j in i+1:n
        if isfinite(distances[i, j])
            degrees[i] += 1
            degrees[j] += 1
        end
    end
    visited = falses(n)
    n == 0 && return visited
    stack = [argmax(degrees)]
    while !isempty(stack)
        i = pop!(stack)
        for j in 1:n
            if !visited[j] && isfinite(distances[min(i, j), max(i, j)])
                visited[j] = true
                push!(stack, j)
            end
        end
    end
    return visited
end

function directlfq_align_reference!(X, buffer)
    n, nruns = size(X)
    for i in 1:n
        if count(isfinite, @view X[i, :]) < 2
            X[i, :] .= NaN
        end
    end
    merged = copy(X)
    distances = fill(Inf, n, n)
    variances = fill(Inf, n, n)
    for i in 1:n, j in i+1:n
        distances[i, j], variances[i, j] = directlfq_distance!(buffer, X, i, j)
    end
    connected = directlfq_connected_mask(distances)
    for i in 1:n
        if !connected[i]
            distances[i, :] .= Inf
            distances[:, i] .= Inf
            variances[i, :] .= Inf
            variances[:, i] .= Inf
        end
    end
    counts = ones(Int, n)
    parents = zeros(Int, n)
    shifts = zeros(Float64, n)
    for _ in 1:n-1
        best_variance, left, right = Inf, 0, 0
        # Row-major tie order matches the reference implementation.
        for i in 1:n, j in i+1:n
            if variances[i, j] < best_variance
                best_variance, left, right = variances[i, j], i, j
            end
        end
        left == 0 && break
        anchor, shifted = counts[left] >= counts[right] ? (left, right) : (right, left)
        shift = distances[left, right] * (anchor == left ? 1 : -1)
        parents[shifted] = anchor
        shifts[shifted] = shift
        for r in 1:nruns
            a, b = merged[anchor, r], X[shifted, r] + shift
            merged[anchor, r] = !isfinite(a) ? b : !isfinite(b) ? a :
                (a * counts[anchor] + b * counts[shifted]) / (counts[anchor] + counts[shifted])
        end
        for other in 1:n
            other == anchor && continue
            i, j = minmax(anchor, other)
            if isfinite(distances[i, j])
                distances[i, j], variances[i, j] = directlfq_distance!(buffer, merged, i, j)
            end
        end
        distances[shifted, :] .= Inf
        distances[:, shifted] .= Inf
        variances[shifted, :] .= Inf
        variances[:, shifted] .= Inf
        counts[anchor] += 1
    end
    for i in 1:n
        total_shift, node = 0.0, i
        while parents[node] != 0
            total_shift += shifts[node]
            node = parents[node]
        end
        X[i, :] .+= total_shift
    end
    return X
end

"""
    solve_directlfq(X; precursor_ids=collect(axes(X, 1)))

Estimate log2 protein abundances from a precursor-by-run matrix of log2 areas.
Returns abundances, selected input row indices, and the number of contributing
aligned precursors per run. Missing/nonfinite entries are not quantified.
Uses the pinned directLFQ defaults: at most 100 precursors and 10 anchor traces.
Between-run normalization is performed separately by Pioneer.
"""
function solve_directlfq(X::AbstractMatrix; precursor_ids=collect(axes(X, 1)))
    n, nruns = size(X)
    length(precursor_ids) == n || throw(DimensionMismatch("precursor IDs must match rows"))
    valid(x) = !ismissing(x) && isfinite(x)
    rows = sortperm(precursor_ids)
    filter!(i -> any(valid, @view X[i, :]), rows)
    if length(rows) > DIRECTLFQ_MAX_PRECURSORS
        sort!(rows; by=i -> (-count(valid, @view X[i, :]),
            -sum(x -> valid(x) ? Float64(x) : 0.0, @view X[i, :])))
        resize!(rows, DIRECTLFQ_MAX_PRECURSORS)
    end
    traces = [valid(X[i, r]) ? Float64(X[i, r]) : NaN for i in rows, r in 1:nruns]
    estimates, contributions = directlfq_profile!(traces)
    return estimates, rows, contributions
end

function directlfq_profile!(traces::Matrix{Float64})
    n, nruns = size(traces)
    estimates = Vector{Union{Missing, Float32}}(missing, nruns)
    contributions = zeros(UInt32, nruns)
    isempty(traces) && return estimates, contributions
    observed_max = maximum(x -> isfinite(x) ? x : -Inf, traces)
    observed_sum = sum(x -> isfinite(x) ? exp2(x - observed_max) : 0.0, traces)
    buffer = Float64[]
    sizehint!(buffer, max(nruns, n))
    if nruns > 1
        if n <= DIRECTLFQ_REFERENCE_PRECURSORS
            directlfq_align_reference!(traces, buffer)
        else
            order = sortperm([count(isfinite, @view traces[i, :]) for i in axes(traces, 1)]; rev=true)
            anchors = order[1:DIRECTLFQ_REFERENCE_PRECURSORS]
            reference_traces = traces[anchors, :]
            directlfq_align_reference!(reference_traces, buffer)
            traces[anchors, :] = reference_traces
            reference = fill(NaN, nruns)
            for r in 1:nruns
                empty!(buffer)
                for i in axes(reference_traces, 1)
                    v = reference_traces[i, r]
                    isfinite(v) && push!(buffer, v)
                end
                isempty(buffer) || (reference[r] = median!(buffer))
            end
            for i in order[DIRECTLFQ_REFERENCE_PRECURSORS+1:end]
                empty!(buffer)
                for r in 1:nruns
                    a, b = reference[r], traces[i, r]
                    isfinite(a) && isfinite(b) && push!(buffer, a - b)
                end
                shift = isempty(buffer) ? NaN : median!(buffer)
                traces[i, :] .+= shift
            end
        end
    end
    profile = fill(NaN, nruns)
    for r in 1:nruns
        empty!(buffer)
        for i in axes(traces, 1)
            v = traces[i, r]
            isfinite(v) && push!(buffer, v)
        end
        contributions[r] = length(buffer)
        isempty(buffer) || (profile[r] = median!(buffer))
    end
    profile_max = maximum(x -> isfinite(x) ? x : -Inf, profile; init=-Inf)
    if isfinite(profile_max)
        profile_sum = sum(x -> isfinite(x) ? exp2(x - profile_max) : 0.0, profile)
        scale = observed_max + log2(observed_sum) - profile_max - log2(profile_sum)
        for r in 1:nruns
            isfinite(profile[r]) && (estimates[r] = Float32(profile[r] + scale))
        end
    end
    return estimates, contributions
end

"""
    directlfq_from_observations(peptides, experiments, abundance, experiments_dict)

Select at most 100 precursor traces before allocating the precursor-by-run
matrix. Input must contain at most one observation per precursor and run.
Returns log2 protein abundances and contributing precursor counts.
"""
function directlfq_from_observations(peptides, experiments, abundance, experiments_dict)
    counts = Dict{UInt32, Int}()
    log_sums = Dict{UInt32, Float64}()
    for i in eachindex(peptides)
        value = abundance[i]
        if !ismissing(value) && isfinite(value) && value > 0
            id = peptides[i]
            counts[id] = get(counts, id, 0) + 1
            log_sums[id] = get(log_sums, id, 0.0) + log2(Float64(value))
        end
    end
    selected = sort!(collect(keys(counts)))
    if length(selected) > DIRECTLFQ_MAX_PRECURSORS
        sort!(selected; by=id -> (-counts[id], -log_sums[id]))
        resize!(selected, DIRECTLFQ_MAX_PRECURSORS)
    end
    indices = Dict(id => i for (i, id) in enumerate(selected))
    traces = fill(NaN, length(selected), length(experiments_dict))
    for i in eachindex(peptides)
        row = get(indices, peptides[i], 0)
        row == 0 && continue
        value = abundance[i]
        if !ismissing(value) && isfinite(value) && value > 0
            col = experiments_dict[experiments[i]]
            isfinite(traces[row, col]) && throw(ArgumentError(
                "Duplicate precursor/run observation in directLFQ input"))
            traces[row, col] = log2(Float64(value))
        end
    end
    return directlfq_profile!(traces)
end
