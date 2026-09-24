# Shared, bounded score ordering and calibration. Temporary records are internal
# to one process; no persistent binary format is exposed.
const SCORE_WORKSPACE_BYTES = 256 * 1024^2
struct ScoreGroup
    score::Float64
    targets::Int64
    decoys::Int64
end
_score_key(score) = isnan(score) ? -Inf : (iszero(score) ? 0.0 : Float64(score))
_group_fdr(t, d, scale) = t == 0 ? Inf32 : Float32(Float64(d) * scale / t)
function _validate_score_options(scale, budget)
    isfinite(scale) && scale > 0 || throw(ArgumentError("FDR scale factor must be positive and finite"))
    budget >= 4096 || throw(ArgumentError("Score workspace must be at least 4096 bytes"))
end
function _read_score_record(io, ::Type{T}) where T
    record = Ref{T}()
    read!(io, record)
    return record[]
end
function _emit_score_arrays(emit, scores, labels)
    length(scores) == length(labels) || throw(DimensionMismatch("Scores and labels differ in length"))
    for i in 1:length(scores)
        emit(scores[i], labels[i])
    end
end

function _visit_sorted_groups(visit, scores, labels)
    order = sortperm(scores; rev=true)
    i = 1
    while i <= length(order)
        score = scores[order[i]]
        t = d = 0
        j = i
        while j <= length(order) && scores[order[j]] == score
            labels[order[j]] ? (t += 1) : (d += 1)
            j += 1
        end
        visit(ScoreGroup(score, t, d))
        i = j
    end
end

function _merge_score_groups(visit, paths)
    streams = IOStream[]
    heap = BinaryMinHeap{Tuple{Float64, Int, Int64, Int64}}()
    try
        for path in paths
            io = open(path, "r")
            push!(streams, io)
            if !eof(io)
                g = _read_score_record(io, ScoreGroup)
                push!(heap, (-g.score, length(streams), g.targets, g.decoys))
            end
        end
        pending = nothing
        processed = 0
        started = last_progress = time()
        while !isempty(heap)
            negative_score, source, t, d = pop!(heap)
            processed += 1
            if processed % 1_000_000 == 0 && time() - last_progress >= 60
                @debug_l1 "Score group merge: inputs=$(length(paths)) records=$processed elapsed=$(round(time()-started, digits=2))s"
                last_progress = time()
            end
            score = -negative_score
            if pending !== nothing && pending.score == score
                pending = ScoreGroup(score, pending.targets + t, pending.decoys + d)
            else
                pending !== nothing && visit(pending)
                pending = ScoreGroup(score, t, d)
            end
            io = streams[source]
            if !eof(io)
                g = _read_score_record(io, ScoreGroup)
                push!(heap, (-g.score, source, g.targets, g.decoys))
            end
        end
        pending !== nothing && visit(pending)
    finally
        foreach(close, streams)
    end
end

"""
    with_sorted_score_groups(consume, produce; memory_budget_bytes, max_fanin)

`produce(emit)` supplies score/target observations. `consume(groups)` receives a
callback iterator: `groups(visit)` emits descending, unique ScoreGroups. Sorting
spills compact grouped chunks when needed, with bounded fan-in and workspace.
NaN scores are treated as -Inf (worst score); signed zero is one score group.
The budget excludes caller-owned input/output and Arrow's input-file buffers.
"""
function with_sorted_score_groups(consume, produce;
    memory_budget_bytes::Int=SCORE_WORKSPACE_BYTES, max_fanin::Int=32,
    temp_parent::AbstractString=tempdir(),
)
    memory_budget_bytes >= 4096 || throw(ArgumentError("Score workspace must be at least 4096 bytes"))
    max_fanin >= 2 || throw(ArgumentError("Score merge fan-in must be at least two"))
    capacity = max(1, memory_budget_bytes ÷ 64)
    mktempdir(temp_parent) do directory
        scores = Float64[]
        labels = Bool[]
        # Level compaction bounds the number of retained paths as well as open files.
        levels = Vector{Vector{String}}()
        counter = 0
        next_path() = joinpath(directory, "$(counter += 1).bin")
        function stage(path, level=1)
            while length(levels) < level
                push!(levels, String[])
            end
            push!(levels[level], path)
            if length(levels[level]) == max_fanin
                merged = next_path()
                open(merged, "w") do io
                    _merge_score_groups(g -> write(io, Ref(g)), levels[level])
                end
                foreach(rm, levels[level])
                empty!(levels[level])
                stage(merged, level + 1)
            end
        end
        function flush_chunk()
            path = next_path()
            open(path, "w") do io
                _visit_sorted_groups(g -> write(io, Ref(g)), scores, labels)
            end
            empty!(scores); empty!(labels)
            stage(path)
        end
        rows = 0
        started = last_log = time()
        produce() do score, target
            push!(scores, _score_key(score)); push!(labels, Bool(target))
            rows += 1
            if length(scores) == capacity
                flush_chunk()
                if time() - last_log >= 60
                    @debug_l1 "Score grouping: rows=$rows chunks=$counter elapsed=$(round(time()-started, digits=2))s"
                    last_log = time()
                end
            end
        end
        if isempty(levels)
            @debug_l1 "Score grouping in memory: rows=$rows budget_bytes=$memory_budget_bytes"
            return consume(visit -> _visit_sorted_groups(visit, scores, labels))
        end
        !isempty(scores) && flush_chunk()
        paths = reduce(vcat, levels; init=String[])
        while length(paths) > max_fanin
            next = String[]
            for first in 1:max_fanin:length(paths)
                batch = paths[first:min(first + max_fanin - 1, end)]
                path = next_path()
                open(path, "w") do io
                    _merge_score_groups(g -> write(io, Ref(g)), batch)
                end
                foreach(rm, batch)
                push!(next, path)
            end
            paths = next
        end
        @debug_l1 "Score grouping final merge: rows=$rows inputs=$(length(paths)) budget_bytes=$memory_budget_bytes"
        return consume(visit -> _merge_score_groups(visit, paths))
    end
end

"""Find exact score and target cutoffs for several q-value thresholds in one
forward scan. A `nothing` cutoff means no qualifying observation.
"""
function qvalue_score_cutoffs(produce, thresholds; fdr_scale_factor=1.0f0,
    memory_budget_bytes=SCORE_WORKSPACE_BYTES, max_fanin=32, temp_parent=tempdir())
    _validate_score_options(fdr_scale_factor, memory_budget_bytes)
    with_sorted_score_groups(produce; memory_budget_bytes, max_fanin, temp_parent) do groups
        targets = decoys = 0
        last_target = nothing
        floors = Union{Nothing,Float64}[nothing for _ in thresholds]
        target_floors = copy(floors)
        groups() do g
            targets += g.targets; decoys += g.decoys
            g.targets > 0 && (last_target = g.score)
            q = _group_fdr(targets, decoys, fdr_scale_factor)
            for i in eachindex(thresholds)
                if q <= thresholds[i]
                    floors[i] = g.score
                    target_floors[i] = last_target
                end
            end
        end
        return (; floors, target_floors)
    end
end

"""Exact lowest passing target score, or `nothing` when no target passes."""
function qvalue_score_cutoff(produce; q_threshold=0.01f0, kwargs...)
    return only(qvalue_score_cutoffs(produce, (q_threshold,); kwargs...).target_floors)
end

struct ScoreFitRecord
    score::Float64
    qvalue::Float64
    pep::Float64
end
struct ScorePAVABlock
    decoys::Float64
    weight::Float64
    groups::Int64
end

# A disk-backed PAVA stack keeps only its top block buffer resident.
mutable struct ScorePAVAStack
    io::IOStream
    disk_blocks::Int
    buffer::Vector{ScorePAVABlock}
    capacity::Int
end
function _stack_top(s::ScorePAVAStack)
    if isempty(s.buffer) && s.disk_blocks > 0
        n = min(s.capacity, s.disk_blocks)
        s.disk_blocks -= n
        seek(s.io, s.disk_blocks * sizeof(ScorePAVABlock))
        resize!(s.buffer, n)
        read!(s.io, s.buffer)
    end
    return isempty(s.buffer) ? nothing : last(s.buffer)
end
function _stack_push!(s::ScorePAVAStack, b::ScorePAVABlock)
    while true
        top = _stack_top(s)
        (top === nothing || top.decoys / top.weight <= b.decoys / b.weight) && break
        pop!(s.buffer)
        b = ScorePAVABlock(top.decoys + b.decoys, top.weight + b.weight, top.groups + b.groups)
    end
    push!(s.buffer, b)
    if length(s.buffer) >= 2 * s.capacity
        seek(s.io, s.disk_blocks * sizeof(ScorePAVABlock))
        write(s.io, view(s.buffer, 1:s.capacity))
        s.disk_blocks += s.capacity
        deleteat!(s.buffer, 1:s.capacity)
    end
end

mutable struct ScoreCalibrationStore
    path::String
    io::IOStream
    n::Int
    cache::Dict{Int, Vector{ScoreFitRecord}}
    lock::ReentrantLock
    page_size::Int
    max_pages::Int
end
struct ScoreCalibration
    store::ScoreCalibrationStore
    column::Symbol
end
function _fit_record(store, index)
    page = (index - 1) ÷ store.page_size
    records = get(store.cache, page, nothing)
    if records === nothing
        length(store.cache) >= store.max_pages && empty!(store.cache)
        start = page * store.page_size
        records = Vector{ScoreFitRecord}(undef, min(store.page_size, store.n - start))
        seek(store.io, start * sizeof(ScoreFitRecord)); read!(store.io, records)
        store.cache[page] = records
    end
    return records[(index - 1) % store.page_size + 1]
end
function (cal::ScoreCalibration)(score::Real)
    store = cal.store
    x = _score_key(score)
    lock(store.lock) do
        lo, hi = 1, store.n
        while lo <= hi
            mid = (lo + hi) >>> 1
            record = _fit_record(store, mid)
            record.score == x && return getfield(record, cal.column)
            if record.score > x
                lo = mid + 1
            else
                hi = mid - 1
            end
        end
        hi < 1 && return getfield(_fit_record(store, 1), cal.column)
        lo > store.n && return getfield(_fit_record(store, store.n), cal.column)
        upper, lower = _fit_record(store, hi), _fit_record(store, lo)
        a, b = getfield(upper, cal.column), getfield(lower, cal.column)
        a == b && return a
        (!isfinite(a) || !isfinite(b) || !isfinite(lower.score)) && return b
        return b + (a-b) * ((x-lower.score)/(upper.score-lower.score))
    end
end

"""Build exact grouped q-values and weighted PAVA in bounded disk-backed workspace.
The returned callable mappings share a bounded page cache; their temporary file
is removed when both mappings become unreachable. Between observed scores the
mapping interpolates linearly and uses flat extrapolation outside the range.
"""
function build_score_calibration(produce; compute_pep=true, fdr_scale_factor=1.0f0,
    memory_budget_bytes=SCORE_WORKSPACE_BYTES, max_fanin=32, temp_parent=tempdir())
    _validate_score_options(fdr_scale_factor, memory_budget_bytes)
    directory = mktempdir(temp_parent)
    path = joinpath(directory, "calibration.bin")
    n = 0
    try
        open(path, "w+") do output
            open(joinpath(directory, "pava.bin"), "w+") do stack_io
                stack = ScorePAVAStack(stack_io, 0, ScorePAVABlock[], max(1, min(1024, memory_budget_bytes ÷ 256)))
                compute_pep && _stack_push!(stack, ScorePAVABlock(0.0, 1.0, 0))
                with_sorted_score_groups(produce; memory_budget_bytes, max_fanin, temp_parent) do groups
                    t = d = 0
                    groups() do g
                        t += g.targets; d += g.decoys; n += 1
                        write(output, Ref(ScoreFitRecord(g.score, _group_fdr(t,d,fdr_scale_factor), NaN)))
                        if compute_pep
                            weighted_d = Float64(g.decoys) * fdr_scale_factor
                            _stack_push!(stack, ScorePAVABlock(weighted_d, g.targets + weighted_d, 1))
                        end
                    end
                end
                @debug_l1 "Score calibration grouped counts complete: groups=$n"
                if compute_pep
                    @debug_l1 "Score calibration PAVA expansion starting: groups=$n"
                    # Flush the remaining stack, then assign fitted blocks sequentially.
                    seek(stack_io, stack.disk_blocks * sizeof(ScorePAVABlock))
                    write(stack_io, stack.buffer)
                    blocks = stack.disk_blocks + length(stack.buffer)
                    flush(stack_io); seekstart(stack_io); seekstart(output)
                    fitted_path = joinpath(directory, "fitted.bin")
                    open(fitted_path, "w") do fitted
                        for _ in 1:blocks
                            b = _read_score_record(stack_io, ScorePAVABlock)
                            pep = clamp(b.decoys / (b.weight - b.decoys), 0.0, 1.0)
                            for _ in 1:b.groups
                                record = _read_score_record(output, ScoreFitRecord)
                                write(fitted, Ref(ScoreFitRecord(record.score, record.qvalue, pep)))
                            end
                        end
                    end
                    seekstart(output)
                    open(fitted_path, "r") do fitted
                        buffer = Vector{UInt8}(undef, max(24, min(65536, memory_budget_bytes ÷ 8)))
                        while !eof(fitted)
                            count = readbytes!(fitted, buffer)
                            write(output, view(buffer, 1:count))
                        end
                    end
                    rm(fitted_path)
                end
                # Reverse by bounded blocks, without reversing the full permutation.
                @debug_l1 "Score calibration reverse q-value scan starting: groups=$n"
                minimum_q = Inf
                block_size = max(1, min(4096, memory_budget_bytes ÷ 192))
                for stop in n:-block_size:1
                    start = max(1, stop - block_size + 1)
                    records = Vector{ScoreFitRecord}(undef, stop-start+1)
                    seek(output, (start-1)*sizeof(ScoreFitRecord)); read!(output, records)
                    for i in length(records):-1:1
                        r = records[i]; minimum_q = min(minimum_q, r.qvalue)
                        records[i] = ScoreFitRecord(r.score, minimum_q, r.pep)
                    end
                    seek(output, (start-1)*sizeof(ScoreFitRecord)); write(output, records)
                end
            end
        end
        rm(joinpath(directory, "pava.bin"))
        if n == 0
            rm(directory; recursive=true)
            return nothing
        end
        @debug_l1 "Score calibration compaction starting: groups=$n"
        compact_path = joinpath(directory, "compact.bin")
        n = _compact_score_calibration(path, compact_path)
        mv(compact_path, path; force=true)
        @debug_l1 "Score calibration complete: knots=$n"
        page_size = max(1, min(1024, memory_budget_bytes ÷ 192))
        store = ScoreCalibrationStore(path, open(path,"r"), n,
            Dict{Int,Vector{ScoreFitRecord}}(), ReentrantLock(), page_size,
            max(1, memory_budget_bytes ÷ (4 * page_size * sizeof(ScoreFitRecord))))
        finalizer(store) do s
            close(s.io)
            rm(dirname(s.path); recursive=true, force=true)
        end
        return (qval_spline=ScoreCalibration(store, :qvalue),
                pep_interp=compute_pep ? ScoreCalibration(store, :pep) : nothing)
    catch
        rm(directory; recursive=true, force=true)
        rethrow()
    end
end

# Retain both endpoints of every constant span: interpolation remains exact at
# observed scores, while long q-value/PEP plateaus require only two records.
function _compact_score_calibration(path, destination)
    n = 0
    open(path, "r") do source
        open(destination, "w") do output
            first = last = _read_score_record(source, ScoreFitRecord)
            write(output, Ref(first)); n += 1
            while !eof(source)
                current = _read_score_record(source, ScoreFitRecord)
                same = isequal(current.qvalue, last.qvalue) && isequal(current.pep, last.pep)
                if !same
                    if last.score != first.score
                        write(output, Ref(last)); n += 1
                    end
                    write(output, Ref(current)); n += 1
                    first = current
                end
                last = current
            end
            if last.score != first.score
                write(output, Ref(last)); n += 1
            end
        end
    end
    return n
end

function _close_score_calibration(fit)
    store = fit.qval_spline.store
    close(store.io)
    empty!(store.cache)
    rm(dirname(store.path); recursive=true, force=true)
end

# Sequentially cached permutation for bulk calibration. The caller owns scores
# and output arrays; only the sort/merge buffers and PAVA stack use workspace.
mutable struct DiskScoreOrder <: AbstractVector{Int64}
    io::IOStream
    n::Int
    buffer::Vector{Int64}
    first::Int
    capacity::Int
end
Base.size(order::DiskScoreOrder) = (order.n,)
Base.IndexStyle(::Type{DiskScoreOrder}) = IndexLinear()
@inline function Base.getindex(order::DiskScoreOrder, i::Int)
    @boundscheck checkbounds(order, i)
    if !(order.first <= i < order.first + length(order.buffer))
        order.first = ((i - 1) ÷ order.capacity) * order.capacity + 1
        resize!(order.buffer, min(order.capacity, order.n - order.first + 1))
        seek(order.io, (order.first - 1) * sizeof(Int64))
        read!(order.io, order.buffer)
    end
    @inbounds return order.buffer[i - order.first + 1]
end

function _merge_score_orders!(destination, paths, scores, budget)
    orders = DiskScoreOrder[]
    heap = BinaryMinHeap{Tuple{Float64, Int64, Int}}()
    positions = ones(Int, length(paths))
    capacity = max(1, min(8192, budget ÷ (32 * (length(paths) + 1))))
    output = Int64[]
    sizehint!(output, capacity)
    try
        for path in paths
            order = DiskScoreOrder(open(path, "r"), filesize(path) ÷ 8, Int64[], 0, capacity)
            push!(orders, order)
            if !isempty(order)
                row = order[1]
                push!(heap, (-_score_key(scores[row]), row, length(orders)))
            end
        end
        open(destination, "w") do io
            while !isempty(heap)
                _, row, source = pop!(heap)
                push!(output, row)
                if length(output) == capacity
                    write(io, output)
                    empty!(output)
                end
                positions[source] += 1
                order = orders[source]
                if positions[source] <= length(order)
                    next_row = order[positions[source]]
                    push!(heap, (-_score_key(scores[next_row]), next_row, source))
                end
            end
            write(io, output)
        end
    finally
        foreach(order -> close(order.io), orders)
    end
    return nothing
end

"""Visit a descending score permutation, spilling indexed runs with bounded fan-in."""
function with_score_order(consume, scores; memory_budget_bytes=SCORE_WORKSPACE_BYTES,
    max_fanin=32, temp_parent=tempdir())
    _validate_score_options(1.0, memory_budget_bytes)
    max_fanin >= 2 || throw(ArgumentError("Score merge fan-in must be at least two"))
    capacity = max(1, memory_budget_bytes ÷ 64)
    length(scores) <= capacity && return consume(_score_order(scores))
    mktempdir(temp_parent) do directory
        levels = Vector{Vector{String}}()
        counter = 0
        next_path() = joinpath(directory, "$(counter += 1).bin")
        function stage(path, level=1)
            while length(levels) < level
                push!(levels, String[])
            end
            push!(levels[level], path)
            if length(levels[level]) == max_fanin
                merged = next_path()
                _merge_score_orders!(merged, levels[level], scores, memory_budget_bytes)
                foreach(rm, levels[level])
                empty!(levels[level])
                stage(merged, level + 1)
            end
        end
        for first in 1:capacity:length(scores)
            rows = collect(Int64, first:min(first + capacity - 1, length(scores)))
            sort!(rows; by=i -> _score_key(scores[i]), rev=true)
            path = next_path()
            open(io -> write(io, rows), path, "w")
            stage(path)
        end
        paths = reduce(vcat, levels; init=String[])
        while length(paths) > 1
            next = String[]
            for first in 1:max_fanin:length(paths)
                batch = paths[first:min(first + max_fanin - 1, end)]
                path = next_path()
                _merge_score_orders!(path, batch, scores, memory_budget_bytes)
                foreach(rm, batch)
                push!(next, path)
            end
            paths = next
        end
        open(only(paths), "r") do io
            order = DiskScoreOrder(io, length(scores), Int64[], 0,
                max(1, min(8192, memory_budget_bytes ÷ 64)))
            return consume(order)
        end
    end
end
