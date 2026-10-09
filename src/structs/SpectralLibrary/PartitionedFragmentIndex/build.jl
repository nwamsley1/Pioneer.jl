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

"Initial partition (1-based) of a precursor m/z: `partition_width` bins from the smallest precursor m/z."
@inline _initial_partition(pmz, min_prec_mz, partition_width, n_initial) =
    clamp(floor(Int, (pmz - min_prec_mz) / partition_width) + 1, 1, n_initial)

"""
    resolve_local_id_type(requested, prec_mzs, partition_width) -> (id_type, n_over, n_bins)

Local precursor ID type of the fragment index. `requested` is `"UInt16"`, `"UInt32"` or `"auto"`. With `"auto"`,
UInt32 is chosen only when some `partition_width` bin of precursor m/z holds more precursors than a UInt16 partition
can (`MAX_LOCAL_PRECS`), i.e. when the ID width rather than `partition_width` would set the partition layout (UInt16
would split those bins, e.g. 5 Da into ~2.5 Da on a 10 M-precursor library). Returns the type, the number of bins
over the limit, and the number of bins.
"""
function resolve_local_id_type(requested::AbstractString, prec_mzs::AbstractVector{<:Real}, partition_width::Real)
    requested in ("auto", "UInt16", "UInt32") ||
        throw(ArgumentError("frag_index_local_id_type must be \"auto\", \"UInt16\" or \"UInt32\" (got \"$requested\")"))
    isempty(prec_mzs) && return (requested == "UInt32" ? UInt32 : UInt16), 0, 0
    lo, hi = extrema(prec_mzs)
    n = max(1, ceil(Int, (hi - lo) / partition_width))
    counts = zeros(Int, n)
    for pmz in prec_mzs
        counts[_initial_partition(pmz, lo, partition_width, n)] += 1
    end
    n_over = count(>(MAX_LOCAL_PRECS), counts)
    id_type = requested == "UInt16" ? UInt16 : requested == "UInt32" ? UInt32 : (n_over > 0 ? UInt32 : UInt16)
    return id_type, n_over, n
end

"""
    IndexFragSelection

Each precursor's fragment-index fragments: its first `max_rank` (8) fragments that pass the index ion-type filters,
in the library's rank order, as `SimpleFrag`s carrying the GLOBAL precursor id and the rank bitmask
`1 << (rank-1)`. Independent of the partition width and of the RT binning, so one selection serves every index a
library is built with (`build_partitioned_index_from_selection`). Precursor `pid`'s fragments are
`frags[starts[pid] : starts[pid] + counts[pid] - 1]` (`index_frag_range`); `frags` may be memory-mapped, and a
precursor's fragments need not follow the previous precursor's (a piece's can be gathered into their own file).
"""
struct IndexFragSelection
    frags::Vector{SimpleFrag{Float32}}
    starts::Vector{Int}
    counts::Vector{UInt8}
    prec_mzs::Vector{Float32}
end

"A selection whose precursors' fragments are consecutive: precursor `pid`'s are `frags[offsets[pid]:offsets[pid+1]-1]`."
function IndexFragSelection(frags::Vector{SimpleFrag{Float32}}, offsets::Vector{Int}, prec_mzs::Vector{Float32})
    n = length(offsets) - 1
    return IndexFragSelection(frags, offsets[1:n], UInt8[offsets[k + 1] - offsets[k] for k in 1:n], prec_mzs)
end

@inline index_frag_range(sel::IndexFragSelection, pid::Integer) =
    sel.starts[pid]:(sel.starts[pid] + Int(sel.counts[pid]) - 1)
@inline n_index_frags(sel::IndexFragSelection, pid::Integer) = Int(sel.counts[pid])

"""
    _visit_index_frags(f, frag_lookup, detailed_frags, pid, (y_start_index, b_start_index, include_p_index)) -> Int

Calls `f(rank, dfrag)` for precursor `pid`'s index fragments in rank order: its fragments that pass the index
ion-type filters, at most 8 (the UInt8 rank bitmask). Returns how many it visited.
"""
@inline _visit_index_frags(f, frag_lookup, detailed_frags, pid, filt::Tuple{UInt8, UInt8, Bool}) =
    _visit_index_frags(f, detailed_frags, getPrecFragRange(frag_lookup, pid), filt)

"The same over the fragments `detailed_frags[frag_range]` (one precursor's, in rank order)."
@inline function _visit_index_frags(f, detailed_frags, frag_range::AbstractUnitRange, filt::Tuple{UInt8, UInt8, Bool})
    y_start_index, b_start_index, include_p_index = filt
    rank = 0
    for fi in frag_range
        dfrag = detailed_frags[fi]
        # Apply fragment index ion-type filters
        if isY(dfrag)
            getIonPosition(dfrag) < y_start_index && continue
        elseif isB(dfrag)
            getIonPosition(dfrag) < b_start_index && continue
        elseif isP(dfrag)
            include_p_index || continue
        end
        isIso(dfrag) && continue
        rank == 8 && break
        rank += 1
        f(rank, dfrag)
    end
    return rank
end

"""
    select_index_fragments(spec_lib; y_start_index=UInt8(4), b_start_index=UInt8(3), include_p_index=false)
        -> IndexFragSelection

The fragments every partitioned index of `spec_lib` is built from (see `IndexFragSelection`). Threaded over
precursors: a counting pass sizes each precursor's slot, a second pass fills it.
"""
function select_index_fragments(
    spec_lib::SpectralLibrary;
    y_start_index::UInt8 = UInt8(4),
    b_start_index::UInt8 = UInt8(3),
    include_p_index::Bool = false,
)
    precursors = getPrecursors(spec_lib)
    frag_lookup = getFragmentLookupTable(spec_lib)
    detailed_frags = getFragments(frag_lookup)
    prec_mzs = Vector{Float32}(getMz(precursors))
    prec_irts = getIrt(precursors)
    n_precursors = length(prec_mzs)
    filt = (y_start_index, b_start_index, include_p_index)

    counts = zeros(Int, n_precursors)
    Threads.@threads :static for pid in UInt32(1):UInt32(n_precursors)
        counts[pid] = _visit_index_frags((_, _) -> nothing, frag_lookup, detailed_frags, pid, filt)
    end
    offsets = Vector{Int}(undef, n_precursors + 1)
    offsets[1] = 1
    for pid in 1:n_precursors
        offsets[pid + 1] = offsets[pid] + counts[pid]
    end
    frags = Vector{SimpleFrag{Float32}}(undef, offsets[end] - 1)
    Threads.@threads :static for pid in UInt32(1):UInt32(n_precursors)
        pmz = prec_mzs[pid]; pirt = prec_irts[pid]; o = offsets[pid] - 1
        _visit_index_frags(frag_lookup, detailed_frags, pid, filt) do rank, dfrag
            frags[o + rank] = SimpleFrag{Float32}(
                getMz(dfrag),
                pid,
                pmz,
                pirt,
                UInt8(0),
                UInt8(1) << UInt8(rank - 1),  # bitmask: rank 1→bit0, rank 2→bit1, ...
            )
        end
    end
    return IndexFragSelection(frags, offsets, prec_mzs)
end

"""
    build_partitioned_index_from_lib(spec_lib; partition_width=5.0f0,
        frag_bin_tol_ppm=2.5f0, rt_bin_tol=3.0f0,
        y_start_index=UInt8(4), b_start_index=UInt8(3),
        include_p_index=false, id_type=UInt16)

Build a LocalPartitionedFragmentIndex from scratch using the spectral library's
DetailedFrag data. Each partition gets its own independently-constructed index
with LocalFragment entries (UInt16 local_id + UInt8 score = 4 bytes).

`id_type = UInt32` builds a `LocalPartitionedFragmentIndex32` instead (`LocalFragment32`, 8 bytes): partitions are
then never split, so they keep the nominal `partition_width`.

Per-fragment score is a bitmask `1 << (rank-1)`, capped at 8 ranks (UInt8).

Precursor IDs are remapped to partition-local values (1..N; N ≤ 65535 for UInt16).
Partitions that would exceed that many unique precursors are automatically split.

One index: `select_index_fragments`, then `build_partitioned_index_from_selection`. To build several indexes of one
library (widths, main and presearch), select once and build each from the selection.
"""
function build_partitioned_index_from_lib(
    spec_lib::SpectralLibrary;
    partition_width::Float32 = 5.0f0,
    frag_bin_tol_ppm::Float32 = 0.0f0,
    frag_bin_tol_mda::Float32 = 2.0f0,
    rt_bin_tol::Float32 = 3.0f0,
    y_start_index::UInt8 = UInt8(4),
    b_start_index::UInt8 = UInt8(3),
    include_p_index::Bool = false,
    id_type::Type{<:Unsigned} = UInt16,
)
    sel = select_index_fragments(spec_lib; y_start_index = y_start_index, b_start_index = b_start_index,
                                 include_p_index = include_p_index)
    return build_partitioned_index_from_selection(sel; partition_width = partition_width,
        frag_bin_tol_ppm = frag_bin_tol_ppm, frag_bin_tol_mda = frag_bin_tol_mda, rt_bin_tol = rt_bin_tol,
        id_type = id_type)
end

"""
    initial_partitions(prec_mzs, partition_width) -> Vector{Vector{UInt32}}

Global precursor IDs per initial `partition_width` bin of precursor m/z, from the smallest precursor m/z (Step 1 of
`build_partitioned_index_from_selection`). Building from a consecutive run of these bins (`build_index_pieces`)
gives exactly the partitions a single index over all bins has.
"""
function initial_partitions(prec_mzs::AbstractVector{<:Real}, partition_width::Real)
    isempty(prec_mzs) && return [UInt32[]]
    min_prec_mz, max_prec_mz = Float32.(extrema(prec_mzs))
    width = Float32(partition_width)
    n_initial = max(1, ceil(Int, (max_prec_mz - min_prec_mz) / width))
    @debug_l2 "build_partitioned_index: prec m/z [$(round(min_prec_mz, digits=2)), $(round(max_prec_mz, digits=2))], $(width) Da → $(n_initial) initial partitions"
    pids = [UInt32[] for _ in 1:n_initial]
    for pid in UInt32(1):UInt32(length(prec_mzs))
        push!(pids[_initial_partition(Float32(prec_mzs[pid]), min_prec_mz, width, n_initial)], pid)
    end
    return pids
end

"""
    build_partitioned_index_from_selection(sel; partition_width, frag_bin_tol_ppm, frag_bin_tol_mda, rt_bin_tol, id_type)

The partitioned index of one width from a precomputed `IndexFragSelection`; see `build_partitioned_index_from_lib`.
Partitions are built in parallel. Each partition's fragments are laid out exactly as the serial builder did
(precursors in partition order, each precursor's fragments in rank order), so the index is identical.
"""
function build_partitioned_index_from_selection(
    sel::IndexFragSelection;
    partition_width::Float32 = 5.0f0,
    frag_bin_tol_ppm::Float32 = 0.0f0,
    frag_bin_tol_mda::Float32 = 2.0f0,
    rt_bin_tol::Float32 = 3.0f0,
    id_type::Type{<:Unsigned} = UInt16,
    initial_partition_pids::Vector{Vector{UInt32}} = initial_partitions(sel.prec_mzs, partition_width),
)
    max_local = max_local_precs(id_type)
    prec_mzs = sel.prec_mzs

    # ── Step 2: Split partitions exceeding max_local (balanced halving) ───────
    final_partition_pids = Vector{UInt32}[]
    function _split_balanced!(out::Vector{Vector{UInt32}}, pids::Vector{UInt32}, prec_mzs)
        if length(pids) <= max_local
            push!(out, pids)
        else
            sort!(pids, by = pid -> prec_mzs[pid])
            mid = length(pids) ÷ 2
            _split_balanced!(out, pids[1:mid], prec_mzs)
            _split_balanced!(out, pids[mid+1:end], prec_mzs)
        end
    end
    for pids in initial_partition_pids
        _split_balanced!(final_partition_pids, pids, prec_mzs)
    end
    n_partitions = length(final_partition_pids)
    n_initial = length(initial_partition_pids)
    if n_partitions != n_initial
        @debug_l2 "build_partitioned_index: split to $(n_partitions) partitions ($(n_partitions - n_initial) extra from UInt16 limit)"
    end

    # ── Steps 3-4, per partition in parallel: the partition's SimpleFrags with local IDs (precursors in
    # partition order, each precursor's fragments in rank order), then its LocalPartition ──────────────
    PT = local_partition_type(id_type){Float32}
    partitions = Vector{PT}(undef, n_partitions)
    Threads.@threads :dynamic for k in 1:n_partitions
        pids = final_partition_pids[k]
        l2g = Vector{UInt32}(pids)                   # local_id i → global prec_id
        n_local = id_type(length(l2g))
        frags_k = Vector{SimpleFrag{Float32}}(undef, sum(pid -> n_index_frags(sel, pid), pids; init = 0))
        j = 0
        for (i, pid) in enumerate(pids)
            for fi in index_frag_range(sel, pid)
                f = sel.frags[fi]
                frags_k[j += 1] = SimpleFrag{Float32}(f.mz, UInt32(i), f.prec_mz, f.prec_irt, f.prec_charge, f.score)
            end
        end

        if isempty(frags_k)
            partitions[k] = PT(
                SoAFragBins{Float32}(Float32[], Float32[], UInt32[], UInt32[]),
                FragIndexBin{Float32}[],
                local_fragment_type(id_type)[],
                l2g,
                n_local,
                UInt16[],
            )
            continue
        end

        partitions[k] = _build_local_partition(frags_k, l2g, n_local,
                                                frag_bin_tol_ppm, frag_bin_tol_mda, rt_bin_tol)
    end
    @debug_l2 "build_partitioned_index: $(length(sel.frags)) total fragments across $(n_partitions) partitions"

    # Compute per-partition prec_mz bounds
    partition_bounds = Vector{Tuple{Float32, Float32}}(undef, n_partitions)
    for k in 1:n_partitions
        pids = final_partition_pids[k]
        if isempty(pids)
            partition_bounds[k] = (Float32(Inf), Float32(-Inf))
        else
            pmin = Float32(Inf)
            pmax = Float32(-Inf)
            for pid in pids
                pmz = prec_mzs[pid]
                pmin = min(pmin, pmz)
                pmax = max(pmax, pmz)
            end
            partition_bounds[k] = (pmin, pmax)
        end
    end

    return local_index_type(id_type){Float32}(partitions, partition_bounds, n_partitions)
end

"""
Build a LocalPartition from SimpleFrags whose prec_id field already contains
local IDs (stored as UInt32). Produces LocalFragment entries (LocalPartition32 /
LocalFragment32 when `n_local` is a UInt32).
"""
function _build_local_partition(
    frag_ions::Vector{SimpleFrag{Float32}},
    local_to_global::Vector{UInt32},
    n_local::I,
    frag_bin_tol_ppm::Float32,
    frag_bin_tol_mda::Float32,
    rt_bin_tol::Float32,
) where {I<:Unsigned}
    scratch = similar(frag_ions)   # one sort buffer for every sort of this partition (each sort! would allocate its own)
    sort!(frag_ions, by = x -> getIRT(x), scratch = scratch)

    n = length(frag_ions)
    local_fragments = Vector{local_fragment_type(I)}(undef, n)
    # bins are appended as they are found (sized per fragment, they cost 32 bytes per fragment in allocation)
    rt_bins = FragIndexBin{Float32}[]
    soa = SoAFragBins{Float32}(Float32[], Float32[], UInt32[], UInt32[])
    rt_bin_idx = 0
    frag_bin_idx = 0

    start_idx = 1
    start_irt = getIRT(frag_ions[1])

    for i in 1:n
        stop_irt = getIRT(frag_ions[i])
        if (stop_irt - start_irt > rt_bin_tol) && (i > start_idx)
            stop_idx = i - 1
            stop_irt_val = getIRT(frag_ions[stop_idx])
            sort!(@view(frag_ions[start_idx:stop_idx]), by = x -> getMZ(x), scratch = scratch)
            first_fb = frag_bin_idx + 1
            frag_bin_idx = _build_local_frag_bins!(local_fragments, soa,
                frag_bin_idx, frag_ions, start_idx, stop_idx, frag_bin_tol_ppm, frag_bin_tol_mda, scratch)
            rt_bin_idx += 1
            push!(rt_bins, FragIndexBin{Float32}(start_irt, stop_irt_val, UInt32(first_fb), UInt32(frag_bin_idx)))
            start_idx = i
            start_irt = getIRT(frag_ions[i])
        end
    end

    # Last RT bin
    stop_idx = n
    stop_irt_val = getIRT(frag_ions[stop_idx])
    sort!(@view(frag_ions[start_idx:stop_idx]), by = x -> getMZ(x), scratch = scratch)
    first_fb = frag_bin_idx + 1
    frag_bin_idx = _build_local_frag_bins!(local_fragments, soa,
        frag_bin_idx, frag_ions, start_idx, stop_idx, frag_bin_tol_ppm, frag_bin_tol_mda, scratch)
    rt_bin_idx += 1
    push!(rt_bins, FragIndexBin{Float32}(start_irt, stop_irt_val, UInt32(first_fb), UInt32(frag_bin_idx)))

    # SIMD padding on highs
    for _ in 1:7
        push!(soa.highs, Float32(Inf))  # safe sentinel for SIMD _vload8
    end
    # push! grows capacity geometrically; release the excess so a built index holds only what it uses
    for v in (soa.lows, soa.highs, soa.first_bins, soa.last_bins)
        sizehint!(v, length(v); shrink = true)
    end
    rb_final = sizehint!(rt_bins, length(rt_bins); shrink = true)
    skip_hints = _compute_skip_hints(soa, rb_final)

    return local_partition_type(I){Float32}(
        soa,
        rb_final,
        local_fragments,
        local_to_global,
        n_local,
        skip_hints,
    )
end

"""
Build fragment m/z bins producing LocalFragment / LocalFragment32 entries.
Writes bin metadata into SoA parallel arrays. Returns updated frag_bin_idx.

If `frag_bin_tol_ppm > 0`, uses ppm-width bins (legacy). Otherwise uses
fixed mDa-width bins via `frag_bin_tol_mda` (default 2.0 mDa).
"""
function _build_local_frag_bins!(
    local_fragments::Vector{F},
    soa::SoAFragBins{Float32},
    frag_bin_idx::Int,
    frag_ions::Vector{SimpleFrag{Float32}},
    start::Int, stop::Int,
    frag_bin_tol_ppm::Float32,
    frag_bin_tol_mda::Float32,
    scratch::Vector{SimpleFrag{Float32}} = similar(frag_ions, 0),
) where {F<:AbstractLocalFragment}
    I = local_id_type(F)
    use_ppm = frag_bin_tol_ppm > 0.0f0
    mda_tol = frag_bin_tol_mda * 0.001f0  # convert mDa to Da
    start_idx = start
    start_mz = getMZ(frag_ions[start])

    for i in start:stop
        stop_mz = getMZ(frag_ions[i])
        diff_mz = stop_mz - start_mz
        exceeds_tol = if use_ppm
            mean_mz = (stop_mz + start_mz) / 2
            diff_mz / (mean_mz / 1.0f6) > frag_bin_tol_ppm
        else
            diff_mz > mda_tol
        end
        if exceeds_tol && (i > start_idx)
            bin_stop = i - 1
            bin_stop_mz = getMZ(frag_ions[bin_stop])
            sort!(@view(frag_ions[start_idx:bin_stop]), by = x -> getPrecMZ(x), scratch = scratch)
            frag_bin_idx += 1
            push!(soa.lows, start_mz); push!(soa.highs, bin_stop_mz)
            push!(soa.first_bins, UInt32(start_idx)); push!(soa.last_bins, UInt32(bin_stop))
            for idx in start_idx:bin_stop
                sf = frag_ions[idx]
                local_fragments[idx] = F(
                    I(getPrecID(sf)), getScore(sf))
            end
            start_idx = i
            start_mz = getMZ(frag_ions[i])
        end
    end

    # Last frag bin
    stop_mz = getMZ(frag_ions[stop])
    sort!(@view(frag_ions[start_idx:stop]), by = x -> getPrecMZ(x), scratch = scratch)
    frag_bin_idx += 1
    push!(soa.lows, start_mz); push!(soa.highs, stop_mz)
    push!(soa.first_bins, UInt32(start_idx)); push!(soa.last_bins, UInt32(stop))
    for idx in start_idx:stop
        sf = frag_ions[idx]
        local_fragments[idx] = F(
            I(getPrecID(sf)), getScore(sf))
    end

    return frag_bin_idx
end

"""
Compute per-frag-bin skip hints for the hinted search.
hints[j] = k where frag_bins.lows[j+k] - frag_bins.lows[j] >= 5.0 Da.
"""
function _compute_skip_hints(
    frag_bins::SoAFragBins{Float32},
    rt_bins::Vector{FragIndexBin{Float32}},
)
    n_fb = length(frag_bins)
    hints = ones(UInt16, n_fb)
    lows = frag_bins.lows

    for rt_bin in rt_bins
        range = getSubBinRange(rt_bin)
        fb_start = Int(first(range))
        fb_end = Int(last(range))
        fb_start > fb_end && continue

        for j in fb_start:fb_end
            target_low = lows[j] + 5.0f0
            max_k = fb_end - j
            max_k <= 0 && continue

            if lows[fb_end] < target_low
                hints[j] = UInt16(max_k)
                continue
            end

            # Binary search for smallest k where lows[j+k] >= target_low
            lo_k = 1
            hi_k = max_k
            result_k = hi_k
            while lo_k <= hi_k
                mid_k = (lo_k + hi_k) >>> 1
                if lows[j + mid_k] >= target_low
                    result_k = mid_k
                    hi_k = mid_k - 1
                else
                    lo_k = mid_k + 1
                end
            end

            hints[j] = UInt16(clamp(result_k, 1, 65535))
        end
    end

    return hints
end
