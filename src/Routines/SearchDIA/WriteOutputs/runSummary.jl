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
Per-run QC summary (`run_summary.tsv`), in the spirit of DIA-NN's
`report.stats.tsv` / alphaDIA's `stats.tsv`: one row per raw file with
identification and quantification counts, signal, calibration and
peptide-property medians.

The per-precursor statistics are accumulated from the chunks that
`ProteinQuantificationSearch` already streams into `precursors_long.arrow`, so the summary
adds no extra pass over the precursor Arrow output. Exact order statistics use
temporary disk partitions, so retained observation data does not grow in RAM.
"""

const SUMMARY_FANOUT = 16

mutable struct RunSummaryStats
    file_name::String
    precursors_identified::Int
    precursors_quantified::Int
    precursors_mbr::Int
    peptides_identified::Int
    protein_groups_identified::Int
    protein_groups_quantified::Int
    total_peak_area::Float64
    medians::NTuple{8, Union{Missing, Float32}}
end

RunSummaryStats(name::String) = RunSummaryStats(name, 0, 0, 0, 0, 0, 0, 0.0, ntuple(_ -> missing, 8))

struct RunSummaryRecord
    run::UInt32
    peptide::UInt32
    valid::UInt32
    values::NTuple{8, Float32}
end

mutable struct RunSummaryAccumulator
    stats::Vector{RunSummaryStats}
    peptide_ids::Dict{String, UInt32}
    streams::Vector{IOStream}
    buffers::Vector{Vector{RunSummaryRecord}}
    capacity::Int
end

function _summary_flush!(io, buffer)
    isempty(buffer) && return
    write(io, buffer)
    empty!(buffer)
end

function _scan_summary_records(f, path, budget)
    capacity = max(1, min(65536, budget ÷ (8sizeof(RunSummaryRecord))))
    buffer = Vector{RunSummaryRecord}(undef, capacity)
    open(path, "r") do io
        while !eof(io)
            n = min(capacity, (filesize(path) - position(io)) ÷ sizeof(RunSummaryRecord))
            resize!(buffer, n)
            read!(io, buffer)
            f(buffer)
        end
    end
end

_summary_float_key(x::Float32) = signbit(x) ? ~reinterpret(UInt32, x) : reinterpret(UInt32, x) ⊻ 0x80000000
_summary_key_float(x::UInt32) = reinterpret(Float32, x & 0x80000000 == 0 ? ~x : x ⊻ 0x80000000)

# Exact Float32 order statistics using four byte-wise passes over an oversized run.
function _large_summary_medians(path, n_peptides, budget)
    seen = falses(n_peptides)
    counts = zeros(Int, 8)
    has_nan = falses(8)
    _scan_summary_records(path, budget) do records
        for record in records
            seen[record.peptide] = true
            for metric in 1:8
                record.valid & (UInt32(1) << (metric-1)) == 0 && continue
                counts[metric] += 1
                has_nan[metric] |= isnan(record.values[metric])
            end
        end
    end
    ranks = [(counts[m] + k) ÷ 2 for m in 1:8, k in 1:2]
    prefixes = zeros(UInt32, 8, 2)
    bins = zeros(Int, 256, 8, 2)
    for shift in (24, 16, 8, 0)
        fill!(bins, 0)
        mask = shift == 24 ? UInt32(0) : typemax(UInt32) << (shift + 8)
        _scan_summary_records(path, budget) do records
            for record in records, metric in 1:8
                (has_nan[metric] || record.valid & (UInt32(1) << (metric-1)) == 0) && continue
                key = _summary_float_key(record.values[metric])
                for k in 1:2
                    key & mask == prefixes[metric, k] || continue
                    bins[Int((key >> shift) & 0xff) + 1, metric, k] += 1
                end
            end
        end
        for metric in 1:8, k in 1:2
            (counts[metric] == 0 || has_nan[metric]) && continue
            for bin in 1:256
                n = bins[bin, metric, k]
                if ranks[metric, k] > n
                    ranks[metric, k] -= n
                else
                    prefixes[metric, k] |= UInt32(bin-1) << shift
                    break
                end
            end
        end
    end
    medians = ntuple(8) do m
        counts[m] == 0 && return missing
        has_nan[m] && return Float32(NaN)
        median(Float32[_summary_key_float(prefixes[m, 1]), _summary_key_float(prefixes[m, 2])])
    end
    return medians, count(seen)
end

function _finish_summary_partition!(stats, path, depth, n_peptides, budget)
    n = filesize(path) ÷ sizeof(RunSummaryRecord)
    n == 0 && return rm(path)
    first_run = open(io -> read(io, UInt32), path)
    if filesize(path) <= budget ÷ 4
        records = Vector{RunSummaryRecord}(undef, n)
        open(io -> read!(io, records), path)
        sort!(records; by=r -> r.run)
        seen = falses(n_peptides)
        values = Float32[]
        first = 1
        while first <= n
            last = first
            while last < n && records[last+1].run == records[first].run
                last += 1
            end
            fill!(seen, false)
            for i in first:last
                seen[records[i].peptide] = true
            end
            s = stats[records[first].run]
            s.peptides_identified = count(seen)
            s.medians = ntuple(8) do metric
                empty!(values)
                for i in first:last
                    r = records[i]
                    r.valid & (UInt32(1) << (metric-1)) == 0 || push!(values, r.values[metric])
                end
                isempty(values) ? missing : median!(values)
            end
            first = last + 1
        end
    elseif stats[first_run].precursors_identified == n
        stats[first_run].medians, stats[first_run].peptides_identified =
            _large_summary_medians(path, n_peptides, budget)
    else
        paths = [path * ".$i" for i in 1:SUMMARY_FANOUT]
        streams = IOStream[]
        buffers = [RunSummaryRecord[] for _ in paths]
        capacity = max(1, budget ÷ (8SUMMARY_FANOUT * sizeof(RunSummaryRecord)))
        try
            for child in paths
                push!(streams, open(child, "w"))
            end
            _scan_summary_records(path, budget) do records
                for record in records
                    bucket = Int(((record.run - 1) >> (4depth)) & 0x0f) + 1
                    push!(buffers[bucket], record)
                    length(buffers[bucket]) >= capacity && _summary_flush!(streams[bucket], buffers[bucket])
                end
            end
            for i in eachindex(streams)
                _summary_flush!(streams[i], buffers[i])
            end
        finally
            foreach(close, streams)
        end
        rm(path)
        for child in paths
            _finish_summary_partition!(stats, child, depth + 1, n_peptides, budget)
        end
        return
    end
    rm(path)
end

"""
    with_run_summary(f, file_names; temp_parent=tempdir(), memory_budget_bytes=64*1024^2)

Accumulate exact summary statistics through bounded disk partitions. The callback
receives an accumulator for `accumulate_run_summary!`; finalized statistics are
returned after all temporary files have been removed. Run counters and a shared
sequence dictionary are metadata outside the numeric workspace budget.
"""
function with_run_summary(f, file_names; temp_parent=tempdir(), memory_budget_bytes=64*1024^2)
    memory_budget_bytes >= 65536 || throw(ArgumentError("Summary workspace must be at least 64 KiB"))
    mktempdir(temp_parent; prefix=".run_summary_") do dir
        stats = RunSummaryStats.(file_names)
        paths = [joinpath(dir, "$i.bin") for i in 1:SUMMARY_FANOUT]
        streams = IOStream[]
        buffers = [RunSummaryRecord[] for _ in paths]
        capacity = max(1, min(16384, memory_budget_bytes ÷ (8SUMMARY_FANOUT * sizeof(RunSummaryRecord))))
        acc = RunSummaryAccumulator(stats, Dict{String,UInt32}(), streams, buffers, capacity)
        try
            for path in paths
                push!(streams, open(path, "w"))
            end
            f(acc)
            for i in eachindex(streams)
                _summary_flush!(streams[i], buffers[i])
            end
        finally
            foreach(close, streams)
        end
        n_peptides = length(acc.peptide_ids)
        empty!(acc.buffers)
        acc.peptide_ids = Dict{String,UInt32}()
        for path in paths
            _finish_summary_partition!(stats, path, 1, n_peptides, memory_budget_bytes)
        end
        return stats
    end
end

function accumulate_run_summary!(acc::RunSummaryAccumulator, tbl)
    cols = Tables.columntable(tbl)
    _accumulate_run_summary!(acc, cols.ms_file_idx, cols.target, cols.sequence, cols.peak_area,
        get(cols, :peak_area_normalized, nothing), get(cols, :mbr_recovered, nothing),
        cols.irt_error, cols.rt_fwhm, cols.points_integrated, cols.charge, cols.missed_cleavage)
    return acc
end

function _accumulate_run_summary!(acc::RunSummaryAccumulator, ms_file_idx, target, sequence,
    peak_area, peak_area_normalized, mbr_recovered, irt_error, rt_fwhm, points_integrated,
    charge, missed_cleavage)
    for i in axes(ms_file_idx, 1)
        target[i] || continue
        run = ms_file_idx[i]
        s = acc.stats[run]
        s.precursors_identified += 1
        peptide = get!(acc.peptide_ids, sequence[i]) do
            UInt32(length(acc.peptide_ids) + 1)
        end
        mbr_recovered !== nothing && mbr_recovered[i] && (s.precursors_mbr += 1)
        area, factor, valid = 0.0f0, 0.0f0, UInt32(0xfc)
        raw_area = peak_area[i]
        if !ismissing(raw_area) && !(raw_area <= 0)
            s.precursors_quantified += 1
            s.total_peak_area += raw_area
            area = Float32(raw_area)
            valid |= 0x01
            if peak_area_normalized !== nothing
                norm = peak_area_normalized[i]
                if !ismissing(norm) && norm > 0
                    factor = Float32(norm / raw_area)
                    valid |= 0x02
                end
            end
        end
        values = (area, factor, Float32(abs(irt_error[i])), Float32(rt_fwhm[i]),
            Float32(points_integrated[i]), Float32(length(sequence[i])), Float32(charge[i]),
            Float32(missed_cleavage[i]))
        bucket = Int((run - 1) % SUMMARY_FANOUT) + 1
        push!(acc.buffers[bucket], RunSummaryRecord(run, peptide, valid, values))
        length(acc.buffers[bucket]) >= acc.capacity && _summary_flush!(acc.streams[bucket], acc.buffers[bucket])
    end
end

"""
    add_protein_group_counts!(stats, protein_groups_long, file_names)

Per-file protein-group counts from the long protein-group table: identified =
target rows for the file, quantified = those with a positive abundance.
"""
function add_protein_group_counts!(stats::Vector{RunSummaryStats},
                                   protein_groups_long::DataFrame,
                                   file_names::Vector{String})
    idx = Dict(name => i for (i, name) in enumerate(file_names))
    for row in eachrow(protein_groups_long)
        (ismissing(row.file_name) || !row.target) && continue
        i = get(idx, row.file_name, 0)
        i == 0 && continue
        stats[i].protein_groups_identified += 1
        a = row.abundance
        !ismissing(a) && a > 0 && (stats[i].protein_groups_quantified += 1)
    end
    return stats
end

function add_protein_group_counts!(stats::Vector{RunSummaryStats},
                                   path::AbstractString,
                                   file_names::Vector{String})
    for batch in Arrow.Stream(path)
        add_protein_group_counts!(stats, DataFrame(batch; copycols=false), file_names)
    end
    return stats
end

_median_or_missing(v::AbstractVector) = isempty(v) ? missing : median(v)

_mass_tol_unit(::SimpleMassErrorModel) = "ppm"
_mass_tol_unit(::LinearBiasPpmTolMassErrorModel) = "ppm"
_mass_tol_unit(::AbstractMassErrorModel) = "Da"

function _mass_tol_columns(models::Dict{Int64, AbstractMassErrorModel}, i::Int)
    haskey(models, i) || return (missing, missing, missing)
    m = models[i]
    return (getLeftTol(m), getRightTol(m), _mass_tol_unit(m))
end

# Scan-level metadata straight from the (memory-mapped) raw file.
function _raw_file_columns(path::String)
    tbl = Arrow.Table(path)
    orders = tbl[:msOrder]
    rts = tbl[:retentionTime]
    n_ms1 = count(==(UInt8(1)), orders)
    n_ms2 = length(orders) - n_ms1
    gradient = isempty(rts) ? missing : Float32(maximum(rts))
    return (gradient, n_ms1, n_ms2)
end

"""
    write_run_summary(path, stats, search_context)

Write `run_summary.tsv`: one row per raw file, in file order.
"""
function write_run_summary(path::String, stats::Vector{RunSummaryStats},
                           search_context::SearchContext)
    ms_data = getMSData(search_context)
    raw_paths = getFilePaths(ms_data)
    n = length(stats)
    ms2_tols = [_mass_tol_columns(search_context.mass_error_model, i) for i in 1:n]
    ms1_tols = [_mass_tol_columns(search_context.ms1_mass_error_model, i) for i in 1:n]
    raw = [_raw_file_columns(raw_paths[i]) for i in 1:n]
    df = DataFrame(
        file_name = [s.file_name for s in stats],
        precursors_identified = [s.precursors_identified for s in stats],
        precursors_quantified = [s.precursors_quantified for s in stats],
        precursors_mbr = [s.precursors_mbr for s in stats],
        peptides_identified = [s.peptides_identified for s in stats],
        protein_groups_identified = [s.protein_groups_identified for s in stats],
        protein_groups_quantified = [s.protein_groups_quantified for s in stats],
        total_peak_area = [s.total_peak_area for s in stats],
        median_peak_area = [s.medians[1] for s in stats],
        median_normalization_factor = [s.medians[2] for s in stats],
        median_irt_error = [s.medians[3] for s in stats],
        median_rt_fwhm = [s.medians[4] for s in stats],
        median_points_integrated = [s.medians[5] for s in stats],
        median_peptide_length = [s.medians[6] for s in stats],
        median_charge = [s.medians[7] for s in stats],
        median_missed_cleavages = [s.medians[8] for s in stats],
        ms2_mass_tol_low = [t[1] for t in ms2_tols],
        ms2_mass_tol_high = [t[2] for t in ms2_tols],
        ms2_mass_tol_unit = [t[3] for t in ms2_tols],
        ms1_mass_tol_low = [t[1] for t in ms1_tols],
        ms1_mass_tol_high = [t[2] for t in ms1_tols],
        ms1_mass_tol_unit = [t[3] for t in ms1_tols],
        gradient_length_min = [r[1] for r in raw],
        n_ms1_scans = [r[2] for r in raw],
        n_ms2_scans = [r[3] for r in raw],
    )
    add_calibration_qc_columns!(df, search_context.calibration_qc, n)
    CSV.write(path, df, delim = '\t')
    return df
end
