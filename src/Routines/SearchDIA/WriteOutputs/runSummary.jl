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
`MaxLFQSearch` already streams into `precursors_long.arrow`, so the summary
adds no extra pass over the data.
"""

"""
    RunSummaryStats(file_name)

Accumulator for one raw file. Counts are updated in place; the vectors hold
the per-precursor values whose medians are reported.
"""
mutable struct RunSummaryStats
    file_name::String
    precursors_identified::Int
    precursors_quantified::Int
    precursors_mbr::Int
    peptides::Set{String}
    protein_groups_identified::Int
    protein_groups_quantified::Int
    total_peak_area::Float64
    peak_areas::Vector{Float32}
    normalization_factors::Vector{Float32}
    irt_errors::Vector{Float32}
    rt_fwhms::Vector{Float32}
    points_integrated::Vector{Float32}
    peptide_lengths::Vector{Float32}
    charges::Vector{Float32}
    missed_cleavages::Vector{Float32}
end

RunSummaryStats(file_name::String) = RunSummaryStats(
    file_name, 0, 0, 0, Set{String}(), 0, 0, 0.0,
    Float32[], Float32[], Float32[], Float32[], Float32[], Float32[], Float32[], Float32[])

"""
    accumulate_run_summary!(stats, tbl)

Fold one chunk of identified precursors (the table `MaxLFQSearch` writes to
`precursors_long.arrow`) into the per-file accumulators. Decoy rows are
skipped; `stats` is indexed by `ms_file_idx`.
"""
function accumulate_run_summary!(stats::Vector{RunSummaryStats}, tbl)
    cols = Tables.columntable(tbl)
    _accumulate_run_summary!(
        stats,
        cols.ms_file_idx, cols.target, cols.sequence, cols.peak_area,
        haskey(cols, :peak_area_normalized) ? cols.peak_area_normalized : nothing,
        haskey(cols, :mbr_recovered) ? cols.mbr_recovered : nothing,
        cols.irt_error, cols.rt_fwhm, cols.points_integrated,
        cols.charge, cols.missed_cleavage)
    return stats
end

# Function barrier: the column types are only known once the table is opened.
function _accumulate_run_summary!(
    stats::Vector{RunSummaryStats},
    ms_file_idx::AbstractVector, target::AbstractVector, sequence::AbstractVector,
    peak_area::AbstractVector, peak_area_normalized, mbr_recovered,
    irt_error::AbstractVector, rt_fwhm::AbstractVector, points_integrated::AbstractVector,
    charge::AbstractVector, missed_cleavage::AbstractVector)
    # ChainedVector eachindex returns indices tied to one column; use row numbers.
    for i in axes(ms_file_idx, 1)
        target[i] || continue
        s = stats[ms_file_idx[i]]
        s.precursors_identified += 1
        push!(s.peptides, String(sequence[i]))
        push!(s.irt_errors, Float32(abs(irt_error[i])))
        push!(s.rt_fwhms, Float32(rt_fwhm[i]))
        push!(s.points_integrated, Float32(points_integrated[i]))
        push!(s.peptide_lengths, Float32(length(sequence[i])))
        push!(s.charges, Float32(charge[i]))
        push!(s.missed_cleavages, Float32(missed_cleavage[i]))
        mbr_recovered !== nothing && mbr_recovered[i] && (s.precursors_mbr += 1)
        # Quantified == the same rule that blanks peak_area in the output tables.
        area = peak_area[i]
        (ismissing(area) || area <= 0) && continue
        s.precursors_quantified += 1
        s.total_peak_area += area
        push!(s.peak_areas, Float32(area))
        if peak_area_normalized !== nothing
            norm = peak_area_normalized[i]
            !ismissing(norm) && norm > 0 && push!(s.normalization_factors, Float32(norm / area))
        end
    end
    return stats
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
        peptides_identified = [length(s.peptides) for s in stats],
        protein_groups_identified = [s.protein_groups_identified for s in stats],
        protein_groups_quantified = [s.protein_groups_quantified for s in stats],
        total_peak_area = [s.total_peak_area for s in stats],
        median_peak_area = [_median_or_missing(s.peak_areas) for s in stats],
        median_normalization_factor = [_median_or_missing(s.normalization_factors) for s in stats],
        median_irt_error = [_median_or_missing(s.irt_errors) for s in stats],
        median_rt_fwhm = [_median_or_missing(s.rt_fwhms) for s in stats],
        median_points_integrated = [_median_or_missing(s.points_integrated) for s in stats],
        median_peptide_length = [_median_or_missing(s.peptide_lengths) for s in stats],
        median_charge = [_median_or_missing(s.charges) for s in stats],
        median_missed_cleavages = [_median_or_missing(s.missed_cleavages) for s in stats],
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
    CSV.write(path, df, delim = '\t')
    return df
end
