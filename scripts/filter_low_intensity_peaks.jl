#!/usr/bin/env julia
#
# Filter low-intensity peaks from MS2 spectra in Pioneer .arrow files.
# For each MS2 scan, removes peaks with intensity < 1e-4 * sum(all intensities in that scan).
# MS1 scans are left unchanged. Original files are not modified.
#

using Arrow, Tables

const SRC_DIR = "/Users/nathanwamsley/Data/For_Figures/OlsenAstralThreeProteome200ng"
const DST_DIR = "/Users/nathanwamsley/Data/For_Figures/OlsenAstralThreeProteome200ng_filtered1e4"
const THRESHOLD_FRACTION = 1e-4

mkpath(DST_DIR)

arrow_files = filter(f -> endswith(f, ".arrow"), readdir(SRC_DIR))
println("Found $(length(arrow_files)) arrow files to process")

for fname in arrow_files
    src_path = joinpath(SRC_DIR, fname)
    dst_path = joinpath(DST_DIR, fname)
    println("\nProcessing: $fname")

    tbl = Arrow.Table(src_path)
    ms_orders = tbl[:msOrder]
    n_scans = length(ms_orders)
    println("  $n_scans total scans")

    # Build filtered mz/intensity columns
    orig_mz = tbl[:mz_array]
    orig_int = tbl[:intensity_array]
    new_mz = Vector{Vector{Union{Missing, Float32}}}(undef, n_scans)
    new_int = Vector{Vector{Union{Missing, Float32}}}(undef, n_scans)

    n_peaks_before = 0
    n_peaks_after = 0
    n_ms2 = 0

    for i in 1:n_scans
        mzs = orig_mz[i]
        ints = orig_int[i]
        n_peaks_before += length(ints)

        if ms_orders[i] >= 2 && length(ints) > 0
            n_ms2 += 1
            total_intensity = sum(skipmissing(ints))
            cutoff = Float32(THRESHOLD_FRACTION * total_intensity)

            keep = [!ismissing(v) && v >= cutoff for v in ints]
            new_mz[i] = collect(mzs[keep])
            new_int[i] = collect(ints[keep])
        else
            new_mz[i] = collect(mzs)
            new_int[i] = collect(ints)
        end

        n_peaks_after += length(new_int[i])
    end

    pct_removed = round(100.0 * (1 - n_peaks_after / max(1, n_peaks_before)), digits=1)
    println("  MS2 scans: $n_ms2, peaks: $n_peaks_before → $n_peaks_after ($pct_removed% removed)")

    # Build columnar output — copy all columns, replacing mz/intensity
    col_names = Tables.columnnames(tbl)
    columns = map(col_names) do cn
        if cn == :mz_array
            new_mz
        elseif cn == :intensity_array
            new_int
        else
            collect(Tables.getcolumn(tbl, cn))
        end
    end
    out_table = Tables.columntable(NamedTuple{Tuple(col_names)}(Tuple(columns)))

    Arrow.write(dst_path, out_table)
    println("  Written: $dst_path")
end

println("\nDone! Filtered files in: $DST_DIR")
