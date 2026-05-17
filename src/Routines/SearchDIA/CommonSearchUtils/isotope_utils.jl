# Isotope capture and transmission utilities.
#
# DataFrame-level wrappers around isotopeSplines.jl functions
# (getPrecursorIsotopeSet, getPrecursorFractionTransmitted!).
# Shared by MainSearch and IntegrateChromatogramsSearch.

#############################################################################
# Per-row isotope computation (shared inner loop)
#############################################################################

"""
    _compute_scan_window(scan_id, centerMz, isolationWidthMz)

Extract scan center m/z and half-width for a given scan, handling missing values.
Returns (low_mz, high_mz, scan_mz, window_width).
"""
@inline function _compute_scan_window(scan_id::UInt32,
                                       centerMz::AbstractVector{Union{Missing, Float32}},
                                       isolationWidthMz::AbstractVector{Union{Missing, Float32}})
    scan_mz = coalesce(centerMz[scan_id], zero(Float32))::Float32
    window_width = coalesce(isolationWidthMz[scan_id], zero(Float32))::Float32
    low_mz = Float32(scan_mz - window_width / 2)
    high_mz = Float32(scan_mz + window_width / 2)
    return low_mz, high_mz, scan_mz, window_width
end

"""
    _compute_fraction_transmitted(iso_splines, quad_transmission_model,
                                  mz, charge, sulfur, scan_mz, window_width,
                                  precursor_transmission)

Compute the fraction of precursor signal transmitted through the quadrupole.
Uses the pre-allocated `precursor_transmission` buffer to avoid allocation.
"""
@inline function _compute_fraction_transmitted(
    iso_splines::IsotopeSplineModel{40, Float32},
    quad_transmission_model::QuadTransmissionModel,
    mz::Float32, charge::UInt8, sulfur::UInt8,
    scan_mz::Float32, window_width::Float32,
    precursor_transmission::Vector{Float32})

    return getPrecursorFractionTransmitted!(
        precursor_transmission,
        iso_splines,
        (1, 5),
        getQuadTransmissionFunction(quad_transmission_model, scan_mz, window_width),
        mz, charge, sulfur)
end

#############################################################################
# Public API
#############################################################################

"""
    quant_isotope_range_mask(isotopes_captured) -> UInt32

Encode one captured precursor-isotope range as a bit in a compact mask.
Valid ranges are `(lo, hi)` with `0 <= lo <= hi <= MAX_PRECURSOR_ISOTOPE`.
Invalid or empty ranges encode to zero.
"""
@inline function quant_isotope_range_mask(isotopes_captured::Tuple{Int8, Int8})::UInt32
    lo = Int(first(isotopes_captured))
    hi = Int(last(isotopes_captured))
    (lo < 0 || hi < lo || hi > MAX_PRECURSOR_ISOTOPE) && return UInt32(0)

    bit_idx = 0
    if lo > 0
        @inbounds for first_iso in 0:(lo - 1)
            bit_idx += MAX_PRECURSOR_ISOTOPE - first_iso + 1
        end
    end
    bit_idx += hi - lo
    return UInt32(1) << bit_idx
end

@inline function quant_isotope_range_allowed(mask::UInt32, isotopes_captured::Tuple{Int8, Int8})::Bool
    bit = quant_isotope_range_mask(isotopes_captured)
    return bit != UInt32(0) && (mask & bit) != UInt32(0)
end

"""
    add_quant_isotope_masks_from_scores!(psms, score_col, pep_threshold)

For each precursor, build the quantification isotope-range allowlist from rows
whose per-run model score passes the current MainSearch PEP threshold. The
allowlist is stored as `:quant_isotope_mask` on every row for that precursor.
"""
function add_quant_isotope_masks_from_scores!(
    psms::DataFrame,
    score_col::Symbol,
    pep_threshold::Float32,
)
    n = nrow(psms)
    if n == 0
        psms[!, :quant_isotope_mask] = UInt32[]
        return psms
    end

    prec_col = psms[!, :precursor_idx]::AbstractVector{UInt32}
    iso_col = psms[!, :isotopes_captured]::AbstractVector{Tuple{Int8, Int8}}

    keep = trues(n)
    if pep_threshold < 1.0f0
        probs = Float32.(psms[!, score_col])
        targets = Vector{Bool}(psms[!, :target])
        peps = Vector{Float32}(undef, n)
        get_PEP!(probs, targets, peps; doSort=true, fdr_scale_factor=1.0f0)
        keep .= peps .<= pep_threshold
    end

    masks = Dict{UInt32, UInt32}()
    @inbounds for i in 1:n
        keep[i] || continue
        bit = quant_isotope_range_mask(iso_col[i])
        bit == UInt32(0) && continue
        pid = prec_col[i]
        masks[pid] = get(masks, pid, UInt32(0)) | bit
    end

    psms[!, :quant_isotope_mask] = UInt32[get(masks, pid, UInt32(0)) for pid in prec_col]
    return psms
end

"""
    apply_quant_isotope_masks_to_chromatograms!(chromatograms, passing_psms)

Apply precursor-level quant isotope allowlists to chromatogram rows by removing
disallowed isotope ranges before integration. Removing rows avoids leaving
zero-transmission duplicate RT points in combined-trace chromatograms.
If `passing_psms` lacks `:quant_isotope_mask`, this is a no-op for backwards
compatibility with older intermediate files.
"""
function apply_quant_isotope_masks_to_chromatograms!(
    chromatograms::DataFrame,
    passing_psms::DataFrame,
)
    empty_summary = (
        chrom_rows = nrow(chromatograms),
        kept_rows = nrow(chromatograms),
        psm_precursors = hasproperty(passing_psms, :precursor_idx) ?
            length(Set(passing_psms[!, :precursor_idx])) : 0,
        zero_mask_precursors = 0,
        disallowed_rows = 0,
        all_disallowed_precursors = 0,
        partially_disallowed_precursors = 0,
    )
    hasproperty(passing_psms, :quant_isotope_mask) || return empty_summary
    (nrow(chromatograms) == 0 || !hasproperty(chromatograms, :isotopes_captured)) &&
        return empty_summary

    psm_prec_col = passing_psms[!, :precursor_idx]::AbstractVector{UInt32}
    psm_mask_col = passing_psms[!, :quant_isotope_mask]::AbstractVector{UInt32}
    mask_by_pid = Dict{UInt32, UInt32}()
    sizehint!(mask_by_pid, length(psm_prec_col))
    @inbounds for i in eachindex(psm_prec_col)
        mask_by_pid[psm_prec_col[i]] = psm_mask_col[i]
    end

    chrom_prec_col = chromatograms[!, :precursor_idx]::AbstractVector{UInt32}
    chrom_iso_col = chromatograms[!, :isotopes_captured]::AbstractVector{Tuple{Int8, Int8}}
    total_rows_by_pid = Dict{UInt32, Int}()
    allowed_rows_by_pid = Dict{UInt32, Int}()
    keep_rows = trues(nrow(chromatograms))
    disallowed_rows = 0

    @inbounds for i in eachindex(chrom_prec_col)
        pid = chrom_prec_col[i]
        total_rows_by_pid[pid] = get(total_rows_by_pid, pid, 0) + 1
        mask = get(mask_by_pid, pid, UInt32(0))
        if quant_isotope_range_allowed(mask, chrom_iso_col[i])
            allowed_rows_by_pid[pid] = get(allowed_rows_by_pid, pid, 0) + 1
        else
            keep_rows[i] = false
            disallowed_rows += 1
        end
    end

    all_disallowed_precursors = 0
    partially_disallowed_precursors = 0
    for (pid, total_rows) in total_rows_by_pid
        allowed_rows = get(allowed_rows_by_pid, pid, 0)
        if allowed_rows == 0
            all_disallowed_precursors += 1
        elseif allowed_rows < total_rows
            partially_disallowed_precursors += 1
        end
    end

    if disallowed_rows > 0
        deleteat!(chromatograms, findall(!, keep_rows))
    end

    return (
        chrom_rows = length(keep_rows),
        kept_rows = length(keep_rows) - disallowed_rows,
        psm_precursors = length(mask_by_pid),
        zero_mask_precursors = count(==(UInt32(0)), values(mask_by_pid)),
        disallowed_rows = disallowed_rows,
        all_disallowed_precursors = all_disallowed_precursors,
        partially_disallowed_precursors = partially_disallowed_precursors,
    )
end

"""
    get_isotopes_captured!(chroms, quad_transmission_model, search_data,
                           scan_idx, prec_charge, prec_mz, sulfur_count,
                           centerMz, isolationWidthMz;
                           compute_isotope_set=true)

Compute precursor fraction transmitted and optionally which isotopes are
captured in each isolation window. Adds columns to `chroms`:
- `:precursor_fraction_transmitted` (always)
- `:isotopes_captured` (when `compute_isotope_set=true`)

Multi-threaded: partitions rows into chunks and spawns tasks.
"""
function get_isotopes_captured!(chroms::DataFrame,
                                quad_transmission_model::QuadTransmissionModel,
                                search_data::Vector{SimpleLibrarySearch{IsotopeSplineModel{40, Float32}}},
                                scan_idx::AbstractVector{UInt32},
                                prec_charge::AbstractArray{UInt8},
                                prec_mz::AbstractArray{Float32},
                                sulfur_count::AbstractArray{UInt8},
                                centerMz::AbstractVector{Union{Missing, Float32}},
                                isolationWidthMz::AbstractVector{Union{Missing, Float32}};
                                compute_isotope_set::Bool = true)
    n = size(chroms, 1)
    precursor_fraction_transmitted = Vector{Float32}(undef, n)
    isotopes_captured = compute_isotope_set ? Vector{Tuple{Int8, Int8}}(undef, n) : nothing

    # Cache column vector (avoid DataFrame row-indexing in hot loop)
    prec_idx_col = chroms[!, :precursor_idx]

    parallel_foreach!(n) do chunk
        thread_id = (first(chunk) % Threads.nthreads()) + 1
        iso_splines = getIsoSplines(search_data[thread_id])
        precursor_transmission = zeros(Float32, 5)

        for i in chunk
            prec_id = prec_idx_col[i]
            mz = prec_mz[prec_id]
            charge = prec_charge[prec_id]
            sulfur = sulfur_count[prec_id]
            scan_id = scan_idx[i]

            low_mz, high_mz, scan_mz, window_width = _compute_scan_window(
                scan_id, centerMz, isolationWidthMz)

            precursor_fraction_transmitted[i] = _compute_fraction_transmitted(
                iso_splines, quad_transmission_model,
                mz, charge, sulfur, scan_mz, window_width,
                precursor_transmission)

            if compute_isotope_set
                isotopes_captured[i] = getPrecursorIsotopeSet(mz, charge, low_mz, high_mz)
            end
        end
    end

    chroms[!, :precursor_fraction_transmitted] = precursor_fraction_transmitted
    if compute_isotope_set
        chroms[!, :isotopes_captured] = isotopes_captured
    end
    return nothing
end
