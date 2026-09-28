# Scanning-quad (ZT) per-file state: whether a file is searched in ZT mode, its Q1 bin lattice,
# and the settings that follow from it. Everything else in SearchDIA reaches ZT through
# `getZTGeometry(search_context, ms_file_idx)`, which is `nothing` for every other file.

"""
   getZTGeometry(s::SearchContext, index::Integer) -> Union{Nothing, ZTGeometry}

Q1 bin lattice for a scanning-quad file, or `nothing` when the file is not a scanning
acquisition (or detection has not run for it yet). Unlike `getQuadTransmissionModel`, absence
is the normal case — most files are not scanning acquisitions — so this does not warn.
"""
getZTGeometry(s::SearchContext, index::I) where {I<:Integer} =
    get(s.zt_geometry, Int64(index), nothing)::Union{Nothing, ZTGeometry}

setZTGeometry!(s::SearchContext, index::I, g::Union{Nothing, ZTGeometry}) where {I<:Integer} =
    (s.zt_geometry[Int64(index)] = g)

"""True when `index` has been checked and is a scanning-quad acquisition."""
isZTFile(s::SearchContext, index::I) where {I<:Integer} =
    getZTGeometry(s, index) !== nothing

"""
    zt_prepare_file!(search_context, ms_file_idx, spectra)

Called by `execute_search` before each file of each stage. Detects the file's Q1 lattice on first
sight (`ensure_zt_geometry!`) and sets the deconvolution convergence tolerance for this file:
`ZT_DECONV_CONVERGENCE_TOL` on scanning-quad files, each stage's own tolerance otherwise.
"""
function zt_prepare_file!(search_context::SearchContext, ms_file_idx::Int64, spectra::MassSpecData)
    ensure_zt_geometry!(search_context, ms_file_idx, spectra)
    tol = getZTGeometry(search_context, ms_file_idx) === nothing ? NaN32 : ZT_DECONV_CONVERGENCE_TOL
    foreach(sd -> sd.deconv_tol = tol, getSearchData(search_context))
    return nothing
end

"""
    zt_mode(spectra) -> Bool

Whether to search a file in scanning-quad (ZT) mode: exactly when its converter recorded
`acquisition_type = "zt_scan_dia"`. SciexWiff writes that when the user declares the run ZT Scan
DIA at conversion (the `.wiff` does not record the scan mode). There is no config override.
"""
function zt_mode(spectra::MassSpecData)
    meta = getAcquisitionMetadata(spectra)
    return meta !== nothing && get(meta, "acquisition_type", "") == "zt_scan_dia"
end

"""
    ensure_zt_geometry!(search_context, ms_file_idx, spectra) -> Nothing

Detect and record the Q1 bin lattice for one file, the first time that file is seen.
Idempotent — later search phases hit the cache and return immediately.

Called from `execute_search`'s per-file loop, so it measures the `spectra` already open for
that iteration: one file at a time, no extra file opens. The lattice is per file: sibling runs
from the same method differ slightly in sweep start and step.

The tuning stages run with the provisional `ZT_METASCAN_K_DEFAULT`; QuadTuningSearch then sets
the search k from the fitted transmission profile (`zt_search_k`).

Every ZT file logs its status once. Non-ZT files log at debug level only, so their output
matches a build without ZT support.
"""
function ensure_zt_geometry!(search_context::SearchContext, ms_file_idx::Int64,
                             spectra::MassSpecData)
    haskey(search_context.zt_geometry, ms_file_idx) && return nothing
    fname = getParsedFileName(getMassSpecData(search_context), ms_file_idx)
    if !zt_mode(spectra)
        setZTGeometry!(search_context, ms_file_idx, nothing)
        @debug_l1 "Scanning-quad (ZT) [$fname]: OFF"
        return nothing
    end
    g = detect_zt_geometry(spectra, ZT_METASCAN_K_DEFAULT)
    if g === nothing
        setZTGeometry!(search_context, ms_file_idx, nothing)
        @user_warn "Scanning-quad (ZT) [$fname]: file is marked ZT Scan DIA, but has no usable MS2 isolation metadata — ZT off for this file"
        return nothing
    end

    setZTGeometry!(search_context, ms_file_idx, g)

    # Install the transmission model for this file up front, so every consumer
    # (BitVecCalibration, HuberTuning, precursor_fraction_transmitted, chromatogram
    # integration, quant) sees ONE coherent model instead of a per-call-site patch.
    #
    # A flat box spanning the meta-scan is the right choice for deconvolution: after expansion
    # each (precursor, scan) gets its OWN weight, so the true transmission is absorbed into that
    # weight rather than needing to be in the model. Measured on the reference ZT data the
    # profile is ~Gaussian with FWHM ~6.3-7.0 Da, and modelling it as flat mis-states only the
    # intra-scan M0/M1 ratio — by <=11% within +/-2 Da, where 53% of transmission lives.
    #
    # The box is derived from the expansion span so k stays the single control over which bins a
    # precursor participates in; a box narrower than the expansion would silently clip bins the
    # expansion deliberately included. QuadTuningSearch re-installs it with the search k.
    setQuadTransmissionModel!(search_context, ms_file_idx,
                              SquareQuadModel(zt_deconv_overhang(g)))

    tiling = zt_tiles_contiguously(g) ? "contiguous" : "NON-CONTIGUOUS (check acquisition!)"
    @user_info "Scanning-quad (ZT) [$fname]: ON (file metadata)  tuning k=$(g.metascan_k) (±$(g.metascan_k) bins, " *
               "±$(round(g.metascan_k * g.bin_step, digits=2)) Da)  S=$(round(g.bin_step, digits=4)) Da  " *
               "nominal width=$(round(g.nominal_width, digits=4)) Da ($tiling)  " *
               "bins/ramp=$(g.bins_per_ramp)  span=$(round(zt_span_mz(g), digits=2)) Da"
    return nothing
end
