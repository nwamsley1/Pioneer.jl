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
    zt_prepare_file!(search_context, params, ms_file_idx, spectra)

Called by `execute_search` before each file of each stage. Detects the file's Q1 lattice on first
sight (`ensure_zt_geometry!`) and sets the deconvolution convergence tolerance for this file:
`ZT_DECONV_CONVERGENCE_TOL` on scanning-quad files, each stage's own tolerance otherwise.
"""
function zt_prepare_file!(search_context::SearchContext, params::PioneerParameters,
                          ms_file_idx::Int64, spectra::MassSpecData)
    ensure_zt_geometry!(search_context, params, ms_file_idx, spectra)
    tol = getZTGeometry(search_context, ms_file_idx) === nothing ? NaN32 : ZT_DECONV_CONVERGENCE_TOL
    foreach(sd -> sd.deconv_tol = tol, getSearchData(search_context))
    return nothing
end

"""
    zt_mode(acq, spectra) -> (on::Bool, source::String, file_is_zt::Bool)

Whether to search a file in scanning-quad (ZT) mode. `acquisition.scanning_quad` decides when the config sets
it; otherwise the file's own metadata does (SciexWiff records `acquisition_type = "zt_scan_dia"` when it
converts a ZT Scan DIA run). Files without metadata (e.g. msConvert → `convertMzML`) stay off unless configured.
"""
function zt_mode(acq, spectra::MassSpecData)
    meta = getAcquisitionMetadata(spectra)
    file_is_zt = meta !== nothing && get(meta, "acquisition_type", "") == "zt_scan_dia"
    hasproperty(acq, :scanning_quad) && return (Bool(acq.scanning_quad), "config", file_is_zt)
    return (file_is_zt, "file metadata", file_is_zt)
end

"""
    ensure_zt_geometry!(search_context, params, ms_file_idx, spectra) -> Nothing

Detect and record the Q1 bin lattice for one file, the first time that file is seen.
Idempotent — later search phases hit the cache and return immediately.

Called from `execute_search`'s per-file loop, so it measures the `spectra` already open for
that iteration: one file at a time, no extra file opens.

Activation is `acquisition.scanning_quad` (default `false`); `acquisition.metascan_k` sets the
expansion half-width in bins (default 6). The physical swept width is NOT recoverable from the
data — `scanHeader` is empty and `lowMz`/`highMz` are the fragment range — so `k` is a config
input, not an inference.

The measured lattice is PER FILE: sibling runs from the same acquisition method differ slightly
in sweep start and step. `metascan_k` is a config value and so is the same for every file.

Every ZT file logs its status once, ON or OFF — a silent no-op of the ZT path is the single most
expensive bug this work has hit. Non-ZT files log OFF at debug level only, so their output
matches a build without ZT support.
"""
function ensure_zt_geometry!(
    search_context::SearchContext,
    params::PioneerParameters,
    ms_file_idx::Int64,
    spectra::MassSpecData,
)
    haskey(search_context.zt_geometry, ms_file_idx) && return nothing

    fname = getParsedFileName(getMassSpecData(search_context), ms_file_idx)
    acq = params.acquisition

    zt_on, source, file_is_zt = zt_mode(acq, spectra)
    if !zt_on
        setZTGeometry!(search_context, ms_file_idx, nothing)
        if file_is_zt
            @user_info "Scanning-quad (ZT) [$fname]: OFF (file is ZT Scan DIA, but acquisition.scanning_quad = false)"
        else
            @debug_l1 "Scanning-quad (ZT) [$fname]: OFF"
        end
        return nothing
    end

    # `acquisition.metascan_k` absent (or 0) means DERIVE it: the tuning stages run with the
    # default 6, then QuadTuningSearch replaces it with `k_implied = round(h / bin_step)` from
    # the fitted transmission profile. An explicit value is authoritative and only warned on.
    _k_cfg = hasproperty(acq, :metascan_k) ? Int(acq.metascan_k) : 0
    k = _k_cfg > 0 ? _k_cfg : ZT_METASCAN_K_DEFAULT
    fwhm = hasproperty(acq, :transmission_fwhm_mz) ? Float32(acq.transmission_fwhm_mz) :
                                                     ZT_TRANSMISSION_FWHM_DEFAULT
    g = detect_zt_geometry(spectra, k, fwhm; metascan_k_derived = _k_cfg <= 0)
    if g === nothing
        setZTGeometry!(search_context, ms_file_idx, nothing)
        @user_warn "Scanning-quad (ZT) [$fname]: ON, but no usable MS2 isolation metadata — ZT off for this file"
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
    # The box is derived from the expansion span so `metascan_k` stays the single control over
    # which bins a precursor participates in; a box narrower than the expansion would silently
    # clip bins the expansion deliberately included.
    setQuadTransmissionModel!(search_context, ms_file_idx,
                              SquareQuadModel(zt_deconv_overhang(g)))

    tiling = zt_tiles_contiguously(g) ? "contiguous" : "NON-CONTIGUOUS (check acquisition!)"
    @user_info "Scanning-quad (ZT) [$fname]: ON ($source)  k=$(g.metascan_k) (±$(g.metascan_k) bins, " *
               "±$(round(g.metascan_k * g.bin_step, digits=2)) Da)  S=$(round(g.bin_step, digits=4)) Da  " *
               "nominal width=$(round(g.nominal_width, digits=4)) Da ($tiling)  " *
               "bins/ramp=$(g.bins_per_ramp)  span=$(round(zt_span_mz(g), digits=2)) Da"
    return nothing
end
