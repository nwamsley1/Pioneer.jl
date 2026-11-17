# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

#==========================================================
Acquisition mode (DDA / DIA) of each MS file, decided once at search start.

A converter may record the decision in the file metadata (`acquisition_mode` = dda | dia | targeted | unknown, with
`acquisition_mode_source` = vendor | inferred and `acquisition_mode_evidence`). Files without a usable label are
classified here from the data: coverage of the MS2 m/z range per MS1 cycle (`infer_acquisition_mode`).
==========================================================#

const ACQUISITION_MODES = ("dda", "dia", "targeted", "unknown")
const ACQUISITION_MODE_VERSION = "1"
# data rule constants
const ACQ_N_PROBES = 400
const ACQ_N_STARTS = 16
const ACQ_COVER_TARGET = 0.95
const ACQ_MIN_CYCLES = 20
const ACQ_MIN_EVENTS = 100      # fewer MS2 isolation events: too little to classify (unknown)
const ACQ_DIA_MAX_RATIO = 2.5

"""
    acquisition_mode_label(meta) -> Union{Nothing, String}

The converter-recorded `acquisition_mode` (dda, dia or targeted) in a file's metadata, or `nothing` when there is no
usable label (absent, `unknown` or unrecognised).
"""
function acquisition_mode_label(meta)
    meta === nothing && return nothing
    m = lowercase(strip(string(get(meta, "acquisition_mode", ""))))
    return m in ("dda", "dia", "targeted") ? m : nothing
end

"""
    coverage_ratio(windows_per_group) -> Float64

How many groups (MS1 cycles, or single isolation events) it takes to cover 95% of the MS2 m/z range, relative to the
ideal number given the median coverage of one group. `ACQ_N_PROBES` probe m/z span the 2nd-98th percentile of the
window bounds; the count is the median over `ACQ_N_STARTS` starting groups (`Inf` when the run ends first). DIA tiles
its range every cycle (about 1, whatever the window shape); DDA covers intensity-chosen spots (much more, or never).
"""
function coverage_ratio(groups::Vector{Vector{NTuple{2, Float32}}})
    lo = Float32[w[1] for g in groups for w in g]; hi = Float32[w[2] for g in groups for w in g]
    probes = collect(range(quantile(lo, 0.02), quantile(hi, 0.98); length = ACQ_N_PROBES))
    masks = map(groups) do g
        m = falses(ACQ_N_PROBES)
        for (l, h) in g
            i1 = searchsortedfirst(probes, l); i2 = searchsortedlast(probes, h)
            i1 <= i2 && (m[i1:i2] .= true)
        end
        m
    end
    per = median(count.(masks))
    per == 0 && return Inf
    cov = falses(ACQ_N_PROBES); ks = Float64[]
    for st in round.(Int, range(1, max(1, length(masks) - 1); length = ACQ_N_STARTS))
        fill!(cov, false); k = 0; done = false
        for i in st:length(masks)
            k += 1; cov .|= masks[i]
            count(cov) >= ACQ_COVER_TARGET * ACQ_N_PROBES && (done = true; break)
        end
        push!(ks, done ? k : Inf)
    end
    return median(ks) / (ACQ_N_PROBES / per)
end

"""
    infer_acquisition_mode(orders, centers, widths) -> (mode, evidence)

The data rule for files without a label. MS1 cycles (the MS2 scans between consecutive MS1 scans) are the coverage
groups; files with fewer than `ACQ_MIN_CYCLES` such cycles use isolation events instead. An event is a run of
consecutive MS2 scans with the same centre (rounded to 0.01 m/z), so the mobility slices of one DDA-PASEF isolation
or one diaPASEF window count once. Coverage per cycle, not per event, because sliding windows (Synchro-PASEF, ZT)
overlap from scan to scan: counted per event they need many steps to cover the range and look like DDA. Ratio <= `ACQ_DIA_MAX_RATIO` -> dia; otherwise dda, or targeted when isolation centres repeat
(fewer than 5% of isolation events are new centres, i.e. a fixed target list). Fewer than `ACQ_MIN_EVENTS` MS2
isolation events, or fewer than `ACQ_MIN_CYCLES` groups -> unknown.
"""
function infer_acquisition_mode(orders::AbstractVector, centers::AbstractVector, widths::AbstractVector)
    cycles = Vector{Vector{NTuple{2, Float32}}}(); cur = NTuple{2, Float32}[]
    events = Vector{Vector{NTuple{2, Float32}}}()
    keys = Set{Int}(); n_events = 0; prev = typemin(Int)
    for i in 1:length(orders)
        o = orders[i]
        if o == 1
            isempty(cur) || (push!(cycles, cur); cur = NTuple{2, Float32}[])
            continue
        end
        o == 2 || continue
        c = centers[i]; w = widths[i]
        (ismissing(c) || ismissing(w) || !isfinite(c) || !(w > 0)) && continue
        win = (Float32(c - w / 2), Float32(c + w / 2))
        push!(cur, win)
        k = round(Int, 100 * c)
        if k != prev
            n_events += 1; push!(keys, k); push!(events, [win])
        end
        prev = k
    end
    isempty(cur) || push!(cycles, cur)
    f = n_events == 0 ? NaN : length(keys) / n_events
    groups, unit = length(cycles) >= ACQ_MIN_CYCLES ? (cycles, "cycles") : (events, "events")
    (n_events < ACQ_MIN_EVENTS || length(groups) < ACQ_MIN_CYCLES) &&
        return "unknown", "too few MS2 isolation events ($(n_events)) to classify"
    r = coverage_ratio(groups)
    mode = r <= ACQ_DIA_MAX_RATIO ? "dia" : f >= 0.05 ? "dda" : "targeted"
    return mode, "coverage ratio $(round(r; digits = 2)) over $(length(groups)) $(unit), distinct-centre fraction $(round(f; digits = 4))"
end

infer_acquisition_mode(spectra::MassSpecData) =
    infer_acquisition_mode(getMsOrders(spectra), getCenterMzs(spectra), getIsolationWidthMzs(spectra))

"""
    file_acquisition_mode(spectra) -> (mode, source)

The acquisition mode of one opened MS file: the converter's label when present, else the data rule. `source` is a
short description for the log.
"""
function file_acquisition_mode(spectra::MassSpecData)
    meta = getAcquisitionMetadata(spectra)
    label = acquisition_mode_label(meta)
    if label !== nothing
        ev = get(meta, "acquisition_mode_evidence", "")
        return label, "recorded at conversion ($(get(meta, "acquisition_mode_source", "?"))$(isempty(ev) ? "" : ": " * ev))"
    end
    mode, ev = infer_acquisition_mode(spectra)
    recorded = meta === nothing ? "" : string(get(meta, "acquisition_mode", ""))
    why = isempty(recorded) ? "no label in the file" :
          "converter recorded \"$recorded\"" * (haskey(meta, "acquisition_mode_evidence") ?
                                                 " ($(meta["acquisition_mode_evidence"]))" : "")
    return mode, "inferred from the data, $why: $ev"
end

"""Configured `search.acquisition_mode`: `auto` (default; per file, from the label or the data), `dda` or `dia`."""
function configured_acquisition_mode(params)
    s = params.search
    m = hasproperty(s, :acquisition_mode) ? lowercase(String(s.acquisition_mode)) : "auto"
    m in ("auto", "dda", "dia") ||
        throw(ArgumentError("search.acquisition_mode must be \"auto\", \"dda\" or \"dia\", not \"$(s.acquisition_mode)\""))
    return m
end

"""
    detect_acquisition_modes!(search_context, params) -> Nothing

Decide every file's acquisition mode before the search (`dda_files`): `search.acquisition_mode` when set to dda/dia, else
the converter label or the data rule per file.
Scanning-quad (ZT) files are DIA. `targeted` files are searched as DDA, with a warning. `unknown` files (too little
data to classify, e.g. near-empty runs) follow the classified files of the search, or DIA when there are none. A search
mixing DDA and DIA files is refused: they are quantified and scored differently.
"""
function detect_acquisition_modes!(search_context::SearchContext, params)
    configured = configured_acquisition_mode(params)
    paths = getMSData(search_context).file_paths
    modes = Dict{Int, String}()
    for (i, path) in enumerate(paths)
        name = basename(rstrip(path, '/'))
        if configured != "auto"
            mode, why = configured, "set by search.acquisition_mode"
        else
            spectra = loadMassSpecData(path)
            mode, why = zt_mode(spectra) ? ("dia", "scanning-quad (ZT) acquisition") : file_acquisition_mode(spectra)
        end
        modes[i] = mode
        @user_info "Acquisition mode $(uppercase(mode)): $name ($why)"
        mode == "targeted" &&
            @user_warn "$name looks like a targeted acquisition (PRM / SureQuant): not supported; searching it as DDA"
    end
    n = length(paths)
    known = [i for i in 1:n if modes[i] != "unknown"]
    is_dda = Dict(i => modes[i] in ("dda", "targeted") for i in known)
    n_dda = count(i -> is_dda[i], known)
    if 0 < n_dda < length(known)
        dda_names = [basename(rstrip(paths[i], '/')) for i in known if is_dda[i]]
        throw(ArgumentError("This search mixes DDA ($(n_dda)) and DIA ($(length(known) - n_dda)) files; search them " *
                            "separately (or set search.acquisition_mode). DDA: " * join(first(dda_names, 5), ", ") *
                            (n_dda > 5 ? ", ..." : "")))
    end
    rest_dda = n_dda > 0                                   # the classified files agree: unknown ones follow them
    for i in 1:n
        if modes[i] == "unknown"
            @user_warn "$(basename(rstrip(paths[i], '/'))): acquisition mode could not be determined; searching it as " *
                       "$(rest_dda ? "DDA" : "DIA") like the rest of the search (set search.acquisition_mode to override)"
            search_context.dda_files[i] = rest_dda
        else
            search_context.dda_files[i] = is_dda[i]
        end
    end
    return nothing
end

"Whether file `ms_file_idx` is DDA (decided at search start by `detect_acquisition_modes!`)."
is_dda_file(search_context::SearchContext, ms_file_idx::Integer) = get(search_context.dda_files, Int64(ms_file_idx), false)
