"""Per-file calibration QC. Statuses and reason flags are stored without strings."""
@enum CalibrationQCStatus::UInt8 QC_NOT_ASSESSED QC_NORMAL QC_SUSPICIOUS QC_WARNING QC_FAILED
const CALIBRATION_QC_STAGES = (:rt, :ms2_mass, :ms1_mass, :quadrupole, :nce, :ion_mobility)
const CALIBRATION_QC_BASELINE_LIMIT = 50
const CALIBRATION_QC_SUSPICIOUS_LIMIT = 50
const QC_LOW_SUPPORT = UInt16(1)
const QC_FALLBACK = UInt16(2)
const QC_INVALID = UInt16(4)
const QC_METRIC_FLAGS = (UInt16(8), UInt16(16), UInt16(32), UInt16(64))

struct CalibrationQCRecord
    status::CalibrationQCStatus
    reasons::UInt16
    n::UInt32
    metrics::NTuple{4,Float32}
end
CalibrationQCRecord() = CalibrationQCRecord(QC_NOT_ASSESSED, 0, 0, (NaN32, NaN32, NaN32, NaN32))

mutable struct CalibrationQCState
    records::NTuple{6,Vector{CalibrationQCRecord}}
    selected::NTuple{6,Set{Int}}
    pages::NTuple{6,Vector{String}}
end
CalibrationQCState() = CalibrationQCState(ntuple(_ -> CalibrationQCRecord[], 6),
    ntuple(_ -> Set{Int}(), 6), ntuple(_ -> String[], 6))
_qc_stage(stage::Symbol) = something(findfirst(==(stage), CALIBRATION_QC_STAGES))
_qc_name(s::CalibrationQCStatus) = ("not_assessed", "normal", "suspicious", "warning", "failed")[Int(s)+1]

# Conservative within-file screening heuristics, not fit acceptance criteria.
# Coverage is a fraction; other units are encoded by the metric names below.
const CALIBRATION_QC_METRICS = (
    (:rt_coverage, :residual_mad_fraction, :residual_trend_fraction, :unused),
    (:unused, :persistent_biased_peptide_fraction, :global_median_bias_over_tolerance, :outside_tolerance_fraction),
    (:unused, :residual_mad_ppm, :residual_trend_ppm, :unused),
    (:edge_coverage, :log_ratio_rmse, :unused, :parameter_at_bound),
    (:fitted_charge_coverage, :unused, :weak_nce_fraction, :endpoint_fraction),
    (:unused, :unused, :unused, :unused),   # ion mobility: support and fallback only
)
const CALIBRATION_QC_LIMITS = (
    (0.5f0, 0.05f0, 0.01f0, Inf32),
    (-Inf32, 0.10f0, 0.25f0, 0.10f0),
    (-Inf32, 10.0f0, 5.0f0, Inf32),
    (0.5f0, 0.5f0, Inf32, 0.5f0),
    (0.5f0, Inf32, 0.5f0, 0.5f0),
    (-Inf32, Inf32, Inf32, Inf32),
)

function assess_calibration_qc(stage, n, metrics; min_support=1, fallback=false, failed=false)
    values = NTuple{4,Float32}(metrics)
    reasons = UInt16(0)
    n < min_support && (reasons |= QC_LOW_SUPPORT)
    fallback && (reasons |= QC_FALLBACK)
    failed && (reasons |= QC_INVALID)
    limits = CALIBRATION_QC_LIMITS[_qc_stage(stage)]
    for i in 1:4
        isfinite(values[i]) || continue
        (i == 1 ? values[i] < limits[i] : values[i] > limits[i]) && (reasons |= QC_METRIC_FLAGS[i])
    end
    status = failed ? QC_FAILED : (reasons & (QC_LOW_SUPPORT | QC_FALLBACK)) != 0 ? QC_WARNING :
        reasons != 0 ? QC_SUSPICIOUS : QC_NORMAL
    return CalibrationQCRecord(status, reasons, UInt32(n), values)
end

function record_calibration_qc!(state::CalibrationQCState, stage, file_idx, record)
    records = state.records[_qc_stage(stage)]
    while length(records) < file_idx
        push!(records, CalibrationQCRecord())
    end
    records[file_idx] = record
    return record
end
function calibration_qc_record(state, stage, file_idx)
    records = state.records[_qc_stage(stage)]
    return file_idx <= length(records) ? records[file_idx] : CalibrationQCRecord()
end

"""Called in input-file order, once assessment for the stage is complete."""
function select_calibration_plot!(state::CalibrationQCState, stage, file_idx)
    selected = state.selected[_qc_stage(stage)]
    file_idx in selected && return true
    record = calibration_qc_record(state, stage, file_idx)
    suspicious = record.status in (QC_SUSPICIOUS, QC_WARNING, QC_FAILED)
    if file_idx <= CALIBRATION_QC_BASELINE_LIMIT ||
       (suspicious && count(>(CALIBRATION_QC_BASELINE_LIMIT), selected) < CALIBRATION_QC_SUSPICIOUS_LIMIT)
        push!(selected, file_idx)
        return true
    end
    return false
end

function calibration_qc_reason(stage, record)
    record.status == QC_NORMAL && return ""
    record.status == QC_NOT_ASSESSED && return "Calibration skipped, unavailable, or not reached"
    reasons = String[]
    record.reasons & QC_LOW_SUPPORT != 0 && push!(reasons, "Insufficient support (n=$(record.n))")
    record.reasons & QC_FALLBACK != 0 && push!(reasons, "Fallback calibration")
    record.reasons & QC_INVALID != 0 && push!(reasons, "Calibration failed or produced invalid diagnostics")
    idx = _qc_stage(stage)
    for i in 1:4
        record.reasons & QC_METRIC_FLAGS[i] == 0 && continue
        comparator = i == 1 ? "<" : ">"
        push!(reasons, "$(CALIBRATION_QC_METRICS[idx][i])=$(round(record.metrics[i], sigdigits=4)) $comparator $(CALIBRATION_QC_LIMITS[idx][i])")
    end
    return join(reasons, "; ")
end

function calibration_qc_title(context, stage, file_idx, name)
    record = calibration_qc_record(context.calibration_qc, stage, file_idx)
    reason = calibration_qc_reason(stage, record)
    return name * "\nQC: " * _qc_name(record.status) * (isempty(reason) ? "" : "\n" * replace(reason, "; " => "\n"))
end

function write_calibration_page!(context, stage, plot)
    pages = context.calibration_qc.pages[_qc_stage(stage)]
    directory = joinpath(getDataOutDir(context), "qc_plots", ".calibration_pages", string(stage))
    mkpath(directory)
    path = joinpath(directory, "$(length(pages)+1).png")
    Plots.savefig(plot, path)
    push!(pages, path)
    return nothing
end

function finish_calibration_report!(context, stage, path)
    pages = context.calibration_qc.pages[_qc_stage(stage)]
    if !isempty(pages)
        try
            mkpath(dirname(path))
            PDFGenerator.write_pngs_to_pdf(pages, path)
            foreach(p -> rm(p; force=true), pages)
            empty!(pages)
        catch err
            @user_warn "Calibration report failed for $stage" exception=(err, catch_backtrace())
        end
    end
    n = length(getFilePaths(getMSData(context)))
    counts = zeros(Int, 5)
    for i in 1:n
        counts[Int(calibration_qc_record(context.calibration_qc, stage, i).status)+1] += 1
    end
    @debug_l1 "Calibration QC $stage: normal=$(counts[2]) suspicious=$(counts[3]) warning=$(counts[4]) failed=$(counts[5]) not_assessed=$(counts[1])"
    return nothing
end

function render_calibration_safely(render, context, stage, file_idx)
    try
        render()
    catch err
        @user_warn "Calibration plot failed for $stage, file $file_idx" exception=(err, catch_backtrace())
    end
    return nothing
end
function calibration_notice!(context, stage, file_idx)
    render_calibration_safely(context, stage, file_idx) do
        title = calibration_qc_title(context, stage, file_idx, getParsedFileName(context, file_idx))
        p = Plots.plot(; title, axis=false, ticks=false, legend=false, size=(900, 450))
        write_calibration_page!(context, stage, p)
    end
end

function rt_calibration_qc(rts, irts, model, scan_rts; min_support)
    n = length(rts)
    if n == 0 || model isa IdentityModel
        return assess_calibration_qc(:rt, n, (NaN,NaN,NaN,NaN); min_support, fallback=true)
    end
    residuals = Float64[irts[i] - model(rts[i]) for i in eachindex(rts)]
    valid = all(isfinite, residuals) && all(isfinite, rts) && all(isfinite, irts)
    valid || return assess_calibration_qc(:rt,n,(NaN,NaN,NaN,NaN); min_support, failed=true)
    rt_span = maximum(scan_rts)-minimum(scan_rts)
    irt_span = maximum(irts)-minimum(irts)
    if !isfinite(rt_span) || !isfinite(irt_span) || rt_span <= 0 || irt_span <= 0
        return assess_calibration_qc(:rt,n,(NaN,NaN,NaN,NaN); min_support, failed=true)
    end
    center = median(residuals)
    trend = _qc_binned_bias(rts, residuals)
    spread = 1.4826 * median!(abs.(residuals .- center))
    metrics = ((maximum(rts)-minimum(rts))/rt_span, spread/irt_span, max(abs(center),trend)/irt_span, NaN)
    return assess_calibration_qc(:rt,n,metrics; min_support)
end

# Thresholds are screening defaults, not validated fit acceptance criteria.
const MS2_QC_MIN_PEPTIDES = 50
const MS2_QC_BIAS_LIMIT = 0.25

# Approximate 95% distribution-free median interval, using binomial ranks.
# Each entry is one peptide median, so repeated fragments do not inflate support.
function _qc_supported_median(values)
    v = sort(Float64.(values))
    n = length(v)
    isempty(v) && return (median = NaN, low = NaN, high = NaN)
    rank = max(1, floor(Int, n / 2 - 1.96 * sqrt(n) / 2))
    return (median = median(v), low = v[rank], high = v[n - rank + 1])
end

_qc_bias_direction(m) = abs(m.median) <= MS2_QC_BIAS_LIMIT ? 0 :
    m.low > 0 ? 1 : m.high < 0 ? -1 : 0

"""
    ms2_calibration_diagnostics(samples, models; group_ids)

Assess robust residual centering on peptide medians. `models[i]` must exclude
peptide `group_ids[i]` from its fit for an independent assessment. Equal-count
peptide bins cover m/z, log2 intensity and RT. Local bias requires two adjacent
bins in the same direction, at least 50 peptides per bin, and more than 10% of
peptides affected on an axis. Report actual-window coverage separately.
"""
function ms2_calibration_diagnostics(samples, models; group_ids, fallback=false)
    n = length(samples)
    length(models) == length(group_ids) == n || throw(DimensionMismatch("MS2 QC input lengths differ"))
    groups = Dict{Any, Vector{Int}}()
    residuals = Vector{Float64}(undef, n)
    errors_mda = similar(residuals)
    outside = falses(n)
    for (i, s) in enumerate(samples)
        corrected, low, high = getCorrectedMzAndBounds(models[i], s.observed_mz, s.intensity, s.rt)
        width = (high - low) / 2
        if !isfinite(width) || width <= 0 || !isfinite(corrected) ||
           !isfinite(s.theoretical_mz) || !isfinite(s.intensity) || s.intensity <= 0 || !isfinite(s.rt)
            return (record = assess_calibration_qc(:ms2_mass,n,(NaN,NaN,NaN,NaN); failed=true),
                    bins = NamedTuple[], n_peptides = 0)
        end
        residuals[i] = (corrected - s.theoretical_mz) / width
        isfinite(residuals[i]) || return (
            record=assess_calibration_qc(:ms2_mass,n,(NaN,NaN,NaN,NaN); failed=true),
            bins=NamedTuple[], n_peptides=0)
        errors_mda[i] = 1000 * (corrected - s.theoretical_mz)
        outside[i] = !(low <= s.theoretical_mz <= high)
        push!(get!(groups, group_ids[i], Int[]), i)
    end
    peptide_rows = collect(values(groups))
    ng = length(peptide_rows)
    peptide_residuals = [median(residuals[rows]) for rows in peptide_rows]
    global_bias = _qc_supported_median(peptide_residuals)
    supported_global_bias = ng >= MS2_QC_MIN_PEPTIDES && _qc_bias_direction(global_bias) != 0 ?
        abs(global_bias.median) : 0.0
    # Give every peptide equal weight in coverage, regardless of fragment count.
    outside_fraction = ng == 0 ? NaN : mean(mean(outside[rows]) for rows in peptide_rows)
    bins = NamedTuple[]
    affected_fraction = 0.0
    for (axis, coordinate) in ((:mz, s -> Float64(s.theoretical_mz)),
                              (:intensity, s -> log2(Float64(s.intensity))),
                              (:rt, s -> Float64(s.rt)))
        coordinates = [median(coordinate(samples[i]) for i in rows) for rows in peptide_rows]
        nbins = min(10, ng ÷ MS2_QC_MIN_PEPTIDES)
        nbins == 0 && continue
        # A constant coordinate has no local trend; global centering still applies.
        minimum(coordinates) == maximum(coordinates) && continue
        order = sortperm(coordinates)
        chunks = [order[(fld((b-1)*ng,nbins)+1):fld(b*ng,nbins)] for b in 1:nbins]
        stats = [_qc_supported_median(peptide_residuals[chunk]) for chunk in chunks]
        directions = _qc_bias_direction.(stats)
        persistent = [directions[b] != 0 &&
            ((b > 1 && directions[b-1] == directions[b]) ||
             (b < nbins && directions[b+1] == directions[b])) for b in 1:nbins]
        fraction = sum(length(chunks[b]) for b in 1:nbins if persistent[b]; init=0) / ng
        affected_fraction = max(affected_fraction, fraction)
        for b in 1:nbins
            chunk = chunks[b]
            before = mean(mean(outside[peptide_rows[g]]) for g in chunk)
            after = mean(mean(abs.(residuals[peptide_rows[g]] .- stats[b].median) .> 1) for g in chunk)
            push!(bins, (axis=String(axis), bin=b, n_peptides=length(chunk),
                coordinate_min=minimum(coordinates[chunk]), coordinate_max=maximum(coordinates[chunk]),
                median_error_mda=median([median(errors_mda[peptide_rows[g]]) for g in chunk]),
                median_error_over_tolerance=stats[b].median,
                median_ci_low=stats[b].low, median_ci_high=stats[b].high,
                persistent_bias=persistent[b], outside_tolerance_fraction=before,
                diagnostic_recentered_outside_fraction=after,
                diagnostic_net_coverage_gain=before-after))
        end
    end
    record = assess_calibration_qc(:ms2_mass,ng,
        (NaN,affected_fraction,supported_global_bias,outside_fraction);
        min_support=MS2_QC_MIN_PEPTIDES, fallback)
    return (; record, bins, n_peptides=ng,
        global_median_error_over_tolerance=global_bias.median,
        global_median_error_mda=ng == 0 ? NaN : median([median(errors_mda[rows]) for rows in peptide_rows]))
end

# Compatibility entry point for simple models and synthetic diagnostics. Production
# uses peptide-held-out models below, rather than scoring the installed fit's training data.
function ms2_calibration_qc(samples, model; fallback=false, group_ids=nothing)
    samples === nothing && (samples = MassErrSample[])
    ids = group_ids === nothing ? collect(eachindex(samples)) : group_ids
    return ms2_calibration_diagnostics(samples, fill(model,length(samples)); group_ids=ids, fallback).record
end

"""Fit five diagnostic models with whole peptide sequences held out; never install them."""
function heldout_ms2_calibration_qc(samples, model, sequences; fallback=false)
    if samples === nothing || isempty(samples)
        return (record=ms2_calibration_qc(samples,model; fallback=true), bins=NamedTuple[], n_peptides=0)
    end
    if any(s -> !(1 <= s.precursor_idx <= length(sequences)), samples)
        return (record=assess_calibration_qc(:ms2_mass,0,(NaN,NaN,NaN,NaN); min_support=50),
                bins=NamedTuple[], n_peptides=0)
    end
    ids = [String(sequences[s.precursor_idx]) for s in samples]
    unique_ids = sort!(unique(ids))
    if length(unique_ids) < 2 * MS2_QC_MIN_PEPTIDES || !(model isa IntensityMassErrorModel)
        # Do not label an in-sample or unsupported diagnostic as independent QC.
        return (record=assess_calibration_qc(:ms2_mass,length(unique_ids),(NaN,NaN,NaN,NaN);
                    min_support=2*MS2_QC_MIN_PEPTIDES, fallback=true), bins=NamedTuple[], n_peptides=length(unique_ids))
    end
    fold_by_id = Dict(id => mod1(i,5) for (i,id) in enumerate(unique_ids))
    folds = [fold_by_id[id] for id in ids]
    models = Vector{AbstractMassErrorModel}(undef,length(samples))
    for fold in 1:5
        train = samples[folds .!= fold]
        diagnostic_model = fit_intensity_mass_error_model(train, model; k=Float32(model.k))
        diagnostic_model isa IntensityMassErrorModel ||
            return (record=assess_calibration_qc(:ms2_mass,length(unique_ids),(NaN,NaN,NaN,NaN); fallback=true),
                    bins=NamedTuple[], n_peptides=length(unique_ids))
        models[folds .== fold] .= Ref(diagnostic_model)
    end
    return ms2_calibration_diagnostics(samples, models; group_ids=ids, fallback)
end

_qc_bin(x, low, high) = high > low ? clamp(floor(Int, 10 * (x-low)/(high-low)) + 1, 1, 10) : 1
function _qc_binned_bias(xs, residuals)
    low, high = extrema(xs)
    sums = zeros(Float64, 10)
    counts = zeros(Int, 10)
    for i in eachindex(xs, residuals)
        bin = _qc_bin(xs[i], low, high)
        sums[bin] += residuals[i]; counts[bin] += 1
    end
    bias = 0.0
    for i in eachindex(sums)
        counts[i] >= 20 && (bias = max(bias, abs(sums[i]/counts[i])))
    end
    return bias
end

function add_calibration_qc_columns!(table, state::CalibrationQCState, n)
    for stage in CALIBRATION_QC_STAGES
        records = [calibration_qc_record(state, stage, i) for i in 1:n]
        # Ion mobility is assessed only on ion-mobility data; other searches keep their columns as they were.
        stage === :ion_mobility && all(r -> r.status == QC_NOT_ASSESSED, records) && continue
        table[!, Symbol(stage, "_calibration_qc")] = [_qc_name(r.status) for r in records]
        table[!, Symbol(stage, "_calibration_qc_reason")] = [calibration_qc_reason(stage, r) for r in records]
    end
    return table
end
