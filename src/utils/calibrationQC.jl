"""Per-file calibration QC. Statuses and reason flags are stored without strings."""
@enum CalibrationQCStatus::UInt8 QC_NOT_ASSESSED QC_NORMAL QC_SUSPICIOUS QC_WARNING QC_FAILED
const CALIBRATION_QC_STAGES = (:rt, :ms2_mass, :ms1_mass, :quadrupole, :nce)
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
    records::NTuple{5,Vector{CalibrationQCRecord}}
    selected::NTuple{5,Set{Int}}
    pages::NTuple{5,Vector{String}}
end
CalibrationQCState() = CalibrationQCState(ntuple(_ -> CalibrationQCRecord[], 5),
    ntuple(_ -> Set{Int}(), 5), ntuple(_ -> String[], 5))
_qc_stage(stage::Symbol) = something(findfirst(==(stage), CALIBRATION_QC_STAGES))
_qc_name(s::CalibrationQCStatus) = ("not_assessed", "normal", "suspicious", "warning", "failed")[Int(s)+1]

# Conservative within-file screening heuristics, not fit acceptance criteria.
# Coverage is a fraction; other units are encoded by the metric names below.
const CALIBRATION_QC_METRICS = (
    (:rt_coverage, :residual_mad_fraction, :residual_trend_fraction, :unused),
    (:unused, :residual_rms_over_tolerance, :max_binned_bias_over_tolerance, :outside_tolerance_fraction),
    (:unused, :residual_mad_ppm, :residual_trend_ppm, :unused),
    (:edge_coverage, :log_ratio_rmse, :unused, :parameter_at_bound),
    (:fitted_charge_coverage, :unused, :weak_nce_fraction, :endpoint_fraction),
)
const CALIBRATION_QC_LIMITS = (
    (0.5f0, 0.05f0, 0.01f0, Inf32),
    (-Inf32, 1.0f0, 0.25f0, 0.10f0),
    (-Inf32, 10.0f0, 5.0f0, Inf32),
    (0.5f0, 0.5f0, Inf32, 0.5f0),
    (0.5f0, Inf32, 0.5f0, 0.5f0),
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

function ms2_calibration_qc(samples, model; fallback=false)
    n = samples === nothing ? 0 : length(samples)
    n == 0 && return assess_calibration_qc(:ms2_mass,0,(NaN,NaN,NaN,NaN); min_support=50, fallback=true)
    total = square = 0.0
    outside = 0
    mz_min, mz_max = extrema(s.theoretical_mz for s in samples)
    intensity_min, intensity_max = extrema(log2(max(s.intensity, eps(Float32))) for s in samples)
    sums = zeros(Float64, 10, 2)
    counts = zeros(Int, 10, 2)
    for s in samples
        corrected, low, high = getCorrectedMzAndBounds(model, s.observed_mz, s.intensity, s.rt)
        width = (high-low)/2
        if !isfinite(width) || width <= 0 || !isfinite(corrected)
            return assess_calibration_qc(:ms2_mass,n,(NaN,NaN,NaN,NaN); failed=true)
        end
        residual = (corrected-s.theoretical_mz)/width
        isfinite(residual) || return assess_calibration_qc(:ms2_mass,n,(NaN,NaN,NaN,NaN); failed=true)
        for (axis, x, low, high) in ((1,s.theoretical_mz,mz_min,mz_max),
                                     (2,log2(max(s.intensity,eps(Float32))),intensity_min,intensity_max))
            bin = _qc_bin(x, low, high)
            sums[bin,axis] += residual; counts[bin,axis] += 1
        end
        total += residual; square += residual^2
        outside += !(low <= s.theoretical_mz <= high)
    end
    bias = abs(total/n)
    for i in eachindex(sums)
        counts[i] >= 50 && (bias = max(bias, abs(sums[i]/counts[i])))
    end
    return assess_calibration_qc(:ms2_mass,n,(NaN,sqrt(square/n),bias,outside/n);
        min_support=50, fallback)
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
        table[!, Symbol(stage, "_calibration_qc")] = [_qc_name(r.status) for r in records]
        table[!, Symbol(stage, "_calibration_qc_reason")] = [calibration_qc_reason(stage, r) for r in records]
    end
    return table
end
