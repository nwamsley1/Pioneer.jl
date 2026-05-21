# Chromatogram-boundary shape-model scoring utilities.
const BOUNDARY_SHAPE_SIGMA_IRT_MAX_BINS = 24
const BOUNDARY_SHAPE_SIGMA_IRT_MIN_BIN_ROWS = 25
const BOUNDARY_SHAPE_FIT_PRIOR_SCORE_WEIGHT = 1.0f0
const BOUNDARY_SHAPE_FIT_PRIOR_NORMALIZATION_FLOOR = 0.02f0

# Temporary hardcoded override: set to false to re-enable shape-model bounds.
const FORCE_FALLBACK_CHROMATOGRAM_BOUNDS = true

const DEBUG_BOUNDARY_CANDIDATE_TARGETS = Ref{Set{Tuple{UInt32,UInt16}}}(
    Set([
        (UInt32(370714), UInt16(1)),
        (UInt32(909488), UInt16(3)),
        (UInt32(1047223), UInt16(3)),
    ])
)

@inline function egh_boundary_value(
    t::Real,
    apex_t::Real,
    height::Real,
    sigma::Real,
    tau::Real,
)
    height_f = Float32(height)
    sigma_f = Float32(sigma)
    sigma_f > 0.0f0 || return 0.0f0

    delta = Float32(t - apex_t)
    denom = 2.0f0 * sigma_f * sigma_f + Float32(tau) * delta
    denom > eps(Float32) || return 0.0f0
    val = height_f * exp(-(delta * delta) / denom)
    return isfinite(val) ? Float32(max(0.0f0, val)) : 0.0f0
end

function _shape_point_widths(x::AbstractVector{<:Real})
    n = length(x)
    widths = Vector{Float32}(undef, n)
    if n == 1
        widths[1] = 1.0f0
        return widths
    end

    @inbounds for i in 1:n
        left = i == 1 ? 0.0f0 : abs(Float32(x[i] - x[i - 1]))
        right = i == n ? 0.0f0 : abs(Float32(x[i + 1] - x[i]))
        widths[i] = max(0.5f0 * (left + right), eps(Float32))
    end
    return widths
end

function _fit_egh_linearized(
    irt::AbstractVector{<:Real},
    signal::AbstractVector{<:Real},
    apex_pos::Int,
    height::Float32,
    widths::Vector{Float32},
)
    floor_ratio = 1.0f-4
    ceiling_ratio = 0.98f0
    apex_t = Float32(irt[apex_pos])

    s00 = 0.0f0
    s01 = 0.0f0
    s11 = 0.0f0
    b0 = 0.0f0
    b1 = 0.0f0

    @inbounds for i in eachindex(signal)
        i == apex_pos && continue
        delta = Float32(irt[i]) - apex_t
        abs(delta) > eps(Float32) || continue
        ratio = clamp(Float32(signal[i]) / height, floor_ratio, ceiling_ratio)
        denom_obs = (delta * delta) / max(-log(ratio), eps(Float32))
        weight = widths[i] * max(0.25f0, sqrt(max(ratio, floor_ratio)))
        s00 += weight
        s01 += weight * delta
        s11 += weight * delta * delta
        b0 += weight * denom_obs
        b1 += weight * delta * denom_obs
    end

    det = s00 * s11 - s01 * s01
    if det <= eps(Float32)
        return nothing
    end

    intercept = (b0 * s11 - b1 * s01) / det
    tau = (s00 * b1 - s01 * b0) / det
    if !isfinite(intercept) || intercept <= eps(Float32) || !isfinite(tau)
        return nothing
    end

    sigma = sqrt(0.5f0 * intercept)
    sigma > eps(Float32) || return nothing
    return (sigma = Float32(sigma), tau = Float32(tau))
end

function _fit_gaussian_moment_fallback(
    irt::AbstractVector{<:Real},
    signal::AbstractVector{<:Real},
    apex_pos::Int,
    height::Float32,
    widths::Vector{Float32},
)
    apex_t = Float32(irt[apex_pos])
    moment = 0.0f0
    total = 0.0f0

    @inbounds for i in eachindex(signal)
        val = max(0.0f0, Float32(signal[i]))
        weight = widths[i] * max(val / height, 0.05f0)
        delta = Float32(irt[i]) - apex_t
        moment += weight * delta * delta
        total += weight
    end

    total > eps(Float32) || return nothing
    sigma = sqrt(moment / total)
    sigma > eps(Float32) || return nothing
    return (sigma = Float32(sigma), tau = 0.0f0)
end

function _egh_fit_loss(
    irt::AbstractVector{<:Real},
    signal::AbstractVector{<:Real},
    apex_pos::Int,
    height::Float32,
    widths::Vector{Float32},
    sigma::Float32,
    tau::Float32,
)
    apex_t = Float32(irt[apex_pos])
    total_weight = 0.0f0
    fit_loss = 0.0f0
    overshoot_penalty = 0.0f0

    @inbounds for i in eachindex(signal)
        observed = max(0.0f0, Float32(signal[i]))
        predicted = egh_boundary_value(irt[i], apex_t, height, sigma, tau)
        residual = (observed - predicted) / height
        point_loss = residual * residual
        weight = widths[i]
        fit_loss += weight * point_loss
        total_weight += weight

        overshoot = max(0.0f0, observed / height - 1.0f0)
        overshoot_penalty += weight * overshoot * overshoot
    end

    total_weight = max(total_weight, eps(Float32))
    edge_fraction = 0.5f0 * (
        egh_boundary_value(first(irt), apex_t, height, sigma, tau) +
        egh_boundary_value(last(irt), apex_t, height, sigma, tau)
    ) / height

    return (
        fit_loss = Float32(fit_loss / total_weight),
        edge_fraction = Float32(edge_fraction),
        overshoot_penalty = Float32(overshoot_penalty / total_weight),
    )
end

function fit_egh_boundary_shape(
    irt::AbstractVector{<:Real},
    signal::AbstractVector{<:Real},
    apex_pos::Integer,
)
    n = length(irt)
    if n != length(signal) || n < 3 || apex_pos < 1 || apex_pos > n
        return (
            valid = false,
            fit_loss = Inf32,
            sigma = 0.0f0,
            tau = 0.0f0,
            edge_fraction = 1.0f0,
            overshoot_penalty = Inf32,
        )
    end

    apex_i = Int(apex_pos)
    height = Float32(signal[apex_i])
    if !isfinite(height) || height <= eps(Float32)
        return (
            valid = false,
            fit_loss = Inf32,
            sigma = 0.0f0,
            tau = 0.0f0,
            edge_fraction = 1.0f0,
            overshoot_penalty = Inf32,
        )
    end

    widths = _shape_point_widths(irt)
    fit = _fit_egh_linearized(irt, signal, apex_i, height, widths)
    fit === nothing && (fit = _fit_gaussian_moment_fallback(irt, signal, apex_i, height, widths))
    fit === nothing && return (
        valid = false,
        fit_loss = Inf32,
        sigma = 0.0f0,
        tau = 0.0f0,
        edge_fraction = 1.0f0,
        overshoot_penalty = Inf32,
    )

    loss = _egh_fit_loss(irt, signal, apex_i, height, widths, fit.sigma, fit.tau)
    return (
        valid = true,
        fit_loss = loss.fit_loss,
        sigma = fit.sigma,
        tau = fit.tau,
        edge_fraction = loss.edge_fraction,
        overshoot_penalty = loss.overshoot_penalty,
    )
end

function _get_column_or_fill(df::AbstractDataFrame, col::Symbol, default)
    hasproperty(df, col) && return df[!, col]
    return fill(default, nrow(df))
end

function _fallback_rows_by_group(
    candidates::AbstractDataFrame,
    group_col::Symbol,
    fallback_col::Symbol,
)
    groups = candidates[!, group_col]
    fallback = _get_column_or_fill(candidates, fallback_col, false)
    rows_by_group = Dict{UInt64, Vector{Int}}()

    for i in 1:nrow(candidates)
        group_id = UInt64(groups[i])
        push!(get!(() -> Int[], rows_by_group, group_id), i)
    end

    fallback_row = Dict{UInt64, Int}()
    for (group_id, rows) in rows_by_group
        chosen = rows[1]
        for row in rows
            if Bool(fallback[row])
                chosen = row
                break
            end
        end
        fallback_row[group_id] = chosen
    end

    return rows_by_group, fallback_row
end

function _median_mad_scale(values::Vector{Float32}; floor::Float32)
    isempty(values) && return (median = 0.0f0, scale = floor)
    med = Float32(median(values))
    deviations = Float32[abs(v - med) for v in values]
    mad = isempty(deviations) ? 0.0f0 : Float32(median(deviations))
    scale = max(1.4826f0 * mad, floor)
    return (median = med, scale = scale)
end

function _smooth_curve_values(values::Vector{Float32})
    length(values) < 3 && return copy(values)
    smoothed = similar(values)
    @inbounds for i in eachindex(values)
        start = max(firstindex(values), i - 1)
        stop = min(lastindex(values), i + 1)
        smoothed[i] = Float32(median(@view values[start:stop]))
    end
    return smoothed
end

function _shape_sigma_irt_curve(
    apex_irt::Vector{Float32},
    log_sigma::Vector{Float32};
    max_bins::Integer = BOUNDARY_SHAPE_SIGMA_IRT_MAX_BINS,
    min_bin_rows::Integer = BOUNDARY_SHAPE_SIGMA_IRT_MIN_BIN_ROWS,
)
    n = length(log_sigma)
    if n == 0 || length(apex_irt) != n
        return (irt = Float32[], log_sigma = Float32[])
    end

    if n < 2 * min_bin_rows
        return (
            irt = Float32[Float32(median(apex_irt))],
            log_sigma = Float32[Float32(median(log_sigma))],
        )
    end

    order = sortperm(apex_irt)
    n_bins = max(2, min(Int(max_bins), n ÷ Int(min_bin_rows)))
    knots_irt = Float32[]
    knots_log_sigma = Float32[]
    sizehint!(knots_irt, n_bins)
    sizehint!(knots_log_sigma, n_bins)

    for bin in 1:n_bins
        start = floor(Int, (bin - 1) * n / n_bins) + 1
        stop = floor(Int, bin * n / n_bins)
        stop >= start || continue
        idxs = @view order[start:stop]
        bin_irt = Float32[apex_irt[idx] for idx in idxs]
        bin_log_sigma = Float32[log_sigma[idx] for idx in idxs]
        push!(knots_irt, Float32(median(bin_irt)))
        push!(knots_log_sigma, Float32(median(bin_log_sigma)))
    end

    smoothed_log_sigma = _smooth_curve_values(knots_log_sigma)
    keep_irt = Float32[]
    keep_log_sigma = Float32[]
    for i in eachindex(knots_irt)
        if isempty(keep_irt) || knots_irt[i] > last(keep_irt) + eps(Float32)
            push!(keep_irt, knots_irt[i])
            push!(keep_log_sigma, smoothed_log_sigma[i])
        else
            keep_log_sigma[end] = 0.5f0 * (keep_log_sigma[end] + smoothed_log_sigma[i])
        end
    end

    return (irt = keep_irt, log_sigma = keep_log_sigma)
end

function _interpolate_shape_log_sigma(
    apex_irt::Float32,
    knots_irt::Vector{Float32},
    knots_log_sigma::Vector{Float32},
    fallback::Float32,
)
    n = min(length(knots_irt), length(knots_log_sigma))
    n == 0 && return fallback
    n == 1 && return knots_log_sigma[1]
    apex_irt <= knots_irt[1] && return knots_log_sigma[1]
    apex_irt >= knots_irt[n] && return knots_log_sigma[n]

    right = searchsortedfirst(knots_irt, apex_irt)
    left = max(1, right - 1)
    denom = max(knots_irt[right] - knots_irt[left], eps(Float32))
    frac = clamp((apex_irt - knots_irt[left]) / denom, 0.0f0, 1.0f0)
    return Float32((1.0f0 - frac) * knots_log_sigma[left] + frac * knots_log_sigma[right])
end

function _shape_sigma_curve_inputs(candidates::AbstractDataFrame, rows::Vector{Int})
    apex_irt = Float32[]
    log_sigma = Float32[]
    hasproperty(candidates, :shape_apex_irt) || return apex_irt, log_sigma
    sizehint!(apex_irt, length(rows))
    sizehint!(log_sigma, length(rows))

    for row in rows
        sigma = Float32(candidates[row, :shape_sigma_irt])
        irt = Float32(candidates[row, :shape_apex_irt])
        if isfinite(sigma) && sigma > eps(Float32) && isfinite(irt)
            push!(apex_irt, irt)
            push!(log_sigma, log2(sigma))
        end
    end

    return apex_irt, log_sigma
end

function _shape_sigma_prior(candidates::AbstractDataFrame, rows::Vector{Int}, fallback)
    log_sigma = Float32[
        log2(Float32(candidates[row, :shape_sigma_irt])) for row in rows
    ]
    base = isempty(log_sigma) ? fallback : _median_mad_scale(log_sigma; floor = 0.5f0)
    apex_irt, curve_log_sigma = _shape_sigma_curve_inputs(candidates, rows)
    if isempty(curve_log_sigma)
        return (
            median = base.median,
            scale = base.scale,
            curve_irt = Float32[],
            curve_log_sigma = Float32[],
        )
    end

    curve = _shape_sigma_irt_curve(apex_irt, curve_log_sigma)
    residuals = Float32[
        curve_log_sigma[i] - _interpolate_shape_log_sigma(
            apex_irt[i],
            curve.irt,
            curve.log_sigma,
            base.median,
        )
        for i in eachindex(curve_log_sigma)
    ]
    residual_scale = _median_mad_scale(residuals; floor = 0.5f0)
    return (
        median = base.median,
        scale = residual_scale.scale,
        curve_irt = curve.irt,
        curve_log_sigma = curve.log_sigma,
    )
end

function boundary_shape_priors(candidates::AbstractDataFrame)
    valid_cols = all(hasproperty(candidates, col) for col in (
        :shape_model_valid,
        :shape_sigma_irt,
        :shape_tau_irt,
        :ms_file_idx,
    ))
    valid_cols || return Dict{UInt16, NamedTuple}()

    fallback = _get_column_or_fill(candidates, :is_fallback, false)
    primary_apex = hasproperty(candidates, :is_primary_apex) ?
        candidates[!, :is_primary_apex] :
        fill(true, nrow(candidates))
    valid = candidates[!, :shape_model_valid]
    files = sort(unique(UInt16.(candidates[!, :ms_file_idx])))
    priors = Dict{UInt16, NamedTuple}()

    all_rows = [
        i for i in 1:nrow(candidates)
        if Bool(fallback[i]) && Bool(primary_apex[i]) && Bool(valid[i]) &&
           isfinite(Float32(candidates[i, :shape_sigma_irt])) &&
           Float32(candidates[i, :shape_sigma_irt]) > eps(Float32)
    ]
    if isempty(all_rows)
        all_rows = [
            i for i in 1:nrow(candidates)
            if Bool(fallback[i]) && Bool(valid[i]) &&
               isfinite(Float32(candidates[i, :shape_sigma_irt])) &&
               Float32(candidates[i, :shape_sigma_irt]) > eps(Float32)
        ]
    end
    isempty(all_rows) && return priors

    global_log_sigma = Float32[
        log2(Float32(candidates[row, :shape_sigma_irt])) for row in all_rows
    ]
    global_tau_ratio = Float32[
        Float32(candidates[row, :shape_tau_irt]) /
            max(Float32(candidates[row, :shape_sigma_irt]), eps(Float32))
        for row in all_rows
    ]
    global_sigma_prior = _median_mad_scale(global_log_sigma; floor = 0.5f0)
    global_sigma_curve_prior = _shape_sigma_prior(candidates, all_rows, global_sigma_prior)
    global_tau_prior = _median_mad_scale(global_tau_ratio; floor = 1.0f0)

    for file_idx in files
        rows = [
            i for i in all_rows
            if UInt16(candidates[i, :ms_file_idx]) == file_idx
        ]
        sigma_prior = isempty(rows) ? global_sigma_curve_prior :
            _shape_sigma_prior(candidates, rows, global_sigma_prior)
        tau_ratio = isempty(rows) ? global_tau_ratio : Float32[
            Float32(candidates[row, :shape_tau_irt]) /
                max(Float32(candidates[row, :shape_sigma_irt]), eps(Float32))
            for row in rows
        ]
        tau_prior = isempty(tau_ratio) ? global_tau_prior :
            _median_mad_scale(tau_ratio; floor = 1.0f0)
        priors[file_idx] = (
            log_sigma_median = sigma_prior.median,
            log_sigma_scale = sigma_prior.scale,
            sigma_irt_knots = sigma_prior.curve_irt,
            log_sigma_knots = sigma_prior.curve_log_sigma,
            tau_ratio_median = tau_prior.median,
            tau_ratio_scale = tau_prior.scale,
        )
    end

    return priors
end

function _boundary_shape_prior_penalty(candidates::AbstractDataFrame, row::Int, priors)
    hasproperty(candidates, :ms_file_idx) || return 0.0f0
    hasproperty(candidates, :shape_sigma_irt) || return 0.0f0
    hasproperty(candidates, :shape_tau_irt) || return 0.0f0

    sigma = Float32(candidates[row, :shape_sigma_irt])
    isfinite(sigma) && sigma > eps(Float32) || return 0.0f0
    file_idx = UInt16(candidates[row, :ms_file_idx])
    haskey(priors, file_idx) || return 0.0f0
    prior = priors[file_idx]

    log_sigma = log2(sigma)
    apex_irt = hasproperty(candidates, :shape_apex_irt) ?
        Float32(candidates[row, :shape_apex_irt]) :
        NaN32
    expected_log_sigma = isfinite(apex_irt) ?
        _interpolate_shape_log_sigma(
            apex_irt,
            prior.sigma_irt_knots,
            prior.log_sigma_knots,
            prior.log_sigma_median,
        ) :
        prior.log_sigma_median
    tau_ratio = Float32(candidates[row, :shape_tau_irt]) / sigma
    z_sigma = clamp(
        (log_sigma - expected_log_sigma) / prior.log_sigma_scale,
        -4.0f0,
        4.0f0,
    )
    z_tau = clamp(
        (tau_ratio - prior.tau_ratio_median) / prior.tau_ratio_scale,
        -4.0f0,
        4.0f0,
    )
    return Float32(0.05f0 * (z_sigma * z_sigma + z_tau * z_tau))
end

function _boundary_shape_unexplained_area_penalty(candidates::AbstractDataFrame, row::Int)
    hasproperty(candidates, :shape_deconvolution_area_fraction) || return 0.0f0
    fraction = Float32(candidates[row, :shape_deconvolution_area_fraction])
    isfinite(fraction) || return 0.0f0
    explained_fraction = clamp(fraction, 0.0f0, 1.0f0)
    return 1.0f0 - explained_fraction
end

function score_boundary_candidates_by_shape!(
    candidates::DataFrame;
    score_col::Symbol = :boundary_model_score,
)
    n = nrow(candidates)
    candidates[!, :boundary_shape_score] = fill(Inf32, n)
    candidates[!, :shape_prior_penalty] = zeros(Float32, n)
    candidates[!, :shape_fit_prior_score] = zeros(Float32, n)
    candidates[!, :shape_normalized_fit_prior_score] = zeros(Float32, n)
    candidates[!, :shape_unexplained_area_penalty] = zeros(Float32, n)
    candidates[!, score_col] = fill(-Inf32, n)
    n == 0 && return candidates

    required = (
        :shape_model_valid,
        :shape_fit_loss,
        :shape_sigma_irt,
        :shape_tau_irt,
    )
    if !all(hasproperty(candidates, col) for col in required)
        fallback = _get_column_or_fill(candidates, :is_fallback, false)
        candidates[!, score_col] = Float32[Bool(x) ? 1.0f0 : 0.0f0 for x in fallback]
        candidates[!, :boundary_shape_score] = Float32[Bool(x) ? 0.0f0 : Inf32 for x in fallback]
        return candidates
    end

    priors = boundary_shape_priors(candidates)
    valid_rows = Int[]
    best_fit_prior_by_group = Dict{UInt64, Tuple{Float32,Float32}}()
    @inbounds for row in 1:n
        if !Bool(candidates[row, :shape_model_valid])
            continue
        end

        fit_loss = Float32(candidates[row, :shape_fit_loss])
        if !isfinite(fit_loss)
            continue
        end

        prior_penalty = _boundary_shape_prior_penalty(candidates, row, priors)
        fit_prior_score = fit_loss + prior_penalty
        group_id = hasproperty(candidates, :boundary_group_id) ?
            UInt64(candidates[row, :boundary_group_id]) :
            UInt64(0)
        best_fit_prior, second_best_fit_prior = get(
            best_fit_prior_by_group,
            group_id,
            (Inf32, Inf32),
        )
        if fit_prior_score <= best_fit_prior
            second_best_fit_prior = best_fit_prior
            best_fit_prior = fit_prior_score
        elseif fit_prior_score < second_best_fit_prior
            second_best_fit_prior = fit_prior_score
        end
        best_fit_prior_by_group[group_id] = (best_fit_prior, second_best_fit_prior)

        candidates[row, :shape_prior_penalty] = prior_penalty
        candidates[row, :shape_fit_prior_score] = fit_prior_score
        candidates[row, :shape_unexplained_area_penalty] =
            _boundary_shape_unexplained_area_penalty(candidates, row)
        push!(valid_rows, row)
    end

    @inbounds for row in valid_rows
        group_id = hasproperty(candidates, :boundary_group_id) ?
            UInt64(candidates[row, :boundary_group_id]) :
            UInt64(0)
        best_fit_prior, second_best_fit_prior = get(
            best_fit_prior_by_group,
            group_id,
            (0.0f0, Inf32),
        )
        normalization_scale = isfinite(second_best_fit_prior) ?
            max(
                second_best_fit_prior - best_fit_prior,
                BOUNDARY_SHAPE_FIT_PRIOR_NORMALIZATION_FLOOR,
            ) :
            BOUNDARY_SHAPE_FIT_PRIOR_NORMALIZATION_FLOOR
        fit_prior_delta = max(
            Float32(candidates[row, :shape_fit_prior_score]) - best_fit_prior,
            0.0f0,
        )
        normalized_fit_prior = clamp(fit_prior_delta / normalization_scale, 0.0f0, 1.0f0)
        shape_score = BOUNDARY_SHAPE_FIT_PRIOR_SCORE_WEIGHT * normalized_fit_prior +
            Float32(candidates[row, :shape_unexplained_area_penalty])

        candidates[row, :shape_normalized_fit_prior_score] = normalized_fit_prior
        candidates[row, :boundary_shape_score] = Float32(shape_score)
        candidates[row, score_col] = -Float32(shape_score)
    end

    return candidates
end

function _boundary_shape_qc_rows(candidates::AbstractDataFrame)
    required = (
        :ms_file_idx,
        :is_fallback,
        :shape_model_valid,
        :shape_sigma_irt,
        :shape_tau_irt,
        :shape_apex_irt,
        :shape_apex_height,
    )
    all(hasproperty(candidates, col) for col in required) || return Int[]

    primary_apex = hasproperty(candidates, :is_primary_apex) ?
        candidates[!, :is_primary_apex] :
        fill(true, nrow(candidates))
    rows = Int[]
    sizehint!(rows, nrow(candidates))
    @inbounds for row in 1:nrow(candidates)
        sigma = Float32(candidates[row, :shape_sigma_irt])
        tau = Float32(candidates[row, :shape_tau_irt])
        apex_irt = Float32(candidates[row, :shape_apex_irt])
        apex_height = Float32(candidates[row, :shape_apex_height])
        if Bool(candidates[row, :is_fallback]) &&
           Bool(primary_apex[row]) &&
           Bool(candidates[row, :shape_model_valid]) &&
           isfinite(sigma) && sigma > eps(Float32) &&
           isfinite(tau) &&
           isfinite(apex_irt) &&
           isfinite(apex_height) && apex_height > eps(Float32)
            push!(rows, row)
        end
    end
    return rows
end

function _downsample_boundary_shape_qc_rows(rows::Vector{Int}, max_points::Integer)
    max_points <= 0 && return rows
    length(rows) <= max_points && return rows

    step = ceil(Int, length(rows) / max_points)
    sampled = rows[1:step:end]
    if length(sampled) > max_points
        sampled = sampled[1:Int(max_points)]
    end
    return collect(sampled)
end

function _boundary_shape_qc_plot(
    sigma::Vector{Float64},
    tau::Vector{Float64},
    tau_ratio::Vector{Float64},
    apex_irt::Vector{Float64},
    log2_apex_height::Vector{Float64},
    title_prefix::AbstractString;
    sigma_curve_irt::Vector{Float64} = Float64[],
    sigma_curve::Vector{Float64} = Float64[],
)
    n = length(sigma)
    bins = max(5, min(60, ceil(Int, sqrt(max(n, 1)))))
    scatter_alpha = n > 1000 ? 0.25 : 0.55

    p_sigma_hist = Plots.histogram(
        sigma;
        bins = bins,
        title = "$(title_prefix) sigma",
        xlabel = "sigma (iRT)",
        ylabel = "count",
        legend = false,
        color = :dodgerblue,
        alpha = 0.75,
    )
    p_tau_hist = Plots.histogram(
        tau;
        bins = bins,
        title = "$(title_prefix) tau",
        xlabel = "tau (iRT)",
        ylabel = "count",
        legend = false,
        color = :darkorange,
        alpha = 0.75,
    )
    p_ratio_hist = Plots.histogram(
        tau_ratio;
        bins = bins,
        title = "$(title_prefix) tau / sigma",
        xlabel = "tau / sigma",
        ylabel = "count",
        legend = false,
        color = :seagreen,
        alpha = 0.75,
    )

    p_sigma_irt = Plots.scatter(
        apex_irt,
        sigma;
        title = "sigma vs apex iRT",
        xlabel = "apex iRT",
        ylabel = "sigma (iRT)",
        legend = false,
        alpha = scatter_alpha,
        markersize = 2,
        markerstrokewidth = 0,
        color = :dodgerblue,
    )
    if !isempty(sigma_curve_irt) && length(sigma_curve_irt) == length(sigma_curve)
        Plots.plot!(
            p_sigma_irt,
            sigma_curve_irt,
            sigma_curve;
            label = "sigma prior",
            legend = :topright,
            color = :black,
            linewidth = 3,
        )
    end
    p_tau_irt = Plots.scatter(
        apex_irt,
        tau;
        title = "tau vs apex iRT",
        xlabel = "apex iRT",
        ylabel = "tau (iRT)",
        legend = false,
        alpha = scatter_alpha,
        markersize = 2,
        markerstrokewidth = 0,
        color = :darkorange,
    )
    p_ratio_irt = Plots.scatter(
        apex_irt,
        tau_ratio;
        title = "tau / sigma vs apex iRT",
        xlabel = "apex iRT",
        ylabel = "tau / sigma",
        legend = false,
        alpha = scatter_alpha,
        markersize = 2,
        markerstrokewidth = 0,
        color = :seagreen,
    )

    p_sigma_intensity = Plots.scatter(
        log2_apex_height,
        sigma;
        title = "sigma vs log2 apex height",
        xlabel = "log2 apex height",
        ylabel = "sigma (iRT)",
        legend = false,
        alpha = scatter_alpha,
        markersize = 2,
        markerstrokewidth = 0,
        color = :dodgerblue,
    )
    p_tau_intensity = Plots.scatter(
        log2_apex_height,
        tau;
        title = "tau vs log2 apex height",
        xlabel = "log2 apex height",
        ylabel = "tau (iRT)",
        legend = false,
        alpha = scatter_alpha,
        markersize = 2,
        markerstrokewidth = 0,
        color = :darkorange,
    )
    p_ratio_intensity = Plots.scatter(
        log2_apex_height,
        tau_ratio;
        title = "tau / sigma vs log2 apex height",
        xlabel = "log2 apex height",
        ylabel = "tau / sigma",
        legend = false,
        alpha = scatter_alpha,
        markersize = 2,
        markerstrokewidth = 0,
        color = :seagreen,
    )

    return Plots.plot(
        p_sigma_hist,
        p_tau_hist,
        p_ratio_hist,
        p_sigma_irt,
        p_tau_irt,
        p_ratio_irt,
        p_sigma_intensity,
        p_tau_intensity,
        p_ratio_intensity;
        layout = (3, 3),
        size = (1500, 1200),
        dpi = 150,
    )
end

function write_boundary_shape_parameter_qc_plots(
    candidates::AbstractDataFrame,
    output_root::AbstractString;
    file_name_by_idx::AbstractDict = Dict{UInt16, String}(),
    max_points_per_file::Integer = 5000,
)
    rows = _boundary_shape_qc_rows(candidates)
    isempty(rows) && return String[]

    out_dir = joinpath(
        String(output_root),
        "qc_plots",
        "chromatogram_boundary_shape_parameters",
    )
    isdir(out_dir) || mkpath(out_dir)

    paths = String[]
    files = sort(unique(UInt16(candidates[row, :ms_file_idx]) for row in rows))
    priors = boundary_shape_priors(candidates)
    withenv("GKSwstype" => "100", "GKS_WSTYPE" => "100") do
        Plots.gr()

        for file_idx in files
            file_rows = [row for row in rows if UInt16(candidates[row, :ms_file_idx]) == file_idx]
            isempty(file_rows) && continue
            plot_rows = _downsample_boundary_shape_qc_rows(file_rows, max_points_per_file)

            sigma = Float64[Float64(candidates[row, :shape_sigma_irt]) for row in plot_rows]
            tau = Float64[Float64(candidates[row, :shape_tau_irt]) for row in plot_rows]
            tau_ratio = Float64[
                Float64(candidates[row, :shape_tau_irt]) /
                    max(Float64(candidates[row, :shape_sigma_irt]), eps(Float64))
                for row in plot_rows
            ]
            apex_irt = Float64[Float64(candidates[row, :shape_apex_irt]) for row in plot_rows]
            log2_apex_height = Float64[
                log2(max(Float64(candidates[row, :shape_apex_height]), eps(Float64)))
                for row in plot_rows
            ]
            prior = get(priors, file_idx, nothing)
            sigma_curve_irt = Float64[]
            sigma_curve = Float64[]
            if prior !== nothing && length(prior.sigma_irt_knots) > 0
                sigma_curve_irt = Float64.(prior.sigma_irt_knots)
                sigma_curve = Float64.(2.0f0 .^ prior.log_sigma_knots)
            end

            file_name = string(get(file_name_by_idx, file_idx, "run"))
            safe_name = debug_sanitize_chromatogram_filename(file_name)
            path = joinpath(
                out_dir,
                "egh_shape_parameters_file_$(Int(file_idx))_$(safe_name).png",
            )
            plot_obj = _boundary_shape_qc_plot(
                sigma,
                tau,
                tau_ratio,
                apex_irt,
                log2_apex_height,
                "file $(Int(file_idx))";
                sigma_curve_irt = sigma_curve_irt,
                sigma_curve = sigma_curve,
            )
            Plots.savefig(plot_obj, path)
            push!(paths, path)
        end
    end

    joined_paths = join(paths, ",")
    isempty(paths) || @debug_l1 "boundary_shape_parameter_qc_plots paths=$(joined_paths)"
    return paths
end

function write_boundary_shape_parameter_qc_plots(
    candidates::AbstractDataFrame,
    search_context::SearchContext;
    kwargs...,
)
    if !hasproperty(candidates, :ms_file_idx)
        return String[]
    end

    ms_data = getMSData(search_context)
    file_name_by_idx = Dict{UInt16, String}()
    for file_idx in sort(unique(UInt16.(candidates[!, :ms_file_idx])))
        file_name_by_idx[file_idx] = string(getFileIdToName(ms_data, Int(file_idx)))
    end
    return write_boundary_shape_parameter_qc_plots(
        candidates,
        getDataOutDir(search_context);
        file_name_by_idx = file_name_by_idx,
        kwargs...,
    )
end

function boundary_candidate_category_counts(selected_candidates::AbstractDataFrame)
    counts = Dict{String, Int}(label => 0 for label in BOUNDARY_CANDIDATE_CATEGORY_LABELS)
    hasproperty(selected_candidates, :candidate_category) || return counts

    for category in selected_candidates[!, :candidate_category]
        label = String(category)
        counts[label] = get(counts, label, 0) + 1
    end
    return counts
end

function boundary_candidate_category_tally_lines(selected_candidates::AbstractDataFrame)
    total = nrow(selected_candidates)
    counts = boundary_candidate_category_counts(selected_candidates)
    lines = ["Chromatogram boundary shape-model selected candidate categories (total $(total)):"]

    for label in BOUNDARY_CANDIDATE_CATEGORY_LABELS
        push!(lines, "    $(rpad(label, 24)) $(counts[label])")
    end

    extra_labels = sort(setdiff(collect(keys(counts)), collect(BOUNDARY_CANDIDATE_CATEGORY_LABELS)))
    for label in extra_labels
        push!(lines, "    $(rpad(label, 24)) $(counts[label])")
    end

    return lines
end

function log_boundary_candidate_category_tally(selected_candidates::AbstractDataFrame)
    lines = boundary_candidate_category_tally_lines(selected_candidates)
    @debug_l1 join(lines, "\n")
    return nothing
end

function quant_boundary_candidate_rows(candidates::AbstractDataFrame)
    hasproperty(candidates, :quant_trace_selected) || return copy(candidates)
    mask = Bool.(candidates[!, :quant_trace_selected])
    return candidates[mask, :]
end

function _boundary_debug_value(df::AbstractDataFrame, row::Integer, col::Symbol)
    hasproperty(df, col) || return "NA"
    value = df[Int(row), col]
    if value isa AbstractFloat
        isfinite(value) || return string(value)
        return string(round(Float64(value), sigdigits = 5))
    end
    return string(value)
end

function _boundary_debug_selected_keys(selected_candidates::AbstractDataFrame)
    keys = Set{Tuple{UInt64, UInt16}}()
    if !hasproperty(selected_candidates, :boundary_group_id) ||
       !hasproperty(selected_candidates, :candidate_index)
        return keys
    end

    for row in eachrow(selected_candidates)
        push!(keys, (UInt64(row.boundary_group_id), UInt16(row.candidate_index)))
    end
    return keys
end

function boundary_candidate_debug_lines(
    candidates::AbstractDataFrame,
    selected_candidates::AbstractDataFrame;
    target_precursor_idx::UInt32 = UInt32(370714),
    target_ms_file_idx::UInt16 = UInt16(1),
)
    hasproperty(candidates, :precursor_idx) ||
        return ["Boundary candidate debug precursor_idx=$(target_precursor_idx) ms_file_idx=$(target_ms_file_idx) rows=0 reason=missing_precursor_idx"]
    hasproperty(candidates, :ms_file_idx) ||
        return ["Boundary candidate debug precursor_idx=$(target_precursor_idx) ms_file_idx=$(target_ms_file_idx) rows=0 reason=missing_ms_file_idx"]

    rows = findall(
        i -> UInt32(candidates[i, :precursor_idx]) == target_precursor_idx &&
             UInt16(candidates[i, :ms_file_idx]) == target_ms_file_idx,
        axes(candidates, 1),
    )
    if hasproperty(candidates, :boundary_group_id) && hasproperty(candidates, :candidate_index)
        sort!(rows, by = i -> (UInt64(candidates[i, :boundary_group_id]), UInt16(candidates[i, :candidate_index])))
    end

    selected_keys = _boundary_debug_selected_keys(selected_candidates)
    lines = [
        "Boundary candidate debug precursor_idx=$(target_precursor_idx) " *
        "ms_file_idx=$(target_ms_file_idx) rows=$(length(rows))",
    ]

    for row in rows
        group_id = hasproperty(candidates, :boundary_group_id) ?
            UInt64(candidates[row, :boundary_group_id]) : UInt64(0)
        candidate_index = hasproperty(candidates, :candidate_index) ?
            UInt16(candidates[row, :candidate_index]) : UInt16(0)
        selected = (group_id, candidate_index) in selected_keys
        range_idx = _boundary_debug_value(candidates, row, :candidate_start_idx) *
            ":" * _boundary_debug_value(candidates, row, :candidate_stop_idx)
        range_scan = _boundary_debug_value(candidates, row, :candidate_start_scan) *
            ":" * _boundary_debug_value(candidates, row, :candidate_stop_scan)

        push!(lines, join((
            "  group=$(group_id)",
            "candidate_index=$(candidate_index)",
            "selected=$(selected)",
            "category=$(_boundary_debug_value(candidates, row, :candidate_category))",
            "range_idx=$(range_idx)",
            "range_scan=$(range_scan)",
            "apex_scan=$(_boundary_debug_value(candidates, row, :new_best_scan))",
            "area=$(_boundary_debug_value(candidates, row, :peak_area))",
            "points=$(_boundary_debug_value(candidates, row, :points_integrated))",
            "score=$(_boundary_debug_value(candidates, row, :boundary_model_score))",
            "shape_score=$(_boundary_debug_value(candidates, row, :boundary_shape_score))",
            "shape_fit=$(_boundary_debug_value(candidates, row, :shape_fit_loss))",
            "shape_prior=$(_boundary_debug_value(candidates, row, :shape_prior_penalty))",
            "shape_fit_prior=$(_boundary_debug_value(candidates, row, :shape_fit_prior_score))",
            "shape_fit_prior_norm=$(_boundary_debug_value(candidates, row, :shape_normalized_fit_prior_score))",
            "deconv_area_explained=$(_boundary_debug_value(candidates, row, :shape_deconvolution_area_fraction))",
            "unexplained_area_penalty=$(_boundary_debug_value(candidates, row, :shape_unexplained_area_penalty))",
            "shape_sigma_tau=$(_boundary_debug_value(candidates, row, :shape_sigma_irt)):$(_boundary_debug_value(candidates, row, :shape_tau_irt))",
            "shape_edge=$(_boundary_debug_value(candidates, row, :shape_edge_fraction))",
            "shape_overshoot=$(_boundary_debug_value(candidates, row, :shape_overshoot_penalty))",
            "fallback=$(_boundary_debug_value(candidates, row, :is_fallback))",
            "primary_apex=$(_boundary_debug_value(candidates, row, :is_primary_apex))",
            "quant_trace_selected=$(_boundary_debug_value(candidates, row, :quant_trace_selected))",
            "isotope_key=$(_boundary_debug_value(candidates, row, :isotope_key))",
        ), " "))
    end

    return lines
end

function log_boundary_candidate_debug(
    candidates::AbstractDataFrame,
    selected_candidates::AbstractDataFrame;
    target_precursor_idx::Union{Nothing,UInt32} = nothing,
    target_ms_file_idx::Union{Nothing,UInt16} = nothing,
)
    DEBUG_CONSOLE_LEVEL[] >= 1 || return nothing
    targets = if target_precursor_idx !== nothing && target_ms_file_idx !== nothing
        [(target_precursor_idx, target_ms_file_idx)]
    else
        sort!(collect(DEBUG_BOUNDARY_CANDIDATE_TARGETS[]))
    end

    for (precursor_idx, ms_file_idx) in targets
        lines = boundary_candidate_debug_lines(
            candidates,
            selected_candidates;
            target_precursor_idx = precursor_idx,
            target_ms_file_idx = ms_file_idx,
        )
        @debug_l1 join(lines, "\n")
    end
    return nothing
end

function select_boundary_candidate_rows_by_shape(
    candidates::DataFrame;
    group_col::Symbol = :boundary_group_id,
    score_col::Symbol = :boundary_model_score,
)
    selection_candidates = quant_boundary_candidate_rows(candidates)
    nrow(selection_candidates) == 0 && return copy(selection_candidates)
    scored = hasproperty(selection_candidates, score_col) ?
        copy(selection_candidates) :
        score_boundary_candidates_by_shape!(copy(selection_candidates); score_col = score_col)

    rows_by_group, _ = _fallback_rows_by_group(scored, group_col, :is_fallback)
    fallback = _get_column_or_fill(scored, :is_fallback, false)
    selected_rows = Int[]
    sizehint!(selected_rows, length(rows_by_group))

    for (_, rows) in rows_by_group
        best_row = rows[1]
        best_score = Float32(scored[best_row, score_col])
        for row in rows
            score = Float32(scored[row, score_col])
            if score > best_score ||
               (score == best_score && Bool(fallback[row]) && !Bool(fallback[best_row]))
                best_row = row
                best_score = score
            end
        end
        push!(selected_rows, best_row)
    end

    sort!(selected_rows)
    return scored[selected_rows, :]
end

function select_fallback_boundary_candidate_rows(
    candidates::AbstractDataFrame;
    group_col::Symbol = :boundary_group_id,
)
    selection_candidates = quant_boundary_candidate_rows(candidates)
    nrow(selection_candidates) == 0 && return copy(selection_candidates)

    _, fallback_row_by_group = _fallback_rows_by_group(selection_candidates, group_col, :is_fallback)
    selected_rows = sort!(collect(values(fallback_row_by_group)))
    return selection_candidates[selected_rows, :]
end

function append_boundary_candidate_rows!(
    rows::Vector{NamedTuple},
    candidates,
    boundary_group_id::UInt64,
    ms_file_idx::Integer,
    precursor_idx::UInt32,
    cv_fold::UInt8,
    protein_key,
    sequence_key,
    isotope_key,
    target::Bool,
    quant_trace_selected::Bool,
)
    candidates === nothing && return rows

    for candidate in candidates
        push!(rows, (;
            boundary_group_id = boundary_group_id,
            ms_file_idx = UInt16(ms_file_idx),
            precursor_idx = precursor_idx,
            cv_fold = cv_fold,
            protein_key = protein_key,
            sequence_key = sequence_key,
            isotope_key = isotope_key,
            target = target,
            quant_trace_selected = quant_trace_selected,
            candidate...
        ))
    end

    return rows
end

function boundary_candidate_table(
    boundary_candidate_data::AbstractVector,
    psms::DataFrame,
    search_context::SearchContext,
    ms_file_idx::Integer;
    extra_rows::Vector{NamedTuple} = NamedTuple[],
)
    precursors = getPrecursors(getSpecLib(search_context))
    has_isotope_key = hasproperty(psms, :isotopes_captured)
    has_target = hasproperty(psms, :target)
    rows = copy(extra_rows)

    for i in eachindex(boundary_candidate_data)
        candidates = boundary_candidate_data[i]
        candidates === nothing && continue

        pid = UInt32(psms[i, :precursor_idx])
        boundary_group_id = (UInt64(ms_file_idx) << 32) | UInt64(i)
        protein_key = getAccessionNumbers(precursors)[pid]
        sequence_key = getSequence(precursors)[pid]
        cv_fold = getCvFold(precursors, pid)
        isotope_key = has_isotope_key ? psms[i, :isotopes_captured] : (Int8(0), Int8(0))
        target = has_target ? Bool(psms[i, :target]) : true

        append_boundary_candidate_rows!(
            rows,
            candidates,
            boundary_group_id,
            ms_file_idx,
            pid,
            cv_fold,
            protein_key,
            sequence_key,
            isotope_key,
            target,
            true,
        )
    end

    return DataFrame(rows)
end

function boundary_candidate_dir(search_context::SearchContext)
    return joinpath(getDataOutDir(search_context), "temp_data", "boundary_candidates")
end

function boundary_candidate_path(search_context::SearchContext, ms_file_idx::Integer)
    return joinpath(boundary_candidate_dir(search_context), "file_$(Int(ms_file_idx)).arrow")
end

function write_boundary_candidate_table!(
    boundary_candidate_data::AbstractVector,
    psms::DataFrame,
    search_context::SearchContext,
    ms_file_idx::Integer,
    extra_rows::Vector{NamedTuple} = NamedTuple[],
)
    candidates = boundary_candidate_table(
        boundary_candidate_data,
        psms,
        search_context,
        ms_file_idx;
        extra_rows = extra_rows,
    )
    nrow(candidates) == 0 && return nothing

    out_dir = boundary_candidate_dir(search_context)
    isdir(out_dir) || mkpath(out_dir)
    writeArrow(boundary_candidate_path(search_context, ms_file_idx), candidates)
    return nothing
end

function collect_isotope_boundary_candidate_rows(
    chromatograms::DataFrame,
    passing_psms::DataFrame,
    selected_quant_trace::Union{Nothing, Dict{UInt32, Tuple{Tuple{Int8, Int8}, Float32}}},
    search_context::SearchContext,
    ms_file_idx::Integer,
    min_fraction_transmitted::Float32,
    λ::Float32;
    max_groups::Integer = typemax(Int),
    rng::AbstractRNG = Random.default_rng(),
)
    nrow(chromatograms) == 0 && return NamedTuple[]
    hasproperty(chromatograms, :isotopes_captured) || return NamedTuple[]

    # The quant path may be combined by precursor only. For learner-only isotope
    # traces, regroup the already extracted points by captured isotope set after
    # the fallback quantification has finished.
    sort_chromatograms_for_integration!(chromatograms, SeperateTraces())

    precursors = getPrecursors(getSpecLib(search_context))
    rt_to_irt_model = getRtIrtModel(search_context, ms_file_idx)
    chrom_index, max_chrom_len = build_chrom_index(chromatograms, SeperateTraces())
    max_chrom_len <= 0 && return NamedTuple[]

    rt_all = chromatograms[!, :rt]::AbstractVector{Float32}
    scan_idx_all = chromatograms[!, :scan_idx]::AbstractVector{UInt32}
    intensity_all = chromatograms[!, :intensity]::AbstractVector{Float32}
    fraction_all = chromatograms[!, :precursor_fraction_transmitted]::AbstractVector{Float32}

    psm_row_by_pid = Dict{UInt32, Int}()
    for i in 1:nrow(passing_psms)
        psm_row_by_pid[UInt32(passing_psms[i, :precursor_idx])] = i
    end

    ws = WHWorkspace(max_chrom_len)
    state = Chromatogram(zeros(Float32, max_chrom_len), zeros(Float32, max_chrom_len), 0)
    rows = NamedTuple[]
    extra_group_idx = 0
    has_target = hasproperty(passing_psms, :target)

    chrom_keys = collect(keys(chrom_index))
    Random.shuffle!(rng, chrom_keys)
    max_iter = min(length(chrom_keys), max(0, Int(max_groups)))

    for chrom_key in @view chrom_keys[1:max_iter]
        chrom_range = chrom_index[chrom_key]
        pid, isotope_key = chrom_key
        psm_row = get(psm_row_by_pid, pid, 0)
        psm_row == 0 && continue

        selected = selected_quant_trace === nothing ? nothing : get(selected_quant_trace, pid, nothing)
        selected !== nothing && isotope_key == selected[1] && continue
        isnothing(findfirst(x -> x > 0.0f0, @view intensity_all[chrom_range])) && continue

        avg_cycle_time = (rt_all[last(chrom_range)] - rt_all[first(chrom_range)]) / length(chrom_range)
        apex_scan = find_nearest_scan(
            @view(scan_idx_all[chrom_range]),
            UInt32(passing_psms[psm_row, :scan_idx]),
        )

        candidate_ref = Ref{Any}(nothing)
        integrate_chrom(
            @view(rt_all[chrom_range]),
            @view(scan_idx_all[chrom_range]),
            @view(intensity_all[chrom_range]),
            @view(fraction_all[chrom_range]),
            apex_scan,
            ws,
            state,
            avg_cycle_time,
            λ;
            min_fraction_transmitted = min_fraction_transmitted,
            mainsearch_1pct_start_scan = hasproperty(passing_psms, :mainsearch_1pct_start_scan) ?
                UInt32(passing_psms[psm_row, :mainsearch_1pct_start_scan]) :
                UInt32(0),
            mainsearch_1pct_stop_scan = hasproperty(passing_psms, :mainsearch_1pct_stop_scan) ?
                UInt32(passing_psms[psm_row, :mainsearch_1pct_stop_scan]) :
                UInt32(0),
            rt_to_irt_model = rt_to_irt_model,
            boundary_candidate_data = candidate_ref,
        )
        reset!(state)
        candidate_ref[] === nothing && continue

        extra_group_idx += 1
        boundary_group_id = (UInt64(ms_file_idx) << 32) |
            UInt64(nrow(passing_psms) + extra_group_idx)
        target = has_target ? Bool(passing_psms[psm_row, :target]) : true
        append_boundary_candidate_rows!(
            rows,
            candidate_ref[],
            boundary_group_id,
            ms_file_idx,
            pid,
            getCvFold(precursors, pid),
            getAccessionNumbers(precursors)[pid],
            getSequence(precursors)[pid],
            isotope_key,
            target,
            false,
        )
    end

    return rows
end

function load_boundary_candidate_tables(search_context::SearchContext)
    candidate_dir = boundary_candidate_dir(search_context)
    isdir(candidate_dir) || return DataFrame()

    paths = filter(path -> endswith(path, ".arrow"), readdir(candidate_dir; join = true))
    isempty(paths) && return DataFrame()

    tables = DataFrame[]
    sizehint!(tables, length(paths))
    for path in paths
        push!(tables, DataFrame(Tables.columntable(Arrow.Table(path))))
    end
    isempty(tables) && return DataFrame()
    return reduce(vcat, tables; cols = :union)
end

function _candidate_lookup_key(precursor_idx, isotope_key)
    return (UInt32(precursor_idx), isotope_key)
end

function apply_selected_boundary_candidates!(
    selected_candidates::DataFrame,
    search_context::SearchContext,
)
    nrow(selected_candidates) == 0 && return 0

    updated = 0
    ms_data = getMSData(search_context)
    for group in groupby(selected_candidates, :ms_file_idx)
        ms_file_idx = Int(first(group.ms_file_idx))
        psm_path = getPassingPsms(ms_data)[ms_file_idx]
        (isempty(psm_path) || !isfile(psm_path)) && continue

        lookup = Dict{Tuple{UInt32, Any}, Tuple{Float32, UInt32, UInt32}}()
        for row in eachrow(group)
            area = Float32(row.peak_area)
            area > 0.0f0 && isfinite(area) || continue
            lookup[_candidate_lookup_key(row.precursor_idx, row.isotope_key)] = (
                area,
                UInt32(row.new_best_scan),
                UInt32(row.points_integrated),
            )
        end
        isempty(lookup) && continue

        psms = DataFrame(Tables.columntable(Arrow.Table(psm_path)))
        nrow(psms) == 0 && continue
        has_isotope_key = hasproperty(psms, :isotopes_captured)

        for i in 1:nrow(psms)
            isotope_key = has_isotope_key ? psms[i, :isotopes_captured] : (Int8(0), Int8(0))
            value = get(lookup, _candidate_lookup_key(psms[i, :precursor_idx], isotope_key), nothing)
            value === nothing && continue

            psms[i, :peak_area] = value[1]
            psms[i, :new_best_scan] = value[2]
            psms[i, :points_integrated] = value[3]
            updated += 1
        end

        writeArrow(psm_path, psms)
    end

    return updated
end
