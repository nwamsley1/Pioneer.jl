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

_split_title(fname::String, desc::String; max_len::Int=45) =
    length(fname) + length(desc) > max_len ? "$fname\n$desc" : "$fname: $desc"

"""
    add_tuning_search_columns!(psms::DataFrame, MS_TABLE::Arrow.Table,
                             prec_is_decoy::Arrow.BoolVector{Bool},
                             prec_irt::Arrow.Primitive{T, Vector{T}},
                             prec_charge::Arrow.Primitive{UInt8, Vector{UInt8}},
                             scan_retention_time::AbstractVector{Float32},
                             tic::AbstractVector{Float32}) where {T<:AbstractFloat}

Adds essential columns to PSM DataFrame for parameter tuning analysis.

# Arguments
- `psms`: DataFrame containing PSMs to modify
- `MS_TABLE`: Mass spectrometry data table
- `prec_is_decoy`: Boolean vector indicating decoy status
- `prec_irt`: Vector of iRT values
- `prec_charge`: Vector of precursor charges
- `scan_retention_time`: Vector of scan retention times
- `tic`: Vector of total ion currents

# Added Columns
- Basic metrics: RT, iRT predicted, charge, TIC
- Analysis columns: target, decoy
- Scoring columns: q_value, prob, intercept

Uses parallel processing for efficiency through data chunking.
"""
function add_tuning_search_columns!(psms::DataFrame, 
                                MS_TABLE::MassSpecData, 
                                prec_is_decoy::Arrow.BoolVector{Bool},
                                prec_irt::Arrow.Primitive{T, Vector{T}},
                                prec_charge::Arrow.Primitive{UInt8, Vector{UInt8}},
                                scan_retention_time::AbstractVector{Float32},
                                tic::AbstractVector{Float32}) where {T<:AbstractFloat}
    
    N = size(psms, 1)
    decoys = zeros(Bool, N);
    targets = zeros(Bool, N);
    TIC = zeros(Float16, N);
    charge = zeros(UInt8, N);
    spectrum_peak_count = zeros(UInt32, N);
    irt_pred = zeros(Float32, N);
    rt = zeros(Float32, N);
    scan_idx::Vector{UInt32} = psms[!,:scan_idx]
    precursor_idx::Vector{UInt32} = psms[!,:precursor_idx]

    parallel_foreach!(size(psms, 1); tasks_per_thread=10) do chunk
        for i in chunk
            decoys[i] = prec_is_decoy[precursor_idx[i]];
            targets[i] = decoys[i] == false
            irt_pred[i] = Float32(prec_irt[precursor_idx[i]]);
            rt[i] = Float32(scan_retention_time[scan_idx[i]]);
            TIC[i] = Float16(log2(tic[scan_idx[i]]));
            charge[i] = UInt8(prec_charge[precursor_idx[i]]);
        end
    end
    psms[!,:decoy] = decoys
    psms[!,:irt_predicted] = irt_pred
    psms[!,:rt] = rt
    psms[!,:TIC] = TIC
    psms[!,:target] = targets
    psms[!,:charge] = charge
    psms[!,:spectrum_peak_count] = spectrum_peak_count
    psms[!,:intercept] = ones(Float16, size(psms, 1))
    psms[!,:q_value] = zeros(Float16, N);
    psms[!,:prob] = zeros(Float16, N);
end


#==========================================================
PSM diagnostics
==========================================================#

#==========================================================
Bitmask filter construction for scout/collect
==========================================================#

# Build a 256-entry LUT requiring the top `n_required` ranked fragments
# (bits 0..n_required-1) to be present AND total count_ones >= min_total.
# Fragment ranks use 1<<(rank-1): rank 1 = bit 0, rank 7 = bit 6.
function make_top_n_required_lut(n_required::Int, min_total::Int)
    required_mask = UInt8((1 << n_required) - 1)  # e.g. n=4 → 0x0F
    lut = Vector{Bool}(undef, 256)
    @inbounds for p in 0:255
        lut[p + 1] = (UInt8(p) & required_mask) == required_mask &&
                     count_ones(UInt8(p)) >= min_total
    end
    return lut
end

#==========================================================
Scout: lightweight target/decoy excess estimation
==========================================================#

# Fast scout: fragment index search only (no spectral matching or scoring).
# Counts target vs decoy precursors that pass the bitmask filter.
# Returns (n_target, n_decoy, excess) where excess = n_target - n_decoy.
"""
    score_presearch!(psms::DataFrame)

Performs initial scoring of PSMs using probit regression.

# Arguments
- `psms`: DataFrame containing PSMs to score

# Process
1. Uses features: gof, fitted_manhattan_distance, fitted_hellinger, weight,
   y_count, error, TIC, intercept
2. Performs probit regression with 20 max iterations
3. Calculates probability scores for each PSM

Adds 'prob' column to psms DataFrame containing probability scores.
Uses parallel processing for efficiency.
"""
function score_presearch!(psms::DataFrame)
    # Select feature set based on which columns are present (deconvolution vs simple path)
    features = if hasproperty(psms, :gof)
        [:gof,:fitted_manhattan_distance,
         :fitted_hellinger,:poisson,
         :y_count,:error,
         :TIC,:intercept]
    else
        [:entropy_score,:city_block,
         :scribe,:spectral_contrast,
         :y_count,:error,
         :TIC,:intercept]
    end
    psms[!,:prob] = zeros(Float32, size(psms, 1))

    M = size(psms, 1)
    if M > 50
        try
            tasks_per_thread = 10
            chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
            data_chunks = partition(1:M, chunk_size)
            β = zeros(Float64, length(features))
            β = ProbitRegression(β, psms[!,features], psms[!,:target], data_chunks, max_iter = 20)
            ModelPredictProbs!(psms[!,:prob], psms[!,features], β, data_chunks)
        catch
            psms[!,:prob] = ones(Float32, size(psms, 1))
        end
    else
        psms[!,:prob] = ones(Float32, size(psms, 1))
    end

end


# Note: get_qvalues! and get_PEP! functions have been moved to src/utils/ML/fdrUtilities.jl
# They are imported through the module loading system and are available in the Pioneer module namespace

#==========================================================
RT modeling
==========================================================#

"""
    filter_top_psms_per_precursor(psms::DataFrame, max_per_precursor::Int=3)

Filter PSMs to keep only the top N PSMs per unique precursor_idx.
Selects PSMs with highest probability scores to avoid over-representation
of abundant precursors in RT model fitting.

# Arguments
- `psms`: DataFrame containing PSMs with precursor_idx and prob columns
- `max_per_precursor`: Maximum number of PSMs to keep per precursor (default: 3)

# Returns
- Filtered DataFrame with at most max_per_precursor PSMs per precursor_idx
"""
function filter_top_psms_per_precursor(psms::DataFrame, max_per_precursor::Int=3)
    if isempty(psms) || !("precursor_idx" in names(psms))
        return psms
    end

    # Group by precursor_idx and keep top N by probability score
    filtered_psms = combine(groupby(psms, :precursor_idx)) do group
        n_to_keep = min(nrow(group), max_per_precursor)
        if n_to_keep == nrow(group)
            return group
        end
        # Sort by prob (descending) and take top N
        sorted_indices = sortperm(group.prob, rev=true)
        return group[sorted_indices[1:n_to_keep], :]
    end

    return filtered_psms
end

"""
    fit_irt_model(params::P, psms::DataFrame) where {P<:ParameterTuningSearchParameters}

Wrapper around shared fit_irt_model from rt_alignment_utils.jl.
Adapts ParameterTuningSearchParameters to keyword arguments.

# Arguments
- `params`: Parameter tuning search parameters
- `psms`: DataFrame containing PSMs with RT information

# Returns
Tuple containing:
- Final RT conversion model (SplineRtConversionModel, LinearRtConversionModel, or IdentityModel)
- Valid RT values
- Valid iRT values
- iRT median absolute deviation
"""
function fit_irt_model(
    params::P,
    psms::DataFrame
) where {P<:ParameterTuningSearchParameters}

    # Call shared fit_irt_model from CommonSearchUtils/rt_alignment_utils.jl
    # with extracted parameters
    return Pioneer.fit_irt_model(
        psms;
        lambda_penalty = Float32(getRtAlignmentLambdaPenalty(params)),
        ransac_threshold = getRtAlignmentRansacThreshold(params),
        min_psms = getRtAlignmentMinPsms(params),
        spline_degree = getSplineDegree(params),
        max_knots = getSplineNKnots(params),
        outlier_threshold = Float32(getOutlierThreshold(params))
    )
end

function calculate_ppm_errors(samples::AbstractVector{MassErrSample})
    return [(s.observed_mz - s.theoretical_mz)/(s.theoretical_mz/1f6) for s in samples]
end

#==========================================================
Mass Error Modeling
==========================================================#

function getScanToPrecIdx(scan_idxs::Vector{UInt32}, n_scans::Int64)
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, n_scans)
    start_idx = stop_idx = 1
    
    for i in 1:length(scan_idxs)
        stop_idx = i
        if scan_idxs[start_idx] == scan_idxs[stop_idx]
            scan_to_prec_idx[scan_idxs[i]] = start_idx:stop_idx
        else
            scan_to_prec_idx[scan_idxs[i]] = i:i
            start_idx = i
        end
    end
    scan_to_prec_idx
end

"""
    mass_error_search(spectra::MassSpecData, scan_idxs::Vector{UInt32},
                     precursor_idxs::Vector{UInt32}, ms_file_idx::UInt32,
                     spec_lib::SpectralLibrary, search_data::AbstractVector{S},
                     mem::M, params::P) where {M<:MassErrorModel, S<:SearchDataStructures, P<:SearchParameters}

Performs mass error-focused fragment matching.

# Arguments
- `spectra`: MS spectra data
- `scan_idxs`: Scan indices to process
- `precursor_idxs`: Precursor indices
- `ms_file_idx`: MS file index
- `spec_lib`: Spectral library
- `search_data`: Search data structures
- `mem`: Mass error model
- `params`: Search parameters

# Process
1. Maps scans to precursor indices
2. Selects transitions for mass error estimation
3. Matches peaks and collects errors
4. Processes matches in parallel across threads

Returns matched fragments for mass error analysis.
"""
# Default method without chromatogram type - calls MS2CHROM version
function mass_error_search(
    spectra::MassSpecData,
    scan_idxs::Vector{UInt32},
    precursor_idxs::Vector{UInt32},
    ms_file_idx::UInt32,
    spec_lib::SpectralLibrary,
    search_data::AbstractVector{S},
    mem::M,
    params::P,
    nce_model::NceModel{Float32}
) where {
        M<:AbstractMassErrorModel,
        S<:SearchDataStructures,
        P<:SearchParameters
        }
    return mass_error_search(spectra, scan_idxs, precursor_idxs, ms_file_idx,
                            spec_lib, search_data, mem, params, nce_model, MS2CHROM())
end

function mass_error_search(
    spectra::MassSpecData,
    scan_idxs::Vector{UInt32},
    precursor_idxs::Vector{UInt32},
    ms_file_idx::UInt32,
    spec_lib::SpectralLibrary,
    search_data::AbstractVector{S},
    mem::M,
    params::P,
    nce_model::NceModel{Float32},
    ::MS2CHROM
) where {
        M<:AbstractMassErrorModel,
        S<:SearchDataStructures,
        P<:SearchParameters
        }

    # Sort scans and setup thread tasks
    sorted_indices = sortperm(scan_idxs)
    scan_idxs = scan_idxs[sorted_indices]
    precursor_idxs = precursor_idxs[sorted_indices]

    scan_to_prec_idx = getScanToPrecIdx(scan_idxs, length(spectra))
    thread_tasks = partition_scans(spectra, Threads.nthreads())

    ion_list = getFragmentLookupTable(spec_lib)
    max_rank = Int64(getMaxFragsForMassErrEstimation(params))

    # Process mass errors in parallel — fused kernel writes MassErrSample
    # triples (theoretical_mz, raw observed_mz, intensity) directly into the
    # per-thread mass_err_samples buffer.
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            sd = search_data[thread_id]
            corr_mz   = getScanCorrectedMz(sd)
            obs_low   = getScanObsLow(sd)
            obs_high  = getScanObsHigh(sd)
            samples   = getMassErrSamples(sd)
            sample_idx = 0

            for scan_idx in last(thread_task)
                (scan_idx == 0 || scan_idx > length(spectra)) && continue
                ismissing(scan_to_prec_idx[scan_idx]) && continue

                scan_mz  = getMzArray(spectra, scan_idx)
                scan_int = getIntensityArray(spectra, scan_idx)
                scan_rt = Float32(getRetentionTime(spectra, scan_idx))
                peak_mz_len = prepare_scan_peaks!(corr_mz, obs_low, obs_high,
                                                  mem, scan_mz, scan_int, scan_rt)

                sample_idx = run_fused_masserr!(
                    samples, sample_idx,
                    corr_mz, obs_low, obs_high, peak_mz_len,
                    ion_list,
                    precursor_idxs, scan_to_prec_idx[scan_idx],
                    mem, scan_mz, scan_int, scan_rt;
                    max_rank = max_rank)
            end

            @view(getMassErrSamples(sd)[1:sample_idx])
        end
    end
    fetch.(tasks)
end

function fit_mass_err_model(
    params::P,
    samples::AbstractVector{MassErrSample}
) where {P<:FragmentIndexSearchParameters}
    
    # Check if we have enough samples after filtering
    if length(samples) == 0
        return MassErrorModel(
            zero(Float32),
            (zero(Float32), zero(Float32))
        ), Float32[]
    end

    # Calculate PPM errors
    ppm_errs = calculate_ppm_errors(samples)
    mass_err = median(ppm_errs)
    ppm_errs .-= mass_err
    
    # Calculate error bounds using configured quantile
    # Exponential scale estimated from median absolute error (robust to outliers):
    #   median = scale * ln(2)  →  scale = median / ln(2)
    frag_err_quantile = getFragErrQuantile(params)
    r_errs = abs.(ppm_errs[ppm_errs .> 0.0f0])
    l_errs = abs.(ppm_errs[ppm_errs .< 0.0f0])
    r_bound, l_bound = nothing, nothing
    r_scale, l_scale = 0.0, 0.0
    try
    r_scale = median(r_errs) / log(2.0)
    l_scale = median(l_errs) / log(2.0)
    r_bound = quantile(Exponential(r_scale), 1 - frag_err_quantile)
    l_bound = quantile(Exponential(l_scale), 1 - frag_err_quantile)
    catch e
        @debug_l1 "length(r_errs) $(length(r_errs)) length(l_errs) $(length(l_errs)) sum(r_errs) $(sum(r_errs)) sum(l_errs) $(sum(l_errs))"
        return MassErrorModel(
            zero(Float32),
            (zero(Float32), zero(Float32))
        ), ppm_errs
    end

    # Empirical quantiles for comparison
    empirical_r_quantile = length(r_errs) > 0 ? quantile(r_errs, 1 - frag_err_quantile) : 0.0
    empirical_l_quantile = length(l_errs) > 0 ? quantile(l_errs, 1 - frag_err_quantile) : 0.0

    @debug_l1 "  SimpleMassErrorModel tolerance fitting:" *
        "\n    offset: $(round(mass_err, digits=3)) ppm" *
        "\n    n_fragments: $(length(ppm_errs)) ($(length(r_errs)) right, $(length(l_errs)) left)" *
        "\n    median |err|: right=$(round(median(r_errs), digits=3)), left=$(round(median(l_errs), digits=3)) ppm" *
        "\n    exponential scale (median/ln2): right=$(round(r_scale, digits=3)), left=$(round(l_scale, digits=3)) ppm" *
        "\n    quantile=$(1 - frag_err_quantile) → bounds: -$(round(l_bound, digits=2)) / +$(round(r_bound, digits=2)) ppm" *
        "\n    empirical $(1 - frag_err_quantile) quantile: -$(round(empirical_l_quantile, digits=2)) / +$(round(empirical_r_quantile, digits=2)) ppm"

    return MassErrorModel(
        Float32(mass_err),
        (Float32(abs(l_bound)), Float32(abs(r_bound)))
    ), ppm_errs
end

"""
    fit_and_install_intensity_model!(spectra, search_context, params, ms_file_idx,
                                     scan_indices, fragments)

After tolerance scanning, fit an IntensityMassErrorModel from the collected
fragment data and install it as the mass error model for this file.
The SimpleMassErrorModel remains available as fallback.
Spread is fitted in mDa space; bias stays in Da.
"""
function fit_and_install_intensity_model!(
    search_context::SearchContext,
    ms_file_idx::Int64,
    samples::AbstractVector{MassErrSample};
    k::Float32 = 1.96f0
)
    old_model = getMassErrorModel(search_context, ms_file_idx)
    new_model = fit_intensity_mass_error_model(samples, old_model; k=k)
    if new_model isa IntensityMassErrorModel
        setMassErrorModel!(search_context, ms_file_idx, new_model)
        @debug_l1 "Installed IntensityMassErrorModel for file $ms_file_idx"
    end
end

# Deprecated: apply_mass_error_buffer removed – buffer applied inline per-file before plotting

function get_matched_fragments(
    spectra::MassSpecData,
    psms::DataFrame,
    search_context::SearchContext,
    params::P,
    ms_file_idx::Int64
) where {P<:FragmentIndexSearchParameters}

    return vcat(mass_error_search(
        spectra,
        psms[!,:scan_idx],
        psms[!,:precursor_idx],
        UInt32(ms_file_idx),
        getSpecLib(search_context),
        getSearchData(search_context),
        getMassErrorModel(search_context, ms_file_idx),
        params,
        getNceModel(search_context, ms_file_idx),
        MS2CHROM()
    )...)
end

# Fit mass error model from pre-collected samples (no sample collection).
function fit_models_from_fragments(
    params::P,
    samples::AbstractVector{MassErrSample}
) where {P<:FragmentIndexSearchParameters}
    if isempty(samples)
        return nothing, nothing, nothing
    end
    mass_err_model, ppm_errs = fit_mass_err_model(params, samples)
    frag_mzs = Float32[s.theoretical_mz for s in samples]
    return mass_err_model, ppm_errs, frag_mzs
end

#==========================================================
Plotting Helpers
==========================================================#

"""
    generate_rt_plot(results::ParameterTuningSearchResults, plot_path::String, title::String)

Generates retention time alignment visualization plot.

# Arguments
- `results`: Parameter tuning results containing RT data
- `plot_path`: Path to save plot
- `title`: Plot title

Creates scatter plot of RT vs iRT with fitted spline curve.
"""
function generate_rt_plot(
    results::ParameterTuningSearchResults,
    title::String;
    irt_tol::Union{Float32, Float64, Nothing} = nothing
)
    n = length(results.rt)
    p = plot(
        results.rt,
        results.irt,
        seriestype=:scatter,
        title = title*"\n n = $n",
        xlabel = "Retention Time RT (min)",
        ylabel = "Indexed Retention Time iRT (min)",
        label = nothing,
        alpha = 0.1,
        size = (900, 900)
    )

    pbins = LinRange(minimum(results.rt), maximum(results.rt), 100)
    spline_vals = getRtToIrtModel(results).(pbins)
    plot!(pbins, spline_vals, lw=3, label=nothing)

    # Add ±iRT tolerance band as dotted lines, or indicate infinite tolerance
    if irt_tol !== nothing && isfinite(irt_tol)
        plot!(pbins, spline_vals .+ irt_tol, lw=1.5, ls=:dot, color=:red, label="±$(round(irt_tol, digits=2)) iRT")
        plot!(pbins, spline_vals .- irt_tol, lw=1.5, ls=:dot, color=:red, label=nothing)
    elseif irt_tol === nothing || !isfinite(irt_tol)
        annotate!([(mean(pbins), maximum(spline_vals),
            Plots.text("iRT tolerance: ∞ (insufficient PSMs)", 10, :red, :center))])
    end

    return p
end


"""
    generate_mass_error_plot_mda(fragments, model, fname)

Generates bias-corrected mDa mass error histogram with tolerance bounds.
Uses the IntensityMassErrorModel to subtract m/z + intensity + RT bias, then shows
both the max observed error clamp and the empirical 99.9% quantile clamp.
"""
function generate_mass_error_plot_mda(
    frag_data::NamedTuple,
    model::IntensityMassErrorModel,
    fname::String
)
    frag_mzs = frag_data.frag_mzs
    log2I = frag_data.log2I
    rts = frag_data.rts
    da_errs = frag_data.da_errs
    n = length(da_errs)
    n < 10 && return nothing

    raw_mda = da_errs .* 1e3
    raw_ppm = [da_errs[i] / (Float64(frag_mzs[i]) / 1e6) for i in 1:n]

    corrected_mda = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        bias_da = Float64(_total_bias_da(
            model,
            Float32(frag_mzs[i]),
            Float32(log2I[i]),
            Float32(rts[i]),
        ))
        corrected_mda[i] = (da_errs[i] - bias_da) * 1e3
    end

    plots = Plots.Plot[]

    # ── Page 1: Raw mDa errors (no bias correction) ──
    raw_median = median(raw_mda)
    raw_range = ceil(max(quantile(abs.(raw_mda), 0.999), abs(raw_median)) * 2.0)
    raw_bins = collect(-raw_range:1.0:raw_range)  # 1 mDa bin width
    p_raw = Plots.histogram(raw_mda,
        orientation = :h, yflip = true, size = (900, 900),
        title = _split_title(fname, "raw mDa") * "\nn=$n, median=$(round(raw_median, digits=2)) mDa",
        xlabel = "Count", ylabel = "Raw error (mDa)",
        label = nothing, bins = raw_bins,
        ylim = (-raw_range, raw_range), topmargin = 15mm)
    Plots.hline!(p_raw, [0.0], label = nothing, color = :black, lw = 1, ls = :dash)
    Plots.hline!(p_raw, [raw_median], label = nothing, color = :red, lw = 2)
    Plots.annotate!(p_raw, last(xlims(p_raw)), raw_median,
        text("median: $(round(raw_median, digits=2)) mDa", :red, :right, :bottom, 10))
    push!(plots, p_raw)

    # ── Page 2: Bias-corrected mDa errors ──
    max_err_mda = Float64(model.max_da_tolerance) * 1e3
    q999_mda = quantile(abs.(corrected_mda), 0.999)
    plot_range = ceil(max(max_err_mda, q999_mda) * 1.5)
    bins = collect(-plot_range:1.0:plot_range)  # 1 mDa bin width
    p_corr = Plots.histogram(corrected_mda,
        orientation = :h, yflip = true, size = (900, 900),
        title = _split_title(fname, "bias-corrected mDa") * "\nn=$n",
        xlabel = "Count", ylabel = "Corrected error (mDa)",
        label = nothing, bins = bins,
        ylim = (-plot_range, plot_range), topmargin = 15mm)
    Plots.hline!(p_corr, [0.0], label = nothing, color = :black, lw = 1, ls = :dash)
    Plots.hline!(p_corr, [max_err_mda, -max_err_mda], label = nothing, color = :black, lw = 2)
    Plots.annotate!(p_corr, last(xlims(p_corr)), max_err_mda,
        text("max obs: $(round(max_err_mda, digits=2)) mDa", :black, :right, :bottom, 10))
    Plots.hline!(p_corr, [q999_mda, -q999_mda], label = nothing, color = :red, lw = 2, ls = :dash)
    Plots.annotate!(p_corr, last(xlims(p_corr)), q999_mda,
        text("99.9%: $(round(q999_mda, digits=2)) mDa", :red, :right, :top, 10))
    push!(plots, p_corr)

    return plots
end

# Fallback for SimpleMassErrorModel or missing fragments
function generate_mass_error_plot_mda(frag_data::NamedTuple, model::SimpleMassErrorModel, fname::String;
)
    mda = frag_data.mda_errs
    mzs = frag_data.frag_mzs
    n = length(mda)
    n < 5 && return nothing
    raw_median = median(mda)
    raw_range = ceil(quantile(abs.(mda), 0.999) * 1.5)
    raw_bins = collect(-raw_range:1.0:raw_range)

    offset_ppm = getMassOffset(model)
    left_tol = getLeftTol(model)
    right_tol = getRightTol(model)
    model_desc = "SimpleMassErrorModel: offset=$(round(offset_ppm, digits=1))ppm, tol=-$(round(left_tol, digits=1))/+$(round(right_tol, digits=1))ppm"

    p_raw = Plots.histogram(mda,
        orientation = :h, yflip = true, size = (900, 900),
        title = _split_title(fname, "raw mDa (SimpleMassErrorModel)") *
                "\nn=$n, median=$(round(raw_median, digits=2)) mDa\n$model_desc",
        xlabel = "Count", ylabel = "Raw error (mDa)",
        label = nothing, bins = raw_bins,
        ylim = (-raw_range, raw_range), topmargin = 20Plots.mm)
    Plots.hline!(p_raw, [0.0], label = nothing, color = :black, lw = 1, ls = :dash)
    Plots.hline!(p_raw, [raw_median], label = nothing, color = :red, lw = 2)

    # mDa scatter vs m/z (same data, shows m/z dependence)
    mz_range = range(minimum(mzs), maximum(mzs), length=100)
    offset_mda = [offset_ppm * m / 1e3 for m in mz_range]
    upper_mda = [(offset_ppm + right_tol) * m / 1e3 for m in mz_range]
    lower_mda = [(offset_ppm - left_tol) * m / 1e3 for m in mz_range]

    p_scatter = Plots.scatter(Float64.(mzs), Float64.(mda),
        alpha = 0.15, markersize = 1.5, color = :steelblue, label = nothing,
        xlabel = "Fragment m/z", ylabel = "Raw error (mDa)",
        title = _split_title(fname, "raw mDa vs m/z (SimpleMassErrorModel)") *
                "\nn=$n",
        size = (900, 900), topmargin = 15Plots.mm)
    Plots.plot!(p_scatter, mz_range, offset_mda, lw = 2, color = :red, label = "offset $(round(offset_ppm, digits=1)) ppm")
    Plots.plot!(p_scatter, mz_range, upper_mda, lw = 1.5, ls = :dash, color = :red,
        label = "tol: -$(round(left_tol, digits=1))/+$(round(right_tol, digits=1)) ppm")
    Plots.plot!(p_scatter, mz_range, lower_mda, lw = 1.5, ls = :dash, color = :red, label = nothing)
    Plots.hline!(p_scatter, [0.0], color = :black, lw = 1, ls = :dot, label = nothing)

    return [p_raw, p_scatter]
end

#==========================================================
# IntensityMassErrorModel Diagnostic Plots
==========================================================#

"""
    generate_intensity_model_plots(fragments, model, fname)

Generate diagnostic plots for an IntensityMassErrorModel showing:
1. Bias fit: corrected residuals vs m/z and vs log2(intensity) with fitted curves
2. RT-dependent bias: residuals after m/z + intensity correction vs RT
3. Spread fit: Laplace scale vs log2(intensity) with fitted spline
4. Corrected scatter: residuals after full bias correction vs m/z (should be centered at 0)

Returns a vector of Plot objects.
"""
# Extract arrays from FragmentMatch vector for plotting. Returns NamedTuple
# with (frag_mzs, log2I, rts, da_errs, mda_errs). Call once, pass to multiple plotters.
function extract_fragment_plot_data(samples::AbstractVector{MassErrSample})
    n = length(samples)
    frag_mzs = Vector{Float32}(undef, n)
    intensities = Vector{Float32}(undef, n)
    rts = Vector{Float32}(undef, n)
    da_errs = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        s = samples[i]
        frag_mzs[i] = s.theoretical_mz
        intensities[i] = s.intensity
        rts[i] = s.rt
        da_errs[i] = Float64(s.observed_mz - s.theoretical_mz)
    end

    valid = (intensities .> 0.0f0) .& isfinite.(rts)
    frag_mzs = frag_mzs[valid]
    rts = rts[valid]
    da_errs = da_errs[valid]
    log2I = log2.(Float64.(intensities[valid]))
    mda_errs = da_errs .* 1e3

    return (frag_mzs=frag_mzs, log2I=log2I, rts=rts, da_errs=da_errs, mda_errs=mda_errs)
end

function _diagnostic_axis_grid(values::AbstractVector, n_bins::Integer)
    n = max(Int(n_bins), 2)
    finite_values = Float64[x for x in values if isfinite(x)]
    isempty(finite_values) && return collect(range(0.0, 1.0; length=n))

    lo, hi = extrema(finite_values)
    if lo == hi
        pad = max(abs(lo) * 0.01, 1.0)
        lo -= pad
        hi += pad
    end
    return collect(range(lo, hi; length=n))
end

function _model_tolerance_heatmap_grid(
    model::IntensityMassErrorModel,
    frag_mzs::AbstractVector,
    log2I::AbstractVector;
    n_mz_bins::Integer=80,
    n_intensity_bins::Integer=80,
)
    mz_grid = _diagnostic_axis_grid(frag_mzs, n_mz_bins)
    log2I_grid = _diagnostic_axis_grid(log2I, n_intensity_bins)
    tolerance_mda = Matrix{Float64}(undef, length(log2I_grid), length(mz_grid))

    @inbounds for (iy, intensity) in pairs(log2I_grid)
        sigma_mda = Float64(_gaussian_sigma_mda(model, Float32(intensity)))
        for (ix, mz) in pairs(mz_grid)
            correction = Float64(_mz_spread_correction(model, Float32(mz)))
            tolerance_mda[iy, ix] = Float64(model.k) * sigma_mda * correction
        end
    end

    return mz_grid, log2I_grid, tolerance_mda
end

function generate_intensity_model_plots(
    frag_data::NamedTuple,
    model::IntensityMassErrorModel,
    fname::String
)
    plots_out = Plots.Plot[]
    frag_mzs = frag_data.frag_mzs
    log2I = frag_data.log2I
    rts = frag_data.rts
    da_errs = frag_data.da_errs
    mda_errs = frag_data.mda_errs
    n = length(da_errs)
    n < 50 && return plots_out

    # ── Compute model predictions for each fragment ──
    # Bias in Da, spread in ppm (converted to Da per-fragment for calibration)
    mz_bias = Vector{Float64}(undef, n)
    int_bias = Vector{Float64}(undef, n)
    rt_bias = Vector{Float64}(undef, n)
    total_bias = Vector{Float64}(undef, n)
    laplace_scale_mda = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        mz_bias[i] = Float64(_mz_bias_da(model, Float32(frag_mzs[i])))
        int_bias[i] = Float64(_intensity_bias_da(model, Float32(log2I[i])))
        rt_bias[i] = Float64(_rt_bias_da(model, Float32(rts[i])))
        total_bias[i] = mz_bias[i] + int_bias[i] + rt_bias[i]
        laplace_scale_mda[i] = Float64(_laplace_scale_mda(model, Float32(log2I[i])))  # now actually mDa
    end

    corrected = da_errs .- total_bias
    corrected_mda = corrected .* 1e3

    # ── Bias diagnostics: two separate pages ──
    diagnostic_bin_sz = TUNING_CALIBRATION_BIN_SIZE
    n_bins = max(n ÷ diagnostic_bin_sz, 1)
    mz_order = sortperm(frag_mzs)
    mz_f64 = Float64.(frag_mzs)

    # Page 1: m/z-dependent bias
    bin_mz_c = Float64[]
    bin_med_raw = Float64[]
    for b in 1:n_bins
        s = (b - 1) * diagnostic_bin_sz + 1
        e = b == n_bins ? n : b * diagnostic_bin_sz
        idx = mz_order[s:e]
        push!(bin_mz_c, Float64(median(frag_mzs[idx])))
        push!(bin_med_raw, median(mda_errs[idx]))
    end

    p_mz_bias = plot(size=(600, 600), bottommargin=10Plots.mm, leftmargin=10Plots.mm)
    scatter!(p_mz_bias, mz_f64, mda_errs,
             alpha=0.08, markersize=1, color=:gray, label=nothing)
    scatter!(p_mz_bias, bin_mz_c, bin_med_raw,
             xlabel="Fragment m/z", ylabel="Error (mDa)",
             label="bin median", marker=:circle, color=:blue, ms=4)
    mz_range = range(minimum(frag_mzs), maximum(frag_mzs), length=200)
    mz_fit = [Float64(_mz_bias_da(model, Float32(m))) * 1e3 for m in mz_range]
    plot!(p_mz_bias, mz_range, mz_fit, lw=2.5, color=:red, label="m/z bias spline")
    hline!(p_mz_bias, [0.0], color=:black, lw=1, ls=:dot, label=nothing)
    vline!(p_mz_bias, [Float64(model.mz_bias_extrap.lo), Float64(model.mz_bias_extrap.hi)],
           color=:orange, lw=1.5, ls=:dot, label="extrap boundary")
    title!(p_mz_bias, _split_title(fname, "m/z-dependent bias (mDa)"))
    push!(plots_out, p_mz_bias)

    # Page 2: Intensity-dependent bias
    int_order = sortperm(log2I)
    bin_log2I_c = Float64[]
    bin_med_resid = Float64[]
    mz_resid_mda = (da_errs .- mz_bias) .* 1e3
    for b in 1:n_bins
        s = (b - 1) * diagnostic_bin_sz + 1
        e = b == n_bins ? n : b * diagnostic_bin_sz
        idx = int_order[s:e]
        push!(bin_log2I_c, median(log2I[idx]))
        push!(bin_med_resid, median(mz_resid_mda[idx]))
    end

    log2I_range = range(minimum(log2I), maximum(log2I), length=200)
    p_int_bias = plot(size=(600, 600), bottommargin=10Plots.mm, leftmargin=10Plots.mm)
    scatter!(p_int_bias, log2I, mz_resid_mda,
             alpha=0.08, markersize=1, color=:gray, label=nothing)
    scatter!(p_int_bias, bin_log2I_c, bin_med_resid,
             xlabel="log2(intensity)", ylabel="Residual after m/z bias (mDa)",
             label="bin median", marker=:circle, color=:blue, ms=4)
    int_fit = [Float64(_intensity_bias_da(model, Float32(x))) * 1e3 for x in log2I_range]
    plot!(p_int_bias, log2I_range, int_fit, lw=2.5, color=:red, label="intensity bias spline")
    hline!(p_int_bias, [0.0], color=:black, lw=1, ls=:dot, label=nothing)
    vline!(p_int_bias, [Float64(model.int_bias_extrap.lo), Float64(model.int_bias_extrap.hi)],
           color=:orange, lw=1.5, ls=:dot, label="extrap boundary")
    title!(p_int_bias, _split_title(fname, "Intensity-dependent bias (mDa)"))
    push!(plots_out, p_int_bias)

    # Page 3: RT-dependent bias
    rt_raw = Float64.(rts)
    rt_order = sortperm(rt_raw)
    bin_rt_c = Float64[]
    bin_med_rt_resid = Float64[]
    mz_int_resid_mda = (da_errs .- mz_bias .- int_bias) .* 1e3
    for b in 1:n_bins
        s = (b - 1) * diagnostic_bin_sz + 1
        e = b == n_bins ? n : b * diagnostic_bin_sz
        idx = rt_order[s:e]
        push!(bin_rt_c, median(rt_raw[idx]))
        push!(bin_med_rt_resid, median(mz_int_resid_mda[idx]))
    end

    rt_range = range(minimum(rt_raw), maximum(rt_raw), length=200)
    rt_fit = [Float64(_rt_bias_da(model, Float32(x))) * 1e3
              for x in rt_range]
    p_rt_bias = plot(size=(600, 600), bottommargin=10Plots.mm, leftmargin=10Plots.mm)
    scatter!(p_rt_bias, rt_raw, mz_int_resid_mda,
             alpha=0.08, markersize=1, color=:gray, label=nothing)
    scatter!(p_rt_bias, bin_rt_c, bin_med_rt_resid,
             xlabel="RT (min)", ylabel="Residual after m/z + intensity bias (mDa)",
             label="bin median", marker=:circle, color=:blue, ms=4)
    plot!(p_rt_bias, rt_range, rt_fit, lw=2.5, color=:red, label="RT bias spline")
    hline!(p_rt_bias, [0.0], color=:black, lw=1, ls=:dot, label=nothing)
    vline!(p_rt_bias, [Float64(model.rt_bias_extrap.lo), Float64(model.rt_bias_extrap.hi)],
           color=:orange, lw=1.5, ls=:dot, label="extrap boundary")
    title!(p_rt_bias, _split_title(fname, "RT-dependent bias (mDa)") *
                      "\nRT $(round(model.rt_min, digits=2))-$(round(model.rt_max, digits=2)) min")
    push!(plots_out, p_rt_bias)

    # ── Shared spread computation ──
    spread_bin_sz = TUNING_CALIBRATION_BIN_SIZE
    k_sel = Float64(model.k)
    scatter_idx = 1:n  # plot all points — PNG rendering handles size

    mz_f64 = Float64.(frag_mzs)

    # ── Plot 2: Spread (mDa) vs log2(I) — single panel ──
    I_order = sortperm(log2I)
    I_centers_mda, I_b_mda = binned_laplace_scale(collect(1:n)[I_order], log2I, corrected_mda, spread_bin_sz)
    log2I_range_local = range(minimum(log2I), maximum(log2I), length=200)
    spread_fit_mda = [Float64(_laplace_scale_mda(model, Float32(x))) for x in log2I_range_local]

    p_spread = plot(size=(600, 600),
        plot_title=_split_title(fname, "Spread (mDa) vs log2(I)"))
    scatter!(p_spread, I_centers_mda, I_b_mda,
        xlabel="log2(intensity)", ylabel="MAD-based spread (mDa)",
        label="b per bin", marker=:circle, color=:steelblue, ms=4)
    plot!(p_spread, log2I_range_local, spread_fit_mda, lw=2.5, color=:green,
        label="spline (mDa)")
    push!(plots_out, p_spread)

    # ── Plot 3: Envelope (mDa) vs log2(I) — single panel with Gaussian σ boundary ──
    eval_I = collect(log2I_range_local)
    laplace_to_gauss = log(2.0) / 0.6744897501960817
    fit_σ_mda = [Float64(_laplace_scale_mda(model, Float32(x))) * laplace_to_gauss for x in eval_I]

    # Compute fraction of fragments within k×σ boundary
    n_within_env = 0
    for i in 1:n
        σ_i = Float64(_laplace_scale_mda(model, Float32(log2I[i]))) * laplace_to_gauss *
              Float64(_mz_spread_correction(model, Float32(mz_f64[i])))
        if abs(corrected_mda[i]) <= k_sel * σ_i
            n_within_env += 1
        end
    end
    frac_within_env = round(100.0 * n_within_env / n, digits=1)

    n_scatter = length(scatter_idx)
    p_envelope = plot(size=(600, 600),
        plot_title=_split_title(fname, "Spread envelope") * "\n(n=$n, $(frac_within_env)% within k×σ boundary)")

    scatter!(p_envelope, log2I[scatter_idx], corrected_mda[scatter_idx],
        alpha=0.08, markersize=1, color=:steelblue, label=nothing,
        xlabel="log2(Intensity)", ylabel="Corrected error (mDa)")
    hline!(p_envelope, [0.0], color=:black, lw=1, ls=:dash, label=nothing)
    plot!(p_envelope, eval_I, k_sel .* fit_σ_mda, lw=2, ls=:dash, color=:green,
          label="k·σ boundary (k=$(round(k_sel, digits=2)), $(round(100*erf(k_sel/sqrt(2)), digits=1))% Gaussian coverage)")
    plot!(p_envelope, eval_I, -k_sel .* fit_σ_mda, lw=2, ls=:dash, color=:green, label=nothing)

    # Intensity quantile lines
    q_pcts = [10, 25, 50, 75, 90]
    q_vals = quantile(log2I, [0.10, 0.25, 0.50, 0.75, 0.90])
    for (pct, qv) in zip(q_pcts, q_vals)
        vline!(p_envelope, [qv]; lw=1, ls=:dot, color=:black, label=nothing)
        annotate!(p_envelope, [(qv, k_sel * fit_σ_mda[end] * 0.9, Plots.text("$(pct)%", 7, :center, :bottom, :black))])
    end
    push!(plots_out, p_envelope)

    # ── Plot 4: m/z-dependent spread correction (top panel only) ──
    mz_bin_sz = TUNING_CALIBRATION_BIN_SIZE
    mz_sort = sortperm(mz_f64)
    n_mz_bins = max(n ÷ mz_bin_sz, 1)

    bin_mz_centers = Float64[]
    bin_k_half     = Float64[]

    for b in 1:n_mz_bins
        s = (b - 1) * mz_bin_sz + 1
        e = b == n_mz_bins ? n : b * mz_bin_sz
        idx = mz_sort[s:e]
        r = abs.(corrected_mda[idx]) ./ max.(laplace_scale_mda[idx], 1e-4)
        push!(bin_mz_centers, median(mz_f64[idx]))
        push!(bin_k_half, median(r))
    end

    bin_correction = bin_k_half ./ log(2)
    n_mz = length(bin_mz_centers)
    X_mz = hcat(ones(n_mz), bin_mz_centers)
    coef_corr = X_mz \ bin_correction
    mz_fit_range = range(minimum(bin_mz_centers), maximum(bin_mz_centers), length=200)
    fit_corr = coef_corr[1] .+ coef_corr[2] .* mz_fit_range

    model_corr = [Float64(_mz_spread_correction(model, Float32(mz))) for mz in mz_fit_range]
    model_label = "model: correction spline ($(length(model.mz_spread_spline.coeffs)) coeffs)"

    p_corr = plot(size=(600, 600),
        plot_title=_split_title(fname, "m/z-dependent spread correction"),
        bottommargin=10Plots.mm, leftmargin=10Plots.mm)
    scatter!(p_corr, bin_mz_centers, bin_correction;
        marker=:circle, ms=4, color=:steelblue, label="c(mz) = k_half / ln(2)",
        xlabel="Fragment m/z", ylabel="Correction factor c(mz)")
    plot!(p_corr, mz_fit_range, fit_corr; lw=1.5, ls=:dot, color=:gray,
        label="OLS: $(round(coef_corr[1], digits=3)) + $(round(coef_corr[2], sigdigits=2))·mz")
    plot!(p_corr, mz_fit_range, model_corr; lw=2.5, color=:green, label=model_label)
    hline!(p_corr, [1.0]; lw=1.5, ls=:dash, color=:black, label="ideal (1.0)")
    push!(plots_out, p_corr)

    # ── Plot 5: fitted tolerance surface over fragment m/z and intensity ──
    mz_heat, log2I_heat, tolerance_heat_mda =
        _model_tolerance_heatmap_grid(model, frag_mzs, log2I)
    p_tolerance_heatmap = heatmap(
        mz_heat,
        log2I_heat,
        tolerance_heat_mda;
        xlabel="Fragment m/z",
        ylabel="log2(intensity)",
        title=_split_title(fname, "Fitted tolerance heatmap") *
              "\ncolor = k*sigma tolerance half-width (mDa)",
        size=(600, 600),
        color=:viridis,
        rightmargin=12Plots.mm,
        bottommargin=10Plots.mm,
        leftmargin=10Plots.mm,
    )
    push!(plots_out, p_tolerance_heatmap)

    return plots_out
end

# Fallback for SimpleMassErrorModel: no intensity plots
generate_intensity_model_plots(::NamedTuple, ::SimpleMassErrorModel, ::String) = Plots.Plot[]

"""
    generate_fallback_rt_plot_in_memory(results, fname, search_context, ms_file_idx)

Generates a diagnostic RT plot when no data is available (fallback/borrowed parameters).
Returns the plot object without saving to disk.
"""
function generate_fallback_rt_plot_in_memory(
    results::ParameterTuningSearchResults,
    fname::String,
    search_context::SearchContext,
    ms_file_idx::Int64
)
    # Create an informative plot showing the fallback/borrowed status
    p = Plots.plot(
        title = fname * "\n⚠️ Using Fallback/Borrowed Parameters",
        xlabel = "Retention Time RT (min)",
        ylabel = "Indexed Retention Time iRT (min)",
        size = (900, 900),
        grid = true
    )
    
    # Add identity line if using identity model
    rt_model = getRtToIrtModel(results)
    if isa(rt_model, IdentityModel)
        rt_range = [0, 120]  # Default RT range for visualization
        Plots.plot!(rt_range, rt_range, 
                   lw=2, ls=:dash, color=:red,
                   label="Identity Model (Fallback)")
    else
        # If borrowed, show the borrowed model
        rt_range = LinRange(0, 120, 100)
        Plots.plot!(rt_range, rt_model.(rt_range),
                   lw=2, color=:blue,
                   label="Borrowed RT Model")
    end
    
    # Add text annotation explaining the situation
    Plots.annotate!(60, 20, 
                   text("Insufficient PSMs for RT calibration\n" *
                        "Using conservative parameters", 
                        :center, 10))
    
    return p  # Return plot for collection
end

"""
    generate_fallback_mass_error_plot_in_memory(results, fname, search_context, ms_file_idx)

Generates a diagnostic mass error plot when no data is available (fallback/borrowed parameters).
Returns the plot object without saving to disk.
"""
function generate_fallback_mass_error_plot_in_memory(
    results::ParameterTuningSearchResults,
    fname::String,
    search_context::SearchContext,
    ms_file_idx::Int64
)
    mem = getMassErrorModel(results)
    mass_offset = getMassOffset(mem)
    left_tol = getLeftTol(mem)
    right_tol = getRightTol(mem)
    
    # Create an informative plot showing the tolerances being used
    p = Plots.plot(
        title = fname * "\n⚠️ Using Fallback/Borrowed Parameters",
        xlabel = "Count",
        ylabel = "Mass Error (ppm)",
        size = (900, 900),
        ylim = (-max(left_tol, right_tol) - 10, max(left_tol, right_tol) + 10),
        grid = true,
        yflip = true,
        orientation = :h
    )
    
    # Draw tolerance boundaries
    Plots.hline!([mass_offset], color=:black, lw=2, label="Mass Offset: $(round(mass_offset, digits=2)) ppm")
    Plots.hline!([mass_offset - left_tol], color=:red, lw=1, ls=:dash, label="Lower Bound: -$(round(left_tol, digits=1)) ppm")
    Plots.hline!([mass_offset + right_tol], color=:red, lw=1, ls=:dash, label="Upper Bound: +$(round(right_tol, digits=1)) ppm")
    
    # Add text annotation
    Plots.annotate!(0, 0,
                   text("Insufficient fragment matches\n" *
                        "Tolerance: ±$(round((left_tol+right_tol)/2, digits=1)) ppm",
                        :center, 10))
    
    return p  # Return plot for collection
end

# Legacy functions that save to disk (kept for backward compatibility)
function generate_fallback_rt_plot(
    results::ParameterTuningSearchResults,
    plot_path::String,
    fname::String,
    search_context::SearchContext,
    ms_file_idx::Int64
)
    p = generate_fallback_rt_plot_in_memory(results, fname, search_context, ms_file_idx)
    savefig(p, plot_path)
    return p
end

function generate_fallback_mass_error_plot(
    results::ParameterTuningSearchResults,
    plot_path::String,
    fname::String,
    search_context::SearchContext,
    ms_file_idx::Int64
)
    p = generate_fallback_mass_error_plot_in_memory(results, fname, search_context, ms_file_idx)
    savefig(p, plot_path)
    return p
end







"""
    get_fallback_parameters(search_context, ms_file_idx)

Get fallback parameters from another file or use defaults.
Returns (mass_err_model, rt_model, borrowed_from_idx)
"""
function get_fallback_parameters(
    search_context::SearchContext,
    ms_file_idx::Int64
)

    # Default mass error model
    mass_err = MassErrorModel(0.0f0, (20.0f0, 20.0f0))
    
    # Default RT model (identity)
    rt_model = IdentityModel()
    
    return mass_err, rt_model, nothing
end

"""
    record_file_status!(diagnostics, ms_file_idx, file_name, converged, used_fallback, iteration_state)

Record the parameter tuning status for a file.
"""
function record_file_status!(
    diagnostics::ParameterTuningDiagnostics,
    ms_file_idx::Int64,
    file_name::String,
    converged::Bool,
    used_fallback::Bool,
    iteration_state::IterationState
)
    status = ParameterTuningStatus(
        ms_file_idx,
        file_name,
        converged,
        used_fallback,
        used_fallback ? "Failed convergence" : "",
        iteration_state.scan_attempt,
        0,  # PSM count - could be tracked if needed
        0.0f0,  # Mass offset - could be extracted if needed
        (0.0f0, 0.0f0)  # Tolerances - could be extracted if needed
    )
    
    diagnostics.file_statuses[ms_file_idx] = status
    
    if converged && !used_fallback
        diagnostics.n_successful += 1
    elseif used_fallback
        diagnostics.n_fallback += 1
    else
        diagnostics.n_failed += 1
    end
end

"""
    record_tuning_status!(diagnostics, status)

Record the parameter tuning status for a file.
"""
function record_tuning_status!(diagnostics::ParameterTuningDiagnostics, status::ParameterTuningStatus)
    diagnostics.file_statuses[status.file_idx] = status
    
    if status.converged && !status.used_fallback
        diagnostics.n_successful += 1
    elseif status.used_fallback
        diagnostics.n_fallback += 1
    else
        diagnostics.n_failed += 1
    end
end

"""
    store_tuning_results!(history::ParameterHistory, file_idx::Int64, results::TuningResults)

Store the parameter tuning results for a file.
"""
function store_tuning_results!(history::ParameterHistory, file_idx::Int64, results::TuningResults)
    history.file_parameters[file_idx] = results
    
    # Update global statistics if this was successful
    if results.converged
        update_global_statistics!(history)
    end
end

"""
    update_global_statistics!(history::ParameterHistory)

Recompute global statistics from all successful files.
"""
function update_global_statistics!(history::ParameterHistory)
    successful_results = [r for r in values(history.file_parameters) if r.converged]
    
    if isempty(successful_results)
        return
    end
    
    # Extract values
    mass_offsets = [r.mass_offset for r in successful_results]
    tolerances = [mean(r.mass_tolerance) for r in successful_results]
    
    # Compute statistics
    stats = history.global_stats
    stats.median_mass_offset = median(mass_offsets)
    stats.mass_offset_mad = mad(mass_offsets, normalize=false)
    stats.median_tolerance = median(tolerances)
    stats.tolerance_mad = mad(tolerances, normalize=false)
    stats.n_successful_files = length(successful_results)
end
