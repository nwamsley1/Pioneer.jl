
using Logging
using DataStructures
using Statistics

# Configuration types with validation
struct PreSearchConfig 
    sample_rate::Float32
    max_presearch_iters::Int
    frag_tol_ppm::Float32
    max_qval::Float32
    min_samples::Int
    frag_err_quantile::Float32

    function PreSearchConfig(sample_rate, max_iters, tol_ppm, qval, min_samples, err_quantile)
        0 < sample_rate ≤ 1 || throw(ArgumentError("Sample rate must be between 0 and 1"))
        max_iters > 0 || throw(ArgumentError("Max iterations must be positive"))
        tol_ppm > 0 || throw(ArgumentError("Fragment tolerance must be positive"))
        
        new(sample_rate, max_iters, tol_ppm, qval, min_samples, err_quantile)
    end
end

struct IRTConfig
    n_sigma_tol::Float32
    spline_degree::Int
    spline_knots::Int
    outlier_threshold::Float32

    function IRTConfig(n_sigma_tol, degree=3, knots=5, thresh=10.0f0)
        n_sigma_tol > 0 || throw(ArgumentError("Sigma tolerance must be positive"))
        1 ≤ degree ≤ 5 || throw(ArgumentError("Spline degree must be between 1 and 5"))
        knots ≥ 3 || throw(ArgumentError("Must have at least 3 knots"))
        
        new(n_sigma_tol, degree, knots, thresh)
    end
end

# Core context and results types
struct SearchContext
    ms_table_paths::Vector{String}
    spec_lib::Dict{String,Any}
    ion_matches::Vector{FragmentMatch{Float32}}
    ion_misses::Vector{FragmentMatch{Float32}}
    id_to_col::Vector{ArrayDict{UInt32, UInt16}}
    iso_splines::IsotopeSplineModel
    output_dirs::NamedTuple{(:rt_align, :mass_err),Tuple{String,String}}
end

mutable struct SearchResults
    rt_to_irt_map::Dict{Int64,Any}
    frag_err_dist::Dict{Int64,MassErrorModel}
    irt_errs::Dict{Int64,Float64}
    metrics::Dict{String,Any}
    cache::Dict{String,Any}
end

# Custom error types
struct PreSearchError <: Exception
    file::String
    message::String
end

struct IRTError <: Exception
    file::String
    message::String
end

# Main search implementation
struct ParameterTuningSearch <: TuningStrategy
    presearch_config::PreSearchConfig
    irt_config::IRTConfig
    use_cache::Bool
end

# Core algorithm implementations using multiple dispatch
function fit_irt_model(::Type{ParameterTuningSearch}, psms::DataFrame, config::IRTConfig)
    # Initial fit
    rt_to_irt_map = UniformSpline(
        psms.irt_predicted,
        psms.rt,
        config.spline_degree,
        config.spline_knots
    )
    
    # Calculate residuals
    psms.irt_observed = rt_to_irt_map.(psms.rt)
    residuals = psms.irt_observed .- psms.irt_predicted
    irt_mad = mad(residuals)
    
    # Remove outliers and refit
    valid_psms = psms[abs.(residuals) .< (irt_mad * config.outlier_threshold), :]
    
    final_map = UniformSpline(
        valid_psms.irt_predicted,
        valid_psms.rt,
        config.spline_degree,
        config.spline_knots
    )
    
    return final_map, irt_mad
end

function execute_presearch!(
    strategy::ParameterTuningSearch,
    context::SearchContext,
    ms_file_idx::Int,
    ms_table::Arrow.Table
)::DataFrame
    
    mass_err_model = MassErrorModel(
        0.0f0,
        (strategy.presearch_config.frag_tol_ppm, strategy.presearch_config.frag_tol_ppm)
    )

    psms = DataFrame()
    n_attempts = 0
    
    while n_attempts ≤ strategy.presearch_config.max_presearch_iters
        try
            results = LibrarySearch(
                ms_table,
                context.spec_lib;
                fragment_lookup_table = context.spec_lib["f_det"],
                rt_to_irt_spline = identity,
                ms_file_idx = UInt32(ms_file_idx),
                mass_err_model = mass_err_model,
                sample_rate = strategy.presearch_config.sample_rate
            )
            
            new_psms = process_search_results(results)
            psms = isempty(psms) ? new_psms : vcat(psms, new_psms)
            
            if check_convergence(psms, strategy.presearch_config)
                break
            end
        catch e
            @warn "Presearch attempt failed" attempt=n_attempts exception=e
        end
        
        n_attempts += 1
    end
    
    validate_presearch_results!(psms, strategy.presearch_config) || 
        throw(PreSearchError("Insufficient high-quality PSMs found"))
        
    return psms
end

function estimate_mass_errors(
    ::Type{ParameterTuningSearch},
    fragments::Vector{<:MatchIon},
    config::PreSearchConfig,
    output_dir::String,
    file_name::String
)
    ppm_errors = [calculate_ppm_error(frag) for frag in fragments]
    
    return ModelMassErrs(
        ppm_errors;
        frag_err_quantile = config.frag_err_quantile,
        out_fdir = output_dir,
        out_fname = file_name
    )
end

# Main execution function
function execute_search(
    strategy::ParameterTuningSearch,
    context::SearchContext
)::SearchResults
    
    @info "Starting parameter tuning search" n_files=length(context.ms_table_paths)
    
    results = SearchResults(
        Dict{Int64,Any}(),
        Dict{Int64,MassErrorModel}(),
        Dict{Int64,Float64}(),
        Dict{String,Any}(),
        Dict{String,Any}()
    )
    
    clear_output_directories!(context.output_dirs)
    
    Threads.@threads for (idx, path) in collect(enumerate(context.ms_table_paths))
        process_file!(strategy, context, idx, path, results)
    end
    
    generate_qc_plots!(results, context.output_dirs)
    
    return results
end

# Helper function for processing individual files
function process_file!(
    strategy::ParameterTuningSearch,
    context::SearchContext,
    file_idx::Int,
    file_path::String,
    results::SearchResults
)
    @debug "Processing file" file_path
    
    if strategy.use_cache && haskey(results.cache, file_path)
        @debug "Using cached results" file_path
        update_results_from_cache!(results, file_idx, file_path)
        return
    end
    
    try
        ms_table = Arrow.Table(file_path)
        psms = execute_presearch!(strategy, context, file_idx, ms_table)
        
        rt_model, irt_mad = fit_irt_model(
            ParameterTuningSearch,
            psms,
            strategy.irt_config
        )
        
        fragments = get_matched_fragments(context, psms)
        mass_errors = estimate_mass_errors(
            ParameterTuningSearch,
            fragments,
            strategy.presearch_config,
            context.output_dirs.mass_err,
            basename(file_path)
        )
        
        # Update results atomically
        lock(results) do
            results.rt_to_irt_map[file_idx] = rt_model
            results.frag_err_dist[file_idx] = mass_errors
            results.irt_errs[file_idx] = strategy.irt_config.n_sigma_tol * irt_mad
            
            if strategy.use_cache
                results.cache[file_path] = (rt_model, mass_errors, irt_mad)
            end
        end
        
    catch e
        handle_processing_error(e, file_path, file_idx, results)
    end
end

# Usage example:
function run_parameter_tuning(
    ms_paths::Vector{String},
    spec_lib::Dict,
    output_dirs::NamedTuple,
    ion_data::NamedTuple;
    sample_rate::Float32 = 0.1f0,
    max_iters::Int = 10
)
    presearch_config = PreSearchConfig(
        sample_rate,
        max_iters,
        20.0f0,  # frag_tol_ppm
        0.01f0,  # max_qval
        1000,    # min_samples
        0.95f0   # frag_err_quantile
    )
    
    irt_config = IRTConfig(3.0f0)
    
    strategy = ParameterTuningSearch(presearch_config, irt_config, true)
    
    context = SearchContext(
        ms_paths,
        spec_lib,
        ion_data.matches,
        ion_data.misses,
        ion_data.id_to_col,
        ion_data.iso_splines,
        output_dirs
    )
    
    return execute_search(strategy, context)
end