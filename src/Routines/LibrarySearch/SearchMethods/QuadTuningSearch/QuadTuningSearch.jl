"""
    QuadTuningSearch

Search method for optimizing quadrupole transmission models.

This search:
1. Collects PSMs with extended precursor isotope patterns
2. Performs deconvolution to estimate relative isotope abundances
3. Fits transmission model based on isotope ratio deviations
4. Stores optimized models in SearchContext for other methods

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :isotope_err_bounds => (0, 0),  # Fixed for quad tuning
    :presearch_params => Dict(
        "frag_tol_ppm" => 30.0,
        "min_index_search_score" => 3,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "n_frag_isotopes" => 2,
        "max_frag_rank" => 10,
        "quad_tuning_sample_rate" => 0.1,
        "min_quad_tuning_fragments" => 3,
        "min_quad_tuning_psms" => 1000,
        "abreviate_precursor_calc" => false
    ),
    :deconvolution_params => Dict(
        "max_iter_newton" => 100,
        "max_iter_bisection" => 100,
        "max_iter_outer" => 100,
        "accuracy_newton" => 1e-5,
        "accuracy_bisection" => 1e-4,
        "max_diff" => 1e-5
    )
)

# Execute search
results = execute_search(QuadTuningSearch(), search_context, params)
```
"""
struct QuadTuningSearch <: TuningMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for quadrupole tuning search.
"""
struct QuadTuningSearchResults <: SearchResults
    tuning_results::Vector{Vector{@NamedTuple{
        precursor_idx::UInt32,
        scan_idx::UInt32,
        weight::Float32,
        iso_idx::UInt8,
        center_mz::Float32,
        n_matches::UInt8
    }}}
    quad_model::Base.Ref{QuadTransmissionModel}
    quad_plot_dir::String
end

"""
Parameters for quadrupole tuning search.
Configures deconvolution and quad model fitting.
"""
struct QuadTuningSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # Search parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    frag_tol_ppm::Float32
    min_index_search_score::UInt8
    min_log2_matched_ratio::Float32
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::UInt8
    max_frag_rank::UInt8
    n_frag_isotopes::Int64
    sample_rate::Float32
    irt_tol::Float32
    spec_order::Set{Int64}
    
    # Deconvolution parameters
    max_iter_newton::Int64
    max_iter_bisection::Int64
    max_iter_outer::Int64
    accuracy_newton::Float32
    accuracy_bisection::Float32
    max_diff::Float32
    
    # Quad tuning specific parameters
    min_quad_tuning_fragments::Int64
    min_quad_tuning_psms::Int64
    prec_estimation::P

    function QuadTuningSearchParameters(params::Any)
        pp = params[:presearch_params]
        dp = params[:deconvolution_params]
        prec_estimation = PartialPrecCapture() #pp["abreviate_precursor_calc"] ? FullPrecCapture() : PartialPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(0), UInt8(0)),  # Fixed for quad tuning
            Float32(pp["frag_tol_ppm"]),
            UInt8(pp["min_index_search_score"]),
            typemin(Float32),
            Int64(pp["min_frag_count"]),
            Float32(pp["min_spectral_contrast"]),
            (Int64(first(pp["min_topn_of_m"])), Int64(last(pp["min_topn_of_m"]))),
            UInt8(pp["max_best_rank"]),
            UInt8(pp["max_frag_rank"]),
            one(Int64),
            Float32(pp["quad_tuning_sample_rate"]),
            typemax(Float32),
            Set(2),

            Int64(dp["max_iter_newton"]),
            Int64(dp["max_iter_bisection"]),
            Int64(dp["max_iter_outer"]),
            Float32(dp["accuracy_newton"]),
            Float32(dp["accuracy_bisection"]),
            Float32(dp["max_diff"]),

            Int64(pp["min_quad_tuning_fragments"]),
            Int64(pp["min_quad_tuning_psms"]),
            prec_estimation
        )
    end
end

#==========================================================
Interface Implementation
==========================================================#

# Getters
getQuadModel(q::QuadTuningSearchResults) = q.quad_model[]

# Setters
function setQuadModel(q::QuadTuningSearchResults, model::Q) where {Q<:QuadTransmissionModel}
    q.quad_model[] = model
end

get_parameters(::QuadTuningSearch, params::Any) = QuadTuningSearchParameters(params)

function init_search_results(::QuadTuningSearchParameters, search_context::SearchContext)
    # Initialize empty tuning results vector in each search data structure
    out_dir = getDataOutDir(search_context)
    qc_dir = joinpath(out_dir, "qc_plots")
    !isdir(qc_dir) && mkdir(qc_dir)

    qpp = joinpath(qc_dir, "quad_transmission_model")
    !isdir(qpp) && mkdir(qpp)
    quad_data = joinpath(qpp, "quad_data")
    !isdir(quad_data) && mkdir(quad_data)
    quad_models = joinpath(qpp, "quad_models")
    !isdir(quad_models) && mkdir(quad_models)
    temp_data = Vector{Vector{@NamedTuple{
            precursor_idx::UInt32,
            scan_idx::UInt32,
            weight::Float32,
            iso_idx::UInt8,
            center_mz::Float32,
            n_matches::UInt8
        }}}()#(undef, length(getSearchData(search_context)))
    for i in range(1, length(getSearchData(search_context)))
        push!(temp_data, Vector{@NamedTuple{
            precursor_idx::UInt32,
            scan_idx::UInt32,
            weight::Float32,
            iso_idx::UInt8,
            center_mz::Float32,
            n_matches::UInt8
        }}())
    end
    return QuadTuningSearchResults(
        temp_data,
        Ref{QuadTransmissionModel}(),
        qpp
    )
end

#==========================================================
Core Processing Methods
==========================================================#

"""
Main file processing method for quad tuning search.
"""
function process_file!(
    results::QuadTuningSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::Arrow.Table) where {P<:QuadTuningSearchParameters}

    setQuadTransmissionModel!(search_context, ms_file_idx, SquareQuadModel(1.0f0))
    
    try
        setNceModel!(
            getFragmentLookupTable(getSpecLib(search_context)), 
            getNceModelModel(search_context, ms_file_idx)
        )
        # Adjust arrays for isotope variants
        adjust_precursor_arrays!(search_context)
        
        # Check window widths
        window_widths = check_window_widths(spectra)
        if length(window_widths) != 1
            @warn "Multiple window sizes detected: $(join(collect(window_widths), ';'))"
            setQuadModel(results, GeneralGaussModel(5.0f0, 0.0f0))
            return results
        end
        
        # Collect and process PSMs
        total_psms = collect_psms(spectra, search_context, results, params, ms_file_idx)
        
        if nrow(total_psms) == 0
            setQuadModel(results, GeneralGaussModel(5.0f0, 0.0f0))
            return results
        end

        # Plot charge states
        plot_charge_distributions(total_psms, results, getFileIdToName(getMSData(search_context), ms_file_idx))
        
        # Fit quad model
        window_width = parse(Float64, first(window_widths))

        fitted_model = fit_quad_model(total_psms, window_width)
        setQuadModel(results, RazoQuadModel(fitted_model))
        
    catch e
        throw(e)
        @warn "Quad transmission function fit failed" exception=e
        setQuadModel(results, GeneralGaussModel(5.0f0, 0.0f0))
    end
    
    return results
end

function process_search_results!(
    results::QuadTuningSearchResults,
    ::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    ::Arrow.Table
) where {P<:QuadTuningSearchParameters}
    
    setQuadTransmissionModel!(search_context, ms_file_idx, getQuadModel(results))
end

function summarize_results!(
    results::QuadTuningSearchResults,
    ::P,
    search_context::SearchContext
) where {P<:QuadTuningSearchParameters}
    
    plot_bins = LinRange(0-3, 0+3, 100)
    for (ms_file_idx, quad_model) in pairs(search_context.quad_transmission_model)
        fname = getFileIdToName(getMSData(search_context), ms_file_idx)
        quad_func = getQuadTransmissionFunction(quad_model, 0.0f0, 2.0f0)
        p = plot(plot_bins, quad_func.(plot_bins), lw = 2, alpha = 0.5, title = "$fname")
        savefig(p, joinpath(results.quad_plot_dir, "quad_models", fname*".pdf"))
    end

    qmp = [x for x in readdir(joinpath(results.quad_plot_dir, "quad_models"), join=true) if endswith(x, ".pdf")]
    if !isempty(qmp)
        merge_pdfs(qmp, 
                  joinpath(results.quad_plot_dir, "quad_models", "quad_model_plots.pdf"), 
                  cleanup=true)
    end

    qmp = [x for x in readdir(joinpath(results.quad_plot_dir, "quad_data"), join=true) if endswith(x, ".pdf")]
    if !isempty(qmp)
        merge_pdfs(qmp, 
                  joinpath(results.quad_plot_dir, "quad_data", "quad_data_plots.pdf"), 
                  cleanup=true)
    end

    reset_precursor_arrays!(search_context)
    return nothing
end

function reset_results!(results::QuadTuningSearchResults)
    for r in results.tuning_results
        empty!(r)
    end
    return nothing
end
