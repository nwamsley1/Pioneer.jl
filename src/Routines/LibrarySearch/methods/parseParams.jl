function parseParams(
    params::Dict{String, Any}
)

return (
    isotope_err_bounds = Tuple([Int64(bound) for bound in params["isotope_err_bounds"]]),
    choose_most_intense = Bool(params["choose_most_intense"]),
    presearch_params = Dict{String, Any}(k => v for (k, v) in params["presearch_params"]),
    first_search_params = Dict{String, Any}(k => v for (k, v) in params["first_search_params"]),
    quant_search_params = Dict{String, Any}(k => v for (k, v) in params["quant_search_params"]),
    irt_mapping_params = Dict{String, Any}(k => v for (k, v) in params["irt_mapping_params"]),
    deconvolution_params = Dict{String, Any}(k => v for (k, v) in params["deconvolution_params"]),
    summarize_first_search_params = Dict{String, Any}(k => v for (k, v) in params["summarize_first_search_params"]),
    qc_plot_params = Dict{String, Any}(k => v for (k, v) in params["qc_plot_params"]),
    normalization_params = Dict{String, Any}(k => v for (k, v) in params["normalization_params"]),
    benchmark_params = Dict{String, Any}(k => v for (k, v) in params["benchmark_params"])
    );
end