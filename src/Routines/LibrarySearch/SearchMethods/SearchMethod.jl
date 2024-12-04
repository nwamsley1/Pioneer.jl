#This method should be generic 
function execute_search(
    search_type::SearchMethod, 
    search_context::SearchContext,
    params::Any)

    msdr = getMassSpecData(search_context)

    @info "Starting parameter tuning search" n_files=length(msdr.file_paths)

    search_parameters = get_parameters(search_type, params)

    for (ms_file_idx, spectra) in ProgressBar(enumerate(msdr))

        #Initialize results
        search_results = init_search_results(search_parameters, search_context, ms_file_idx)

        #Analysis of MS data 
        process_file!(search_results, search_parameters, search_context, ms_file_idx, spectra)

        #Modifies `search_results` and/or `search_context` in place 
        process_search_results!(search_results, search_parameters, search_context, ms_file_idx)

        #Reset search results for the next file
        reset_results!(search_results)
    end

    summarize_results!(search_results, search_parameters, search_context)
return search_results
#error("execute_search not implemented for $(typeof(strategy))")
end

function get_parameters(search_type::SearchMethod, params::Any)
 error("get_parameters not implemented for search method of type $(typeof(search_type))")
end
function init_search_results(search_parameters::SearchParameters, search_context::SearchContext, ms_file_idx::Int64) 
error("init_search_results not implemented for params of type $(typeof(search_parameters))")
end
function process_file!(results::SearchResults, params::SearchParameters, search_context::SearchContext, ms_file_idx::Int64, spectra::Arrow.Table)
error("process_file! not implemented for params of type $(typeof(paarams)) and results of type$(typeof(results)) ")
end
function process_file!(results::SearchResults, params::SearchParameters, search_context::SearchContext, ms_file_idx::Int64, spectra::Arrow.Table)
error("process_file! not implemented for params of type $(typeof(paarams)) and results of type$(typeof(results)) ")
end
function summarize_results!(results::SearchResults, params::SearchParameters, search_context::SearchContext)
    error("summarize_results! not implemented for params of type $(typeof(paarams)) and results of type$(typeof(results)) ")

end
   
# Helper functions
function partition_scans(ms_table, n_threads)
    thread_tasks, total_peaks = partitionScansToThreads(
        ms_table[:mz_array],
        ms_table[:retentionTime],
        ms_table[:centerMz],
        ms_table[:msOrder],
        n_threads,
        1
    )
    return thread_tasks
end

