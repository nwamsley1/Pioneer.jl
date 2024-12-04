# Define abstract types and traits
abstract type SearchMethod end
abstract type TuningMethod <: SearchMethod end

abstract type SearchResults end 
include("src/Routines/LibrarySearch/selectTransitions/selectTransitions.jl")

# Common interface for containers of pre-allocated
# data structures used for conducting searches. Does not hold
# the search results, rather, intermediate data used during the search. 
abstract type SearchDataStructures end
getIonMatches(s::SearchDataStructures) = s.ion_matches
getIonMisses(s::SearchDataStructures) = s.ion_misses
getIdToCol(s::SearchDataStructures) = s.id_to_col
getPrecursorScores(s::SearchDataStructures) = s.prec_count
getIonTemplates(s::SearchDataStructures) = s.ion_templates
getIsoSplines(s::SearchDataStructures) = s.iso_splines
getScoredPsms(s::SearchDataStructures) = s.scored_psms
getUnscoredPsms(s::SearchDataStructures) = s.unscored_psms
getSpectralScores(s::SearchDataStructures) = s.spectral_scores

struct SimpleLibrarySearch{I<:IsotopeSplineModel} <: SearchDataStructures
    ion_matches::Vector{FragmentMatch{Float32}}
    ion_misses::Vector{FragmentMatch{Float32}}
    id_to_col::ArrayDict{UInt32, UInt16}
    prec_count::Counter{UInt32, UInt8}
    ion_templates::Vector{DetailedFrag{Float32}}
    iso_splines::I
    scored_psms::Vector{SimpleScoredPSM{Float32, Float16}}
    unscored_psms::Vector{SimpleUnscoredPSM{Float32}}
    spectral_scores::Vector{SpectralScoresSimple{Float16}}
end

#Common interface for references to mass spec data used in searches
abstract type MassSpecDataReference end
getMSData(msdr::MassSpecDataReference, ms_file_idx::I) where {I<:Integer} = Arrow.Table(msdr.file_paths[ms_file_idx])
struct ArrowTableReference <: MassSpecDataReference
    file_paths::Vector{AbstractString}
end
import Base: enumerate
function enumerate(msdr::ArrowTableReference)
    return zip(1:length(msdr.file_paths), (getMSData(msdr, i) for i in 1:length(msdr.file_paths)))
end

# Common interface for containers for search data
struct SearchContext{N,L<:FragmentIndexLibrary,S<:SearchDataStructures,M<:MassSpecDataReference}
    spec_lib::L
    temp_structures::AbstractVector{S}
    mass_spec_data_reference::M
    data_out_dir::Base.Ref{String}
    quad_transmission_model::Dict{Int64, QuadTransmissionModel}
    mass_error_model::Dict{Int64, MassErrorModel}
    rt_to_irt_model::Dict{Int64, Any}

    function SearchContext(
        spec_lib::L,
        temp_structures::AbstractVector{S},
        mass_spec_data_reference::M
    ) where {L<:FragmentIndexLibrary,S<:SearchDataStructures,M<:MassSpecDataReference}
        N = length(temp_structures)
        new{N,L,S,M}(
            spec_lib,
            temp_structures,
            mass_spec_data_reference,
            Ref{Sring}()
            Dict{Int64, QuadTransmissionModel}(),
            Dict{Int64, MassErrorModel}(),
            Dict{Int64, Any}()
        )
    end

end

getMassSpecData(s::SearchContext) = s.mass_spec_data_reference
getSpecLib(s::SearchContext) = s.spec_lib
getSearchData(s::SearchContext) = s.temp_structures

function getQuadTransmissionModel(s::SearchContext, index::Int64) 
    if haskey(s.quad_transmission_model, index)
        return s.quad_transmission_model[index]
    else
        #Return a sensible default 
        return GeneralGaussModel(5.0f0, 0.0f0)
    end
end

function setQuadTransmissionModel!(s::SearchContext, index::Int64, model::QuadTransmissionModel)
    s.quad_transmission_model[index] = model
end

function getMassErrorModel(s::SearchContext, index::Int64) 
    if haskey(s.quad_transmission_model, index)
        return s.quad_transmission_model[index]
    else
        #Return a sensible default 
        return MassErrorModel(zero(Float32), #(getFragTolPpm(params), getFragTolPpm(params)),
        (30.0f0, 30.0f0))
    end
end

function setMassErrorModel!(s::SearchContext, index::Int64, model::MassErrorModel)
    s.mass_error_model[index] = model
end

function getRtIrtModel(s::SearchContext, index::Int64) 
    if haskey(s.rt_to_irt_model, index)
        return s.rt_to_irt_model[index]
    else
        #Return a sensible default 
        return x::Float32->x::Float32
    end
end

function setRtIrtModel!(s::SearchContext, index::Int64, model::Any)
    s.rt_to_irt_model[index] = model
end

#Common interface for tunable search parameters
abstract type SearchParameters end 
getIsotopeErrBounds(fsp::SearchParameters) = fsp.isotope_err_bounds
getMinIndexSearchScore(fsp::SearchParameters) = fsp.min_index_search_score
getMinFragCount(fsp::SearchParameters) = fsp.min_frag_count
getMinSpectralContrast(fsp::SearchParameters) = fsp.min_spectral_contrast
getMinLog2MatchedRatio(fsp::SearchParameters) = fsp.min_log2_matched_ratio
getMinTopNofM(fsp::SearchParameters) = fsp.min_topn_of_m
getMaxBestRank(fsp::SearchParameters) = fsp.max_best_rank
getNFragIsotopes(fsp::SearchParameters) = fsp.n_frag_isotopes
getMaxFragRank(fsp::SearchParameters) = fsp.max_frag_rank
getSampleRate(fsp::SearchParameters) = fsp.sample_rate
getIRTTol(fsp::SearchParameters) = fsp.irt_tol
getSpecOrder(fsp::SearchParameters) = fsp.spec_order
getPrecEstimation(fsp::SearchParameters) = fsp.prec_estimation
getFragTolPpm(fsp::SearchParameters) = fsp.prec_estimation
getSplineDegree(fsp::SearchParameters) = fsp.spline_degree
getSplineNKnots(fsp::SearchParameters) = fsp.spline_n_knots
getOutlierThreshold(fsp::SearchParameters) = fsp.spline_n_knots

abstract type SpectralLibrary end
struct FragmentIndexLibrary <: SpectralLibrary
    presearch_fragment_index::FragmentIndex{Float32}
    fragment_index::FragmentIndex{Float32}
    precursors::Arrow.Table
    fragment_lookup_table::SplineFragmentLookup
end
getPresearchFragmentIndex(sl::SpectralLibrary) = sl.presearch_fragment_index
getFragmentIndex(sl::SpectralLibrary) = sl.fragment_index
getPrecursors(sl::SpectralLibrary) = sl.precursors
getFragmentLookupTable(sl::SpectralLibrary) = sl.fragment_lookup_table 


#=
struct FileSpecificModels{Q<:QuadTransmissionModel, M<:MassErrorModel, N<:NceModel, U<:UniformBasisCubicSpline}
    quad_transmission_model::Dict{Int64, Q}
    mass_err_model::Dict{Int64, M}
    nce_model::Dict{Int64, N}
    rt_to_irt_model::Dict{Int64, U}
end
=#
function initSimpleSearchContext(
    iso_splines::IsotopeSplineModel,
    n_precursors::Int64,
    M::Int64
    )
    SimpleLibrarySearch(
        [FragmentMatch{Float32}() for _ in range(1, M)],
        [FragmentMatch{Float32}() for _ in range(1, M)],
        ArrayDict(UInt32, UInt16, n_precursors),
        Counter(UInt32, UInt8,n_precursors ),
        [DetailedFrag{Float32}() for _ in range(1, M)],
        iso_splines,
        Vector{SimpleScoredPSM{Float32, Float16}}(undef, 5000),
        [SimpleUnscoredPSM{Float32}() for _ in range(1, 5000)],
        Vector{SpectralScoresSimple{Float16}}(undef, 5000)
    )
end
function initSimpleSearchContexts(
    iso_splines::IsotopeSplineModel,
    n_precursors::Int64,
    N::Int64,
    M::Int64)
    [initSimpleSearchContext(iso_splines, n_precursors, M) for _ in 1:N]
end

#This method should be generic 
function execute_search(
        search_type::SearchMethod, 
        search_context::SearchContext,
        params::Any)

        msdr = getMassSpecData(search_context)

        @info "Starting parameter tuning search" n_files=length(msdr.file_paths)
    
        search_parameters = get_parameters(search_type, params)

        for (ms_file_idx, spectra) in ProgressBar(enumerate(msdr))

            #Initialize results that will be returned. 
            #Some results are instaed written to a file specified in `search_context``: SearchContext.data_out_dir::Base.Ref{String}
            #And some results are simply modification of `search_context` data. 
            search_results = init_search_results(search_parameters, search_context, ms_file_idx)

            #Analysis of MS data 
            process_file!(search_results, search_parameters, search_context, ms_file_idx, spectra)

            #Modifies `search_results` and/or `search_context` in place 
            process_search_results!(search_results, search_parameters, search_context, ms_file_idx)

        end

    return search_results
    #error("execute_search not implemented for $(typeof(strategy))")
end




function execute_presearch!(
    search_type::ParameterTuningSearch, 
    ms_file_idx::Int64,
    spectra::Arrow.Table,
    spec_lib::SpectralLibrary,
    search_data::NTuple{N, S},
    qtm::Q,
    mem::M,
    rt_to_irt_spline::Any,
    params::P) where {N,  M<:MassErrorModel, 
    Q<:QuadTransmissionModel,
    S<:SearchDataStructures,
    P<:FragmentIndexSearchParameters}
        
        LibrarySearch(
            spectra,
            UInt32(ms_file_idx),
            getPresearchFragmentIndex(spec_lib),
            spec_lib,
            search_data,
            qtm,
            mem,
            rt_to_irt_spline,
            params,
        )
end




function process_results(strategy::SearchStrategy, context::SearchContext, results::Any)
    error("process_results not implemented for $(typeof(strategy))")
end

# Example implementation for ParameterTuningSearch
function prepare_search(::ParameterTuningSearch, context::MSSearchContext)
    # Clear directories, initialize models
    [rm(joinpath(rt_alignment_folder, x)) for x in readdir(rt_alignment_folder)]
    [rm(joinpath(mass_err_estimation_folder, x)) for x in readdir(mass_err_estimation_folder)]
    
    rt_to_irt_map_dict = Dict{Int64, Any}()
    frag_err_dist_dict = Dict{Int64, MassErrorModel}()
    irt_errs = Dict{Int64, Float64}()
    
    return (rt_to_irt_map_dict, frag_err_dist_dict, irt_errs)
end

# Common utilities
struct SearchResults
    rt_to_irt_map::Dict
    frag_err_dist::Dict
    irt_errs::Dict
    quad_model::Union{Dict, Nothing}
    nce_model::Union{Dict, Nothing}
    peak_fwhms::Union{Dict, Nothing}
    psms_paths::Union{Dict, Nothing}
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

# Unified search execution
function run_search(strategy::SearchStrategy, context::MSSearchContext)
    prepared_data = prepare_search(strategy, context)
    results = execute_search(strategy, context, prepared_data)
    return process_results(strategy, context, results)
end


# Initialize search context
context = MSSearchContext(
    MS_TABLE_PATHS,
    params_,
    spec_lib,
    ionMatches,
    ionMisses,
    IDtoCOL,
    ionTemplates,
    iso_splines,
    scored_PSMs,
    unscored_PSMs,
    spectral_scores
)

# Run searches sequentially
parameter_results = run_search(ParameterTuningSearch(), context)
nce_results = run_search(NCETuningSearch(), context)
quad_results = params_[:presearch_params]["estimate_quad_transmission"] ? 
    run_search(QuadTuningSearch(), context) : 
    create_default_quad_model()
first_search_results = run_search(FirstSearch(), context)
quant_results = run_search(QuantSearch(), context)
