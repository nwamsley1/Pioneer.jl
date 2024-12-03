# Define abstract types and traits
abstract type SearchMethod end
abstract type TuningMethod <: SearchMethod end
include("src/Routines/LibrarySearch/selectTransitions/selectTransitions.jl")
# Concrete strategies
struct ParameterTuningSearch <: TuningMethod end
struct NCETuningSearch <: TuningMethod end
struct QuadTuningSearch <: TuningMethod end
struct FirstSearch <: SearchMethod end
struct SecondSearch <: SearchMethod end
struct ChromatogramSearch <: SearchMethod end

# Common interface for containers of pre-allocated
# data structures used for conducting searches. Does not hold
# the search results, rather, intermediate date used during the search. 
abstract type SearchDataStructures end
getIonMatches(s::SearchDataStructures) = s.ion_matches
getIonMisses(s::SearchDataStructures) = s.ion_matches
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

#Common interface for references to mass spec data used in seraches
abstract type MassSpecDataReference end
getMSData(msdr::MassSpecDataReference, ms_file_idx::I) where {I<:Integer} = Arrow.Table(msdr.file_paths[ms_file_idx])
struct ArrowTableReference <: MassSpecDataReference
    file_paths::Vector{AbstractString}
end
import Base: enumerate
function enumerate(msdr::ArrowTableReference)
    return zip(1:length(msdr.file_paths), (getMSData(msdr, i) for i in 1:length(msdr.file_paths)))
end

abstract type SearchParameters end 

struct FragmentIndexSearchParameters{P<:PrecEstimation} <: SearchParameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_index_search_score::UInt8
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::UInt8
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    sample_rate::Float32
    irt_tol::Float32
    spec_order::Set{Int64}
    prec_estimation::P
end 
getIsotopeErrBounds(fsp::FragmentIndexSearchParameters) = fsp.isotope_err_bounds
getMinIndexSearchScore(fsp::FragmentIndexSearchParameters) = fsp.min_index_search_score
getMinFragCount(fsp::FragmentIndexSearchParameters) = fsp.min_frag_count
getMinSpectralContrast(fsp::FragmentIndexSearchParameters) = fsp.min_spectral_contrast
getMinLog2MatchedRatio(fsp::FragmentIndexSearchParameters) = fsp.min_log2_matched_ratio
getMinTopNofM(fsp::FragmentIndexSearchParameters) = fsp.min_topn_of_m
getMaxBestRank(fsp::FragmentIndexSearchParameters) = fsp.max_best_rank 
getNFragIsotopes(fsp::FragmentIndexSearchParameters) = fsp.n_frag_isotopes
getMaxFragRank(fsp::FragmentIndexSearchParameters) = fsp.max_frag_rank
getSampleRate(fsp::FragmentIndexSearchParameters) = fsp.sample_rate
getIRTTol(fsp::FragmentIndexSearchParameters) = fsp.irt_tol
getSpecOrder(fsp::FragmentIndexSearchParameters) = fsp.spec_order
getPrecEstimation(fsp::FragmentIndexSearchParameters) = fsp.prec_estimation

abstract type SpectralLibrary end
struct FragmentIndexLibrary <: SpectralLibrary
    presearch_fragment_index::FragmentIndex{Float32}
    fragment_index::FragmentIndex{Float32}
    precursors::Arrow.Table
    fragment_lookup_table::SplineFragmentLookup
end
getPresearchFragmentIndex(sl::SpectralLibrary) = sl.fragment_index
getFragmentIndex(sl::SpectralLibrary) = sl.fragment_index
getPrecursors(sl::SpectralLibrary) = sl.precursors
getFragmentLookupTable(sl::SpectralLibrary) = sl.fragment_lookup_table 

function getParameters(s::ParameterTuningSearch, params::Any)
    pp = params[:presearch_params]
    prec_estimation = pp["abreviate_precursor_calc"] ? FullPrecCapture() : PartialPrecCapture()
    return FragmentIndexSearchParameters(
        (zero(UInt8), zero(UInt8)),
        UInt8(pp["min_index_search_score"]),
        Int64(pp["min_frag_count"]),
        Float32(pp["min_spectral_contrast"]),
        Float32(pp["min_log2_matched_ratio"]),
        (
            Int64(first(pp["min_topn_of_m"])),
            Int64(last(pp["min_topn_of_m"]))
        ),
        UInt8(pp["max_best_rank"]),
        Int64(pp["n_frag_isotopes"]),
        UInt8(pp["max_frag_rank"]),
        Float32(pp["sample_rate"]),
        typemax(Float32),
        Set(2),
        prec_estimation
    )
end
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
    Tuple([initSimpleSearchContext(iso_splines, n_precursors, M) for _ in 1:N])
end

# Trait functions
#=
function prepare_search(
        search_type::ParameterTuningSearch, 
        context::Vector{SearchContext},
        params::Any) where {S<:SearchContext}
    search_parameters = getParameters(search_type, params)
    error("prepare_search not implemented for $(typeof(strategy))")
end
=#
function execute_search(
        search_type::ParameterTuningSearch, 
        msdr::MassSpecDataReference,
        spec_lib::SpectralLibrary,
        search_data::NTuple{N, S},
        params::Any) where {N, S<:SearchDataStructures}
        search_parameters = getParameters(search_type, params)
        mem = MassErrorModel(
                zero(Float32),
                #(getFragTolPpm(params), getFragTolPpm(params)),
                (30.0f0, 30.0f0)
        )
        for (ms_file_idx, spectra) in ProgressBar(enumerate(msdr))
            #println("ms_file_io")
            return LibrarySearch(
                spectra,
                UInt32(ms_file_idx),
                getFragmentIndex(spec_lib),
                spec_lib,
                search_data,
                search_parameters,
                GeneralGaussModel(5.0f0, 0.0f0),
                mem,
                x::Float32->x::Float32,
            )
        end
    #error("execute_search not implemented for $(typeof(strategy))")
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
