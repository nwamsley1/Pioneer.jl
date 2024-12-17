# Define abstract types and traits
#=
using Arrow, ArrowTypes, ArgParse
using LaTeXStrings
#using BSplineKit Don't need this imports anymore?
using Base64
using Base.Order
using Base.Iterators: partition
using CSV, CategoricalArrays, Combinatorics, CodecZlib
using DataFrames, DataStructures, Dictionaries, Distributions 
using FASTX
using Interpolations
using JSON, JLD2
using LinearAlgebra, LoopVectorization, LinearSolve, LightXML
using Measures
using NumericalIntegration
using Optim
using Plots, PrettyPrinting, Polynomials, PDFmerger, ProgressBars, Pkg
using Tables, Test
using StatsPlots
using Random
using StaticArrays, StatsBase, SpecialFunctions, Statistics
using XGBoost
using KernelDensity
using FastGaussQuadrature
include("structs/Ion.jl")
include("structs/LibraryFragmentIndex.jl")
include("structs/LibraryIon.jl")
include("structs/LibraryFragmentIndex.jl")
include("structs/LibraryIon.jl")
include("utils/ML/uniformBasisCubicSpline.jl")
 include("structs/RetentionTimeConversionModel.jl")
include("utils/quadTransmissionModeling/quadTransmissionModel.jl")
include("utils/isotopeSplines.jl")
include(joinpath(@__DIR__, "Routines","LibrarySearch","importScripts.jl"))
importScripts()
include("Routines/LibrarySearch/SearchMethods/SearchTypes.jl")
importScripts()
include("Routines/LibrarySearch/SearchMethods/NceTuningSearch.jl")
include("Routines/LibrarySearch/SearchMethods/SearchMethods.jl")
include("Routines/LibrarySearch/SearchMethods/ParameterTuningSearch.jl")
include("Routines/LibrarySearch/SearchMethods/FirstPassSearch.jl")
include("structs/Ion.jl")
include("structs/LibraryFragmentIndex.jl")
include("structs/LibraryIon.jl")
include("utils/ML/uniformBasisCubicSpline.jl")
 include("structs/RetentionTimeConversionModel.jl")
include("utils/quadTransmissionModeling/quadTransmissionModel.jl")
include("utils/isotopeSplines.jl")
include(joinpath(@__DIR__, "Routines","LibrarySearch","importScripts.jl"))
importScripts()
include("Routines/LibrarySearch/SearchMethods/SearchTypes.jl")
importScripts()
include("Routines/LibrarySearch/SearchMethods/SearchMethods.jl")
include("Routines/LibrarySearch/SearchMethods/ParameterTuningSearch.jl")
include("Routines/LibrarySearch/SearchMethods/FirstPassSearch.jl")

#include(joinpath(@__DIR__, "Routines","LibrarySearch","method"s,"loadSpectralLibrary.jl"))
const methods_path = joinpath(@__DIR__, "Routines","LibrarySearch")       
include(joinpath(@__DIR__, "Routines","SearchDIA.jl"))
include(joinpath(@__DIR__, "Routines","ThreeProteomeAnalysis.jl"))
const charge_facs = Float64[1, 0.9, 0.85, 0.8, 0.75]
const prediction_model_options =  Set(["unispec","prosit_2020_hcd","AlphaPeptDeep"])
const prediction_model_to_annotation_type = Dict(
    "unispec" => UniSpecFragAnnotation("y1^1"),
    "prosit_2020_hcd" => GenericFragAnnotation("y1+1"),
    "AlphaPeptDeep" => GenericFragAnnotation("y1+1")
)
const prediction_model_to_model_type = Dict(
    "unispec" => InstrumentSpecificModel("unispec"),
    "prosit_2020_hcd" => InstrumentAgnosticModel("prosit_2020_hcd"),
    "AlphaPeptDeep" => InstrumentSpecificModel("AlphaPeptDeep")
)
const H2O::Float64 = Float64(18.010565)
const PROTON::Float64 = Float64(1.0072764)
const NEUTRON::Float64 = Float64(1.00335)
const NCE_MODEL_BREAKPOINT::Float32 = Float32(500.0f0)

=#
#include("Routines/LibrarySearch/selectTransitions/selectTransitions.jl")
#include("src/Routine")
abstract type SearchMethod end
abstract type TuningMethod <: SearchMethod end

abstract type SearchResults end 
#include("src/Routines/LibrarySearch/selectTransitions/selectTransitions.jl")

# Common interface for containers of pre-allocated
# data structures used for conducting searches. Does not hold
# the search results, rather, intermediate data used during the search. 
abstract type SearchDataStructures end
getIonMatches(s::SearchDataStructures) = s.ion_matches
getIonMisses(s::SearchDataStructures) = s.ion_misses
getMassErrMatches(s::SearchDataStructures) = s.mass_err_matches
getIdToCol(s::SearchDataStructures) = s.id_to_col
getPrecursorScores(s::SearchDataStructures) = s.prec_count
getIonTemplates(s::SearchDataStructures) = s.ion_templates
getIsoSplines(s::SearchDataStructures) = s.iso_splines
getScoredPsms(s::SearchDataStructures) = s.scored_psms
getUnscoredPsms(s::SearchDataStructures) = s.unscored_psms
getSpectralScores(s::SearchDataStructures) = s.spectral_scores

mutable struct SimpleLibrarySearch{I<:IsotopeSplineModel} <: SearchDataStructures
    ion_matches::Vector{FragmentMatch{Float32}}
    ion_misses::Vector{FragmentMatch{Float32}}
    mass_err_matches::Vector{FragmentMatch{Float32}}
    id_to_col::ArrayDict{UInt32, UInt16}
    prec_count::Counter{UInt32, UInt8}
    ion_templates::Vector{DetailedFrag{Float32}}
    iso_splines::I
    scored_psms::Vector{SimpleScoredPSM{Float32, Float16}}
    unscored_psms::Vector{SimpleUnscoredPSM{Float32}}
    spectral_scores::Vector{SpectralScoresSimple{Float16}}
    #
    Hs::SparseArray           # For design matrix
    precursor_weights::Vector{Float32}
    temp_weights::Vector{Float32}        # For deconvolution weights
    residuals::Vector{Float32}      # For deconvolution residuals
    isotopes::Vector{Float32}       # For isotope calculations
    precursor_transmission::Vector{Float32}  # For transmission calculations
    tuning_results::Vector{@NamedTuple{precursor_idx::UInt32, scan_idx::UInt32, weight::Float32, iso_idx::UInt8, center_mz::Float32, n_matches::UInt8}}  # Store results per thread
end

getTempWeights(s::SimpleLibrarySearch) = s.temp_weights
getPrecursorWeights(s::SimpleLibrarySearch) = s.precursor_weights
#Common interface for references to mass spec data used in searches
abstract type MassSpecDataReference end
getMSData(msdr::MassSpecDataReference, ms_file_idx::I) where {I<:Integer} = Arrow.Table(msdr.file_paths[ms_file_idx])

struct ArrowTableReference <: MassSpecDataReference
    file_paths::Vector{AbstractString}
    file_id_to_name::Dict{Int64, String}
end

getParsedFileName(s::ArrowTableReference, ms_file_idx::Int64) = s.file_id_to_name[ms_file_idx]

import Base: enumerate
function enumerate(msdr::ArrowTableReference)
    return zip(1:length(msdr.file_paths), (getMSData(msdr, i) for i in 1:length(msdr.file_paths)))
end


# Common interface for containers for search data
mutable struct SearchContext{N,L<:FragmentIndexLibrary,M<:MassSpecDataReference}
    spec_lib::L
    temp_structures::AbstractVector{<:SearchDataStructures}
    mass_spec_data_reference::M
    data_out_dir::Base.Ref{String}
    qc_plot_folder::Base.Ref{String}
    rt_alignment_plot_folder::Base.Ref{String}
    mass_err_plot_folder::Base.Ref{String}
    quad_transmission_model::Dict{Int64, QuadTransmissionModel}
    mass_error_model::Dict{Int64, MassErrorModel}
    rt_to_irt_model::Dict{Int64, RtConversionModel}
    nce_model::Dict{Int64, NceModel}
    irt_rt_map::Base.Ref{Any}
    rt_irt_map::Base.Ref{Any}
    precursor_dict::Base.Ref{Dict}
    rt_index_paths::Base.Ref{Vector{String}}
    irt_errors::Dict{String, Float32}
    # New fields
    n_threads::Int64
    n_precursors::Int64
    buffer_size::Int64  # This is M in your notation

    function SearchContext(
        spec_lib::L,
        temp_structures::AbstractVector{<:SearchDataStructures},
        mass_spec_data_reference::M,
        n_threads::Int64,
        n_precursors::Int64,
        buffer_size::Int64
    ) where {L<:FragmentIndexLibrary,M<:MassSpecDataReference}
        N = length(temp_structures)
        new{N,L,M}(
            spec_lib,
            temp_structures,
            mass_spec_data_reference,
            Ref{String}(),
            Ref{String}(),
            Ref{String}(),
            Ref{String}(),
            Dict{Int64, QuadTransmissionModel}(),
            Dict{Int64, MassErrorModel}(),
            Dict{Int64, RtConversionModel}(),
            Dict{Int64, NceModel}(),
            Ref{Any}(),
            Ref{Any}(),
            Ref{Dict}(),
            Ref{Vector{String}}(),
            Dict{String, Float32}(),
            n_threads,
            n_precursors,
            buffer_size
        )
    end
end

# Add getters for new fields
getNThreads(s::SearchContext) = s.n_threads
getNPrecursors(s::SearchContext) = s.n_precursors
getBufferSize(s::SearchContext) = s.buffer_size
getMassSpecData(s::SearchContext) = s.mass_spec_data_reference
getSpecLib(s::SearchContext) = s.spec_lib
getSearchData(s::SearchContext) = s.temp_structures
getDataOutDir(s::SearchContext) = s.data_out_dir[]
getQcPlotfolder(s::SearchContext) = s.qc_plot_folder[]
getRtAlignPlotFolder(s::SearchContext) = s.rt_alignment_plot_folder[]
getMassErrPlotFolder(s::SearchContext) = s.mass_err_plot_folder[]
getParsedFileName(s::SearchContext, ms_file_idx::Int64) = getParsedFileName(s.mass_spec_data_reference, ms_file_idx)
getIrtRtMap(s::SearchContext) = s.irt_rt_map[]
getRtIrtMap(s::SearchContext) = s.rt_irt_map[]
getPrecursorDict(s::SearchContext) = s.precursor_dict[]
getRtIndexPaths(s::SearchContext) = s.rt_index_paths[]
getIrtErrors(s::SearchContext) = s.irt_errors

function getQuadTransmissionModel(s::SearchContext, index::I) where {I<:Integer}  
    if haskey(s.quad_transmission_model,index)
        return s.quad_transmission_model[index]
    else
        #Return a sensible default 
         @warn "Mass error model not found for ms_file_idx $index. Returning default GeneralGaussModel(5.0f0, 0.0f0)"
        return GeneralGaussModel(5.0f0, 0.0f0)
    end
end

function setQuadTransmissionModel!(s::SearchContext, index::I, model::QuadTransmissionModel)where {I<:Integer}
    s.quad_transmission_model[index] = model
end

function getMassErrorModel(s::SearchContext, index::I) where {I<:Integer} 
    if haskey(s.mass_error_model, index)
        return s.mass_error_model[index]
    else
        #Return a sensible default 
        @warn "Mass error model not found for ms_file_idx $index. Returning default +/- 30ppm"
        return MassErrorModel(zero(Float32), #(getFragTolPpm(params), getFragTolPpm(params)),
        (30.0f0, 30.0f0))
    end
end

function setMassErrorModel!(s::SearchContext, index::I, model::MassErrorModel) where {I<:Integer} 
    s.mass_error_model[index] = model
end

function getRtIrtModel(s::SearchContext, index::I) where {I<:Integer} 
    if haskey(s.rt_to_irt_model, index)
        return s.rt_to_irt_model[index]
    else
        #Return a sensible default 
        return IdentityModel()
    end
end

function setRtIrtModel!(s::SearchContext, index::I, model::Any) where {I<:Integer}
    s.rt_to_irt_model[index] = model
end

function getNceModelModel(s::SearchContext, index::I) where {I<:Integer} 
    if haskey(s.nce_model, index)
        return s.nce_model[index]
    else
        #Return a sensible default 
        return PiecewiseNceModel(30.0f0)
    end
end

function setNceModel!(s::SearchContext, index::I, model::NceModel) where {I<:Integer}
    s.nce_model[index] = model
end

function setDataOutDir!(s::SearchContext, dir::String)
    s.data_out_dir[] = dir
    qc_plot_folder = joinpath(dir, "qc_plots")
    !isdir(qc_plot_folder) && mkdir(qc_plot_folder)
    rt_alignment_folder = joinpath(qc_plot_folder, "rt_alignment_plots")
    mass_error_folder = joinpath(qc_plot_folder, "mass_error_plots")
    !isdir(rt_alignment_folder) && mkdir(rt_alignment_folder)
    !isdir(mass_error_folder) && mkdir(mass_error_folder)
    s.qc_plot_folder[] = qc_plot_folder
    s.rt_alignment_plot_folder[] = rt_alignment_folder
    s.mass_err_plot_folder[] = mass_error_folder
end

function setIrtRtMap!(s::SearchContext, map::Any)
    s.irt_rt_map[] = map
end

function setRtIrtMap!(s::SearchContext, map::Any)
    s.rt_irt_map[] = map
end

function setPrecursorDict!(s::SearchContext, dict::Dict)
    s.precursor_dict[] = dict
end

function setRtIndexPaths!(s::SearchContext, paths::Vector{String})
    s.rt_index_paths[] = paths
end

#function setIrtErrors!(s::SearchContext, errors::Dict{String, Float32})
#    s.irt_errors = errors
#end

# Add getter/setter for irt_errs
getIrtErrs(s::SearchContext) = s.irt_errors
function setIrtErrs!(s::SearchContext, errs::Dict{String, Float32})
    for (k,v) in pairs(errs)
        s.irt_errors[k] = v
    end
end

#Common interface for tunable search parameters
abstract type SearchParameters end 
abstract type FragmentIndexSearchParameters <: SearchParameters end
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
getMaxPresearchIters(fsp::SearchParameters) = fsp.max_presearch_iters
getIRTTol(fsp::SearchParameters) = fsp.irt_tol
getSpecOrder(fsp::SearchParameters) = fsp.spec_order
getPrecEstimation(fsp::SearchParameters) = fsp.prec_estimation
getFragTolPpm(fsp::SearchParameters) = fsp.frag_tol_ppm
getSplineDegree(fsp::SearchParameters) = fsp.spline_degree
getSplineNKnots(fsp::SearchParameters) = fsp.spline_n_knots
getOutlierThreshold(fsp::SearchParameters) = fsp.spline_fit_outlier_sd
getMinPsms(fsp::SearchParameters) = fsp.min_psms
getMaxQVal(fsp::SearchParameters) = fsp.max_q_val
getFragErrQuantile(fsp::SearchParameters) = fsp.frag_err_quantile


function initSimpleSearchContext(
    iso_splines::IsotopeSplineModel,
    n_precursors::Int64,
    M::Int64
    )
    SimpleLibrarySearch(
        [FragmentMatch{Float32}() for _ in range(1, M)],
        [FragmentMatch{Float32}() for _ in range(1, M)],
        [FragmentMatch{Float32}() for _ in range(1, M)],
        ArrayDict(UInt32, UInt16, n_precursors),
        Counter(UInt32, UInt8,n_precursors ),
        [DetailedFrag{Float32}() for _ in range(1, M)],
        iso_splines,
        Vector{SimpleScoredPSM{Float32, Float16}}(undef, 5000),
        [SimpleUnscoredPSM{Float32}() for _ in range(1, 5000)],
        Vector{SpectralScoresSimple{Float16}}(undef, 5000),
        SparseArray(UInt32(5000)),
        zeros(Float32, n_precursors),
        zeros(Float32, 5000),
        zeros(Float32, 5000),
        zeros(Float32, 5),
        zeros(Float32, 5),
        Vector{@NamedTuple{precursor_idx::UInt32, scan_idx::UInt32, weight::Float32, iso_idx::UInt8, center_mz::Float32, n_matches::UInt8}}()
    )
end



getHs(s::SearchDataStructures) = s.Hs
getWeights(s::SearchDataStructures) = s.weights
getResiduals(s::SearchDataStructures) = s.residuals
getIsotopes(s::SearchDataStructures) = s.isotopes
getPrecursorTransmission(s::SearchDataStructures) = s.precursor_transmission
getTuningResults(s::SearchDataStructures) = s.tuning_results

function initSimpleSearchContexts(
    iso_splines::IsotopeSplineModel,
    n_precursors::Int64,
    N::Int64,
    M::Int64)
    [initSimpleSearchContext(iso_splines, n_precursors, M) for _ in 1:N]
end

# Update conversion back to simple search
function convertToSimpleContext!(search_context::SearchContext)
    search_context.temp_structures = initSimpleSearchContexts(
        getIsoSplines(search_context.temp_structures[1]),
        getNPrecursors(search_context),
        getNThreads(search_context),
        getBufferSize(search_context)
    )
    return search_context
end


# Update the initialization of SearchContext
function initSearchContext(
    spec_lib::FragmentIndexLibrary,
    iso_splines::IsotopeSplineModel,
    ms_data_reference::MassSpecDataReference,
    n_threads::Int64,
    n_precursors::Int64,
    buffer_size::Int64
)
    temp_structures = initSimpleSearchContexts(
        iso_splines,
        n_precursors,
        n_threads,
        buffer_size
    )
    
    return SearchContext(
        spec_lib,
        temp_structures,
        ms_data_reference,
        n_threads,
        n_precursors,
        buffer_size
    )
end