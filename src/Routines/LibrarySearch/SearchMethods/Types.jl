# Define abstract types and traits
#=
include("structs/Ion.jl")
include("structs/LibraryFragmentIndex.jl")
include("structs/LibraryIon.jl")
include("utils/ML/uniformBasisCubicSpline.jl")
 include("structs/RetentionTimeConversionModel.jl")
include("utils/quadTransmissionModeling/quadTransmissionModel.jl")
include("utils/isotopeSplines.jl")
include(joinpath(@__DIR__, "Routines","LibrarySearch","importScripts.jl"))
importScripts()
include("Routines/LibrarySearch/SearchMethods/Types.jl")
importScripts()
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

struct SimpleLibrarySearch{I<:IsotopeSplineModel} <: SearchDataStructures
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
end

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
struct SearchContext{N,L<:FragmentIndexLibrary,S<:SearchDataStructures,M<:MassSpecDataReference}
    spec_lib::L
    temp_structures::AbstractVector{S}
    mass_spec_data_reference::M
    data_out_dir::Base.Ref{String}
    qc_plot_folder::Base.Ref{String}
    rt_alignment_plot_folder::Base.Ref{String}
    mass_err_plot_folder::Base.Ref{String}
    quad_transmission_model::Dict{Int64, QuadTransmissionModel}
    mass_error_model::Dict{Int64, MassErrorModel}
    rt_to_irt_model::Dict{Int64, RtConversionModel}
    nce_model::Dict{Int64, NceModel}
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
            Ref{String}(),
            Ref{String}(),
            Ref{String}(),
            Ref{String}(),
            Dict{Int64, QuadTransmissionModel}(),
            Dict{Int64, MassErrorModel}(),
            Dict{Int64, RtConversionModel}(),
            Dict{Int64, NceModel}()
        )
    end

end

getMassSpecData(s::SearchContext) = s.mass_spec_data_reference
getSpecLib(s::SearchContext) = s.spec_lib
getSearchData(s::SearchContext) = s.temp_structures
getDataOutDir(s::SearchContext) = s.data_out_dir[]
getQcPlotfolder(s::SearchContext) = s.qc_plot_folder[]
getRtAlignPlotFolder(s::SearchContext) = s.rt_alignment_plot_folder[]
getMassErrPlotFolder(s::SearchContext) = s.mass_err_plot_folder[]
getParsedFileName(s::SearchContext, ms_file_idx::Int64) = getParsedFileName(s.mass_spec_data_reference, ms_file_idx)

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
        return IdentityModel()
    end
end

function setRtIrtModel!(s::SearchContext, index::Int64, model::Any)
    s.rt_to_irt_model[index] = model
end

function getNceModelModel(s::SearchContext, index::Int64) 
    if haskey(s.nce_model, index)
        return s.nce_model[index]
    else
        #Return a sensible default 
        return PiecewiseNceModel(30.0f0)
    end
end

function setNceModel!(s::SearchContext, index::Int64, model::NceModel)
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