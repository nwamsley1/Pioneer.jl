"""
Core types and interfaces for implementing search methods in mass spectrometry analysis.

# Search Method Implementation Guide

To implement a new search method:

1. Create types:
   - Subtype SearchMethod or TuningMethod 
   - Subtype SearchResults to store search outputs
   - Subtype SearchParameters to define search configuration
   - (Optional) Subtype SearchDataStructures if custom intermediate data needed

2. Implement required interface methods:
   - get_parameters(::YourSearchMethod, params::Any)
   - init_search_results(::YourParameters, ::SearchContext, ::Int64)
   - process_file!(::YourResults, ::YourParameters, ::SearchContext, ::Int64, ::Arrow.Table)
   - process_search_results!(::YourResults, ::YourParameters, ::SearchContext, ::Int64)
   - summarize_results!(::YourResults, ::YourParameters, ::SearchContext)

See ParameterTuningSearch or FirstPassSearch for example implementations.
"""

#==========================================================
Abstract Types and Interfaces
==========================================================#

"""
Base type for all search methods. Subtypes should implement the core search interface.
"""
abstract type SearchMethod end

"""
Specialized search method for parameter tuning operations.
"""
abstract type TuningMethod <: SearchMethod end

"""
Base type for search results. Subtypes define the structure of search outputs.
"""
abstract type SearchResults end

"""
Base type for search data structures used during execution.
Holds intermediate data and working arrays, not final results.
"""
abstract type SearchDataStructures end

"""
Base type for mass spec data references. Defines how to access MS data files.
"""
abstract type MassSpecDataReference end

"""
Base type for search parameters. Subtypes define configuration for specific searches.
"""
abstract type SearchParameters end

"""
Parameters specifically for fragment index-based searches.
"""
abstract type FragmentIndexSearchParameters <: SearchParameters end

#==========================================================
Concrete Types
==========================================================#

"""
Reference to MS data stored in Arrow files.
"""
struct ArrowTableReference <: MassSpecDataReference
    file_paths::Vector{AbstractString}
    file_id_to_name::Dict{Int64, String}
end

"""
Basic search data structure for library searches.
Contains pre-allocated arrays and intermediate data structures.
"""
mutable struct SimpleLibrarySearch{I<:IsotopeSplineModel} <: SearchDataStructures
    # Match data
    ion_matches::Vector{FragmentMatch{Float32}}
    ion_misses::Vector{FragmentMatch{Float32}}
    mass_err_matches::Vector{FragmentMatch{Float32}}
    
    # Indexing and scoring
    id_to_col::ArrayDict{UInt32, UInt16}
    prec_count::Counter{UInt32, UInt8}
    ion_templates::Vector{DetailedFrag{Float32}}
    iso_splines::I
    
    # PSM scoring
    scored_psms::Vector{SimpleScoredPSM{Float32, Float16}}
    unscored_psms::Vector{SimpleUnscoredPSM{Float32}}
    spectral_scores::Vector{SpectralScoresSimple{Float16}}
    
    # Working arrays
    Hs::SparseArray
    prec_ids::Vector{UInt32}
    precursor_weights::Vector{Float32}
    temp_weights::Vector{Float32}
    residuals::Vector{Float32}
    isotopes::Vector{Float32}
    precursor_transmission::Vector{Float32}
    tuning_results::Vector{@NamedTuple{
        precursor_idx::UInt32,
        scan_idx::UInt32,
        weight::Float32,
        iso_idx::UInt8,
        center_mz::Float32,
        n_matches::UInt8
    }}
end

"""
Primary search context holding all data structures and state for search execution.
"""
mutable struct SearchContext{N,L<:FragmentIndexLibrary,M<:MassSpecDataReference}
    # Core components
    spec_lib::L
    temp_structures::AbstractVector{<:SearchDataStructures}
    mass_spec_data_reference::M
    
    # Output directories
    data_out_dir::Base.Ref{String}
    qc_plot_folder::Base.Ref{String}
    rt_alignment_plot_folder::Base.Ref{String}
    mass_err_plot_folder::Base.Ref{String}
    
    # Models and mappings
    quad_transmission_model::Dict{Int64, QuadTransmissionModel}
    mass_error_model::Dict{Int64, MassErrorModel}
    rt_to_irt_model::Dict{Int64, RtConversionModel}
    nce_model::Dict{Int64, NceModel}
    huber_delta::Base.Ref{Float32}

    # Results and paths
    irt_rt_map::Dict{Int64, RtConversionModel}
    rt_irt_map::Dict{Int64, RtConversionModel}
    precursor_dict::Base.Ref{Dictionary}
    rt_index_paths::Base.Ref{Vector{String}}
    irt_errors::Dict{Int64, Float32}
    
    # Configuration
    n_threads::Int64
    n_precursors::Int64
    buffer_size::Int64

    # Constructor
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
            spec_lib, temp_structures, mass_spec_data_reference,
            Ref{String}(), Ref{String}(), Ref{String}(), Ref{String}(),
            Dict{Int64, QuadTransmissionModel}(),
            Dict{Int64, MassErrorModel}(),
            Dict{Int64, RtConversionModel}(),
            Dict{Int64, NceModel}(), Ref(100000.0f0),
            Dict{Int64, RtConversionModel}(), 
            Dict{Int64, RtConversionModel}(), 
            Ref{Dictionary}(), 
            Ref{Vector{String}}(),
            Dict{Int64, Float32}(),
            n_threads, n_precursors, buffer_size
        )
    end
end

#==========================================================
Interface Methods for Parameter Access
==========================================================#
#MassSpecDataReference interface getters 
getMSData(msdr::MassSpecDataReference, ms_file_idx::I) where {I<:Integer} = Arrow.Table(msdr.file_paths[ms_file_idx])
getParsedFileName(s::ArrowTableReference, ms_file_idx::Int64) = s.file_id_to_name[ms_file_idx]

import Base: enumerate
function enumerate(msdr::ArrowTableReference)
    return zip(1:length(msdr.file_paths), (getMSData(msdr, i) for i in 1:length(msdr.file_paths)))
end

# SearchParameters interface getters
getFragErrQuantile(fsp::SearchParameters)      = fsp.frag_err_quantile
getFragTolPpm(fsp::SearchParameters)           = fsp.frag_tol_ppm
getIRTTol(fsp::SearchParameters)               = fsp.irt_tol
getIsotopeErrBounds(fsp::SearchParameters)      = fsp.isotope_err_bounds
getMaxBestRank(fsp::SearchParameters)           = fsp.max_best_rank
getMaxFragRank(fsp::SearchParameters)           = fsp.max_frag_rank
getMaxPresearchIters(fsp::SearchParameters)    = fsp.max_presearch_iters
getMaxQVal(fsp::SearchParameters)              = fsp.max_q_val
getMinFragCount(fsp::SearchParameters)          = fsp.min_frag_count
getMinIndexSearchScore(fsp::SearchParameters)   = fsp.min_index_search_score
getMinLog2MatchedRatio(fsp::SearchParameters)   = fsp.min_log2_matched_ratio
getMinPsms(fsp::SearchParameters)              = fsp.min_psms
getMinSpectralContrast(fsp::SearchParameters)   = fsp.min_spectral_contrast
getMinTopNofM(fsp::SearchParameters)            = fsp.min_topn_of_m
getNFragIsotopes(fsp::SearchParameters)         = fsp.n_frag_isotopes
getOutlierThreshold(fsp::SearchParameters)     = fsp.spline_fit_outlier_sd
getPrecEstimation(fsp::SearchParameters)       = fsp.prec_estimation
getSampleRate(fsp::SearchParameters)           = fsp.sample_rate
getSpecOrder(fsp::SearchParameters)            = fsp.spec_order
getSplineDegree(fsp::SearchParameters)         = fsp.spline_degree
getSplineNKnots(fsp::SearchParameters)         = fsp.spline_n_knots
# ... rest of parameter getters ...

# SearchDataStructures interface getters

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
getHs(s::SearchDataStructures) = s.Hs
getPrecIds(s::SearchDataStructures) = s.prec_ids
getWeights(s::SearchDataStructures) = s.weights
getResiduals(s::SearchDataStructures) = s.residuals
getIsotopes(s::SearchDataStructures) = s.isotopes
getPrecursorTransmission(s::SearchDataStructures) = s.precursor_transmission
getTuningResults(s::SearchDataStructures) = s.tuning_results
getTempWeights(s::SimpleLibrarySearch) = s.temp_weights
getPrecursorWeights(s::SimpleLibrarySearch) = s.precursor_weights

#==========================================================
SearchContext Getters and Setters
==========================================================#


# Simple getters
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
getIrtRtMap(s::SearchContext) = s.irt_rt_map
getRtIrtMap(s::SearchContext) = s.rt_irt_map
getPrecursorDict(s::SearchContext) = s.precursor_dict[]
getRtIndexPaths(s::SearchContext) = s.rt_index_paths[]
getIrtErrors(s::SearchContext) = s.irt_errors
getHuberDelta(s::SearchContext) = s.huber_delta[]

"""
   getQuadTransmissionModel(s::SearchContext, index::Integer)

Get quad transmission model for MS file index. Returns default model if not found.
"""
function getQuadTransmissionModel(s::SearchContext, index::I) where {I<:Integer}  
   if haskey(s.quad_transmission_model,index)
       return s.quad_transmission_model[index]
   else
       @warn "Quad Transmission model not found for ms_file_idx $index. Returning default GeneralGaussModel(5.0f0, 0.0f0)"
       return GeneralGaussModel(5.0f0, 0.0f0)
   end
end

# Simple setters
setQuadTransmissionModel!(s::SearchContext, index::I, model::QuadTransmissionModel) where {I<:Integer} = (s.quad_transmission_model[index] = model)

"""
   getMassErrorModel(s::SearchContext, index::Integer)

Get mass error model for MS file index. Returns default ±30ppm model if not found.
"""
function getMassErrorModel(s::SearchContext, index::I) where {I<:Integer} 
   if haskey(s.mass_error_model, index)
       return s.mass_error_model[index]
   else
       @warn "Mass error model not found for ms_file_idx $index. Returning default +/- 30ppm"
       return MassErrorModel(zero(Float32), (30.0f0, 30.0f0))
   end
end

setMassErrorModel!(s::SearchContext, index::I, model::MassErrorModel) where {I<:Integer} = (s.mass_error_model[index] = model)

"""
   getRtIrtModel(s::SearchContext, index::Integer)

Get RT to iRT conversion model for MS file index. Returns identity model if not found.
"""
function getRtIrtModel(s::SearchContext, index::I) where {I<:Integer} 
   if haskey(s.rt_to_irt_model, index)
       return s.rt_to_irt_model[index]
   else
       return IdentityModel()
   end
end

setRtIrtModel!(s::SearchContext, index::I, model::Any) where {I<:Integer} = (s.rt_to_irt_model[index] = model)

"""
   getNceModelModel(s::SearchContext, index::Integer)

Get NCE model for MS file index. Returns default 30 NCE model if not found.
"""
function getNceModelModel(s::SearchContext, index::I) where {I<:Integer} 
   if haskey(s.nce_model, index)
       return s.nce_model[index]
   else
       return PiecewiseNceModel(30.0f0)
   end
end

# More simple setters
setNceModel!(s::SearchContext, index::I, model::NceModel) where {I<:Integer} = (s.nce_model[index] = model)
function setIrtRtMap!(s::SearchContext, rcm::RtConversionModel, index::I) where {I<:Integer }
    s.irt_rt_map[index] = rcm
end
function setRtIrtMap!(s::SearchContext, rcm::RtConversionModel, index::I) where {I<:Integer}
    s.rt_irt_map[index] = rcm
end
function setPrecursorDict!(s::SearchContext, dict::Dictionary{UInt32, @NamedTuple{best_prob::Float32, best_ms_file_idx::UInt32, best_scan_idx::UInt32, best_irt::Float32, mean_irt::Union{Missing, Float32}, var_irt::Union{Missing, Float32}, n::Union{Missing, UInt16}, mz::Float32}})
    s.precursor_dict[] = dict
end
setRtIndexPaths!(s::SearchContext, paths::Vector{String}) = (s.rt_index_paths[] = paths)
setHuberDelta!(s::SearchContext, delta::Float32) = (s.huber_delta[] = delta)
function setIrtErrors!(s::SearchContext, errs::Dictionary{Int64, Float32})
    for (k,v) in pairs(errs)
        s.irt_errors[k] = v
    end
end

# ... rest of data structure getters ...