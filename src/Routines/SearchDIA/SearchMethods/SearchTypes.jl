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

"""
Base type for chromatogram extraction methods.
Used to distinguish between MS1 and MS2 chromatogram searches.
"""
abstract type CHROMATOGRAM end

"""
MS2 chromatogram type for fragment-based searches.
"""
struct MS2CHROM <: CHROMATOGRAM end

"""
MS1 chromatogram type for precursor-based searches.
"""
struct MS1CHROM <: CHROMATOGRAM end

#==========================================================
Concrete Types
==========================================================#

"""
Reference to MS data stored in Arrow files.
"""
struct ArrowTableReference{N} <: MassSpecDataReference
    file_paths::NTuple{N, String}
    file_id_to_name::NTuple{N, String}
    first_pass_psms::Vector{String}
    second_pass_psms::Vector{String}
    passing_psms::Vector{String}
    passing_proteins::Vector{String}
    rt_index_paths::Vector{String}
    failed_search_indicator::Vector{Bool}

    # Internal constructor
    function ArrowTableReference(file_paths::Vector{String})
        file_paths = [arrow_path for arrow_path in file_paths if endswith(arrow_path, ".arrow")]
        file_id_to_name = parseFileNames(file_paths)
        if length(file_id_to_name) != length(file_paths)
            file_id_to_name = ["" for x in 1:length(file_id_to_name)]
            @user_warn "Improper File Names Parsing. "
        elseif length(file_paths) == 0
            @user_warn "Could not find any files ending in `arrow` in the paths supplied: $file_paths"
        end
        n = length(file_paths)
        new{n}(
            NTuple{n, String}(file_paths), 
            NTuple{n, String}(file_id_to_name),
            fill("", n),
            fill("", n),
            fill("", n),
            fill("", n),
            fill("", n),
            fill(false, n)
        )
    end

    # Internal constructor
    function ArrowTableReference(file_dir::String)
        file_paths = [arrow_path for arrow_path in readdir(file_dir, join=true) if endswith(arrow_path, ".arrow")]
        if length(file_paths) == 0
            @user_warn "Could not find any files ending in `arrow` in the directory: $file_dir"
        end
        n = length(file_paths)
        new{n}(
            NTuple{n, String}(file_paths...), 
            fill("", n),
            fill("", n),
            fill("", n),
            fill("", n),
            fill("", n),
            fill("", n),
            fill(false, n)
        )
    end

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
    complex_scored_psms::Vector{ComplexScoredPSM{Float32, Float16}}
    complex_unscored_psms::Vector{ComplexUnscoredPSM{Float32}}
    complex_spectral_scores::Vector{SpectralScoresComplex{Float16}}
    ms1_scored_psms::Vector{Ms1ScoredPSM{Float32, Float16}}
    ms1_unscored_psms::Vector{Ms1UnscoredPSM{Float32}}
    ms1_spectral_scores::Vector{SpectralScoresMs1{Float16}}

    # Working arrays
    Hs::SparseArray
    prec_ids::Vector{UInt32}
    precursor_weights::Vector{Float32}
    temp_weights::Vector{Float32}
    residuals::Vector{Float32}
    isotopes::Vector{Float32}
    precursor_transmission::Vector{Float32}
end

"""
Primary search context holding all data structures and state for search execution.
"""
mutable struct SearchContext{N,L<:SpectralLibrary,M<:MassSpecDataReference}
    # Core components
    spec_lib::L
    temp_structures::AbstractVector{<:SearchDataStructures}
    mass_spec_data_reference::M
    
    # Output directories
    data_out_dir::Base.Ref{String}
    qc_plot_folder::Base.Ref{String}
    rt_alignment_plot_folder::Base.Ref{String}
    mass_err_plot_folder::Base.Ref{String}
    ms1_mass_err_plot_folder::Base.Ref{String}
    
    # Models and mappings
    quad_transmission_model::Dict{Int64, QuadTransmissionModel}
    mass_error_model::Dict{Int64, MassErrorModel}
    ms1_mass_error_model::Dict{Int64, MassErrorModel}
    #rt_to_irt_model::Dict{Int64, RtConversionModel}
    nce_model::Dict{Int64, NceModel}
    huber_delta::Base.Ref{Float32}
    deconvolution_stop_tolerance::Base.Ref{Float32}
    # Results and paths
    irt_rt_map::Dict{Int64, RtConversionModel}
    rt_irt_map::Dict{Int64, RtConversionModel}
    rt_to_refined_irt_map::Dict{Int64, RtConversionModel}
    refined_irt_to_rt_map::Dict{Int64, RtConversionModel}
    irt_refinement_models::Dict{Int64, Union{IrtRefinementModel, Nothing}}
    precursor_dict::Base.Ref{Dictionary}
    rt_index_paths::Base.Ref{Vector{String}}
    irt_errors::Dict{Int64, Float32}
    irt_obs::Dict{UInt32, Float32}
    pg_score_to_qval::Ref{Any}
    pg_name_to_global_pg_score::Ref{Dict{ProteinKey, Float32}}
    global_pg_score_to_qval_dict::Ref{Dict{Tuple{String,Bool,UInt8}, Float32}}
    pg_score_to_pep::Ref{Any}
    
    # Method results storage
    method_results::Dict{Type{<:SearchMethod}, Any}
    
    # Configuration
    n_threads::Int64
    n_precursors::Int64
    buffer_size::Int64
    
    # Library target/decoy statistics for FDR calculation
    n_library_targets::Int64
    n_library_decoys::Int64
    library_fdr_scale_factor::Float32

    # Failed file tracking
    failed_files::Set{Int64}
    file_failure_reasons::Dict{Int64, String}

    # Constructor
    function SearchContext(
        spec_lib::L,
        temp_structures::AbstractVector{<:SearchDataStructures},
        mass_spec_data_reference::M,
        n_threads::Int64,
        n_precursors::Int64,
        buffer_size::Int64
    ) where {L<:SpectralLibrary,M<:MassSpecDataReference}
        N = length(temp_structures)
        new{N,L,M}(
            spec_lib, temp_structures, mass_spec_data_reference,
            Ref{String}(), Ref{String}(), Ref{String}(), Ref{String}(),Ref{String}(),
            Dict{Int64, QuadTransmissionModel}(),
            Dict{Int64, MassErrorModel}(),
            Dict{Int64, MassErrorModel}(),
            Dict{Int64, NceModel}(), Ref(100000.0f0), 10.0f0,
            Dict{Int64, RtConversionModel}(),
            Dict{Int64, RtConversionModel}(),
            Dict{Int64, RtConversionModel}(),
            Dict{Int64, RtConversionModel}(),
            Dict{Int64, Union{IrtRefinementModel, Nothing}}(),
            Ref{Dictionary}(), 
            Ref{Vector{String}}(),
            Dict{Int64, Float32}(),
            Dict{UInt32, Float32}(),
            Ref{Any}(), Ref(Dict{ProteinKey, Float32}()), Ref(Dict{Tuple{String,Bool,UInt8}, Float32}()), Ref{Any}(),
            Dict{Type{<:SearchMethod}, Any}(),  # Initialize method_results
            n_threads, n_precursors, buffer_size,
            0, 0, 1.0f0,  # Initialize library stats with defaults
            Set{Int64}(),  # Initialize failed_files
            Dict{Int64, String}()  # Initialize file_failure_reasons
        )
    end
end

#==========================================================
Failed File Tracking Functions
==========================================================#

"""
    markFileFailed!(ctx::SearchContext, ms_file_idx::Int64, reason::String)

Mark a file as failed with the given reason.
"""
function markFileFailed!(ctx::SearchContext, ms_file_idx::Int64, reason::String)
    push!(ctx.failed_files, ms_file_idx)
    ctx.file_failure_reasons[ms_file_idx] = reason
end



#==========================================================
Interface Methods for Parameter Access
==========================================================#
#MassSpecDataReference interface getters 
getMSData(msdr::MassSpecDataReference, ms_file_idx::I) where {I<:Integer} = BasicMassSpecData(msdr.file_paths[ms_file_idx])
getMSData(sc::SearchContext) = sc.mass_spec_data_reference
getParsedFileName(s::ArrowTableReference, ms_file_idx::Int64) = s.file_id_to_name[ms_file_idx]

# Add length method for ArrowTableReference
Base.length(::ArrowTableReference{N}) where N = N

import Base: enumerate
function enumerate(msdr::ArrowTableReference)
    return zip(1:length(msdr.file_paths), (getMSData(msdr, i) for i in 1:length(msdr.file_paths)))
end

# Getter methods
getFileIdToName(ref::ArrowTableReference, index::Int) = ref.file_id_to_name[index]
getFirstPassPsms(ref::ArrowTableReference, index::Int) = ref.first_pass_psms[index]
getSecondPassPsms(ref::ArrowTableReference, index::Int) = ref.second_pass_psms[index]
getPassingPsms(ref::ArrowTableReference, index::Int) = ref.passing_psms[index]
getPassingProteins(ref::ArrowTableReference, index::Int) = ref.passing_proteins[index]
getRtIndex(ref::ArrowTableReference, index::Int) = ref.rt_index_paths[index]
getFailedIndicator(ref::ArrowTableReference, index::Int) = ref.failed_search_indicator[index]
getParsedFileNames(ref::ArrowTableReference) = ref.file_id_to_name

getFilePaths(ref::ArrowTableReference) = ref.file_paths

getFileIdToName(ref::ArrowTableReference) = ref.file_id_to_name
getFirstPassPsms(ref::ArrowTableReference) = ref.first_pass_psms
getSecondPassPsms(ref::ArrowTableReference) = ref.second_pass_psms
getPassingPsms(ref::ArrowTableReference) = ref.passing_psms
getPassingProteins(ref::ArrowTableReference) = ref.passing_proteins
getRtIndex(ref::ArrowTableReference) = ref.rt_index_paths


# Setter methods
setFileIdToName!(ref::ArrowTableReference, index::Int, value::String) = ref.file_id_to_name[index] = value
setFirstPassPsms!(ref::ArrowTableReference, index::Int, value::String) = ref.first_pass_psms[index] = value
setSecondPassPsms!(ref::ArrowTableReference, index::Int, value::String) = ref.second_pass_psms[index] = value
setPassingPsms!(ref::ArrowTableReference, index::Int, value::String) = ref.passing_psms[index] = value
setPassingProteins!(ref::ArrowTableReference, index::Int, value::String) = ref.passing_proteins[index] = value
setRtIndex!(ref::ArrowTableReference, index::Int, value::String) = ref.rt_index_paths[index] = value
setFailedIndicator!(ref::ArrowTableReference, index::Int, value::Bool) = ref.failed_search_indicator[index] = value

# SearchParameters interface getters
getFragErrQuantile(fsp::SearchParameters)      = fsp.frag_err_quantile
getFragTolPpm(fsp::SearchParameters)           = fsp.frag_tol_ppm
getIRTTol(fsp::SearchParameters)               = fsp.irt_tol
getIsotopeErrBounds(fsp::SearchParameters)      = fsp.isotope_err_bounds
getMinFractionTransmitted(fsp::SearchParameters) = fsp.min_fraction_transmitted
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

getComplexScoredPsms(s::SearchDataStructures) = s.complex_scored_psms
getComplexUnscoredPsms(s::SearchDataStructures) = s.complex_unscored_psms
getComplexSpectralScores(s::SearchDataStructures) = s.complex_spectral_scores

getMs1ScoredPsms(s::SearchDataStructures) = s.ms1_scored_psms
getMs1UnscoredPsms(s::SearchDataStructures) = s.ms1_unscored_psms
getMs1SpectralScores(s::SearchDataStructures) = s.ms1_spectral_scores


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
getMs1MassErrPlotFolder(s::SearchContext) = s.ms1_mass_err_plot_folder[]
getParsedFileName(s::SearchContext, ms_file_idx::Int64) = getParsedFileName(s.mass_spec_data_reference, ms_file_idx)
getIrtRtMap(s::SearchContext) = s.irt_rt_map
getRtIrtMap(s::SearchContext) = s.rt_irt_map
getPrecursorDict(s::SearchContext) = s.precursor_dict[]
getRtIndexPaths(s::SearchContext) = s.rt_index_paths[]
getIrtErrors(s::SearchContext) = s.irt_errors
getPredIrt(s::SearchContext) = s.irt_obs
getPredIrt(s::SearchContext, prec_idx::Int64) = s.irt_obs[prec_idx]
getPredIrt(s::SearchContext, prec_idx::UInt32) = s.irt_obs[prec_idx]
getHuberDelta(s::SearchContext) = s.huber_delta[]
setPredIrt!(s::SearchContext, prec_idx::Int64, irt::Float32) = s.irt_obs[prec_idx] = irt
setPredIrt!(s::SearchContext, prec_idx::UInt32, irt::Float32) = s.irt_obs[prec_idx] = irt
getLibraryTargetCount(s::SearchContext) = s.n_library_targets
getLibraryDecoyCount(s::SearchContext) = s.n_library_decoys
getLibraryFdrScaleFactor(s::SearchContext) = s.library_fdr_scale_factor
"""
   getQuadTransmissionModel(s::SearchContext, index::Integer)

Get quad transmission model for MS file index. Returns default model if not found.
"""
function getQuadTransmissionModel(s::SearchContext, index::I) where {I<:Integer}  
   if haskey(s.quad_transmission_model,index)
       return s.quad_transmission_model[index]
   else
       @user_warn "Quad Transmission model not found for ms_file_idx $index. Returning default GeneralGaussModel(5.0f0, 0.0f0)"
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
       @user_warn "Mass error model not found for ms_file_idx $index. Returning default +/- 30ppm"
       return MassErrorModel(zero(Float32), (30.0f0, 30.0f0))
   end
end
setMassErrorModel!(s::SearchContext, index::I, model::MassErrorModel) where {I<:Integer} = (s.mass_error_model[index] = model)

"""
   getMassErrorModel(s::SearchContext, index::Integer)

Get mass error model for MS file index. Returns default ±30ppm model if not found.
"""
function getMs1MassErrorModel(s::SearchContext, index::I) where {I<:Integer} 
   if haskey(s.ms1_mass_error_model, index)
       return s.ms1_mass_error_model[index]
   else
       @user_warn "Mass error model not found for ms_file_idx $index. Returning default +/- 30ppm"
       return MassErrorModel(zero(Float32), (30.0f0, 30.0f0))
   end
end


setMs1MassErrorModel!(s::SearchContext, index::I, model::MassErrorModel) where {I<:Integer} = (s.ms1_mass_error_model[index] = model)

"""
   getRtIrtModel(s::SearchContext, index::Integer)

Get RT to iRT conversion model for MS file index. Returns identity model if not found.
"""
function getRtIrtModel(s::SearchContext, index::I) where {I<:Integer} 
   if haskey(s.rt_irt_map, index)
       return s.rt_irt_map[index]
   else
       return IdentityModel()
   end
end

"""
   getRtIrtModel(s::SearchContext, index::Integer)

Get RT to iRT conversion model for MS file index. Returns identity model if not found.
"""
function getRtIrtModel(s::SearchContext)
   return s.rt_irt_map
end

"""
    getRtToRefinedIrtModel(s::SearchContext, index::Integer)

Get RT → refined_iRT model for file index.
Falls back to library iRT if refined unavailable.
Returns identity model if neither exists.

# Usage
observed_refined_irt = getRtToRefinedIrtModel(context, file_idx)(scan_rt)
"""
function getRtToRefinedIrtModel(s::SearchContext, index::I) where {I<:Integer}
    if haskey(s.rt_to_refined_irt_map, index)
        return s.rt_to_refined_irt_map[index]
    elseif haskey(s.rt_irt_map, index)
        @debug "Refined iRT model not found for file $index, falling back to library iRT model"
        return s.rt_irt_map[index]
    else
        return IdentityModel()
    end
end

"""
    getRefinedIrtToRtModel(s::SearchContext, index::Integer)

Get refined_iRT → RT model for file index.
Falls back to library iRT if refined unavailable.
Returns identity model if neither exists.

# Usage
predicted_rt = getRefinedIrtToRtModel(context, file_idx)(refined_irt)
"""
function getRefinedIrtToRtModel(s::SearchContext, index::I) where {I<:Integer}
    if haskey(s.refined_irt_to_rt_map, index)
        return s.refined_irt_to_rt_map[index]
    elseif haskey(s.irt_rt_map, index)
        @debug "Refined iRT model not found for file $index, falling back to library iRT model"
        return s.irt_rt_map[index]
    else
        return IdentityModel()
    end
end

"""
   getNceModel(s::SearchContext, index::Integer)

Get NCE model for MS file index. Returns default 30 NCE model if not found.
"""
function getNceModel(s::SearchContext, index::I) where {I<:Integer} 
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

"""
    getRtToRefinedIrtMap(s::SearchContext)

Get dictionary of RT → refined_iRT conversion models for all files.
"""
getRtToRefinedIrtMap(s::SearchContext) = s.rt_to_refined_irt_map

"""
    getRefinedIrtToRtMap(s::SearchContext)

Get dictionary of refined_iRT → RT conversion models for all files.
"""
getRefinedIrtToRtMap(s::SearchContext) = s.refined_irt_to_rt_map

"""
    setRtToRefinedIrtMap!(s::SearchContext, rcm::RtConversionModel, index::Integer)

Store RT → refined_iRT conversion model for MS file index.
"""
function setRtToRefinedIrtMap!(s::SearchContext, rcm::RtConversionModel, index::I) where {I<:Integer}
    s.rt_to_refined_irt_map[index] = rcm
end

"""
    setRefinedIrtToRtMap!(s::SearchContext, rcm::RtConversionModel, index::Integer)

Store refined_iRT → RT conversion model for MS file index.
"""
function setRefinedIrtToRtMap!(s::SearchContext, rcm::RtConversionModel, index::I) where {I<:Integer}
    s.refined_irt_to_rt_map[index] = rcm
end

"""
    getIrtRefinementModel(s::SearchContext, index::Integer)

Get iRT refinement model for MS file index. Returns nothing if not found.
"""
function getIrtRefinementModel(s::SearchContext, index::I) where {I<:Integer}
    return get(s.irt_refinement_models, index, nothing)
end

"""
    setIrtRefinementModel!(s::SearchContext, model::Union{IrtRefinementModel, Nothing}, index::Integer)

Store iRT refinement model for MS file index.
"""
function setIrtRefinementModel!(s::SearchContext, model::Union{IrtRefinementModel, Nothing}, index::I) where {I<:Integer}
    s.irt_refinement_models[index] = model
end

function setPrecursorDict!(s::SearchContext, dict::Dictionary{UInt32, @NamedTuple{best_prob::Float32, best_ms_file_idx::UInt32, best_scan_idx::UInt32, best_refined_irt::Float32, mean_refined_irt::Union{Missing, Float32}, var_refined_irt::Union{Missing, Float32}, n::Union{Missing, UInt16}, mz::Float32}})
    s.precursor_dict[] = dict
end
setRtIndexPaths!(s::SearchContext, paths::Vector{String}) = (s.rt_index_paths[] = paths)
setHuberDelta!(s::SearchContext, delta::Float32) = (s.huber_delta[] = delta)
function setIrtErrors!(s::SearchContext, errs::Dictionary{Int64, Float32})
    for (k,v) in pairs(errs)
        s.irt_errors[k] = v
    end
end

#==========================================================
Method Results Storage Accessors
==========================================================#

"""
    store_results!(ctx::SearchContext, ::Type{T}, results) where T<:SearchMethod

Store results from a search method in the context.
"""
function store_results!(ctx::SearchContext, ::Type{T}, results) where T<:SearchMethod
    ctx.method_results[T] = results
    return nothing
end

"""
    get_results(ctx::SearchContext, ::Type{T}) where T<:SearchMethod

Retrieve stored results for a search method. Returns nothing if not found.
"""
function get_results(ctx::SearchContext, ::Type{T}) where T<:SearchMethod
    return get(ctx.method_results, T, nothing)
end

"""
    has_results(ctx::SearchContext, ::Type{T}) where T<:SearchMethod

Check if results exist for a search method.
"""
function has_results(ctx::SearchContext, ::Type{T}) where T<:SearchMethod
    return haskey(ctx.method_results, T)
end
