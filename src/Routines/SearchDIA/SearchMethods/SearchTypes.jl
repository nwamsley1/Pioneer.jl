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

See ParameterTuningSearch or MainSearch for example implementations.
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

"""
RT-bin-adaptive chromatographic tolerance.

Stores per-RT-bin tolerance values (in empirical RT space) computed from
high-confidence precursor peak widths, so that every precursor gets enough
chromatogram points regardless of gradient steepness.
"""
struct RTBinnedTolerance
    rt_edges::Vector{Float32}   # N+1 bin edges in empirical RT space
    rt_tols::Vector{Float32}    # N tolerance values (half-width, one per bin)
end

"""
    get_rt_tol(tol::RTBinnedTolerance, rt::Float32) -> Float32

Look up the RT tolerance (half-width) for the given empirical RT value.
Uses binary search on bin edges.
"""
function get_rt_tol(tol::RTBinnedTolerance, rt::Float32)::Float32
    idx = searchsortedlast(tol.rt_edges, rt)
    return tol.rt_tols[clamp(idx, 1, length(tol.rt_tols))]
end

#==========================================================
Concrete Types
==========================================================#

"""
Reference to MS data stored in Arrow files.
"""
struct ArrowTableReference{N} <: MassSpecDataReference
    file_paths::NTuple{N, String}
    file_id_to_name::NTuple{N, String}
    main_search_psms::Vector{String}
    fragment_index_matches::Vector{String}
    filtered_fragment_matches::Vector{String}
    second_pass_psms::Vector{String}
    passing_psms::Vector{String}
    passing_proteins::Vector{String}
    rt_index_paths::Vector{String}

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
            fill("", n),
            fill("", n)
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
            fill("", n),
            fill("", n)
        )
    end

end



"""
Basic search data structure for library searches.
Contains pre-allocated arrays and intermediate data structures.
"""
mutable struct SimpleLibrarySearch{I<:IsotopeSplineModel} <: SearchDataStructures
    # Mass-error sample buffer — output of run_fused_masserr! during
    # ParameterTuning, consumed by fit_mass_err_model + friends.
    mass_err_samples::Vector{MassErrSample}

    # Indexing and scoring
    id_to_col::SparsePrecMap{UInt16}
    iso_splines::I
    
    # PSM scoring.
    # MainSearch buffers — `main_unscored_psms` (per-(scan, precursor)
    # accumulator written by apply_main_scoring!) and `main_search_scored_psms`
    # (per-scan rows produced by Score!) carry the full MainSearch feature set
    # including fragment-chromatogram captures.
    main_unscored_psms::Vector{MainUnscoredPSM{Float32}}
    main_search_scored_psms::Vector{MainSearchScoredPSM{Float32, Float16}}
    main_search_spectral_scores::Vector{SpectralScoresMainSearch{Float16, Float32}}
    # Tuning buffers — slim variants used by ParameterTuning, QuadTuning, and
    # IntegrateChromatograms paths. Same structural shape minus the
    # MainSearch-only fragment-chromatogram fields (frag1..8_int,
    # matched_rank_mask, rank1/top3/top5_matched).
    tuning_unscored_psms::Vector{TuningUnscoredPSM{Float32}}
    tuning_scored_psms::Vector{TuningScoredPSM{Float32, Float16}}

    # Working arrays
    prec_ids::Vector{UInt32}
    precursor_weights::SparsePrecMap{Float32}
    temp_weights::Vector{Float32}
    colnorm2::Vector{Float32}
    residuals::Vector{Float32}
    mu::Vector{Float32}
    observed::Vector{Float32}
    isotopes::Vector{Float32}
    precursor_transmission::Vector{Float32}

    # Per-thread scratch for the fused scan path (used when
    # `MainSearchParameters.use_fused_scan == true`).
    fused_scratch::FusedScratch
    Hs_fused::SparseArrayFused{UInt32, Float32}
    # Per-scan precomputed peak window (SoA — each array SIMD-loadable):
    # `scan_corrected_mz[j]` = bias-corrected center m/z (used as bsearch anchor
    # and for ppm_err on matches). `scan_obs_low[j]` / `scan_obs_high[j]` =
    # intensity-aware Da tolerance bounds for matching (classic semantics:
    # theoretical_mz must satisfy `obs_low[j] ≤ theoretical ≤ obs_high[j]`).
    scan_corrected_mz::Vector{Float32}
    scan_obs_low::Vector{Float32}
    scan_obs_high::Vector{Float32}
end

"""
Primary search context holding all data structures and state for search execution.
"""
# No thread-count type parameter. It used to carry `N = length(temp_structures)`, which typed
# nothing -- `temp_structures` is an AbstractVector, not an NTuple{N} -- and was never dispatched on
# or read. Its only effect was to give every method touching a SearchContext a separate
# specialisation per thread count, which no precompile statement file can cover because the value
# comes from the user's machine. See the note on ArrowTableReference below.
mutable struct SearchContext{L<:SpectralLibrary,M<:MassSpecDataReference}
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
    mass_error_model::Dict{Int64, AbstractMassErrorModel}
    ms1_mass_error_model::Dict{Int64, AbstractMassErrorModel}
    #rt_to_irt_model::Dict{Int64, RtConversionModel}
    nce_model::Dict{Int64, NceModel}
    huber_delta::Base.Ref{Float32}
    deconvolution_stop_tolerance::Base.Ref{Float32}
    # Results and paths
    irt_rt_map::Dict{Int64, RtConversionModel}
    rt_irt_map::Dict{Int64, RtConversionModel}
    precursor_dict::Base.Ref{Dictionary}
    rt_index_paths::Base.Ref{Vector{String}}
    irt_errors::Dict{Int64, Float32}
    rt_tolerances::Dict{Int64, RTBinnedTolerance}
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

    # BitVec pre-filter (256-entry Bool table per file, trained by BitVecCalibration)
    bitvec_filter::Dict{Int64, Vector{Bool}}
    # Per-file bitmask-pattern rank by target/decoy excess. Index is
    # Int(mask)+1, lower rank is better.
    bitvec_excess_rank::Dict{Int64, Vector{UInt16}}

    # Reusable per-thread scratch buffers for searchFragmentIndexPartitionMajorHinted
    # (see FragIndexScratch in src/structs/Counter.jl). Lifted out of the per-call
    # path so BitVec's adaptive while loop doesn't re-allocate ~50–100 MB per batch.
    frag_index_scratch::FragIndexScratch

    # Constructor
    function SearchContext(
        spec_lib::L,
        temp_structures::AbstractVector{<:SearchDataStructures},
        mass_spec_data_reference::M,
        n_threads::Int64,
        n_precursors::Int64,
        buffer_size::Int64
    ) where {L<:SpectralLibrary,M<:MassSpecDataReference}
        new{L,M}(
            spec_lib, temp_structures, mass_spec_data_reference,
            Ref{String}(), Ref{String}(), Ref{String}(), Ref{String}(),Ref{String}(),
            Dict{Int64, QuadTransmissionModel}(),
            Dict{Int64, AbstractMassErrorModel}(),
            Dict{Int64, AbstractMassErrorModel}(),
            Dict{Int64, NceModel}(),
            Ref(300.0f0),
            10.0f0,
            Dict{Int64, RtConversionModel}(), 
            Dict{Int64, RtConversionModel}(), 
            Ref{Dictionary}(), 
            Ref{Vector{String}}(),
            Dict{Int64, Float32}(),
            Dict{Int64, RTBinnedTolerance}(),
            Dict{UInt32, Float32}(),
            Ref{Any}(), Ref(Dict{ProteinKey, Float32}()), Ref(Dict{Tuple{String,Bool,UInt8}, Float32}()), Ref{Any}(),
            Dict{Type{<:SearchMethod}, Any}(),  # Initialize method_results
            n_threads, n_precursors, buffer_size,
            0, 0, 1.0f0,  # Initialize library stats with defaults
            Dict{Int64, Vector{Bool}}(),  # Initialize bitvec_filter
            Dict{Int64, Vector{UInt16}}(), # Initialize bitvec_excess_rank
            FragIndexScratch(n_threads),  # reusable frag-index scratch buffers
        )
    end
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
getMainSearchPsms(ref::ArrowTableReference, index::Int) = ref.main_search_psms[index]
getFragmentIndexMatches(ref::ArrowTableReference, index::Int) = ref.fragment_index_matches[index]
getFilteredFragmentMatches(ref::ArrowTableReference, index::Int) = ref.filtered_fragment_matches[index]
getSecondPassPsms(ref::ArrowTableReference, index::Int) = ref.second_pass_psms[index]

"""
    getSecondPassPsmsFold(ref::ArrowTableReference, index::Int, fold::UInt8) -> String

Get the path to a specific CV fold's second pass PSM file.

The second_pass_psms field stores the base path (without fold suffix).
This function constructs the full path for a specific fold.

# Arguments
- `ref`: ArrowTableReference containing file paths
- `index`: MS file index
- `fold`: CV fold number (0 or 1)

# Returns
- Full path to the fold-specific Arrow file, or empty string if base path is empty
"""
function getSecondPassPsmsFold(ref::ArrowTableReference, index::Int, fold::UInt8)
    base = ref.second_pass_psms[index]
    return isempty(base) ? "" : "$(base)_fold$(fold).arrow"
end
getPassingPsms(ref::ArrowTableReference, index::Int) = ref.passing_psms[index]
getPassingProteins(ref::ArrowTableReference, index::Int) = ref.passing_proteins[index]
getRtIndex(ref::ArrowTableReference, index::Int) = ref.rt_index_paths[index]
getParsedFileNames(ref::ArrowTableReference) = ref.file_id_to_name

getFilePaths(ref::ArrowTableReference) = ref.file_paths

getFileIdToName(ref::ArrowTableReference) = ref.file_id_to_name
getMainSearchPsms(ref::ArrowTableReference) = ref.main_search_psms
getFragmentIndexMatches(ref::ArrowTableReference) = ref.fragment_index_matches
getFilteredFragmentMatches(ref::ArrowTableReference) = ref.filtered_fragment_matches
getSecondPassPsms(ref::ArrowTableReference) = ref.second_pass_psms
getPassingPsms(ref::ArrowTableReference) = ref.passing_psms
getPassingProteins(ref::ArrowTableReference) = ref.passing_proteins
getRtIndex(ref::ArrowTableReference) = ref.rt_index_paths


# Setter methods
setFileIdToName!(ref::ArrowTableReference, index::Int, value::String) = ref.file_id_to_name[index] = value
setMainSearchPsms!(ref::ArrowTableReference, index::Int, value::String) = ref.main_search_psms[index] = value
setFragmentIndexMatches!(ref::ArrowTableReference, index::Int, value::String) = ref.fragment_index_matches[index] = value
setFilteredFragmentMatches!(ref::ArrowTableReference, index::Int, value::String) = ref.filtered_fragment_matches[index] = value
setSecondPassPsms!(ref::ArrowTableReference, index::Int, value::String) = ref.second_pass_psms[index] = value
setPassingPsms!(ref::ArrowTableReference, index::Int, value::String) = ref.passing_psms[index] = value
setPassingProteins!(ref::ArrowTableReference, index::Int, value::String) = ref.passing_proteins[index] = value
setRtIndex!(ref::ArrowTableReference, index::Int, value::String) = ref.rt_index_paths[index] = value

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
getMinIndexSearchScore(fsp::SearchParameters)   = fsp.min_index_search_score
getPrefilterMinScanCount(::FragmentIndexSearchParameters) = 0
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

getMassErrSamples(s::SearchDataStructures) = s.mass_err_samples
getIdToCol(s::SearchDataStructures) = s.id_to_col
getIsoSplines(s::SearchDataStructures) = s.iso_splines

getMainUnscoredPsms(s::SearchDataStructures) = s.main_unscored_psms

getMainSearchScoredPsms(s::SearchDataStructures) = s.main_search_scored_psms
getMainSearchSpectralScores(s::SearchDataStructures) = s.main_search_spectral_scores

# Tuning-path buffer getters (slim PSM types).
getTuningUnscoredPsms(s::SearchDataStructures) = s.tuning_unscored_psms
getTuningScoredPsms(s::SearchDataStructures) = s.tuning_scored_psms

getPrecIds(s::SearchDataStructures) = s.prec_ids
getWeights(s::SearchDataStructures) = s.weights
getResiduals(s::SearchDataStructures) = s.residuals
getIsotopes(s::SearchDataStructures) = s.isotopes
getPrecursorTransmission(s::SearchDataStructures) = s.precursor_transmission
getFusedScratch(s::SearchDataStructures) = s.fused_scratch
getHsFused(s::SearchDataStructures) = s.Hs_fused
getScanCorrectedMz(s::SearchDataStructures) = s.scan_corrected_mz
getScanObsLow(s::SearchDataStructures) = s.scan_obs_low
getScanObsHigh(s::SearchDataStructures) = s.scan_obs_high
getTuningResults(s::SearchDataStructures) = s.tuning_results
getTempWeights(s::SimpleLibrarySearch) = s.temp_weights
getColNorm2(s::SimpleLibrarySearch) = s.colnorm2
getPrecursorWeights(s::SimpleLibrarySearch) = s.precursor_weights
getMu(s::SimpleLibrarySearch) = s.mu
getObserved(s::SimpleLibrarySearch) = s.observed

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
getRtTolerances(s::SearchContext) = s.rt_tolerances
getRtTolerance(s::SearchContext, ms_file_idx::Int64) = s.rt_tolerances[ms_file_idx]
getHuberDelta(s::SearchContext) = s.huber_delta[]
# Use library iRT array directly — O(1) indexing, no Dict overhead
getPredIrt(s::SearchContext) = getIrt(getPrecursors(getSpecLib(s)))
getPredIrt(s::SearchContext, prec_idx::Int64) = getIrt(getPrecursors(getSpecLib(s)))[prec_idx]
getPredIrt(s::SearchContext, prec_idx::UInt32) = getIrt(getPrecursors(getSpecLib(s)))[prec_idx]
# irt_obs Dict on SearchContext is now unused; setPredIrt! removed
getLibraryTargetCount(s::SearchContext) = s.n_library_targets
getLibraryDecoyCount(s::SearchContext) = s.n_library_decoys
getBitVecFilter(s::SearchContext, ms_file_idx::Int64) = get(s.bitvec_filter, ms_file_idx, nothing)
getBitVecExcessRanks(s::SearchContext, ms_file_idx::Int64) = get(s.bitvec_excess_rank, ms_file_idx, nothing)
getFragIndexScratch(s::SearchContext) = s.frag_index_scratch
setBitVecFilter!(s::SearchContext, ms_file_idx::Int64, filter::Vector{Bool}) = (s.bitvec_filter[ms_file_idx] = filter)
setBitVecExcessRanks!(s::SearchContext, ms_file_idx::Int64, ranks::Vector{UInt16}) = (s.bitvec_excess_rank[ms_file_idx] = ranks)
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
setMassErrorModel!(s::SearchContext, index::I, model::AbstractMassErrorModel) where {I<:Integer} = (s.mass_error_model[index] = model)

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


setMs1MassErrorModel!(s::SearchContext, index::I, model::AbstractMassErrorModel) where {I<:Integer} = (s.ms1_mass_error_model[index] = model)

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
