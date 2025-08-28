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

using Test
using Pioneer  # Load the Pioneer module for exported functions

# Import specific internal types and functions needed for testing
# These are NOT exported by Pioneer but we need them for unit tests
using Pioneer: SparseArray, ArrayDict, FragmentMatch, parseIsoXML
using Pioneer: UniSpecFragAnnotation, GenericFragAnnotation
using Pioneer: InstrumentSpecificModel, InstrumentAgnosticModel
using Pioneer: SplineCoefficientModel, RetentionTimeModel
using Pioneer: H2O, PROTON, NEUTRON, NCE_MODEL_BREAKPOINT
using Pioneer: InterpolationTypeAlias
using Pioneer: DetailedFrag, SimpleFrag, LibraryFragmentLookup
using Pioneer: DEBUG_CONSOLE_LEVEL
using Pioneer: buildDesignMatrix!  # For buildDesignMatrix.jl test
using Pioneer: IsotopeSplineModel   # For isotopeSplines.jl test
using Pioneer: getPrecursorIsotopeSet, getFragIsotopes!  # For isotopeSplines.jl test
using Pioneer: MassErrorModel, matchPeaks!, reset!  # For matchPeaks.jl test
using Pioneer: FragIndexBin, findFirstFragmentBin  # For queryFragmentIndex.jl test
using Pioneer: exponentialFragmentBinSearch, searchFragmentBin!  # For queryFragmentIndex.jl test
using Pioneer: getHigh, Counter, IndexFragment, queryFragment!  # For queryFragmentIndex.jl test
using Pioneer: UniformSpline  # For uniformBasisCubicSpline.jl test
using Pioneer: ProteinKey, PeptideKey, InferenceResult  # For test_protein_inference.jl test
using Pioneer: adjustNCE  # For ChronologerPrepTests.jl
using Pioneer: PeptideMod, matchVarMods, add_pair_indices!  # For FastaDigestTests.jl
using Pioneer: digest_sequence, getFixedMods!, countVarModCombinations  # For FastaDigestTests.jl
using Pioneer: FastaEntry, parse_fasta, PeptideSequenceSet  # For FastaDigestTests.jl
using Pioneer: buildFragmentIndex!, FragBoundModel, cleanUpLibrary  # For BuildPionLibTest.jl
using Pioneer: RazoQuadParams, simmulateQuad  # For RazoQuadModel.jl
using Pioneer: buildPionLib  # For BuildPionLibTest.jl
using Pioneer: digest_fasta, combine_shared_peptides  # For FastaDigestTests.jl  
using Pioneer: add_decoy_sequences, add_entrapment_sequences  # For FastaDigestTests.jl
using Pioneer: fillVarModStrings!, fragFilter  # For BuildPionLibTest.jl
using Pioneer: get_base_pep_id, get_base_prec_id, get_charge  # For FastaDigestTests.jl
using Pioneer: get_proteome, get_sequence, get_structural_mods  # For FastaDigestTests.jl
using Pioneer: get_isotopic_mods, get_description, get_id  # For FastaDigestTests.jl
using Pioneer: get_entrapment_pair_id  # For FastaDigestTests.jl
using Pioneer: filter_by_threshold, filter_by_multiple_thresholds  # For file operations tests
using Pioneer: getDetailedFrags, getSeqSet, getSimpleFrags, getMZ  # For BuildPionLibTest.jl
using Pioneer: getIRT, getPrecCharge, getPrecID, getPrecMZ, getScore  # For BuildPionLibTest.jl
using Pioneer: is_decoy, SplineDetailedFrag  # For FastaDigestTests.jl and BuildPionLibTest.jl
using Pioneer: PSMFileReference, TransformPipeline, add_column, sort_by, apply_pipeline!  # For file operations tests
# Note: iso_splines is loaded dynamically via parseIsoXML in RazoQuadModel.jl

# Package dependencies that tests use directly
using Arrow, ArrowTypes, ArgParse
using Base64, Base.Order
using Base.Iterators: partition
using CSV, Combinatorics, CodecZlib
using DataFrames, DataStructures, Dictionaries, Distributions
using FASTX, Interpolations, JSON, JLD2
using LinearAlgebra, LoopVectorization, LinearSolve, LightXML
using Measures, NumericalIntegration, Optim
using Plots, Polynomials, ProgressBars
using Tables, StatsPlots, SentinelArrays
using Random, StaticArrays, StatsBase, SpecialFunctions, Statistics
using EvoTrees
using MLJModelInterface: fit, predict
using KernelDensity, FastGaussQuadrature
using LaTeXStrings, Printf
using SparseArrays, Dates
using HTTP
using Interpolations: Gridded, Linear, Throw, OnGrid, Line

# Test configuration - suppress debug output during tests
Pioneer.DEBUG_CONSOLE_LEVEL[] = 0

# Test-specific directory paths
main_dir = joinpath(@__DIR__, "../src")

# InterpolationTypeAlias is imported from Pioneer module, no need to redefine

const methods_path = joinpath(@__DIR__, "Routines","LibrarySearch")
const CHARGE_ADJUSTMENT_FACTORS = Float64[1, 0.9, 0.85, 0.8, 0.75]

# H2O, PROTON, NEUTRON, NCE_MODEL_BREAKPOINT are imported from Pioneer module


const MODEL_CONFIGS = Dict(
    "unispec" => (
        annotation_type = UniSpecFragAnnotation("y1^1"),
        model_type = InstrumentSpecificModel("unispec"),
        instruments = Set(["QE","QEHFX","LUMOS","ELITE","VELOS","NONE"])
    ),
    "prosit_2020_hcd" => (
        annotation_type = GenericFragAnnotation("y1+1"), 
        model_type = InstrumentAgnosticModel("prosit_2020_hcd"),
        instruments = Set([])
    ),
    "AlphaPeptDeep" => (
        annotation_type = GenericFragAnnotation("y1+1"),
        model_type = InstrumentSpecificModel("AlphaPeptDeep"),
        instruments = Set(["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"])
    )
)


const KOINA_URLS = Dict(
    "unispec" => "https://koina.wilhelmlab.org:443/v2/models/UniSpec/infer",
    "prosit_2020_hcd" => "https://koina.wilhelmlab.org:443/v2/models/Prosit_2020_intensity_HCD/infer",
    "AlphaPeptDeep" => "https://koina.wilhelmlab.org:443/v2/models/AlphaPeptDeep_ms2_generic/infer",
    "chronologer" => "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer"
)

results_dir = joinpath(@__DIR__, "../data/ecoli_test/ecoli_test_results")
if isdir(results_dir)
    # Delete all files and subdirectories within the directory
    #for item in readdir(results_dir, join=true)
    #    rm(item, force=true, recursive=true)
    #end
end
@testset "Pioneer.jl" begin
    println("dir ", @__DIR__)
    
    #@testset "process_test_speclib" begin 
    #    @test size(ParseSpecLib(joinpath(@__DIR__, "./../data/library_test/defaultParseEmpiricalLibParams2.json")).libdf, 1)==120
    #end
    #include("./UnitTests/empiricalLibTests.jl")
    #=
    @testset "process_test" begin 
        @test SearchDIA(joinpath(@__DIR__, "../data/ecoli_test/ecoli_test_params.json"))===nothing
    end
    
    
    #Test FASTA parameter enhancement
    include("./Routines/BuildSpecLib/params/test_fasta_params.jl")
    
    # Test BuildSpecLib functionality
    include("./Routines/BuildSpecLib/test_build_spec_lib.jl")
    =#
    
    include("./UnitTests/buildDesignMatrix.jl")
    include("./UnitTests/isotopeSplines.jl")
    include("./UnitTests/matchPeaks.jl")
    include("./UnitTests/queryFragmentIndex.jl")
    include("./UnitTests/testIsotopesJun13.jl")
    include("./UnitTests/uniformBassisCubicSpline.jl")
    include("./UnitTests/test_protein_inference.jl")
    include("./UnitTests/ChronologerPrepTests.jl")
    include("./UnitTests/FastaDigestTests.jl")
    include("./UnitTests/BuildPionLibTest.jl")
    include("./utils/FileOperations/test_file_operations_suite.jl")
    include("./UnitTests/RazoQuadModel.jl")
    
end
