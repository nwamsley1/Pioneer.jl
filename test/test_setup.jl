# Shared test setup — imports + configuration.
#
# Each runtests_part*.jl includes this file so its subprocess has the
# same names available. Keeping the imports in one place avoids drift
# between parts.

using Test
using Pioneer  # Load the Pioneer module for exported functions

# Import specific internal types and functions needed for testing.
# These are NOT exported by Pioneer but we need them for unit tests.
using Pioneer: parseIsoXML
using Pioneer: UniSpecFragAnnotation, GenericFragAnnotation
using Pioneer: SplineCoefficientModel, RetentionTimeModel
using Pioneer: H2O, PROTON, C13_C12_MASS_DIFF
using Pioneer: InterpolationTypeAlias
using Pioneer: DetailedFrag, SimpleFrag, LibraryFragmentLookup
using Pioneer: DEBUG_CONSOLE_LEVEL
using Pioneer: IsotopeSplineModel
using Pioneer: getPrecursorIsotopeSet, getFragIsotopes!
using Pioneer: getHigh, Counter, IndexFragment
using Pioneer: UniformSpline
using Pioneer: ProteinKey, PeptideKey, InferenceResult
using Pioneer: PeptideMod, matchVarMods, add_pair_indices!
using Pioneer: digest_sequence, getFixedMods!, countVarModCombinations
using Pioneer: FastaEntry, parse_fasta, PeptideSequenceSet
using Pioneer: buildFragmentIndex!, FragBoundModel, cleanUpLibrary, get_fragment_bounds
using Pioneer: RazoQuadParams, simmulateQuad, fitRazoQuadModel, MergeBins
using Pioneer: buildPionLib
using Pioneer: digest_fasta, combine_shared_peptides
using Pioneer: add_decoy_sequences, add_entrapment_sequences
using Pioneer: fillVarModStrings!, fragFilter
using Pioneer: get_base_pep_id, get_charge
using Pioneer: make_koina_request, prepare_koina_batch, parse_koina_batch
using Pioneer: KoinaBatchResult
using Pioneer: get_proteome, get_sequence, get_structural_mods
using Pioneer: get_isotopic_mods, get_description, get_id
using Pioneer: get_entrapment_pair_id
using Pioneer: get_gene, get_protein, get_organism
using Pioneer: filter_by_threshold, filter_by_multiple_thresholds
using Pioneer: compute_diann_mutation_shifts, create_diann_decoy_fragments
using Pioneer: extract_mod_positions, DIANN_MUTATION_TABLE, AA_to_mass
using Pioneer: solvePoissonMM_fast!, initMu!
using Pioneer: logodds
using Pioneer: getDetailedFrags, getSeqSet, getSimpleFrags, getMZ
using Pioneer: getIRT, getPrecCharge, getPrecID, getPrecMZ, getScore
using Pioneer: is_decoy, SplineDetailedFrag
using Pioneer: PSMFileReference, TransformPipeline, add_column, sort_by, apply_pipeline!

# Package dependencies that tests use directly.
using Arrow, ArrowTypes, ArgParse
using Base64, Base.Order
using Base.Iterators: partition
using CSV, Combinatorics, CodecZlib
using DataFrames, Dictionaries, Distributions
import DataStructures  # qualified to avoid reset! conflict
using FASTX, Interpolations, JSON, JLD2
using LinearAlgebra, LoopVectorization, LinearSolve, LightXML
using Measures, NumericalIntegration, Optim
using Plots, Polynomials, ProgressBars
using Tables, StatsPlots, SentinelArrays
using Random, StaticArrays, StatsBase, SpecialFunctions, Statistics
using LightGBM
using MLJModelInterface: fit, predict
using KernelDensity, FastGaussQuadrature
using LaTeXStrings, Printf
using SparseArrays, Dates
using HTTP
using Interpolations: Gridded, Linear, Throw, OnGrid, Line

# Test configuration — suppress debug output during tests.
Pioneer.DEBUG_CONSOLE_LEVEL[] = 0

# Test-specific directory paths.
const main_dir = joinpath(@__DIR__, "../src")

# Windows-robust recursive removal for test cleanup. BuildSpecLib tests leave
# memory-mapped Arrow/.jls files (inside .poin dirs) whose OS handles Windows
# keeps open until GC runs the finalizers, so a plain `rm(...; recursive=true)`
# throws ENOTEMPTY/EACCES and turns into a spurious test *error* that aborts the
# whole test part. GC to release the mmaps, retry with backoff, and on final
# failure warn (never throw) so cleanup can never fail the suite. On POSIX the
# first rm succeeds immediately.
function safe_rmdir(path::AbstractString; max_tries::Int=5)
    (isdir(path) || isfile(path)) || return nothing
    for i in 1:max_tries
        try
            GC.gc(); GC.gc()
            rm(path; recursive=true, force=true)
            return nothing
        catch err
            if i == max_tries
                @warn "safe_rmdir: could not remove $path after $max_tries attempts" err
                return nothing
            end
            sleep(0.2 * i)
        end
    end
    return nothing
end
