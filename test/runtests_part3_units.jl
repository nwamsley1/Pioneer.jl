# Pioneer test suite — part 3 (unit tests, group B: containers + algorithms).
#
# Runs in its own Julia subprocess (separate from parts 1 and 2) so it
# has a clean compiler and heap state. See test/runtests.jl for the
# driver rationale.
#
# This file holds the "later" units — PSM containers, scoring,
# quantification, fused-scan pipeline, file-operations, partitioned
# fragment index, parallel utilities, and BuildSpecLib decoy / synthetic
# Koina pieces. The split point (between part 2 and part 3) is
# deliberately just before `test_psm_container.jl`, which kept being
# the next victim of upstream-poisoned state in the single-process
# suite. With a fresh subprocess at this boundary, those crashes go
# away.

include(joinpath(@__DIR__, "test_setup.jl"))

# No outer @testset wrapper — see runtests_part1_heavy.jl for why.

println("dir ", @__DIR__)

include("./UnitTests/test_psm_container.jl")
include("./UnitTests/test_scoring_workspace.jl")
include("./UnitTests/test_fast_df_sort.jl")

# Quantification tests
include("./UnitTests/test_maxLFQ.jl")
include("./UnitTests/test_normalizeQuant.jl")
include("./UnitTests/test_download_speclib.jl")

# FDR/q-value utilities
include("./UnitTests/test_fdr_utilities.jl")

# Fused per-precursor scan pipeline
include("./UnitTests/test_fused_prec_filters.jl")
include("./UnitTests/test_run_fused.jl")
include("./UnitTests/test_fragment_isotope_trace_intensities.jl")
include("./UnitTests/test_sort_detailed_fragments.jl")
include("./UnitTests/test_feature_finiteness.jl")
include("./UnitTests/test_bitvec_rank_features.jl")
include("./UnitTests/test_fragment_correlation_summaries.jl")
include("./UnitTests/test_chromatogram_integration_trace_order.jl")
include("./UnitTests/test_chromatogram_integration_bounds.jl")
include("./UnitTests/test_huber_tuning_global.jl")

# FileOperations focused tests
include("./utils/FileOperations/io/test_arrow_operations_basic.jl")
include("./utils/FileOperations/io/test_safe_file_ops.jl")
include("./utils/FileOperations/core/test_core_references_basic.jl")
include("./utils/FileOperations/streaming/test_stream_sorted_merge_basic.jl")

# Partitioned fragment index tests
include("./UnitTests/partitionedFragmentIndex.jl")
include("./UnitTests/buildPartitionedIndex.jl")
include("./UnitTests/test_build_fragment_index_exact.jl")
include("./UnitTests/test_synthetic_koina.jl")
include("./UnitTests/test_koina_http_retry.jl")
include("./Routines/BuildSpecLib/fragments/test_get_frag_bounds.jl")
include("./integration/test_buildspeclib_synthetic.jl")
include("./Routines/BuildSpecLib/test_build_determinism.jl")
include("./UnitTests/test_buildspeclib_filter_equivalence.jl")

# Parallel utilities
include("./UnitTests/test_parallel_utils.jl")

# Spectral deconvolution - Poisson MM
# NOTE: test_poissonMM.jl is excluded from the suite. It's a
# developer-local diagnostic that hardcodes a path under
# /Users/nathanwamsley/Data/For_Figures/... for serialized .dat
# fixtures, AND it imports `Pioneer.SparseArray`, which was removed
# when the deconvolution path was unified onto SparseArrayFused
# (the file src/structs/SparseArray.jl now defines only the
# AbstractSparseDesignMatrix supertype, not a concrete SparseArray).
# Re-enable if either the type is reintroduced or the test is
# ported to SparseArrayFused with self-contained fixtures.
# include("./UnitTests/test_poissonMM.jl")

# BuildSpecLib - decoy generation
include("./UnitTests/test_diann_decoys.jl")

# QuadTuning utilities
include("./UnitTests/test_summarize_precursor.jl")

# PrecursorScoringSearch
include("./UnitTests/test_run_similarity.jl")
include("./Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/test_global_precursor_scoring.jl")
include("./Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/test_wide_window_features.jl")
include("./Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/test_mbr_file_aware_decoys.jl")
include("./Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/test_mbr_empirical_spectra_hellinger.jl")
include("./Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/test_mbr_postintegration_pipeline.jl")
include("./UnitTests/test_mainsearch_irt_refinement.jl")
include("./UnitTests/test_scoring_semisupervised.jl")

# FileOperations pipeline
include("./utils/FileOperations/pipeline/test_pipeline_filtering.jl")

# Coverage data is normally flushed to .cov files by an atexit hook
# that `_exit` skips. Flush it explicitly before terminating.
Base.JLOptions().code_coverage != 0 &&
    ccall(:jl_write_coverage_data, Cvoid, (Cstring,), C_NULL)

# Bypass Julia's normal shutdown — see runtests_part2_units.jl for why.
ccall(:_exit, Cvoid, (Cint,), 0)
