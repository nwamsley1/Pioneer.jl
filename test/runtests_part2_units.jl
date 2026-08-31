# Pioneer test suite — part 2 (unit tests, group A: foundations).
#
# Runs in its own Julia subprocess (separate from part 1 and part 3) so
# it has a clean compiler and heap state. See test/runtests.jl for the
# driver rationale.
#
# This file holds the "early" units — types, splines, FASTA, protein
# inference, mass-spec/logging utilities, BuildPionLib. Some of these
# ARE the suspected poisoners, so isolating later tests into part 3 is
# the whole point of the split.

include(joinpath(@__DIR__, "test_setup.jl"))

include("./UnitTests/test_portable_cpu_boundary.jl")

# No outer @testset wrapper — see runtests_part1_heavy.jl for why.
# Per-file @testsets each print their own summary; the driver decides
# pass/fail from this subprocess's exit code.

println("dir ", @__DIR__)

# NOTE: the four hoisted includes below stay first as a
# belt-and-suspenders defense alongside the runtests.jl subprocess
# split. See git log around 78fbd55a / 80c65e22 / 2ec91877 / c56e4b86
# for backstory.
include("./UnitTests/test_whittaker_henderson.jl")
include("./UnitTests/test_retention_time_index.jl")
include("./UnitTests/test_arrow_psm_container.jl")
include("./UnitTests/test_empty_search_results.jl")
include("./UnitTests/test_quant_trace_selection.jl")
include("./UnitTests/test_buildpionlib_index_filters.jl")
include("./UnitTests/test_enzymatic_termini_metadata.jl")

include("./UnitTests/isotopeSplines.jl")
include("./UnitTests/testIsotopesJun13.jl")
include("./UnitTests/test_uniform_spline.jl")
include("./UnitTests/test_bspline.jl")
include("./UnitTests/test_protein_inference.jl")
include("./UnitTests/test_global_protein_inference.jl")
include("./UnitTests/ChronologerPrepTests.jl")
include("./UnitTests/FastaDigestTests.jl")
include("./UnitTests/FastaEntryConstructorsTests.jl")
include("./UnitTests/PrecursorPositionOutputTests.jl")
include("./UnitTests/BuildPionLibTest.jl")
include("./UnitTests/RazoQuadModel.jl")

include("./UnitTests/MassSpecAndFilteredDataTests.jl")
include("./UnitTests/LoggingTests.jl")
include("./UnitTests/LogTruncationTests.jl")

include("./UnitTests/test_counter.jl")

# Coverage: small modules with previously-poor unit coverage.
include("./UnitTests/test_mass_error_model.jl")
include("./UnitTests/test_ms1_feature_helpers.jl")
include("./UnitTests/test_parse_mods.jl")
include("./UnitTests/test_quad_models_basic.jl")
include("./UnitTests/test_solve_poisson_mm.jl")
include("./UnitTests/test_scoring_workspace_coverage.jl")
include("./UnitTests/test_rt_alignment_utils.jl")
include("./UnitTests/test_protein_model_fit.jl")
include("./UnitTests/test_protein_mbr_features.jl")
include("./UnitTests/test_protein_ambiguous_scoring.jl")
include("./UnitTests/test_model_limits.jl")
include("./UnitTests/test_quant_bin_occupancy.jl")
include("./Routines/SearchDIA/SearchMethods/ProteinScoringSearch/test_global_protein_scoring.jl")

# Coverage data is normally flushed to .cov files by an atexit hook
# that `_exit` skips. Flush it explicitly before terminating, otherwise
# `--code-coverage=...` produces zero output and codecov reports 0%.
Base.JLOptions().code_coverage != 0 &&
    ccall(:jl_write_coverage_data, Cvoid, (Cstring,), C_NULL)

# Bypass Julia's normal shutdown via POSIX _exit(2) — atexit hooks +
# the final GC sweep both trigger more inference / compile work, and
# Julia 1.12.6 segfaults in that path once enough code has been
# compiled in this subprocess. If any @testset above failed, its
# default `finish` would have thrown a TestSetException and we'd
# never reach this line, so the only way control gets here is
# all-passed → exit code 0.
ccall(:_exit, Cvoid, (Cint,), 0)
