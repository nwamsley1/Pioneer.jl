# Restore Staged SearchDIA Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make current `develop` use a staged SearchDIA identification pipeline again while preserving the current tuning stack: current parameter/NCE-related tuning behavior, current quadrupole tuning, BitVec calibration, first-pass search, Huber tuning, second-pass search, precursor scoring, chromatogram integration, protein inference/scoring, and MaxLFQ output.

**Architecture:** Start from current `develop` and replace the executed fused `MainSearch` identification path with the staged first-pass/Huber/second-pass stack from `a33f96466d2d38538a5da919fe8cc28effa76ecd`, while preserving current `ParameterTuningSearch`, current NCE-related tuning behavior, current `QuadTuningSearch`, and current `BitVecCalibrationSearch`. Do not add a runtime "legacy" switch; the fused version remains available by checking out its existing branch/history. Keep current tuning infrastructure, including fused matching utilities used by tuning, and avoid library-prediction/schema changes unless the restored staged search cannot read current libraries without a compatibility fix.

**Tech Stack:** Julia, Pioneer.jl SearchDIA internals, Arrow/DataFrames/JLD2 library artifacts, LightGBM/percolator-style scoring, probit regression, Huber-loss deconvolution, `Pkg.test`, targeted unit tests, and a local compatible SearchDIA run using `/private/tmp/search.json`.

---

## Design Decision

This plan supersedes the earlier "roll `develop` back to `a33f964`" plan.

The new strategy is:

- Keep current `develop` as the base.
- Keep the current parameter/NCE/quad tuning stack.
- Restore the staged first-pass, Huber tuning, second-pass, and old PSM scoring implementation removed in PR #375.
- Make the staged pipeline canonical in `SearchDIA(params_path)`.
- Remove or stop loading the fused-only `MainSearch` identification path.
- Keep `BitVecCalibrationSearch` as a canonical SearchDIA step.
- Do not provide a config flag or runtime switch for fused search.
- Test staged search behavior one layer at a time before running a local compatible end-to-end search.

This is better for the current goal because the suspected regression is concentrated in the fused first/second search replacement, while other post-`a33f964` work may still be worth keeping.

## Current Versus Target Pipeline

Current `develop` pipeline:

```julia
[
    ("Parameter Tuning", ParameterTuningSearch()),
    ("Quadrupole Tuning", QuadTuningSearch()),
    ("BitVec Calibration", BitVecCalibrationSearch()),
    ("Main Search", MainSearch()),
    ("Precursor Scoring", PrecursorScoringSearch()),
    ("Chromatogram Integration", IntegrateChromatogramSearch()),
    ("Protein Inference", ProteinInferenceSearch()),
    ("Protein Scoring", ProteinScoringSearch()),
    ("Quantification & Output", MaxLFQSearch())
]
```

Target canonical pipeline:

```julia
[
    ("Parameter Tuning", ParameterTuningSearch()),
    ("Quadrupole Tuning", QuadTuningSearch()),
    ("BitVec Calibration", BitVecCalibrationSearch()),
    ("First Pass Search", FirstPassSearch()),
    ("Huber Tuning", HuberTuningSearch()),
    ("Second Pass Search", SecondPassSearch()),
    ("Precursor Scoring", PrecursorScoringSearch()),
    ("Chromatogram Integration", IntegrateChromatogramSearch()),
    ("Protein Inference", ProteinInferenceSearch()),
    ("Protein Scoring", ProteinScoringSearch()),
    ("Quantification & Output", MaxLFQSearch())
]
```

## Scope

In scope:

- Restore `FirstPassSearch`, `HuberTuningSearch`, and `SecondPassSearch`.
- Restore old fragment-index search helpers: `matchPeaks`, `queryFragmentIndex`, `buildDesignMatrix`, and `selectTransitions/*`.
- Restore old probit/percolator-style scoring support: `percolatorSortOf`, `percolator_generic`, `scoring_workspace`, `scoring_config`, `pairing`, `training_selection`, `iteration_scheme`, `mbr_update`, `spectralLinearRegression`, and `ftrUtilities`.
- Rewire `SearchDIA.jl` so the staged pipeline is the only canonical executed pipeline.
- Preserve `BitVecCalibrationSearch` in the canonical pipeline, unless implementation proves it depends directly on fused `MainSearch` internals.
- Reconcile `SearchTypes.jl`, `importScripts.jl`, and downstream scoring/integration expectations with the staged outputs.
- Add tests that fail while `MainSearch` is canonical and pass once staged search is canonical.
- Run fast compile/unit tests before slow regression.

Out of scope:

- A broad rollback of every post-`a33f964` commit.
- A runtime switch between fused and staged search.
- Restoring old `ParameterTuningSearch`, old `NceTuningSearch`, or old `QuadTuningSearch`.
- Removing fused matching utilities that current tuning depends on.
- Changing library prediction, BuildSpecLib output schema, or `.poin` library format unless required for staged-search compatibility.
- Removing current protein-inference/species/decoy-output changes unless they are incompatible with staged search.
- Reworking docs, release scripts, or unrelated dependency bounds.

## File Map

Restore from `a33f96466d2d38538a5da919fe8cc28effa76ecd`:

```text
src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl
src/Routines/SearchDIA/SearchMethods/FirstPassSearch/getBestPrecursorsAccrossRuns.jl
src/Routines/SearchDIA/SearchMethods/FirstPassSearch/utils.jl
src/Routines/SearchDIA/SearchMethods/HuberTuningSearch/HuberTuningSearch.jl
src/Routines/SearchDIA/SearchMethods/HuberTuningSearch/utils.jl
src/Routines/SearchDIA/SearchMethods/SecondPassSearch/SecondPassSearch.jl
src/Routines/SearchDIA/SearchMethods/SecondPassSearch/utils.jl
src/Routines/SearchDIA/CommonSearchUtils/buildDesignMatrix.jl
src/Routines/SearchDIA/CommonSearchUtils/matchPeaks.jl
src/Routines/SearchDIA/CommonSearchUtils/queryFragmentIndex.jl
src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/fillTransitionList.jl
src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/massErrEstimationStrategy.jl
src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/quadEstimationSelection.jl
src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/rtIndexTransitionSelection.jl
src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/selectTransitions.jl
src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/standardTransitionSelection.jl
src/utils/ML/ftrUtilities.jl
src/utils/ML/mbr_update.jl
src/utils/ML/pairing.jl
src/utils/ML/percolatorSortOf.jl
src/utils/ML/percolator_generic.jl
src/utils/ML/scoring_config.jl
src/utils/ML/scoring_workspace.jl
src/utils/ML/spectralLinearRegression.jl
src/utils/ML/training_selection.jl
```

Reconcile manually:

```text
src/Routines/SearchDIA.jl
src/importScripts.jl
src/Routines/SearchDIA/SearchMethods/SearchMethods.jl
src/Routines/SearchDIA/SearchMethods/SearchTypes.jl
src/Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/PrecursorScoringSearch.jl
src/Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/model_config.jl
src/Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/score_psms.jl
src/Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/scoring_interface.jl
src/Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/utils.jl
src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/IntegrateChromatogramsSearch.jl
src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/utils.jl
src/Routines/SearchDIA/ParseInputs/loadSpectralLibrary.jl
```

Remove or stop loading after staged pipeline compiles:

```text
src/Routines/SearchDIA/SearchMethods/MainSearch/MainSearch.jl
src/Routines/SearchDIA/SearchMethods/MainSearch/deconvolution.jl
src/Routines/SearchDIA/SearchMethods/MainSearch/features.jl
src/Routines/SearchDIA/SearchMethods/MainSearch/prescore_aggregation.jl
src/Routines/SearchDIA/SearchMethods/MainSearch/scoring.jl
src/Routines/SearchDIA/SearchMethods/MainSearch/types.jl
src/Routines/SearchDIA/SearchMethods/MainSearch/utils.jl
src/Routines/SearchDIA/process_scans.jl
src/Routines/SearchDIA/process_scans_fused.jl
```

Keep current tuning stack and dependencies; modify only for narrow staged-search compatibility:

```text
src/Routines/SearchDIA/CommonSearchUtils/fusedMatch.jl
src/Routines/SearchDIA/CommonSearchUtils/fusedScan.jl
src/Routines/SearchDIA/CommonSearchUtils/getFragIsotopes.jl
src/Routines/SearchDIA/CommonSearchUtils/isotope_utils.jl
src/utils/SpectralDeconvolution/*
src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/*
src/Routines/SearchDIA/SearchMethods/QuadTuningSearch/*
src/Routines/SearchDIA/SearchMethods/BitVecCalibrationSearch/*
```

Keep unless directly incompatible:

```text
src/Routines/SearchDIA/SearchMethods/ProteinInferenceSearch/*
src/Routines/SearchDIA/SearchMethods/ProteinScoringSearch/*
src/Routines/SearchDIA/SearchMethods/MaxLFQSearch/*
src/utils/maxLFQ.jl
src/utils/proteinInference.jl
.github/workflows/regression.yml
test/runtests.jl
test/runtests_part1_heavy.jl
test/runtests_part2_units.jl
test/runtests_part3_units.jl
```

## Task 1: Create A Dedicated Restoration Branch

**Files:**
- No source edits

- [ ] **Step 1: Confirm the starting branch**

Run:

```bash
git status --short --branch
```

Expected:

```bash
## develop...origin/develop
?? docs/superpowers/
```

- [ ] **Step 2: Create the working branch**

Run:

```bash
git switch -c hotfix/restore-staged-searchdia-pipeline
```

Expected:

```bash
Switched to a new branch 'hotfix/restore-staged-searchdia-pipeline'
```

- [ ] **Step 3: Save a machine-readable changed-file inventory**

Run:

```bash
git diff --name-status a33f96466d2d38538a5da919fe8cc28effa76ecd..HEAD -- src/Routines/SearchDIA src/utils/ML src/structs > /tmp/searchdia-overhaul-name-status.txt
```

Expected: `/tmp/searchdia-overhaul-name-status.txt` exists and lists the search/scoring removals and additions.

## Task 2: Add A Failing Pipeline-Order Test

**Files:**
- Create: `test/UnitTests/test_searchdia_pipeline_order.jl`
- Modify: `test/runtests_part2_units.jl`

- [ ] **Step 1: Add the test file**

Create `test/UnitTests/test_searchdia_pipeline_order.jl`:

```julia
using Test
using Pioneer

@testset "canonical SearchDIA pipeline is staged" begin
    @test isdefined(Pioneer, :canonical_search_pipeline)
    names = first.(Pioneer.canonical_search_pipeline())

    @test names == [
        "Parameter Tuning",
        "Quadrupole Tuning",
        "BitVec Calibration",
        "First Pass Search",
        "Huber Tuning",
        "Second Pass Search",
        "Precursor Scoring",
        "Chromatogram Integration",
        "Protein Inference",
        "Protein Scoring",
        "Quantification & Output",
    ]

    method_types = typeof.(last.(Pioneer.canonical_search_pipeline()))
    @test Pioneer.MainSearch ∉ method_types
    @test Pioneer.BitVecCalibrationSearch ∈ method_types
    @test Pioneer.FirstPassSearch ∈ method_types
    @test Pioneer.HuberTuningSearch ∈ method_types
    @test Pioneer.SecondPassSearch ∈ method_types
end
```

- [ ] **Step 2: Include it in the part-2 unit suite**

Add this line to `test/runtests_part2_units.jl` after the basic utility tests:

```julia
include("./UnitTests/test_searchdia_pipeline_order.jl")
```

- [ ] **Step 3: Run the test and verify it fails for the right reason**

Run:

```bash
julia --startup-file=no --project=. test/UnitTests/test_searchdia_pipeline_order.jl
```

Expected: failure because `canonical_search_pipeline` and/or staged search types are not defined on current `develop`.

- [ ] **Step 4: Commit the failing test**

Run:

```bash
git add test/UnitTests/test_searchdia_pipeline_order.jl test/runtests_part2_units.jl
git commit -m "test: assert staged SearchDIA pipeline order"
```

## Task 3: Restore Deleted Staged-Search Source Files

**Files:**
- Restore the files listed in "Restore from `a33f964...`"

- [ ] **Step 1: Restore search-method directories**

Run:

```bash
git restore --source=a33f96466d2d38538a5da919fe8cc28effa76ecd -- \
  src/Routines/SearchDIA/SearchMethods/FirstPassSearch \
  src/Routines/SearchDIA/SearchMethods/HuberTuningSearch \
  src/Routines/SearchDIA/SearchMethods/SecondPassSearch
```

Expected: the three restored staged-identification directories exist in the working tree. Do not restore `NceTuningSearch`; current NCE-related tuning stays with the current tuning stack.

- [ ] **Step 2: Restore legacy common search utilities**

Run:

```bash
git restore --source=a33f96466d2d38538a5da919fe8cc28effa76ecd -- \
  src/Routines/SearchDIA/CommonSearchUtils/buildDesignMatrix.jl \
  src/Routines/SearchDIA/CommonSearchUtils/matchPeaks.jl \
  src/Routines/SearchDIA/CommonSearchUtils/queryFragmentIndex.jl \
  src/Routines/SearchDIA/CommonSearchUtils/selectTransitions
```

Expected: the restored files exist and are tracked as additions relative to current `develop`.

- [ ] **Step 3: Restore legacy scoring utilities**

Run:

```bash
git restore --source=a33f96466d2d38538a5da919fe8cc28effa76ecd -- \
  src/utils/ML/ftrUtilities.jl \
  src/utils/ML/mbr_update.jl \
  src/utils/ML/pairing.jl \
  src/utils/ML/percolatorSortOf.jl \
  src/utils/ML/percolator_generic.jl \
  src/utils/ML/scoring_config.jl \
  src/utils/ML/scoring_workspace.jl \
  src/utils/ML/spectralLinearRegression.jl \
  src/utils/ML/training_selection.jl
```

Expected: all listed legacy scoring files exist in `src/utils/ML/`.

- [ ] **Step 4: Commit the raw restore**

Run:

```bash
git add src/Routines/SearchDIA/SearchMethods/FirstPassSearch src/Routines/SearchDIA/SearchMethods/HuberTuningSearch src/Routines/SearchDIA/SearchMethods/SecondPassSearch src/Routines/SearchDIA/CommonSearchUtils src/utils/ML
git commit -m "refactor: restore staged SearchDIA source files"
```

## Task 4: Rewire Includes For Staged Search

**Files:**
- Modify: `src/importScripts.jl`
- Modify: `src/Routines/SearchDIA/SearchMethods/SearchMethods.jl`
- Modify: `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl`

- [ ] **Step 1: Load legacy ML files before scoring callers**

In `src/importScripts.jl`, expand the `src/utils/ML` include list so it includes:

```julia
[
    "fdrUtilities.jl",
    "ftrUtilities.jl",
    "probitRegression.jl",
    "piecewiseLinearFunction.jl",
    "spectralLinearRegression.jl",
    "wittakerHendersonSmoothing.jl",
    "scoring_config.jl",
    "pairing.jl",
    "training_selection.jl",
    "mbr_update.jl",
    "scoring_workspace.jl",
    "percolator_generic.jl",
    "percolatorSortOf.jl",
]
```

Keep `src/utils/ML/PSMScoring/*` loaded only if current retained code still needs it after PrecursorScoringSearch is reconciled.

- [ ] **Step 2: Load `selectTransitions.jl` before the common utility directory**

In `src/importScripts.jl`, add:

```julia
safe_include!(joinpath(package_root, "src", "Routines", "SearchDIA", "CommonSearchUtils", "selectTransitions", "selectTransitions.jl"))
```

Place it before:

```julia
safe_include_directory!(joinpath(package_root, "src", "Routines", "SearchDIA", "CommonSearchUtils"))
```

- [ ] **Step 3: Load staged search directories in dependency order**

In `src/importScripts.jl`, add explicit include blocks for:

```julia
include_files!(joinpath(search_methods_dir, "FirstPassSearch"), [
    "getBestPrecursorsAccrossRuns.jl",
    "utils.jl",
    "FirstPassSearch.jl",
])

include_files!(joinpath(search_methods_dir, "HuberTuningSearch"), [
    "utils.jl",
    "HuberTuningSearch.jl",
])

include_files!(joinpath(search_methods_dir, "SecondPassSearch"), [
    "utils.jl",
    "SecondPassSearch.jl",
])
```

- [ ] **Step 4: Stop explicitly loading fused execution**

Remove these explicit includes from `src/importScripts.jl`:

```julia
include_files!(joinpath(search_methods_dir, "MainSearch"), ...)
safe_include!(joinpath(package_root, "src", "Routines", "SearchDIA", "process_scans.jl"))
safe_include!(joinpath(package_root, "src", "Routines", "SearchDIA", "process_scans_fused.jl"))
```

Also exclude `MainSearch` in the catch-all `walkdir(search_methods_dir)` loop so the unused fused search types are not loaded accidentally. Do not exclude `BitVecCalibration`; it should continue to load and execute.

- [ ] **Step 5: Compile-check module loading**

Run:

```bash
julia --startup-file=no --project=. -e 'using Pioneer; println("loaded")'
```

Expected: this may fail with type/interface errors in `SearchTypes.jl`; those errors are expected before Task 5.

- [ ] **Step 6: Commit include rewiring**

Run:

```bash
git add src/importScripts.jl
git commit -m "refactor: load staged SearchDIA modules"
```

## Task 5: Reconcile Search Types And Search Context

**Files:**
- Modify: `src/Routines/SearchDIA/SearchMethods/SearchTypes.jl`
- Modify: `src/Routines/SearchDIA/SearchMethods/SearchMethods.jl`
- Modify if required: `src/structs/Counter.jl`
- Modify if required: `src/structs/SpectralLibrary/*`

- [ ] **Step 1: Restore staged path fields on `ArrowTableReference`**

Ensure `ArrowTableReference` has fields and accessors for:

```julia
first_pass_psms::Vector{String}
second_pass_psms::Vector{String}
passing_psms::Vector{String}
passing_proteins::Vector{String}
rt_index_paths::Vector{String}
```

Ensure these functions exist:

```julia
getFirstPassPsms(ref::ArrowTableReference, index::Int)
getSecondPassPsms(ref::ArrowTableReference, index::Int)
getPassingPsms(ref::ArrowTableReference, index::Int)
getPassingProteins(ref::ArrowTableReference, index::Int)
setFirstPassPsms!(ref::ArrowTableReference, index::Int, value::String)
setSecondPassPsms!(ref::ArrowTableReference, index::Int, value::String)
setPassingPsms!(ref::ArrowTableReference, index::Int, value::String)
setPassingProteins!(ref::ArrowTableReference, index::Int, value::String)
```

Remove `main_search_psms`, `fragment_index_matches`, and `filtered_fragment_matches` fields unless another retained non-fused stage still reads them.

- [ ] **Step 2: Restore staged `SimpleLibrarySearch` working arrays**

Ensure `SimpleLibrarySearch` supports the legacy first/second pass functions:

```julia
ion_matches::Vector{FragmentMatch{Float32}}
ion_misses::Vector{FragmentMatch{Float32}}
mass_err_matches::Vector{FragmentMatch{Float32}}
id_to_col
prec_count
ion_templates
scored_psms
unscored_psms
spectral_scores
complex_scored_psms
complex_unscored_psms
complex_spectral_scores
ms1_scored_psms
ms1_unscored_psms
ms1_spectral_scores
Hs
prec_ids
precursor_weights
temp_weights
residuals
isotopes
precursor_transmission
```

Use current concrete map/index types only when they satisfy the old API. If a current type does not support the required old methods, restore the old container type or write a narrow compatibility method with the old method name.

- [ ] **Step 3: Restore Huber state on `SearchContext`**

Ensure `SearchContext` has:

```julia
huber_delta::Base.Ref{Float32}
deconvolution_stop_tolerance::Base.Ref{Float32}
```

Ensure these accessors exist:

```julia
getHuberDelta(s::SearchContext)
setHuberDelta!(s::SearchContext, value::Float32)
getDeconvolutionStopTolerance(s::SearchContext)
setDeconvolutionStopTolerance!(s::SearchContext, value::Float32)
```

- [ ] **Step 4: Restore failed-file helpers only if staged code calls them**

If restored staged files call `check_and_skip_failed_file`, `markFileFailed!`, or `is_file_failed`, keep the current failed-file helpers from current `develop`.

If staged code does not call them after reconciliation, do not reintroduce extra failure-state behavior.

- [ ] **Step 5: Compile-check type reconciliation**

Run:

```bash
julia --startup-file=no --project=. -e 'using Pioneer; @assert isdefined(Pioneer, :FirstPassSearch); @assert isdefined(Pioneer, :SecondPassSearch); @assert isdefined(Pioneer, :HuberTuningSearch); println("staged types loaded")'
```

Expected:

```text
staged types loaded
```

- [ ] **Step 6: Commit type/context reconciliation**

Run:

```bash
git add src/Routines/SearchDIA/SearchMethods/SearchTypes.jl src/Routines/SearchDIA/SearchMethods/SearchMethods.jl src/structs
git commit -m "refactor: reconcile search context for staged SearchDIA"
```

## Task 6: Make The Staged Pipeline Canonical

**Files:**
- Modify: `src/Routines/SearchDIA.jl`
- Test: `test/UnitTests/test_searchdia_pipeline_order.jl`

- [ ] **Step 1: Add a testable canonical pipeline function**

Add this function near the top-level SearchDIA helpers in `src/Routines/SearchDIA.jl`:

```julia
function canonical_search_pipeline()
    return [
        ("Parameter Tuning", ParameterTuningSearch()),
        ("Quadrupole Tuning", QuadTuningSearch()),
        ("BitVec Calibration", BitVecCalibrationSearch()),
        ("First Pass Search", FirstPassSearch()),
        ("Huber Tuning", HuberTuningSearch()),
        ("Second Pass Search", SecondPassSearch()),
        ("Precursor Scoring", PrecursorScoringSearch()),
        ("Chromatogram Integration", IntegrateChromatogramSearch()),
        ("Protein Inference", ProteinInferenceSearch()),
        ("Protein Scoring", ProteinScoringSearch()),
        ("Quantification & Output", MaxLFQSearch()),
    ]
end
```

- [ ] **Step 2: Use the canonical function in `SearchDIA`**

Replace the local `searches = [...]` block with:

```julia
searches = canonical_search_pipeline()
```

- [ ] **Step 3: Update performance-report order**

In `print_performance_report`, replace the execution-order list so it contains:

```julia
execution_order = [
    "Parameter Loading", "Spectral Library Loading", "Search Context Initialization",
    "Parameter Tuning", "Quadrupole Tuning", "BitVec Calibration",
    "First Pass Search", "Huber Tuning", "Second Pass Search",
    "Precursor Scoring", "Chromatogram Integration",
    "Protein Inference", "Protein Scoring", "Quantification & Output",
]
```

- [ ] **Step 4: Run the pipeline-order test**

Run:

```bash
julia --startup-file=no --project=. test/UnitTests/test_searchdia_pipeline_order.jl
```

Expected: pass.

- [ ] **Step 5: Commit canonical pipeline rewiring**

Run:

```bash
git add src/Routines/SearchDIA.jl test/UnitTests/test_searchdia_pipeline_order.jl
git commit -m "refactor: restore staged SearchDIA pipeline order"
```

## Task 7: Preserve Current Parameter/Quad Tuning And Reconcile Huber Tuning

**Files:**
- Modify only if compatibility requires it: `src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/*`
- Modify only if compatibility requires it: `src/Routines/SearchDIA/SearchMethods/QuadTuningSearch/*`
- Modify: `src/Routines/SearchDIA/SearchMethods/HuberTuningSearch/*`
- Test: existing parameter/quad tuning unit tests

- [ ] **Step 1: Compile-check tuning constructors**

Run:

```bash
julia --startup-file=no --project=. -e 'using Pioneer; p = Pioneer.GetSearchParams(); println("params ok")'
```

Expected: params generation succeeds. If `GetSearchParams` is not callable in this form, run:

```bash
julia --startup-file=no --project=. -e 'using Pioneer; println("Pioneer loaded for tuning constructor check")'
```

- [ ] **Step 2: Preserve current parameter and quad tuning**

Do not restore the old `ParameterTuningSearch`, old `NceTuningSearch`, or old `QuadTuningSearch` from `a33f964`.

Keep these current files as the tuning source of truth:

```text
src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/constants.jl
src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/types.jl
src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/ParameterTuningSearch.jl
src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/utils.jl
src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/fit_intensity_mass_error.jl
src/Routines/SearchDIA/SearchMethods/QuadTuningSearch/QuadTuningSearch.jl
src/Routines/SearchDIA/SearchMethods/QuadTuningSearch/utils.jl
```

If staged first/second pass code needs an `nce_model`, use the current tuning output/accessor that populates `SearchContext.nce_model`; do not add the old `NceTuningSearch` stage.

- [ ] **Step 3: Restore Huber deconvolution dispatch**

In `HuberTuningSearch` and `SecondPassSearch`, preserve calls to the old Huber-loss deconvolution utilities from `a33f964`. Keep current `solvePoissonMM`/`solveOLS` only if they are unused by canonical staged search or if they have exact equivalent behavior for old Huber calls.

- [ ] **Step 4: Run focused tuning-related tests**

Run:

```bash
julia --startup-file=no --project=. test/UnitTests/RazoQuadModel.jl
julia --startup-file=no --project=. test/Routines/SearchDIA/SearchMethods/ParameterTuningSearch/test_intensity_spline_lipschitz.jl
```

Expected: pass, or the second command is skipped if that test file remains tied to the fused `fit_intensity_mass_error.jl` split and has not been ported.

- [ ] **Step 5: Commit tuning reconciliation**

Run:

```bash
git add src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch src/Routines/SearchDIA/SearchMethods/QuadTuningSearch src/Routines/SearchDIA/SearchMethods/HuberTuningSearch
git commit -m "fix: reconcile Huber tuning with current tuning stack"
```

## Task 8: Reconcile First And Second Pass Search With Current Library Structures

**Files:**
- Modify: `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/*`
- Modify: `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/*`
- Modify: `src/Routines/SearchDIA/CommonSearchUtils/*`
- Modify if required: `src/Routines/SearchDIA/ParseInputs/loadSpectralLibrary.jl`
- Modify if required: `src/structs/SpectralLibrary/*`

- [ ] **Step 1: Compile-check first/second pass methods**

Run:

```bash
julia --startup-file=no --project=. -e 'using Pioneer; @assert methods(Pioneer.process_file!).ms !== nothing; println("process_file methods loaded")'
```

Expected:

```text
process_file methods loaded
```

- [ ] **Step 2: Validate old query helpers against current library accessors**

Run:

```bash
rg -n "getFragments|getPrecursors|getProteins|LibraryFragmentIndex|DetailedFrag|SplineCompactFrag|PartitionedFragmentIndex" src/Routines/SearchDIA/CommonSearchUtils src/Routines/SearchDIA/SearchMethods/FirstPassSearch src/Routines/SearchDIA/SearchMethods/SecondPassSearch src/structs/SpectralLibrary
```

Expected: every old helper call has either an existing current accessor or a compatibility method.

- [ ] **Step 3: Add compatibility methods instead of reverting library prediction**

If old code expects `DetailedFrag`/`LibraryFragmentIndex` but current libraries expose `SpectralLibrary`/`PartitionedFragmentIndex`, add narrow compatibility accessors in the library/search boundary. Keep these methods read-only and scoped to search execution.

Do not change BuildSpecLib prediction or committed `.poin` fixture schema in this task.

- [ ] **Step 4: Run a minimal SearchDIA smoke test**

Run:

```bash
julia --startup-file=no --project=. -e 'using Pioneer; Pioneer.SearchDIA("/private/tmp/search.json")'
```

Expected: SearchDIA reaches at least the end of `Second Pass Search` in the compatible local dataset. If this fails before first-pass query execution due library shape, fix only the search/library compatibility layer. A later-stage failure is acceptable in this task and should be handled by Tasks 9-10.

- [ ] **Step 5: Commit first/second pass reconciliation**

Run:

```bash
git add src/Routines/SearchDIA/SearchMethods/FirstPassSearch src/Routines/SearchDIA/SearchMethods/SecondPassSearch src/Routines/SearchDIA/CommonSearchUtils src/Routines/SearchDIA/ParseInputs/loadSpectralLibrary.jl src/structs/SpectralLibrary
git commit -m "fix: restore first and second pass search"
```

## Task 9: Restore Probit And Percolator-Style Scoring

**Files:**
- Modify: `src/utils/ML/*`
- Modify: `src/Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/*`
- Test: `test/Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/test_scoring_interface.jl`
- Test: `test/UnitTests/test_psm_container.jl` if retained

- [ ] **Step 1: Restore old PrecursorScoringSearch behavior**

Use the `a33f964` version of `PrecursorScoringSearch` as the reference for second-pass PSM input, q-value/PEP computation, MBR feature handling, and output path registration.

Avoid comments or helper names that describe MainSearch as the producer of scored PSMs.

- [ ] **Step 2: Ensure `sort_of_percolator!` is the canonical PSM scorer where old staged code called it**

Search:

```bash
rg -n "sort_of_percolator!|percolator_scoring!|PSMScoring|ScoringConfig" src/Routines/SearchDIA/SearchMethods src/utils/ML
```

Expected: staged first/second/precursor scoring either calls `sort_of_percolator!` directly or calls the same old helper path from `a33f964`.

- [ ] **Step 3: Keep probit regression available**

Confirm:

```bash
julia --startup-file=no --project=. -e 'using Pioneer; @assert isdefined(Pioneer, :ProbitRegression); println("probit available")'
```

Expected:

```text
probit available
```

If the concrete exported name is different, use `rg -n "struct .*Probit|function .*probit|probitRegression" src/utils/ML src/Routines/SearchDIA` and check the symbol that old staged code actually calls.

- [ ] **Step 4: Run scoring tests**

Run:

```bash
julia --startup-file=no --project=. test/Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/test_scoring_interface.jl
```

Expected: pass after reconciliation.

- [ ] **Step 5: Commit scoring restoration**

Run:

```bash
git add src/utils/ML src/Routines/SearchDIA/SearchMethods/PrecursorScoringSearch test/Routines/SearchDIA/SearchMethods/PrecursorScoringSearch
git commit -m "fix: restore staged precursor scoring"
```

## Task 10: Reconcile Chromatogram Integration And Downstream Protein Stages

**Files:**
- Modify: `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/*`
- Modify if required: `src/Routines/SearchDIA/SearchMethods/ProteinInferenceSearch/*`
- Modify if required: `src/Routines/SearchDIA/SearchMethods/ProteinScoringSearch/*`
- Modify if required: `src/Routines/SearchDIA/SearchMethods/MaxLFQSearch/*`

- [ ] **Step 1: Ensure integration consumes staged passing PSMs**

Search:

```bash
rg -n "getPassingPsms|getSecondPassPsms|getMainSearchPsms|getFilteredFragmentMatches|getFragmentIndexMatches" src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch src/Routines/SearchDIA/SearchMethods/ProteinInferenceSearch src/Routines/SearchDIA/SearchMethods/ProteinScoringSearch src/Routines/SearchDIA/SearchMethods/MaxLFQSearch
```

Expected: downstream stages consume `passing_psms`, `second_pass_psms`, or protein-stage outputs, not `main_search_psms` or fused-only fragment match files.

- [ ] **Step 2: Preserve current decoy-through-protein-scoring behavior if compatible**

Keep the current post-#377 decoy behavior only if the staged integration output includes the target/decoy columns expected by protein scoring.

If staged integration cannot support that behavior cleanly, revert only the decoy-through-protein-scoring patch and document that it needs a separate reintroduction PR after staged search is stable.

- [ ] **Step 3: Run downstream unit tests**

Run:

```bash
julia --startup-file=no --project=. test/UnitTests/test_global_protein_inference.jl
julia --startup-file=no --project=. test/UnitTests/test_maxLFQ.jl
julia --startup-file=no --project=. test/UnitTests/test_normalizeQuant.jl
```

Expected: pass, or only fail because a test depends on fused-specific PSM fixture names; fix those fixtures to staged names.

- [ ] **Step 4: Commit downstream reconciliation**

Run:

```bash
git add src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch src/Routines/SearchDIA/SearchMethods/ProteinInferenceSearch src/Routines/SearchDIA/SearchMethods/ProteinScoringSearch src/Routines/SearchDIA/SearchMethods/MaxLFQSearch src/utils/maxLFQ.jl src/utils/proteinInference.jl test/UnitTests
git commit -m "fix: align downstream stages with staged SearchDIA"
```

## Task 11: Remove Fused-Only MainSearch Canonical Code

**Files:**
- Delete or stop tracking: `src/Routines/SearchDIA/SearchMethods/MainSearch/*`
- Delete or stop tracking: `src/Routines/SearchDIA/process_scans.jl`
- Delete or stop tracking: `src/Routines/SearchDIA/process_scans_fused.jl`
- Modify: `src/importScripts.jl`
- Modify: `test/runtests_part3_units.jl`

- [ ] **Step 1: Remove fused-only tests from the default suite**

Remove these includes from `test/runtests_part3_units.jl`:

```julia
include("./UnitTests/test_fused_prec_filters.jl")
include("./UnitTests/test_run_fused.jl")
include("./UnitTests/test_feature_finiteness.jl")
```

Keep the files on disk only if they are useful as branch-local diagnostics; do not run them in canonical `develop` tests.

- [ ] **Step 2: Delete unused fused execution files**

Run:

```bash
git rm -r src/Routines/SearchDIA/SearchMethods/MainSearch
git rm src/Routines/SearchDIA/process_scans.jl src/Routines/SearchDIA/process_scans_fused.jl
```

Expected: files are staged for deletion. If another retained stage still imports a helper from these files, move that helper to a neutral staged-search location before deletion.

- [ ] **Step 3: Confirm no canonical references remain**

Run:

```bash
rg -n "MainSearch|process_scans_fused|run_fused" src test/runtests*.jl
rg -n "fusedScan|fusedMatch" src/Routines/SearchDIA/SearchMethods/ParameterTuningSearch src/Routines/SearchDIA/SearchMethods/QuadTuningSearch src/Routines/SearchDIA/CommonSearchUtils
```

Expected: no `MainSearch`, `process_scans_fused`, or `run_fused` references in canonical source or default tests. `fusedMatch` and `fusedScan` references may remain only in current tuning code or neutral common utilities they require.

- [ ] **Step 4: Commit fused cleanup**

Run:

```bash
git add -u src test/runtests_part3_units.jl
git commit -m "refactor: remove fused SearchDIA from canonical pipeline"
```

## Task 12: Run Local Unit Verification

**Files:**
- Execute only

- [ ] **Step 1: Load the package**

Run:

```bash
julia --startup-file=no --project=. -e 'using Pioneer; println("loaded")'
```

Expected:

```text
loaded
```

- [ ] **Step 2: Run targeted unit tests**

Run:

```bash
julia --startup-file=no --project=. test/UnitTests/test_searchdia_pipeline_order.jl
julia --startup-file=no --project=. test/Routines/SearchDIA/SearchMethods/PrecursorScoringSearch/test_scoring_interface.jl
julia --startup-file=no --project=. test/UnitTests/test_maxLFQ.jl
```

Expected: all pass.

- [ ] **Step 3: Run the full local test driver**

Run:

```bash
julia --startup-file=no --project=. test/runtests.jl
```

Expected: all test parts pass.

- [ ] **Step 4: Do not use incompatible regression fixtures as the acceptance gate**

Do not run the dedicated GitHub regression workflow fixtures for this restoration branch. Their predicted libraries are not compatible with this staged-search restoration target. The end-to-end acceptance gate is the local compatible SearchDIA run in Task 13.

## Task 13: Run Local Compatible SearchDIA Acceptance

**Files:**
- Execute only

- [ ] **Step 1: Confirm the local params file exists**

Run:

```bash
test -f /private/tmp/search.json
```

Expected: exit code `0`. If this fails, stop and ask for the updated local params path.

- [ ] **Step 2: Run the local end-to-end search**

Run:

```bash
julia --startup-file=no --project=. -e 'using Pioneer; Pioneer.SearchDIA("/private/tmp/search.json")'
```

Expected: exits `0`. The search log should show the staged canonical step order: `Parameter Tuning`, `Quadrupole Tuning`, `BitVec Calibration`, `First Pass Search`, `Huber Tuning`, `Second Pass Search`, `Precursor Scoring`, `Chromatogram Integration`, `Protein Inference`, `Protein Scoring`, and `Quantification & Output`.

- [ ] **Step 3: Confirm outputs were produced**

Run:

```bash
julia --startup-file=no --project=. -e 'using Pioneer; params = Pioneer.parse_pioneer_parameters("/private/tmp/search.json"); results = params.paths[:results]; required = ("precursors_long.arrow", "precursors_wide.arrow", "protein_groups_long.arrow", "protein_groups_wide.arrow", "pioneer_search_report.txt", "pioneer_search_log.log", "pioneer_warnings.log"); missing = [name for name in required if !isfile(joinpath(results, name))]; println("results=", results); isempty(missing) || error("Missing outputs: " * join(missing, ", ")); println("all required outputs present")'
```

Expected: prints the results directory and `all required outputs present`.

- [ ] **Step 4: Commit the accepted local verification notes if any test fixtures changed**

Run:

```bash
git status --short
```

Expected: no untracked or modified result files under the repository. If local verification required fixture updates, commit only those fixture/test changes with:

```bash
git add test
git commit -m "test: update staged SearchDIA local fixtures"
```

## Task 14: Review, PR, And Merge

**Files:**
- Execute only

- [ ] **Step 1: Confirm no library-prediction files changed unintentionally**

Run:

```bash
git diff --name-only origin/develop...HEAD -- src/Routines/BuildSpecLib data/test_build_spec_lib assets/example_config/defaultBuildLibParams.json
```

Expected: no output, unless a staged-search compatibility issue explicitly required a library reader-only change. If output appears for BuildSpecLib prediction/build files, split those changes out.

- [ ] **Step 2: Open PR**

Run:

```bash
gh pr create --base develop --head hotfix/restore-staged-searchdia-pipeline --title "hotfix: restore staged SearchDIA pipeline" --body "Restores the staged first-pass/probit/Huber/second-pass SearchDIA pipeline as the canonical develop behavior while preserving current parameter/NCE-related tuning, current quad tuning, and BitVecCalibration. Removes the fused MainSearch execution path from develop; fused search remains available in branch history. Keeps library prediction/schema changes out of scope. Acceptance uses local compatible SearchDIA params instead of the incompatible dedicated regression library fixtures."
```

- [ ] **Step 3: Merge after local tests and the local SearchDIA acceptance run are reviewed**

Use the repository's normal merge policy.

## Compatibility Rules

- Prefer restoring old search behavior over preserving fused abstractions.
- Prefer narrow compatibility methods over BuildSpecLib/library prediction changes.
- Do not keep both `MainSearch` and staged search as user-facing runtime paths.
- Do not carry forward comments that describe MainSearch as the canonical producer.
- Keep current protein/scoring/output fixes only when they consume staged outputs cleanly.
- If a change requires re-predicting libraries, stop and split it into a separate final-stage plan.
- Do not use the dedicated regression workflow as the acceptance gate for this branch because those predicted libraries are not compatible with the restored staged-search path.

## Self-Review

- The plan starts from current `develop`, not from `a33f964`.
- The plan restores staged SearchDIA as canonical behavior without a runtime switch.
- The plan keeps fused search available only through branch checkout/history.
- The plan keeps current parameter/NCE-related tuning, current quad tuning, and BitVec calibration.
- The plan restores first-pass, probit/percolator scoring, Huber tuning, and second-pass search as one coherent pipeline.
- The plan avoids library prediction/schema changes unless staged search compatibility requires a read-only adapter.
- The plan includes failing tests first, focused compile checks, local unit tests, and a local compatible SearchDIA acceptance run.
