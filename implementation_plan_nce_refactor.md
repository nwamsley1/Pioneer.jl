# Implementation Plan: NCE Model Refactor (Solution 1)

## Objective

Remove the mutable `nce_model::Base.Ref` from `SplineFragmentLookup` and pass NCE model as an explicit parameter throughout the call chain. This eliminates the root cause of the GC crash (Error 2) by removing mutable shared state.

## Benefits

1. **Thread-Safe by Design**: No mutable shared state
2. **Explicit Dependencies**: NCE model dependency visible in function signatures
3. **Better Semantics**: NCE is context, not library property
4. **No GC Issues**: No mutable Ref for GC to corrupt
5. **Testability**: Easy to test with different NCE models
6. **Performance**: No Ref dereferencing overhead

## Implementation Steps

### Phase 1: Core Library Changes

#### Step 1.1: Remove nce_model from SplineFragmentLookup

**File:** `src/structs/LibraryIon.jl`

**Current (lines ~550-556):**
```julia
struct SplineFragmentLookup{N,M,T<:AbstractFloat} <: LibraryFragmentLookup
    frags::Vector{SplineDetailedFrag{N,T}}
    prec_frag_ranges::Vector{UInt64}
    knots::NTuple{M, T}
    nce_model::Base.Ref{<:NceModel{T}}  # ← REMOVE THIS LINE
    degree::Int64
end
```

**New:**
```julia
struct SplineFragmentLookup{N,M,T<:AbstractFloat} <: LibraryFragmentLookup
    frags::Vector{SplineDetailedFrag{N,T}}
    prec_frag_ranges::Vector{UInt64}
    knots::NTuple{M, T}
    # nce_model field removed!
    degree::Int64
end
```

**Constructors to Update:**
- Update any constructors that initialize `nce_model` field
- Remove `nce_model` parameter from constructor calls

#### Step 1.2: Remove setNceModel! functions

**File:** `src/structs/LibraryIon.jl`

**Lines to Remove (~565-573):**
```julia
# DELETE these functions entirely:
function setNceModel!(lookup::SplineFragmentLookup{N,M,T}, new_nce_model::NceModel{T}) where {N,M,T<:AbstractFloat}
    lookup.nce_model[] = new_nce_model
end

function setNceModel!(lookup::StandardFragmentLookup, new_nce_model::NceModel{T}) where {T<:AbstractFloat}
    return nothing
end
```

#### Step 1.3: Update getNCE to accept model parameter

**File:** `src/structs/LibraryIon.jl`

**Current (lines ~591-596):**
```julia
function getNCE(lfp::SplineFragmentLookup, prec_charge::UInt8, prec_mz::T) where {T<:AbstractFloat}
    return lfp.nce_model[](prec_mz, prec_charge)
end

function getNCE(lfp::SplineFragmentLookup)
    return lfp.nce_model[]()
end
```

**New:**
```julia
# Delete old getNCE functions that use nce_model field
# These are no longer needed since we pass NCE explicitly
```

#### Step 1.4: Update getSplineData to accept NCE parameter

**File:** `src/structs/LibraryIon.jl`

**Current (lines ~575-589):**
```julia
function getSplineData(lfp::SplineFragmentLookup{N,M,T}, prec_charge::UInt8, prec_mz::T) where {N,M,T<:AbstractFloat}
    return SplineType(
        getKnots(lfp),
        getNCE(lfp, prec_charge, prec_mz),
        getDegree(lfp)
    )
end

function getSplineData(lfp::SplineFragmentLookup{N,M,T}) where {N,M,T<:AbstractFloat}
    return SplineType(
        getKnots(lfp),
        getNCE(lfp),
        getDegree(lfp)
    )
end
```

**New:**
```julia
# Version with precursor-specific NCE
function getSplineData(
    lfp::SplineFragmentLookup{N,M,T},
    nce_model::NceModel{T},
    prec_charge::UInt8,
    prec_mz::T
) where {N,M,T<:AbstractFloat}
    return SplineType(
        getKnots(lfp),
        nce_model(prec_mz, prec_charge),  # Direct call to model
        getDegree(lfp)
    )
end

# Version with default NCE (for cases without specific precursor)
function getSplineData(
    lfp::SplineFragmentLookup{N,M,T},
    nce_model::NceModel{T}
) where {N,M,T<:AbstractFloat}
    return SplineType(
        getKnots(lfp),
        nce_model(),  # Call with no arguments
        getDegree(lfp)
    )
end

# Keep StandardFragmentLookup versions unchanged (no NCE)
function getSplineData(lfp::StandardFragmentLookup, nce_model::NceModel{T}, prec_charge::UInt8, prec_mz::T) where {T<:AbstractFloat}
    return ConstantType()
end

function getSplineData(lfp::StandardFragmentLookup, nce_model::NceModel{T}) where {T<:AbstractFloat}
    return ConstantType()
end
```

### Phase 2: Update Transition Selection Functions

#### Step 2.1: Update fillTransitionListPrecomputed!

**File:** `src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/fillTransitionList.jl`

**Current signature (line ~95):**
```julia
function fillTransitionListPrecomputed!(
    transitions::Vector{DetailedFrag{Float32}},
    prec_estimation_type::PrecEstimation,
    precursor_fragment_range::UnitRange{UInt64},
    fragment_ions::Vector{F},
    spline_data::G,  # ← This is SplineType created from getSplineData
    prec_mz::Float32,
    prec_charge::UInt8,
    prec_sulfur_count::UInt8,
    transition_idx::Int64,
    precursor_transmission::Vector{Float32},
    prec_isotope_set::Tuple{Int64, Int64},
    isotopes::Vector{Float32},
    n_frag_isotopes::Int64,
    max_frag_rank::UInt8,
    iso_splines::IsotopeSplineModel,
    frag_mz_bounds::Tuple{Float32, Float32},
    block_size::Int64
) where {G<:IntensityDataType, F <: AltimeterFragment}
```

**Changes:**
- `spline_data` parameter already contains the NCE value (from SplineType)
- No signature change needed here!
- The work is in the callers that create spline_data

#### Step 2.2: Update rtIndexTransitionSelection

**File:** `src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/rtIndexTransitionSelection.jl`

**Current (lines 29-53):**
```julia
function _select_transitions_impl!(
    transitions::Vector{DetailedFrag{Float32}},
    ::RTIndexedTransitionSelection,
    prec_estimation_type::PrecEstimation,
    transition_idx::Int64,
    lookup::LibraryFragmentLookup,
    precs_temp::Vector{UInt32},
    prec_mzs::AbstractArray{Float32},
    prec_charges::AbstractArray{UInt8},
    prec_sulfur_counts::AbstractArray{UInt8},
    iso_splines::IsotopeSplineModel,
    quad_transmission_func::QuadTransmissionFunction,
    precursor_transmission::Vector{Float32},
    isotopes::Vector{Float32},
    n_frag_isotopes::Int64,
    max_frag_rank::UInt8,
    rt_index::retentionTimeIndex{Float32, Float32},
    rt_start_idx::Int64,
    rt_stop_idx::Int64,
    frag_mz_bounds::Tuple{Float32, Float32};
    precursors_passing::Union{Set{UInt32}, Nothing} = nothing,
    isotope_err_bounds::Tuple{I, I} = (3, 1),
    block_size::Int64 = 10000,
    min_fraction_transmitted::Float32 = 0.0f0
) where {I<:Integer}
```

**New signature:**
```julia
function _select_transitions_impl!(
    transitions::Vector{DetailedFrag{Float32}},
    ::RTIndexedTransitionSelection,
    prec_estimation_type::PrecEstimation,
    transition_idx::Int64,
    lookup::LibraryFragmentLookup,
    nce_model::NceModel{Float32},  # ← NEW PARAMETER
    precs_temp::Vector{UInt32},
    prec_mzs::AbstractArray{Float32},
    prec_charges::AbstractArray{UInt8},
    prec_sulfur_counts::AbstractArray{UInt8},
    iso_splines::IsotopeSplineModel,
    quad_transmission_func::QuadTransmissionFunction,
    precursor_transmission::Vector{Float32},
    isotopes::Vector{Float32},
    n_frag_isotopes::Int64,
    max_frag_rank::UInt8,
    rt_index::retentionTimeIndex{Float32, Float32},
    rt_start_idx::Int64,
    rt_stop_idx::Int64,
    frag_mz_bounds::Tuple{Float32, Float32};
    precursors_passing::Union{Set{UInt32}, Nothing} = nothing,
    isotope_err_bounds::Tuple{I, I} = (3, 1),
    block_size::Int64 = 10000,
    min_fraction_transmitted::Float32 = 0.0f0
) where {I<:Integer}
```

**Update call to getSplineData (line ~121):**
```julia
# OLD:
transition_idx = @inline fillTransitionListPrecomputed!(
    transitions,
    prec_estimation_type,
    getPrecFragRange(lookup, prec_idx),
    getFragments(lookup),
    getSplineData(lookup, prec_charge, prec_mz),  # ← OLD
    prec_mz,
    prec_charge,
    prec_sulfur_count,
    transition_idx,
    precursor_transmission,
    prec_isotope_set,
    isotopes,
    n_frag_isotopes,
    max_frag_rank,
    iso_splines,
    frag_mz_bounds,
    block_size,
)

# NEW:
transition_idx = @inline fillTransitionListPrecomputed!(
    transitions,
    prec_estimation_type,
    getPrecFragRange(lookup, prec_idx),
    getFragments(lookup),
    getSplineData(lookup, nce_model, prec_charge, prec_mz),  # ← Pass nce_model
    prec_mz,
    prec_charge,
    prec_sulfur_count,
    transition_idx,
    precursor_transmission,
    prec_isotope_set,
    isotopes,
    n_frag_isotopes,
    max_frag_rank,
    iso_splines,
    frag_mz_bounds,
    block_size,
)
```

#### Step 2.3: Update standardTransitionSelection

**File:** `src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/standardTransitionSelection.jl`

**Add nce_model parameter to signature** (similar to Step 2.2)

**Update getSplineData call** (line ~69):
```julia
# OLD:
getSplineData(lookup, prec_charge, prec_mz),

# NEW:
getSplineData(lookup, nce_model, prec_charge, prec_mz),
```

#### Step 2.4: Update quadEstimationSelection

**File:** `src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/quadEstimationSelection.jl`

**Add nce_model parameter and update call** (line ~66)

#### Step 2.5: Update massErrEstimationStrategy

**File:** `src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/massErrEstimationStrategy.jl`

**Update getSplineData call** (line ~57):
```julia
# OLD:
transitions[transition_idx] = convert_to_detailed(frag, getSplineData(lookup))

# NEW:
transitions[transition_idx] = convert_to_detailed(frag, getSplineData(lookup, nce_model))
```

#### Step 2.6: Update selectTransitions! wrapper

**File:** `src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/selectTransitions.jl`

**Update signature** (line ~87):
```julia
# OLD:
function selectTransitions!(
    transitions::Vector{DetailedFrag{Float32}},
    strategy::TransitionSelectionStrategy,
    prec_estimation_type::PrecEstimation,
    lookup::LibraryFragmentLookup,
    precs_temp::Vector{UInt32},
    # ... other parameters
)

# NEW:
function selectTransitions!(
    transitions::Vector{DetailedFrag{Float32}},
    strategy::TransitionSelectionStrategy,
    prec_estimation_type::PrecEstimation,
    lookup::LibraryFragmentLookup,
    nce_model::NceModel,  # ← NEW PARAMETER (after lookup)
    precs_temp::Vector{UInt32},
    # ... other parameters
)
```

**Update call to _select_transitions_impl!:**
```julia
# Pass nce_model to strategy-specific implementation
transition_idx, precs_temp_size = _select_transitions_impl!(
    transitions,
    strategy,
    prec_estimation_type,
    0,  # transition_idx
    lookup,
    nce_model,  # ← NEW
    precs_temp,
    # ... rest
)
```

### Phase 3: Update Search Method Callers

#### Step 3.1: Update IntegrateChromatogramsSearch

**File:** `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/IntegrateChromatogramsSearch.jl`

**Remove setNceModel! call** (lines 188-191):
```julia
# DELETE these lines:
setNceModel!(
    getFragmentLookupTable(getSpecLib(search_context)),
    getNceModelModel(search_context, ms_file_idx)
)
```

**File:** `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/utils.jl`

**Update extract_chromatograms function:**

Add nce_model parameter and pass to build_chromatograms:
```julia
function extract_chromatograms(
    spectra::MassSpecData,
    passing_psms::DataFrame,
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64,
    chrom_type::ChromType
)
    # Get NCE model for this file
    nce_model = getNceModelModel(search_context, ms_file_idx)

    # ... existing setup code ...

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            search_data = getSearchData(search_context)[thread_id]

            return build_chromatograms(
                spectra,
                last(thread_task),
                precursor_sets[thread_id],
                rt_index,
                search_context,
                search_data,
                params,
                ms_file_idx,
                nce_model,  # ← NEW PARAMETER
                chrom_type
            )
        end
    end
    return vcat(fetch.(tasks)...)
end
```

**Update build_chromatograms signatures** (lines ~216 and ~430):
```julia
function build_chromatograms(
    spectra::MassSpecData,
    scan_range::Vector{Int64},
    precursors_passing::Set{UInt32},
    rt_index::retentionTimeIndex,
    search_context::SearchContext,
    search_data::SearchDataStructures,
    params::IntegrateChromatogramSearchParameters,
    ms_file_idx::Int64,
    nce_model::NceModel{Float32},  # ← NEW PARAMETER
    ::MS2CHROM
)
```

**Update selectTransitions! call** (line ~270):
```julia
ion_idx, prec_temp_size = selectTransitions!(
    getIonTemplates(search_data),
    RTIndexedTransitionSelection(),
    params.prec_estimation,
    getFragmentLookupTable(getSpecLib(search_context)),
    nce_model,  # ← NEW PARAMETER (after lookup)
    precs_temp,
    getMz(getPrecursors(getSpecLib(search_context))),
    getCharge(getPrecursors(getSpecLib(search_context))),
    getSulfurCount(getPrecursors(getSpecLib(search_context))),
    getIsoSplines(search_data),
    quad_func,
    getPrecursorTransmission(search_data),
    getIsotopes(search_data),
    params.n_frag_isotopes,
    params.max_frag_rank,
    rt_index,
    irt_start,
    irt_stop,
    (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx));
    precursors_passing = precursors_passing,
    isotope_err_bounds = params.isotope_err_bounds,
    block_size = 10000,
    min_fraction_transmitted = params.min_fraction_transmitted
)
```

#### Step 3.2: Update SecondPassSearch

**File:** `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/SecondPassSearch.jl`

**Remove setNceModel! call** (lines 311-314):
```julia
# DELETE:
setNceModel!(
    getFragmentLookupTable(getSpecLib(search_context)),
    getNceModelModel(search_context, ms_file_idx)
)
```

**Get NCE model and pass to perform_second_pass_search:**
```julia
# ADD after line 314:
nce_model = getNceModelModel(search_context, ms_file_idx)

# UPDATE perform_second_pass_search call to pass nce_model:
psms = perform_second_pass_search(
    spectra,
    rt_index,
    search_context,
    params,
    ms_file_idx,
    nce_model  # ← NEW PARAMETER
)
```

**Update perform_second_pass_search signature and selectTransitions! calls inside it**

#### Step 3.3: Update FirstPassSearch

**File:** `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl`

**Remove setNceModel! call** (lines 234-237):
```julia
# DELETE:
setNceModel!(
    getFragmentLookupTable(getSpecLib(search_context)),
    getNceModelModel(search_context, ms_file_idx)
)
```

**Pass nce_model to library_search:**
```julia
nce_model = getNceModelModel(search_context, ms_file_idx)
return library_search(spectra, search_context, params, ms_file_idx, nce_model)
```

**Update library_search and any selectTransitions! calls**

#### Step 3.4: Update LibrarySearch (NCE tuning)

**File:** `src/Routines/SearchDIA/LibrarySearch.jl`

**Update NCE grid search** (lines 385-398):
```julia
# OLD:
all_results = map(nce_grid) do nce
    # Update NCE model in fragment lookup table
    setNceModel!(
        getFragmentLookupTable(spec_lib),
        PiecewiseNceModel(nce)
    )

    # Run getPSMS with updated NCE model
    # ...
end

# NEW:
all_results = map(nce_grid) do nce
    # Create NCE model for this iteration
    nce_model = PiecewiseNceModel(nce)

    # Run getPSMS with explicit NCE model parameter
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            try
                psms = getPSMS(
                    getThreadSpecData(search_context, thread_id),
                    # ... other parameters ...
                    nce_model  # ← Pass as parameter
                )
            catch e
                # error handling
            end
        end
    end
    # ...
end
```

**Update getPSMS and related functions to accept nce_model parameter**

### Phase 4: Update Library Loading

#### Step 4.1: Update library construction

**Files that create SplineFragmentLookup:**
- Check library building code
- Remove nce_model initialization from constructors
- Update any code that expects nce_model field

**Likely files:**
- `src/Routines/BuildSpecLib/*.jl`
- Anywhere `SplineFragmentLookup` is constructed

### Phase 5: Testing and Validation

#### Test 1: Compilation
```bash
julia --project=. -e 'using Pkg; Pkg.build()'
```

Should compile without errors.

#### Test 2: Unit Tests
```bash
julia --project=. -e 'using Pkg; Pkg.test("Pioneer")'
```

All tests should pass.

#### Test 3: Integration Test - BuildSpecLib
```bash
julia --project=. -e 'using Pkg; Pkg.test("Pioneer"; test_args=["BuildSpecLib"])'
```

Should build library without nce_model field.

#### Test 4: Integration Test - Chromatogram Integration
Run the workload that previously crashed:
- 312 files
- IntegrateChromatograms step
- Monitor for EXCEPTION_ACCESS_VIOLATION

#### Test 5: Verify NCE Values
Add logging to verify correct NCE model is used for each file:
```julia
@info "Processing file $ms_file_idx with NCE model: $(nce_model)"
```

### Phase 6: Cleanup

#### Remove dead code:
- `setNceModel!` functions (already done in Step 1.2)
- `getNCE` functions that use nce_model field (already done in Step 1.3)
- Any other NCE-model-Ref-related code

#### Update documentation:
- Update docstrings to reflect NCE as parameter
- Update CLAUDE.md files with new pattern
- Note that NCE model is file-specific context, not library state

## File Checklist

### Files to Modify (Core):
- [ ] `src/structs/LibraryIon.jl` - Remove nce_model field, update getSplineData
- [ ] `src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/selectTransitions.jl` - Add nce_model parameter
- [ ] `src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/rtIndexTransitionSelection.jl` - Add nce_model parameter
- [ ] `src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/standardTransitionSelection.jl` - Add nce_model parameter
- [ ] `src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/quadEstimationSelection.jl` - Add nce_model parameter
- [ ] `src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/massErrEstimationStrategy.jl` - Add nce_model parameter

### Files to Modify (Search Methods):
- [ ] `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/IntegrateChromatogramsSearch.jl` - Remove setNceModel!, pass to utils
- [ ] `src/Routines/SearchDIA/SearchMethods/IntegrateChromatogramsSearch/utils.jl` - Add nce_model to functions, pass to selectTransitions!
- [ ] `src/Routines/SearchDIA/SearchMethods/SecondPassSearch/SecondPassSearch.jl` - Remove setNceModel!, pass nce_model
- [ ] `src/Routines/SearchDIA/SearchMethods/FirstPassSearch/FirstPassSearch.jl` - Remove setNceModel!, pass nce_model
- [ ] `src/Routines/SearchDIA/LibrarySearch.jl` - Update NCE grid search to pass nce_model

### Files to Check (Library Building):
- [ ] `src/Routines/BuildSpecLib/*.jl` - Update SplineFragmentLookup construction

### Files to Update (Documentation):
- [ ] `src/Routines/SearchDIA/CLAUDE.md`
- [ ] `src/Routines/SearchDIA/CommonSearchUtils/selectTransitions/CLAUDE.md`

## Estimated Effort

- **Core changes**: 2-3 hours
- **Search method updates**: 2-3 hours
- **Library building updates**: 1-2 hours
- **Testing**: 2-3 hours
- **Documentation**: 1 hour

**Total**: ~8-12 hours

## Risk Assessment

**Low Risk:**
- Type system will catch missing parameters at compile time
- No runtime behavior changes (same NCE values used)
- Each step is independently testable

**Potential Issues:**
- Missing a call site (will fail at compile time)
- Library building code needs nce_model removed from construction
- Performance overhead of parameter passing (negligible)

## Success Criteria

1. ✅ Code compiles without errors
2. ✅ All tests pass
3. ✅ No EXCEPTION_ACCESS_VIOLATION crashes
4. ✅ Correct NCE model used for each file (verified by logging)
5. ✅ No performance regression (<1% expected)
6. ✅ Clean separation: NCE is context, not library state

## Rollback Plan

If issues arise:
1. Git revert to previous commit
2. Use Solution 4 (GC.@preserve) as temporary workaround
3. Debug specific issues
4. Retry with fixes

## Post-Implementation

After successful implementation:
1. Monitor Windows runs for stability
2. Performance profiling
3. Consider additional refactoring opportunities
4. Document pattern for future file-specific context needs
