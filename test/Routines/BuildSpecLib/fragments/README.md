# Fragment Processing Unit Tests

This directory contains comprehensive unit tests for the fragment processing functionality in Pioneer.jl's spectral library building pipeline.

## Test Files

### test_fragment_annotation.jl
Tests for fragment ion annotation parsing and indexing:
- `get_ion_annotation_set` - Parsing regular, internal, and immonium ions
- `parse_internal_ion` - Extracting modification info from internal ions
- `create_ion_annotation_index` - Creating deterministic ion type indices
- `parse_fragment_annotation` - Parsing full fragment annotations with modifications
- `getAnnotationToID` - Bidirectional mapping between annotations and IDs
- `get_altimeter_ion_dict` - Loading Altimeter ion dictionaries

### test_fragment_parse.jl
Tests for fragment data parsing and processing:
- `parseInternalIon` - Internal ion string parsing
- `getIonAnnotationSet` - Creating unique ion type sets
- `getIonAnnotationDict` - Mapping annotations to unique IDs
- `countSulfurLoss` - Counting sulfurs in neutral losses/gains
- `count_sulfurs!` - Tracking sulfur content in sequences and modifications
- `fill_isotope_mods!` - Handling isotope-labeled modifications
- `get_fragment_indices` - Calculating sequence bounds for different ion types
- `apply_isotope_mod` - Applying isotope mass shifts
- `getNumeric` - Extracting numeric values from strings
- `get_immonium_sulfur_dict` - Loading immonium ion sulfur counts
- Integration test for `parse_koina_fragments`

### test_fragment_predict.jl
Tests for fragment intensity prediction via Koina API:
- Model type validation (InstrumentSpecific, InstrumentAgnostic, SplineCoefficient)
- `filter_fragments!` - Filtering by intensity and m/z thresholds
- `sort_fragments!` - Sorting by precursor and intensity
- `predict_fragments_batch` - Batch prediction for different model types
- `predict_fragments` - Main prediction dispatcher
- Batch processing calculations

### test_get_frag_bounds.jl
Tests for fragment m/z boundary detection:
- `get_fragment_bounds` from MS data - Linear model fitting
- Auto-detection from raw files with fallback to defaults
- `FragBoundModel` polynomial evaluation
- Edge cases and error handling

### test_fragments_suite.jl
Main test runner that executes all fragment tests in sequence.

## Running Tests

To run all fragment tests:
```julia
julia> include("test/Routines/BuildSpecLib/fragments/test_fragments_suite.jl")
```

To run individual test files:
```julia
julia> include("test/Routines/BuildSpecLib/fragments/test_fragment_annotation.jl")
```

## Test Coverage

The tests cover:
- **Data structures**: PioneerFrag, PioneerSplineFrag, PioneerFragAnnotation, FragBoundModel
- **Parsing logic**: Regular ions, internal ions, immonium ions, modifications
- **Chemical tracking**: Sulfur counting, isotope modifications
- **API integration**: Koina model types and batch processing
- **Edge cases**: Missing data, invalid formats, empty inputs
- **Performance**: Large dataset handling, batch size calculations

## Mock Data

Tests use mock data to avoid external dependencies:
- Mock Arrow tables for precursor and fragment data
- Temporary files for I/O testing
- Synthetic MS data for fragment bound detection

## Integration Points

These tests verify integration with:
- Koina API models (UniSpec, Prosit, Altimeter, AlphaPeptDeep)
- Arrow file I/O
- JLD2 serialization
- DataFrames processing
- Polynomial fitting for fragment bounds

## Known Limitations

- Actual Koina API calls are not tested (would require live server)
- Some complex regex patterns may have edge cases not covered
- Performance benchmarks not included (focus on correctness)