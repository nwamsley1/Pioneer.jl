# BuildSpecLib Integration Test Implementation Plan

## Overview
Create comprehensive integration tests for BuildSpecLib using test FASTA files with the Altimeter model for fragment prediction via Koina API. All test scenarios will include fragment prediction.

## Test Structure

### 1. Create new test file: `test/Routines/BuildSpecLib/test_build_spec_lib.jl`

### 2. Test Scenarios (All with Altimeter Fragment Prediction)

#### Scenario A: Missed Cleavage Verification with Minimal Protein
- **FASTA source**: minimal_protein.fasta (12 AA protein designed for easy verification)
- **Test Cases**:
  - A1: 0 missed cleavages → expect 0 peptides (all too short)
  - A2: 1 missed cleavage → expect 2 peptides (MARALYK, ALYKDER)
  - A3: 2 missed cleavages → expect 4 peptides (adds MARALYKDER, ALYKDERPK)
  - A4: 3 missed cleavages → expect 5 peptides (adds full sequence MARALYKDERPK)
- **Common Parameters**:
  - `predict_fragments: true`
  - `prediction_model: "altimeter"`
  - `min_length: 7`, `max_length: 15`
  - `min_charge: 2`, `max_charge: 2`
  - `add_decoys: false` (for easier counting)
  - `max_koina_batch: 10` (very small batch)
- **Verification**: Count peptides exactly matches expected for each missed cleavage setting

#### Scenario B: Standard Library Build with Altimeter
- **FASTA source**: Single small FASTA file (single_protein.fasta with 1 protein)
- **Parameters**:
  - `predict_fragments: true` 
  - `prediction_model: "altimeter"`
  - `min_length: 7`, `max_length: 12` (keep peptide count low)
  - `missed_cleavages: 0`
  - `min_charge: 2`, `max_charge: 2`
  - `add_decoys: true`
  - `max_koina_batch: 100` (smaller batches for testing)
- **Verification**:
  - Library directory created (.poin extension)
  - Peptides CSV file exists with expected peptides
  - Fragment index files created (*.arrow files)
  - Fragments have Altimeter spline coefficients
  - Verify tryptic digestion produces expected peptides
  - Check decoys are properly generated

#### Scenario C: Multi-file with Different Parameters + Altimeter
- **FASTA source**: dir1 (2 proteins) 
- **Parameters**:
  - `predict_fragments: true`
  - `prediction_model: "altimeter"`
  - `min_length: 8`, `max_length: 15` (moderate peptide count)
  - `missed_cleavages: 1`
  - `min_charge: 2`, `max_charge: 3`
  - `variable_mods: ["M"]` with oxidation
  - `max_koina_batch: 200`
- **Verification**:
  - Correct number of peptides with modifications
  - Fragment predictions for all charge states
  - Altimeter ion dictionary properly loaded
  - Fragment index contains spline coefficients
  - Missed cleavages handled properly

#### Scenario D: Comprehensive Test with All Features
- **FASTA source**: All test proteins (dir1 + dir2)
- **Parameters**:
  - `predict_fragments: true`
  - `prediction_model: "altimeter"`
  - `min_length: 7`, `max_length: 20`
  - `missed_cleavages: 2`
  - `min_charge: 2`, `max_charge: 4`
  - `variable_mods: ["M"]`
  - `add_decoys: true`
  - `entrapment_r: 0.1` (add entrapment sequences)
- **Verification**:
  - All proteins processed correctly
  - Entrapment sequences added
  - Fragment predictions complete
  - Library contains all expected components

### 3. Test Data Setup

#### Test FASTA Files to Create:

1. **minimal_protein.fasta** - For missed cleavage verification:
```fasta
>TEST001|MINI Mini protein for cleavage testing
MARALYKDERPK
```

2. **Keep existing files**:
   - single_protein.fasta
   - proteins1.fasta, proteins2.fasta (in dir1)
   - proteins3.fasta (in dir2)

#### Directory Structure:
```
data/
├── test_fasta_files/
│   ├── minimal_protein.fasta (new)
│   ├── single_protein.fasta (existing)
│   ├── fastas_dir1/ (existing)
│   └── fastas_dir2/ (existing)
├── test_build_spec_lib/
│   ├── scenario_a_missed_cleavages/
│   │   ├── params_0mc.json
│   │   ├── params_1mc.json
│   │   ├── params_2mc.json
│   │   ├── params_3mc.json
│   │   └── output/
│   ├── scenario_b_standard/
│   │   ├── params_altimeter.json
│   │   └── output/
│   ├── scenario_c_multifile/
│   │   ├── params_altimeter.json  
│   │   └── output/
│   └── scenario_d_comprehensive/
│       ├── params_altimeter.json
│       └── output/
```

### 4. Parameter JSON Configuration for Altimeter

Key settings for all scenarios:
```json
{
  "library_params": {
    "prediction_model": "altimeter",
    "auto_detect_frag_bounds": false,
    "frag_mz_min": 150.0,
    "frag_mz_max": 2000.0,
    "max_frag_rank": 50,
    "include_neutral_diff": true
  },
  "predict_fragments": true,
  "max_koina_requests": 12,
  "max_koina_batch": 100
}
```

### 5. Verification Functions

#### `verify_altimeter_library(lib_dir)`
- Check for Altimeter-specific files
- Verify fragment index format
- Validate spline coefficients present

#### `verify_fragment_predictions(lib_dir)`
- Load fragment index Arrow files
- Check fragments have intensities
- Verify Altimeter ion annotations
- Validate spline coefficient structure

#### `verify_peptide_generation(peptides_file, expected_params)`
- Parse peptides CSV
- Check peptide lengths and charges
- Verify modifications applied
- Count targets vs decoys

#### `verify_koina_integration(lib_dir)`
- Check build_log.txt for Koina API calls
- Verify successful fragment predictions
- No excessive retry warnings (using our new debug level)

### 6. Integration with runtests.jl
```julia
# Add after FASTA params test
@testset "BuildSpecLib Integration Tests" begin
    include("./Routines/BuildSpecLib/test_build_spec_lib.jl")
end
```

### 7. Handling Koina API Calls

- Use small peptide sets to minimize API calls
- Implement retry logic awareness
- Set reasonable timeouts
- Check for network availability (skip tests if offline)

### 8. Expected Outputs with Altimeter

For single_protein.fasta (MEPVDPRLEPWKHPGSQPKTACTNCYCKKCCFHCQVCFITKALGISYGRKKRRQRRRRAHQNSY - 64 AA):
- Tryptic sites at K and R positions
- With min_length=7, max_length=12:
  - Peptides: LEPWKHPGSQPK, TACTNCYCK, ALGISYGR, etc.
  - ~5-7 peptides meeting criteria
- With decoys: double the count
- Each peptide gets Altimeter fragment predictions:
  - b and y ions with spline coefficients
  - Neutral losses if enabled
  - Multiple charge states for fragments

### 9. Test Execution Flow

```julia
@testset "Scenario A - Minimal Altimeter" begin
    # 1. Generate params JSON
    # 2. Run BuildSpecLib
    # 3. Verify outputs
    # 4. Check fragment predictions
end
```

### 10. Performance Expectations

- Koina API calls: ~10-50 per scenario (depending on peptide count)
- Total test time: 30-60 seconds (including API calls)
- Network dependency: Required for Altimeter predictions
- Disk usage: ~10-50 MB per library

## Implementation Steps

1. Create test file with imports
2. Implement parameter JSON generators with Altimeter config
3. Add network connectivity check
4. Write Scenario A test (minimal)
5. Write Scenario B test (modifications)  
6. Write Scenario C test (comprehensive)
7. Add verification functions for Altimeter outputs
8. Integrate with runtests.jl
9. Document Altimeter-specific validations

## Key Validation Points for Altimeter

- Fragment index contains spline coefficients (not just intensities)
- Ion dictionary properly loaded from assets
- Fragment annotations follow Altimeter format
- API responses successfully parsed
- Coefficient arrays have correct dimensions
- Fragment m/z values match theoretical calculations

## Error Handling

- Network timeout handling
- Koina API rate limiting
- Invalid peptide sequences
- Missing ion dictionary file
- Malformed API responses

## Success Criteria

- All three scenarios complete without errors
- Fragment predictions present for all peptides
- Altimeter spline coefficients properly stored
- Decoys have predicted fragments
- Modified peptides handled correctly
- Test execution time < 2 minutes

## Test Peptide Examples

### Special Test Case: Ultra-short protein for missed cleavage verification
Create a new test file: `minimal_protein.fasta`
```
>TEST001|MINI Mini protein for cleavage testing
MARKTKAAR
```

This 9 AA sequence has 3 tryptic cleavage sites (after R3, K5, K7):
- **0 missed cleavages**: 
  - `MAR` (3 AA) - too short if min_length=7
  - `K` (1 AA) - too short
  - `TK` (2 AA) - too short  
  - `AAR` (3 AA) - too short
  - No valid peptides with min_length=7!

Let's use a better test sequence:
```
>TEST001|MINI Mini protein for cleavage testing  
MARALYKDERPK
```

This 12 AA sequence analysis:
- **0 missed cleavages** (cleave after R4, K7, R10, K12):
  - `MAR` (3 AA) - too short
  - `ALYK` (4 AA) - too short
  - `DER` (3 AA) - too short
  - `PK` (2 AA) - too short
  
- **1 missed cleavage**:
  - `MARALYK` (7 AA) - VALID ✓
  - `ALYKDER` (7 AA) - VALID ✓
  - `DERPK` (5 AA) - too short

- **2 missed cleavages**:
  - `MARALYKDER` (10 AA) - VALID ✓
  - `ALYKDERPK` (9 AA) - VALID ✓
  
- **3 missed cleavages**:
  - `MARALYKDERPK` (12 AA) - VALID ✓

This makes it very easy to verify:
- 0 missed cleavages = 0 peptides (all too short)
- 1 missed cleavage = 2 peptides
- 2 missed cleavages = 2 peptides (different ones)
- 3+ missed cleavages = 1 peptide (full sequence)

### From single_protein.fasta (Custom protein)
Sequence: `MEPVDPRLEPWKHPGSQPKTACTNCYCKKCCFHCQVCFITKALGISYGRKKRRQRRRRAHQNSY`

Expected tryptic peptides (0 missed cleavages, 7-12 AA):
1. `LEPWKHPGSQPK` (12 AA)
2. `TACTNCYCK` (9 AA) 
3. `ALGISYGR` (8 AA)

### From proteins1.fasta (P12345|PROT1)
Sequence: `MKLLSSIEQACDICRLKKLKCSKEKPKCAKCLKNNWECRYSPKTKRSPLTRAHLTEVESRLERL`

Expected tryptic peptides (0 missed cleavages, 7-12 AA):
1. `LLSSIEQACDICR` (12 AA, after initial MK cleavage)
2. `YSPKTKR` (7 AA)
3. `AHLTEVESRLER` (12 AA)

### From proteins2.fasta (P67890|PROT2)  
Sequence: `MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFL`

Expected tryptic peptides (0 missed cleavages, 7-12 AA):
1. `VLRDNIQGITK` (11 AA)
2. `ISGLIYEETR` (10 AA)

## Notes on Implementation

- The test should be self-contained and not rely on external test data beyond what we create
- Use try-catch blocks around Koina API calls to handle network issues gracefully
- Consider implementing a mock Koina server for offline testing (future enhancement)
- Fragment predictions should be validated against known Altimeter output characteristics
- The test should clean up after itself but leave artifacts for debugging if tests fail