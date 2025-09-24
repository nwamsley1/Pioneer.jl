# Parameter Configuration

Pioneer.jl uses JSON configuration files to control analysis. This guide explains the parameters for both SearchDIA and BuildSpecLib functions.

## SearchDIA Configuration

Pioneer.jl uses JSON configuration files to control analysis. This guide explains the parameters for both SearchDIA and BuildSpecLib functions.

## SearchDIA Configuration

### Frequently Modified Parameters


Most parameters should not be changed, but the following may need adjustement. 

* `first_search.fragment_settings.min_score`: The minimum score determines which fragments must match in the fragment-index search in order for the precursor to pass. Each precursor is awarded a score based on which fragments match the spectrum. The score assigned to each fragment depends on its intensity rank. The default scheme is 8,4,4,2,2,1,1. That is, if the 1st, 3rd, and 7th ranking fragments matched the spectrum, the precursor would be awarded a score of 8+4+1=13. If all 7 of the fragments matched, the precursor would be awarded a score of 22. For normal instrument settings on an Orbitrap or Astral mass analyzer, the mass tolerance is about +/- 5-15 ppm and 15 is a reasonable default score threshold. However, for instruments with less mass accuracy (Sciex ZenoTOF 7600 or different Orbitrap scan settings), the score threshold may need to be set higher, perhaps to 20. It may be worthwile to test different values when searching data from a new instrument or sample type. In order to pass the first search, a precursor need only pass the threshold and score sufficiently well in at least one of the MS data files.

* `first_search.fragment_settings.max_rank`: Search against only the n'th most abundant fragment for each precursor. Including more fragments can improve performance but increase memory consumption, and the search could take longer. From experience, there are diminishing returns after 25-50 fragments. 

* `quant_search.fragment_settings.max_rank`: See above 

* `quant_search.fragment_settings.n_isotopes`: If searching with non-Altimeter libraries (not recommended), such as Prosit or UniSpec, this should be set to 1 as the second fragment isotopes will not be calculated accurately.

* `acquisition.nce`: This is the initial guess for the normalized collision energy that will best align the Altimeter Library with the empirical data. Altimeter values should agree with those from Thermo Instruments manufactured in Bremen Germany. If upon inspection of the quality control plots the initial guess is far from the estimated value, it might be possible to improve search results slightly by re-searching with a better initial guess.

* `acquisition.quad_transmission.fit_from_data`: Estimate the quad transmission function from the data. Otherwise defaults to symmetric, smooth function. 

* `optimization.machine_learning.max_samples`: This is the maximum number of PSMs to use for training the LightGBM model. These PSMs need to comfortably fit in memory in addition to the spectral library. As a rule of thumb, 7M rows is about 1GB. At the default maximum of 50M rows, the PSMs table will consume 7GB of memory.
* During LightGBM training, any missing feature values are replaced with the column median. If a column is entirely missing, the values are filled with zero of the appropriate type.

* `global.isotope_settings.combine_traces`: Some precursors may be split accross different acquisition windows. Pioneer refers to these as seperate isotope traces. When set to true, Pioneer does not distinguish between a precursor's isotope traces. They are combined for scoring and quantitation. With a clever acquisition scheme this can increase the number of data points accross chromatographic peaks. This is recomended only for acquisition windows 2-4 m/z. It should also be combined with `aquisition.quad_transmission.fit_from_data` = true. 


### Global Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `isotope_settings.err_bounds_first_pass` | [Int, Int] | Precursor monoisotope may lie NEUTRON/charge Thompsons (left, right) outside the quadrupole isolation window (default: [1, 0]) |
| `isotope_settings.err_bounds_second_pass` | [Int, Int] | Precursor monoisotope may lie NEUTRON/charge Thompsons (left, right) outside the quadrupole isolation window (default: [3, 1]) |
| `isotope_settings.combine_traces` | Boolean | Whether to combine precursor isotope traces in quantification. Experimental, so set to false (default: false) |
| `isotope_settings.partial_capture` | Boolean | Whether to estimate the conditional fragment isotope distribution (true) or assume complete transmission the entire precursor isotopic envelope (default: true) |
| `isotope_settings.min_fraction_transmitted` | Float | Minimum fraction of the precursor isotope distribution that must be isolated for scoring and quantitation (default: 0.25) |
| `scoring.q_value_threshold` | Float | Global q-value threshold for filtering results. Also controls false transfer rate of MBR (default: 0.01) |
| `normalization.n_rt_bins` | Int | Number of retention time bins for quant normalization (default: 100) |
| `normalization.spline_n_knots` | Int | Number of knots in quant normalization spline (default: 7) |
| `match_between_runs` | Boolean | Whether to attempt to transfer peptide identifications across runs. Turning this on will add additional features to the LightGBM model (default: true) |

### Parameter Tuning Settings

| Parameter | Type | Description |
|-----------|------|-------------|
| `fragment_settings.min_count` | Int | Minimum number of matching fragment ions (default: 7) |
| `fragment_settings.max_rank` | Int | Maximum rank of fragments to consider (default: 25, means 26th-last most abundant fragments per precursor are filtered out) |
| `fragment_settings.tol_ppm` | Float | Initial tragment mass tolerance guess in parts per million (default: 20.0, should be set lower for some TOF instruments) |
| `fragment_settings.min_score` | Int | Minimum fragment-index score threshold for fragment matches (default: 22) |
| `fragment_settings.min_spectral_contrast` | Float | Minimum cosine simmilarity score (default: 0.9) |
| `fragment_settings.relative_improvement_threshold` | Float | Minimum relative Scribe score improvement needed to ignore an interferring peak (default: 1.25) |
| `fragment_settings.min_log2_ratio` | Float | Minimum log2 ratio of matched library fragment intensities to unmatched library fragment intensities (default: 1.5) |
| `fragment_settings.min_top_n` | [Int, Int] | Minimum number of top N matches - [requirement, denominator]. Default: `[3, 3]` |
| `fragment_settings.n_isotopes` | Int | Number of fragment isotopes to consider in matching (default: 1, mono only) |
| `search_settings.min_samples` | Int | Minimum number of samples required for tuning (default: 3500) |
| `search_settings.min_quad_tuning_psms` | Int | Minimum number of psms required for estimating quad transmission (default: 5000) |
| `search_settings.min_quad_tuning_fragments` | Int | Must match at least n fragments to each quad tuning psm (default: 3) |
| `search_settings.max_presearch_iters` | Int | Maximum number of parameter tuning iterations (default: 10) |
| `search_settings.frag_err_quantile` | Float | Quantile for fragment error estimation (default: 0.01) |

### First Search Parameters 

| Parameter | Type | Description |
|-----------|------|-------------|
| `fragment_settings.min_count` | Int | Minimum number of matching fragments (default: 4) |
| `fragment_settings.max_rank` | Int | Maximum fragment rank to consider (default: 50 means 50th-last most abundant fragments per precursor are filtered out) |
| `fragment_settings.min_score` | Int | Minimum score for fragment matches (default: 15) |
| `fragment_settings.min_spectral_contrast` | Float | Minimum cosine simmilarity required (default: 0.5) |
| `fragment_settings.relative_improvement_threshold` | Float | Minimum relative Scribe score improvement needed to ignore an interferring peak (default: 1.25) |
| `fragment_settings.min_log2_ratio` | Float | Minimum log2 ratio of matched library fragment intensities to unmatched library fragment intensities (default: 0.0, means sum of matched library fragment intensities is equal to the sum of unmatched library fragment intensities for the precursor ) |
| `fragment_settings.min_top_n` | [Int, Int] | Minimum top N matches - [requirement, denominator]. Default: `[2, 3]` |
| `fragment_settings.n_isotopes` | Int | Number of isotopes to consider (default: 1) |
| `scoring_settings.n_train_rounds` | Int | Number of training rounds for scoring model (default: 2) |
| `scoring_settings.max_iterations` | Int | Maximum iterations for scoring optimization (default: 20) |
| `scoring_settings.max_q_value_probit_rescore` | Float | Maximum q-value threshold for semi-supervised learning durning probit regression (default: 0.05) |
| `scoring_settings.max_PEP` | Int | Maximum local FDR threshold for passing the first search (default: 0.9) |
| `irt_mapping.max_prob_to_impute_irt` | Int | If probability of the psm is less then x in the first-pass search, then impute irt for the precursor with globably determined value from the other runs (default: 0.75) |
| `irt_mapping.fwhm_nstd` | Float | Number of standard deviations of the fwhm to add to the retention time tolerance (default: 4) |
| `irt_mapping.irt_nstd` | Int | Number of standard deviations of run-to-run irt tolerance to add to the retention time tolerance (default: 4) |

### Quantification Search Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `fragment_settings.min_count` | Int | Minimum fragment count for quantification (default: 3) |
| `fragment_settings.min_y_count` | Int | Minimum number of y-ions required (default: 2) |
| `fragment_settings.max_rank` | Int | Maximum fragment rank (default: 255) |
| `fragment_settings.min_spectral_contrast` | Float | Minimum spectral contrast score (default: 0.0) |
| `fragment_settings.min_log2_ratio` | Float | Minimum log2 ratio of intensities (default: -1.7) |
| `fragment_settings.min_top_n` | [Int, Int] | Minimum top N matches - [requirement, denominator]. Default: `[2, 3]` |
| `fragment_settings.n_isotopes` | Int | Number of isotopes for quantification (default: 2, include the M1 and M2 isotopes) |
| `chromatogram.smoothing_strength` | Float | Strength of chromatogram smoothing (default: 1e-6) |
| `chromatogram.padding` | Int | Number of zeros to pad chromatograms on either side (default: 0) |
| `chromatogram.max_apex_offset` | Int | Maximum allowed apex offset in #scans where the precursor could have been detected between the second-pass search and re-integration with 1 percent FDR precursors (default: 2) |

### Acquisition Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `nce` | Int | Normalized collision energy initial guess (used in pre-search before NCE tuning) (default: 25) |
| `quad_transmission.fit_from_data` | Boolean | Whether to fit quadrupole transmission from data (default: false)|
| `quad_transmission.overhang` | Float | deprecated (default: 0.25) |
| `quad_transmission.smoothness` | Float | Smoothness parameter for transmission curve. Higher value means more "box-like" shape. (default: 5.0) |

### RT Alignment Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `sigma_tolerance` | Int | Number of standard deviations for irt tolerance after pre-search (default: 4) |
| `min_probability` | Float | Minimum probability for alignment psms in pre-search (default: 0.95) |

### Optimization Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `deconvolution.lambda` | Float | L2 regularization parameter for deconvolution (default: 0.0 for NoNorm regularization) |
| `deconvolution.reg_type` | String | Regularization type: "none" (NoNorm), "l1" (L1Norm), or "l2" (L2Norm) (default: "none") |
| `deconvolution.huber_delta` | Float | Delta parameter for Huber loss function (default: 300) |
| `deconvolution.huber_exp` | Float | Exponent for Huber delta progression (default: 1.5) |
| `deconvolution.huber_iters` | Int | Number of Huber outer iterations (default: 15) |
| `deconvolution.newton_iters` | Int | Maximum Newton iterations per outer iteration (recommended: 25, default: 50) |
| `deconvolution.bisection_iters` | Int | Maximum bisection iterations when Newton fails (recommended: 100, default: 100) |
| `deconvolution.outer_iters` | Int | Maximum outer iterations for convergence (recommended: max(1000, n_variables*5), default: 1000) |
| `deconvolution.newton_accuracy` | Float | Absolute convergence threshold for Newton method (default: 10) |
| `deconvolution.bisection_accuracy` | Float | Absolute convergence threshold for bisection method (default: 10) |
| `deconvolution.max_diff` | Float | Relative convergence threshold - maximum relative change in weights between iterations. Also used as relative tolerance for Newton's method (default: 0.01) |
| `machine_learning.max_samples` | Int | Maximum number of samples for LightGBM training (default: 5000000) |
| `machine_learning.min_trace_prob` | Float | Minimum trace probability threshold (default: 0.75) |
| `machine_learning.max_q_value_mbr_rescore` | Float | q-value threshold for match-between-runs candidates during semi-supervised learning with LightGBM (default: 0.20) |
| `machine_learning.min_PEP_neg_threshold_rescore` | Float | Minimum posterior error probabilility threshold for poor scoring targets to be relabeled as negative examples during semi-supervised learning with LightGBM (default: 0.20) |
| `machine_learning.spline_points` | Int | Number of points for probability spline (default: 500) |
| `machine_learning.interpolation_points` | Int | Number of interpolation points (default: 10) |

### Protein Inference Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `min_peptides` | Int | Minimum number of peptides required for a protein group (default: 1) |

### MaxLFQ Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `run_to_run_normalization` | Boolean | Whether to use run-to-run normalized abundances for precursor and protein quantification (default: true) |

### Output Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `write_csv` | Boolean | Whether to write results to CSV |
| `write_decoys` | Boolean | Where to quantify and include decoys in the output files |
| `delete_temp` | Boolean | Whether to delete temporary files |
| `plots_per_page` | Int | Number of plots per page in reports (default: 12) |

### Path Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `library` | String | Path to spectral library file |
| `ms_data` | String | Path to mass spectrometry data directory |
| `results` | String | Path to output results directory |

## BuildSpecLib Configuration

### FASTA Input and Regex Mapping

Pioneer.jl now supports flexible FASTA input through `GetBuildLibParams`:

#### Input Options
1. **Single directory**: Scans for all `.fasta` and `.fasta.gz` files
2. **Single file**: Directly uses the specified FASTA file
3. **Mixed array**: Any combination of directories and files

#### Regex Code Mapping
The regex patterns for parsing FASTA headers can be configured in three ways:

1. **Single regex set for all files** (default):
   ```julia
   GetBuildLibParams(out_dir, lib_name, [dir1, dir2, file1])
   # All FASTA files use the same default regex patterns
   ```

2. **Custom single regex set**:
   ```julia
   GetBuildLibParams(out_dir, lib_name, [dir1, file1],
       regex_codes = Dict(
           "accessions" => "^>(\\S+)",
           "genes" => "GN=(\\S+)",
           "proteins" => "\\s+(.+?)\\s+OS=",
           "organisms" => "OS=(.+?)\\s+GN="
       ))
   # All files use these custom patterns
   ```

3. **Positional mapping** (one regex set per input):
   ```julia
   GetBuildLibParams(out_dir, lib_name, [uniprot_dir, custom_file],
       regex_codes = [
           Dict("accessions" => "^\\w+\\|(\\w+)\\|", ...),  # For uniprot_dir files
           Dict("accessions" => "^>(\\S+)", ...)             # For custom_file
       ])
   ```

### FASTA Digest Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `min_length` | Int | Minimum peptide length (default: 7) |
| `max_length` | Int | Maximum peptide length (default: 30) |
| `min_charge` | Int | Minimum charge state (default: 2) |
| `max_charge` | Int | Maximum charge state (default: 4) |
| `cleavage_regex` | String | Regular expression for cleavage sites (default: "[KR][^_\|$]", to exclude cleavage after proline: "[KR][^P|$]") |
| `missed_cleavages` | Int | Maximum allowed missed cleavages (default: 1) |
| `max_var_mods` | Int | Maximum variable modifications per peptide (default: 1) |
| `add_decoys` | Boolean | Generate decoy sequences (default: true) |
| `entrapment_r` | Float | Ratio of entrapment sequences (default: 0) |
| `fasta_header_regex_accessions` | [String] | Regex with a capture group for the accession, one per FASTA file |
| `fasta_header_regex_genes` | [String] | Regex with a capture group for the gene name, one per FASTA file |
| `fasta_header_regex_proteins` | [String] | Regex with a capture group for the protein name, one per FASTA file |
| `fasta_header_regex_organisms` | [String] | Regex with a capture group for the organism, one per FASTA file |

### NCE Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `nce` | Float | Base normalized collision energy (default: 25.0) |
| `default_charge` | Int | Default charge state for NCE calculations (default: 2) |
| `dynamic_nce` | Boolean | Use charge-dependent NCE adjustments (default: true) |

### Library Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `rt_bin_tol` | Float | Retention time binning tolerance in minutes (default: 1.0) |
| `frag_bin_tol_ppm` | Float | Fragment mass tolerance in PPM (default: 10.0) |
| `rank_to_score` | [Int] | Intensity multipliers for ranked peaks (default: [8,4,4,2,2,1,1]) |
| `y_start_index` | Int | Starting index for y-ion annotation (default: 4) |
| `b_start_index` | Int | Starting index for b-ion annotation (default: 3) |
| `y_start` | Int | Minimum y-ion to consider (default: 3) |
| `b_start` | Int | Minimum b-ion to consider (default: 2) |
| `include_p_index` | Boolean | Include proline-containing index fragments (default: false) |
| `include_p` | Boolean | Include proline-containing fragments (default: false) |
| `auto_detect_frag_bounds` | Boolean | Auto-detect fragment mass bounds (default: true) |
| `calibration_raw_file` | String | Path to calibration raw file |
| `frag_mz_min` | Float | Minimum fragment m/z (default: 150.0) |
| `frag_mz_max` | Float | Maximum fragment m/z (default: 2020.0) |
| `prec_mz_min` | Float | Minimum precursor m/z (default: 390.0) |
| `prec_mz_max` | Float | Maximum precursor m/z (default: 1010.0) |
| `max_frag_charge` | Int | Maximum fragment ion charge (default: 3) |
| `max_frag_rank` | Int | Maximum fragment rank (default: 50) |
| `min_frag_intensity` | Float | Minimum relative fragment intensity (default: 0.00) |
| `include_isotope` | Boolean | Include isotope peak annotations (default: false) |
| `include_internal` | Boolean | Include internal fragment annotations (default: false) |
| `include_immonium` | Boolean | Include immonium ion annotations (default: false) |
| `include_neutral_diff` | Boolean | Include neutral loss annotations (default: true) |
| `instrument_type` | String | Instrument type for predictions (default: "NONE") |
| `prediction_model` | String | Model for fragment predictions (default: "altimeter") |

### Modification Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `variable_mods.pattern` | [String] | Amino acids to modify (default: ["M"]) |
| `variable_mods.mass` | [Float] | Modification masses (default: [15.99491]) |
| `variable_mods.name` | [String] | Modification identifiers (default: ["Unimod:35"]) |
| `fixed_mods.pattern` | [String] | Amino acids to modify (default: ["C"]) |
| `fixed_mods.mass` | [Float] | Modification masses (default: [57.021464]) |
| `fixed_mods.name` | [String] | Modification identifiers (default: ["Unimod:4"]) |

### Processing Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `max_koina_requests` | Int | Maximum concurrent Prosit API requests (default: 24) |
| `max_koina_batch` | Int | Maximum batch size for API requests (default: 1000) |
| `match_lib_build_batch` | Int | Batch size for library building (default: 100000) |

!!! note "Koina API Retry Behavior"
    As of version 0.1.13, Koina API retry warnings are now logged at debug level 2 instead of being shown to users by default. To see retry attempts during debugging, set `debug_console_level: 2` in your SearchDIA parameters. The library build will only fail if all retry attempts are exhausted.

### Path Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `fasta_paths` | [String] | List of FASTA file paths |
| `fasta_names` | [String] | Names for each FASTA file |
| `include_contaminants` | Boolean | Append a contaminants FASTA to the build (default: true) |
| `out_dir` | String | Output directory path |
| `lib_name` | String | Base name for library files |
| `new_lib_name` | String | Name for updated library files |
| `out_name` | String | Output filename |
| `predict_fragments` | Boolean | Predict fragment intensities (default: true) |

