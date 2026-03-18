#!/usr/bin/env julia
"""
Compare ScoringSearch results using FirstPassSearch PSMs vs SecondPassSearch PSMs.

Question: Is SecondPassSearch necessary? Or can we get comparable FDR-controlled IDs
directly from FirstPassSearch PSMs scored with a global LightGBM model?

Approach:
  1. Load first pass PSM tables (all features, one best PSM per precursor per file)
  2. Load second pass PSM tables (all features, re-deconvoluted with filtered precursor set)
  3. Run the same AdvancedLightGBM scoring (percolator_scoring!) on both
  4. Compare target IDs at 1% FDR

Conditions:
  1. FirstPass ALL precursors, shared features       (9M PSMs)
  2. FirstPass FILTERED to globally-passing precs    (~700K PSMs, apples-to-apples)
  3. SecondPass, shared features                     (702K PSMs)
  4. SecondPass, full features                       (702K PSMs)

Comparing conditions 2 vs 3 isolates the effect of re-deconvolution on the same precursor set.
"""

using Pioneer
using Arrow, DataFrames, Tables
using LightGBM

# ─── Configuration ────────────────────────────────────────────────────────────
const DATA_DIR = "/Users/nathanwamsley/Data/For_Figures/OlsenAstralThreeProteome200ng/OlsenAstralThreeProteome200ng_Standard/temp_data"
const FILES = ["E20H50Y30", "E45H50Y5", "E5H50Y45"]

# ─── Load PSMs ────────────────────────────────────────────────────────────────

function load_psms(subdir::String)
    dfs = DataFrame[]
    for (i, f) in enumerate(FILES)
        path = joinpath(DATA_DIR, subdir, "$f.arrow")
        if !isfile(path)
            @warn "Missing file: $path"
            continue
        end
        df = DataFrame(Tables.columntable(Arrow.Table(path)))
        # Ensure ms_file_idx is present and correct type
        if !hasproperty(df, :ms_file_idx)
            df[!, :ms_file_idx] = fill(UInt32(i), nrow(df))
        else
            df[!, :ms_file_idx] = UInt32.(df[!, :ms_file_idx])
        end
        # Fix broken ms_file_idx (always 0 in older first pass outputs)
        if all(==(UInt32(0)), df[!, :ms_file_idx])
            df[!, :ms_file_idx] .= UInt32(i)
        end
        push!(dfs, df)
    end
    return vcat(dfs...)
end

println("Loading first pass PSMs...")
fp_psms = load_psms("first_pass_psms")
println("  $(nrow(fp_psms)) PSMs, $(ncol(fp_psms)) columns")
println("  Columns: $(names(fp_psms))")

println("\nLoading second pass PSMs...")
sp_psms = load_psms("second_pass_all_features")
println("  $(nrow(sp_psms)) PSMs, $(ncol(sp_psms)) columns")
println("  Columns: $(names(sp_psms))")

# ─── Define feature sets ─────────────────────────────────────────────────────

# Features available in BOTH first and second pass
const SHARED_FEATURES = [
    # Deconv quality
    :spectral_contrast, :gof, :fitted_manhattan_distance,
    :max_matched_residual, :max_unmatched_residual,
    :matched_ratio, :poisson, :log2_intensity_explained,
    :err_norm, :weight,
    # Ion counts
    :y_count, :b_count, :isotope_count, :total_ions,
    :best_rank, :best_rank_iso, :topn, :topn_iso,
    :longest_y, :longest_b,
    # Peptide properties
    :missed_cleavage, :Mox, :sequence_length, :charge,
    # RT
    :irt_error, :irt_pred,
    # Other
    :spectrum_peak_count,
]

# Additional features only in second pass
const SECOND_PASS_EXTRA = [
    :fitted_spectral_contrast, :percent_theoretical_ignored, :scribe,
    :irt_diff, :tic, :prec_mz, :adjusted_intensity_explained,
    # Multi-scan aggregates
    :num_scans, :max_y_ions, :y_ions_sum,
    :max_gof, :max_fitted_manhattan_distance, :max_fitted_spectral_contrast,
    :max_matched_ratio, :max_scribe, :max_weight,
    # Amino acids
    :aa_H, :aa_P, :aa_L,
    # MS1
    :ms1_ms2_rt_diff, :ms1_features_missing,
]

# Filter to features actually present in each DataFrame
function available_features(df::DataFrame, features::Vector{Symbol})
    return [f for f in features if hasproperty(df, f)]
end

const SCRIBE_FEATURES = Set([:scribe, :max_scribe])

fp_features = available_features(fp_psms, SHARED_FEATURES)
sp_features_shared = available_features(sp_psms, SHARED_FEATURES)
sp_features_full = available_features(sp_psms, vcat(SHARED_FEATURES, SECOND_PASS_EXTRA))
sp_features_no_scribe = [f for f in sp_features_full if f ∉ SCRIBE_FEATURES]

println("\n=== Feature availability ===")
println("First pass:  $(length(fp_features)) features")
println("Second pass (shared only): $(length(sp_features_shared)) features")
println("Second pass (full): $(length(sp_features_full)) features")
println("Second pass (no scribe): $(length(sp_features_no_scribe)) features")

# Show which shared features are missing from first pass
missing_fp = setdiff(SHARED_FEATURES, fp_features)
if !isempty(missing_fp)
    println("\nFeatures in SHARED_FEATURES missing from first pass: $missing_fp")
end

# ─── Build filtered first pass (apples-to-apples) ───────────────────────────

# Get the set of precursors present in second pass (globally-passing set)
sp_precursors = Set(sp_psms.precursor_idx)
println("\nSecondPass precursor set: $(length(sp_precursors)) unique precursors")

# Filter first pass to only those precursors
fp_filtered = fp_psms[in.(fp_psms.precursor_idx, Ref(sp_precursors)), :]
println("FirstPass filtered to passing precursors: $(nrow(fp_filtered)) PSMs (was $(nrow(fp_psms)))")

# ─── Scoring function ────────────────────────────────────────────────────────

function score_psms_lgbm(psms::DataFrame, features::Vector{Symbol}; label::String="")
    println("\n━━━ Scoring: $label ━━━")
    println("  $(nrow(psms)) PSMs, $(length(features)) features")

    # Ensure required columns
    if !hasproperty(psms, :target)
        error("PSMs must have :target column")
    end
    if !hasproperty(psms, :cv_fold)
        error("PSMs must have :cv_fold column")
    end

    n_targets = count(psms.target)
    n_decoys = nrow(psms) - n_targets
    println("  Targets: $n_targets, Decoys: $n_decoys")

    # Initialize scoring columns
    psms[!, :trace_prob] = zeros(Float32, nrow(psms))
    psms[!, :q_value] = zeros(Float64, nrow(psms))
    if !hasproperty(psms, :decoy)
        psms[!, :decoy] = .!psms.target
    end

    # Need pair_id and irt_bin_idx for percolator_scoring!
    # We'll use NoPairing equivalent: each precursor gets its own pair_id
    if !hasproperty(psms, :pair_id)
        psms[!, :pair_id] = UInt32.(1:nrow(psms))
    end
    if !hasproperty(psms, :isotopes_captured)
        psms[!, :isotopes_captured] = fill((Int8(0), Int8(0)), nrow(psms))
    end

    # Wrap in PSMContainer for sort_of_percolator!
    psm_container = Pioneer.DataFramePSMContainer(psms)

    # Use sort_of_percolator! which handles all the scoring logic
    models = Pioneer.sort_of_percolator!(
        psm_container,
        features,
        false;  # match_between_runs = false (no MBR for this comparison)
        feature_fraction = 0.5,
        learning_rate = 0.05,
        min_data_in_leaf = 500,
        bagging_fraction = 0.25,
        min_gain_to_split = 0.5,
        max_depth = 10,
        num_leaves = 63,
        iter_scheme = [200],
        print_importance = true,
        show_progress = true
    )

    # Compute q-values from trace_prob
    qvals = Vector{Float64}(undef, nrow(psms))
    Pioneer.get_qvalues!(psms.trace_prob, psms.target, qvals)
    psms[!, :q_value] = qvals

    # Count IDs at various FDR thresholds
    for thresh in [0.01, 0.05, 0.10]
        n_pass = count((qvals .<= thresh) .& psms.target)
        println("  Targets at $(Int(thresh*100))% FDR: $n_pass")
    end

    return psms, models
end

# ─── Run comparisons ─────────────────────────────────────────────────────────

# 1. First pass PSMs with shared features (ALL precursors)
fp_scored, fp_models = score_psms_lgbm(
    copy(fp_psms), fp_features;
    label="FirstPass ALL precursors (shared features)"
)

# 2. First pass PSMs FILTERED to globally-passing precursors (apples-to-apples)
fp_filt_scored, fp_filt_models = score_psms_lgbm(
    copy(fp_filtered), available_features(fp_filtered, SHARED_FEATURES);
    label="FirstPass FILTERED to passing precs (shared features)"
)

# 3. Second pass PSMs with shared features only
sp_scored_shared, sp_models_shared = score_psms_lgbm(
    copy(sp_psms), sp_features_shared;
    label="SecondPass PSMs (shared features only)"
)

# 4. Second pass PSMs with full feature set
sp_scored_full, sp_models_full = score_psms_lgbm(
    copy(sp_psms), sp_features_full;
    label="SecondPass PSMs (full features)"
)

# 5. Second pass PSMs with full features minus scribe
sp_scored_no_scribe, sp_models_no_scribe = score_psms_lgbm(
    copy(sp_psms), sp_features_no_scribe;
    label="SecondPass PSMs (full features, no scribe)"
)

# ─── Summary comparison ──────────────────────────────────────────────────────

println("\n" * "="^80)
println("COMPARISON SUMMARY")
println("="^80)

function count_at_fdr(psms::DataFrame, thresh::Float64)
    return count((psms.q_value .<= thresh) .& psms.target)
end

println("\n$(rpad("Condition", 55)) $(rpad("1% FDR", 10)) $(rpad("5% FDR", 10)) 10% FDR")
println("-"^85)

conditions = [
    ("1. FirstPass ALL precs (shared features)", fp_scored),
    ("2. FirstPass FILTERED to passing precs (shared feat)", fp_filt_scored),
    ("3. SecondPass (shared features)", sp_scored_shared),
    ("4. SecondPass (full features)", sp_scored_full),
    ("5. SecondPass (full features, no scribe)", sp_scored_no_scribe),
]

for (label, scored) in conditions
    n1 = count_at_fdr(scored, 0.01)
    n5 = count_at_fdr(scored, 0.05)
    n10 = count_at_fdr(scored, 0.10)
    println("$(rpad(label, 55)) $(rpad(n1, 10)) $(rpad(n5, 10)) $n10")
end

# Per-file breakdown using actual ms_file_idx values
println("\n--- Per-file breakdown at 1% FDR ---")

# Get unique file indices from each dataset
all_file_idxs = sort(unique(vcat(
    unique(fp_scored.ms_file_idx),
    unique(fp_filt_scored.ms_file_idx),
    unique(sp_scored_shared.ms_file_idx),
    unique(sp_scored_full.ms_file_idx),
    unique(sp_scored_no_scribe.ms_file_idx),
)))

println(rpad("ms_file_idx", 15) * rpad("FP_all", 18) * rpad("FP_filt", 18) * rpad("SP_shared", 18) * rpad("SP_full", 18) * rpad("SP_no_scribe", 18))
println("-"^105)

for fidx in all_file_idxs
    row = rpad(string(fidx), 15)
    for (_, scored) in conditions
        mask = scored.ms_file_idx .== fidx
        n = count((scored.q_value[mask] .<= 0.01) .& scored.target[mask])
        row *= rpad(n, 18)
    end
    println(row)
end

# Unique precursors across all files
println("\n--- Unique precursors at 1% FDR ---")
fp_passing = unique(fp_scored[fp_scored.q_value .<= 0.01 .&& fp_scored.target, :precursor_idx])
fp_f_passing = unique(fp_filt_scored[fp_filt_scored.q_value .<= 0.01 .&& fp_filt_scored.target, :precursor_idx])
sp_s_passing = unique(sp_scored_shared[sp_scored_shared.q_value .<= 0.01 .&& sp_scored_shared.target, :precursor_idx])
sp_f_passing = unique(sp_scored_full[sp_scored_full.q_value .<= 0.01 .&& sp_scored_full.target, :precursor_idx])
sp_ns_passing = unique(sp_scored_no_scribe[sp_scored_no_scribe.q_value .<= 0.01 .&& sp_scored_no_scribe.target, :precursor_idx])

println("1. FirstPass ALL (shared):        $(length(fp_passing)) unique precursors")
println("2. FirstPass FILTERED (shared):   $(length(fp_f_passing)) unique precursors")
println("3. SecondPass (shared):           $(length(sp_s_passing)) unique precursors")
println("4. SecondPass (full):             $(length(sp_f_passing)) unique precursors")
println("5. SecondPass (full, no scribe):  $(length(sp_ns_passing)) unique precursors")

println("\n--- Key comparisons ---")
println("Effect of noise reduction (1 vs 2):  $(length(fp_f_passing) - length(fp_passing)) precursors")
println("Effect of re-deconv (2 vs 3):        $(length(sp_s_passing) - length(fp_f_passing)) precursors")
println("Effect of extra features (3 vs 4):   $(length(sp_f_passing) - length(sp_s_passing)) precursors")
println("Effect of scribe (5 vs 4):           $(length(sp_f_passing) - length(sp_ns_passing)) precursors")

println("\n--- Overlap analysis (2 vs 3, apples-to-apples) ---")
println("Overlap FP_filt ∩ SP(shared):  $(length(intersect(fp_f_passing, sp_s_passing))) precursors")
println("Only in FP_filt:               $(length(setdiff(fp_f_passing, sp_s_passing))) precursors")
println("Only in SP(shared):            $(length(setdiff(sp_s_passing, fp_f_passing))) precursors")

println("\n--- Overlap analysis (1 vs 4, overall) ---")
println("Overlap FP_all ∩ SP(full):     $(length(intersect(fp_passing, sp_f_passing))) precursors")
println("Only in FP_all:                $(length(setdiff(fp_passing, sp_f_passing))) precursors")
println("Only in SP(full):              $(length(setdiff(sp_f_passing, fp_passing))) precursors")
