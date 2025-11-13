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

"""
iRT Refinement Utilities

Provides functionality to refine library iRT predictions using linear regression
models trained on amino acid composition and observed retention time errors.
"""

using DataFrames, StatsModels, GLM, Statistics

# Standard amino acids for feature extraction (sorted for consistent ordering)
const STANDARD_AAS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Create mapping from amino acid to index for fast lookup
const AA_TO_IDX = Dict{Char, Int}(aa => i for (i, aa) in enumerate(STANDARD_AAS))

"""
    IrtRefinementModel

Callable linear model for iRT refinement with zero-allocation application.

# Fields
- `use_refinement::Bool` - Whether to apply refinement (based on validation performance)
- `aa_weights::Dict{Char, Float32}` - Direct AA character → weight mapping for O(1) lookup
- `intercept::Float32` - Model intercept
- `irt_coef::Float32` - Coefficient for library_irt feature
- `mae_original::Float32` - MAE before refinement
- `mae_refined::Float32` - MAE after refinement
- `r2_train::Float32` - R² on training set
- `r2_val::Float32` - R² on validation set

# Callable Interface
The model can be called directly on sequences:
```julia
refined_irt = model(sequence::String, library_irt::Float32)
```

This zero-allocation design eliminates the need for amino acid count buffers.
"""
struct IrtRefinementModel
    use_refinement::Bool
    aa_weights::Dict{Char, Float32}  # AA → weight mapping (zero-allocation lookup)
    intercept::Float32
    irt_coef::Float32
    mae_original::Float32
    mae_refined::Float32
    r2_train::Float32
    r2_val::Float32
end

"""
    (model::IrtRefinementModel)(sequence::String, library_irt::Float32) -> Float32

Apply iRT refinement to a peptide sequence. Zero-allocation implementation.

# Algorithm
1. If refinement disabled, return library_irt
2. Start with intercept + (irt_coef × library_irt)
3. Loop through sequence, accumulate AA weights (single pass!)
4. Return: library_irt - predicted_error

# Performance
- Zero allocations (no counts buffer needed)
- Single loop through sequence
- O(1) Dict lookup per amino acid
- ~2x faster than array-based implementation

# Arguments
- `sequence`: Peptide sequence string
- `library_irt`: Library iRT prediction

# Returns
Refined iRT value (or library_irt if refinement disabled)
"""
function (model::IrtRefinementModel)(sequence::String, library_irt::Float32)::Float32
    # Early return if refinement disabled
    if !model.use_refinement
        return library_irt
    end

    # Start with intercept + library contribution
    predicted_error = model.intercept + model.irt_coef * library_irt

    # Single loop: accumulate AA weights directly (zero allocation!)
    for aa in sequence
        # Dict lookup returns 0.0f0 if AA not in dict (handles non-standard AAs)
        predicted_error += get(model.aa_weights, aa, 0.0f0)
    end

    # Refined iRT = library - predicted_error
    return library_irt - predicted_error
end

"""
    count_amino_acids!(counts::Vector{Int}, sequence::String)

Count occurrences of each standard amino acid in a peptide sequence.
Updates the pre-allocated counts vector in-place (indexed by STANDARD_AAS order).
Non-standard amino acids are ignored.

**Note**: This function is used for model training (feature extraction).
For model application (refinement), use the callable interface which has
zero allocations: `refined_irt = model(sequence, library_irt)`
"""
function count_amino_acids!(counts::Vector{Int}, sequence::String)
    # Reset counts to zero
    fill!(counts, 0)

    # Count each amino acid
    for char in sequence
        idx = get(AA_TO_IDX, char, 0)
        if idx > 0
            counts[idx] += 1
        end
    end

    return nothing
end

"""
    prepare_features_dataframe(sequences::Vector{String},
                               irt_predicted::Vector{Float32},
                               irt_errors::Vector{Float32}) -> DataFrame

Prepare feature matrix for iRT refinement model training.

# Arguments
- `sequences` - Peptide sequences
- `irt_predicted` - Library iRT predictions
- `irt_errors` - Prediction errors (irt_predicted - irt_observed)

# Returns
DataFrame with columns:
- `irt_predicted` - Library iRT values
- `count_A`, `count_C`, ... - Counts for each of 20 standard amino acids
- `error` - Target variable (iRT error to predict)
"""
function prepare_features_dataframe(sequences::Vector{String},
                                   irt_predicted::Vector{Float32},
                                   irt_errors::Vector{Float32})
    n = length(sequences)
    n_aas = length(STANDARD_AAS)

    # Initialize feature matrix with pre-allocated arrays
    # Convert to Float64 for GLM compatibility
    features = Dict{Symbol, Vector}()
    features[:irt_predicted] = Float64.(irt_predicted)
    features[:error] = Float64.(irt_errors)

    # Pre-allocate amino acid count columns
    aa_count_arrays = [zeros(Int, n) for _ in 1:n_aas]
    for (i, aa) in enumerate(STANDARD_AAS)
        col_name = Symbol("count_", aa)
        features[col_name] = aa_count_arrays[i]
    end

    # Reusable counts vector
    counts = zeros(Int, n_aas)

    # Count amino acids for each sequence using in-place counting
    for (seq_idx, sequence) in enumerate(sequences)
        count_amino_acids!(counts, sequence)
        for (aa_idx, count) in enumerate(counts)
            aa_count_arrays[aa_idx][seq_idx] = count
        end
    end

    return DataFrame(features)
end

"""
    fit_irt_refinement_model(sequences::Vector{String},
                            irt_predicted::Vector{Float32},
                            irt_observed::Vector{Float32};
                            ms_file_idx::Int=1,
                            min_psms::Int=20,
                            train_fraction::Float64=0.67) -> Union{IrtRefinementModel, Nothing}

Train a linear regression model to predict library iRT prediction errors based on amino acid composition.

Model predicts: error = irt_predicted - irt_observed (positive when library overestimates)
Corrected iRT = irt_predicted - predicted_error

# Workflow
1. Filter to sequences with sufficient PSMs (min_psms)
2. Split into training (train_fraction) and validation sets
3. Train linear model on training set
4. Evaluate on validation set
5. If validation MAE improves, retrain on full dataset for final model
6. Return model if refinement helps, nothing otherwise

# Arguments
- `sequences` - Peptide sequences
- `irt_predicted` - Library iRT predictions
- `irt_observed` - Observed iRT values from RT alignment
- `ms_file_idx` - File index for logging (default: 1)
- `min_psms` - Minimum PSMs required per feature (default: 20)
- `train_fraction` - Fraction of data for training (default: 0.67)

# Returns
`IrtRefinementModel` if refinement improves MAE, `nothing` otherwise
"""
function fit_irt_refinement_model(sequences::Vector{String},
                                  irt_predicted::Vector{Float32},
                                  irt_observed::Vector{Float32};
                                  ms_file_idx::Int=1,
                                  min_psms::Int=20,
                                  train_fraction::Float64=0.67)
    n = length(sequences)

    # Check minimum data requirement (21 features: 20 AAs + irt_predicted)
    if n < min_psms
        @debug_l1 "File $ms_file_idx: Insufficient PSMs ($n < $min_psms), skipping iRT refinement"
        return nothing
    end

    # Calculate errors (library - observed; positive when library overestimates)
    irt_errors = irt_predicted .- irt_observed

    # Prepare features
    features_df = prepare_features_dataframe(sequences, irt_predicted, irt_errors)

    # Train/validation split
    n_train = round(Int, n * train_fraction)
    n_val = n - n_train

    if n_val < 5  # Need reasonable validation set
        @debug_l1 "File $ms_file_idx: Insufficient validation data ($n_val < 5), skipping iRT refinement"
        return nothing
    end

    # Random shuffle indices
    indices = randperm(n)
    train_idx = indices[1:n_train]
    val_idx = indices[(n_train+1):end]

    train_df = features_df[train_idx, :]
    val_df = features_df[val_idx, :]

    # Build formula (all 20 amino acids + irt_predicted)
    aa_terms = [Symbol("count_", aa) for aa in STANDARD_AAS]
    formula_str = "error ~ irt_predicted + " * join(string.(aa_terms), " + ")
    formula = @eval @formula($(Meta.parse(formula_str)))

    # Train model on training set
    model = lm(formula, train_df)

    # Extract coefficients
    coef_values = coef(model)

    intercept = Float32(coef_values[1])
    irt_coef = Float32(coef_values[2])  # irt_predicted is always second

    # Extract AA coefficients in STANDARD_AAS order (indices 3-22)
    aa_coefficients = Float32.(coef_values[3:end])

    # Calculate training R²
    r2_train = Float32(r2(model))

    # Predict on validation set manually using coefficients
    # Build design matrix for validation set
    val_matrix = hcat(ones(n_val), val_df[:, :irt_predicted], Matrix(val_df[:, aa_terms]))
    val_predictions = val_matrix * coef_values

    # Calculate validation R²
    val_errors = val_df.error
    ss_res = sum((val_errors .- val_predictions).^2)
    ss_tot = sum((val_errors .- mean(val_errors)).^2)
    r2_val = Float32(1 - ss_res / ss_tot)

    # Calculate MAEs
    mae_original = Float32(mean(abs.(val_errors)))
    mae_refined = Float32(mean(abs.(val_errors .- val_predictions)))

    # Decision: use refinement if it improves validation MAE
    use_refinement = mae_refined < mae_original

    # Calculate improvement metrics
    mae_improvement = mae_original - mae_refined
    mae_improvement_pct = (mae_improvement / mae_original) * 100.0

    if use_refinement
        @user_info "File $ms_file_idx iRT Refinement ENABLED:" *
                  " Training R²=$(round(r2_train, digits=4))," *
                  " Validation R²=$(round(r2_val, digits=4))," *
                  " MAE: $(round(mae_original, digits=4)) → $(round(mae_refined, digits=4))" *
                  " (Δ=$(round(mae_improvement, digits=4)), $(round(mae_improvement_pct, digits=2))% improvement)\n"
    else
        @user_info "File $ms_file_idx iRT Refinement DISABLED:" *
                  " MAE: $(round(mae_original, digits=4)) → $(round(mae_refined, digits=4))" *
                  " (Δ=$(round(mae_improvement, digits=4)), $(round(mae_improvement_pct, digits=2))%)" *
                  " - No improvement, using library iRT\n"
    end

    # If refinement improves MAE, retrain on FULL dataset for final model
    if use_refinement
        @debug_l1 "File $ms_file_idx: Retraining on full dataset (train + validation combined)"

        # Retrain on combined train + validation data
        final_model = lm(formula, features_df)

        # Extract coefficients from retrained model
        final_coef_values = coef(final_model)

        final_intercept = Float32(final_coef_values[1])
        final_irt_coef = Float32(final_coef_values[2])

        # Create Dict mapping AA character to weight for zero-allocation lookup
        aa_weights = Dict{Char, Float32}()
        for (i, aa) in enumerate(STANDARD_AAS)
            aa_weights[aa] = Float32(final_coef_values[2 + i])  # Indices 3-22 in coef vector
        end

        return IrtRefinementModel(
            true,
            aa_weights,      # Dict instead of Vector
            final_intercept,
            final_irt_coef,
            mae_original,
            mae_refined,
            r2_train,
            r2_val
        )
    else
        # Return nothing if refinement doesn't help
        return nothing
    end
end
