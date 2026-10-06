# Ion-mobility prediction residuals, in library 1/K0 units. The observed
# mobility is never corrected: donor/receiver comparisons need a common scale.
struct ImCorrectionModel
    charges::Vector{Int}
    scale::Vector{Float64}
    beta::Vector{Float64}
end

"""
    im_sequence_fold(sequence; n_folds=5)

Stable sequence-only fold assignment. All charge states, modifications and
replicates of a peptide stay together; assignments do not depend on Julia's hash seed.
"""
function im_sequence_fold(sequence::AbstractString; n_folds::Int = 5)
    n_folds > 1 || throw(ArgumentError("n_folds must exceed one"))
    h = UInt64(0xcbf29ce484222325)
    for b in codeunits(sequence)
        h = (h ⊻ UInt64(b)) * UInt64(0x100000001b3)
    end
    return Int(mod(h, UInt64(n_folds))) + 1
end

function _im_design(pred, charge, tokens, charges)
    nz = length(charges)
    nt = size(tokens, 2)
    X = zeros(Float64, length(pred), (2 + nt) * nz)
    for i in eachindex(pred)
        k = findfirst(==(Int(charge[i])), charges)
        k === nothing && continue
        X[i, 2k - 1] = 1
        X[i, 2k] = pred[i]
        lo = 2nz + (k - 1) * nt + 1
        X[i, lo:(lo + nt - 1)] .= view(tokens, i, :)
    end
    return X
end

"""
    fit_im_correction(pred, obs, charge, tokens; min_charge=100, ridge=10.0)

Fit `obs - pred` using charge-specific offset/slopes and
residue/modification/terminal counts and squared peptide length. Ridge regularization
acts on standardized composition features. Unsupported charges remain unchanged.
Inputs must be finite, distinct high-confidence precursor anchors.
"""
function fit_im_correction(pred, obs, charge, tokens;
                           min_charge::Int = 100,
                           ridge::Float64 = 10.0)
    ridge >= 0 || throw(ArgumentError("ridge must be nonnegative"))
    n = length(pred)
    length(obs) == length(charge) == size(tokens, 1) == n || throw(DimensionMismatch("IM anchors"))
    all(isfinite, pred) && all(isfinite, obs) && all(isfinite, tokens) ||
        throw(ArgumentError("IM anchors must be finite"))
    charges = sort([Int(z) for z in unique(charge) if count(==(z), charge) >= min_charge])
    isempty(charges) && return nothing
    rows = findall(z -> Int(z) in charges, charge)
    X = _im_design(pred[rows], charge[rows], tokens[rows, :], charges)
    # Scale without centering to retain charge indicator intercepts. Drop
    # unsupported/constant composition terms through their ridge penalty.
    scale = [max(sqrt(sum(abs2, view(X, :, j)) / length(rows)), 1e-8) for j in axes(X, 2)]
    X ./= transpose(scale)
    penalty = fill(ridge, size(X, 2))
    penalty[1:(2length(charges))] .= 1e-6
    beta = (X' * X + Diagonal(penalty)) \ (X' * (Float64.(obs[rows]) .- pred[rows]))
    all(isfinite, beta) || return nothing
    return ImCorrectionModel(charges, scale, beta)
end

function predict_im_correction(model::ImCorrectionModel, pred, charge, tokens)
    X = _im_design(pred, charge, tokens, model.charges)
    return (X ./ transpose(model.scale)) * model.beta
end

"""
    crossfit_im_correction(pred, obs, charge, tokens, sequences; min_charge=100)

Return predictions corrected with charge-specific composition models on five
sequence folds. Each fold trains only on other peptide sequences. Charges with
insufficient training anchors retain the original predictions.
"""
function crossfit_im_correction(pred, obs, charge, tokens, sequences;
                                min_charge::Int = 100,
                                training_mask = trues(length(pred)), n_folds::Int = 5)
    n = length(pred)
    length(obs) == length(charge) == length(sequences) == length(training_mask) == size(tokens, 1) == n ||
        throw(DimensionMismatch("IM crossfit inputs"))
    folds = im_sequence_fold.(sequences; n_folds = n_folds)
    refined = Float64.(pred)
    valid = training_mask .& isfinite.(pred) .& isfinite.(obs)
    for fold in 1:n_folds
        train = findall(valid .& (folds .!= fold))
        test = findall((folds .== fold) .& isfinite.(pred))
        isempty(train) && continue
        model = fit_im_correction(pred[train], obs[train], charge[train], tokens[train, :];
                                  min_charge = min_charge)
        model === nothing && continue
        correction = predict_im_correction(model, pred[test], charge[test], tokens[test, :])
        for (j, row) in enumerate(test)
            p = pred[row] + correction[j]
            isfinite(p) && p > 0 && (refined[row] = p)
        end
    end
    return (; refined, folds)
end

function im_composition_tokens(sequences, modifications)
    tokens = zeros(Float64, length(sequences), N_IRT_TOKENS + 1)
    scratch = IrtCountScratch()
    for i in eachindex(sequences)
        # Residue counts already encode linear length; this adds curvature.
        tokens[i, end] = (length(sequences[i]) / 20)^2
        count_token_ids!(scratch, sequences[i], modifications[i])
        for id in scratch.touched
            tokens[i, id] = scratch.counts[id]
        end
    end
    return tokens
end

"""
    refine_im_error!(psms, precursors, sigma)

Correct library mobility predictions on sequence-held-out folds using target
anchors with first-pass probability > 0.9 and q-value <= 0.01. Preserve `im_obs`
and the existing scan calibration/gate. Both ordinary and partitioned searches
use this best-per-precursor table before experiment-wide and MBR scoring.
"""
function refine_im_error!(psms::DataFrame, precursors, sigma::Real)
    im_lib = getInvIonMobility(precursors)
    (im_lib === nothing || nrow(psms) == 0) && return nothing
    ids = psms.precursor_idx
    sequences = getSequence(precursors)[ids]
    modifications = getStructuralMods(precursors)[ids]
    pred = Float64.(im_lib[ids])
    tokens = im_composition_tokens(sequences, modifications)
    qvals = Vector{Float32}(undef, nrow(psms))
    get_qvalues!(psms.lgbm_prob, collect(Bool, psms.target), qvals; doSort = true)
    anchors = psms.target .& (psms.lgbm_prob .> 0.9f0) .& (qvals .<= 0.01f0)
    # Isotope trace variants must not weight one precursor multiple times.
    seen = Set{UInt32}()
    for i in sortperm(psms.lgbm_prob; rev = true)
        anchors[i] || continue
        if ids[i] in seen
            anchors[i] = false
        else
            push!(seen, ids[i])
        end
    end
    result = crossfit_im_correction(pred, psms.im_obs, psms.charge, tokens, sequences;
                                    training_mask = anchors)
    psms[!, :im_pred] = Float32.(pred)
    psms[!, :im_pred_refined] = Float32.(result.refined)
    psms[!, :im_error_uncorrected] = copy(psms.im_error)
    psms[!, :im_error] = Float32.((result.refined .- psms.im_obs) ./ max(Float64(sigma), 1e-6))
    @debug_l1 "IM prediction refinement: anchors=$(count(anchors))"
    return result
end
