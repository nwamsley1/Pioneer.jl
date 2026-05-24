# Predicted iRT refinement for MainSearch.
#
# MainSearch first trains a per-run classifier on the library iRT predictions.
# High-confidence target precursors from that pass provide observed residuals
# (irt_obs - irt_pred). This helper fits a small peptide-composition correction
# model out-of-fold and applies it to every candidate row before the per-run
# classifier is reapplied to the refreshed iRT-dependent features.

struct MainSearchIrtCorrectionModel
    intercept::Float32
    irt_pred_coefficient::Float32
    irt_pred_squared_coefficient::Float32
    coefficients::Dict{String, Float32}
end

struct MainSearchIrtRefinement{P<:LibraryPrecursors}
    precursors::P
    q_value_threshold::Float32
    min_precursors::Int
end

function MainSearchIrtRefinement(
    precursors::LibraryPrecursors;
    q_value_threshold::Float32 = PRESCORE_QVALUE_THRESHOLD,
    min_precursors::Int = 250,
)
    return MainSearchIrtRefinement(precursors, q_value_threshold, min_precursors)
end

@inline function _irt_pred_basis(current_irt_pred::Float32)
    x = Float64(current_irt_pred)
    return x, x * x
end

function _precursor_token_counts(
    sequence::AbstractString,
    structural_mods::Union{Missing, AbstractString},
)
    counts = Dict{String, Float32}()
    position_mods = Dict{Int, Vector{String}}()
    residue_tokens = String[]

    if !ismissing(structural_mods) && !isempty(structural_mods)
        mod_regex = r"\((\d+),([A-Z]|[nc]),([^,\)]+)\)"
        for m in eachmatch(mod_regex, structural_mods)
            position = parse(Int, m.captures[1])
            site = only(m.captures[2])
            mod_name = m.captures[3]

            if site == 'n' || site == 'c'
                token = string(site, "|", mod_name)
                counts[token] = get(counts, token, 0.0f0) + 1.0f0
            elseif 1 <= position <= length(sequence)
                push!(get!(position_mods, position, String[]), mod_name)
            end
        end
    end

    for (position, aa) in enumerate(sequence)
        mods_here = get(position_mods, position, nothing)
        token = if isnothing(mods_here) || isempty(mods_here)
            string(aa)
        else
            string(aa, "|", join(sort(mods_here), "&"))
        end
        push!(residue_tokens, token)
        counts[token] = get(counts, token, 0.0f0) + 1.0f0
    end

    if !isempty(residue_tokens)
        counts["NTERM|" * first(residue_tokens)] = get(counts, "NTERM|" * first(residue_tokens), 0.0f0) + 1.0f0
        counts["CTERM|" * last(residue_tokens)] = get(counts, "CTERM|" * last(residue_tokens), 0.0f0) + 1.0f0
    end

    return counts
end

function precursor_token_counts(
    strategy::MainSearchIrtRefinement,
    precursor_idx::Integer,
)
    pid = UInt32(precursor_idx)
    sequence = String(getSequence(strategy.precursors)[pid])
    structural_mods = getStructuralMods(strategy.precursors)[pid]
    return _precursor_token_counts(sequence, structural_mods)
end

function _unique_sorted_precursor_ids(precursor_idx::AbstractVector)
    precursor_ids = unique(UInt32.(precursor_idx))
    sort!(precursor_ids)
    return precursor_ids
end

function _passing_precursor_targets(
    precursor_idx::AbstractVector,
    targets::AbstractVector{Bool},
    trace_prob::AbstractVector,
    irt_pred::AbstractVector,
    irt_obs::AbstractVector,
    q_value_threshold::Float32,
)
    n = length(trace_prob)
    n == 0 && return UInt32[], Float32[], Float32[]

    q_values = Vector{Float32}(undef, n)
    get_qvalues!(Float32.(trace_prob), collect(Bool, targets), q_values; doSort = true)
    pass_mask = collect(Bool, targets) .& (q_values .<= q_value_threshold)

    passing_precursor_ids = _unique_sorted_precursor_ids(precursor_idx[pass_mask])
    sum_pred = Dict{UInt32, Float64}()
    sum_irt = Dict{UInt32, Float64}()
    count_irt = Dict{UInt32, Int}()
    for i in eachindex(pass_mask)
        pass_mask[i] || continue
        pid = UInt32(precursor_idx[i])
        sum_pred[pid] = get(sum_pred, pid, 0.0) + Float64(irt_pred[i])
        sum_irt[pid] = get(sum_irt, pid, 0.0) + Float64(irt_obs[i])
        count_irt[pid] = get(count_irt, pid, 0) + 1
    end

    irt_pred_inputs = Float32[sum_pred[pid] / count_irt[pid] for pid in passing_precursor_ids]
    irt_corrections = Float32[
        (sum_irt[pid] / count_irt[pid]) - (sum_pred[pid] / count_irt[pid])
        for pid in passing_precursor_ids
    ]
    return passing_precursor_ids, irt_pred_inputs, irt_corrections
end

function fit_irt_refinement_model(
    strategy::MainSearchIrtRefinement,
    precursor_ids::Vector{UInt32},
    irt_pred_inputs::Vector{Float32},
    irt_corrections::Vector{Float32},
)
    n = length(precursor_ids)
    if n < strategy.min_precursors ||
       n != length(irt_pred_inputs) ||
       n != length(irt_corrections)
        return nothing
    end

    token_to_col = Dict{String, Int}()
    precursor_counts = Vector{Dict{String, Float32}}(undef, n)

    for i in eachindex(precursor_ids)
        counts = precursor_token_counts(strategy, precursor_ids[i])
        precursor_counts[i] = counts
        for token in keys(counts)
            if !haskey(token_to_col, token)
                token_to_col[token] = length(token_to_col) + 1
            end
        end
    end
    isempty(token_to_col) && return nothing

    use_quadratic_basis = n >= 3
    irt_basis_cols = use_quadratic_basis ? 2 : 1
    X = zeros(Float64, n, length(token_to_col) + 1 + irt_basis_cols)
    X[:, 1] .= 1.0

    for i in eachindex(precursor_counts)
        linear_term, quadratic_term = _irt_pred_basis(irt_pred_inputs[i])
        X[i, 2] = linear_term
        if use_quadratic_basis
            X[i, 3] = quadratic_term
        end
        for (token, count) in precursor_counts[i]
            X[i, token_to_col[token] + 1 + irt_basis_cols] = Float64(count)
        end
    end

    y = Float64.(irt_corrections)
    beta = try
        X \ y
    catch
        return nothing
    end

    coefficients = Dict{String, Float32}()
    for (token, col_idx) in token_to_col
        coefficients[token] = Float32(beta[col_idx + 1 + irt_basis_cols])
    end

    return MainSearchIrtCorrectionModel(
        Float32(beta[1]),
        Float32(beta[2]),
        use_quadratic_basis ? Float32(beta[3]) : 0.0f0,
        coefficients,
    )
end

function predict_irt_refinement(
    model::MainSearchIrtCorrectionModel,
    counts::Dict{String, Float32},
    current_irt_pred::Float32,
)
    linear_term, quadratic_term = _irt_pred_basis(current_irt_pred)
    prediction = model.intercept +
                 model.irt_pred_coefficient * Float32(linear_term) +
                 model.irt_pred_squared_coefficient * Float32(quadratic_term)
    for (token, count) in counts
        prediction += get(model.coefficients, token, 0.0f0) * count
    end
    return prediction
end

function predict_irt_refinement(
    strategy::MainSearchIrtRefinement,
    model::MainSearchIrtCorrectionModel,
    precursor_idx::Integer,
    current_irt_pred::Float32,
)
    return predict_irt_refinement(
        model,
        precursor_token_counts(strategy, precursor_idx),
        current_irt_pred,
    )
end

function _refresh_predicted_irt_dependent_features!(psms::DataFrame)
    if hasproperty(psms, :irt_error) && hasproperty(psms, :irt_obs) && hasproperty(psms, :irt_pred)
        psms[!, :irt_error] = abs.(Float32.(psms[!, :irt_obs]) .- Float32.(psms[!, :irt_pred]))
    end
    if hasproperty(psms, :irt_diff) && hasproperty(psms, :irt_obs) && hasproperty(psms, :irt_pred)
        psms[!, :irt_diff] = abs.(Float32.(psms[!, :irt_obs]) .- Float32.(psms[!, :irt_pred]))
    end
    return nothing
end

function apply_mainsearch_irt_refinement_model!(
    psms::DataFrame,
    strategy::MainSearchIrtRefinement,
    model::MainSearchIrtCorrectionModel,
)
    nrow(psms) == 0 && return nothing
    !(hasproperty(psms, :irt_pred) && hasproperty(psms, :precursor_idx)) && return nothing

    current_pred = Float32.(psms[!, :irt_pred])
    refined_pred = copy(current_pred)
    prec_idx = psms.precursor_idx
    # PSMs are sorted contiguous-by-precursor_idx upstream (see
    # `permute_psms_by_precursor_idx!`). current_pred is the library iRT,
    # constant per precursor. Each chunk walks its row range with a last-pid
    # sentinel — correction is recomputed only when pid changes. Chunks write
    # to disjoint rows so the threading is allocation-free.
    n = nrow(psms)
    nt = Threads.nthreads()
    chunk = max(1, cld(n, nt))
    @sync for t in 1:nt
        c_start = (t - 1) * chunk + 1
        c_start > n && break
        c_end = min(t * chunk, n)
        Threads.@spawn begin
            last_pid = UInt32(0)
            last_corr = 0f0
            have_corr = false
            @inbounds for row_idx in c_start:c_end
                pid = UInt32(prec_idx[row_idx])
                if pid != last_pid || !have_corr
                    corr = predict_irt_refinement(strategy, model, pid, current_pred[row_idx])
                    last_corr = isfinite(corr) ? Float32(corr) : 0f0
                    last_pid = pid
                    have_corr = true
                end
                refined_pred[row_idx] = current_pred[row_idx] + last_corr
            end
        end
    end
    psms[!, :irt_pred] = refined_pred
    _refresh_predicted_irt_dependent_features!(psms)
    return nothing
end

function apply_mainsearch_irt_refinement_model!(
    psms::DataFrame,
    strategy::MainSearchIrtRefinement,
    models::Dict{UInt8, MainSearchIrtCorrectionModel},
)
    nrow(psms) == 0 && return nothing
    !(hasproperty(psms, :irt_pred) && hasproperty(psms, :precursor_idx) && hasproperty(psms, :cv_fold)) &&
        return nothing
    isempty(models) && return nothing

    current_pred = Float32.(psms[!, :irt_pred])
    refined_pred = copy(current_pred)
    folds = UInt8.(psms[!, :cv_fold])
    prec_idx = psms.precursor_idx
    # cv_fold is precursor-keyed; PSMs are sorted contiguous-by-precursor_idx
    # upstream. Each chunk walks with a last-pid sentinel — predict_irt_refinement
    # runs once per pid change (vs once per PSM previously). Chunks write disjoint
    # rows so threading is allocation-free.
    n = nrow(psms)
    nt = Threads.nthreads()
    chunk = max(1, cld(n, nt))
    t_loop = time()
    @sync for t in 1:nt
        c_start = (t - 1) * chunk + 1
        c_start > n && break
        c_end = min(t * chunk, n)
        Threads.@spawn begin
            last_pid = UInt32(0)
            last_corr = 0f0
            have_corr = false
            @inbounds for row_idx in c_start:c_end
                pid = UInt32(prec_idx[row_idx])
                if pid != last_pid || !have_corr
                    model = get(models, folds[row_idx], nothing)
                    if isnothing(model)
                        last_corr = 0f0
                    else
                        corr = predict_irt_refinement(strategy, model, pid, current_pred[row_idx])
                        last_corr = isfinite(corr) ? Float32(corr) : 0f0
                    end
                    last_pid = pid
                    have_corr = true
                end
                refined_pred[row_idx] = current_pred[row_idx] + last_corr
            end
        end
    end
    t_predict_loop = time() - t_loop
    psms[!, :irt_pred] = refined_pred
    t_r = time()
    _refresh_predicted_irt_dependent_features!(psms)
    t_refresh = time() - t_r
    @debug_l1 "    apply (cv_fold dict, n=$(nrow(psms))): predict_loop=$(round(t_predict_loop, digits=2))s refresh=$(round(t_refresh, digits=2))s"
    return nothing
end

function refine_mainsearch_irt_predictions!(
    psms::DataFrame,
    best_psms::DataFrame,
    scores::Vector{Float32},
    strategy::MainSearchIrtRefinement,
)
    required = (:precursor_idx, :target, :irt_pred, :irt_obs)
    if !all(c -> hasproperty(best_psms, c), required) ||
       !all(c -> hasproperty(psms, c), (:precursor_idx, :irt_pred)) ||
       length(scores) != nrow(best_psms)
        return (refined = false, training_target_precursors = UInt32[], model = nothing)
    end

    if hasproperty(best_psms, :cv_fold) && hasproperty(psms, :cv_fold)
        best_folds = UInt8.(best_psms[!, :cv_fold])
        fold_values = sort!(unique(best_folds))
        models = Dict{UInt8, MainSearchIrtCorrectionModel}()
        training_target_precursors = UInt32[]

        t_qval = 0.0
        t_fit  = 0.0
        for fold in fold_values
            train_rows = findall(!=(fold), best_folds)
            isempty(train_rows) && continue

            tq = time()
            precursor_ids, irt_pred_inputs, irt_corrections = _passing_precursor_targets(
                best_psms.precursor_idx[train_rows],
                best_psms.target[train_rows],
                scores[train_rows],
                best_psms.irt_pred[train_rows],
                best_psms.irt_obs[train_rows],
                strategy.q_value_threshold,
            )
            t_qval += time() - tq
            append!(training_target_precursors, precursor_ids)

            tf = time()
            model = fit_irt_refinement_model(
                strategy,
                precursor_ids,
                irt_pred_inputs,
                irt_corrections,
            )
            t_fit += time() - tf
            isnothing(model) && continue
            models[fold] = model
        end

        training_target_precursors = sort!(unique(training_target_precursors))
        if isempty(models)
            return (
                refined = false,
                training_target_precursors = training_target_precursors,
                model = nothing,
            )
        end

        t_ap = time()
        apply_mainsearch_irt_refinement_model!(psms, strategy, models)
        t_apply_psms = time() - t_ap
        t_ab = time()
        apply_mainsearch_irt_refinement_model!(best_psms, strategy, models)
        t_apply_best = time() - t_ab
        tokens_str = join([length(m.coefficients) for m in values(models)], ",")
        @debug_l1 "  iRT refine breakdown: qval+agg=$(round(t_qval, digits=2))s " *
                   "fit=$(round(t_fit, digits=2))s " *
                   "apply_psms=$(round(t_apply_psms, digits=2))s " *
                   "apply_best=$(round(t_apply_best, digits=2))s " *
                   "(n_psms=$(nrow(psms))  n_best=$(nrow(best_psms))  " *
                   "n_train_total=$(length(training_target_precursors))  " *
                   "n_tokens_per_fold=$(tokens_str))"

        return (
            refined = true,
            training_target_precursors = training_target_precursors,
            model = models,
        )
    end

    precursor_ids, irt_pred_inputs, irt_corrections = _passing_precursor_targets(
        best_psms.precursor_idx,
        best_psms.target,
        scores,
        best_psms.irt_pred,
        best_psms.irt_obs,
        strategy.q_value_threshold,
    )
    model = fit_irt_refinement_model(strategy, precursor_ids, irt_pred_inputs, irt_corrections)
    if isnothing(model)
        return (refined = false, training_target_precursors = precursor_ids, model = nothing)
    end

    apply_mainsearch_irt_refinement_model!(psms, strategy, model)
    apply_mainsearch_irt_refinement_model!(best_psms, strategy, model)

    return (refined = true, training_target_precursors = precursor_ids, model = model)
end
