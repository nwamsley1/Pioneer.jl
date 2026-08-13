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
    coefficients::Vector{Float32}   # indexed by token id (1..N_IRT_TOKENS); 0 if absent from the fit
end

struct MainSearchIrtRefinement{P<:LibraryPrecursors, S, M}
    precursors::P
    sequences::S          # getSequence(precursors), fetched once (vs per-precursor)
    structural_mods::M    # getStructuralMods(precursors), fetched once
    q_value_threshold::Float32
    min_precursors::Int
end

function MainSearchIrtRefinement(
    precursors::LibraryPrecursors;
    q_value_threshold::Float32 = PRESCORE_QVALUE_THRESHOLD,
    min_precursors::Int = 250,
)
    return MainSearchIrtRefinement(
        precursors,
        getSequence(precursors),
        getStructuralMods(precursors),
        q_value_threshold,
        min_precursors,
    )
end

@inline function _irt_pred_basis(current_irt_pred::Float32)
    x = Float64(current_irt_pred)
    return x, x * x
end

# ============================================================================
# Hardcoded peptide-composition token table for iRT refinement.
#
# `sequence` holds the 20 canonical amino acids; the only modifications the
# library admits are carbamidomethyl on C (Unimod:4) and oxidation on M
# (Unimod:35). The entire token universe is therefore a fixed set of 66 ids, so
# token counting is pure integer arithmetic — no per-precursor Dict, no token
# strings (previously the single largest allocation cluster in SearchDIA). The
# id scheme was verified to reproduce the previous string-keyed counts for every
# precursor in the benchmark library.
#
# Layout (k = residue-token index, 1..22):
#   1..20  plain residue, alphabetical: A C D E F G H I K L M N P Q R S T V W Y
#   21     C|Unimod:4          22     M|Unimod:35
#   NTERM token id = 22 + k         CTERM token id = 44 + k
#
# Column order in the design matrix is the ascending id of tokens present in a
# fold — fixed by this table, not by Dict hash order, so the fit is
# deterministic (no hash-layout sensitivity). If the library ever admits new
# residues/mods, extend _IRT_AA_K / _residue_token_k and N_IRT_TOKENS;
# unrecognized tokens are dropped (and warned once) rather than misbinned.
# ============================================================================
const _IRT_AA_LETTERS = collect("ACDEFGHIKLMNPQRSTVWY")
const _IRT_AA_K = let v = zeros(Int, 128)
    for (i, c) in enumerate(_IRT_AA_LETTERS)
        v[Int(c) + 1] = i
    end
    v
end
const _IRT_K_MODC = 21          # C|Unimod:4
const _IRT_K_MODM = 22          # M|Unimod:35
const N_IRT_TOKENS = 66         # 22 residue tokens + 22 NTERM + 22 CTERM

# Reusable per-thread scratch: token counts indexed by id, the list of touched
# ids (so reset/iteration is O(peptide length), not O(N_IRT_TOKENS)), and a
# residue-position -> mod-code buffer filled by the structural_mods parser.
struct IrtCountScratch
    counts::Vector{Float32}
    touched::Vector{Int}
    mod_at_pos::Vector{Int8}    # 0 none, 1 Unimod:4, 2 Unimod:35, -1 unrecognized
end
IrtCountScratch() = IrtCountScratch(zeros(Float32, N_IRT_TOKENS), Int[], Int8[])

@inline function _reset!(s::IrtCountScratch)
    @inbounds for id in s.touched
        s.counts[id] = 0.0f0
    end
    empty!(s.touched)
    return s
end

@inline function _add!(s::IrtCountScratch, id::Int)
    (id < 1 || id > N_IRT_TOKENS) && return
    @inbounds begin
        s.counts[id] == 0.0f0 && push!(s.touched, id)
        s.counts[id] += 1.0f0
    end
    return
end

# Compare the codeunit range str[lo:hi] against a literal, allocation-free.
@inline function _name_is(str, lo::Int, hi::Int, lit::String)
    ncodeunits(lit) == (hi - lo + 1) || return false
    @inbounds for k in 1:ncodeunits(lit)
        codeunit(str, lo + k - 1) != codeunit(lit, k) && return false
    end
    return true
end

# Mod-name codeunit range -> code (1 Unimod:4, 2 Unimod:35, -1 otherwise).
@inline function _mod_name_code(str, lo::Int, hi::Int)
    _name_is(str, lo, hi, "Unimod:4")  && return Int8(1)
    _name_is(str, lo, hi, "Unimod:35") && return Int8(2)
    return Int8(-1)
end

# Parse `structural_mods` ("(pos,site,modname)...") into scratch.mod_at_pos[1:L]
# by scanning codeunits — no regex, no per-mod String/Dict allocation. Terminal
# (n/c) mods are ignored (warned). A second mod at one position -> -1 (matches the
# old multi-mod path: no combined token exists, so it is treated as unrecognized).
function _parse_residue_mods!(s::IrtCountScratch, structural_mods, L::Int)
    length(s.mod_at_pos) < L && resize!(s.mod_at_pos, L)
    @inbounds for i in 1:L
        s.mod_at_pos[i] = Int8(0)
    end
    (ismissing(structural_mods) || isempty(structural_mods)) && return
    str = structural_mods
    n = ncodeunits(str)
    i = 1
    @inbounds while i <= n
        if codeunit(str, i) != UInt8('(')
            i += 1
            continue
        end
        i += 1
        pos = 0
        while i <= n && (b = codeunit(str, i)) >= UInt8('0') && b <= UInt8('9')
            pos = pos * 10 + Int(b - UInt8('0'))
            i += 1
        end
        (i <= n && codeunit(str, i) == UInt8(',')) || break
        i += 1
        site = codeunit(str, i); i += 1
        (i <= n && codeunit(str, i) == UInt8(',')) || break
        i += 1
        mstart = i
        while i <= n && codeunit(str, i) != UInt8(')')
            i += 1
        end
        mend = i - 1
        i += 1   # past ')'
        if site == UInt8('n') || site == UInt8('c')
            @debug_l2 "iRT refinement: unhandled terminal modification (ignored)"
        elseif 1 <= pos <= L
            code = _mod_name_code(str, mstart, mend)
            s.mod_at_pos[pos] = (s.mod_at_pos[pos] == 0) ? code : Int8(-1)
        end
    end
    return
end

# Residue-token index k (1..22), or 0 if the residue/modification is unrecognized.
@inline function _residue_token_k(aa::Char, code::Integer)
    if code == 0
        b = UInt32(aa)
        return b <= 0x7f ? @inbounds(_IRT_AA_K[b + 1]) : 0
    elseif code == 1
        return aa == 'C' ? _IRT_K_MODC : 0
    elseif code == 2
        return aa == 'M' ? _IRT_K_MODM : 0
    end
    return 0
end

# Fill `scratch` with token-id counts for one peptide; returns `scratch`.
function count_token_ids!(
    scratch::IrtCountScratch,
    sequence::AbstractString,
    structural_mods::Union{Missing, AbstractString},
)
    _reset!(scratch)
    L = length(sequence)
    _parse_residue_mods!(scratch, structural_mods, L)

    first_k = 0
    last_k = 0
    have_any = false
    for (position, aa) in enumerate(sequence)
        code = position <= length(scratch.mod_at_pos) ? scratch.mod_at_pos[position] : Int8(0)
        k = _residue_token_k(aa, code)
        if k == 0
            @debug_l2 "iRT refinement: unhandled residue token for '$(aa)' (ignored)"
            continue
        end
        _add!(scratch, k)
        have_any || (first_k = k; have_any = true)
        last_k = k
    end
    if have_any
        _add!(scratch, 22 + first_k)   # NTERM|<first residue token>
        _add!(scratch, 44 + last_k)    # CTERM|<last residue token>
    end
    return scratch
end

# Fill `scratch` for a library precursor (columns cached on the strategy).
@inline function count_token_ids!(
    scratch::IrtCountScratch,
    strategy::MainSearchIrtRefinement,
    precursor_idx::Integer,
)
    pid = UInt32(precursor_idx)
    return count_token_ids!(scratch, strategy.sequences[pid], strategy.structural_mods[pid])
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

    scratch = IrtCountScratch()

    # Pass 1: which token ids occur in this fold. Column order is the ascending
    # id (capacity-independent, deterministic) of the tokens that actually
    # appear, so X has no all-zero columns and the solve stays well-posed.
    present = falses(N_IRT_TOKENS)
    for pid in precursor_ids
        count_token_ids!(scratch, strategy, pid)
        for id in scratch.touched
            present[id] = true
        end
    end
    present_ids = findall(present)
    isempty(present_ids) && return nothing
    col_of = zeros(Int, N_IRT_TOKENS)
    for (j, id) in enumerate(present_ids)
        col_of[id] = j
    end

    # Pass 2: build the design matrix (counting is allocation-free, so recomputing
    # per precursor is cheap relative to the least-squares solve).
    use_quadratic_basis = n >= 3
    irt_basis_cols = use_quadratic_basis ? 2 : 1
    X = zeros(Float64, n, length(present_ids) + 1 + irt_basis_cols)
    X[:, 1] .= 1.0
    for i in eachindex(precursor_ids)
        linear_term, quadratic_term = _irt_pred_basis(irt_pred_inputs[i])
        X[i, 2] = linear_term
        if use_quadratic_basis
            X[i, 3] = quadratic_term
        end
        count_token_ids!(scratch, strategy, precursor_ids[i])
        for id in scratch.touched
            X[i, col_of[id] + 1 + irt_basis_cols] = Float64(scratch.counts[id])
        end
    end

    y = Float64.(irt_corrections)
    beta = try
        X \ y
    catch
        return nothing
    end

    # Dense coefficient vector indexed by token id; 0 for ids absent from this fold.
    coefficients = zeros(Float32, N_IRT_TOKENS)
    for (j, id) in enumerate(present_ids)
        coefficients[id] = Float32(beta[j + 1 + irt_basis_cols])
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
    scratch::IrtCountScratch,
    current_irt_pred::Float32,
)
    linear_term, quadratic_term = _irt_pred_basis(current_irt_pred)
    prediction = model.intercept +
                 model.irt_pred_coefficient * Float32(linear_term) +
                 model.irt_pred_squared_coefficient * Float32(quadratic_term)
    coef = model.coefficients
    @inbounds for id in scratch.touched
        prediction += coef[id] * scratch.counts[id]
    end
    return prediction
end

function predict_irt_refinement(
    scratch::IrtCountScratch,
    strategy::MainSearchIrtRefinement,
    model::MainSearchIrtCorrectionModel,
    precursor_idx::Integer,
    current_irt_pred::Float32,
)
    count_token_ids!(scratch, strategy, precursor_idx)
    return predict_irt_refinement(model, scratch, current_irt_pred)
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
    scratches = [IrtCountScratch() for _ in 1:nt]   # one reusable buffer per task
    @sync for t in 1:nt
        c_start = (t - 1) * chunk + 1
        c_start > n && break
        c_end = min(t * chunk, n)
        scratch = scratches[t]
        Threads.@spawn begin
            last_pid = UInt32(0)
            last_corr = 0f0
            have_corr = false
            @inbounds for row_idx in c_start:c_end
                pid = UInt32(prec_idx[row_idx])
                if pid != last_pid || !have_corr
                    corr = predict_irt_refinement(scratch, strategy, model, pid, current_pred[row_idx])
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
    scratches = [IrtCountScratch() for _ in 1:nt]   # one reusable buffer per task
    t_loop = time()
    @sync for t in 1:nt
        c_start = (t - 1) * chunk + 1
        c_start > n && break
        c_end = min(t * chunk, n)
        scratch = scratches[t]
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
                        corr = predict_irt_refinement(scratch, strategy, model, pid, current_pred[row_idx])
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
