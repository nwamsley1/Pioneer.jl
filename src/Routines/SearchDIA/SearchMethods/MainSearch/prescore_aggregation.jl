#==========================================================
Prescore Aggregation — Raw Log-Odds Combination

Combines per-file LightGBM prescore probabilities into a
single global score per precursor via log-odds averaging.

Uses raw LightGBM probabilities directly (no PEP calibration).
Simulation showed isotonic regression PEP calibration introduces
systematic anti-conservative FDR bias (mean inflation 1.29x at
1% FDR across 100 seeds) because each sample's own target/decoy
label influences its own PEP estimate. Raw log-odds of CV-scored
LightGBM probabilities is perfectly calibrated (mean inflation
0.999).
==========================================================#

"""
    _logodds_combine(probs, top_n, floor) -> Float32

Take top-n values from `probs`, clamp to [floor, 1-ε],
convert to log-odds, average, and convert back to probability.

Optimized to avoid per-call allocations:
- Fast path for single-element vectors (scalar ops, no sort)
- Fast path for n==1 (find max, scalar ops)
- General case uses partialsort! (in-place) + scalar loop
"""
function _logodds_combine(probs::Vector{Float32}, top_n::Int, floor::Float32)::Float32
    isempty(probs) && return 0.0f0
    eps = 1f-6
    n = min(length(probs), top_n)
    if n == 1
        # Single value needed: find max without sorting
        p = length(probs) == 1 ? @inbounds(probs[1]) : maximum(probs)
        c = clamp(p, floor, 1 - eps)
        # logit → average(1 element) → sigmoid is identity
        return c
    end
    # In-place partial sort — safe because probs is a temporary consumed only once
    partialsort!(probs, 1:n; rev=true)
    lo_sum = 0.0f0
    @inbounds for i in 1:n
        c = clamp(probs[i], floor, 1 - eps)
        lo_sum += log(c / (1 - c))
    end
    avg = lo_sum / n
    return 1.0f0 / (1 + exp(-avg))
end

"""
Typed inner function for the per-file aggregation loop in `aggregate_prescore_globally!`.
Accepts Arrow columns as `AbstractVector` so the compiler specializes on concrete column
types at each call site, eliminating dynamic dispatch inside the hot loop.
"""
function _aggregate_file_scores!(
        prec_idx_col::AbstractVector{UInt32},
        calibrated_probs::AbstractVector{Float32},
        n_above_hm_col::AbstractVector{UInt16},
        weight_col::AbstractVector{Float32},
        rt_fwhm_col::AbstractVector{Float32},
        best_rt_col::AbstractVector{Float32},
        prec_probs_by_run::Dictionary{UInt32, Vector{Float32}},
        prec_best_prob::Dictionary{UInt32, Float32},
        prec_best_n_above_hm::Dictionary{UInt32, UInt16},
        prec_best_weight::Dictionary{UInt32, Float32},
        prec_best_rt_fwhm::Dictionary{UInt32, Float32},
        prec_best_rt::Dictionary{UInt32, Float32})

    @inbounds for i in eachindex(prec_idx_col)
        pid = prec_idx_col[i]
        p = calibrated_probs[i]
        if haskey(prec_probs_by_run, pid)
            push!(prec_probs_by_run[pid], p)
            if p > prec_best_prob[pid]
                prec_best_prob[pid] = p
                prec_best_n_above_hm[pid] = n_above_hm_col[i]
                prec_best_weight[pid] = weight_col[i]
                prec_best_rt_fwhm[pid] = rt_fwhm_col[i]
                prec_best_rt[pid] = best_rt_col[i]
            end
        else
            insert!(prec_probs_by_run, pid, Float32[p])
            insert!(prec_best_prob, pid, p)
            insert!(prec_best_n_above_hm, pid, n_above_hm_col[i])
            insert!(prec_best_weight, pid, weight_col[i])
            insert!(prec_best_rt_fwhm, pid, rt_fwhm_col[i])
            insert!(prec_best_rt, pid, best_rt_col[i])
        end
    end
    return nothing
end

# Hardcoded prescore q-value threshold. The active filter cutoff between
# MainSearch's per-file LightGBM and ScoringSearch's experiment-wide LightGBM.
# Tuned via per-file ID sweep (see commit 5a24ccdf): Astral peaks at 1.5%,
# Eclipse near-tied with 1%; 0.015 is the cross-dataset compromise.
const PRESCORE_QVALUE_THRESHOLD = 0.015f0

# Floor for the per-file logit clamp inside _logodds_combine. Caps the
# negative contribution of a single bad file (worst-case single-file logit
# = log(floor / (1-floor))).
const LOGODDS_CLAMP_FLOOR = 0.02f0

"""
    aggregate_prescore_globally!(search_context; fold_suffix="") -> (rt_binned_tol = ...,)

Load per-file LightGBM prescore arrow files, aggregate raw LightGBM probabilities
via log-odds averaging of the top ⌊√N⌋ values per precursor, compute global
target/decoy q-values, and use the high-confidence (q ≤ 0.001) targets to fit an
RT-binned chromatographic tolerance. Returns the `RTBinnedTolerance` (or `nothing`
if insufficient data). The cross-file q-value is used only as a selector for the
RT-tolerance fit; it is **not** propagated as a filter to ScoringSearch.
"""
function aggregate_prescore_globally!(search_context::SearchContext;
                                      fold_suffix::String="")
    r(t) = round(t; digits=2)
    t_total_start = time()

    ms_data = getMSData(search_context)
    n_files = length(ms_data)
    prescore_dir = joinpath(getDataOutDir(search_context), "temp_data", "main_search_psms")

    # Collect per-file probs and metrics per precursor
    prec_probs_by_run = Dictionary{UInt32, Vector{Float32}}()
    prec_best_prob = Dictionary{UInt32, Float32}()
    prec_best_n_above_hm = Dictionary{UInt32, UInt16}()
    prec_best_weight = Dictionary{UInt32, Float32}()
    prec_best_rt_fwhm = Dictionary{UInt32, Float32}()
    prec_best_rt = Dictionary{UInt32, Float32}()
    n_valid_files = 0

    t_reads = 0.0
    t_loop = 0.0

    for ms_file_idx in 1:n_files
        file_name = getParsedFileName(ms_data, ms_file_idx)
        score_path = joinpath(prescore_dir, "$(file_name)$(fold_suffix).arrow")
        !isfile(score_path) && continue

        t_r = time()
        scores = Arrow.Table(score_path)
        t_reads += time() - t_r
        n_valid_files += 1

        t_l = time()
        _aggregate_file_scores!(
            scores[:precursor_idx], scores[:lgbm_prob],
            scores[:n_above_hm], scores[:weight],
            scores[:rt_fwhm], scores[:best_rt],
            prec_probs_by_run, prec_best_prob,
            prec_best_n_above_hm, prec_best_weight,
            prec_best_rt_fwhm, prec_best_rt
        )
        t_loop += time() - t_l
    end

    # Aggregate via log-odds; look up target/decoy from spectral library
    t_qval_start = time()
    is_decoy = getIsDecoy(getPrecursors(getSpecLib(search_context)))
    sqrt_n = max(1, floor(Int, sqrt(n_valid_files)))
    n_unique = length(prec_probs_by_run)
    global_prec_idxs = Vector{UInt32}(undef, n_unique)
    global_probs = Vector{Float32}(undef, n_unique)
    global_targets = Vector{Bool}(undef, n_unique)

    for (i, (pid, probs)) in enumerate(pairs(prec_probs_by_run))
        global_prec_idxs[i] = pid
        global_probs[i] = _logodds_combine(probs, sqrt_n, LOGODDS_CLAMP_FLOOR)
        global_targets[i] = !is_decoy[pid]
    end

    # Q-values on global scores — used only as a selector for the RT-tolerance fit
    # below. Not propagated to ScoringSearch.
    global_qvals = Vector{Float32}(undef, n_unique)
    get_qvalues!(global_probs, global_targets, global_qvals; doSort=true)

    # --- Compute per-RT-bin tolerances ---
    n_rt_bins = 20
    min_n_scans = 10  # floor: at least this many scans worth of RT half-width tolerance
    rt_qval_cutoff = 0.001f0

    # Collect (best_rt, rt_fwhm, weight, n_above_hm) for q<=0.1% targets
    rt_obs_hc = Float32[]
    rt_fwhm_hc = Float32[]
    rt_weight_hc = Float32[]
    rt_n_hm_hc = UInt16[]
    for i in eachindex(global_prec_idxs)
        if global_qvals[i] <= rt_qval_cutoff && global_targets[i]
            pid = global_prec_idxs[i]
            haskey(prec_best_rt, pid) || continue
            haskey(prec_best_rt_fwhm, pid) || continue
            haskey(prec_best_weight, pid) || continue
            haskey(prec_best_n_above_hm, pid) || continue
            push!(rt_obs_hc, prec_best_rt[pid])
            push!(rt_fwhm_hc, prec_best_rt_fwhm[pid])
            push!(rt_weight_hc, prec_best_weight[pid])
            push!(rt_n_hm_hc, prec_best_n_above_hm[pid])
        end
    end

    rt_binned_tol = nothing
    if length(rt_obs_hc) >= n_rt_bins * 5
        # Compute quantile bin edges
        quantiles = range(0f0, 1f0, length=n_rt_bins+1)
        rt_bin_edges = Float32[Float32(quantile(rt_obs_hc, q)) for q in quantiles]
        # Ensure edges are strictly increasing
        for i in 2:length(rt_bin_edges)
            if rt_bin_edges[i] <= rt_bin_edges[i-1]
                rt_bin_edges[i] = rt_bin_edges[i-1] + 0.001f0
            end
        end

        rt_tols = Vector{Float32}(undef, n_rt_bins)
        n_top_per_bin = 100

        for bin in 1:n_rt_bins
            lo = rt_bin_edges[bin]
            hi = rt_bin_edges[bin + 1]

            # Select precursors in this bin
            bin_mask = (rt_obs_hc .>= lo) .& (rt_obs_hc .< hi)
            if bin == n_rt_bins  # last bin includes right edge
                bin_mask .|= (rt_obs_hc .== hi)
            end

            bin_fwhm = rt_fwhm_hc[bin_mask]
            bin_weight = rt_weight_hc[bin_mask]
            bin_n_hm = rt_n_hm_hc[bin_mask]

            if length(bin_fwhm) >= 3
                # Take top-K by weight
                n_use = min(n_top_per_bin, length(bin_weight))
                top_idx = partialsortperm(bin_weight, 1:n_use; rev=true)
                selected_fwhm = bin_fwhm[top_idx]
                selected_n_hm = bin_n_hm[top_idx]

                bin_cycle_times = Float32[]
                for j in eachindex(selected_fwhm)
                    n_hm = selected_n_hm[j]
                    n_hm >= 2 || continue
                    ct = selected_fwhm[j] / Float32(n_hm - 1)
                    ct > 0.001f0 || continue
                    push!(bin_cycle_times, ct)
                end
                selected_fwhm_median = Float32(median(selected_fwhm))
                rt_cycle_time = !isempty(bin_cycle_times) ? Float32(quantile(bin_cycle_times, 0.75)) : selected_fwhm_median

                fwhm_half = 2.0f0 * selected_fwhm_median
                cycle_half = Float32(min_n_scans) * rt_cycle_time
                rt_tols[bin] = max(fwhm_half, cycle_half)
            else
                rt_tols[bin] = Float32(median(rt_fwhm_hc)) / 2.0f0
            end
        end

        rt_binned_tol = RTBinnedTolerance(rt_bin_edges, rt_tols)

        @debug_l1 "RT-binned tolerance ($n_rt_bins bins): RT $(round(rt_bin_edges[1], digits=1))-$(round(rt_bin_edges[end], digits=1)), tol +/-$(round(minimum(rt_tols), digits=3)) to +/-$(round(maximum(rt_tols), digits=3)) min"
    else
        @debug_l1 "Not enough precursors for RT-binned tolerance ($(length(rt_obs_hc)) < $(n_rt_bins * 5))"
    end

    t_qval = time() - t_qval_start
    t_total = time() - t_total_start
    @debug_l1 "Prescore aggregation: file_reads=$(r(t_reads))s, loop=$(r(t_loop))s, combination+qvalues=$(r(t_qval))s, total=$(r(t_total))s"

    return (rt_binned_tol = rt_binned_tol,)
end
