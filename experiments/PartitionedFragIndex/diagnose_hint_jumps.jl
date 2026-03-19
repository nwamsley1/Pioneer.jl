#!/usr/bin/env julia
#
# Diagnostic v2: 5-Da direct hint + lb=first_matching_bin approach
#
# Tests the new hint design:
#   hint[j] = k where getLow(frag_bins[j+k]) - getLow(frag_bins[j]) >= 5.0 Da
#   est_step = hint[lb] * (delta_mz / 5.0)
#   lb advances to first_matching_bin (provably safe by binary search guarantee)
#   Hint-based lb advancement before exponential search (provably safe for delta > 5 Da,
#     approximate with conservative factor for delta <= 5 Da)
#
# Usage:
#   julia --project=. experiments/PartitionedFragIndex/diagnose_hint_jumps.jl /path/to/params.json [n_scans]

using Pioneer
using Statistics

include(joinpath(@__DIR__, "partitioned_types.jl"))
include(joinpath(@__DIR__, "build_partitioned_index.jl"))
include(joinpath(@__DIR__, "search_partitioned_index.jl"))

# ── 5-Da direct hint computation ────────────────────────────────────────

"""
    _compute_skip_hints_5da(frag_bins, rt_bins) -> Vector{UInt16}

For each frag bin j within an RT bin, binary search for the smallest k where
`getLow(frag_bins[j+k]) >= getLow(frag_bins[j]) + 5.0f0`.
Store as UInt16 (dense regions may exceed 255). Clamped to RT bin range bounds.
"""
function _compute_skip_hints_5da(
    frag_bins::Vector{Pioneer.FragIndexBin{Float32}},
    rt_bins::Vector{Pioneer.FragIndexBin{Float32}},
)
    n_fb = length(frag_bins)
    hints = ones(UInt16, n_fb)  # default 1 (minimum skip)

    for rt_bin in rt_bins
        range = Pioneer.getSubBinRange(rt_bin)
        fb_start = Int(first(range))
        fb_end = Int(last(range))
        fb_start > fb_end && continue

        for j in fb_start:fb_end
            target_low = Pioneer.getLow(frag_bins[j]) + 5.0f0
            # Binary search for smallest k where getLow(frag_bins[j+k]) >= target_low
            lo_k = 1
            hi_k = fb_end - j
            hi_k <= 0 && continue  # at end of RT bin

            # If even the last bin doesn't reach +5 Da, set hint to max available
            if Pioneer.getLow(frag_bins[fb_end]) < target_low
                hints[j] = UInt16(fb_end - j)
                continue
            end

            # Binary search
            result_k = hi_k
            while lo_k <= hi_k
                mid_k = (lo_k + hi_k) >>> 1
                if Pioneer.getLow(frag_bins[j + mid_k]) >= target_low
                    result_k = mid_k
                    hi_k = mid_k - 1
                else
                    lo_k = mid_k + 1
                end
            end

            hints[j] = UInt16(clamp(result_k, 1, 65535))
        end
    end

    return hints
end

# ── Record type ──────────────────────────────────────────────────────────

struct JumpRecord2
    delta_mz::Float32           # frag_mz_min difference between consecutive peaks
    hint_5da::UInt16            # 5-Da hint value at lb
    est_step::Float32           # hint * (delta_mz / 5.0)
    per_peak_jump::UInt32       # first_match_current - first_match_previous
    new_lb::UInt32              # lb after hint-based advancement
    first_match::UInt32         # actual first matching bin (from findFirstFragmentBin)
    ub_guess::UInt32            # hint-based UB guess
    ub_final::UInt32            # final UB (after exponential doubling if needed)
    range_size::UInt32          # ub_final - new_lb + 1
    lb_overshoot::Bool          # did new_lb > first_match? (bad if true)
    ub_guess_sufficient::Bool   # was ub_guess >= frag_mz_max before exponential?
    bin_width::Float32          # getHigh - getLow for the first_match bin
    lb_advance_factor::Float32  # which factor was used (for sweeping)
end

# ── Instrumented query function ──────────────────────────────────────────

"""
    queryFragmentInstrumented_v2!(records, counter, ..., hints_5da, prev_first_match,
                                   lb_factor, ub_overshoot_factor)

New approach:
  1. Hint-based lb advancement (provably safe for delta > 5 Da)
  2. Hint-based UB guess
  3. Exponential doubling only if UB guess insufficient
  4. findFirstFragmentBin for exact first match
  5. Return (first_match, ub) so lb advances for next peak
"""
@inline function queryFragmentInstrumented_v2!(
        records::Vector{JumpRecord2},
        counter,
        frag_bin_max_idx::UInt32,
        lower_bound_guess::UInt32,
        upper_bound_guess::UInt32,
        frag_bins::Vector{Pioneer.FragIndexBin{T}},
        fragments::AbstractVector,
        frag_mz_min::Float32,
        frag_mz_max::Float32,
        hints_5da::Vector{UInt16},
        prev_first_match::UInt32,
        lb_factor::Float32,
        ub_overshoot_factor::Float32) where {T<:AbstractFloat}

    # ── Step 1: Hint-based lb advancement ────────────────────────────
    new_lb = lower_bound_guess
    hint_val = UInt16(1)
    est_step = 0.0f0
    delta_mz = 0.0f0

    if prev_first_match > zero(UInt32) && lower_bound_guess <= frag_bin_max_idx
        delta_mz = frag_mz_min - Pioneer.getLow(frag_bins[lower_bound_guess])
        if delta_mz > 0.0f0
            hint_val = @inbounds hints_5da[lower_bound_guess]
            est_step = Float32(hint_val) * (delta_mz / 5.0f0)

            if delta_mz > 5.0f0
                # Provably safe: hint says k bins = 5 Da in getLow.
                # At lb+k, getLow is only 5 Da higher — still below +delta_mz target.
                advance = UInt32(hint_val)
                new_lb = min(lower_bound_guess + advance, frag_bin_max_idx)
            else
                # Approximate: assume roughly linear density within 5 Da
                raw_advance = floor(Int, Float32(hint_val) * (delta_mz / 5.0f0) * lb_factor)
                advance = UInt32(max(raw_advance, 1))
                new_lb = min(lower_bound_guess + advance, frag_bin_max_idx)
            end
        end
    end

    # ── Step 2: Hint-based UB guess ──────────────────────────────────
    ub_guess = new_lb + UInt32(max(1, ceil(Int, est_step * ub_overshoot_factor)))
    ub_guess = min(ub_guess, frag_bin_max_idx)

    # ── Step 3: Check if UB guess is sufficient ──────────────────────
    ub_guess_sufficient = @inbounds Pioneer.getHigh(frag_bins[ub_guess]) >= frag_mz_max
    ub_final = ub_guess

    if !ub_guess_sufficient
        # Exponential doubling from ub_guess until sufficient
        step = one(UInt32)
        n = zero(UInt8)
        while @inbounds Pioneer.getHigh(frag_bins[ub_final]) < frag_mz_max
            ub_final += step
            step = step << one(UInt8)
            if ub_final > frag_bin_max_idx
                ub_final = frag_bin_max_idx
                break
            end
            n += one(UInt8)
        end
    end

    # ── Step 4: findFirstFragmentBin (provably safe) ─────────────────
    first_match = Pioneer.findFirstFragmentBin(
        frag_bins, new_lb, ub_final, frag_mz_min)

    # ── Step 5: Check lb safety ──────────────────────────────────────
    lb_overshoot = new_lb > first_match && first_match > zero(UInt32)

    # ── Step 6: Record diagnostics ───────────────────────────────────
    per_peak_jump = if prev_first_match > zero(UInt32) && first_match > zero(UInt32)
        first_match > prev_first_match ? first_match - prev_first_match : zero(UInt32)
    else
        zero(UInt32)
    end

    range_size = ub_final >= new_lb ? ub_final - new_lb + one(UInt32) : one(UInt32)

    bin_width = if first_match > zero(UInt32) && first_match <= frag_bin_max_idx
        @inbounds Pioneer.getHigh(frag_bins[first_match]) - Pioneer.getLow(frag_bins[first_match])
    else
        0.0f0
    end

    if prev_first_match > zero(UInt32) && delta_mz > 0.0f0
        push!(records, JumpRecord2(
            delta_mz, hint_val, est_step, per_peak_jump,
            new_lb, first_match, ub_guess, ub_final, range_size,
            lb_overshoot, ub_guess_sufficient, bin_width, lb_factor))
    end

    # ── Step 7: Score matching bins ──────────────────────────────────
    frag_bin_idx = first_match
    @inbounds @fastmath begin
        while frag_bin_idx <= frag_bin_max_idx
            frag_bin = frag_bins[frag_bin_idx]
            if Pioneer.getLow(frag_bin) > frag_mz_max
                break
            else
                if frag_bin_max_idx === frag_bin_idx
                    if Pioneer.getHigh(frag_bin) < frag_mz_min
                        break
                    end
                end
                frag_id_range = Pioneer.getSubBinRange(frag_bin)
                searchFragmentBinUnconditional!(counter, fragments, frag_id_range)
                frag_bin_idx += one(UInt32)
            end
        end
    end

    # Return (first_match, ub_final) — lb advances for next peak!
    return first_match, ub_final
end

"""
Instrumented partition scoring with the new v2 approach.
"""
function _score_partition_instrumented_v2!(
        records::Vector{JumpRecord2},
        local_counter::LocalCounter{UInt16, UInt8},
        partition::LocalPartition{T},
        irt_low::Float32,
        irt_high::Float32,
        masses::AbstractArray{Union{Missing, U}},
        mass_err_model::Pioneer.MassErrorModel,
        hints_5da::Vector{UInt16};
        lb_factor::Float32 = 0.8f0,
        ub_overshoot_factor::Float32 = 1.5f0,
        ) where {T<:AbstractFloat, U<:AbstractFloat}

    p_rt_bins   = getRTBins(partition)
    p_frag_bins = getFragBins(partition)
    p_fragments = getFragments(partition)

    isempty(p_frag_bins) && return nothing
    n_rt = length(p_rt_bins)

    local_rt_bin_idx = _find_rt_bin_start(p_rt_bins, irt_low)
    local_rt_bin_idx > n_rt && return nothing

    @inbounds @fastmath while Pioneer.getLow(p_rt_bins[local_rt_bin_idx]) < irt_high
        sub_bin_range = Pioneer.getSubBinRange(p_rt_bins[local_rt_bin_idx])
        min_frag_bin = first(sub_bin_range)
        max_frag_bin = last(sub_bin_range)

        if min_frag_bin <= max_frag_bin
            lower_bound_guess = min_frag_bin
            upper_bound_guess = min_frag_bin
            prev_first_match = zero(UInt32)

            for mass in masses
                corrected_mz = Pioneer.getCorrectedMz(mass_err_model, mass)
                frag_min, frag_max = Pioneer.getMzBoundsReverse(mass_err_model, corrected_mz)

                lower_bound_guess, upper_bound_guess = queryFragmentInstrumented_v2!(
                    records, local_counter, max_frag_bin,
                    lower_bound_guess, upper_bound_guess,
                    p_frag_bins, p_fragments,
                    frag_min, frag_max,
                    hints_5da, prev_first_match,
                    lb_factor, ub_overshoot_factor)
                prev_first_match = lower_bound_guess  # first_match was returned as lb
            end
        end

        local_rt_bin_idx += 1
        if local_rt_bin_idx > n_rt
            break
        end
    end

    return nothing
end

# ── Printing helpers ─────────────────────────────────────────────────────

function percentiles(v::Vector{<:Real}, ps=[0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.0])
    isempty(v) && return fill(NaN, length(ps))
    sv = sort(v)
    return [sv[max(1, ceil(Int, p * length(sv)))] for p in ps]
end

function print_percentile_table(name::String, values::Vector{<:Real};
                                 ps=[0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.0])
    pvals = percentiles(values, ps)
    labels = ["p$(Int(p*100))" for p in ps]
    labels[end] = "max"
    println("  $name (n=$(length(values)), mean=$(round(mean(values), digits=2))):")
    print("    ")
    for (l, v) in zip(labels, pvals)
        print("$l=$(round(v, digits=2))  ")
    end
    println()
end

function print_histogram(name::String, values::Vector{<:Real}, edges::Vector{<:Real})
    isempty(values) && return
    counts = zeros(Int, length(edges))
    for v in values
        for (i, e) in enumerate(edges)
            if v <= e
                counts[i] += 1
                break
            end
        end
    end
    println("  $name histogram (n=$(length(values))):")
    prev = 0.0
    for (i, e) in enumerate(edges)
        pct = round(100.0 * counts[i] / length(values), digits=1)
        label = i == length(edges) ? ">$(edges[end-1])" : "$(prev)-$(e)"
        bar = repeat("█", max(1, round(Int, pct / 2)))
        println("    $(rpad(label, 12)) $(lpad(counts[i], 8))  $(lpad(pct, 5))%  $bar")
        prev = e
    end
end

# ── Main ─────────────────────────────────────────────────────────────────

if isempty(ARGS)
    error("Usage: julia --project=. diagnose_hint_jumps.jl <params.json> [n_scans]")
end

const PARAMS_PATH = ARGS[1]
isfile(PARAMS_PATH) || error("Params file not found: $PARAMS_PATH")
const N_SCANS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 100

# ── Load & calibrate ────────────────────────────────────────────────────
println("Loading parameters from $PARAMS_PATH ...")
Pioneer.checkParams(PARAMS_PATH)
params = Pioneer.parse_pioneer_parameters(PARAMS_PATH)
mkpath(params.paths[:results])

MS_DATA_DIR = params.paths[:ms_data]
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, f) for f in readdir(MS_DATA_DIR)
                  if isfile(joinpath(MS_DATA_DIR, f)) && endswith(f, ".arrow")]
SPEC_LIB_DIR = params.paths[:library]
SPEC_LIB = Pioneer.loadSpectralLibrary(SPEC_LIB_DIR, params)

SC = Pioneer.initSearchContext(SPEC_LIB, Pioneer.parseIsoXML(Pioneer.isotope_spline_path()),
    Pioneer.ArrowTableReference(MS_TABLE_PATHS), Threads.nthreads(), 250000)
Pioneer.setDataOutDir!(SC, params.paths[:results])

Pioneer.execute_search(Pioneer.ParameterTuningSearch(), SC, params)
Pioneer.execute_search(Pioneer.NceTuningSearch(), SC, params)
Pioneer.execute_search(Pioneer.QuadTuningSearch(), SC, params)

search_parameters = Pioneer.FirstPassSearchParameters(params)
msdr = Pioneer.getMassSpecData(SC)
ms_file_idx, spectra = first(enumerate(msdr))
qtm = Pioneer.getQuadTransmissionModel(SC, ms_file_idx)
mem = Pioneer.getMassErrorModel(SC, ms_file_idx)
rt_to_irt_spline = Pioneer.getRtIrtModel(SC, ms_file_idx)
irt_tol = Pioneer.getIrtErrors(SC)[ms_file_idx]
precursor_mzs = Pioneer.getMz(Pioneer.getPrecursors(Pioneer.getSpecLib(SC)))

# ── Build partitioned index ─────────────────────────────────────────────
println("\nBuilding partitioned index (7 frags, weighted) ...")
pi = build_partitioned_index_from_lib(Pioneer.getSpecLib(SC);
    partition_width=5.0f0, frag_bin_tol_ppm=2.5f0, rt_bin_tol=3.0f0,
    rank_to_score=UInt8[8, 4, 4, 2, 2, 1, 1])

# ── Compute 5-Da hints for each partition ───────────────────────────────
println("\nComputing 5-Da direct hints ...")
partition_hints_5da = Vector{Vector{UInt16}}(undef, getNPartitions(pi))
for k in 1:getNPartitions(pi)
    p = getPartition(pi, k)
    partition_hints_5da[k] = _compute_skip_hints_5da(getFragBins(p), getRTBins(p))
end
total_hint_mem = sum(sizeof, partition_hints_5da)
println("  5-Da hints memory: $(round(total_hint_mem/1024^2, digits=2)) MB (UInt16)")

# Print hint statistics
all_hints = UInt16[]
for h in partition_hints_5da
    append!(all_hints, h)
end
println("  Hint value stats: min=$(minimum(all_hints)), median=$(Int(median(Float64.(all_hints)))), " *
        "mean=$(round(mean(Float64.(all_hints)), digits=1)), p95=$(Int(percentiles(Float64.(all_hints), [0.95])[1])), " *
        "max=$(maximum(all_hints))")

# ── Collect valid MS2 scans ─────────────────────────────────────────────
thread_tasks = Pioneer.partition_scans(spectra, Threads.nthreads())
all_scans = Int[]
for (_, ts) in thread_tasks
    for s in ts
        (s <= 0 || s > length(spectra)) && continue
        Pioneer.getMsOrder(spectra, s) ∉ Pioneer.getSpecOrder(search_parameters) && continue
        push!(all_scans, s)
    end
end
sort!(all_scans)

n_scans = min(N_SCANS, length(all_scans))
scan_subset = all_scans[1:n_scans]
println("\nDiagnosing on $n_scans / $(length(all_scans)) MS2 scans (single-threaded)")

# ── Run instrumented search with multiple lb_factor values ─────────────
iso_bounds = Pioneer.getIsotopeErrBounds(search_parameters)
max_local = maximum(p -> Int(p.n_local_precs), getPartitions(pi); init=0)
lc = LocalCounter(UInt16, UInt8, max_local + 1)

# Test multiple lb advancement factors
lb_factors = Float32[0.5, 0.7, 0.8, 0.9, 1.0]
ub_overshoot_factor = 1.5f0

all_records = Dict{Float32, Vector{JumpRecord2}}()
for lbf in lb_factors
    all_records[lbf] = JumpRecord2[]
end

println("\nRunning instrumented searches with lb_factors = $lb_factors ...")
for (fi, lbf) in enumerate(lb_factors)
    records = all_records[lbf]
    for scan_idx in scan_subset
        irt_lo, irt_hi = Pioneer.getRTWindow(
            rt_to_irt_spline(Pioneer.getRetentionTime(spectra, scan_idx)), irt_tol)
        quad_func = Pioneer.getQuadTransmissionFunction(
            qtm, Pioneer.getCenterMz(spectra, scan_idx),
            Pioneer.getIsolationWidthMz(spectra, scan_idx))
        prec_min = Float32(Pioneer.getPrecMinBound(quad_func) - Pioneer.NEUTRON * first(iso_bounds) / 2)
        prec_max = Float32(Pioneer.getPrecMaxBound(quad_func) + Pioneer.NEUTRON * last(iso_bounds) / 2)
        masses = Pioneer.getMzArray(spectra, scan_idx)

        first_k, last_k = get_partition_range(pi, prec_min, prec_max)
        for k in first_k:last_k
            partition = getPartition(pi, k)
            _score_partition_instrumented_v2!(records, lc, partition,
                irt_lo, irt_hi, masses, mem,
                partition_hints_5da[k];
                lb_factor=lbf, ub_overshoot_factor=ub_overshoot_factor)
            reset!(lc)
        end
    end
    println("  lb_factor=$lbf: $(length(records)) records")
end

# ── Analyze the default factor (0.8) in detail ──────────────────────────
println("\n" * "="^70)
println("5-Da HINT DIAGNOSTIC RESULTS (v2)")
println("="^70)

default_factor = 0.8f0
records = all_records[default_factor]
println("\nUsing lb_factor=$default_factor, ub_overshoot_factor=$ub_overshoot_factor")
println("Total records: $(length(records))")

if isempty(records)
    println("No records collected — nothing to analyze.")
    exit(0)
end

# Extract vectors
delta_mzs = Float32[r.delta_mz for r in records]
hints = UInt16[r.hint_5da for r in records]
est_steps = Float32[r.est_step for r in records]
per_peak_jumps = UInt32[r.per_peak_jump for r in records]
range_sizes = UInt32[r.range_size for r in records]
lb_overshoots = Bool[r.lb_overshoot for r in records]
ub_sufficient = Bool[r.ub_guess_sufficient for r in records]
bin_widths = Float32[r.bin_width for r in records]
first_matches = UInt32[r.first_match for r in records]

# ── 1. Per-peak jump (the TRUE metric with advancing lb) ────────────────
println("\n── 1. Per-Peak Jump (first_match_curr - first_match_prev) ──")
valid_jumps = per_peak_jumps[per_peak_jumps .> 0]
if !isempty(valid_jumps)
    print_percentile_table("per_peak_jump", Float64.(valid_jumps))
    print_histogram("per_peak_jump",
        Float64.(valid_jumps),
        Float64[1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, Inf])
else
    println("  No valid per-peak jumps recorded")
end

# ── 2. Delta m/z distribution ────────────────────────────────────────────
println("\n── 2. Delta m/z Between Consecutive Peaks ──")
print_percentile_table("delta_mz", Float64.(delta_mzs))

# ── 3. Hint accuracy: est_step vs per_peak_jump ─────────────────────────
println("\n── 3. Hint Accuracy (est_step vs per_peak_jump) ──")
valid_mask = per_peak_jumps .> 0
if sum(valid_mask) > 0
    ratios = Float64.(est_steps[valid_mask]) ./ Float64.(per_peak_jumps[valid_mask])
    print_percentile_table("est/actual ratio", ratios)
    print_histogram("est/actual ratio",
        ratios,
        Float64[0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0, Inf])

    n_overshoot = sum(est_steps[valid_mask] .> Float64.(per_peak_jumps[valid_mask]))
    n_undershoot = sum(est_steps[valid_mask] .< Float64.(per_peak_jumps[valid_mask]))
    n_valid = sum(valid_mask)
    println("\n  Overshoot (est > actual): $n_overshoot / $n_valid ($(round(100*n_overshoot/n_valid, digits=1))%)")
    println("  Undershoot (est < actual): $n_undershoot / $n_valid ($(round(100*n_undershoot/n_valid, digits=1))%)")
end

# ── 4. Search range sizes ────────────────────────────────────────────────
println("\n── 4. Search Range Sizes (ub_final - new_lb + 1) ──")
print_percentile_table("range_size", Float64.(range_sizes))
print_histogram("range_size",
    Float64.(range_sizes),
    Float64[1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, Inf])

# ── 5. UB guess tightness ───────────────────────────────────────────────
println("\n── 5. UB Guess Tightness ──")
n_sufficient = sum(ub_sufficient)
n_total = length(ub_sufficient)
println("  UB guess sufficient (no exponential needed): $n_sufficient / $n_total ($(round(100*n_sufficient/n_total, digits=1))%)")

# When ub_guess was sufficient, how tight was it?
has_match = first_matches .> 0
sufficient_and_match = ub_sufficient .& has_match
if sum(sufficient_and_match) > 0
    ub_guesses = UInt32[r.ub_guess for r in records]
    ub_slack = Int64.(ub_guesses[sufficient_and_match]) .- Int64.(first_matches[sufficient_and_match])
    print_percentile_table("ub_guess - first_match (when sufficient)", Float64.(ub_slack))
end

# ── 6. LB safety across all factors ─────────────────────────────────────
println("\n── 6. LB Safety Across Factors ──")
println("  lb_factor  |  n_records  |  overshoots  |  overshoot_rate  |  median_range")
for lbf in lb_factors
    recs = all_records[lbf]
    n = length(recs)
    n_over = sum(r.lb_overshoot for r in recs)
    rate = n > 0 ? round(100 * n_over / n, digits=2) : 0.0
    med_range = n > 0 ? round(median(Float64[r.range_size for r in recs]), digits=1) : NaN
    println("    $(rpad(lbf, 11)) | $(lpad(n, 9)) | $(lpad(n_over, 10)) | $(lpad(rate, 14))%  | $(med_range)")
end

# ── 7. Bin width distribution ────────────────────────────────────────────
println("\n── 7. Bin Width Distribution (getHigh - getLow) ──")
valid_bw = bin_widths[bin_widths .> 0.0f0]
if !isempty(valid_bw)
    print_percentile_table("bin_width", Float64.(valid_bw))
    print_histogram("bin_width (Da)",
        Float64.(valid_bw),
        Float64[0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, Inf])
end

# ── 8. 5-Da hint value distribution ─────────────────────────────────────
println("\n── 8. 5-Da Hint Values (bins to +5 Da) ──")
print_percentile_table("hint_5da", Float64.(hints))
print_histogram("hint_5da",
    Float64.(hints),
    Float64[10, 25, 50, 100, 200, 500, 1000, 2000, 5000, Inf])

# ── 9. Detailed factor sweep: per-peak jump coverage ────────────────────
println("\n── 9. LB Advancement Factor Sweep (detailed) ──")
println("  For each factor, what fraction of lb advancements are safe AND useful?")
for lbf in lb_factors
    recs = all_records[lbf]
    n = length(recs)
    n == 0 && continue

    n_over = sum(r.lb_overshoot for r in recs)
    n_safe = n - n_over

    # Among safe cases, how much did lb advance vs per_peak_jump?
    safe_recs = [r for r in recs if !r.lb_overshoot && r.per_peak_jump > 0]
    if !isempty(safe_recs)
        # How much of the per-peak jump did the lb advancement capture?
        advance_fracs = Float64[]
        for r in safe_recs
            lb_advance = Float64(r.new_lb) - Float64(r.new_lb) + Float64(r.per_peak_jump) > 0 ?
                (Float64(r.first_match) > Float64(r.new_lb) ?
                 Float64(r.first_match - r.new_lb) / Float64(r.per_peak_jump) : 0.0) : 0.0
            # Actually: remaining_gap = first_match - new_lb. Smaller is better (lb got closer).
            push!(advance_fracs, Float64(r.first_match) - Float64(r.new_lb))
        end
        med_gap = round(median(advance_fracs), digits=1)
        p90_gap = round(percentiles(advance_fracs, [0.90])[1], digits=1)
        println("  factor=$lbf: overshoots=$n_over/$n ($(round(100*n_over/n, digits=1))%), " *
                "median remaining_gap=$med_gap, p90_gap=$p90_gap")
    else
        println("  factor=$lbf: overshoots=$n_over/$n ($(round(100*n_over/n, digits=1))%), no safe records with jump>0")
    end
end

# ── 10. Comparison: old approach vs new approach range sizes ─────────────
println("\n── 10. Range Size Comparison ──")
# The old approach used saved_lb (never advancing), causing cumulative growth.
# Simulate: for each RT bin crossing in the records, sum range sizes.
println("  New approach (lb=first_match) range sizes already shown in section 4.")
println("  Key question: are ranges now small enough to skip binary search?")
small_range = sum(range_sizes .<= 16)
med_range = sum(16 .< range_sizes .<= 64)
large_range = sum(range_sizes .> 64)
n_total = length(range_sizes)
println("  Ranges <= 16 (linear scan): $small_range / $n_total ($(round(100*small_range/n_total, digits=1))%)")
println("  Ranges 17-64 (small binary): $med_range / $n_total ($(round(100*med_range/n_total, digits=1))%)")
println("  Ranges > 64 (large): $large_range / $n_total ($(round(100*large_range/n_total, digits=1))%)")

println("\n" * "="^70)
println("DONE")
println("="^70)
