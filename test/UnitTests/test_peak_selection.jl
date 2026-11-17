# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Unit tests for fragment peak selection (PeakSelection: NearestPeak / MostIntensePeak).
#
# Target functions (src/Routines/SearchDIA/CommonSearchUtils/fusedMatch.jl):
#   scan_for_nearest_in_window, scan_for_most_intense_in_window, scan_for_match
#   run_fused!(...; peak_select), run_fused_masserr!(...; peak_select)
#
# Run: julia --project=. test/UnitTests/test_peak_selection.jl

using Test
using Pioneer
using Random

# Reuses make_fused_fixture from test_run_fused.jl (already loaded when run from runtests).
isdefined(@__MODULE__, :make_fused_fixture) || include(joinpath(@__DIR__, "test_run_fused.jl"))

const NearestPeak     = Pioneer.NearestPeak
const MostIntensePeak = Pioneer.MostIntensePeak
const scan_for_nearest_in_window      = Pioneer.scan_for_nearest_in_window
const scan_for_most_intense_in_window = Pioneer.scan_for_most_intense_in_window
const scan_for_match                  = Pioneer.scan_for_match
const run_fused_masserr!              = Pioneer.run_fused_masserr!
const MassErrSample                   = Pioneer.MassErrSample

# Plain scalar references: the definitions the SIMD scans must reproduce exactly.
function ref_scan(corrected_mz, obs_low, obs_high, start_idx, n_peaks, iso_mz, conservative_high, scan_int, most_intense)
    best_peak, best_key = 0, Inf32
    j = start_idx
    while j <= n_peaks && corrected_mz[j] <= conservative_high
        if obs_low[j] <= iso_mz <= obs_high[j]
            key = most_intense ? -(ismissing(scan_int[j]) ? 0f0 : Float32(scan_int[j])) : abs(corrected_mz[j] - iso_mz)
            if best_peak == 0 || key < best_key
                best_peak, best_key = j, key
            end
        end
        j += 1
    end
    best_peak == 0 && return (0, typemax(Float32), 0f0)
    return (best_peak, abs(corrected_mz[best_peak] - iso_mz), corrected_mz[best_peak])
end

# A scan with peaks spread over [lo, hi] m/z, per-peak ppm windows, some missing and some tied intensities.
function random_scan(rng, n; lo = 300f0, hi = 301f0, ppm_range = (5f0, 120f0), p_missing = 0.05, p_tie = 0.2)
    mz = sort!(lo .+ (hi - lo) .* rand(rng, Float32, n))
    tol = ppm_range[1] .+ (ppm_range[2] - ppm_range[1]) .* rand(rng, Float32, n)
    obs_low  = mz .- mz .* tol .* 1f-6
    obs_high = mz .+ mz .* tol .* 1f-6
    ints = Vector{Union{Missing, Float32}}(undef, n)
    for i in 1:n
        ints[i] = rand(rng) < p_missing ? missing :
                  rand(rng) < p_tie     ? 1000f0 : Float32(round(1000 * rand(rng)))
    end
    return mz, obs_low, obs_high, ints
end

@testset "peak selection" begin

    @testset "SIMD scans match the scalar references (random scans, all block offsets)" begin
        rng = MersenneTwister(20261114)
        n_checked = 0; n_matched = 0; n_differ = 0
        for trial in 1:400
            n = rand(rng, 0:41)                 # covers empty, < 8, exact blocks and tails
            mz, lo, hi, ints = random_scan(rng, n)
            for _ in 1:10
                iso = 300f0 + rand(rng, Float32)
                chigh = iso + iso * Float32(rand(rng, (20, 60, 150))) * 1f-6
                start = n == 0 ? 1 : rand(rng, 1:n)
                for (most, f) in ((false, scan_for_nearest_in_window), (true, scan_for_most_intense_in_window))
                    got = most ? f(mz, lo, hi, start, n, iso, chigh, ints) : f(mz, lo, hi, start, n, iso, chigh)
                    @test got === ref_scan(mz, lo, hi, start, n, iso, chigh, ints, most)
                end
                a = scan_for_nearest_in_window(mz, lo, hi, start, n, iso, chigh)
                b = scan_for_most_intense_in_window(mz, lo, hi, start, n, iso, chigh, ints)
                n_checked += 1; n_matched += a[1] > 0; n_differ += a[1] != b[1]
            end
        end
        # The random scans must actually exercise matches and disagreements, or the test proves nothing.
        @test n_matched > 500
        @test n_differ > 100
    end

    @testset "scan_for_match dispatches to the right scan" begin
        rng = MersenneTwister(7)
        mz, lo, hi, ints = random_scan(rng, 30)
        for _ in 1:200
            iso = 300f0 + rand(rng, Float32); chigh = iso + iso * 150f-6; start = rand(rng, 1:30)
            @test scan_for_match(NearestPeak(), mz, lo, hi, start, 30, iso, chigh, ints) ===
                  scan_for_nearest_in_window(mz, lo, hi, start, 30, iso, chigh)
            @test scan_for_match(MostIntensePeak(), mz, lo, hi, start, 30, iso, chigh, ints) ===
                  scan_for_most_intense_in_window(mz, lo, hi, start, 30, iso, chigh, ints)
        end
    end

    @testset "one candidate per window: both rules agree" begin
        # Peaks 1 m/z apart with 20 ppm windows: at most one window covers any target.
        mz = Float32.(300:315); n = length(mz)
        lo = mz .- mz .* 20f-6; hi = mz .+ mz .* 20f-6
        ints = Vector{Union{Missing, Float32}}(Float32.(100:100:1600))
        for k in 1:n, off_ppm in (-15f0, 0f0, 12f0)
            iso = mz[k] * (1 + off_ppm * 1f-6); chigh = iso + iso * 40f-6
            @test scan_for_nearest_in_window(mz, lo, hi, 1, n, iso, chigh) ===
                  scan_for_most_intense_in_window(mz, lo, hi, 1, n, iso, chigh, ints)
            @test scan_for_nearest_in_window(mz, lo, hi, 1, n, iso, chigh)[1] == k
        end
    end

    @testset "satellite next to a shifted strong peak (the scout failure)" begin
        # Theoretical 500; the real peak reads -38 ppm, a weak satellite +2 ppm. Wide ±100 ppm windows.
        iso = 500f0
        for pad in (0, 3, 9)   # put the pair in the scalar tail and inside SIMD blocks
            mz = vcat(Float32[490 + i for i in 1:pad], Float32[iso * (1 - 38f-6), iso * (1 + 2f-6)])
            n = length(mz)
            lo = mz .- mz .* 100f-6; hi = mz .+ mz .* 100f-6
            ints = Vector{Union{Missing, Float32}}(vcat(fill(5f0, pad), Float32[1.0f5, 2.0f3]))
            chigh = iso + iso * 100f-6
            near = scan_for_nearest_in_window(mz, lo, hi, 1, n, iso, chigh)
            most = scan_for_most_intense_in_window(mz, lo, hi, 1, n, iso, chigh, ints)
            @test near[1] == pad + 2                       # satellite
            @test most[1] == pad + 1                       # real peak
            @test most[3] == mz[pad + 1]
            @test most[2] ≈ abs(mz[pad + 1] - iso)
        end
    end

    @testset "edge cases" begin
        mz = Float32[100, 200, 300]; lo = mz .- 1f0; hi = mz .+ 1f0
        ints = Union{Missing, Float32}[1f0, 2f0, 3f0]
        # no match
        @test scan_for_most_intense_in_window(mz, lo, hi, 1, 3, 150f0, 400f0, ints) === (0, typemax(Float32), 0f0)
        @test scan_for_most_intense_in_window(Float32[], Float32[], Float32[], 1, 0, 150f0, 400f0,
                                              Union{Missing, Float32}[]) === (0, typemax(Float32), 0f0)
        # early stop: a more intense covering peak past conservative_high is not taken
        mz = Float32[200.0, 200.1]; lo = Float32[199, 199]; hi = Float32[201, 201]
        ints = Union{Missing, Float32}[1f0, 9f0]
        @test scan_for_most_intense_in_window(mz, lo, hi, 1, 2, 200f0, 200.05f0, ints)[1] == 1
        @test scan_for_most_intense_in_window(mz, lo, hi, 1, 2, 200f0, 200.2f0, ints)[1] == 2
        # start_idx skips earlier peaks
        @test scan_for_most_intense_in_window(Float32[200, 200.01], Float32[199, 199], Float32[201, 201], 2, 2,
                                              200f0, 201f0, Union{Missing, Float32}[9f0, 1f0])[1] == 2
        # ties: first (lowest m/z) peak; missing intensity counts as 0 but still matches
        mz = Float32[200.0, 200.01, 200.02]; lo = mz .- 1f0; hi = mz .+ 1f0
        @test scan_for_most_intense_in_window(mz, lo, hi, 1, 3, 200.02f0, 201f0,
                                              Union{Missing, Float32}[5f0, 5f0, 1f0])[1] == 1
        @test scan_for_most_intense_in_window(mz, lo, hi, 1, 3, 200.02f0, 201f0,
                                              Union{Missing, Float32}[missing, missing, missing])[1] == 1
        @test scan_for_most_intense_in_window(mz, lo, hi, 1, 3, 200.02f0, 201f0,
                                              Union{Missing, Float32}[missing, 0.5f0, missing])[1] == 2
        # plain Vector{Float32} intensities work too
        @test scan_for_most_intense_in_window(mz, lo, hi, 1, 3, 200.02f0, 201f0, Float32[1, 3, 2])[1] == 2
    end

    # run_fused! with an explicit peak_select (same argument order as call_run_fused!).
    function call_run_fused_select!(fx, peak_select)
        intensity_model = Pioneer.prepare_fragment_intensity_model(fx.ion_list, fx.nce_model)
        return Pioneer.run_fused!(
            fx.kind, fx.Hs, fx.unscored_psms, fx.id_to_col, fx.scratch,
            fx.scan_corrected_mz, fx.scan_obs_low, fx.scan_obs_high, fx.peak_mz_len,
            fx.isotopes_buf, fx.prec_trans_buf, fx.ion_list, intensity_model,
            fx.precursors_passed, fx.prec_range,
            fx.prec_mzs, fx.prec_charges, fx.prec_sulfur_counts, fx.prec_irts,
            fx.iso_splines, fx.qfunc, fx.mem, fx.scan_int,
            fx.scan_irt, fx.irt_tol, fx.frag_mz_bounds, fx.n_frag_isotopes,
            fx.isotope_err_bounds; peak_select = peak_select)
    end
    hs_state(fx) = (fx.Hs.n, fx.Hs.n_vals, copy(fx.Hs.rowval[1:fx.Hs.n_vals]), copy(fx.Hs.x[1:fx.Hs.n_vals]),
                    copy(fx.Hs.nzval[1:fx.Hs.n_vals]))

    # Fragments at 200 and 300; each has a strong peak -38 ppm and a weak one +2 ppm (±100 ppm windows).
    sat_peaks = Float32[200 * (1 - 38f-6), 200 * (1 + 2f-6), 300 * (1 - 38f-6), 300 * (1 + 2f-6)]
    sat_ints  = Float32[9000, 300, 8000, 200]

    @testset "run_fused!: default is NearestPeak, MostIntensePeak picks the strong peaks" begin
        fx_default = make_fused_fixture(peak_mz = sat_peaks, peak_int = sat_ints, ppm_tol = 100f0)
        fx_near    = make_fused_fixture(peak_mz = sat_peaks, peak_int = sat_ints, ppm_tol = 100f0)
        fx_most    = make_fused_fixture(peak_mz = sat_peaks, peak_int = sat_ints, ppm_tol = 100f0)
        r_default = call_run_fused!(fx_default)
        r_near    = call_run_fused_select!(fx_near, NearestPeak())
        r_most    = call_run_fused_select!(fx_most, MostIntensePeak())
        @test r_default == r_near == r_most == (2, 0)
        @test hs_state(fx_default) == hs_state(fx_near)
        @test fx_near.Hs.rowval[1:2] == UInt32[2, 4]      # satellites
        @test fx_most.Hs.rowval[1:2] == UInt32[1, 3]      # strong peaks
        @test fx_most.Hs.x[1:2] == Float32[9000, 8000]    # observed intensities of the strong peaks
    end

    @testset "run_fused!: one candidate per window gives identical output for both rules" begin
        fx_near = make_fused_fixture()
        fx_most = make_fused_fixture()
        @test call_run_fused_select!(fx_near, NearestPeak()) == call_run_fused_select!(fx_most, MostIntensePeak())
        @test hs_state(fx_near) == hs_state(fx_most)
    end

    @testset "run_fused_masserr!: default is NearestPeak, MostIntensePeak samples the strong peaks" begin
        # y ions, rank 1, ion position 5 (passes the mass-error filter).
        frags = [(UInt32(1), 200.0f0, Float32(1000), UInt8(1), :y, UInt8(5)),
                 (UInt32(1), 300.0f0, Float32( 800), UInt8(1), :y, UInt8(5))]
        fx = make_fused_fixture(peak_mz = sat_peaks, peak_int = sat_ints, ppm_tol = 100f0, frags = frags)
        scan_mz = Vector{Union{Missing, Float32}}(sat_peaks)
        run(sel) = begin
            samples = [MassErrSample() for _ in 1:4]
            n = sel === nothing ?
                run_fused_masserr!(samples, 0, fx.scan_corrected_mz, fx.scan_obs_low, fx.scan_obs_high,
                                   fx.peak_mz_len, fx.ion_list, fx.precursors_passed, 1:1, fx.mem,
                                   scan_mz, fx.scan_int, 10f0) :
                run_fused_masserr!(samples, 0, fx.scan_corrected_mz, fx.scan_obs_low, fx.scan_obs_high,
                                   fx.peak_mz_len, fx.ion_list, fx.precursors_passed, 1:1, fx.mem,
                                   scan_mz, fx.scan_int, 10f0; peak_select = sel)
            [(s.theoretical_mz, s.observed_mz, s.intensity) for s in samples[1:n]]
        end
        d = run(nothing); near = run(NearestPeak()); most = run(MostIntensePeak())
        @test d == near
        @test [x[2] for x in near] == sat_peaks[[2, 4]]
        @test [x[2] for x in most] == sat_peaks[[1, 3]]
        @test [x[3] for x in most] == Float32[9000, 8000]
    end
end
