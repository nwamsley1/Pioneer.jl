# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# DDA MS1 quantification: per-file window / smoothing calibration (ms1_quant.jl) and the fixed WH smoothing
# scale (`rt_scale`) of integrate_chrom.

# Gaussian MS1 trace on a uniform RT grid (minutes); FWHM in minutes.
_gauss(rt, mu, fwhm, h) = h .* exp.(-0.5 .* ((rt .- mu) ./ (fwhm / 2.3548f0)) .^ 2)

@testset "DDA MS1 calibration: local peak measurement" begin
    rt = Float32.(collect(0:0.01:2))              # 0.6 s spacing
    y = Float32.(_gauss(rt, 1.02f0, 0.05f0, 1f4))  # FWHM 3 s, apex 1.2 s after the MS2 scan at 1.00
    o, w, dt, cut = Pioneer._ms1_local_peak(rt, y, 1.0f0)
    @test o ≈ 0.02f0 atol = 0.005f0
    @test w ≈ 0.05f0 rtol = 0.10                  # 3-point pre-smoothing broadens slightly
    @test dt ≈ 0.01f0 atol = 1f-4
    @test !cut

    # climbs to the local peak at the MS2 scan, not the larger neighbour 18 s away
    y2 = y .+ Float32.(_gauss(rt, 1.32f0, 0.05f0, 5f4))
    o2, w2, _, _ = Pioneer._ms1_local_peak(rt, y2, 1.0f0)
    @test o2 ≈ 0.02f0 atol = 0.005f0
    @test w2 ≈ 0.05f0 rtol = 0.10

    # trace ending before the peak's far half-height crossing: flagged as cut off
    k = findall(r -> 0.9f0 <= r <= 1.03f0, rt)
    _, _, _, cut3 = Pioneer._ms1_local_peak(rt[k], y[k], 1.0f0)
    @test cut3

    # no signal / too short
    @test isnan(Pioneer._ms1_local_peak(rt, zeros(Float32, length(rt)), 1.0f0)[1])
    @test isnan(Pioneer._ms1_local_peak(rt[1:4], y[1:4], 0.01f0)[1])
end

@testset "DDA MS1 calibration: window, retry and smoothing rules" begin
    n = Pioneer.MS1_WINDOW_CALIB_MIN_PEAKS
    @test Pioneer.ms1_window_need(fill(0.01f0, n - 1), fill(0.05f0, n - 1)) == Inf32
    @test Pioneer.ms1_window_need(fill(0.01f0, n), fill(0.05f0, n)) ≈
          0.01f0 + Pioneer.MS1_WINDOW_FWHM_MULT * 0.05f0

    # too narrow when the peaks call for more than half the calibration window, or too many are cut off
    @test Pioneer.ms1_calib_too_narrow(0.11, 0.20, 0.0)
    @test !Pioneer.ms1_calib_too_narrow(0.11, 0.30, 0.0)
    @test Pioneer.ms1_calib_too_narrow(0.11, 0.30, 2 * Pioneer.MS1_WINDOW_CALIB_MAX_CUT_FRAC)
    @test Pioneer.ms1_calib_too_narrow(Inf32, 5.0, 0.0)

    # accept / widen / fallback
    act = Pioneer.ms1_calib_action
    @test act(n, 3000, 0.10, 0.60, 0.0, true) === :accept              # enough peaks, ample window
    @test act(n, 3000, 0.40, 0.60, 0.0, true) === :widen               # peaks call for > h0/2
    @test act(n, 3000, 0.40, 0.60, 0.0, false) === :accept             # cannot widen further: use what we have
    @test act(n - 1, 3000, Inf32, 0.05, 1.0, true) === :widen          # plenty of candidates, few peaks: too narrow
    @test act(n - 1, 3000, Inf32, 5.0, 0.0, false) === :fallback       # still too few after widening
    @test act(100, n - 1, Inf32, 0.60, 0.0, true) === :fallback        # few confident IDs: never widen

    # S = max(k * median FWHM, floor * median spacing)
    k = Pioneer.MS1_SMOOTH_FWHM_MULT; c = Pioneer.MS1_SMOOTH_MIN_CYCLES
    @test Pioneer.ms1_smoothing_scale(fill(0.05f0, 100), fill(0.01f0, 100)) ≈ max(k * 0.05f0, c * 0.01f0)
    @test Pioneer.ms1_smoothing_scale(fill(0.05f0, 100), fill(0.03f0, 100)) ≈ max(k * 0.05f0, c * 0.03f0)
    # medians: a few wild FWHM values do not move S
    w = vcat(fill(0.05f0, 95), fill(2.0f0, 5))
    @test Pioneer.ms1_smoothing_scale(w, fill(0.01f0, 100)) ≈ k * 0.05f0

    # calibrated value used; fallback when calibration failed
    @test Pioneer.ms1_smoothing_rt_scale(0.3f0) == 0.3f0
    @test Pioneer.ms1_smoothing_rt_scale(0f0) == Pioneer.MS1_SMOOTH_SCALE_DEFAULT
end

# synthetic scan tables for the acquisition-mode rule: per cycle one MS1 scan, then the cycle's MS2 isolation
# windows as (centre, width), each repeated `slices` times (mobility slices of one isolation event)
function _acq_scans(cycles; slices = 1)
    o = UInt8[]; c = Union{Missing, Float32}[]; w = Union{Missing, Float32}[]
    for cyc in cycles
        push!(o, 0x01); push!(c, missing); push!(w, missing)
        for (x, wd) in cyc, _ in 1:slices
            push!(o, 0x02); push!(c, Float32(x)); push!(w, Float32(wd))
        end
    end
    o, c, w
end
_acq_mode(cycles; kw...) = first(Pioneer.infer_acquisition_mode(_acq_scans(cycles; kw...)...))

@testset "acquisition mode: coverage per MS1 cycle" begin
    rng = Random.MersenneTwister(7)
    # Orbitrap DIA: 300 x 2 m/z windows tiling 380-980 every cycle
    @test _acq_mode([[(381 + 2.0k, 2.0) for k in 0:299] for _ in 1:60]) == "dia"
    # diaPASEF: 24 x 25 m/z windows, mobility slices
    @test _acq_mode([[(412.5 + 25.0k, 25.0) for k in 0:23] for _ in 1:60]; slices = 30) == "dia"
    # alternating windows: odd cycles offset by half a window
    @test _acq_mode([[(381 + 2.0k + (iseven(n) ? 0 : 1), 2.0) for k in 0:299] for n in 1:60]) == "dia"
    # diagonal / sliding windows (Synchro-PASEF-like): a 25 m/z window moving 1 m/z per slice over the cycle;
    # counted per isolation event this looked like DDA
    @test _acq_mode([[(400 + 1.0k, 25.0) for k in 0:599] for _ in 1:40]) == "dia"
    # DDA: top-20 of intensity-weighted random precursors (clustered at mid m/z), 2 m/z windows
    dda = [[(clamp(650 + 120randn(rng), 380, 1500), 2.0) for _ in 1:20] for _ in 1:600]
    @test _acq_mode(dda) == "dda"
    # short DDA run (old rule: < 1,000 distinct centres -> wrongly DIA)
    @test _acq_mode(dda[1:30]) == "dda"
    # targeted: 40 fixed targets re-isolated every cycle -> neither tiles the range nor samples it
    targets = sort(380 .+ 900rand(rng, 40))
    @test _acq_mode([[(t, 1.0) for t in targets] for _ in 1:200]) == "targeted"
    # nothing to classify; near-empty runs (fewer than ACQ_MIN_EVENTS isolation events) are not called
    @test _acq_mode([Tuple{Float64, Float64}[] for _ in 1:100]) == "unknown"
    @test _acq_mode(dda[1:3]) == "unknown"
    o, c, w = _acq_scans([[(500.0, 2.0)] for _ in 1:5])
    @test first(Pioneer.infer_acquisition_mode(o, c, w)) == "unknown"

    # labels written at conversion win; "unknown" (or none) falls back to the data rule
    @test Pioneer.acquisition_mode_label(Dict("acquisition_mode" => "DDA")) == "dda"
    @test Pioneer.acquisition_mode_label(Dict("acquisition_mode" => "unknown")) === nothing
    @test Pioneer.acquisition_mode_label(Dict("acquisition_type" => "swath")) === nothing
    @test Pioneer.acquisition_mode_label(nothing) === nothing
    mktempdir() do d
        o, c, w = _acq_scans([[(381 + 2.0k, 2.0) for k in 0:299] for _ in 1:30])
        n = length(o)
        df = DataFrame(mz_array = [Union{Missing, Float32}[500f0] for _ in 1:n], intensity_array = [Union{Missing, Float32}[1f0] for _ in 1:n],
                       scanHeader = fill("", n), scanNumber = Int32.(1:n), packetType = zeros(Int32, n),
                       lowMz = fill(100f0, n), highMz = fill(1700f0, n), TIC = ones(Float32, n),
                       collisionEnergyField = Vector{Union{Missing, Float32}}(fill(30f0, n)), collisionEnergyEvField = zeros(Float32, n),
                       retentionTime = Float32.(range(0, 10, length = n)), msOrder = o, centerMz = c, isolationWidthMz = w)
        Arrow.write(joinpath(d, "labelled.arrow"), df; metadata = Dict("acquisition_mode" => "dda", "acquisition_mode_source" => "vendor"))
        Arrow.write(joinpath(d, "unlabelled.arrow"), df)
        Arrow.write(joinpath(d, "unknown.arrow"), df; metadata = Dict("acquisition_mode" => "unknown"))
        @test first(Pioneer.file_acquisition_mode(Pioneer.BasicMassSpecData(joinpath(d, "labelled.arrow")))) == "dda"
        @test first(Pioneer.file_acquisition_mode(Pioneer.BasicMassSpecData(joinpath(d, "unlabelled.arrow")))) == "dia"
        mode, why = Pioneer.file_acquisition_mode(Pioneer.BasicMassSpecData(joinpath(d, "unknown.arrow")))
        @test mode == "dia" && occursin("converter recorded", why)
    end
end

@testset "DDA MS1 solve defaults are the validated configuration" begin
    # Huber (quadratic) on raw columns with ridge lambda = 0.01
    @test Pioneer.MS1_LAMBDA_DEFAULT == 0.01f0
    @test Pioneer.ms1_solver() isa Pioneer.HuberSolver
end

@testset "MS1 isotope matching: nearest peak to the offset-corrected target" begin
    iso = Pioneer.parseIsoXML(joinpath(pkgdir(Pioneer), "assets", "IsotopeSplines_10kDa_10isotopes.xml"))
    prec_mzs = Float32[500.0]; charges = UInt8[2]; sulfurs = UInt8[0]
    groups = Pioneer.Ms1MzGroups()
    Pioneer.build_ms1_mz_groups!(groups, UInt32[1], 1, prec_mzs, charges)
    # M0 window: a strong peak 9 ppm below the theoretical m/z, a weak one at the +6 ppm offset; tolerance +/-10 ppm
    t = 500f0
    mz = Float32[t * (1 - 9f-6), t * (1 + 6f-6)]; it = Float32[1000, 100]
    function m0(mem)
        Hs = Pioneer.SparseArrayFused(UInt32(100)); obs = zeros(Float32, 1)
        Pioneer.build_ms1_design_matrix!(Hs, groups, 1, mz, it, length(mz), prec_mzs, charges, sulfurs, iso, mem, obs,
                                         zeros(Float32, 1), UInt8[0])
        return obs[1]
    end
    @test m0(Pioneer.MassErrorModel(6f0, (10f0, 10f0))) == 100f0   # nearest to the corrected target (+6 ppm)
    @test m0(Pioneer.MassErrorModel(0f0, (10f0, 10f0))) == 100f0   # nearest to the theoretical m/z (6 vs 9 ppm)
end

@testset "MS1 design matrix: unmatched isotopes get zero-observation rows (as MS2)" begin
    iso = Pioneer.parseIsoXML(joinpath(pkgdir(Pioneer), "assets", "IsotopeSplines_10kDa_10isotopes.xml"))
    d = Pioneer.C13_C12_MASS_DIFF_F32 / 2
    # precursor 1 (500, 2+): M0 and M2 observed, M1 missing; precursor 2 (700, 2+): nothing observed
    prec_mzs = Float32[500.0, 700.0]; charges = UInt8[2, 2]; sulfurs = UInt8[0, 0]
    groups = Pioneer.Ms1MzGroups()
    Pioneer.build_ms1_mz_groups!(groups, UInt32[1, 2], 2, prec_mzs, charges)
    mz = Float32[500, 500 + 2d, 650]; it = Float32[1000, 200, 50]
    mem = Pioneer.MassErrorModel(0f0, (10f0, 10f0))
    Hs = Pioneer.SparseArrayFused(UInt32(100)); obs = zeros(Float32, 2); niso = zeros(UInt8, 2)
    Pioneer.build_ms1_design_matrix!(Hs, groups, 2, mz, it, length(mz), prec_mzs, charges, sulfurs, iso, mem, obs,
                                     zeros(Float32, 2), niso)
    c1 = Int(Hs.colptr[1]):Int(Hs.colptr[2]) - 1
    @test Int.(Hs.rowval[c1]) == [1, 2, 4]          # M0 -> peak 1, M2 -> peak 2, M1 -> synthetic row n_peaks + 1
    @test Hs.x[c1] == Float32[1000, 200, 0]
    @test Hs.nzval[last(c1)] ≈ iso(0, 1, (500f0 - Float32(Pioneer.PROTON)) * 2)
    @test niso[1] == 2
    @test Hs.colptr[3] == Hs.colptr[2]              # precursor 2: no match, empty column (its isotopes are out of range anyway)
    @test Hs.m == 4
    # with two isotopes modelled, M1 is still a miss and M2 is not modelled
    Pioneer.build_ms1_design_matrix!(Hs, groups, 2, mz, it, length(mz), prec_mzs, charges, sulfurs, iso, mem, obs,
                                     zeros(Float32, 2), niso, 0f0, 2)
    c1 = Int(Hs.colptr[1]):Int(Hs.colptr[2]) - 1
    @test Int.(Hs.rowval[c1]) == [1, 4]
    @test Pioneer.MS1_N_ISOTOPES == 3
end

@testset "DDA MS1 quant mass model: intensity-aware fit at a coverage target, and matching" begin
    # synthetic MS1 residuals: m/z-dependent bias (ppm) and a spread that narrows with intensity
    rng = Random.MersenneTwister(11)
    function synth(n)
        out = Pioneer.MassErrSample[]
        for _ in 1:n
            mz = 400f0 + 800f0 * rand(rng, Float32); l2 = 8f0 + 10f0 * rand(rng, Float32); rt = 20f0 * rand(rng, Float32)
            bias = 3f0 + 0.008f0 * (mz - 700f0)
            sd = 7f0 - 0.4f0 * (l2 - 8f0)                    # 7 ppm (weak) -> 3 ppm (strong)
            ppm = bias + sd * randn(rng, Float32)
            push!(out, Pioneer.MassErrSample(mz, mz * (1f0 + ppm * 1f-6), 2f0^l2, rt, UInt32(1)))
        end
        return out
    end
    train = synth(9000); test = synth(20000)
    fallback = Pioneer.MassErrorModel(3f0, (25f0, 25f0))
    m = Pioneer.fit_intensity_mass_error_model(train, fallback; coverage = 0.85)
    @test m isa Pioneer.IntensityMassErrorModel
    cov = Pioneer.ms1_model_coverage(test, m, fallback)
    @test 0.80 <= cov.coverage <= 0.90
    weak = [s.intensity < 2f0^10 for s in test]; strong = [s.intensity > 2f0^16 for s in test]
    @test median(cov.hw_ppm[weak]) > 1.5 * median(cov.hw_ppm[strong])
    @test abs(median(cov.err_ppm)) < 0.5                    # bias removed
    # the default coverage (fragments) is unchanged: wider windows than the 85% fit
    m95 = Pioneer.fit_intensity_mass_error_model(train, fallback)
    @test Pioneer.ms1_model_coverage(test, m95, fallback).median_hw > cov.median_hw

    # matching: candidates corrected with their own intensity; nearest corrected peak inside its window wins
    t = 700f0; I = 2f0^16; rtv = 10f0
    c_of(ppm) = t * (1f0 + (3f0 + ppm) * 1f-6)               # observed m/z with the 3 ppm bias plus `ppm`
    mz = Float32[c_of(-4f0), c_of(0.5f0), c_of(40f0)]; it = Float32[I, I / 4, I]
    p, v = Pioneer._ms1_match_peak(m, t, mz, it, 3, rtv)
    @test p == 2 && v == I / 4                               # nearest after correction, though weaker
    p, _ = Pioneer._ms1_match_peak(m, t, mz[3:3], it[3:3], 1, rtv)
    @test p == 0                                             # +40 ppm is outside the window
    @test Pioneer.MS1_QUANT_MODEL_COVERAGE == Pioneer.EMPK_QUANTILE      # same coverage rule as fragments
end

# integrate one trace; returns (area, apex RT, smoothed trace)
function _integrate(rt, y, seed_rt; rt_scale = 0f0)
    n = length(rt)
    sc = UInt32.(1:n)
    ws = Pioneer.WHWorkspace(n); st = Pioneer.Chromatogram(zeros(Float32, n), zeros(Float32, n), 0)
    dbg = Ref{Any}(nothing)
    seed = argmin(abs.(rt .- seed_rt))
    area, best, _ = Pioneer.integrate_chrom(rt, sc, y, ones(Float32, n), seed, ws, st,
                                            Float32((rt[end] - rt[1]) / n), 1f-6;
                                            debug_plot_data = dbg, rt_scale = rt_scale)
    return area, rt[Int(best)], copy(dbg[].wh_smoothed)
end

@testset "integrate_chrom: fixed smoothing scale is window-independent" begin
    rt = Float32.(collect(0:0.01:1))                       # 0.6 s spacing, 60 s
    mu = 0.5f0
    y = Float32.(_gauss(rt, mu, 0.05f0, 1f4) .+ _gauss(rt, mu - 0.2f0, 0.05f0, 3f4))  # bigger peak 12 s earlier
    near = findall(r -> abs(r - mu) <= 10f0 / 60f0, rt)    # +/-10 s window
    S = 20f0 / 60f0
    a10, apex10, _ = _integrate(rt[near], y[near], mu; rt_scale = S)
    a30, apex30, z30 = _integrate(rt, y, mu; rt_scale = S)
    @test apex10 ≈ mu atol = 0.011f0
    @test apex30 ≈ mu atol = 0.011f0
    @test a30 ≈ a10 rtol = 0.05

    # the old window-relative scale smooths a +/-30 s trace far more (penalty ~ width^6): lower smoothed apex
    _, _, z30_old = _integrate(rt, y, mu)
    i = argmin(abs.(rt .- mu))
    @test z30_old[i] < 0.8f0 * z30[i]
end

@testset "WHSmooth!: rt_scale = 0 keeps the window-relative behaviour (DIA path)" begin
    rt = Float32.(collect(0:0.01:0.4))
    y = Float32.(_gauss(rt, 0.2f0, 0.05f0, 1f4)) .+ Float32.(100 .* sin.(37 .* rt) .^ 2)
    n = length(rt); fr = ones(Float32, n)
    ws1 = Pioneer.WHWorkspace(n); ws2 = Pioneer.WHWorkspace(n)
    Pioneer.WHSmooth!(ws1, y, fr, rt, n, 0, 0f0, 1f-6)
    Pioneer.WHSmooth!(ws2, y, fr, rt, n, 0, 0f0, 1f-6; rt_scale = rt[end] - rt[1])
    @test ws1.z[1:n] == ws2.z[1:n]
end
