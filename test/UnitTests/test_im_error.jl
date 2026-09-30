# Tests for the per-file ion-mobility calibration: the z2 line (fit_im_line), the median line for files without
# one (fill_missing_im_lines!), and MainSearch's im_error / im_obs features (add_im_error!).

using Test
using Arrow
using DataFrames
using Pioneer
using Pioneer: fit_im_line, fill_missing_im_lines!, add_im_error!, BasicMassSpecData

# Fake library precursor table with / without inv_ion_mobility (struct must be top-level).
struct _ImTestPrecs; data::NamedTuple; end
Pioneer.getInvIonMobility(p::_ImTestPrecs) = hasproperty(p.data, :inv_ion_mobility) ? p.data.inv_ion_mobility : nothing

@testset "fit_im_line — least squares on the calib rows, MAD sigma, nothing when short" begin
    scan = Float32[]; pred = Float32[]; calib = Bool[]
    for k in 1:300
        s = 200f0 + 600f0 * (k / 300); noise = 0.01f0 * (isodd(k) ? 1 : -1)
        push!(scan, s); push!(pred, 1.45f0 - 0.0009f0 * s + noise); push!(calib, true)
    end
    # rows outside calib must be ignored even with wild values
    for k in 1:50
        push!(scan, 100f0); push!(pred, 5f0); push!(calib, false)
    end
    a, b, s = fit_im_line(scan, pred, calib; min_calib = 100)
    @test isapprox(a, 1.45f0; atol = 0.005) && isapprox(b, -0.0009f0; atol = 1e-5)
    @test isapprox(s, 1.4826f0 * 0.01f0; rtol = 0.15)         # MAD of ±0.01 noise
    @test fit_im_line(scan[1:50], pred[1:50], calib[1:50]; min_calib = 100) === nothing
end

@testset "fill_missing_im_lines! — median line in 1/K0 units, mapped to each file's scans" begin
    # Library 1/K0 = 0.02 + 1.1 * k0 in every file (k0 = instrument 1/K0). Files 1-2 are one scan scale, file 3
    # another; the line in scan units is pred = (0.02 + 1.1 * c0) + (1.1 * m) * scan.
    true_line(c0, m) = (0.02f0 + 1.1f0 * c0, 1.1f0 * m)
    cal_a = (1.45f0, -0.000865f0); cal_b = (1.45f0, -0.00085f0)
    models = Dict{Int64, Dict{Int, NTuple{3, Float32}}}(
        1 => Dict(2 => (true_line(cal_a...)..., 0.013f0)),
        2 => Dict(2 => (true_line(cal_a...)..., 0.015f0)),
        # 3: no line
    )
    cals = Dict{Int64, NTuple{2, Float32}}(1 => cal_a, 2 => cal_a, 3 => cal_b)
    fill_missing_im_lines!(models, cals, 3)
    a3, b3, s3 = models[3][2]
    ea, eb = true_line(cal_b...)
    @test isapprox(a3, ea; atol = 1e-4) && isapprox(b3, eb; rtol = 1e-4)   # file 3's own scan scale
    @test isapprox(s3, 0.014f0; atol = 1e-6)                                # median sigma
    @test models[1][2] == (true_line(cal_a...)..., 0.013f0)                 # donors untouched

    # without an instrument calibration for every file: median in scan units
    m2 = Dict{Int64, Dict{Int, NTuple{3, Float32}}}(1 => Dict(2 => (1.5f0, -0.0009f0, 0.01f0)),
                                                    2 => Dict(2 => (1.6f0, -0.0011f0, 0.03f0)),
                                                    3 => Dict(2 => (1.7f0, -0.0010f0, 0.02f0)))
    fill_missing_im_lines!(m2, Dict{Int64, NTuple{2, Float32}}(1 => cal_a), 4)
    @test m2[4][2] == (1.6f0, -0.0010f0, 0.02f0)

    # nobody has a line: nothing is filled (no gate, im_error 0)
    m3 = Dict{Int64, Dict{Int, NTuple{3, Float32}}}(1 => Dict{Int, NTuple{3, Float32}}())
    fill_missing_im_lines!(m3, cals, 3)
    @test !haskey(m3[1], 2) && !haskey(m3, 2) && !haskey(m3, 3)
end

@testset "add_im_error! — signed z2-sigma residuals for every charge; fallback line; no data" begin
    mktempdir() do d
        # spectra: 400 packet rows with an imScan column
        n = 400
        df = DataFrame(mz_array = [Union{Missing,Float32}[500f0] for _ in 1:n], intensity_array = [Union{Missing,Float32}[1f0] for _ in 1:n],
                       scanHeader = fill("", n), scanNumber = Int32.(1:n), packetType = zeros(Int32, n),
                       retentionTime = Float32.(1:n), lowMz = fill(100f0, n), highMz = fill(1700f0, n), TIC = ones(Float32, n),
                       centerMz = Vector{Union{Missing,Float32}}(fill(612.5f0, n)), isolationWidthMz = Vector{Union{Missing,Float32}}(fill(25f0, n)),
                       collisionEnergyField = Vector{Union{Missing,Float32}}(fill(30f0, n)), collisionEnergyEvField = zeros(Float32, n),
                       msOrder = fill(0x02, n), cycle_idx = Int32.(1:n), imScan = UInt16.(200 .+ (1:n)))
        p = joinpath(d, "packets.arrow"); Arrow.write(p, df)
        spectra = BasicMassSpecData(p)

        # precursors 1:n on the line pred = 1.45 - 0.0009*scan with spread deterministic noise
        # (|noise| <= 0.005, so sigma is ~0.004 rather than the floor), plus two 0.1 outliers
        scans = Float32.(200 .+ (1:n))
        im_pred = 1.45f0 .- 0.0009f0 .* scans .+ Float32[0.005f0 * sin(1.7 * k) for k in 1:n]
        im_pred[10] += 0.1f0; im_pred[20] -= 0.1f0          # outliers: large |im_error|, opposite signs
        precs = _ImTestPrecs((inv_ion_mobility = im_pred,))
        charges = fill(0x02, n); charges[301:end] .= 0x03    # 3+ rows: scored on the z2 line, not fitted
        im_pred[301:end] .+= 0.2f0                           # ... and far off it, so a 3+ fit would be wrong
        mkbest() = DataFrame(precursor_idx = UInt32.(1:n), scan_idx = UInt32.(1:n), charge = copy(charges),
                             target = trues(n), lgbm_prob = fill(0.99f0, n))
        best = mkbest()
        best.target[10] = false                              # outlier 10 is a decoy: excluded from the fit, still scored
        m = add_im_error!(best, best.lgbm_prob, spectra, precs, 1; min_prob = 0.9f0, min_calib = 100)
        @test collect(keys(m)) == [2]
        a, b, s = m[2]
        @test isapprox(a, 1.45f0; atol = 0.005) && isapprox(b, -0.0009f0; atol = 2e-5)
        @test eltype(best.im_error) == Float32
        z2 = setdiff(1:300, [10, 20])
        @test maximum(abs, best.im_error[z2]) < 2f0
        @test best.im_error[10] > 8f0 && best.im_error[20] < -8f0          # signed
        @test all(>(20f0), best.im_error[301:end])                         # 3+ measured from the z2 line
        @test best.im_obs ≈ a .+ b .* scans

        # too few z2 PSMs for its own line: the fallback line is used as given
        fb = (1.40f0, -0.0008f0, 0.02f0)
        best_fb = mkbest()
        m_fb = add_im_error!(best_fb, best_fb.lgbm_prob, spectra, precs, 1; fallback = fb, min_calib = 1000)
        @test m_fb[2] == fb
        @test best_fb.im_obs ≈ fb[1] .+ fb[2] .* scans
        @test best_fb.im_error ≈ (im_pred .- best_fb.im_obs) ./ fb[3]

        # neither its own line nor a fallback: im_error 0, im_obs NaN (MBR: missing)
        best_none = mkbest()
        @test isempty(add_im_error!(best_none, best_none.lgbm_prob, spectra, precs, 1; min_calib = 1000))
        @test all(iszero, best_none.im_error) && all(isnan, best_none.im_obs)

        # library without mobility -> zeros (im_obs too: unchanged for data without ion mobility)
        best2 = mkbest()
        add_im_error!(best2, best2.lgbm_prob, spectra, _ImTestPrecs((mz = scans,)), 1)
        @test all(iszero, best2.im_error) && all(iszero, best2.im_obs)

        # spectra without imScan -> zeros
        q = joinpath(d, "noim.arrow"); Arrow.write(q, select(df, Not(:imScan)))
        best3 = mkbest()
        add_im_error!(best3, best3.lgbm_prob, BasicMassSpecData(q), precs, 1)
        @test all(iszero, best3.im_error) && all(iszero, best3.im_obs)
    end
end

@testset "MBR observed-IM difference treats a NaN im_obs as missing" begin
    donor(im) = Pioneer._MBRDonorEntry(0.99f0, UInt32(1), 100f0, -1f0, 0f0, 10f0, im, 5f0,
                                       ntuple(_ -> 0f0, 8), UInt8(0x03), UInt16(2), UInt32(2))
    @test Pioneer._mbr_observed_im_diff(1.0f0, donor(1.01f0)) ≈ 0.01f0
    @test Pioneer._mbr_observed_im_diff(NaN32, donor(1.01f0)) == -1f0
    @test Pioneer._mbr_observed_im_diff(1.0f0, donor(NaN32)) == -1f0
end
