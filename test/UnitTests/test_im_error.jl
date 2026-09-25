# Tests for MainSearch's per-file ion-mobility calibration feature (fit_im_lines / add_im_error!).

using Test
using Arrow
using DataFrames
using Pioneer
using Pioneer: fit_im_lines, add_im_error!, BasicMassSpecData

# Fake library precursor table with / without inv_ion_mobility (struct must be top-level).
struct _ImTestPrecs; data::NamedTuple; end
Pioneer.getInvIonMobility(p::_ImTestPrecs) = hasproperty(p.data, :inv_ion_mobility) ? p.data.inv_ion_mobility : nothing

@testset "fit_im_lines — per-charge lines with pooled fallback" begin
    # 2+ : pred = 1.45 - 0.0009*scan ; 3+ : pred = 1.30 - 0.0007*scan ; 4+ : too few rows
    scan = Float32[]; pred = Float32[]; charge = UInt8[]; calib = Bool[]
    for k in 1:300
        s = 200f0 + 600f0 * (k / 300); noise = 0.01f0 * (isodd(k) ? 1 : -1)
        push!(scan, s); push!(pred, 1.45f0 - 0.0009f0 * s + noise); push!(charge, 0x02); push!(calib, true)
    end
    for k in 1:150
        s = 300f0 + 500f0 * (k / 150); noise = 0.02f0 * (isodd(k) ? 1 : -1)
        push!(scan, s); push!(pred, 1.30f0 - 0.0007f0 * s + noise); push!(charge, 0x03); push!(calib, true)
    end
    for k in 1:20
        push!(scan, 500f0); push!(pred, 0.9f0); push!(charge, 0x04); push!(calib, true)
    end
    # decoys / low-score rows must be ignored even with wild values
    for k in 1:50
        push!(scan, 100f0); push!(pred, 5f0); push!(charge, 0x02); push!(calib, false)
    end
    models = fit_im_lines(scan, pred, charge, calib; min_calib = 100)
    @test haskey(models, 0) && haskey(models, 2) && haskey(models, 3) && !haskey(models, 4)
    a2, b2, s2 = models[2]
    @test isapprox(a2, 1.45f0; atol = 0.005) && isapprox(b2, -0.0009f0; atol = 1e-5)
    @test isapprox(s2, 1.4826f0 * 0.01f0; rtol = 0.15)        # MAD of ±0.01 noise
    a3, b3, s3 = models[3]
    @test isapprox(a3, 1.30f0; atol = 0.01) && isapprox(b3, -0.0007f0; atol = 2e-5)
    @test s3 > s2
    # too few calibration rows overall -> no model
    @test isempty(fit_im_lines(scan[1:50], pred[1:50], charge[1:50], calib[1:50]; min_calib = 100))
end

@testset "add_im_error! — sigma-scaled residuals; zeros without mobility data" begin
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
        im_pred[10] += 0.1f0; im_pred[20] -= 0.1f0          # outliers: large im_error expected
        precs = _ImTestPrecs((inv_ion_mobility = im_pred,))
        best = DataFrame(precursor_idx = UInt32.(1:n), scan_idx = UInt32.(1:n), charge = fill(0x02, n),
                         target = trues(n), lgbm_prob = fill(0.99f0, n))
        best.target[10] = false                              # outlier 10 is a decoy: excluded from the fit, still scored
        add_im_error!(best, best.lgbm_prob, spectra, precs, 1; min_prob = 0.9f0, min_calib = 100)
        @test hasproperty(best, :im_error) && eltype(best.im_error) == Float32
        # in-line points sit within ~1 sigma; the two 0.1 outliers are >10 sigma out
        @test maximum(best.im_error[Not([10, 20])]) < 2f0
        @test best.im_error[10] > 8f0 && best.im_error[20] > 8f0

        # library without mobility -> zeros
        best2 = DataFrame(precursor_idx = UInt32.(1:n), scan_idx = UInt32.(1:n), charge = fill(0x02, n),
                          target = trues(n), lgbm_prob = fill(0.99f0, n))
        add_im_error!(best2, best2.lgbm_prob, spectra, _ImTestPrecs((mz = scans,)), 1)
        @test all(iszero, best2.im_error)

        # spectra without imScan -> zeros
        q = joinpath(d, "noim.arrow"); Arrow.write(q, select(df, Not(:imScan)))
        best3 = DataFrame(precursor_idx = UInt32.(1:n), scan_idx = UInt32.(1:n), charge = fill(0x02, n),
                          target = trues(n), lgbm_prob = fill(0.99f0, n))
        add_im_error!(best3, best3.lgbm_prob, BasicMassSpecData(q), precs, 1)
        @test all(iszero, best3.im_error)
    end
end
