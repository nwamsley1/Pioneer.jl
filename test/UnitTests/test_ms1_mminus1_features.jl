using Test
using DataFrames
import Pioneer: _ms1_lookup_chunk!, PRESCORE_FEATURES, ADVANCED_FEATURE_SET,
                feature_matrix, getMzArray, getIntensityArray

struct FakeMS1Spectra
    mzs::Vector{Vector{Union{Missing, Float32}}}
    intensities::Vector{Vector{Union{Missing, Float32}}}
end

getMzArray(spectra::FakeMS1Spectra, scan_idx::Integer) = spectra.mzs[scan_idx]
getIntensityArray(spectra::FakeMS1Spectra, scan_idx::Integer) = spectra.intensities[scan_idx]

@testset "MS1 M-1 lookup features" begin
    neutron = 1.00335f0
    target_m0 = 500f0
    target_mminus1 = target_m0 - neutron / 2f0
    observed_mminus1 = target_mminus1 * (1f0 + 2f-6)

    spectra = FakeMS1Spectra(
        [
            Union{Missing, Float32}[
                observed_mminus1,
                target_m0,
                target_m0 + neutron / 2f0,
                600f0,
                600f0 + neutron / 2f0,
            ],
        ],
        [
            Union{Missing, Float32}[80f0, 200f0, 60f0, 100f0, 30f0],
        ],
    )
    psms = DataFrame(
        scan_idx = Int32[2, 2],
        precursor_idx = UInt32[1, 2],
        ms1_m0_mass_err_ppm = zeros(Float32, 2),
        ms1_m0_intensity = zeros(Float32, 2),
        ms1_m1_intensity = zeros(Float32, 2),
        ms1_m1_to_m0_ratio = zeros(Float32, 2),
        ms1_m1_to_m0_pred = zeros(Float32, 2),
        ms1_mminus1_present = zeros(Float32, 2),
        ms1_mminus1_intensity = zeros(Float32, 2),
        ms1_mminus1_mass_err_ppm = fill(-1f0, 2),
        ms1_mminus1_to_m0_ratio = fill(-1f0, 2),
    )

    iso_splines = (_, isotope, _) -> isotope == 0 ? 1f0 : 0.3f0

    _ms1_lookup_chunk!(
        psms,
        spectra,
        Float32[500f0, 600f0],
        UInt8[2, 2],
        UInt8[0, 0],
        iso_splines,
        Int32[1, 1],
        10f0,
        0f0,
        neutron,
        1,
        nrow(psms),
    )

    @test psms.ms1_mminus1_present[1] == 1f0
    @test psms.ms1_mminus1_intensity[1] == 80f0
    expected_mminus1_err = abs(observed_mminus1 - target_mminus1) / target_mminus1 * 1f6
    @test psms.ms1_mminus1_mass_err_ppm[1] ≈ expected_mminus1_err rtol = 1f-6
    @test psms.ms1_mminus1_to_m0_ratio[1] == 0.4f0

    @test psms.ms1_mminus1_present[2] == 0f0
    @test psms.ms1_mminus1_intensity[2] == 0f0
    @test psms.ms1_mminus1_mass_err_ppm[2] == -1f0
    @test psms.ms1_mminus1_to_m0_ratio[2] == -1f0
end

@testset "LightGBM M-1 feature registration" begin
    mminus1_features = [
        :ms1_mminus1_present,
        :ms1_mminus1_intensity,
        :ms1_mminus1_mass_err_ppm,
        :ms1_mminus1_to_m0_ratio,
    ]

    @test all(feature -> feature in PRESCORE_FEATURES, mminus1_features)
    @test all(feature -> feature in ADVANCED_FEATURE_SET, mminus1_features)

    X = feature_matrix(
        DataFrame(
            ms1_mminus1_present = Float32[1],
            ms1_mminus1_intensity = Float32[80],
            ms1_mminus1_mass_err_ppm = Float32[2],
            ms1_mminus1_to_m0_ratio = Float32[0.4],
        ),
        mminus1_features,
    )
    @test vec(X) == Float32[1, 80, 2, 0.4]
end
