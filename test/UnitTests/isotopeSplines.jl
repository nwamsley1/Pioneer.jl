
iso_splines = parseIsoXML(joinpath(dirname(dirname(@__DIR__)),"data", "IsotopeSplines", "IsotopeSplines_10kDa_21isotopes-1.xml"));
SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/ASMS_2024/ASTRAL_THREE_PROTEOME/unispec_chronologer_1mc_1var_by_052724/spec_lib/pioneer_lib/"
 
#Use exact isotope composition from expasy.org 
#https://education.expasy.org/student_projects/isotopident/cgi-bin/iso.pl
@testset "isotope_spline_accuracy" begin
    pep_a = "PEPTIDEPEPTIDEMMM"
    pep_mass_a = 1973.831f0
    pep_sulfur_count = 4 #1 means zero for for 3 sulfurs must be 4
    exact_iso_abundances = Float32[0.953933484,1.000000000,0.707123653,0.373977793,0.162691088,0.060424297]
    spline_iso_abundances = [iso_splines.splines[pep_sulfur_count][x](pep_mass_a) for x in range(1, 6)]
    spline_iso_abundances = spline_iso_abundances./maximum(spline_iso_abundances)
    percent_deviance = 100.0*(exact_iso_abundances .- spline_iso_abundances)./exact_iso_abundances
    #within 3% on first three isotopes
    @test all(abs.(percent_deviance[1:3]) .< 3)

    pep_a = "PEPTIDEPEPTIDE"
    pep_mass_a = 1580.709f0
    pep_sulfur_count = 1
    exact_iso_abundances = Float32[1.000000000,0.840226254,0.406760020,0.144049251,0.041139503,0.009968661]
    spline_iso_abundances = [iso_splines.splines[pep_sulfur_count][x](pep_mass_a) for x in range(1, 6)]
    spline_iso_abundances = spline_iso_abundances./maximum(spline_iso_abundances)
    percent_deviance = 100.0*(exact_iso_abundances .- spline_iso_abundances)./exact_iso_abundances
    #within 3% on first three isotopes
    @test all(abs.(percent_deviance[1:3]) .< 3)

    pep_a = "VLSMSNLR"
    pep_mass_a = 918.496f0
    pep_sulfur_count = 2
    exact_iso_abundances = Float32[1.000000000,0.493068635,0.187476610,0.052138016,0.011652472,0.002172227]
    spline_iso_abundances = [iso_splines.splines[pep_sulfur_count][x](pep_mass_a) for x in range(1, 6)]
    spline_iso_abundances = spline_iso_abundances./maximum(spline_iso_abundances)
    percent_deviance = 100.0*(exact_iso_abundances .- spline_iso_abundances)./exact_iso_abundances
    #within 3% on first three isotopes
    @test all(abs.(percent_deviance[1:3]) .< 3)

    pep_a = "NNSPDGNTDSSALDCYNPMTNQWSPCAPMSVPRNR"
    pep_mass_a = 3952.662f0
    pep_sulfur_count = 5
    exact_iso_abundances = Float32[0.411367694,0.852177261,1.000000000,0.851991166,0.581595404,0.335087683]
    spline_iso_abundances = [iso_splines.splines[pep_sulfur_count][x](pep_mass_a) for x in range(1, 6)]
    spline_iso_abundances = spline_iso_abundances./maximum(spline_iso_abundances)
    percent_deviance = 100.0*(exact_iso_abundances .- spline_iso_abundances)./exact_iso_abundances
    #within 2% on first three isotopes
    @test all(abs.(percent_deviance[1:3]) .< 3)
end
#########
#Precursor isotope set testing 
#########
@testset "precursor_isotope_sets" begin

    @test (1, 1) == getPrecursorIsotopeSet(
        400.0f0, #precursor m/z
        UInt8(2), #precursor charge state
        400.25f0, #minimum quadrupole m/z
        400.75f0, #maximum quadrupole m/z
        max_iso = 5
    )

    @test (0, 1) == getPrecursorIsotopeSet(
        400.0f0, #precursor m/z
        UInt8(2), #precursor charge state
        399.75f0, #minimum quadrupole m/z
        400.75f0, #maximum quadrupole m/z
        max_iso = 5
    )

    @test (0, 0) == getPrecursorIsotopeSet(
        400.0f0, #precursor m/z
        UInt8(2), #precursor charge state
        399.75f0, #minimum quadrupole m/z
        400.25f0, #maximum quadrupole m/z
        max_iso = 5
    )

    @test (0, 2) == getPrecursorIsotopeSet(
        400.0f0, #precursor m/z
        UInt8(3), #precursor charge state
        399.75f0, #minimum quadrupole m/z
        400.75f0, #maximum quadrupole m/z
        max_iso = 5
    )

    @test (3, 3) == getPrecursorIsotopeSet(
        400.0f0, #precursor m/z
        UInt8(3), #precursor charge state
        400.75f0, #minimum quadrupole m/z
        401.25f0, #maximum quadrupole m/z
        max_iso = 5
    )

    #No isotopes in the range given 
    @test (-1, 0) == getPrecursorIsotopeSet(
        400.0f0, #precursor m/z
        UInt8(3), #precursor charge state
        400.75f0, #minimum quadrupole m/z
        400.80f0, #maximum quadrupole m/z
        max_iso = 5
    )

end

#########
#Fragment Isotope Correction
#########

function plotIsotopes(
    iso_a,
    iso_b,
    label_a,
    label_b,
    test_frag;
    title::String = "title here"
    )

    frag_charge = test_frag.frag_charge
    p = plot(title = title)
    iso_idx = 0
    for i in range(1, length(iso_a))
        mz = test_frag.mz + iso_idx*NEUTRON/frag_charge
        plot!(p, [mz, mz], [0.0, iso_a[i]], color = 1, alpha = 0.5, lw = 5, label = nothing)
        iso_idx += 1
    end
    hline!(p,[0.0], lw = 4, color = 1, labels = label_a)
    iso_idx = 0
    for i in range(1, length(iso_b))
        mz = test_frag.mz + iso_idx*NEUTRON/frag_charge
        plot!(p, [mz, mz], [0.0, iso_b[i]], color = 2, alpha = 0.5, lw = 5, label = nothing)
        iso_idx += 1
    end
    hline!(p,[0.0], lw = 4, color = 2, labels = label_b)
    plot!(p, show = true)
end

test_precursor = (mz = 981.5947f0,
 sequence = "TFTLKTVLMIAIQLITR",
 prec_charge = 0x02,
 sulfur_count = 0x01
)
test_fragments = [
    DetailedFrag{Float32}(0x0087b1fe, 1371.8392f0, Float16(1.0), 0x02, false, 0x01, 0x0c, 0x02, 0x01, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 350.171f0, Float16(0.919), 0x01, false, 0x01, 0x03, 0x02, 0x02, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 1171.7231f0, Float16(0.7295), 0x02, false, 0x01, 0x0a, 0x02, 0x03, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 1058.6392f0, Float16(0.531), 0x02, false, 0x01, 0x09, 0x02, 0x04, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 814.5145f0, Float16(0.4707), 0x02, false, 0x01, 0x07, 0x02, 0x05, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 1270.7915f0, Float16(0.4026), 0x02, false, 0x01, 0x0b, 0x02, 0x06, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 249.1234f0, Float16(0.3818), 0x01, false, 0x01, 0x02, 0x02, 0x07, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 927.5986f0, Float16(0.357), 0x02, false, 0x01, 0x08, 0x02, 0x08, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 692.3978f0, Float16(0.3445), 0x01, false, 0x01, 0x06, 0x02, 0x09, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 1499.9342f0, Float16(0.2551), 0x02, false, 0x01, 0x0d, 0x02, 0x0a, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 791.4662f0, Float16(0.1786), 0x01, false, 0x01, 0x07, 0x02, 0x0b, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 463.2551f0, Float16(0.1737), 0x01, false, 0x01, 0x04, 0x02, 0x0c, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 743.4774f0, Float16(0.1628), 0x02, false, 0x01, 0x06, 0x02, 0x0d, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 630.3933f0, Float16(0.1622), 0x02, false, 0x01, 0x05, 0x02, 0x0e, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 591.3501f0, Float16(0.1403), 0x01, false, 0x01, 0x05, 0x02, 0x0f, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 750.4707f0, Float16(0.1305), 0x02, false, 0x02, 0x0d, 0x02, 0x10, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 904.5502f0, Float16(0.0885), 0x01, false, 0x01, 0x08, 0x02, 0x11, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 502.3348f0, Float16(0.0733), 0x02, false, 0x01, 0x04, 0x02, 0x12, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 389.2507f0, Float16(0.0686), 0x02, false, 0x01, 0x03, 0x02, 0x13, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 857.5366f0, Float16(0.057), 0x02, false, 0x02, 0x0f, 0x02, 0x14, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 1035.5907f0, Float16(0.0509), 0x01, false, 0x01, 0x09, 0x02, 0x15, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 807.0128f0, Float16(0.0375), 0x02, false, 0x02, 0x0e, 0x02, 0x16, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 1148.6748f0, Float16(0.0337), 0x01, false, 0x01, 0x0a, 0x02, 0x17, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 1613.0183f0, Float16(0.0222), 0x02, false, 0x01, 0x0e, 0x02, 0x18, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 1219.7119f0, Float16(0.0195), 0x01, false, 0x01, 0x0b, 0x02, 0x19, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 1714.0659f0, Float16(0.0133), 0x02, false, 0x01, 0x0f, 0x02, 0x1a, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 1332.796f0, Float16(0.0109), 0x01, false, 0x01, 0x0c, 0x02, 0x1b, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 931.0708f0, Float16(0.0095), 0x02, false, 0x02, 0x10, 0x02, 0x1c, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 610.3596f0, Float16(0.0082), 0x01, false, 0x02, 0x0b, 0x02, 0x1d, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 346.7025f0, Float16(0.0063), 0x01, false, 0x02, 0x06, 0x02, 0x1e, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 452.7788f0, Float16(0.0062), 0x01, false, 0x02, 0x08, 0x02, 0x1f, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 1573.9386f0, Float16(0.006), 0x01, false, 0x01, 0x0e, 0x02, 0x20, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 666.9016f0, Float16(0.0059), 0x01, false, 0x02, 0x0c, 0x02, 0x21, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 1460.8545f0, Float16(0.0055), 0x01, false, 0x01, 0x0d, 0x02, 0x22, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 396.2367f0, Float16(0.0055), 0x01, false, 0x02, 0x07, 0x02, 0x23, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 730.9309f0, Float16(0.0054), 0x01, false, 0x02, 0x0d, 0x02, 0x24, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 1687.0227f0, Float16(0.0051), 0x01, false, 0x01, 0x0f, 0x02, 0x25, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 787.4729f0, Float16(0.0051), 0x01, false, 0x02, 0x0e, 0x02, 0x26, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 574.841f0, Float16(0.005), 0x01, false, 0x02, 0x0a, 0x02, 0x27, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 894.5388f0, Float16(0.0047), 0x01, false, 0x02, 0x10, 0x02, 0x28, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 518.299f0, Float16(0.0046), 0x01, false, 0x02, 0x09, 0x02, 0x29, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 844.015f0, Float16(0.0044), 0x01, false, 0x02, 0x0f, 0x02, 0x2a, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 686.4233f0, Float16(0.0043), 0x02, false, 0x02, 0x0c, 0x02, 0x2b, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 529.8232f0, Float16(0.0038), 0x02, false, 0x02, 0x09, 0x02, 0x2c, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 407.7609f0, Float16(0.0038), 0x02, false, 0x02, 0x07, 0x02, 0x2d, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 296.1787f0, Float16(0.0037), 0x01, false, 0x02, 0x05, 0x02, 0x2e, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 464.3029f0, Float16(0.0036), 0x02, false, 0x02, 0x08, 0x02, 0x2f, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 1788.0703f0, Float16(0.0036), 0x01, false, 0x01, 0x10, 0x02, 0x30, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 586.3652f0, Float16(0.0035), 0x02, false, 0x02, 0x0a, 0x02, 0x31, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 372.2423f0, Float16(0.0035), 0x02, false, 0x02, 0x06, 0x02, 0x32, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 195.129f0, Float16(0.0034), 0x02, false, 0x02, 0x03, 0x02, 0x33, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 232.1312f0, Float16(0.0031), 0x01, false, 0x02, 0x04, 0x02, 0x34, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 315.7003f0, Float16(0.0031), 0x02, false, 0x02, 0x05, 0x02, 0x35, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 1861.1343f0, Float16(0.0031), 0x02, false, 0x01, 0x10, 0x02, 0x36, 0x01)
    DetailedFrag{Float32}(0x0087b1fe, 251.671f0, Float16(0.0027), 0x02, false, 0x02, 0x04, 0x02, 0x37, 0x00)
    DetailedFrag{Float32}(0x0087b1fe, 635.8994f0, Float16(0.0026), 0x02, false, 0x02, 0x0b, 0x02, 0x38, 0x01)
]

@testset "fragment_isotope_correction" begin 
    #If only the precursor monoisotope is sampled,
    #then the fragment should be all mono
    isotopes = zeros(Float32, 5)
    test_frag = test_fragments[1]
    getFragIsotopes!(
        isotopes,
        iso_splines,
        test_precursor[:mz],
        test_precursor[:prec_charge],
        test_precursor[:sulfur_count],
        test_frag,
        (0, 0)
    )
    isotopes = isotopes./sum(isotopes)
    @test all((isotopes .- [1.0, 0.0, 0.0, 0.0, 0.0]).<1e-6)
    #Repeat with quad transmission function 
    fill!(isotopes, zero(Float32))
    precursor_transmission = Float32[1, 0, 0, 0, 0]
    getFragIsotopes!(
        isotopes,
        precursor_transmission,
        iso_splines,
        test_precursor[:mz],
        test_precursor[:prec_charge],
        test_precursor[:sulfur_count],
        test_frag
    )
    isotopes = isotopes./sum(isotopes)
    @test all((isotopes .- [1.0, 0.0, 0.0, 0.0, 0.0]).<1e-6)

    isotopes = zeros(Float32, 5)
    test_frag = DetailedFrag{Float32}(0x0087b1fe, 1371.8392f0, Float16(1.0), 0x02, false, 0x01, 0x0c, 0x02, 0x01, 0x01)
    getFragIsotopes!(
        isotopes,
        iso_splines,
        test_precursor[:mz],
        test_precursor[:prec_charge],
        test_precursor[:sulfur_count],
        test_frag,
        (0, 5)
    )
    isotopes05 = isotopes./sum(isotopes)

    fill!(isotopes, zero(Float32))
    precursor_transmission = Float32[1, 1, 1, 1, 1]
    getFragIsotopes!(
        isotopes,
        precursor_transmission,
        iso_splines,
        test_precursor[:mz],
        test_precursor[:prec_charge],
        test_precursor[:sulfur_count],
        test_frag
    )
    isotopes = isotopes./sum(isotopes)
    @test all(abs.(isotopes .- isotopes05) .< 1e-2)
    #@test all((isotopes .- [1.0, 0.0, 0.0, 0.0, 0.0]).<1e-6)
    isotopes = zeros(Float32, 5)
    getFragIsotopes!(
        isotopes,
        iso_splines,
        test_precursor[:mz],
        test_precursor[:prec_charge],
        test_precursor[:sulfur_count],
        test_frag,
        (1, 5)
    )
    isotopes15 = isotopes./sum(isotopes)


    fill!(isotopes, zero(Float32))
    precursor_transmission = Float32[0, 1, 1, 1, 1]
    getFragIsotopes!(
        isotopes,
        precursor_transmission,
        iso_splines,
        test_precursor[:mz],
        test_precursor[:prec_charge],
        test_precursor[:sulfur_count],
        test_frag
    )
    isotopes = isotopes./sum(isotopes)
    @test all(abs.(isotopes .- isotopes15) .< 1e-2)

    @test isotopes05[1]>isotopes15[1]
    @test all(isotopes15[2:end].>isotopes05[2:end])

    plotIsotopes(isotopes05, isotopes15, "(0, 5)", "(1, 5)", test_frag; title = "y12+1 of TFTLKTVLMIAIQLITR")


    isotopes = zeros(Float32, 5)
    test_frag = DetailedFrag{Float32}(0x0087b1fe, 350.171f0, Float16(0.919), 0x01, false, 0x01, 0x03, 0x02, 0x02, 0x00)
    getFragIsotopes!(
        isotopes,
        iso_splines,
        test_precursor[:mz],
        test_precursor[:prec_charge],
        test_precursor[:sulfur_count],
        test_frag,
        (0, 5)
    )
    isotopes05 = isotopes./sum(isotopes)


    fill!(isotopes, zero(Float32))
    precursor_transmission = Float32[1, 1, 1, 1, 1]
    getFragIsotopes!(
        isotopes,
        precursor_transmission,
        iso_splines,
        test_precursor[:mz],
        test_precursor[:prec_charge],
        test_precursor[:sulfur_count],
        test_frag
    )
    isotopes = isotopes./sum(isotopes)
    @test all(abs.(isotopes .- isotopes05) .< 1e-2)

    #@test all((isotopes .- [1.0, 0.0, 0.0, 0.0, 0.0]).<1e-6)
    isotopes = zeros(Float32, 5)
    getFragIsotopes!(
        isotopes,
        iso_splines,
        test_precursor[:mz],
        test_precursor[:prec_charge],
        test_precursor[:sulfur_count],
        test_frag,
        (1, 5)
    )
    isotopes15 = isotopes./sum(isotopes)


    fill!(isotopes, zero(Float32))
    precursor_transmission = Float32[0, 1, 1, 1, 1]
    getFragIsotopes!(
        isotopes,
        precursor_transmission,
        iso_splines,
        test_precursor[:mz],
        test_precursor[:prec_charge],
        test_precursor[:sulfur_count],
        test_frag
    )
    isotopes = isotopes./sum(isotopes)
    @test all(abs.(isotopes .- isotopes15) .< 1e-2)

    @test isotopes05[1]>isotopes15[1]
    @test all(isotopes15[2:end].>isotopes05[2:end])

    plotIsotopes(isotopes05, isotopes15, "(0, 5)", "(1, 5)", test_frag; title = "b3+1 of TFTLKTVLMIAIQLITR")
    b3_isotopes15 = copy(isotopes15)

    isotopes = zeros(Float32, 5)
    test_frag = DetailedFrag{Float32}(0x0087b1fe, 1861.1343f0, Float16(0.0031), 0x02, false, 0x01, 0x10, 0x02, 0x36, 0x01)
    getFragIsotopes!(
        isotopes,
        iso_splines,
        test_precursor[:mz],
        test_precursor[:prec_charge],
        test_precursor[:sulfur_count],
        test_frag,
        (0, 5)
    )
    isotopes05 = isotopes./sum(isotopes)


    fill!(isotopes, zero(Float32))
    precursor_transmission = Float32[1, 1, 1, 1, 1]
    getFragIsotopes!(
        isotopes,
        precursor_transmission,
        iso_splines,
        test_precursor[:mz],
        test_precursor[:prec_charge],
        test_precursor[:sulfur_count],
        test_frag
    )
    isotopes = isotopes./sum(isotopes)
    @test all(abs.(isotopes .- isotopes05) .< 1e-2)

    #@test all((isotopes .- [1.0, 0.0, 0.0, 0.0, 0.0]).<1e-6)
    isotopes = zeros(Float32, 5)
    getFragIsotopes!(
        isotopes,
        iso_splines,
        test_precursor[:mz],
        test_precursor[:prec_charge],
        test_precursor[:sulfur_count],
        test_frag,
        (1, 5)
    )
    isotopes15 = isotopes./sum(isotopes)

    fill!(isotopes, zero(Float32))
    precursor_transmission = Float32[0, 1, 1, 1, 1]
    getFragIsotopes!(
        isotopes,
        precursor_transmission,
        iso_splines,
        test_precursor[:mz],
        test_precursor[:prec_charge],
        test_precursor[:sulfur_count],
        test_frag
    )
    isotopes = isotopes./sum(isotopes)
    @test all(abs.(isotopes .- isotopes15) .< 1e-2)

    @test isotopes05[1]>isotopes15[1]
    @test all(isotopes15[2:end].>isotopes05[2:end])
    @test isotopes15[1] < b3_isotopes15[1]
    plotIsotopes(isotopes05, isotopes15, "(0, 5)", "(1, 5)", test_frag; title = "y16+1 of TFTLKTVLMIAIQLITR")
end

