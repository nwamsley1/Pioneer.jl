
# plotIsotopes function is defined in isotopeSplines.jl to avoid duplicate definitions

test_precursor = (mz = 981.5947f0,
 sequence = "TFTLKTVLMIAIQLITR",
 prec_charge = 0x02,
 sulfur_count = 0x01
)

isotopes = zeros(Float32, 5)
test_frag = DetailedFrag{Float32}(0x0087b1fe, 1371.8392f0, Float16(1.0), 0x02, false, false, false, false, 0x01, 0x0c, 0x02, 0x01, 0x01)

getFragIsotopes!(
    isotopes,
    iso_splines,
    test_precursor[:mz],
    test_precursor[:prec_charge],
    test_precursor[:sulfur_count],
    test_frag,
    (0, 0)
)
isotopes05 = isotopes./sum(isotopes)
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

@test isotopes05[1]>isotopes15[1]
@test all(isotopes15[2:end].>isotopes05[2:end])

plotIsotopes(isotopes05, isotopes15, "(0, 5)", "(1, 5)", test_frag; title = "y12+1 of TFTLKTVLMIAIQLITR")


isotopes = zeros(Float32, 5)
test_frag = DetailedFrag{Float32}(0x0087b1fe, 350.171f0, Float16(0.919), 0x01, false, false, false, false, 0x01, 0x03, 0x02, 0x02, 0x00)
getFragIsotopes!(
    isotopes,
    iso_splines,
    test_precursor[:mz],
    test_precursor[:prec_charge],
    test_precursor[:sulfur_count],
    test_frag,
    (0, 0)
)
isotopes05 = isotopes./sum(isotopes)
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

@test isotopes05[1]>isotopes15[1]
@test all(isotopes15[2:end].>isotopes05[2:end])

plotIsotopes(isotopes05, isotopes15, "(0, 5)", "(1, 5)", test_frag; title = "b3+1 of TFTLKTVLMIAIQLITR")
b3_isotopes15 = copy(isotopes15)


#correctPrecursorAbundance(100.0f0, iso_splines, (0, 5), test_precursor[:mz]*test_precursor[:prec_charge], test_precursor[:sulfur_count])