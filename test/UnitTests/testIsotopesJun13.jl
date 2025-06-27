# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.


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