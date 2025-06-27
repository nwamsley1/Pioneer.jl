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

#=

Three peptides, a, b, and c. Sould have some fragments matching ot the same peak, 
Should look something like this 
3×3 Matrix{Float64}:
 1.0  1.0  0.0
 0.0  1.0  1.0
 2.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0
=#
@testset "build_design_matrix.jl" begin
    H = SparseArray(UInt32(10))
    test_id_to_col = ArrayDict(UInt32, UInt16, 3)
    test_matches = sort(
    [
    #Fragments for peptide A
    FragmentMatch{Float32}(1.0, 1.0, 100.0, 100.0, 1, 0x01, 0x01, 0x00, 0x01, true, 0x00000001, 0x01, 0x00000001, 0x00000000, 0x01)
    FragmentMatch{Float32}(1.0, 1.0, 100.0, 100.0, 3, 0x01, 0x01, 0x00, 0x01, true, 0x00000001, 0x01, 0x00000001, 0x00000000, 0x01)
    FragmentMatch{Float32}(1.0, 1.0, 100.0, 100.0, 3, 0x01, 0x01, 0x00, 0x01, true, 0x00000001, 0x01, 0x00000001, 0x00000000, 0x01)
    #Fragments for peptide B
    FragmentMatch{Float32}(1.0, 1.0, 100.0, 100.0, 1, 0x01, 0x01, 0x00, 0x01, true, 0x00000002, 0x01, 0x00000001, 0x00000000, 0x01)
    FragmentMatch{Float32}(1.0, 1.0, 100.0, 100.0, 2, 0x01, 0x01, 0x00, 0x01, true, 0x00000002, 0x01, 0x00000001, 0x00000000, 0x01)
    FragmentMatch{Float32}(1.0, 1.0, 100.0, 100.0, 4, 0x01, 0x01, 0x00, 0x01, true, 0x00000002, 0x01, 0x00000001, 0x00000000, 0x01)
    #Fragments for peptide C 
    FragmentMatch{Float32}(1.0, 1.0, 100.0, 100.0, 2, 0x01, 0x01, 0x00, 0x01, true, 0x00000003, 0x01, 0x00000001, 0x00000000, 0x01)
    FragmentMatch{Float32}(1.0, 1.0, 100.0, 100.0, 5, 0x01, 0x01, 0x00, 0x01, true, 0x00000003, 0x01, 0x00000001, 0x00000000, 0x01)
    ],
    by = x->x.peak_ind)
    test_misses = [
    FragmentMatch{Float32}(1.0, 1.0, 100.0, 100.0, 6, 0x01, 0x01, 0x00, 0x01, true, 0x00000001, 0x01, 0x00000001, 0x00000000, 0x01)
    FragmentMatch{Float32}(1.0, 1.0, 100.0, 100.0, 7, 0x01, 0x01, 0x00, 0x01, true, 0x00000002, 0x01, 0x00000001, 0x00000000, 0x01)
    FragmentMatch{Float32}(1.0, 1.0, 100.0, 100.0, 8, 0x01, 0x01, 0x00, 0x01, true, 0x00000002, 0x01, 0x00000001, 0x00000000, 0x01)
    ]##Vector{FragmentMatch{Float32}}()
    buildDesignMatrix!(
        H,
        test_matches,
        test_misses,
        length(test_matches),
        length(test_misses),
        test_id_to_col 
    )
    @test H.n_vals == 10
    @test H.n == 3
    @test H.m == 8
    H = Matrix(sparse(
        H.rowval[1:H.n_vals],
        H.colval[1:H.n_vals],
        H.nzval[1:H.n_vals]
    ))
    H_test = [
        1.0  1.0  0.0;
        0.0  1.0  1.0;
        2.0  0.0  0.0;
        0.0  1.0  0.0;
        0.0  0.0  1.0;
        1.0  0.0  0.0;
        0.0  1.0  0.0;
        0.0  1.0  0.0
    ]

    @test all(abs.(H.-H_test).<1e-6)
end