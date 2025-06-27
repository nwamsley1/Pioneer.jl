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

using SparseArrays

@testset "buildDesignMatrix.jl" begin

    #Test example. There are three precursors.
    #Each has three matched peaks and one missed peak. 
    #The matched peaks
    matches = [FragmentMatch() for x in 1:9]
    misses = [FragmentMatch() for x in 1:3]

    #Pairs of precursor id's and peak_indices of matched fragments
    matched_frags = [(UInt32(1), 1,1.0, 100.0),
     (UInt32(1), 2, 1.0, 100.0),
     (UInt32(1), 3, 1.0, 100.0), #shared
     (UInt32(2), 4, 1.0, 100.0), #shared
     (UInt32(2), 5, 1.0, 100.0),
     (UInt32(2), 3, 1.0, 100.0), #shared
     (UInt32(3), 4, 1.0, 100.0), #shared
     (UInt32(3), 6, 1.0, 100.0),
     (UInt32(3), 7, 1.0, 100.0)]

    missed_frags = [(UInt32(1), 8, 100.0),
                    (UInt32(2), 9, 100.0),
                    (UInt32(3), 10, 100.0)]

    [matches[i].prec_id = first(x) for (i, x) in enumerate(matched_frags)]
    [matches[i].peak_ind = x[2] for (i, x) in enumerate(matched_frags)]
    [matches[i].intensity = x[3] for (i, x) in enumerate(matched_frags)]
    [matches[i].predicted_intensity = last(x) for (i, x) in enumerate(matched_frags)]

    [misses[i].prec_id = first(x) for (i, x) in enumerate(missed_frags)]
    [misses[i].peak_ind = x[2] for (i, x) in enumerate(missed_frags)]
    [misses[i].predicted_intensity = last(x) for (i, x) in enumerate(missed_frags)]

    sort!(matches, by = x->getPeakInd(x))
    sort!(misses, by = x->getPeakInd(x))
    X, Hs, Hst, IDtoROW, matched_cols = buildDesignMatrix(matches, misses)

    #=
    julia> Matrix(Hs)
    10×3 Matrix{Float64}:
    100.0    0.0    0.0
    100.0    0.0    0.0
    100.0  100.0    0.0 #Third fragment matched by precursor 1 and 2
      0.0  100.0  100.0 #Fourgh fragment matched by precursor 2 and 3
      0.0  100.0    0.0
      0.0    0.0  100.0
      0.0    0.0  100.0
    100.0    0.0    0.0 #Missed fragments
      0.0  100.0    0.0
      0.0    0.0  100.0
    =#
    @test length(X) == length(unique([getPeakInd(x) for x in matches])) + length(unique([getPeakInd(x) for x in misses]))
    Hs_mat = Matrix(Hs);

    @test sum(iszero.(Hs_mat).==false, dims = 2)[:] == [1, 1, 2, 2, 1, 1, 1, 1, 1, 1]

    @test (iszero.(Hs_mat).==false) == BitMatrix([1 0 0; 1 0 0; 1 1 0; 0 1 1; 0 1 0; 0 0 1; 0 0 1; 1 0 0; 0 1 0; 0 0 1])

    @test transpose(Matrix(Hst)) == Hs_mat

    @test all([IDtoROW[i] == UInt32(i) for i in UInt32[1, 2, 3]])

    @test matched_cols == 7

    @test (X.>0.0) == BitVector([1, 1, 1, 1, 1, 1, 1, 0, 0, 0])

end