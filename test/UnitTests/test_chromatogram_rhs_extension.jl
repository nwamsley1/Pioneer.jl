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

using Test, Random
using Pioneer: DensePrecMap, SparsePrecMap, accumulate_max_weight!, select_rhs_precursors,
               rt_and_rhs, CHROM_RHS_MULT

@testset "chromatogram RHS extension selection" begin
    rng = MersenneTwister(7)
    n = 20_000
    files = [(rand(rng, UInt32(1):UInt32(n), 2000),
              [rand(rng) < 0.02 ? 0f0 : Float32(exp(10 + 3randn(rng))) for _ in 1:2000]) for _ in 1:8]
    ref = Dict{UInt32, Float32}()   # max weight per precursor over every row of every file
    for (p, w) in files, i in eachindex(p)
        ref[p[i]] = max(get(ref, p[i], -Inf32), w[i])
    end
    vals = collect(values(ref))
    thr = partialsort!(copy(vals), ceil(Int, 0.95 * length(vals)))
    want = Set(p for (p, v) in ref if v >= thr)
    for m in (DensePrecMap{Float32}(n), SparsePrecMap{Float32}())
        for (p, w) in files
            accumulate_max_weight!(m, p, w)
        end
        selected, n_obs, _ = select_rhs_precursors(m, 0.05f0)
        @test n_obs == length(ref)          # zero-weight precursors stay in the pool
        @test selected == want
    end

    @test rt_and_rhs(Dict{UInt32, Float32}(UInt32(3) => 1.5f0), UInt32(3)) == (1.5f0, 1f0)
    m2 = Dict{UInt32, NTuple{2, Float32}}(UInt32(3) => (1.5f0, CHROM_RHS_MULT))
    @test rt_and_rhs(m2, UInt32(3)) == (1.5f0, CHROM_RHS_MULT)
    @test isnan(first(rt_and_rhs(m2, UInt32(4))))
end
