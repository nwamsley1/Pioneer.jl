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

@testset "Scoring Interface Tests" begin

    @testset "Logodds Function" begin
        # Test single probability
        @test logodds([0.8], 1) == 0.8

        # Test multiple probabilities
        probs = [0.9, 0.8, 0.7, 0.6]
        result = logodds(probs, 2)  # Take top 2
        @test result isa Float64
        @test 0.0 <= result <= 1.0

        # Test edge cases
        @test 0.0 < logodds([0.0, 0.0], 2) < 1.0  # Should handle zeros
        @test 0.0 < logodds([1.0, 1.0], 2) < 1.0  # Should handle ones

        # Test empty case
        @test logodds(Float64[], 1) == 0.0f0
    end
end
