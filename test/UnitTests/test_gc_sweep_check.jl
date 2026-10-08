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

# BuildSpecLib warns when Julia runs a concurrent GC sweep thread (see gc_concurrent_sweep_enabled).

@testset "gc_concurrent_sweep_enabled" begin
    check(; nmark = 0, nsweep = 0, env = Dict{String, String}()) =
        Pioneer.gc_concurrent_sweep_enabled(; nmark, nsweep, env)
    # --gcthreads given: it decides, and the env var is ignored (as Julia does)
    @test check(nmark = 6, nsweep = 1)
    @test !check(nmark = 6, nsweep = 0)
    @test !check(nmark = 6, nsweep = 0, env = Dict("JULIA_NUM_GC_THREADS" => "6,1"))
    # no flag: JULIA_NUM_GC_THREADS decides
    @test check(env = Dict("JULIA_NUM_GC_THREADS" => "6,1"))
    @test !check(env = Dict("JULIA_NUM_GC_THREADS" => "6,0"))
    @test !check(env = Dict("JULIA_NUM_GC_THREADS" => "6"))
    @test !check()
    @test Pioneer.gc_concurrent_sweep_enabled() isa Bool   # this process, real defaults
end
