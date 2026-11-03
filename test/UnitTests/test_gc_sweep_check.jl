# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

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
