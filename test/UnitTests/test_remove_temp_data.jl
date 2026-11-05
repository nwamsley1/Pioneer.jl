# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

using Test
using Pioneer

@testset "End-of-search temp_data removal" begin
    function make_temp_data(root)
        temp = joinpath(root, "temp_data")
        for d in ("jl_AbCdEf", "main_search_psms", "passing_psms", "rt_indices")
            mkpath(joinpath(temp, d))
            write(joinpath(temp, d, "part.arrow"), rand(UInt8, 64))
        end
        write(joinpath(temp, "bitvec_model_global.arrow"), rand(UInt8, 16))
        return temp
    end

    @testset "Open score calibrations are closed and their workspace deleted" begin
        mktempdir() do root
            fit = Pioneer.build_score_calibration(; memory_budget_bytes=4096, temp_parent=root) do emit
                for i in 1:1000
                    emit(Float32(i), isodd(i))
                end
            end
            store = fit.qval_spline.store
            @test haskey(Pioneer.OPEN_SCORE_CALIBRATIONS, store)
            @test isopen(store.io) && isdir(dirname(store.path))
            Pioneer.close_open_score_calibrations!()
            @test !isopen(store.io)
            @test !ispath(dirname(store.path))
            @test !haskey(Pioneer.OPEN_SCORE_CALIBRATIONS, store)
            # Closing again, directly or through the finalizer, is harmless.
            Pioneer._close_score_calibration(fit)
            finalize(store)
        end
    end

    @testset "Removes everything, including a calibration workspace inside temp_data" begin
        mktempdir() do root
            temp = make_temp_data(root)
            fit = Pioneer.build_score_calibration(; memory_budget_bytes=4096, temp_parent=temp) do emit
                for i in 1:1000
                    emit(Float32(i), isodd(i))
                end
            end
            @test Pioneer.remove_temp_data(temp) == String[]
            @test !ispath(temp)
            @test !isopen(fit.qval_spline.store.io)
        end
    end

    @testset "One entry that cannot be deleted does not stop the rest" begin
        # Stands in for NFS, where an open file inside jl_* leaves the directory non-empty.
        mktempdir() do root
            temp = make_temp_data(root)
            stuck = joinpath(temp, "jl_AbCdEf")
            attempts = Ref(0)
            rm_entry = p -> begin
                p == stuck && (attempts[] += 1; throw(Base.IOError("rm(\"$p\"): directory not empty (ENOTEMPTY)", -39)))
                rm(p; recursive=true, force=true)
            end
            @test Pioneer.remove_temp_data(temp; rm_entry) == [stuck]
            @test attempts[] == 2                      # retried once after GC
            @test readdir(temp) == ["jl_AbCdEf"]       # every other entry is gone
        end
    end

    @testset "An entry that frees up on the retry is removed" begin
        mktempdir() do root
            temp = make_temp_data(root)
            stuck = joinpath(temp, "jl_AbCdEf")
            attempts = Ref(0)
            rm_entry = p -> begin
                p == stuck && (attempts[] += 1) == 1 && throw(Base.IOError("busy", -16))
                rm(p; recursive=true, force=true)
            end
            @test Pioneer.remove_temp_data(temp; rm_entry) == String[]
            @test !ispath(temp)
        end
    end

    @testset "Missing folder" begin
        @test Pioneer.remove_temp_data(joinpath(tempdir(), "pioneer_no_such_temp_data")) == String[]
    end
end
