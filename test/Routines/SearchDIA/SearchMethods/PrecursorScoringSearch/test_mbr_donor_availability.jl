using Arrow
using DataFrames
using Pioneer
using SentinelArrays
using Test

@testset "Bounded pre-integration MBR donor availability" begin
    function write_donor_batches(path, rows)
        midpoint = nrow(rows) ÷ 2
        open(Arrow.Writer, path; file=true) do writer
            Arrow.write(writer, rows[1:midpoint, :])
            Arrow.write(writer, rows[midpoint+1:end, :])
        end
        @test Arrow.Table(path).precursor_idx isa SentinelArrays.ChainedVector
        return path
    end

    function qualifying_donor_runs(rows, score_floor)
        reference = Dict{UInt32, Set{UInt32}}()
        for row in eachrow(rows)
            row.trace_prob_prepass >= score_floor || continue
            push!(get!(Set{UInt32}, reference, row.precursor_idx), row.ms_file_idx)
        end
        return reference
    end

    @testset "Eligibility matches complete donor sets" begin
        mktempdir() do dir
            rows = DataFrame(
                precursor_idx = UInt32[1, 1, 1, 2, 2, 2, 3, 3, 4, 5, 5, 5, 6, 1, 2, 2, 3],
                ms_file_idx = UInt16[10, 10, 10, 10, 20, 30, 10, 20, 40, 50, 60, 70, 60, 10, 10, 40, 30],
                trace_prob_prepass = Float32[
                    0.4, 0.5, 0.9, 0.7, 0.6, 0.8, 0.49, NaN, 0.5,
                    Inf, -Inf, NaN, 0.95, 0.8, 0.99, 0.75, -Inf,
                ],
                target = Bool[true, true, true, true, true, true, true, true, true,
                    true, true, true, false, true, true, true, true],
            )
            paths = [
                write_donor_batches(joinpath(dir, "first.arrow"), rows[1:8, :]),
                write_donor_batches(joinpath(dir, "second.arrow"), rows[9:13, :]),
                write_donor_batches(joinpath(dir, "without_target.arrow"),
                    select(rows[14:end, :], Not(:target))),
            ]
            empty_path = joinpath(dir, "empty.arrow")
            Arrow.write(empty_path, rows[1:0, :])
            push!(paths, empty_path)

            for floor in (0.5f0, Inf32, -Inf32, NaN32)
                donors = Pioneer._mbr_preintegration_donor_files(paths, floor)
                reference = qualifying_donor_runs(rows, floor)
                @test donors isa Dict{UInt32, Tuple{UInt32, UInt32}}
                @test Set(keys(donors)) == Set(keys(reference))
                @test all(donors) do (pid, pair)
                    pair[1] in reference[pid] && pair[2] in reference[pid] &&
                        (pair[1] == pair[2]) == (length(reference[pid]) == 1)
                end
                receivers = UInt32[0, 10, 20, 30, 40, 50, 60, 70, 99]
                actual = [Pioneer._mbr_has_cross_run_donor(donors, pid, receiver)
                    for pid in UInt32.(1:7), receiver in receivers]
                expected = [any(!=(receiver), get(reference, pid, Set{UInt32}()))
                    for pid in UInt32.(1:7), receiver in receivers]
                @test actual == expected
            end

            donors = Pioneer._mbr_preintegration_donor_files(paths, 0.5f0)
            @test donors[UInt32(1)] == (UInt32(10), UInt32(10))
            @test !haskey(donors, UInt32(3))
            @test donors[UInt32(4)] == (UInt32(40), UInt32(40))
            @test donors[UInt32(5)] == (UInt32(50), UInt32(50))
            @test Pioneer._mbr_has_cross_run_donor(donors, UInt32(6), UInt32(99))
            @test !Pioneer._mbr_has_cross_run_donor(donors, UInt32(7), UInt32(99))
            @test isempty(Pioneer._mbr_preintegration_donor_files(String[], 0.5f0))
        end
    end

    @testset "6,000 donor runs retain only two run IDs per precursor" begin
        mktempdir() do dir
            rows = DataFrame(
                precursor_idx = fill(UInt32(1), 6_000),
                ms_file_idx = UInt32.(1:6_000),
                trace_prob_prepass = fill(0.8f0, 6_000),
            )
            paths = [
                write_donor_batches(joinpath(dir, "first.arrow"), rows[1:3_000, :]),
                write_donor_batches(joinpath(dir, "second.arrow"), rows[3_001:end, :]),
            ]
            donors = Pioneer._mbr_preintegration_donor_files(paths, 0.5f0)

            @test donors isa Dict{UInt32, Tuple{UInt32, UInt32}}
            @test length(donors) == 1
            pair = donors[UInt32(1)]
            @test length(pair) == 2
            @test pair[1] != pair[2]
            @test all(run -> UInt32(1) <= run <= UInt32(6_000), pair)
            @test all(UInt32.(1:6_001)) do receiver
                Pioneer._mbr_has_cross_run_donor(donors, UInt32(1), receiver)
            end
            @test !Pioneer._mbr_has_cross_run_donor(donors, UInt32(2), UInt32(1))
        end
    end
end
