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

@testset "run summary accumulation" begin
    # Two files; file 2 also receives a decoy row and an unquantified row.
    chunk = (
        ms_file_idx = UInt32[1, 1, 1, 2, 2, 2, 2],
        target = Bool[true, true, true, true, true, true, false],
        sequence = ["PEPTIDEK", "PEPTIDEK", "LLSEQVVDR", "AAAK", "AAAK", "AAAAAAK", "DECOYR"],
        peak_area = Union{Missing, Float32}[10, 30, 20, 5, 0, missing, 99],
        peak_area_normalized = Union{Missing, Float32}[20, 30, 40, 5, 0, missing, 99],
        mbr_recovered = Bool[false, true, false, false, false, true, false],
        irt_error = Float16[-1.0, 2.0, 3.0, 0.5, 0.5, 0.5, 9.0],
        rt_fwhm = Float16[0.1, 0.2, 0.3, 0.4, 0.4, 0.4, 9.0],
        points_integrated = UInt32[3, 5, 7, 4, 4, 4, 1],
        charge = UInt8[2, 3, 2, 1, 1, 1, 2],
        missed_cleavage = UInt8[0, 0, 1, 0, 0, 0, 0],
    )
    @testset "Arrow record batches" begin
        mktempdir() do dir
            for normalized in (false, true), mbr in (false, true)
                rows = DataFrame(chunk)
                normalized || select!(rows, Not(:peak_area_normalized))
                mbr || select!(rows, Not(:mbr_recovered))
                path = joinpath(dir, "summary_$(normalized)_$(mbr).arrow")
                open(Arrow.Writer, path; file=true) do writer
                    Arrow.write(writer, rows[1:3, :])
                    Arrow.write(writer, rows[4:end, :])
                end
                @test length(collect(Arrow.Stream(path))) == 2
                expected = Pioneer.with_run_summary(["a", "b"]) do acc
                    accumulate_run_summary!(acc, rows)
                end
                actual = Pioneer.with_run_summary(["a", "b"]) do acc
                    for batch in Arrow.Stream(path)
                        accumulate_run_summary!(acc, batch)
                    end
                end
                @test all(isequal(getfield(actual[i], field), getfield(expected[i], field))
                          for i in eachindex(expected), field in fieldnames(RunSummaryStats))
            end
        end
    end
    stats = Pioneer.with_run_summary(["a", "b"]) do acc
        accumulate_run_summary!(acc, chunk)
        # A second chunk must fold into the same accumulators.
        accumulate_run_summary!(acc, (
            ms_file_idx = UInt32[1], target = Bool[true], sequence = ["NEWPEPK"],
            peak_area = Union{Missing, Float32}[40], peak_area_normalized = Union{Missing, Float32}[20],
            mbr_recovered = Bool[false], irt_error = Float16[4.0], rt_fwhm = Float16[0.4],
            points_integrated = UInt32[9], charge = UInt8[2], missed_cleavage = UInt8[1]))
    end
    a, b = stats
    @test a.precursors_identified == 4
    @test a.precursors_quantified == 4
    @test a.precursors_mbr == 1
    @test a.peptides_identified == 3            # PEPTIDEK counted once
    @test a.total_peak_area == 100.0
    @test a.medians[1] == 25.0f0
    @test a.medians[2] == 1.5f0   # ratios 2, 1, 2, 0.5
    @test a.medians[3] == 2.5f0             # |−1|,2,3,4
    @test a.medians[7] == 2.0f0
    @test a.medians[8] == 0.5f0
    @test a.medians[6] == 8.0f0        # 8,8,9,7

    @test b.precursors_identified == 3               # decoy skipped
    @test b.precursors_quantified == 1               # zero and missing areas are unquantified
    @test b.precursors_mbr == 1
    @test b.peptides_identified == 2
    @test b.total_peak_area == 5.0
    @test b.medians[2] == 1f0

    pg = DataFrame(
        file_name = ["a", "a", "a", "b", "zzz"],
        target = [true, true, false, true, true],
        abundance = Union{Missing, Float64}[1.0, missing, 3.0, 0.0, 1.0],
    )
    add_protein_group_counts!(stats, pg, ["a", "b"])
    @test a.protein_groups_identified == 2 && a.protein_groups_quantified == 1
    @test b.protein_groups_identified == 1 && b.protein_groups_quantified == 0

    @test Pioneer._median_or_missing(Float32[]) === missing
    @test Pioneer._mass_tol_unit(Pioneer.MassErrorModel(0.0f0, (10.0f0, 12.0f0))) == "ppm"
    @test Pioneer._mass_tol_columns(Dict{Int64, Pioneer.AbstractMassErrorModel}(), 1) === (missing, missing, missing)
    d = Dict{Int64, Pioneer.AbstractMassErrorModel}(1 => Pioneer.MassErrorModel(0.0f0, (10.0f0, 12.0f0)))
    @test Pioneer._mass_tol_columns(d, 1) == (10.0f0, 12.0f0, "ppm")
end

@testset "bounded exact summary partitions" begin
    mktempdir() do dir
        n = 12000
        rows = (
            ms_file_idx = UInt32[isodd(i) ? 1 : 17 for i in 1:n],
            target = trues(n),
            sequence = [string("PEPTIDE", i % 31) for i in 1:n],
            peak_area = Union{Missing,Float32}[i % 11 == 0 ? missing : i % 97 for i in 1:n],
            peak_area_normalized = Float32[i % 107 for i in 1:n],
            irt_error = Float32[sin(i) for i in 1:n],
            rt_fwhm = Float32[isodd(i) ? -0.0 : 0.0 for i in 1:n],
            points_integrated = UInt32[i % 19 for i in 1:n],
            charge = UInt8[i % 4 for i in 1:n],
            missed_cleavage = UInt8[i % 3 for i in 1:n],
        )
        names = string.(1:17)
        small = Pioneer.with_run_summary(names; temp_parent=dir, memory_budget_bytes=65536) do acc
            accumulate_run_summary!(acc, rows)
            @test acc.directory !== nothing
            @test isempty(acc.records)
        end
        large = Pioneer.with_run_summary(names; temp_parent=dir) do acc
            accumulate_run_summary!(acc, rows)
            @test acc.directory === nothing
            @test isempty(acc.streams)
            @test isempty(readdir(dir))
        end
        @test all(isequal(getfield(small[i], f), getfield(large[i], f))
                  for i in 1:17, f in fieldnames(RunSummaryStats))
        @test isempty(readdir(dir))
        for budget in (65536, 64*1024^2)
            @test_throws ErrorException Pioneer.with_run_summary(names; temp_parent=dir, memory_budget_bytes=budget) do acc
                accumulate_run_summary!(acc, rows)
                error("interrupted")
            end
            @test isempty(readdir(dir))
        end
        boundary = Pioneer.with_run_summary(names; temp_parent=dir, memory_budget_bytes=65536) do acc
            limit = acc.record_limit
            table = DataFrame(rows)
            accumulate_run_summary!(acc, view(table, 1:limit, :))
            @test acc.directory === nothing
            @test length(acc.records) == limit
            @test isempty(readdir(dir))
            accumulate_run_summary!(acc, view(table, limit+1:limit+1, :))
            @test acc.directory !== nothing
            @test isempty(acc.records)
            accumulate_run_summary!(acc, view(table, limit+2:n, :))
        end
        @test all(isequal(getfield(boundary[i], f), getfield(large[i], f))
                  for i in 1:17, f in fieldnames(RunSummaryStats))
        @test isempty(readdir(dir))
        empty_stats = Pioneer.with_run_summary(_ -> nothing, names; temp_parent=dir)
        @test all(s -> s.precursors_identified == 0 && all(ismissing, s.medians), empty_stats)
        @test isempty(readdir(dir))
        # Compare disk selection with Julia's median, including nonfinite values.
        for values in (Float32[-Inf, -2, -0.0, 0.0, 2, Inf],
                       Float32[-Inf, Inf], Float32[1, NaN, 2], Float32[1, 3, 9])
            path = joinpath(dir, "records.bin")
            records = [Pioneer.RunSummaryRecord(1, 1, 0xff, ntuple(_ -> v, 8)) for v in values]
            open(io -> write(io, records), path, "w")
            medians, peptides = Pioneer._large_summary_medians(path, 1, 65536)
            @test all(x -> isequal(x, median(values)), medians)
            @test peptides == 1
            rm(path)
        end
    end
end
