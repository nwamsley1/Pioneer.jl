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
    stats = [RunSummaryStats("a"), RunSummaryStats("b")]
    accumulate_run_summary!(stats, chunk)
    # A second chunk must fold into the same accumulators.
    accumulate_run_summary!(stats, (
        ms_file_idx = UInt32[1], target = Bool[true], sequence = ["NEWPEPK"],
        peak_area = Union{Missing, Float32}[40], peak_area_normalized = Union{Missing, Float32}[20],
        mbr_recovered = Bool[false], irt_error = Float16[4.0], rt_fwhm = Float16[0.4],
        points_integrated = UInt32[9], charge = UInt8[2], missed_cleavage = UInt8[1]))

    a, b = stats
    @test a.precursors_identified == 4
    @test a.precursors_quantified == 4
    @test a.precursors_mbr == 1
    @test length(a.peptides) == 3            # PEPTIDEK counted once
    @test a.total_peak_area == 100.0
    @test median(a.peak_areas) == 25.0f0
    @test median(a.normalization_factors) == 1.5f0   # ratios 2, 1, 2, 0.5
    @test median(a.irt_errors) == 2.5f0             # |−1|,2,3,4
    @test median(a.charges) == 2.0f0
    @test median(a.missed_cleavages) == 0.5f0
    @test median(a.peptide_lengths) == 8.0f0        # 8,8,9,7

    @test b.precursors_identified == 3               # decoy skipped
    @test b.precursors_quantified == 1               # zero and missing areas are unquantified
    @test b.precursors_mbr == 1
    @test length(b.peptides) == 2
    @test b.total_peak_area == 5.0
    @test isempty(b.normalization_factors) || b.normalization_factors == Float32[1.0]

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
