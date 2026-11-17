# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# DIA MS1 quantification (maxLFQ.ms1_quantification): parameter parsing and validation.

@testset "DIA MS1 quantification flag" begin
    defaults = JSON.parsefile(joinpath(@__DIR__, "..", "..", "assets", "example_config", "defaultSearchParams.json"))
    mktempdir() do dir
        path = joinpath(dir, "params.json")
        # omitted -> off
        delete!(defaults["maxLFQ"], "ms1_quantification")
        write(path, JSON.json(defaults))
        @test !Pioneer.ms1_quantification_enabled(Pioneer.parse_pioneer_parameters(path))
        for flag in (true, false)
            defaults["maxLFQ"]["ms1_quantification"] = flag
            write(path, JSON.json(defaults))
            @test Pioneer.ms1_quantification_enabled(Pioneer.parse_pioneer_parameters(path)) == flag
        end
        defaults["maxLFQ"]["ms1_quantification"] = "yes"
        write(path, JSON.json(defaults))
        @test_throws Pioneer.InvalidParametersError Pioneer.checkParams(path)
    end
end

@testset "DIA MS1: normalization keeps missing MS1 areas" begin
    mktempdir() do dir
        path = joinpath(dir, "psms.arrow")
        Arrow.write(path, DataFrame(irt_obs = Float32[1, 2, 3],
                                    peak_area = Union{Missing, Float32}[10, 20, 40],
                                    peak_area_ms1 = Union{Missing, Float32}[8, missing, 32]))
        corr = Pioneer.Dictionary{String, Any}()
        Pioneer.insert!(corr, path, x -> 1.0)          # log2 offset 1: halves every area
        Pioneer.applyNormalization!([path], :peak_area_ms1, corr)
        Pioneer.applyNormalization!([path], :peak_area, corr)
        d = DataFrame(Arrow.Table(path))
        @test ismissing(d.peak_area_ms1_normalized[2])
        @test d.peak_area_ms1_normalized[[1, 3]] ≈ [4.0, 16.0]
        @test d.peak_area_normalized ≈ [5.0, 10.0, 20.0]
        @test eltype(d.peak_area_normalized) == Float64           # MS2 column type unchanged
        # no correction for the file: copied, missing kept
        Pioneer.applyNormalization!([path], :peak_area_ms1, Pioneer.Dictionary{String, Any}())
        @test isequal(DataFrame(Arrow.Table(path)).peak_area_ms1_normalized, [8.0, missing, 32.0])
    end
end

@testset "DIA MS1: protein abundance_ms1 next to abundance" begin
    mktempdir() do dir
        n = 5
        t = DataFrame(inferred_protein_group = fill("A", 3n), precursor_idx = repeat(UInt32.(1:n); inner = 3),
                      ms_file_idx = repeat(UInt32[1, 2, 3]; outer = n),
                      peak_area = repeat(Union{Missing, Float32}[10, 20, 40]; outer = n),
                      peak_area_ms1 = repeat(Union{Missing, Float32}[1, 4, missing]; outer = n),
                      use_for_protein_quant = fill(true, 3n), target = fill(true, 3n),
                      entrapment_group_id = zeros(UInt8, 3n), global_pg_qval = fill(0.001f0, 3n),
                      pg_qval = fill(0.001f0, 3n), pg_pep = fill(0.002f0, 3n), pg_score = fill(4f0, 3n),
                      global_pg_score = fill(5f0, 3n))
        input = joinpath(dir, "in.arrow"); Arrow.write(input, t)
        lfq(out; kw...) = (Pioneer.LFQ_chunked([Pioneer.PSMFileReference(input)], out, :peak_area, ["r1", "r2", "r3"],
                                               ["PEPTIDE$i" for i in 1:n], fill(missing, n), fill(missing, n), 0.01f0,
                                               Dict("A" => "human"); kw...); DataFrame(Arrow.Table(out)))
        base = lfq(joinpath(dir, "base.arrow"))
        ms1 = lfq(joinpath(dir, "ms1.arrow"); quant_col_ms1 = :peak_area_ms1)
        @test !hasproperty(base, :abundance_ms1)
        @test ms1.abundance == base.abundance && ms1.file_name == base.file_name   # MS2 output unchanged
        @test ms1.abundance_ms1[1:2] ≈ [5, 20] rtol = 2e-6                          # 5 precursors x (1, 4)
        @test ismissing(ms1.abundance_ms1[3])                                       # no MS1 in run 3
        @test ms1.n_precursors_quantified_ms1 == UInt32[5, 5, 0]
        @test isequal(ms1.total_peak_area_ms1, Union{Missing, Float32}[5, 20, missing])
    end
end

@testset "MS1 coverage: filter-string boxes, segmented ranges, M0/M1 candidacy" begin
    lo = Float32[]; hi = Float32[]
    @test Pioneer.parse_filter_ranges!(lo, hi, "FTMS + p NSI Full msx ms [500.0000-510.0000, 530.0000-540.0000]") == 2
    @test lo == Float32[500.0, 530.0] && hi == Float32[510.0, 540.0]
    @test Pioneer.parse_filter_ranges!(empty!(lo), empty!(hi), "FTMS + p NSI Full ms [300.0000-1650.0000]") == 1
    @test Pioneer.parse_filter_ranges!(empty!(lo), empty!(hi), "") == 0

    mktempdir() do d
        # cycle: full MS1 300-1650, BoxCar MS1 (two boxes), segmented MS1 650-1400 (header without boxes), one MS2
        n = 4
        df = DataFrame(mz_array = [Union{Missing,Float32}[500f0] for _ in 1:n],
                       intensity_array = [Union{Missing,Float32}[1f0] for _ in 1:n],
                       scanHeader = ["FTMS + p NSI Full ms [300.0000-1650.0000]",
                                     "FTMS + p NSI Full msx ms [500.0000-510.0000, 530.0000-540.0000]",
                                     "FTMS + p NSI Full ms [650.0000-1400.0000]",
                                     "FTMS + p NSI d Full ms2 473.2155@hcd27.00 [65.6667-985.0000]"],
                       scanNumber = Int32.(1:n), packetType = zeros(Int32, n),
                       retentionTime = Float32[1.0, 1.01, 1.02, 1.03],
                       lowMz = Float32[300, 500, 650, 65.67], highMz = Float32[1650, 540, 1400, 985],
                       TIC = ones(Float32, n), msOrder = UInt8[1, 1, 1, 2],
                       centerMz = Vector{Union{Missing,Float32}}([missing, missing, missing, 473.2f0]),
                       isolationWidthMz = Vector{Union{Missing,Float32}}([missing, missing, missing, 1.4f0]),
                       collisionEnergyField = Vector{Union{Missing,Float32}}(fill(27f0, n)),
                       collisionEnergyEvField = zeros(Float32, n))
        Arrow.write(joinpath(d, "cov.arrow"), df)
        cov = Pioneer.Ms1Coverage(Pioneer.BasicMassSpecData(joinpath(d, "cov.arrow")))
        @test !cov.uniform
        @test Pioneer.covers(cov, 1, 1600f0) && Pioneer.covers(cov, 2, 505f0) && !Pioneer.covers(cov, 2, 520f0)
        @test Pioneer.covers(cov, 3, 700f0) && !Pioneer.covers(cov, 3, 600f0)
        @test !Pioneer.covers(cov, 4, 500f0)                                     # MS2 scans measure no MS1 range
        mzs = Float32[505.0, 509.5, 520.0, 700.0]; charges = UInt8[2, 1, 2, 2]
        # BoxCar scan: 505 (M0 in box 1), 509.5 z=1 (M0 in box 1, M1 510.5 outside) kept; 520 (both outside) dropped
        precs = UInt32[1, 2, 3, 4]
        k = Pioneer.keep_ms1_covered!(precs, 4, mzs, charges, cov, 2)
        @test precs[1:k] == UInt32[1, 2]
        # 529.6 z=2: M0 outside box 2 (530-540), M1 530.10 inside -> kept
        precs = UInt32[1]
        @test Pioneer.keep_ms1_covered!(precs, 1, Float32[529.6], UInt8[2], cov, 2) == 1
        # a file whose MS1 scans all measure the same range is uniform
        df2 = df[[1, 4], :]
        Arrow.write(joinpath(d, "plain.arrow"), df2)
        @test Pioneer.Ms1Coverage(Pioneer.BasicMassSpecData(joinpath(d, "plain.arrow"))).uniform
    end
end
