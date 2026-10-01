# Tests for the ion-mobility (CCS) Koina model path in BuildSpecLib:
# request preparation, response parsing, the offline synthetic client, the
# CCS -> 1/K0 conversion, and parameter validation of `im_model` and the
# fragment index's partition widths.

using Test
using JSON
using DataFrames
using Pioneer

const _IonMobilityModel      = Pioneer.IonMobilityModel
const _SyntheticKoinaClient  = Pioneer.SyntheticKoinaClient

@testset "IonMobilityModel — prepare / parse" begin
    df = DataFrame(koina_sequence = ["AAAAAKPK", "C[UNIMOD:4]AAAK"],
                   precursor_charge = UInt8[2, 3])
    batches = Pioneer.prepare_koina_batch(_IonMobilityModel("im2deep"), df; batch_size = 1000)
    @test length(batches) == 1
    req = JSON.parse(batches[1])
    @test Set(i["name"] for i in req["inputs"]) == Set(["peptide_sequences", "precursor_charges"])

    resp = Dict{String,Any}("outputs" => Any[Dict{String,Any}(
        "name" => "ccs", "shape" => Any[2, 1], "datatype" => "FP32", "data" => Any[316.2, 513.2])])
    r = Pioneer.parse_koina_batch(_IonMobilityModel("alphapept_ccs"), resp)
    @test r.fragments.ccs == Float32[316.2, 513.2]
    @test_throws ArgumentError Pioneer.parse_koina_batch(
        _IonMobilityModel("alphapept_ccs"), Dict{String,Any}("outputs" => Any[]))
end

@testset "ccs_to_inv_ion_mobility" begin
    # AAAAAKPK 2+ (m/z 393.24) with AlphaPept CCS 316.2 Å² is ~0.776 Vs/cm²
    # by AlphaPeptDeep's Mason-Schamp conversion.
    im = Pioneer.ccs_to_inv_ion_mobility(316.2115f0, 2, 393.2401)
    @test im isa Float32
    @test isapprox(im, 0.7759; atol = 1e-3)
    # Same CCS at higher charge gives lower 1/K0.
    @test Pioneer.ccs_to_inv_ion_mobility(400.0, 3, 500.0) <
          Pioneer.ccs_to_inv_ion_mobility(400.0, 2, 750.0)
end

@testset "SyntheticKoinaClient — CCS endpoints + predict_ion_mobility" begin
    client = _SyntheticKoinaClient()
    df = DataFrame(koina_sequence = ["AAAAAKPK", "LGEHNIDVLEGNEQFINAAK"],
                   precursor_charge = UInt8[2, 2], mz = Float32[393.24, 1085.55])
    req = Pioneer.prepare_koina_batch(_IonMobilityModel("alphapept_ccs"), df)[1]
    for key in ("alphapept_ccs", "im2deep")
        out = Pioneer.koina_request(client, req, Pioneer.KOINA_URLS[key])
        @test out["outputs"][1]["name"] == "ccs"
        @test length(out["outputs"][1]["data"]) == 2
    end

    mktempdir() do d
        inp = joinpath(d, "in.arrow"); outp = joinpath(d, "out.arrow")
        Pioneer.Arrow.write(inp, df)
        Pioneer.with_koina_client(client) do
            Pioneer.predict_ion_mobility(inp, outp, "alphapept_ccs")
        end
        t = DataFrame(Pioneer.Arrow.Table(outp))
        @test hasproperty(t, :ccs) && hasproperty(t, :inv_ion_mobility)
        @test eltype(t.inv_ion_mobility) == Float32
        @test all(0.5 .< t.inv_ion_mobility .< 2.0)
    end
end

@testset "check_params_bsp — im_model / fragment index widths" begin
    defaults_path = Pioneer.asset_path("example_config", "defaultBuildLibParams.json")
    base = JSON.parsefile(defaults_path)
    base["fasta_paths"] = ["dummy.fasta"]; base["fasta_names"] = ["DUMMY"]
    base["library_path"] = "dummy_lib"; base["calibration_raw_file"] = ""

    ok = deepcopy(base)
    ok["library_params"]["im_model"] = "alphapept_ccs"
    ok["library_params"]["prec_partition_width"] = 10.0
    p = Pioneer.check_params_bsp(JSON.json(ok))
    @test p["library_params"]["im_model"] == "alphapept_ccs"
    @test p["library_params"]["prec_partition_width"] == 10.0

    # Defaults: empty im_model (skip); a 5 and a 10 Da fragment index (SearchDIA picks one from the data), for
    # timsTOF libraries too; the hidden prec_partition_width override builds that one width. Local ID type "auto".
    p0 = Pioneer.check_params_bsp(JSON.json(base))
    @test p0["library_params"]["im_model"] == ""
    @test !haskey(p0["library_params"], "prec_partition_width")
    @test !haskey(p0["library_params"], "isolation_window_width")
    @test Pioneer.fragment_index_widths(p0["library_params"]) == (5.0f0, 10.0f0)
    @test Pioneer.fragment_index_widths(Dict{String, Any}()) == (5.0f0, 10.0f0)
    @test p0["library_params"]["frag_index_local_id_type"] == "auto"
    @test Pioneer.frag_index_local_id_request(p0["library_params"]) == "auto"
    @test Pioneer.frag_index_local_id_request(Dict{String, Any}()) == "auto"
    tims = deepcopy(base); tims["library_params"]["im_model"] = "alphapept_ccs"
    @test Pioneer.fragment_index_widths(Pioneer.check_params_bsp(JSON.json(tims))["library_params"]) == (5.0f0, 10.0f0)
    @test Pioneer.fragment_index_widths(p["library_params"]) == (10.0f0,)        # explicit prec_partition_width = 10
    @test Pioneer.fragment_index_files(5.0f0, true) == ("partitioned_fragment_index.jls", "presearch_partitioned_fragment_index.jls")
    @test Pioneer.fragment_index_files(10.0f0, false) == ("partitioned_fragment_index_w10.jls", "presearch_partitioned_fragment_index_w10.jls")
    @test Pioneer.fragment_index_files(2.5f0, false)[1] == "partitioned_fragment_index_w2.5.jls"
    for v in ("auto", "UInt16", "UInt32")
        c = deepcopy(base); c["library_params"]["frag_index_local_id_type"] = v
        @test Pioneer.check_params_bsp(JSON.json(c))["library_params"]["frag_index_local_id_type"] == v
    end
    bad3 = deepcopy(base); bad3["library_params"]["frag_index_local_id_type"] = "UInt64"
    @test_throws Exception Pioneer.check_params_bsp(JSON.json(bad3))

    bad = deepcopy(base); bad["library_params"]["im_model"] = "nope"
    @test_throws Exception Pioneer.check_params_bsp(JSON.json(bad))
    bad2 = deepcopy(base); bad2["library_params"]["prec_partition_width"] = -1
    @test_throws Exception Pioneer.check_params_bsp(JSON.json(bad2))
end
