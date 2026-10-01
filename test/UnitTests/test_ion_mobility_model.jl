# Tests for the ion-mobility (CCS) Koina model path in BuildSpecLib:
# request preparation, response parsing, the offline synthetic client, the
# CCS -> 1/K0 conversion, and parameter validation of `im_model` /
# `prec_partition_width`.

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

@testset "im_koina_sequence — modifications the CCS model cannot encode are left off" begin
    ap = Pioneer.IM_MODEL_CONFIGS["alphapept_ccs"]; im2 = Pioneer.IM_MODEL_CONFIGS["im2deep"]
    dropped = Dict{String, Int}()
    @test Pioneer.im_koina_sequence("PEPTIDEK", missing, ap, dropped) == "PEPTIDEK"
    # AlphaPeptDeep: N-terminal acetyl as a ProForma prefix, oxidation and carbamidomethyl on their residues
    mods = "(1,n,Unimod:1)(2,M,Unimod:35)(4,C,Unimod:4)"
    @test Pioneer.im_koina_sequence("AMGCK", mods, ap, dropped) == "[UNIMOD:1]-AM[UNIMOD:35]GC[UNIMOD:4]K"
    # IM2Deep drops terminal groups, so the N-terminal acetyl rides on the first residue
    @test Pioneer.im_koina_sequence("AMGCK", mods, im2, dropped) == "A[UNIMOD:1]M[UNIMOD:35]GC[UNIMOD:4]K"
    @test isempty(dropped)
    # A peptide with some encodable and some not: only the unencodable ones go, their residues plain
    tmt = "(1,n,Unimod:737)(3,M,Unimod:35)(5,K,Unimod:737)"
    @test Pioneer.im_koina_sequence("PEMTK", tmt, im2, dropped) == "PEM[UNIMOD:35]TK"
    @test dropped == Dict("Unimod:737 on N-term" => 1, "Unimod:737 on K" => 1)
    @test Pioneer.im_koina_sequence("PEMTK", tmt, ap, dropped) == "[UNIMOD:737]-PEM[UNIMOD:35]TK[UNIMOD:737]"
    @test Pioneer.im_koina_sequence("ACK", "(2,C,Unimod:2062)", ap, dropped) == "ACK"     # DBIA: not in AlphaPeptDeep's table
    @test Pioneer.im_koina_sequence("ACK", "(2,C,Unimod:2062)", im2, dropped) == "AC[UNIMOD:2062]K"
    @test dropped["Unimod:2062 on C"] == 1
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

        # with sequence + mods the request is rebuilt per model; an unencodable modification only warns
        dfm = DataFrame(sequence = ["AAAAAKPK", "PEMTK"], mods = Union{Missing, String}[missing, "(5,K,Unimod:737)"],
                        koina_sequence = ["AAAAAKPK", "PEMTK[UNIMOD:737]"],
                        precursor_charge = UInt8[2, 2], mz = Float32[393.24, 400.0])
        Pioneer.Arrow.write(inp, dfm)
        Pioneer.with_koina_client(client) do
            Pioneer.predict_ion_mobility(inp, outp, "im2deep")
        end
        @test nrow(DataFrame(Pioneer.Arrow.Table(outp))) == 2
    end
end

@testset "check_params_bsp — im_model / isolation_window_width / prec_partition_width" begin
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

    # Defaults: empty im_model (skip); the template's isolation_window_width of 5 gives 5 Da partitions, and so does a
    # config without the key; the width is the isolation window snapped to 2.5, 5 or 10 Da, for timsTOF libraries too;
    # an explicit prec_partition_width always wins. Local ID type defaults to "auto".
    p0 = Pioneer.check_params_bsp(JSON.json(base))
    @test p0["library_params"]["im_model"] == ""
    @test p0["library_params"]["isolation_window_width"] == 5.0
    @test !haskey(p0["library_params"], "prec_partition_width")
    @test Pioneer.prec_partition_width(p0["library_params"]) == 5.0f0
    @test Pioneer.prec_partition_width(Dict{String, Any}()) == 5.0f0
    for (w, expect) in ((25.0, 10.0f0), (14.7, 10.0f0), (7.5, 10.0f0), (7.4, 5.0f0), (4.4, 5.0f0), (3.75, 5.0f0),
                        (3.7, 2.5f0), (2.9, 2.5f0), (2.0, 2.5f0), (1, 2.5f0))
        c = deepcopy(base); c["library_params"]["isolation_window_width"] = w
        @test Pioneer.prec_partition_width(Pioneer.check_params_bsp(JSON.json(c))["library_params"]) == expect
    end
    bad4 = deepcopy(base); bad4["library_params"]["isolation_window_width"] = 0
    @test_throws Exception Pioneer.check_params_bsp(JSON.json(bad4))
    @test p0["library_params"]["frag_index_local_id_type"] == "auto"
    @test Pioneer.frag_index_local_id_request(p0["library_params"]) == "auto"
    @test Pioneer.frag_index_local_id_request(Dict{String, Any}()) == "auto"
    tims = deepcopy(base); tims["library_params"]["im_model"] = "alphapept_ccs"
    @test Pioneer.prec_partition_width(Pioneer.check_params_bsp(JSON.json(tims))["library_params"]) == 5.0f0
    tims["library_params"]["isolation_window_width"] = 25.0     # diaPASEF windows
    @test Pioneer.prec_partition_width(Pioneer.check_params_bsp(JSON.json(tims))["library_params"]) == 10.0f0
    @test Pioneer.prec_partition_width(p["library_params"]) == 10.0f0            # explicit prec_partition_width = 10
    tims5 = deepcopy(tims); tims5["library_params"]["prec_partition_width"] = 5.0  # explicit width beats the window
    @test Pioneer.prec_partition_width(Pioneer.check_params_bsp(JSON.json(tims5))["library_params"]) == 5.0f0
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
