# Tests for the offline SyntheticKoinaClient + lookup-by-name parsing.
#
# The synthetic client lets BuildSpecLib unit tests run without network
# access by returning deterministic Koina-shaped responses for the
# Altimeter (SplineCoefficientModel) and Chronologer (RetentionTimeModel)
# endpoints. These tests pin the contract between the client, parse_koina_batch,
# and downstream consumers.

using Test
using JSON
using DataFrames
using Pioneer

const _SyntheticKoinaClient = Pioneer.SyntheticKoinaClient
const _HttpKoinaClient      = Pioneer.HttpKoinaClient
const _with_koina_client    = Pioneer.with_koina_client
const _current_koina_client = Pioneer.current_koina_client
const _koina_request        = Pioneer.koina_request
const _make_koina_request   = Pioneer.make_koina_request
const _SplineCoefficientModel = Pioneer.SplineCoefficientModel
const _RetentionTimeModel   = Pioneer.RetentionTimeModel
const _parse_koina_batch    = Pioneer.parse_koina_batch

@testset "SyntheticKoinaClient — Chronologer" begin
    client = _SyntheticKoinaClient()

    request = JSON.json(Dict(
        "id" => "0",
        "inputs" => Any[
            Dict("name" => "peptide_sequences",
                 "shape" => Any[3, 1],
                 "datatype" => "BYTES",
                 "data" => ["PEPTIDE", "ANOTHER", "TRIPLE"])
        ]
    ))
    response = _koina_request(client, request,
        "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer")

    @test haskey(response, "outputs")
    @test length(response["outputs"]) == 1
    rt_output = response["outputs"][1]
    @test rt_output["name"] == "rt"
    @test length(rt_output["data"]) == 3
    # Determinism: same input → same output
    response2 = _koina_request(client, request,
        "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer")
    @test response["outputs"][1]["data"] == response2["outputs"][1]["data"]
end

@testset "SyntheticKoinaClient — Altimeter" begin
    client = _SyntheticKoinaClient(n_coef_per_frag=4, n_frags_per_prec=10, n_knots=5)

    request = JSON.json(Dict(
        "id" => "0",
        "inputs" => Any[
            Dict("name" => "peptide_sequences",
                 "shape" => Any[2, 1],
                 "datatype" => "BYTES",
                 "data" => ["PEPTIDE", "TESTPEP"]),
            Dict("name" => "precursor_charges",
                 "shape" => Any[2, 1],
                 "datatype" => "INT32",
                 "data" => Int32[2, 3]),
        ]
    ))
    response = _koina_request(client, request,
        "https://koina.wilhelmlab.org:443/v2/models/Altimeter_2024/infer")

    @test haskey(response, "outputs")
    output_names = Set(o["name"] for o in response["outputs"])
    @test output_names == Set(["mz", "annotations", "knots", "coefficients"])

    # Look up coefficients by name (not by index). Validate shape.
    coefs = first(o for o in response["outputs"] if o["name"] == "coefficients")
    @test coefs["shape"] == [2, 4, 10]   # n_precs, K, F
    @test length(coefs["data"]) == 2 * 4 * 10

    knots = first(o for o in response["outputs"] if o["name"] == "knots")
    @test length(knots["data"]) == 5

    # Parse via parse_koina_batch and verify the DataFrame is well-formed.
    model = _SplineCoefficientModel("altimeter")
    batch_result = _parse_koina_batch(model, response)
    @test batch_result.frags_per_precursor == 10
    @test length(batch_result.extra_data) == 5
    @test nrow(batch_result.fragments) == 2 * 10  # n_precs * n_frags_per_prec
    @test :coefficients in propertynames(batch_result.fragments)
    @test :mz in propertynames(batch_result.fragments)
    @test :annotation in propertynames(batch_result.fragments)
    # Coefficient tuple width matches n_coef_per_frag.
    @test length(batch_result.fragments.coefficients[1]) == 4
end

@testset "SyntheticKoinaClient — rejects unknown endpoints" begin
    client = _SyntheticKoinaClient()
    request = JSON.json(Dict(
        "id" => "0",
        "inputs" => Any[
            Dict("name" => "peptide_sequences",
                 "shape" => Any[1, 1],
                 "datatype" => "BYTES",
                 "data" => ["PEPTIDE"])
        ]
    ))
    @test_throws ArgumentError _koina_request(client, request,
        "https://koina.wilhelmlab.org:443/v2/models/UnknownModel/infer")
end

@testset "with_koina_client scoping" begin
    # Default = HttpKoinaClient
    @test _current_koina_client() isa _HttpKoinaClient

    synth = _SyntheticKoinaClient()
    _with_koina_client(synth) do
        @test _current_koina_client() === synth
        # make_koina_request goes through the active client, so within the
        # scope it returns synthetic results without hitting the network.
        request = JSON.json(Dict(
            "id" => "0",
            "inputs" => Any[
                Dict("name" => "peptide_sequences",
                     "shape" => Any[1, 1],
                     "datatype" => "BYTES",
                     "data" => ["PEPTIDE"])
            ]
        ))
        response = _make_koina_request(request,
            "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer")
        @test response["outputs"][1]["name"] == "rt"
    end

    # After scope: back to HttpKoinaClient (default)
    @test _current_koina_client() isa _HttpKoinaClient
end

@testset "Spline parse looks up 'coefficients' by name (not by index 4)" begin
    # Construct a response where coefficients is the FIRST output.
    response = Dict{String,Any}(
        "outputs" => Any[
            Dict{String,Any}("name" => "coefficients",
                             "shape" => Any[1, 2, 3],
                             "datatype" => "FP32",
                             "data" => Float64[
                                 # idx(i,k,j) = (i-1)*K*F + (k-1)*F + j with K=2, F=3
                                 1.0, 2.0, 3.0,   # k=1: j=1,2,3
                                 4.0, 5.0, 6.0,   # k=2: j=1,2,3
                             ]),
            Dict{String,Any}("name" => "mz",
                             "shape" => Any[1, 3],
                             "datatype" => "FP32",
                             "data" => Float64[100.0, 200.0, 300.0]),
            Dict{String,Any}("name" => "annotations",
                             "shape" => Any[1, 3],
                             "datatype" => "INT32",
                             "data" => Int32[10, 20, 30]),
            Dict{String,Any}("name" => "knots",
                             "shape" => Any[3],
                             "datatype" => "FP32",
                             "data" => Float64[0.0, 0.5, 1.0]),
        ]
    )
    model = _SplineCoefficientModel("altimeter")
    batch_result = _parse_koina_batch(model, response)
    @test batch_result.frags_per_precursor == 3
    @test length(batch_result.extra_data) == 3
    # Coefficient tuple per fragment, K=2 wide
    @test batch_result.fragments.coefficients[1] == (1.0f0, 4.0f0)
    @test batch_result.fragments.coefficients[2] == (2.0f0, 5.0f0)
    @test batch_result.fragments.coefficients[3] == (3.0f0, 6.0f0)
end
