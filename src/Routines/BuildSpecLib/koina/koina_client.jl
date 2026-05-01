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

"""
    AbstractKoinaClient

Pluggable backend for Koina inference requests. Concrete clients implement
`koina_request(client, json::String, model_url::String)::Dict{String,Any}`.

The default client (`HttpKoinaClient`) makes real HTTP calls to the Koina
service. Test code can swap in `SyntheticKoinaClient` to get deterministic,
offline responses for the Altimeter and Chronologer endpoints — letting
BuildSpecLib unit tests run without network access.
"""
abstract type AbstractKoinaClient end

"""
    HttpKoinaClient

Production client. Forwards each request to `make_koina_http_request` (the
existing HTTP-with-retries implementation defined in koina_api.jl).
"""
struct HttpKoinaClient <: AbstractKoinaClient
    max_attempts::Int
    retry_delay::Float64
end
HttpKoinaClient(; max_attempts::Int = 100, retry_delay::Float64 = 1.0) =
    HttpKoinaClient(max_attempts, retry_delay)

"""
    SyntheticKoinaClient

Offline client returning deterministic synthetic responses for Altimeter
(SplineCoefficientModel) fragment prediction and Chronologer
(RetentionTimeModel) iRT prediction. The output shapes and dtypes match
real Koina responses so downstream parse functions exercise the same
code paths.

Synthesis is deterministic in the request payload (peptide_sequences /
precursor_charges) so test fixtures are reproducible across runs.
"""
struct SyntheticKoinaClient <: AbstractKoinaClient
    n_coef_per_frag::Int   # spline degree+1 for Altimeter (typically 4)
    n_frags_per_prec::Int  # number of fragments per precursor (typically 174 for Altimeter)
    n_knots::Int           # length of the global knot vector
end
SyntheticKoinaClient(; n_coef_per_frag::Int = 4, n_frags_per_prec::Int = 174,
                     n_knots::Int = 5) =
    SyntheticKoinaClient(n_coef_per_frag, n_frags_per_prec, n_knots)

"""
    koina_request(client, json::String, model_url::String) -> Dict{String,Any}

Dispatch entry point. The currently active client is selected by
`current_koina_client()`; tests use `with_koina_client(f, client)` to
temporarily install a different one.
"""
function koina_request end

# --- Active client: scoped via task-local storage so concurrent BuildSpecLib
# runs in different tasks don't trample on each other. ---
const _KOINA_CLIENT_TLS_KEY = :pioneer_koina_client

"""
    current_koina_client() -> AbstractKoinaClient

Returns the active Koina client for the current task. Defaults to a fresh
`HttpKoinaClient()` when no override has been installed.
"""
current_koina_client() = get(task_local_storage(), _KOINA_CLIENT_TLS_KEY,
                             HttpKoinaClient())

"""
    with_koina_client(f, client::AbstractKoinaClient)

Run `f()` with `client` as the active Koina client for the current task.
Cleans up afterwards.

Example (offline test):
    with_koina_client(SyntheticKoinaClient()) do
        BuildSpecLib("path/to/build_params.json")
    end
"""
function with_koina_client(f, client::AbstractKoinaClient)
    prev = get(task_local_storage(), _KOINA_CLIENT_TLS_KEY, nothing)
    task_local_storage(_KOINA_CLIENT_TLS_KEY, client)
    try
        return f()
    finally
        if prev === nothing
            delete!(task_local_storage(), _KOINA_CLIENT_TLS_KEY)
        else
            task_local_storage(_KOINA_CLIENT_TLS_KEY, prev)
        end
    end
end

# --- HttpKoinaClient implementation: just forward to make_koina_http_request. ---

function koina_request(client::HttpKoinaClient, json::String, model_url::String)
    make_koina_http_request(json, model_url;
                            max_attempts=client.max_attempts,
                            retry_delay=client.retry_delay)
end

# --- Fixture-based clients (record + replay). -------------------------------
#
# `RecordingKoinaClient` wraps another client and captures every
# (request_payload, model_url, response) tuple it sees, so a one-time live
# Koina run can produce a deterministic on-disk fixture. `FixtureKoinaClient`
# replays that fixture, refusing requests that weren't recorded — used in
# CI / GitHub Actions where Koina is unreachable. Both clients key fixtures
# by a stable hash of the canonicalized request body + model URL.

"""
    _koina_fixture_key(json::String, model_url::String) -> String

Stable key for a Koina request used by both the recorder and the replayer.

Parses the JSON, walks the structure sorting every Dict's keys, then
re-serializes — so the resulting string is deterministic regardless of
the request-side Dict iteration order. Without this, a fixture
recorded on one platform / Julia version would miss when replayed on
another (Julia's `Dict{String, Any}` doesn't guarantee insertion
order, and its iteration order can differ across platforms).
"""
function _koina_fixture_key(json::String, model_url::String)
    parsed = JSON.parse(json; dicttype = Dict{String,Any})
    io = IOBuffer()
    _write_canonical_json(io, parsed)
    canonical = String(take!(io))
    return string(hash((canonical, model_url)); base = 16)
end

# Order-independent JSON serializer — sorts dict keys at every level
# of nesting so the byte string we hash is content-determined.
function _write_canonical_json(io::IO, d::AbstractDict)
    write(io, '{')
    first = true
    for k in sort!(collect(keys(d)))
        first || write(io, ',')
        first = false
        JSON.print(io, k)
        write(io, ':')
        _write_canonical_json(io, d[k])
    end
    write(io, '}')
    return nothing
end

function _write_canonical_json(io::IO, v::AbstractVector)
    write(io, '[')
    first = true
    for x in v
        first || write(io, ',')
        first = false
        _write_canonical_json(io, x)
    end
    write(io, ']')
    return nothing
end

_write_canonical_json(io::IO, x) = JSON.print(io, x)

"""
    RecordingKoinaClient(inner = HttpKoinaClient())

Wraps an inner client (defaults to `HttpKoinaClient`) and stores every
(request, model_url, response) triple that flows through it. Use
`save_koina_fixtures(client, path)` after a run to dump the captures to
a JSON file the `FixtureKoinaClient` can consume.

Example:
    rec = RecordingKoinaClient()
    with_koina_client(rec) do
        BuildSpecLib("path/to/build_params.json")
    end
    save_koina_fixtures(rec, "test/.../koina_fixtures.json")
"""
struct RecordingKoinaClient <: AbstractKoinaClient
    inner::AbstractKoinaClient
    captures::Vector{Tuple{String,String,Dict{String,Any}}}
end
RecordingKoinaClient(inner::AbstractKoinaClient = HttpKoinaClient()) =
    RecordingKoinaClient(inner, Tuple{String,String,Dict{String,Any}}[])

function koina_request(client::RecordingKoinaClient, json::String, model_url::String)
    response = koina_request(client.inner, json, model_url)
    push!(client.captures, (json, model_url, response))
    return response
end

"""
    save_koina_fixtures(client::RecordingKoinaClient, path::AbstractString)

Write the recorder's captures to a JSON file keyed by `_koina_fixture_key`.
Last-write-wins on duplicate keys (which shouldn't happen in practice since
the same canonical payload always maps to the same response).
"""
function save_koina_fixtures(client::RecordingKoinaClient, path::AbstractString)
    fixtures = Dict{String,Dict{String,Any}}()
    for (json, model_url, response) in client.captures
        key = _koina_fixture_key(json, model_url)
        fixtures[key] = response
    end
    open(path, "w") do io
        JSON.print(io, fixtures, 2)
    end
    return path
end

"""
    FixtureKoinaClient(fixture_path::AbstractString)

Offline client that replays Koina responses recorded by
`RecordingKoinaClient` + `save_koina_fixtures`. Errors loudly on a cache
miss so test failures surface as obviously-wrong rather than silently
quirky responses. Used in `test_build_spec_lib.jl` so the BuildSpecLib
unit tests work in CI / GitHub Actions without network access.
"""
struct FixtureKoinaClient <: AbstractKoinaClient
    fixtures::Dict{String,Dict{String,Any}}
    fixture_path::String
end
function FixtureKoinaClient(fixture_path::AbstractString)
    isfile(fixture_path) ||
        error("FixtureKoinaClient: fixture file not found at $fixture_path. " *
              "Record one with `Pioneer.RecordingKoinaClient` + " *
              "`Pioneer.save_koina_fixtures`.")
    raw = JSON.parsefile(fixture_path; dicttype = Dict{String,Any})
    fixtures = Dict{String,Dict{String,Any}}()
    for (k, v) in raw
        fixtures[k] = v
    end
    return FixtureKoinaClient(fixtures, String(fixture_path))
end

function koina_request(client::FixtureKoinaClient, json::String, model_url::String)
    key = _koina_fixture_key(json, model_url)
    haskey(client.fixtures, key) ||
        error("FixtureKoinaClient: no fixture for key=$key (model=$model_url, " *
              "payload bytes=$(length(json))). Re-record fixtures via " *
              "`RecordingKoinaClient` against live Koina, or check whether the " *
              "test scenarios changed since the fixture at $(client.fixture_path) " *
              "was captured.")
    return client.fixtures[key]
end

# --- SyntheticKoinaClient implementation. ---

function koina_request(client::SyntheticKoinaClient, json::String, model_url::String)
    request = JSON.parse(json; dicttype=Dict{String,Any})
    inputs = request["inputs"]
    # Find the peptide_sequences input to determine batch size + driver text.
    seq_input = first(i for i in inputs if i["name"] == "peptide_sequences")
    sequences = String.(seq_input["data"])
    n_precs = length(sequences)
    if occursin("Chronologer", model_url) || occursin("chronologer", model_url)
        return _synthetic_chronologer_response(client, sequences)
    elseif occursin("Altimeter", model_url) || occursin("altimeter", model_url)
        return _synthetic_altimeter_response(client, sequences)
    else
        throw(ArgumentError("SyntheticKoinaClient only supports Altimeter and " *
                            "Chronologer endpoints; got model_url=$model_url"))
    end
end

# Hash a string to a stable Float32 in [0, 1). Used for deterministic synthesis.
@inline function _stable_unit(s::AbstractString, salt::Int = 0)
    h = UInt64(hash(s, UInt64(salt)))
    Float32((h % UInt64(10_000_000)) / UInt64(10_000_000))
end

function _synthetic_chronologer_response(::SyntheticKoinaClient, sequences::Vector{String})
    # iRT scaled into a plausible range (~0..200 minutes) deterministically.
    rt_data = Float64[Float64(_stable_unit(s, 1)) * 200.0 for s in sequences]
    return Dict{String,Any}(
        "outputs" => Any[
            Dict{String,Any}(
                "name" => "rt",
                "shape" => Any[length(sequences), 1],
                "datatype" => "FP32",
                "data" => rt_data,
            )
        ]
    )
end

function _synthetic_altimeter_response(client::SyntheticKoinaClient, sequences::Vector{String})
    n_precs = length(sequences)
    F = client.n_frags_per_prec
    K = client.n_coef_per_frag

    # Coefficients shape: (n_precs, K, F) flattened in the same row-major
    # order Altimeter uses: idx(i,k,j) = (i-1)*K*F + (k-1)*F + j  (1-based).
    coefficients = Vector{Float64}(undef, n_precs * K * F)
    mz = Vector{Float64}(undef, n_precs * F)
    annotations = Vector{Int32}(undef, n_precs * F)

    @inbounds for i in 1:n_precs
        seq = sequences[i]
        for j in 1:F
            # Synthetic mz: distinct per (seq, j) but in a plausible range.
            mz[(i-1)*F + j] = 100.0 + 1900.0 * _stable_unit(seq, 100 + j)
            # Synthetic annotation: small positive integer code.
            annotations[(i-1)*F + j] = Int32((j - 1) % 256)
            for k in 1:K
                coefficients[(i-1)*K*F + (k-1)*F + j] =
                    Float64(_stable_unit(seq, 1000 + k * 1000 + j))
            end
        end
    end

    knots = collect(range(0.0, 1.0; length=client.n_knots))

    return Dict{String,Any}(
        "outputs" => Any[
            Dict{String,Any}("name" => "mz",
                             "shape" => Any[n_precs, F],
                             "datatype" => "FP32",
                             "data" => mz),
            Dict{String,Any}("name" => "annotations",
                             "shape" => Any[n_precs, F],
                             "datatype" => "INT32",
                             "data" => annotations),
            Dict{String,Any}("name" => "knots",
                             "shape" => Any[client.n_knots],
                             "datatype" => "FP32",
                             "data" => knots),
            Dict{String,Any}("name" => "coefficients",
                             "shape" => Any[n_precs, K, F],
                             "datatype" => "FP32",
                             "data" => coefficients),
        ]
    )
end
