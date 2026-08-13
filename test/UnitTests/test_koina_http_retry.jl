# Unit tests for make_koina_http_request retry / error-handling logic.
#
# The function is the production HTTP path that lives behind the
# AbstractKoinaClient indirection. We exercise it offline by passing a
# mock `poster` callable in place of HTTP.post — the only contract is
# that `poster(url; body=...)` returns an object with a `.body` field
# containing JSON-encoded bytes (or throws to simulate a network error).
# All tests use very small `retry_delay` so total runtime stays sub-second.

using Test
using JSON
using Pioneer

const _make_koina_http_request = Pioneer.make_koina_http_request
const _KoinaRequestError = Pioneer.KoinaRequestError

# Minimal fake-response stand-in. The real HTTP.post returns
# HTTP.Response, but the function only reads `.body`, so a NamedTuple
# is sufficient.
_resp(payload::Dict) = (body = Vector{UInt8}(JSON.json(payload)),)

@testset "make_koina_http_request — happy path returns parsed response" begin
    n_calls = Ref(0)
    poster = function (url; body)
        n_calls[] += 1
        return _resp(Dict("outputs" => Any[Dict("name" => "rt", "data" => Float64[1.0, 2.0])]))
    end
    response = _make_koina_http_request("ignored", "http://mock";
                                        max_attempts=3, retry_delay=0.0001,
                                        poster=poster)
    @test n_calls[] == 1
    @test haskey(response, "outputs")
    @test response["outputs"][1]["name"] == "rt"
end

@testset "make_koina_http_request — retries on response with 'error' key" begin
    n_calls = Ref(0)
    poster = function (url; body)
        n_calls[] += 1
        if n_calls[] == 1
            return _resp(Dict("error" => "transient backend hiccup"))
        else
            return _resp(Dict("outputs" => Any[]))
        end
    end
    response = _make_koina_http_request("ignored", "http://mock";
                                        max_attempts=3, retry_delay=0.0001,
                                        poster=poster)
    @test n_calls[] == 2
    @test !haskey(response, "error")
    @test haskey(response, "outputs")
end

@testset "make_koina_http_request — retries on thrown exception" begin
    n_calls = Ref(0)
    poster = function (url; body)
        n_calls[] += 1
        n_calls[] == 1 && error("simulated network failure")
        return _resp(Dict("outputs" => Any[]))
    end
    response = _make_koina_http_request("ignored", "http://mock";
                                        max_attempts=3, retry_delay=0.0001,
                                        poster=poster)
    @test n_calls[] == 2
    @test haskey(response, "outputs")
end

@testset "make_koina_http_request — rethrows after max_attempts" begin
    n_calls = Ref(0)
    poster = function (url; body)
        n_calls[] += 1
        error("always fails")
    end
    @test_throws ErrorException _make_koina_http_request("ignored", "http://mock";
                                                         max_attempts=3, retry_delay=0.0001,
                                                         poster=poster)
    @test n_calls[] == 3   # tried 3 times, then rethrew
end

@testset "make_koina_http_request — KoinaRequestError on persistent error response" begin
    poster = function (url; body)
        return _resp(Dict("error" => "permanent failure"))
    end
    @test_throws _KoinaRequestError _make_koina_http_request("ignored", "http://mock";
                                                             max_attempts=2, retry_delay=0.0001,
                                                             poster=poster)
end
