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
Make a real HTTP request to the Koina API with retries. Production path.
Tests typically don't call this directly — they call `make_koina_request`
which dispatches through the active `AbstractKoinaClient` (see
koina_client.jl) and may resolve to a `SyntheticKoinaClient`.

The `poster` kwarg is a dependency-injection seam for unit tests: the
default is `HTTP.post`, but tests can pass a callable returning a
mock response object (anything with a `.body::Vector{UInt8}` field
holding JSON-encoded bytes) to exercise the retry / error-handling
logic without network access.
"""
function make_koina_http_request(json_data::String,
                          model_url::String;
                          max_attempts::Int = 100,
                          retry_delay::Float64 = 1.0,
                          poster = HTTP.post)
    attempt = 1


    while attempt <= max_attempts
        try
            response = JSON.parse(String(poster(model_url, body=json_data).body), dicttype=Dict{String,Any})
            if !haskey(response, "error")
                return response
            end
            throw(KoinaRequestError("Error in response", attempt, response))
        catch e
            if attempt == max_attempts
                @user_error "Failed after $max_attempts attempts: $e"
                rethrow(e)
            end
            if e isa KoinaRequestError
                @debug_l2 "Request failed (attempt $attempt): $(e.message)"
            else
                @debug_l2 "Request failed (attempt $attempt): $(sprint(showerror, e))"
            end
        end
        sleep(retry_delay)
        attempt += 1
    end
    error("Koina API request failed after $max_attempts attempts. This may indicate network issues, API rate limiting, or service unavailability. Check your internet connection and try again later.")
end

"""
    make_koina_request(json_data, model_url) -> Dict{String,Any}

Public entry point. Routes through the currently active `AbstractKoinaClient`
(`HttpKoinaClient` by default, `SyntheticKoinaClient` when installed via
`with_koina_client(...)` for offline tests).
"""
make_koina_request(json_data::String, model_url::String; kwargs...) =
    koina_request(current_koina_client(), json_data, model_url)



"""
    make_koina_batch_requests(json_vec, model_url;
                              concurrency = 24)

Send every JSON string in `json_vec` to `model_url` while keeping at most
`concurrency` HTTP requests active at the same time.
Returns a vector of responses in the **original order**.
"""
function make_koina_batch_requests(json_vec::Vector{String},
                                   model_url::String;
                                   concurrency::Int = 24)

    n          = length(json_vec)
    results    = Vector{Any}(undef, n)          # output buffer
    job_chan   = Channel{Tuple{Int,String}}(n)  # work queue

    # Capture the active Koina client from the calling task. Each @async
    # worker starts with a fresh task_local_storage, so without this the
    # workers fall back to the default HttpKoinaClient — bypassing any
    # RecordingKoinaClient / FixtureKoinaClient / SyntheticKoinaClient
    # installed by the caller via `with_koina_client`.
    parent_client = current_koina_client()

    # ───────────────── fill the channel (producer) ────────────────────
    @async begin
        for (idx, js) in pairs(json_vec)
            put!(job_chan, (idx, js))
        end
        close(job_chan)                         # signal "no more jobs"
    end

    # ───────────────── worker task definition ─────────────────────────
    function worker()
        with_koina_client(parent_client) do
            for (idx, js) in job_chan           # blocks until work available
                results[idx] = make_koina_request(js, model_url)
            end
        end
        return nothing
    end

    # ───────────────── launch the pool of workers ─────────────────────
    tasks = [@async worker() for _ in 1:concurrency]
    fetch.(tasks)                               # wait for all workers
    return results
end
