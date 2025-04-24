"""
Make a request to the Koina API with retries.
"""
function make_koina_request(json_data::String, 
                          model_url::String; 
                          max_attempts::Int = 10,
                          retry_delay::Float64 = 1.0)
    attempt = 1
    

    while attempt <= max_attempts
        try
            response = JSON.parse(String(HTTP.post(model_url, body=json_data).body))
            if !haskey(response, "error")
                return response
            end
            throw(KoinaRequestError("Error in response", attempt, response))
        catch e
            if attempt == max_attempts
                @error "Failed after $max_attempts attempts" exception=e
                rethrow(e)
            end
            if e isa KoinaRequestError
                @warn "Request failed (attempt $attempt): $(e.message)"
            else
                @warn "Request failed (attempt $attempt): $(sprint(showerror, e))"
            end
        end
        sleep(retry_delay)
        attempt += 1
    end
    error("Failed after $max_attempts attempts")
end



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

    # ───────────────── fill the channel (producer) ────────────────────
    @async begin
        for (idx, js) in pairs(json_vec)
            put!(job_chan, (idx, js))
        end
        close(job_chan)                         # signal "no more jobs"
    end

    # ───────────────── worker task definition ─────────────────────────
    function worker()
        for (idx, js) in job_chan               # blocks until work available
            results[idx] = make_koina_request(js, model_url)
        end
        return nothing
    end

    # ───────────────── launch the pool of workers ─────────────────────
    tasks = [@async worker() for _ in 1:concurrency]
    fetch.(tasks)                               # wait for all workers
    return results
end