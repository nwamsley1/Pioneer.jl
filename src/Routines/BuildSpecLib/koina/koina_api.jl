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
Make multiple Koina requests asynchronously.
"""
function make_koina_batch_requests(json_data_list::Vector{String}, model_url::String)
    tasks = [@async make_koina_request(json_data, model_url) for json_data in json_data_list]
    return fetch.(tasks)
end