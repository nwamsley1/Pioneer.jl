
function predictRetentionTimes(
    peptide_table_path::String,
    model_type::KoinaModelType,
    max_koina_batches::Int,
    batch_size::Int,
    model_name::String
)
    if model_name ∉ ["chronologer"]
        error("For retention time prediction, pioneer only suppoer `chronologer`. Model name given was: `$model_name`")
    end
    models_dict = Dict(
        "chronologer" => "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer"
    )
    model_url = nothing
    if haskey(models_dict, model_name)
        model_url = models_dict[model_name]
    else
        throw(ErrorException("The retention time model name '$model_name' isn't valid. Valid models are "*join(keys(models_dict),';')))
    end

    function prepareBatch(
        modified_sequence::AbstractVector{String},
        model_type::RetentionTimeModel)
        nrows = length(modified_sequence)
        data_batch = Dict(
            "id" => "0",
            "inputs" => [
                Dict("name" => "peptide_sequences", "shape" => [nrows,1], "datatype" => "BYTES", "data" => strip.(modified_sequence))
            ]
        )

        return JSON.json(data_batch)
    end

    function getKoinaBatches(
        peptides_df::DataFrame,
        model_type::KoinaModelType,
        model_url::String,
        batch_size::Int)

        batch_size = min(batch_size, 1000)
        nprecs = size(peptides_df, 1)
        batch_start_idxs = collect(one(UInt32):UInt32(batch_size):UInt32(nprecs))
        batch_json_inputs = [prepareBatch(peptides_df[start_idx:min((start_idx + batch_size - 1), nprecs),:chronologer_sequence], 
        model_type) for start_idx in batch_start_idxs]
        #return batch_json_inputs[1]
        out_requests = multi_request(batch_json_inputs, model_url)
        #return out_requestsd
        return vcat([first(parseBatchToTable(request["outputs"], model_type)) for request in out_requests]...)[!,:rt];
    end
        
    #test on first 5K for now
    peptides_df = DataFrame(Arrow.Table(peptide_table_path))#[1:5000,:]
    retention_times = zeros(Float32, size(peptides_df, 1))
    nprecs = size(peptides_df, 1)
    batch_size = min(batch_size, 1000)
    batch_start_idxs = collect(one(UInt32):UInt32(batch_size*max_koina_batches):UInt32(nprecs))
    for start_idx in ProgressBar(batch_start_idxs)
        stop_idx = min(start_idx + batch_size*max_koina_batches - 1, size(peptides_df, 1))
        retention_times[start_idx:stop_idx] = getKoinaBatches(peptides_df[start_idx:stop_idx,:],
                        model_type,
                        model_url,
                        batch_size)
    end
    peptides_df[!,:rt] = retention_times
    precursors_path = joinpath(dirname(peptide_table_path), "precursors.arrow")
    Arrow.write(precursors_path, peptides_df)
end



#=
 cmd = `curl "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer" \
    -d '{
        "id": "0", "inputs": [
            {"name": "peptide_sequences", "shape": [1,1], "datatype": "BYTES", "data": [\"AEVTPSQHGNR\"]}]
        }'`


        batch_json_inputs[1]

model_url = "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer"
json_data = "{\"id\":\"0\",\"inputs\":[{\"name\":\"peptide_sequences\",\"shape\":[5,1],\"data\":[\"AEVTPSQHGNR\",\"AEVTPSQHGNR\",\"AEVTPSQHGNRTFSYTLEDHTK\",\"AEVTPSQHGNRTFSYTLEDHTK\",\"AM[UNIMOD:35]FTNGLR\"],\"datatype\":\"BYTES\"}]}";
cmd = `curl -s $model_url -d $json_data`;
first(JSON.parse(read(cmd, String))["outputs"])
json_data = read(cmd, String)
cmd = `curl -s $model_url -d $json_data`;
read(cmd, String)

parseBatchToTable(JSON.parse(read(cmd, String))["outputs"], 
RetentionTimeModel(""))


model_url = "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer"

json_data =  predictRetentionTimes(
            chronologer_out_path,
            RetentionTimeModel("chronologer"),
            24,
            1000,
            "chronologer"
        )
cmd = `curl -s $model_url -d $json_data`;
@time parseBatchToTable(JSON.parse(read(cmd, String))["outputs"], 
RetentionTimeModel(""))


model_url = "https://koina.wilhelmlab.org:443/v2/models/Chronologer_RT/infer"

json_data =  predictRetentionTimes(
            chronologer_out_path,
            RetentionTimeModel("chronologer"),
            24,
            1000,
            "chronologer"
        )
cmd = `curl -s $model_url -d $json_data`;
@time parseBatchToTable(JSON.parse(read(cmd, String))["outputs"], 
RetentionTimeModel(""))

=#