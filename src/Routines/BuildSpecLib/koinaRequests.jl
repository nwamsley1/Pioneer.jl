"""
    make_request(json_data::String,model_url::String;max_attempts::Int64 = 100)

Send a request to the koina server to predict fragment intensities

### Input
-'json_data::String' json formatted input with peptide sequences, charges, collision energies and instrument types
-'model_url::String' url for a koina model. See https://koina.wilhelmlab.org/docs
-'max_attempts::Int64' maximum number of times to attempt the request before abandoning 

### Output
-'::Vector{Dict{String,Any}}'

### Examples
From the Koina Website

curl "https://koina.wilhelmlab.org:443/v2/models/UniSpec/infer" \
    -d '{
        "id": "0", "inputs": [
            {"name": "peptide_sequences", "shape": [2,1], "datatype": "BYTES", "data": ["YYHTLFTHSLPK", "AGC[UNIMOD:4]FS[UNIMOD:27]PK"]},
            {"name": "precursor_charges", "shape": [2,1], "datatype": "INT32", "data": [1,2]},
            {"name": "collision_energies", "shape": [2,1], "datatype": "FP32", "data": [25, 25]},
            {"name": "instrument_types", "shape": [2,1], "datatype": "BYTES", "data": ["LUMOS", "QE"]}]
        }'

To format and input in Julia, do the following 

julia> data_batch = Dict(
    "id" => "0",
    "inputs" => [
        Dict("name" => "peptide_sequences", "shape" => [nrows,1], "datatype" => "BYTES", "data" => getModifiedSequence.(peptide_table[!,:sequence], "", peptide_table[!,:mods])),
        Dict("name" => "precursor_charges", "shape" => [nrows,1], "datatype" => "INT32", "data" => peptide_table[!,:precursor_charge]),
        Dict("name" => "collision_energies", "shape" => [nrows,1], "datatype" => "FP32", "data" =>  peptide_table[!,:collision_energy]),
        Dict("name" => "instrument_types", "shape" => [nrows,1], "datatype" => "BYTES", "data" => ["QE" for _ in range(1, nrows)])
    ]
)
julia> json_data = JSON.json(data_batch)
julia> model_url = "=>"https://koina.wilhelmlab.org:443/v2/models/UniSpec/infer"
julia> predicted_fragments_json = make_request(json_data, model_url)
"""
function make_request(json_data::String,
                      model_url::String;
                      max_attempts::Int64 = 100
                      )
    attempt_counter = 1
    cmd = `curl -s $model_url -d $json_data`
    #Request and read output into a dictionary
    output = Dict(
        "error" => nothing
    )
    try 
        output = JSON.parse(read(cmd, String))
    catch
        @warn "Koina server failed to process a request. Reatempting... Attempts counter: $attempt_counter"
    end
    while haskey(output, "error")
        attempt_counter += 1
        try
            output = JSON.parse(read(cmd, String))
        catch
            @warn "Koina server failed to process a request. Reatempting... Attempts counter: $attempt_counter"
        end
        if attempt_counter > max_attempts
            throw(ErrorException("Exceeded maximum attempted server requests"))
        end
    end
    return output
end

"""
multi_request(json_data_list::Vector{String},
              model_url::String;
              max_attempts::Int64 = 100)

Asyncrhonously make a series of requests to the koina server

### Input
-'json_data_lst::Vector{String}' list of json formatted inputs with peptide sequences, charges, collision energies and instrument types
-'model_url::String' url for a koina model. See https://koina.wilhelmlab.org/docs
-'max_attempts::Int64' maximum number of times to attempt the request before abandoning 

### Output
-'::Vector{Vector{Dict{String,Any}}}'

"""
function multi_request(json_data_list::Vector{String},
                       model_url::String)
    tasks = [
        @async make_request(json_data, model_url)
        for json_data in json_data_list
    ]
    return fetch.(tasks)
end



function insert_at_indices(original::String, insertions::Vector{Tuple{String, UInt8}})
    # Convert the original string into an array of characters for easier manipulation
    char_array = collect(original)

    # Sort the insertions by index in ascending order
    sorted_insertions = sort(insertions, by = x -> x[2])

    # Adjust the index for each insertion
    offset = 0
    for (substr, idx) in sorted_insertions
        # Adjust the index with the current offset
        insertion_index = idx + offset
        # Insert each character of the substring at the specified index
        for (i, char) in enumerate(substr)
            insert!(char_array, insertion_index + i, char)
        end
        # Update the offset by the length of the inserted substring
        offset += length(substr)
    end

    # Join the array of characters back into a single string
    return join(char_array)
end
function getModifiedSequence(
    sequence::String,
    isotope_mods::String,
    structural_mods::String)

    mods = structural_mods*isotope_mods
    mods = [("["*uppercase(getModName(mod.match))*"]", getModIndex(mod.match)) for mod in parseMods(mods)]
    return insert_at_indices(sequence, mods)
end

function getModifiedSequence(
    sequence::String,
    isotope_mods::Missing,
    structural_mods::String)

    mods = structural_mods
    mods = [("["*uppercase(getModName(mod.match))*"]", getModIndex(mod.match)) for mod in parseMods(mods)]
    return insert_at_indices(sequence, mods)
end

function getModifiedSequence(
    sequence::String,
    isotope_mods::String,
    structural_mods::Missing)

    mods = isotope_mods
    mods = [("["*uppercase(getModName(mod.match))*"]", getModIndex(mod.match)) for mod in parseMods(mods)]
    return insert_at_indices(sequence, mods)
end


"""
parseBatchToTable(json_out::Vector{Any})

Parse text output from koina into a dataframe

### Input
-'json_out::Vector{Any}' list of dictionaries with the predicted fragment intensities, masses, and names 

### Output
-'::DataFrame'
"""
function parseBatchToTable(
    json_out::Vector{Any})::Tuple{DataFrame, Int64}
    df = DataFrame()
    n_precs, n_frags_per_prec = first(json_out)["shape"]
    n_rows = n_precs*n_frags_per_prec
    for col in json_out

        col_name = Symbol(col["name"])
        if col_name ∈[:intensities, :mz]
            df[!,col_name] = Float32.(col["data"])
        else
            #if length(col["data"]) == size(df, 1)
                df[!,:annotation] = string.(col["data"])
            #=
            else
               
                df[!,:annotation] = Vector{String}(undef, first(col["shape"])*last(col["shape"]))
                i = 1
                for peptide in col["data"]
                    for frag in peptide
                        df[i,:annotation] = string.(frag)
                        i += 1
                    end
                end
                #df[!,:annotation] = string.(reduce(vcat, col["data"]))
            end
            =#
        end
    end
    return (df, Int64(last(first(json_out)["shape"])))
end

"""
addPrecursorID!(df::DataFrame, start_prec_idx::UInt32, frags_per_prec::Int64)
adds precursor id's to the koina output formated as a dataframe 

### Input
-'df::DataFrame' has columns :intensities, :mz, and :annotation
-'start_prec_idx::UInt32' precursor_idx of the first precursor
-'frags_per_prec::Int64' number of fragments predicted per precursor. Must be the same for all precursors

### Output
-'::DataFrame'
"""
function addPrecursorID!(df::DataFrame, start_prec_idx::UInt32, frags_per_prec::Int64)
    n_precursors_in_batch = UInt32(fld(size(df, 1), frags_per_prec))
    df[!,:precursor_idx] = repeat(start_prec_idx:(start_prec_idx + n_precursors_in_batch - one(UInt32)), inner=frags_per_prec)
    return nothing
end


function filterEachPrecursor!(
    df::DataFrame;
    intensity_threshold::Float32 = 0.001f0
    )
    ctype = eltype(df[!,:intensities])
    filter!(x->x.intensities::ctype>intensity_threshold, df) #Filter invalid intensities
    filter!(x->x.mz::ctype>zero(Float32), df) #Filter invalid masses
    filter!(x->!occursin('i', x.annotation::String), df) #Filter out M1+ isotopes
end

function sortEachPrecursor!(df::DataFrame, col::Symbol)
    nrows = nrow(df)
    row_idx = 1
    while row_idx < nrows
        start_idx = row_idx
        prec_idx = df[row_idx,:precursor_idx]
        while (prec_idx == df[row_idx,:precursor_idx])
            row_idx += 1
            if row_idx > nrows
                break
            end
        end
        stop_idx = row_idx-1 
        block_view = @view df[start_idx:stop_idx, :]
        sort!(block_view, col, rev = true)
    end
end


"""
predictFragments(
    peptide_table_path::String,
    frags_out_path::String,
    max_koina_batches::Int,
    batch_size::Int,
    model_name::String;
    intensity_threshold::Float32 = 0.001f0)

Asyncrhonously make a series of requests to the koina server

### Input
-'peptide_table_path::String' path to peptide table in .arrow format. At minimum needs columns :sequence, :precursor_charge, and :collision_energy
-'frags_out_path::String' path to write the table of predicted fragment  ions
-'max_attempts::Int64' maximum number of requests to make to koina at once (1-24 are reasonable choices)
-'batch_size::Int64' number of precursors to include in a batch (should not exceed what can easily be stored in memory)
-'model_name::String' name of the model to use. 

### Output
-writes the predicted fragment ions to a .arrow table. 
```
julia> test_frags
Arrow.Table with 10217000 rows, 4 columns, and schema:
 :intensities    Float32
 :mz             Float32
 :annotation     String
 :precursor_idx  UInt32

julia> last(DataFrame(test_frags), 5)
5×4 DataFrame
 Row │ intensities  mz        annotation  precursor_idx 
     │ Float32      Float32   String      UInt32        
─────┼──────────────────────────────────────────────────
   1 │  0.00160101  1587.27   y28-NH3^2          165600
   2 │  0.001587    1651.3    y29-NH3^2          165600
   3 │  0.00140379  1354.67   y24-H2O^2          165600
   4 │  0.00136275   966.155  y26^3              165600
   5 │  0.0012414   1009.84   y27^3              165600

```
### Example
"""
function predictFragments(
    peptide_table_path::String, #"data/for_unispec/keap1_forunispec2.arrow"
    frags_out_path::String, #"/Users/n.t.wamsley/Desktop/test_out.arrow",
    model_type::KoinaModelType,
    instrument_type::String,
    max_koina_batches::Int,
    batch_size::Int,
    model_name::String,
    intensity_threshold::Float32 = 0.001f0)

    if model_name == "unispec'"
        if instrument_type ∉ ["QE","QEHFX","LUMOS","ELITE","VELOS","NONE"]
            error("Instrument type $instrument_type not found for model $model_name")
        end
    elseif model_name == "AlphaPeptDeep"
        if instrument_type ∉ ["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"]
            error("Instrument type $instrument_type not found for model $model_name")
        end
    end
    #Model urls
    models_dict = Dict(
        "unispec"=>"https://koina.wilhelmlab.org:443/v2/models/UniSpec/infer",
        "prosit_2020_hcd"=>"https://koina.wilhelmlab.org:443/v2/models/Prosit_2020_intensity_HCD/infer",
        "AlphaPeptDeep"=>"https://koina.wilhelmlab.org:443/v2/models/AlphaPeptDeep_ms2_generic/infer"
    )
    model_url = nothing
    #Get command for the url and data 
    if haskey(models_dict, model_name)
        model_url = models_dict[model_name]
    else
        throw(ErrorException("The model name '$model_name' isn't valid. Valid models are "*join(keys(models_dict),';')))
    end

    #Format request json string
    function prepareBatch(
        peptide_table::DataFrame,
        model_type::InstrumentSpecificModel)
        nrows = size(peptide_table, 1)
        data_batch = Dict(
            "id" => "0",
            "inputs" => [
                Dict("name" => "peptide_sequences", "shape" => [nrows,1], "datatype" => "BYTES", "data" => getModifiedSequence.(peptide_table[!,:sequence], "", peptide_table[!,:mods])),
                Dict("name" => "precursor_charges", "shape" => [nrows,1], "datatype" => "INT32", "data" => peptide_table[!,:precursor_charge]),
                Dict("name" => "collision_energies", "shape" => [nrows,1], "datatype" => "FP32", "data" =>  peptide_table[!,:collision_energy]),
                Dict("name" => "instrument_types", "shape" => [nrows,1], "datatype" => "BYTES", "data" => [instrument_type for _ in range(1, nrows)])
            ]
        )

        return JSON.json(data_batch)
    end

    function prepareBatch(
        peptide_table::DataFrame,
        model_type::InstrumentAgnosticModel)
        nrows = size(peptide_table, 1)
        data_batch = Dict(
            "id" => "0",
            "inputs" => [
                Dict("name" => "peptide_sequences", "shape" => [nrows,1], "datatype" => "BYTES", "data" => getModifiedSequence.(peptide_table[!,:sequence], "", peptide_table[!,:mods])),
                Dict("name" => "precursor_charges", "shape" => [nrows,1], "datatype" => "INT32", "data" => peptide_table[!,:precursor_charge]),
                Dict("name" => "collision_energies", "shape" => [nrows,1], "datatype" => "FP32", "data" =>  peptide_table[!,:collision_energy]),
            ]
        )

        return JSON.json(data_batch)
    end


    function getKoinaBatches(
        peptides_df::DataFrame,
        model_type::KoinaModelType,
        model_url::String,
        batch_size::Int,
        first_prec_idx::UInt32)
        batch_size = min(batch_size, 1000)
        nprecs = size(peptides_df, 1)
        batch_start_idxs = collect(one(UInt32):UInt32(batch_size):UInt32(nprecs))
        batch_json_inputs = [prepareBatch(peptides_df[start_idx:min((start_idx + batch_size - 1), nprecs),:], 
        model_type) for start_idx in batch_start_idxs]
        #return batch_json_inputs[1]
        out_requests = multi_request(batch_json_inputs, model_url)
        #return out_requestsd
        output_dfs = [parseBatchToTable(request["outputs"]) for request in out_requests]
        for (i, df_frags_per_prec_pair) in enumerate(output_dfs)
            df, frags_per_prec = df_frags_per_prec_pair
            addPrecursorID!(
                df,
                UInt32(first_prec_idx + batch_start_idxs[i] - 1), #First precursor idx
                frags_per_prec
            )
            filterEachPrecursor!(df, intensity_threshold = intensity_threshold)
            sortEachPrecursor!(
                    df, 
                    :intensities,
                    )
        end
        return vcat([first(x) for x in output_dfs]...)
    end

    peptides_df = DataFrame(Arrow.Table(peptide_table_path))
    nprecs = size(peptides_df, 1)
    batch_size = min(batch_size, 1000)
    batch_start_idxs = collect(one(UInt32):UInt32(batch_size*max_koina_batches):UInt32(nprecs))
    for start_idx in ProgressBar(batch_start_idxs)
        stop_idx = min(start_idx + batch_size*max_koina_batches - 1, size(peptides_df, 1))
        frags_out = getKoinaBatches(peptides_df[start_idx:stop_idx,:],
                        model_type,
                        model_url,
                        batch_size,
                        start_idx)
        #return nothing
        Arrow.append(
           frags_out_path,
            frags_out
        )
    end
end

function predictFragments(
    frags_out_path::String, #"/Users/n.t.wamsley/Desktop/test_out.arrow",
    altimeter_dir::String,
    intensity_threshold::Float32 = 0.001f0)

    println("Getting json fpaths...")
    ordered_altimeter_json_paths = joinpath.(
        altimeter_dir,
        sort(
            readdir(altimeter_dir), 
            by = x->parse(Int64, String(split(split(x,'_')[end], '.')[1])) #Names in format of "Altimiter_####.json"
            )
    )
    precursor_idx = one(UInt32)
    println("Reading Altimeter Outputs...")
    batch_dfs = DataFrame()
    batch_counter = 0
    for (fid, batch_path) in ProgressBar(enumerate(ordered_altimeter_json_paths))
        df, frags_per_prec = parseBatchToTable(JSON.parse(read(batch_path, String))["outputs"])
        addPrecursorID!(
            df, 
            UInt32(precursor_idx),
            frags_per_prec
        ) 
        precursor_idx = df[end,:precursor_idx] + one(UInt32)
        filterEachPrecursor!(df, intensity_threshold = intensity_threshold)
        sortEachPrecursor!(
                df, 
                :intensities,
                )
        batch_dfs = vcat(batch_dfs, df)
        batch_counter += 1
        if batch_counter > 500
            println("Batch finished...")
            Arrow.append(frags_out_path,
                    batch_dfs)
            batch_dfs = nothing
            df = nothing
            GC.gc()
            batch_dfs = DataFrame()
            batch_counter = 0
        end
    end
    Arrow.append(frags_out_path,
    batch_dfs)

end