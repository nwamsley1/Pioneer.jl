using JSON

# Prepare the data
data = Dict(
    "id" => "0",
    "inputs" => [
        Dict("name" => "peptide_sequences", "shape" => [2,1], "datatype" => "BYTES", "data" => ["YYHTLFTHSLPK", "AGC[UNIMOD:4]FS[UNIMOD:27]PK"]),
        Dict("name" => "precursor_charges", "shape" => [2,1], "datatype" => "INT32", "data" => [1,2]),
        Dict("name" => "collision_energies", "shape" => [2,1], "datatype" => "FP32", "data" => [25, 25]),
        Dict("name" => "instrument_types", "shape" => [2,1], "datatype" => "BYTES", "data" => ["LUMOS", "QE"])
    ]
)

# Convert the data to JSON
json_data = JSON.json(data)

# Construct the curl command
cmd = `curl -s "https://koina.wilhelmlab.org:443/v2/models/UniSpec/infer" -d $json_data`

# Run the command and capture the output
output = JSON.parse(read(cmd, String))
