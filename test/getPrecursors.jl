
using Dictionaries 
fixed_mods = [(p=r"C", r="C[Carb]")]
var_mods = [(p=r"(K$)", r="[Hlys]"), (p=r"(R$)", r="[Harg]")]
input_string = "PEPTICKDECK"

test_mods::Dict{String, Float32} = 
Dict{String, Float32}(
    "Carb" => Float32(57.021464),
    "Harg" => Float32(10),
    "Hlys" => Float32(8)
)
precursor_table = PrecursorTable()
using Combinatorics



