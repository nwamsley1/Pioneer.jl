using CSV, DataFrames, Dictionaries, ProgressBars, JLD2, Arrow, JSON
using Base.Iterators: partition

include("src/structs/Ion.jl")
include("src/structs/LibraryIon.jl")
include("src/structs/LibraryFragmentIndex.jl")
include("src/Routines/buildFragmentIndex/buildFragmentIndex.jl")

ARGS_ = Dict(
    "params_json" =>  joinpath(pwd(), "data", "example_config","buildFragmentIndex.json"),
    "pioneer_lib_dir" => joinpath("/Users","n.t.wamsley","TEST_DATA",
    "SPEC_LIBS","HUMAN",
    "STANDARD_NCE33_DefCharge2_DYNAMIC","PIONEER","LIBA",
    "UP000005640_9606_NCE29dynamic_DefCharge2_1missedCleavage_1varMod_8to30_wDecoys_fromChronologer_hi.arrow"),
    "out_dir" => joinpath("/Users","n.t.wamsley","TEST_DATA",
    "SPEC_LIBS","HUMAN",
    "STANDARD_NCE33_DefCharge2_DYNAMIC","PIONEER","LIBA","UP000005640_9606_Apr20_24")

)


ARGS_ = Dict(
    "params_json" =>  joinpath("/Users","n.t.wamsley","RIS_temp","ASMS_2024",
                            "UNISPEC_LIBS","THREE_PROTEOME_052224","ASTRAL_NCE25",
                            "buildFragmentIndex.json"),

    "pioneer_lib_dir" => joinpath("/Users","n.t.wamsley","RIS_temp","ASMS_2024",
                            "UNISPEC_LIBS","THREE_PROTEOME_052224","ASTRAL_NCE25",
                            "test_lib_pioneer_lib.arrow"),
    "out_dir" => joinpath("/Users","n.t.wamsley","RIS_temp","ASMS_2024",
    "UNISPEC_LIBS","THREE_PROTEOME_052224","ASTRAL_NCE25","THREE_PROTEOME_052324")

)


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "params_json"
            help = "Path to a .json file with the parameters"
            required = true
        "pioneer_lib_dir"
            help = "path to a pioneer-formated spectral library (.arrow format)"
            required = true
        "out_dir"
            help = "directory in which to place output. Will make if it does not alredy exist"
            required = true
    end

    return parse_args(s)
end

println("Parsing Arguments...")
ARGS = parse_commandline();
params = JSON.parse(read(ARGS_["params_json"], String));
params_ = (
    rt_bin_tol = Float32(params["rt_bin_tol"]),
    frag_bin_tol_ppm = Float32(params["frag_bin_tol_ppm"]),
    frag_mz_min = Float32(params["frag_mz_min"]),
    frag_mz_max = Float32(params["frag_mz_max"]),
    prec_mz_min = Float32(params["prec_mz_min"]),
    prec_mz_max = Float32(params["prec_mz_max"]),
    y_start_index = Int64(params["y_start_index"]),
    b_start_index = Int64(params["b_start_index"]),
    y_start = Int64(params["y_start"]),
    b_start = Int64(params["b_start"]),
    missed_cleavage_regex = Regex(params["missed_cleavage_regex"]),
    rank_to_score = UInt8.(params["rank_to_score"])
)

println("make output folder...")
folder_out = joinpath(ARGS_["out_dir"],"pioneer_lib")
if !isdir(folder_out)
    mkpath(folder_out)
end

#Write params to output folder for record-keeping 
open(joinpath(folder_out, "config.json"),"w") do f 
    JSON.print(f, params)
end

println("reading pioneer lib...")

ARGS_["pioneer_lib_dir"]
#@time pioneer_lib = Arrow.Table("/Users/n.t.wamsley/Desktop/test_lib.arrow")
@time pioneer_lib = Arrow.Table(ARGS_["pioneer_lib_dir"])
#@time pioneer_lib = Arrow.Table(ARGS_["pioneer_lib_dir"])
@time prec_frag_ranges = load(joinpath("/Users","n.t.wamsley","RIS_temp","ASMS_2024",
                            "UNISPEC_LIBS","THREE_PROTEOME_052224","ASTRAL_NCE25",
                            "test_lib_prec_frag_ranges.jld2"))["prec_frag_ranges"];

@time prec_frags = load(joinpath("/Users","n.t.wamsley","RIS_temp","ASMS_2024",
                            "UNISPEC_LIBS","THREE_PROTEOME_052224","ASTRAL_NCE25",
                            "test_lib_prec_frags.jld2"))["prec_frags"];

@time sorted_indices = sortPrositLib(
    pioneer_lib[:mz],
    pioneer_lib[:irt],
    params_[:rt_bin_tol]
)
#using PProf, Profile
#Profile.clear()

#frags_vec, prec_to_frag_range = getFragRanges(pioneer_lib[:frags])
#@profile 
simple_fragments, folder_out = parsePioneerLib(
    folder_out,
    pioneer_lib[:irt],
    pioneer_lib[:mz],
    pioneer_lib[:decoy],

    pioneer_lib[:proteome_name],
    pioneer_lib[:protein_name],
    pioneer_lib[:sequence],
    pioneer_lib[:structural_mods],
    pioneer_lib[:isotope_mods],

    pioneer_lib[:charge],
    pioneer_lib[:sulfur_count],
    pioneer_lib[:seq_length],
    prec_frag_ranges,
    prec_frags,
    params_[:frag_mz_min], 
    params_[:frag_mz_max],
    params_[:prec_mz_min],
    params_[:prec_mz_max],
    rank_to_score=params_[:rank_to_score],
    rt_bin_tol = params_[:rt_bin_tol],
    y_start_index = params_[:y_start_index],
    b_start_index = params_[:b_start_index],
    y_start = params_[:y_start],
    b_start = params_[:b_start],
    missed_cleavage_regex = params_[:missed_cleavage_regex]
);
#pprof(;webport=58603)

#println("building fragment index...")
buildFragmentIndex!(
                    folder_out,
                    simple_fragments,
                    params_[:frag_bin_tol_ppm],
                    params_[:rt_bin_tol],
                    )
#Index with maximally sized RT bin for presearch 
buildFragmentIndex!(
    folder_out,
    simple_fragments,
    params_[:frag_bin_tol_ppm],
    typemax(Float32),
    index_name = "presearch_"
    )
    


"/Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/UP000005640_9606_Apr4_24/pioneer_lib"