using CSV, DataFrames, Dictionaries, ProgressBars, JLD2, Arrow, JSON
using Base.Iterators: partition

include("src/structs/Ion.jl")
include("src/structs/LibraryIon.jl")
include("src/structs/LibraryFragmentIndex.jl")
include("src/Routines/buildFragmentIndex/fragBounds.jl")
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
    "params_json" =>  joinpath("../","spec_lib",
                            "buildFragmentIndex.json"),

    "pioneer_lib_dir" => joinpath("../","spec_lib"),
    "out_dir" => joinpath("../","spec_lib")

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
    rt_bin_tol = Float64(params["rt_bin_tol"]),
    frag_bin_tol_ppm = Float32(params["frag_bin_tol_ppm"]),
    frag_mz_min = Float32(params["frag_mz_min"]),
    frag_mz_max = Float32(params["frag_mz_max"]),
    auto_detect_frag_bounds = parse(Bool, params["auto_detect_frag_bounds"]),
    frag_bound_detection_raw_file = params["frag_bounds_detection_raw_file"],
    prec_mz_min = Float32(params["prec_mz_min"]),
    prec_mz_max = Float32(params["prec_mz_max"]),
    exclude_types_from_index = Set([Char(first(x)) for x in params["exclude_from_index"]]),
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

PIONEER_LIB_DIR = ARGS_["pioneer_lib_dir"]
pioneer_lib_path = [joinpath(PIONEER_LIB_DIR , file) for file in filter(file -> isfile(joinpath(PIONEER_LIB_DIR , file)) && match(r"prec_table.arrow", file) != nothing, readdir(PIONEER_LIB_DIR ))][1];
id_to_annotation_path = [joinpath(PIONEER_LIB_DIR , file) for file in filter(file -> isfile(joinpath(PIONEER_LIB_DIR , file)) && match(r"id_to_annotation.arrow", file) != nothing, readdir(PIONEER_LIB_DIR ))][1];
frags_table_path = [joinpath(PIONEER_LIB_DIR , file) for file in filter(file -> isfile(joinpath(PIONEER_LIB_DIR , file)) && match(r"frag_table.arrow", file) != nothing, readdir(PIONEER_LIB_DIR ))][1];
frag_ranges_path = [joinpath(PIONEER_LIB_DIR , file) for file in filter(file -> isfile(joinpath(PIONEER_LIB_DIR , file)) && match(r"frag_ranges.arrow", file) != nothing, readdir(PIONEER_LIB_DIR ))][1];


@time pioneer_lib = Arrow.Table(pioneer_lib_path)
@time id_to_annotation = Arrow.Table(id_to_annotation_path)
@time frags_table = Arrow.Table(frags_table_path)
@time frag_ranges_table = Arrow.Table(frag_ranges_path)
frag_ranges = Vector{UnitRange{UInt32}}(undef, length(frag_ranges_table[:start]))
for i in range(1, length(frag_ranges))
    frag_ranges[i] = range(UInt32(frag_ranges_table[:start][i]), UInt32(frag_ranges_table[:stop][i]))
end
id_to_annotation = [x for x in id_to_annotation[:id_to_annotation]]
#@time pioneer_lib = Arrow.Table(ARGS_["pioneer_lib_dir"])
#=
@time prec_frag_ranges = load(joinpath("/Users","n.t.wamsley","RIS_temp","ASMS_2024",
                            "UNISPEC_LIBS","THREE_PROTEOME_052224","ASTRAL_NCE25",
                            "test_lib_prec_frag_ranges.jld2"))["prec_frag_ranges"];

@time prec_frags = load(joinpath("/Users","n.t.wamsley","RIS_temp","ASMS_2024",
                            "UNISPEC_LIBS","THREE_PROTEOME_052224","ASTRAL_NCE25",
                            "test_lib_prec_frags.jld2"))["prec_frags"];
=#
#=
@time sorted_indices = sortPrositLib(
    pioneer_lib[:mz],
    pioneer_lib[:irt],
    params_[:rt_bin_tol]
)
=#
frag_bounds, prec_mz_min, prec_mz_max = nothing, nothing, nothing
if params_[:auto_detect_frag_bounds]

    MS_TABLE = Arrow.Table(params_[:frag_bound_detection_raw_file])
    frag_bounds, prec_mz_min, prec_mz_max = getFragBounds(
        MS_TABLE[:centerMass],
        MS_TABLE[:isolationWidth],
        MS_TABLE[:msOrder],
        MS_TABLE[:lowMass],
        MS_TABLE[:highMass])
    prec_mz_min -= 1.0f0
    prec_mz_max += 1.0f0
else
    frag_bounds = FragBoundModel(
        ImmutablePolynomial([params_[:frag_mz_min]]),
        ImmutablePolynomial(params_[:frag_mz_max]) 
    )
    prec_mz_min, prec_mz_max = params_[:prec_mz_min], params_[:prec_mz_max]
end



simple_fragments, folder_out = parsePioneerLib(
    folder_out,
    id_to_annotation,
    pioneer_lib[:irt],
    pioneer_lib[:mz],
    pioneer_lib[:decoy],

    collect(frags_table[:mz]),
    collect(frags_table[:intensity]),
    collect(frags_table[:ion_type]),
    collect(frags_table[:is_y]),
    collect(frags_table[:frag_index]),
    collect(frags_table[:charge]),
    collect(frags_table[:isotope]),
    collect( frags_table[:internal]),
    collect(frags_table[:immonium]),
    collect(frags_table[:internal_ind]),
    collect(frags_table[:sulfur_count]),

    pioneer_lib[:proteome_name],
    pioneer_lib[:protein_name],
    pioneer_lib[:sequence],
    pioneer_lib[:structural_mods],
    pioneer_lib[:isotope_mods],

    pioneer_lib[:charge],
    pioneer_lib[:sulfur_count],
    pioneer_lib[:seq_length],
    frag_ranges,
    #prec_frags,
    frag_bounds,
    prec_mz_min, 
    prec_mz_max,
    rank_to_score=params_[:rank_to_score],
    rt_bin_tol = params_[:rt_bin_tol],
    exclude_from_index = params_[:exclude_types_from_index],
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

prec_table_second= Arrow.Table("/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/THREE_PROTEOME_052224/ASTRAL_NCE25/THREE_PROTEOME_ASTRAL_052624/pioneer_lib/precursor_table.arrow")
prec_frags_second = load("/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/THREE_PROTEOME_052224/ASTRAL_NCE25/THREE_PROTEOME_ASTRAL_052624/pioneer_lib/detailed_fragments.jld2")["detailed_fragments"]
prec_to_frag_ranges_second = load("/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/THREE_PROTEOME_052224/ASTRAL_NCE25/THREE_PROTEOME_ASTRAL_052624/pioneer_lib/precursor_to_fragment_indices.jld2")["precursor_to_fragment_indices"]

prec_table_first = Arrow.Table("/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/THREE_PROTEOME_052224/ASTRAL_NCE25/THREE_PROTEOME_ASTRAL_by_052424/pioneer_lib/precursor_table.arrow")
prec_frags_first = load("/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/THREE_PROTEOME_052224/ASTRAL_NCE25/THREE_PROTEOME_ASTRAL_by_052424/pioneer_lib/detailed_fragments.jld2")["detailed_fragments"]
prec_to_frag_ranges_first = load("/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/THREE_PROTEOME_052224/ASTRAL_NCE25/THREE_PROTEOME_ASTRAL_by_052424/pioneer_lib/precursor_to_fragment_indices.jld2")["precursor_to_fragment_indices"]


N = 1
size(prec_table_first[:sequence])
prec_table_first[:sequence][N]
size(prec_table_second[:sequence])
prec_table_second[:sequence][N]

prec_frags_first[prec_to_frag_ranges_first[N]]
prec_frags_second[prec_to_frag_ranges_second[N]]

N += 1


prec_table = Arrow.Table("/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/THREE_PROTEOME_052224/ASTRAL_NCE25/THREE_PROTEOME_ASTRAL_052624/prec_table.arrow")
frag_range_table = Arrow.Table("/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/THREE_PROTEOME_052224/ASTRAL_NCE25/THREE_PROTEOME_ASTRAL_052624/frag_range_table.arrow")
frag_table = Arrow.Table("/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/THREE_PROTEOME_052224/ASTRAL_NCE25/THREE_PROTEOME_ASTRAL_052624/frag_table.arrow")


frag_table[frag_range_table[:start][N]:frag_range_table[:stop][N],:]
