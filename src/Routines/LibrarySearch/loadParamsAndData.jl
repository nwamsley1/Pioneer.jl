println("Parsing Arguments...")
ARGS = Dict(
    "params_json"=>"data/example_config/LibrarySearch.json"
)
#ARGS = parse_commandline();
#=
library_fragment_lookup_table.prec_frag_ranges[end] = 0x29004baf:(0x29004be8 - 1)
3-protome PC test March 4th 2024 for hupo
julia --threads 15 ./src/Routines/LibrarySearch/MAIN.jl ./data/example_config/LibrarySearch.json /Users/n.t.wamsley/TEST_DATA/PXD046444/arrow/exploris_test /Users/n.t.wamsley/TEST_DATA/SPEC_LIBS/HUMAN/STANDARD_NCE33_DefCharge2_DYNAMIC/PIONEER/LIBA/UP000005640_9606_Apr20_24/pioneer_lib -s true 

#params = JSON.parse(read("/Users/n.t.wamsley/RIS_temp/ASMS_2024/ASTRAL_THREE_PROTEOME/unispec_chronologer_1mc_1var_by_exploris27_02_052724/OAT_103ISO/LibrarySearch.json", String));
params = JSON.parse(read("data/example_config/LibrarySearch.json", String));

#SPEC_LIB_DIR ="/Users/n.t.wamsley/RIS_temp/ASMS_2024/ASTRAL_THREE_PROTEOME/unispec_chronologer_1mc_1var_by_exploris27_02_052724/spec_lib/pioneer_lib"
SPEC_LIB_DIR ="/Users/n.t.wamsley/RIS_temp/ASMS_2024/ASTRAL_THREE_PROTEOME/unispec_chronologer_1mc_1var_by_052724/spec_lib/pioneer_lib"
MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/PXD046444/arrow/astral_test"


MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
EXPERIMENT_NAME = "EXPERIMENT"
=#

params = JSON.parse(read(ARGS["params_json"], String));
EXPERIMENT_NAME = "EXPERIMENT"
params = JSON.parse(read("data/example_config/LibrarySearch.json", String));
MS_DATA_DIR = params["ms_data_dir"];
SPEC_LIB_DIR = params["library_folder"];

#Get all files ending in ".arrow" that are in the MS_DATA_DIR folder. 
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
MS_TABLE_PATH_TO_ID = Dictionary(MS_TABLE_PATHS, UInt32.(collect(range(1,length(MS_TABLE_PATHS)))))
MS_TABLE_ID_TO_PATH = Dictionary(UInt32.(collect(range(1,length(MS_TABLE_PATHS)))), MS_TABLE_PATHS)

params_ = (
    expected_matches = Int64(params["expected_matches"]),
    isotope_err_bounds = Tuple([Int64(bound) for bound in params["isotope_err_bounds"]]),
    choose_most_intense = Bool(params["choose_most_intense"]),
    quadrupole_isolation_width = Float64(params["quadrupole_isolation_width"]),
    irt_err_sigma = params["irt_err_sigma"],


    presearch_params = Dict{String, Any}(k => v for (k, v) in params["presearch_params"]),
    first_search_params = Dict{String, Any}(k => v for (k, v) in params["first_search_params"]),
    quant_search_params = Dict{String, Any}(k => v for (k, v) in params["quant_search_params"]),
    frag_tol_params = Dict{String, Any}(k => v for (k, v) in params["frag_tol_params"]),
    irt_mapping_params = Dict{String, Any}(k => v for (k, v) in params["irt_mapping_params"]),
    integration_params = Dict{String, Any}(k => v for (k, v) in params["integration_params"]),
    deconvolution_params = Dict{String, Any}(k => v for (k, v) in params["deconvolution_params"]),
    summarize_first_search_params = Dict{String, Any}(k => v for (k, v) in params["summarize_first_search_params"]),
    qc_plot_params = Dict{String, Any}(k => v for (k, v) in params["qc_plot_params"]),
    normalization_params = Dict{String, Any}(k => v for (k, v) in params["normalization_params"]),
    benchmark_params = Dict{String, Any}(k => v for (k, v) in params["benchmark_params"])
    );


MS_DATA_DIR = joinpath(params_[:benchmark_params]["results_folder"], EXPERIMENT_NAME);
out_folder = joinpath(MS_DATA_DIR, "Search")
RESULTS_DIR = MS_DATA_DIR
if !isdir(out_folder)
    mkpath(out_folder)
end
qc_plot_folder = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS")
if !isdir(qc_plot_folder)
    mkpath(qc_plot_folder)
else
    [rm(joinpath(qc_plot_folder, x)) for x in readdir(qc_plot_folder) if endswith(x, ".pdf")]
end
const rt_alignment_folder = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS","rt_alignment")
if !isdir(rt_alignment_folder)
    mkpath(rt_alignment_folder)
end
const mass_err_estimation_folder = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS","mass_error_estimation")
if !isdir(mass_err_estimation_folder)
    mkpath(mass_err_estimation_folder)
end
const results_folder = joinpath(MS_DATA_DIR, "Search", "RESULTS")
if !isdir(results_folder)
    mkpath(results_folder)
end
const params_folder = joinpath(MS_DATA_DIR, "Search", "PARAMS")
if !isdir(params_folder )
    mkpath(params_folder)
end
#Write params to folder 
open(joinpath(params_folder, "config.json"),"w") do f 
    JSON.print(f, params)
end



##########
#Load Dependencies 
##########
#Fragment Library Parsing
#=
[include(joinpath(pwd(), "src", "Structs", jl_file)) for jl_file in [
                                                                    "ChromObject.jl",
                                                                    "ArrayDict.jl",
                                                                    "Counter.jl",
                                                                    "Ion.jl",
                                                                    "LibraryIon.jl",
                                                                    "MatchIon.jl",
                                                                    "LibraryFragmentIndex.jl",
                                                                    "SparseArray.jl"]];

#Utilities
[include(joinpath(pwd(), "src", "Utils", jl_file)) for jl_file in [

                                                                    "isotopes.jl",
                                                                    "globalConstants.jl",
                                                                    "uniformBasisCubicSpline.jl",
                                                                    "isotopeSplines.jl",
                                                                    "normalizeQuant.jl",
                                                                    "massErrorEstimation.jl",
                                                                    "SpectralDeconvolution.jl",
                                                                    "percolatorSortOf.jl",
                                                                    "plotRTAlignment.jl",
                                                                    "probitRegression.jl",
                                                                    "partitionThreadTasks.jl",
                                                                    "LFQ.jl",
                                                                    "scoreProteinGroups.jl",
                                                                    "wittakerHendersonSmoothing.jl",
                                                                    "getBestTrace.jl"]];

[include(joinpath(pwd(), "src", "PSM_TYPES", jl_file)) for jl_file in ["PSM.jl","spectralDistanceMetrics.jl","UnscoredPSMs.jl","ScoredPSMs.jl"]];

#Files needed for PRM routines
[include(joinpath(pwd(), "src", "Routines","LibrarySearch","methods",jl_file)) for jl_file in [
                                                                                    "matchPeaks.jl",
                                                                                    "buildDesignMatrix.jl",
                                                                                    "manipulateDataFrames.jl",
                                                                                    "buildRTIndex.jl",
                                                                                    "searchRAW.jl",
                                                                                    "selectTransitions.jl",
                                                                                    "integrateChroms.jl",
                                                                                    "queryFragmentIndex.jl"]];
                                             

=#
###########
#Load Spectral Libraries
###########
f_index_fragments = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_fragments.arrow"))
f_index_rt_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_rt_bins.arrow"))
f_index_frag_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "f_index_fragment_bins.arrow"))


presearch_f_index_fragments = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_fragments.arrow"))
presearch_f_index_rt_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_rt_bins.arrow"))
presearch_f_index_frag_bins = Arrow.Table(joinpath(SPEC_LIB_DIR, "presearch_f_index_fragment_bins.arrow"))


println("Loading spectral libraries into main memory...")
prosit_lib = Dict{String, Any}()
detailed_frags = load(joinpath(SPEC_LIB_DIR,"detailed_fragments.jld2"))["detailed_fragments"]
prec_frag_ranges = load(joinpath(SPEC_LIB_DIR,"precursor_to_fragment_indices.jld2"))["precursor_to_fragment_indices"]
const library_fragment_lookup_table = LibraryFragmentLookup(detailed_frags, prec_frag_ranges)
last_range = library_fragment_lookup_table.prec_frag_ranges[end] #0x29004baf:(0x29004be8 - 1)
last_range = range(first(last_range), last(last_range) - 1)
library_fragment_lookup_table.prec_frag_ranges[end] = last_range
prosit_lib["f_det"] = library_fragment_lookup_table

precursors = Arrow.Table(joinpath(SPEC_LIB_DIR, "precursor_table.arrow"))#DataFrame(precursors)
f_index = FragmentIndex(
    f_index_frag_bins[:FragIndexBin],
    f_index_rt_bins[:FragIndexBin],
    f_index_fragments[:IndexFragment],
);
presearch_f_index = FragmentIndex(
    presearch_f_index_frag_bins[:FragIndexBin],
    presearch_f_index_rt_bins[:FragIndexBin],
    presearch_f_index_fragments[:IndexFragment],
);
prosit_lib["f_index"] = f_index;
prosit_lib["presearch_f_index"] = presearch_f_index;
prosit_lib["precursors"] = precursors;

###########
#Set CV Folds 
###########

const pid_to_cv_fold = getCVFolds(
    collect(range(UInt32(1), UInt32(length(precursors[:sequence])))),#precursor id's, 
    precursors[:accession_numbers]
    )

###########
#Load Pre-Allocated Data Structures. One of each for each thread. 
###########
N = Threads.nthreads()
M = 250000
n_precursors = length(precursors[:mz])
ionMatches = [[FragmentMatch{Float32}() for _ in range(1, M)] for _ in range(1, N)];
ionMisses = [[FragmentMatch{Float32}() for _ in range(1, M)] for _ in range(1, N)];
all_fmatches = [[FragmentMatch{Float32}() for _ in range(1, 1000000)] for _ in range(1, N)];
IDtoCOL = [ArrayDict(UInt32, UInt16, n_precursors ) for _ in range(1, N)];
ionTemplates = [[DetailedFrag{Float32}() for _ in range(1, M)] for _ in range(1, N)];
iso_splines = parseIsoXML("./data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml");
scored_PSMs = [Vector{SimpleScoredPSM{Float32, Float16}}(undef, 5000) for _ in range(1, N)];
unscored_PSMs = [[SimpleUnscoredPSM{Float32}() for _ in range(1, 5000)] for _ in range(1, N)];
spectral_scores = [Vector{SpectralScoresSimple{Float16}}(undef, 5000) for _ in range(1, N)];
precursor_weights = [zeros(Float32, n_precursors ) for _ in range(1, N)];
precs = [Counter(UInt32, UInt8,n_precursors ) for _ in range(1, N)];
chromatograms = [Vector{ChromObject}(undef, 5000) for _ in range(1, N)];
complex_scored_PSMs = [Vector{ComplexScoredPSM{Float32, Float16}}(undef, 5000) for _ in range(1, N)];
complex_unscored_PSMs = [[ComplexUnscoredPSM{Float32}() for _ in range(1, 5000)] for _ in range(1, N)];
complex_spectral_scores = [Vector{SpectralScoresComplex{Float16}}(undef, 5000) for _ in range(1, N)];

###########
#File Names Parsing 
###########
file_names = first.(split.(basename.(MS_TABLE_PATHS), '.'))
split_file_names = split.(file_names, "_")



#If file names have inconsistnat number of delimiters, give up on parsing and use the entire file name
unique_split_file_name_lengths = unique(length.(split_file_names))
parsed_file_names = copy(split_file_names)
M = length(parsed_file_names)

if length(unique_split_file_name_lengths) == 1
    N = first(unique_split_file_name_lengths)
    M = length(file_names)
    split_strings = Array{String}(undef, (M, N))
    for i in 1:N
        for j in 1:M
            split_strings[j, i] = split_file_names[j][i]
        end
    end
    cols_to_keep = zeros(Bool, N)
    for (col_idx, col) in enumerate(eachcol(split_strings))
        if length(unique(col))>1
            cols_to_keep[col_idx] = true
        end
    end
    parsed_file_names = Vector{String}(undef, M)
    split_strings = split_strings[:, cols_to_keep]
    for i in 1:M
        parsed_file_names[i] = join(split_strings[i,:], "_")
    end
else
    parsed_file_names = file_names
end

const file_id_to_parsed_name = Dict(zip(1:M, [string(x) for x in parsed_file_names]))
#Parsed file names 
const parsed_fnames = sort(collect(values(file_id_to_parsed_name)))


