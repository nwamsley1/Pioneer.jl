function SearchDIA(params_path::String)
    #println("JLD2 version is: ", Pkg.installed()["JLD2"])
    total_time = @timed begin
    #params_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/OlsenMixedSpeciesAstral200ng/OlsenMixedAltimeterParamsAltimeter111824_SeperateTracesNoMax.json"
    #params_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/THREE_PROTEOME_5MIN/COMBINE_TRACES.json"
    params_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/DATASETS_ARROW/ASTRAL_MTAC/THREE_PROTEOME_5MIN/COMBINE_TRACES_2.json"
    if !isabspath(params_path)
        params_path = joinpath(@__DIR__, "../../", params_path)
    end
    params = JSON.parse(read(params_path, String));
    MS_DATA_DIR = params["ms_data_dir"];
    #MS_DATA_DIR = "/Users/n.t.wamsley/Desktop/FIRST_TRY_ASTRAL/arrow_out"
    SPEC_LIB_DIR = params["library_folder"];
    if !isabspath(SPEC_LIB_DIR)
        SPEC_LIB_DIR =  joinpath(@__DIR__, "../../", SPEC_LIB_DIR)
    end

    #Get all files ending in ".arrow" that are in the MS_DATA_DIR folder. 
    if isabspath(MS_DATA_DIR)
        MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
    else
        MS_DATA_DIR = joinpath(@__DIR__, "../../", MS_DATA_DIR)
        MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
    end

    params_ = parseParams(params)

    ###########
    #Load Spectral Libraries
    ###########
    #params_[:presearch_params]["nce_guess"] = 30.0f0
    spec_lib = loadSpectralLibrary(SPEC_LIB_DIR, params_);
    SPEC_LIB = FragmentIndexLibrary(
        spec_lib["presearch_f_index"], 
        spec_lib["f_index"], 
        BasicLibraryPrecursors(spec_lib["precursors"]), 
        spec_lib["f_det"])
    ###########
    #Load Pre-Allocated Data Structures. One of each for each thread. 
    ###########
      file_id_to_parsed_name, parsed_fnames,file_path_to_parsed_name = parseFileNames(MS_TABLE_PATHS)


    #MS_TABLE_PATHS = MS_TABLE_PATHS[end-2:end]

    SEARCH_CONTEXT = initSearchContext(
        SPEC_LIB,
        parseIsoXML(joinpath(@__DIR__,"../data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml")),
        ArrowTableReference(MS_TABLE_PATHS),
        Threads.nthreads(),
        250000  # buffer size
    );

    setDataOutDir!(SEARCH_CONTEXT, params_[:benchmark_params]["results_folder"]);

    include("Routines/LibrarySearch/SearchMethods/QuadTuningSearch.jl")
    include("Routines/LibrarySearch/SearchMethods/NceTuningSearch.jl")

    #==========================================================
    1) Estimate empirical (rt) to library retention time conversion (irt) 
        Results stored in SEARCH_CONTEXT.rt_to_irt_model::Dict{Int64, RtConversionModel}
    2) Estimate fragment mass tolerance (ppm) on left and right hand SimpleScoredPSM
        Results stored in SEARCH_CONTEXT.mass_error_model::Dict{Int64, MassErrorModel}
    ==========================================================#
    test = @timed begin
        @time execute_search(
            ParameterTuningSearch(),SEARCH_CONTEXT , params_
        );

        #==========================================================
        Estimate normalized collision energy for each raw file
            Normalized collision energy is modeled as a function of peptide m/z and Charge
            Results stored in SEARCH_CONTEXT.nce_model::Dict{Int64, NceModel}
        ==========================================================#
        include("Routines/LibrarySearch/SearchMethods/NceTuningSearch.jl")
        @time execute_search(
            NceTuningSearch(), SEARCH_CONTEXT, params_
        );

        #==========================================================
        Estimate quadrupole transmission function for each raw file
            Quadrupole transmision function is modeled as an assymtric
            generalized bell function, with a sensible default. 
            See src/utils/quadTransmissionModeling.
            Results stored in SEARCH_CONTEXT.quad_transmission_model::Dict{Int64, QuadTransmissionModel}
        ==========================================================#
        params_[:presearch_params]["quad_tuning_sample_rate"] = 0.02
        #params_[:presearch_params]["min_log2_matched_ratio"] = zero(Float32)#typemin(Float32)
        #params_[:presearch_params]["min_quad_tuning_psms"] = 3000
        include("utils/quadTransmissionModeling/binIsotopeRatioData.jl")
        include("Routines/LibrarySearch/SearchMethods/QuadTuningSearch.jl")
        @time execute_search(
            QuadTuningSearch(), SEARCH_CONTEXT, params_
        );
    
        quad_model_dict = SEARCH_CONTEXT.quad_transmission_model
        p = plot()
        plot_bins = LinRange(400-3, 400+3, 100)
        for (key, value) in pairs(quad_model_dict)
        
           quad_func = getQuadTransmissionFunction(value, 400.0f0, 2.0f0)
           plot!(p, plot_bins, quad_func.(plot_bins), lw = 2, alpha = 0.5,
           #label =  getFileIdToName(getMSData(SEARCH_CONTEXT),key), 
           legend = :outertopleft)
        end
        p
        
        #==========================================================
        1) Estimate empirical (rt) to library retention time conversion (irt) 
            Results stored in SEARCH_CONTEXT.rt_to_irt_model::Dict{Int64, RtConversionModel}
        2) Revised estimate of rt->irt mapping. Also models the inverse irt->rt. Uses 
        psms at 1%fdr. Results stored in SEARCH_CONTEXT.irt_rt_map::Dict{Int64, RtConversionModel},
        and SEARCH_CONTEXT.irt_rt_map::Dict{Int64, RtConversionModel}
        3) Writes top N peptide spectrum matches for each raw file in SEARCH_CONTEXT.mass_spec_data_reference, 
            to the `first_pass_psms` folder as .arrow tables. 
        3) Summarize search results. Accross all raw files in SEARCH_CONTEXT.mass_spec_data_reference, 
            identify the top
        ==========================================================#
        @time execute_search(
            FirstPassSearch(), SEARCH_CONTEXT, params_
        );

        #==========================================================
        Estimate huber (Charboneir) loss smoothing parameter δ. Increasing δ
            behaves like squared error and decreasing δ behaves like absolute error. 
            Used for linear regression of MS/MS spectra onto library spectra. 
            Results stored in SEARCH_CONTEXT.huber_delta::Base.Ref{Float32}
        ==========================================================#
        include("Routines/LibrarySearch/SearchMethods/HuberTuningSearch.jl")
        @time execute_search(
            HuberTuningSearch(), SEARCH_CONTEXT, params_
        );

        #==========================================================
        Estimate huber (Charboneir) loss smoothing parameter δ. Increasing δ
            behaves like squared error and decreasing δ behaves like absolute error. 
            Used for linear regression of MS/MS spectra onto library spectra. 
            Results stored in SEARCH_CONTEXT.huber_delta::Base.Ref{Float32}
        ==========================================================#
        include("Routines/LibrarySearch/SearchMethods/SecondPassSearch.jl")
        @time execute_search(
            SecondPassSearch(), SEARCH_CONTEXT, params_
        );

        include("Routines/LibrarySearch/SearchMethods/ScoringSearch.jl")
        @time execute_search(
            ScoringSearch(), SEARCH_CONTEXT, params_
        );

        include("Routines/LibrarySearch/SearchMethods/IntegrateChromatogramsSearch.jl")
        @time execute_search(
            IntegrateChromatogramSearch(), SEARCH_CONTEXT, params_
        );

        include("Routines/LibrarySearch/SearchMethods/MaxLFQSearch.jl")
        @time execute_search(
            MaxLFQSearch(), SEARCH_CONTEXT, params_
        );
    end
    ###########
    #Normalize Quant 
    ###########
    println("Normalization...")
    normalizeQuant(
            second_quant_folder,
            :peak_area,
            N = params_[:normalization_params]["n_rt_bins"],
            spline_n_knots = params_[:normalization_params]["spline_n_knots"]
        )

    merged_second_quant_path = joinpath(temp_folder, "joined_second_quant.arrow")
    if isfile(merged_second_quant_path)
        rm(merged_second_quant_path, force = true)
    end
    if isfile(joinpath(results_folder, "precursors_long.arrow"))
        rm(joinpath(results_folder, "precursors_long.arrow"), force=true)
    end
    mergeSortedArrowTables(
        second_quant_folder,
        joinpath(results_folder, "precursors_long.arrow"),
        (:protein_idx,:precursor_idx),
        N = 1000000
    )
    #if params_[:output_params]["delete_temp"]
    #    rm(second_quant_folder,recursive=true)
    #end
    precursors_wide_arrow = writePrecursorCSV(
        joinpath(results_folder, "precursors_long.arrow"),
        sort(collect(values(file_id_to_parsed_name))),
        false,#normalized
        write_csv = params_[:output_params]["write_csv"],
        )
     #Summarize Precursor ID's
    println("Max LFQ...")
    if isfile( joinpath(results_folder, "protein_groups_long.arrow"))
        rm( joinpath(results_folder, "protein_groups_long.arrow"), force=true)
    end
    LFQ(
        DataFrame(Arrow.Table(joinpath(results_folder, "precursors_long.arrow"))),
        joinpath(results_folder, "protein_groups_long.arrow"),
        :peak_area,
        file_id_to_parsed_name,
        0.01f0,
        pg_qval_interp,
        batch_size = 100000
    )
    protein_groups_wide_arrow = writeProteinGroupsCSV(
        joinpath(results_folder, "protein_groups_long.arrow"),
        precursors[:sequence],
        precursors[:isotopic_mods],
        precursors[:structural_mods],
        precursors[:prec_charge],
        sort(collect(values(file_id_to_parsed_name))),
        write_csv = params_[:output_params]["write_csv"],
        )

    println("QC Plots")
    if isfile(joinpath(qc_plot_folder, "QC_PLOTS.pdf"))
        rm(joinpath(qc_plot_folder, "QC_PLOTS.pdf"))
    end
    println("parsed_fnames $parsed_fnames")
    println("file_id_to_parsed_name $file_id_to_parsed_name")
    qcPlots(
        precursors_wide_arrow,
        protein_groups_wide_arrow,
        params_,
        precursors,
        parsed_fnames,
        qc_plot_folder,
        file_id_to_parsed_name,
        MS_TABLE_PATHS,
        irt_rt,
        frag_err_dist_dict
    )
    if params_[:output_params]["delete_temp"]
        rm(temp_folder,recursive=true)
    end
end
    return 10
end