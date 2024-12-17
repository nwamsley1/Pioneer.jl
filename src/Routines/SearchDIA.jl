function SearchDIA(params_path::String)
    @info "Loading Parameters..."
    if !isabspath(params_path)
        params_path = joinpath(@__DIR__, "../../", params_path)
    end
    params = JSON.parse(read(params_path, String));
    MS_DATA_DIR = params["ms_data_dir"];
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
    @info "Loading Spetral Library..."
    spec_lib = loadSpectralLibrary(SPEC_LIB_DIR, params_);
    SPEC_LIB = FragmentIndexLibrary(
        spec_lib["presearch_f_index"], 
        spec_lib["f_index"], 
        BasicLibraryPrecursors(spec_lib["precursors"]), 
        spec_lib["f_det"])

    @info "Initializing Search Context..."
    SEARCH_CONTEXT = initSearchContext(
        SPEC_LIB,
        parseIsoXML(joinpath(@__DIR__,"../data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml")),
        ArrowTableReference(MS_TABLE_PATHS),
        Threads.nthreads(),
        250000  # buffer size
    );

    setDataOutDir!(SEARCH_CONTEXT, params_[:benchmark_params]["results_folder"]);

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
        @time execute_search(
            QuadTuningSearch(), SEARCH_CONTEXT, params_
        );
    
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
        @time execute_search(
            HuberTuningSearch(), SEARCH_CONTEXT, params_
        );

        #==========================================================
        Estimate huber (Charboneir) loss smoothing parameter δ. Increasing δ
            behaves like squared error and decreasing δ behaves like absolute error. 
            Used for linear regression of MS/MS spectra onto library spectra. 
            Results stored in SEARCH_CONTEXT.huber_delta::Base.Ref{Float32}
        ==========================================================#
        @time execute_search(
            SecondPassSearch(), SEARCH_CONTEXT, params_
        );

        @time execute_search(
            ScoringSearch(), SEARCH_CONTEXT, params_
        );

        @time execute_search(
            IntegrateChromatogramSearch(), SEARCH_CONTEXT, params_
        );

        @time execute_search(
            MaxLFQSearch(), SEARCH_CONTEXT, params_
        );
    end
end