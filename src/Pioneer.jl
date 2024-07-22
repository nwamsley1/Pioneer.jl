module Pioneer

    using Arrow, ArgParse
    using BSplineKit, Base64
    using Base.Order
    using Base.Iterators: partition
    using CSV, CategoricalArrays, Combinatorics, CodecZlib
    using DataFrames, DataStructures, Dictionaries, Distributions 
    using ExpectationMaximization
    using FASTX
    using Interpolations
    using JSON, JLD2
    using LinearAlgebra, LoopVectorization, LinearSolve, LightXML
    using Measures
    using NumericalIntegration
    using Plots, PrettyPrinting, Polynomials, PDFmerger, ProgressBars
    using Tables
    using StatsPlots
    using Random
    using StaticArrays, StatsBase, SpecialFunctions, Statistics
    using XGBoost
    #create_app("../Pioneer","../Pioneer_Compiled", force = true)


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
    [include(joinpath("Utils", jl_file)) for jl_file in [

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

    [include(joinpath("PSM_TYPES", jl_file)) for jl_file in ["PSM.jl","spectralDistanceMetrics.jl","UnscoredPSMs.jl","ScoredPSMs.jl"]];

    #Files needed for PRM routines
    [include(joinpath("Routines","LibrarySearch","methods",jl_file)) for jl_file in [
                                                                                    "matchPeaks.jl",
                                                                                    "buildDesignMatrix.jl",
                                                                                    "manipulateDataFrames.jl",
                                                                                    "buildRTIndex.jl",
                                                                                    "searchRAW.jl",
                                                                                    "selectTransitions.jl",
                                                                                    "integrateChroms.jl",
                                                                                    "queryFragmentIndex.jl"]];
    include("Routines/LibrarySearch/MAIN.jl")     
    include("Routines/LibrarySearch/parameterTuningSearch.jl")
    include("Routines/LibrarySearch/firstSearch.jl")
    include("Routines/LibrarySearch/quantitativeSearch.jl")
    include("Routines/LibrarySearch/scoreTraces.jl")
    include("Routines/LibrarySearch/secondQuant.jl")
    include("Routines/LibrarySearch/proteinQuant.jl")
    include("Routines/LibrarySearch/qcPlots.jl")
    const methods_path = joinpath(pwd(), "src","Routines","LibrarySearch")                  
    function julia_main()::Cint
        ARGS = Dict(
                "params_json"=>"data/example_config/LibrarySearch.json"
                )

        params = JSON.parse(read(ARGS["params_json"], String));
        params = JSON.parse(read("data/example_config/LibrarySearch.json", String));
        MS_DATA_DIR = params["ms_data_dir"];
        SPEC_LIB_DIR = params["library_folder"];

        #Get all files ending in ".arrow" that are in the MS_DATA_DIR folder. 
        MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
 
        params_ = parseParams(params)


        OUT_DIR = joinpath(params_[:benchmark_params]["results_folder"], "RESULTS");
        if !isdir(out_folder)
            mkpath(out_folder)
        end

        search_folder = joinpath(OUT_DIR, "Search")
        if !isdir(search_folder)
            mkpath(search_folder)
        end
        qc_plot_folder = joinpath(search_folder, "QC_PLOTS")
        if !isdir(qc_plot_folder)
            mkpath(qc_plot_folder)
        else
            [rm(joinpath(qc_plot_folder, x)) for x in readdir(qc_plot_folder) if endswith(x, ".pdf")]
        end
        rt_alignment_folder = joinpath(search_folder, "QC_PLOTS","rt_alignment")
        if !isdir(rt_alignment_folder)
            mkpath(rt_alignment_folder)
        end
        mass_err_estimation_folder = joinpath(search_folder, "QC_PLOTS","mass_error_estimation")
        if !isdir(mass_err_estimation_folder)
            mkpath(mass_err_estimation_folder)
        end
        results_folder = joinpath(results_folder, "RESULTS")
        if !isdir(results_folder)
            mkpath(results_folder)
        end
        params_folder = joinpath(params_folder, "PARAMS")
        if !isdir(params_folder )
            mkpath(params_folder)
        end
        #Write params to folder 
        open(joinpath(params_folder, "config.json"),"w") do f 
            JSON.print(f, params)
        end

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
        library_fragment_lookup_table = LibraryFragmentLookup(detailed_frags, prec_frag_ranges)
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
        pid_to_cv_fold = getCVFolds(
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
        file_id_to_parsed_name, parsed_fnames = parseFileNames(MS_TABLE_PATHS)

        ###########
        #Tune Parameters
        println("Parameter Tuning Search...")
        RT_to_iRT_map_dict, frag_err_dist_dict, irt_errs = parameterTuningSearch(rt_alignment_folder,
                                                                                mass_err_estimation_folder,
                                                                                MS_TABLE_PATHS,
                                                                                params_,
                                                                                prosit_lib,
                                                                                library_fragment_lookup_table,
                                                                                ionMatches,
                                                                                ionMisses,
                                                                                all_fmatches,
                                                                                IDtoCOL,
                                                                                ionTemplates,
                                                                                iso_splines,
                                                                                scored_PSMs,
                                                                                unscored_PSMs,
                                                                                spectral_scores,
                                                                                precs)
        println("Parameter Tuning Search...")
        PSMs_Dict = firstSearch(
            RT_to_iRT_map_dict,
            frag_err_dist_dict,
            irt_errs,
            file_id_to_parsed_name,
            MS_TABLE_PATHS,
            params_,
            prosit_lib,
            library_fragment_lookup_table,
            ionMatches,
            ionMisses,
            all_fmatches,
            IDtoCOL,
            ionTemplates,
            iso_splines,
            scored_PSMs,
            unscored_PSMs,
            spectral_scores,
            precs
        )

        ##########
        #Combine First Search Results
        ##########
        println("Combining First Search Results...")
        iRT_RT, RT_iRT = mapRTandiRT(PSMs_Dict,
        min_prob = Float64(params_[:irt_mapping_params]["min_prob"])
        )
        precID_to_iRT = getPrecIDtoiRT(PSMs_Dict, 
                RT_iRT, 
                max_precursors = Int64(params_[:summarize_first_search_params]["max_precursors"]))
        irt_err = getRTErr(PSMs_Dict, RT_iRT)
        #irt_err = 0.4
        RT_INDICES = makeRTIndices(PSMs_Dict,
            precID_to_iRT,
            RT_iRT,
            bin_rt_size = 0.1,
            min_prob = params_[:summarize_first_search_params]["max_prob_to_impute"])
        
        ############
        #Quantitative Search
        println("Begining Quantitative Search...")
        best_psms = vcat(values(quantSearch(
            frag_err_dist_dict,
            pid_to_cv_fold,
            precID_to_iRT,
            RT_INDICES,
            RT_iRT,
            irt_err,
            chromatograms,
            file_id_to_parsed_name,
            MS_TABLE_PATHS,
            params_,
            precursors,
            prosit_lib,
            library_fragment_lookup_table,
            ionMatches,
            ionMisses,
            IDtoCOL,
            ionTemplates,
            iso_splines,
            complex_scored_PSMs,
            complex_unscored_PSMs,
            complex_spectral_scores,
            precursor_weights
            ))...)
        println("Traning Target-Decoy Model...")
        best_psms = scoreTraces!(best_psms, precursors)

        #Get protein names 
        transform!(best_psms, AsTable(:) => ByRow(psm -> 
        prosit_lib["precursors"][:accession_numbers][psm[:precursor_idx]]
        ) => :accession_numbers
        );


        traces_passing = Set(best_psms[(best_psms[!,:q_value].<=0.01).&(best_psms[!,:target]),:precursor_idx])

        ###########
        #Score Protein Groups
        scored_proteins = scoreProteinGroups!(best_psms)
        protein_to_q_value = Dict{Tuple{UInt32, String}, Float32}()
        for i in range(1, size(scored_proteins, 1))
            protein_to_q_value[
                (
                    UInt32(scored_proteins[i,:ms_file_idx]),
                    scored_proteins[i,:accession_numbers]
                )
            ] = scored_proteins[i,:q_value]
        end

        ###########
        #Re-quantify with 1% fdr precursors 
        best_psms[!,:peak_area] = zeros(Float32, size(best_psms, 1))
        best_psms[!,:new_best_scan] = zeros(UInt32, size(best_psms, 1))
        grouped_best_psms = groupby(best_psms, :file_name)
        println("Re-Quantifying Precursors at FDR Threshold...")
        secondQuant!( 
                        grouped_best_psms,   
                        frag_err_dist_dict,
                        traces_passing,
                        RT_INDICES,
                        RT_iRT,
                        irt_err,
                        chromatograms,
                        file_id_to_parsed_name,
                        MS_TABLE_PATHS,
                        params_,
                        precursors,
                        prosit_lib,
                        library_fragment_lookup_table,
                        ionMatches,
                        ionMisses,
                        IDtoCOL,
                        ionTemplates,
                        iso_splines,
                        complex_scored_PSMs,
                        complex_unscored_PSMs,
                        complex_spectral_scores,
                        precursor_weights)

        ###########
        #Ungroup precursors and get the best trace for each 
        best_psms = DataFrame(grouped_best_psms)
        #Must be based on weight and not peak_area because peak_area was not calculated for decoys
        getBestTrace!(best_psms, 0.01, :weight)
        filter!(x->x.best_trace, best_psms)
        #precursor level q_value 
        getQvalues!(best_psms[!,:prob], best_psms[:,:target], best_psms[!,:q_value]);
        #filter out decoys and low-scoring targets
        filter!(x->(x.q_value<=0.01)&(x.target)&(!isnan(x.peak_area)), best_psms)
        #IDs_PER_FILE = value_counts(best_psms, [:file_name])
        best_psms[!,:species] = [precursors[:proteome_identifiers][pid] for pid in best_psms[!,:precursor_idx]]
        ###########
        #Normalize Quant 
        ###########
        println("Normalization...")
        normalizeQuant(
                best_psms,
                :peak_area,
                N = params_[:normalization_params]["n_rt_bins"],
                spline_n_knots = params_[:normalization_params]["spline_n_knots"],
                max_q_value = params_[:normalization_params]["max_q_value"],
                min_points_above_FWHM = params_[:normalization_params]["min_points_above_FWHM"]
            )

        ############
        #Prep for Protein Inference 
        best_psms[!,:species] = [precursors[:proteome_identifiers][pid] for pid in best_psms[!,:precursor_idx]]
        best_psms[!,:ms_file_idx] =  UInt32.(best_psms[!,:ms_file_idx])
        best_psms[!,:peak_area] =  allowmissing(best_psms[!,:peak_area])
        best_psms[!,:peak_area_normalized] =  allowmissing(best_psms[!,:peak_area_normalized])

        gbpsms = groupby(best_psms,:ms_file_idx)
        file_idx_dicts = [Dict{UInt32, @NamedTuple{prob::Float32, qvalue::Float32, peak_area::Float32}}() for _ in range(1, length(gbpsms))]
        for (file_idx, bpsms) in pairs(gbpsms)
            for i in range(1, size(bpsms, 1))
                    file_idx_dicts[file_idx[:ms_file_idx]][bpsms[i,:precursor_idx]] = (
                    prob = bpsms[i,:prob],
                    qvalue = bpsms[i,:q_value],
                    peak_area = bpsms[i,:peak_area])
            end
        end
        best_psms[!,:structural_mods] = [precursors[:structural_mods][pid] for pid in best_psms[!,:precursor_idx]]
        best_psms[!,:isotopic_mods] = [precursors[:isotopic_mods][pid] for pid in best_psms[!,:precursor_idx]]

        features = [:species,:accession_numbers,:sequence,:structural_mods,
                    :isotopic_mods,:charge,:precursor_idx,:target,:weight,:peak_area,:peak_area_normalized,
                    :scan_idx,:prob,:q_value,:prec_mz,:RT,:irt_obs,:irt_pred,:best_rank,:best_rank_iso,:topn,:topn_iso,:longest_y,:longest_b,:b_count,
                    :y_count,:p_count,:non_cannonical_count,:isotope_count
                    ,:matched_ratio]
        sort!(best_psms,[:species,:accession_numbers,:sequence,:target])
        Arrow.write(joinpath(results_folder,"best_psms.arrow"),best_psms[!,features])
        wide_psms_quant = unstack(best_psms,[:species,:accession_numbers,:sequence,:structural_mods,:isotopic_mods,:precursor_idx,:target],:file_name,:peak_area_normalized)
        sort!(wide_psms_quant,[:species,:accession_numbers,:sequence,:target])
        CSV.write(joinpath(results_folder,"best_psms_wide.arrow"),wide_psms_quant)

        #Summarize Precursor ID's
        value_counts(df, col) = combine(groupby(df, col), nrow)

        precursor_id_table = value_counts(best_psms,:file_name)
        #CSV.write("/Users/n.t.wamsley/Desktop/precursor_ids_table.csv")
        println("Max LFQ...")
        protein_quant = proteinQuant(
            best_psms,
            protein_to_q_value
        )

        protein_quant[!,:experiments] = UInt32.(protein_quant[!,:experiments])
        protein_quant[!,:file_name] = [file_id_to_parsed_name[ms_file_idx] for ms_file_idx in protein_quant[!,:experiments]]

        #Wide DataFormat 
        wide_protein_quant = unstack(protein_quant,[:species,:protein,:target],:experiments,:log2_abundance)
        sort!(wide_protein_quant,[:species,:protein,:target])
        CSV.write(joinpath(results_folder,"proteins_wide.csv"),wide_protein_quant)

        qcPlots(
            best_psms,
            protein_quant,
            params_,
            parsed_fnames,
            qc_plot_folder,
            file_id_to_parsed_name,
            MS_TABLE_PATHS,
            iRT_RT,
            frag_err_dist_dict
        )
        return 0
    end

end

#=
using Revise
using Pkg
@timed begin
Pkg.activate(".")
using Pioneer
Pioneer.julia_main()
end
println("urmom")
=#