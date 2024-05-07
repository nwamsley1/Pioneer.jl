function parameterTuningSearch(
    #Mandatory Args
    spectra::Arrow.Table, 
    frag_index::FragmentIndex{Float32},
    precursors::Arrow.Table,
    ion_list::LibraryFragmentLookup{Float32},
    iRT_to_RT_spline::Any,
    ms_file_idx::UInt32,
    min_max_ppm::Tuple{Float32, Float32},
    params::NamedTuple,
    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
    collect_fmatches::Bool,
    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
    ionTemplates::Vector{Vector{DetailedFrag{Float32}}},
    iso_splines::IsotopeSplineModel{Float32},
    scored_PSMs::Vector{Vector{S}},
    unscored_PSMs::Vector{Vector{Q}},
    spectral_scores::Vector{Vector{R}},
    precursor_weights::Vector{Vector{Float32}},
    precs::Vector{Counter{UInt32, UInt8}}) where {S<:ScoredPSM{Float32, Float16},
                                                            Q<:UnscoredPSM{Float32},
                                                            R<:SpectralScores{Float16}}

    #constant_spline =  approximate(x->0.0f0, BSplineBasis(BSplineOrder(1), range(0.0, 1.0)), MinimiseL2Error())
    N = 200
    t = collect(LinRange(0.0, 4*π, N))
    u = sin.(t) 
    u .+= randn(N)./50
    placeholder_spline = UniformSpline(u, t, 3, 3)
    
    err_dist = MassErrorModel(
                placeholder_spline, 
                placeholder_spline,
                last(min_max_ppm),
                last(min_max_ppm)
    )                                
    return parameterTuningSearch(
        spectra, 
        frag_index,
        precursors, 
        ion_list,
        iRT_to_RT_spline,
        ms_file_idx,
        err_dist,
        searchScan!,
        ionMatches,
        ionMisses,
        all_fmatches,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        scored_PSMs,
        unscored_PSMs,
        spectral_scores,
        precursor_weights,
        precs,
        collect_fmatches = collect_fmatches,
        expected_matches = params[:expected_matches],
        isotope_err_bounds = Tuple([Int64(x) for x in params[:isotope_err_bounds]]),
        min_index_search_score = UInt8(params[:presearch_params]["min_index_search_score"]),
        min_frag_count = Int64(params[:presearch_params]["min_frag_count"]),
        min_log2_matched_ratio = Float32(params[:presearch_params]["min_log2_matched_ratio"]),
        min_spectral_contrast = Float32(params[:presearch_params]["min_spectral_contrast"]),
        min_topn_of_m = Tuple([Int64(x) for x in params[:presearch_params]["min_topn_of_m"]]),
        min_max_ppm = min_max_ppm,
        #min_max_ppm = ( Float32(params[:presearch_params]["frag_tol_ppm"]),  
        #                Float32(params[:presearch_params]["frag_tol_ppm"])
        #             ),
        max_best_rank = Int64(params[:presearch_params]["max_best_rank"]),
        quadrupole_isolation_width = params[:quadrupole_isolation_width],
        sample_rate = Float64(params[:presearch_params]["sample_rate"]),
    )
end
function parameterTuningSearch(
                    #Mandatory Args
                    spectra::Arrow.Table, 
                    frag_index::Union{FragmentIndex{Float32}, Missing},
                    precursors::Union{Arrow.Table, Missing},
                    ion_list::Union{LibraryFragmentLookup{Float32}, Missing},
                    rt_to_irt_spline::Any,
                    ms_file_idx::UInt32,
                    mass_err_model::MassErrorModel,
                    searchScan!::Union{Function, Missing},
                    ionMatches::Vector{Vector{FragmentMatch{Float32}}},
                    ionMisses::Vector{Vector{FragmentMatch{Float32}}},
                    all_fmatches::Vector{Vector{FragmentMatch{Float32}}},
                    IDtoCOL::Vector{ArrayDict{UInt32, UInt16}},
                    ionTemplates::Vector{Vector{L}},
                    iso_splines::IsotopeSplineModel{Float32},
                    scored_PSMs::Vector{Vector{S}},
                    unscored_PSMs::Vector{Vector{Q}},
                    spectral_scores::Vector{Vector{R}},
                    precursor_weights::Vector{Vector{Float32}},
                    precs::Union{Missing, Vector{Counter{UInt32, UInt8}}};
                    #keyword args
                    collect_fmatches = false,
                    expected_matches::Int64 = 100000,
                    frag_ppm_err::Float32 = 0.0f0,

                    δ::Float32 = 10000f0,
                    λ::Float32 = 0f0,
                    max_iter_newton::Int64 = 100,
                    max_iter_bisection::Int64 = 100,
                    max_iter_outer::Int64 = 100,
                    accuracy_newton::Float32 = 100f0,
                    accuracy_bisection::Float32 = 100000f0,
                    max_diff::Float32 = 0.01f0,

                    isotope_dict::Union{UnorderedDictionary{UInt32, Vector{Isotope{Float32}}}, Missing} = missing,
                    isotope_err_bounds::Tuple{Int64, Int64} = (3, 1),
                    min_frag_count::Int64 = 1,
                    min_spectral_contrast::Float32  = 0f0,
                    min_log2_matched_ratio::Float32 = -Inf32,
                    min_index_search_score::UInt8 = zero(UInt8),
                    min_topn_of_m::Tuple{Int64, Int64} = (2, 3),
                    min_max_ppm::Tuple{Float32, Float32} = (-Inf, Inf),
                    filter_by_rank::Bool = false, 
                    filter_by_count::Bool = true,
                    max_best_rank::Int64 = one(Int64),
                    n_frag_isotopes::Int64 = 1,
                    quadrupole_isolation_width::Float64 = 8.5,
                    rt_index::Union{retentionTimeIndex{T, Float32}, Vector{Tuple{Union{U, Missing}, UInt32}}, Missing} = missing,
                    irt_tol::Float64 = Inf,
                    sample_rate::Float64 = Inf,
                    spec_order::Set{Int64} = Set(2)
                    ) where {T,U<:AbstractFloat, 
                                                            L<:LibraryIon{Float32}, 
                                                            S<:ScoredPSM{Float32, Float16},
                                                            Q<:UnscoredPSM{Float32},
                                                            R<:SpectralScores{Float16}}


    ########
    #Each thread needs to handle a similair number of peaks. 
    #For example if there are 10,000 scans and two threads, choose n so that
    #thread 1 handles (0, n) and thread 2 handls (n+1, 10,000) and both seriestype
    #of scans have an equal number of fragment peaks in the spectra
    thread_tasks, total_peaks = partitionScansToThreads(spectra[:masses],
                                                        spectra[:retentionTime],
                                                         spectra[:centerMass],
                                                         spectra[:msOrder],
                                                        Threads.nthreads(),
                                                        1)

    if ismissing(precs)
        precs = [missing for _ in range(1, Threads.nthreads())]
    end
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra[:msOrder]))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return searchFragmentIndex(
                                spectra,
                                last(thread_task), #getRange(thread_task),
                                frag_index,
                                scan_to_prec_idx,
                                rt_to_irt_spline,
                                mass_err_model,
                                searchScan!,
                                frag_ppm_err,
                                precs[thread_id],
                                isotope_err_bounds,
                                min_index_search_score,
                                min_max_ppm,#(5.0f0, 5.0f0),#min_max_ppm,
                                quadrupole_isolation_width,
                                irt_tol,
                                sample_rate,
                                spec_order
                            )
        end
    end
    precursors_passed_scoring = fetch.(tasks)
    #return precursors_passed_scoring, scan_to_prec_idx
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return getPSMS(
                                spectra,
                                last(thread_task), #getRange(thread_task),
                                precursors,
                                scan_to_prec_idx,
                                precursors_passed_scoring[thread_id],
                                ion_list, 
                                rt_to_irt_spline,
                                ms_file_idx,
                                mass_err_model,
                                collect_fmatches,
                                frag_ppm_err,
                                δ,
                                λ,
                                max_iter_newton,
                                max_iter_bisection,
                                max_iter_outer,
                                accuracy_newton,
                                accuracy_bisection,
                                max_diff,
                                ionMatches[thread_id],
                                ionMisses[thread_id],
                                all_fmatches[thread_id],
                                IDtoCOL[thread_id],
                                ionTemplates[thread_id],
                                scored_PSMs[thread_id],
                                unscored_PSMs[thread_id],
                                spectral_scores[thread_id],
                                precursor_weights[thread_id],
                                precs[thread_id],
                                isotope_dict,
                                isotope_err_bounds,
                                min_frag_count,
                                min_spectral_contrast,
                                min_log2_matched_ratio,
                                min_topn_of_m,
                                min_max_ppm,
                                filter_by_rank,
                                filter_by_count,
                                max_best_rank,
                                n_frag_isotopes,
                                quadrupole_isolation_width,
                                irt_tol,
                                spec_order
                            )
        end
    end
    psms = fetch.(tasks)
    return psms
end

#Clear directory of previous QC Plots
[rm(joinpath(rt_alignment_folder, x)) for x in readdir(rt_alignment_folder)]
[rm(joinpath(mass_err_estimation_folder, x)) for x in readdir(mass_err_estimation_folder)]

test_time = @time begin
RT_to_iRT_map_dict = Dict{Int64, Any}()
frag_err_dist_dict = Dict{Int64,MassErrorModel}()
irt_errs = Dict{Int64, Float64}()
#MS_TABLE_PATH = "/Users/n.t.wamsley/TEST_DATA/PXD046444/arrow/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E30H50Y20_30SPD_DIA_4.arrow"
#MS_TABLE_PATH = "/Users/n.t.wamsley/TEST_DATA/PXD046444/arrow/test/20230324_OLEP08_200ng_30min_E20H50Y30_180K_2Th3p5ms_02.arrow"
for (ms_file_idx, MS_TABLE_PATH) in ProgressBar(collect(enumerate(MS_TABLE_PATHS)))
    MS_TABLE = Arrow.Table(MS_TABLE_PATH)
    out_fname = String(first(split(splitpath(MS_TABLE_PATH)[end],".")));
    #Randomly sample spectra to search and retain only the 
    #most probable psms as specified in "first_seach_params"
    start_ppm = 20.0f0
    max_ppm = 30.0f0
    params_[:presearch_params]["sample_rate"] = 0.02
    data_points = 0
    n = 0
    rtPSMs = nothing
    while n < 10
        RESULT =  parameterTuningSearch(
                                                MS_TABLE,
                                                prosit_lib["f_index"],
                                                prosit_lib["precursors"],
                                                library_fragment_lookup_table,
                                                x->x, #RT to iRT map'
                                                UInt32(ms_file_idx), #MS_FILE_IDX
                                                (start_ppm, start_ppm),
                                                params_,
                                                ionMatches,
                                                ionMisses,
                                                false,
                                                all_fmatches,
                                                IDtoCOL,
                                                ionTemplates,
                                                iso_splines,
                                                scored_PSMs,
                                                unscored_PSMs,
                                                spectral_scores,
                                                precursor_weights,
                                                precs,
                                                );

        rtPSMs = vcat([result for result in RESULT]...)
        #all_matches = vcat([last(result) for result in RESULT]...)
        #@time begin
        addPreSearchColumns!(rtPSMs, 
                                    MS_TABLE, 
                                    prosit_lib["precursors"][:is_decoy],
                                    prosit_lib["precursors"][:irt],
                                    prosit_lib["precursors"][:prec_charge],
                                    min_prob = params_[:presearch_params]["min_prob"])
        mean(rtPSMs[!,:error]./(rtPSMs[!,:b_count] .+ rtPSMs[!,:y_count]))
        if size(rtPSMs, 1) > 2000
            println("a ", size(rtPSMs, 1))
            break
        else# size(rtPSMs, 1) > 500
            println("b ", size(rtPSMs, 1))
            test = Float32(quantile(rtPSMs[!,:error]./(rtPSMs[!,:b_count] .+ rtPSMs[!,:y_count]), 0.99))
            #start_ppm = min(max_ppm, start_ppm*Float32(2))#Float32(quantile(rtPSMs[!,:error]./(rtPSMs[!,:b_count] .+ rtPSMs[!,:y_count]), 0.9))
            start_ppm = min(max_ppm, test*2)
            params_[:presearch_params]["sample_rate"] *= 1.5
        end
            #break
        #else
        #    println("c ", size(rtPSMs, 1))
        #    start_ppm *= min(max_ppm, start_ppm*Float32(2))
        #    #params_[:presearch_params]["sample_rate"] *= 1.5
        #end
        println("start_ppm $start_ppm")
        println("sample rate ", params_[:presearch_params]["sample_rate"])
        n += 1
    end

    RT_to_iRT_map = KDEmapping(rtPSMs[1:end,:RT], 
    rtPSMs[1:end,:iRT_predicted], 
    n = params_[:irt_mapping_params]["n_bins"], 
    bandwidth = params_[:irt_mapping_params]["bandwidth"]);

    plotRTAlign(rtPSMs[:,:RT], 
                rtPSMs[:,:iRT_predicted], 
                RT_to_iRT_map, 
                out_fdir = rt_alignment_folder,
                out_fname = out_fname
                );
    rtPSMs[!,:iRT_observed] = RT_to_iRT_map.(rtPSMs[!,:RT])
    irt_MAD = mad(rtPSMs[!,:iRT_observed] .- rtPSMs[!,:iRT_predicted])

    irt_σ = irt_MAD #Robust estimate of standard deviation

    irt_errs[ms_file_idx] = params_[:irt_err_sigma]*irt_σ
    RT_to_iRT_map_dict[ms_file_idx] = RT_to_iRT_map
    #=
    rtPSMs[!,:target] = [precursors[:is_decoy][x] for x in rtPSMs[!,:precursor_idx]]
    bins = LinRange(0, 2, 100)
    histogram(rtPSMs[rtPSMs[!,:target], :entropy_score], alpha = 0.5, bins = bins, normalize = :pdf)
    histogram!(rtPSMs[rtPSMs[!,:target].==false, :entropy_score], alpha = 0.5, bins = bins, normalize=:pdf)
    plot(rtPSMs[!,:RT],  rtPSMs[!,:iRT_predicted], seriestype=:scatter, alpha = 0.1)
     test_spline = UniformSpline( rtPSMs[!,:iRT_predicted], rtPSMs[!,:RT], 3, 5);
     plot!(LinRange(0, 45, 100),  test_spline.(LinRange(0, 45, 100)))
    =#
    function _getPPM(a::T, b::T) where {T<:AbstractFloat}
        (a-b)/(a/1e6)
    end
    #Get Retention Times and Target/Decoy Status 
    ####################
    #Use best_psms to estimate 
    #1) RT to iRT curve and 
    #2) mass error (ppm) distribution 

    RESULT =  parameterTuningSearch(
        MS_TABLE,
        prosit_lib["f_index"],
        prosit_lib["precursors"],
        library_fragment_lookup_table,
        RT_to_iRT_map, #RT to iRT map'
        UInt32(ms_file_idx), #MS_FILE_IDX
        (20.0f0, 20.0f0),#(max_ppm, max_ppm),
        params_,
        ionMatches,
        ionMisses,
        true,
        all_fmatches,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        scored_PSMs,
        unscored_PSMs,
        spectral_scores,
        precursor_weights,
        precs,
        );

    rtPSMs = vcat([first(result) for result in RESULT]...)
    all_matches = vcat([last(result) for result in RESULT]...)
    Int64(maximum([x.predicted_rank for x in all_matches]))
    addPreSearchColumns!(rtPSMs, 
                        MS_TABLE, 
                        prosit_lib["precursors"][:is_decoy],
                        prosit_lib["precursors"][:irt],
                        prosit_lib["precursors"][:prec_charge],
                        min_prob = 0.99)#params_[:presearch_params]["min_prob"])

    best_precursors = Set(rtPSMs[:,:precursor_idx]);
    best_matches = [match for match in all_matches if match.prec_id ∈ best_precursors];
    frag_ppm_errs = [_getPPM(match.theoretical_mz, match.match_mz) for match in best_matches];
    #frag_ppm_errs = [Float64(match.theoretical_mz - match.match_mz) for match in best_matches];
    frag_ppm_intensities = [match.intensity for match in best_matches];
    #frag_ppm_intensities = [match.theoretical_mz for match in best_matches];
    mass_err_model = ModelMassErrs(
        frag_ppm_intensities,
        frag_ppm_errs,
        Float64(max_ppm),
        max_n_bins = 10,
        min_bin_size = 400,
        frag_err_quantile = 0.99,
        out_fdir = mass_err_estimation_folder,
        out_fname = out_fname
    )
#=
    err_df = mass_err_model
    bins = LinRange(-30, 30, 50)
    histogram( groupby(  err_df,:bins)[end][!,:ppm_errs], alpha = 0.5, bins = bins, normalize=:pdf)
    histogram!( groupby(  err_df,:bins)[end-1][!,:ppm_errs], alpha = 0.5, bins = bins, normalize=:pdf)
    histogram!( groupby(  err_df,:bins)[end-2][!,:ppm_errs], alpha = 0.5, bins = bins, normalize=:pdf)
    histogram!( err_df[!,:ppm_errs], alpha = 0.5, bins = bins, normalize=:pdf)
    quantile(frag_ppm_errs .- median(frag_ppm_errs), 0.99)
    quantile(frag_ppm_errs .- median(frag_ppm_errs), 0.01)
=#

    #=

    ities = [match.theoretical_mz for match in best_matches];
    test_intensities, test_shapes = ModelMassErrs(
        frag_ppm_intensities,
        frag_ppm_errs,
        Float64(max_ppm),
        max_n_bins = 30,
        min_bin_size = 300,
        frag_err_quantile = 0.99,
        out_fdir = mass_err_estimation_folder,
        out_fname = out_fname
    )
    test_spline = BSplineApprox(test_shapes, 
                                test_intensities, 3, 4, :ArcLen, :Uniform, extrapolate = true)
    eval_points = LinRange(minimum(test_intensities)/1.2, maximum(test_intensities)*1.2, 100)
    plot(test_intensities, test_shapes, seriestype=:scatter)
    plot!(eval_points, test_spline.(eval_points))


        test_spline = BSplineApprox(test_shapes, 
                                test_intensities, 3, 4, :ArcLen, :Uniform, extrapolate = true)
        @btime test_spline(10.0f0)   
                test_spline = BSplineApprox(test_shapes, 
                                test_intensities, 3, 4, :Uniform, :Uniform, extrapolate = true)
        @btime test_spline(10.0f0)          
                        test_spline = BSplineApprox(test_shapes, 
                                test_intensities, 3, 4, :Uniform, :Average, extrapolate = true)
        @btime test_spline(10.0f0)              
        ppm_diffs = Int64[]
        for i in range(200000, 220000) 
            if MS_TABLE[:msOrder][i] == 2
                push!(ppm_diffs, sum(diff(MS_TABLE[:masses][i])./(MS_TABLE[:masses][i][1:end - 1]./1e6).<40.0))
            end
        end
            err_df = ModelMassErrs(
            frag_ppm_intensities,
            frag_ppm_errs,
            params_[:presearch_params]["frag_tol_ppm"],
            n_intensity_bins = length(frag_ppm_errs)÷250,#Int64(params_[:presearch_params]["samples_per_mass_err_bin"]),
            frag_err_quantile = params_[:frag_tol_params]["frag_tol_quantile"],
            out_fdir = mass_err_estimation_folder,
            out_fname = out_fname
        )
        bins = LinRange(-12, 12, 50)
        histogram(groupby(  err_df,:bins)[end-2][!,:ppm_errs], alpha = 0.5, bins = bins)
        histogram!(groupby(  err_df,:bins)[end-1][!,:ppm_errs], alpha = 0.5, bins = bins)
        histogram!(groupby(  err_df,:bins)[end-15][!,:ppm_errs], alpha = 0.5, bins = bins)
        #histogram!(groupby(  err_df,:bins)[end - 2][!,:ppm_errs], alpha = 0.5, bins = bins)

        plot(frag_ppm_intensities, frag_ppm_errs, alpha = 0.1, seriestype=:scatter)


        plot(log2.(frag_ppm_intensities), frag_ppm_errs, alpha = 0.02, seriestype=:scatter)
        test_loess = loess(log2.(frag_ppm_intensities), frag_ppm_errs, span = 0.5)
        plot!(LinRange(8, 21, 100), Loess.predict(test_loess, LinRange(8, 21, 100)))
        intensities = Float64[]
        offsets = Float64[]
        for (int_bin, subdf) in pairs(groupby(err_df, :bins))
            push!(intensities, median(2 .^subdf[!,:log2_intensities]))
            push!(offsets, (median((subdf[!,:ppm_errs]))))
        end


        plot(intensities, offsets, seriestype=:scatter)
        #plot(intensities, frag_ppm_errs, alpha = 0.02, seriestype=:scatter)
        test_loess = loess((intensities), offsets, span = 0.5)
        loess_range = LinRange(minimum(intensities), maximum(intensities), 100)
        plot!(loess_range, 
                Loess.predict(test_loess, 
                                loess_range)
            )


        test_interp = LinearInterpolation(intensities, offsets,extrapolation_bc=Line()) 
        plot(intensities, offsets, seriestype=:scatter)
        poly = Polynomials.fit(intensities, offsets, 10) 
        plot!(collect(LinRange(0, 2e5, 100)), poly.(LinRange(0, 2e5, 100)))
        #plot!(LinRange(0, 2e5, 100), test_interp.(LinRange(0, 2e5, 100)))

        bins = LinRange(-12, 12, 100)
        histogram( groupby(  err_df,:bins)[end][!,:ppm_errs], alpha = 0.5, bins = bins)
        histogram!( groupby(  err_df,:bins)[end-1][!,:ppm_errs], alpha = 0.5, bins = bins)
        histogram!( groupby(  err_df,:bins)[end-2][!,:ppm_errs], alpha = 0.5, bins = bins)
        histogram!( err_df[!,:ppm_errs], alpha = 0.5, bins = bins)

        histogram!(groupby(  err_df,:bins)[end-1][!,:ppm_errs], alpha = 0.5, bins = bins)


    =#
    PLOT_PATH = joinpath(MS_DATA_DIR, "Search", "QC_PLOTS", split(splitpath(MS_TABLE_PATH)[end],".")[1])
    #File name but remove file type
    frag_err_dist_dict[ms_file_idx] = mass_err_model
end
end

#Merge Quality Control PDFs 
merge_pdfs([joinpath(rt_alignment_folder, x) for x in readdir(rt_alignment_folder) if endswith(x, ".pdf")], 
                joinpath(rt_alignment_folder, "rt_alignment_plots.pdf"), cleanup=true)
merge_pdfs([joinpath(mass_err_estimation_folder, x) for x in readdir(mass_err_estimation_folder) if endswith(x, ".pdf")], 
    joinpath(mass_err_estimation_folder, "mass_err_estimation_folder.pdf"), cleanup=true)
