function SearchDIA(params_path::String)
    #println("JLD2 version is: ", Pkg.installed()["JLD2"])
    total_time = @timed begin
    #params_path = "/Users/n.t.wamsley/RIS_temp/PIONEER_PAPER/SPEC_LIBS/Altimeter111124_MixedSpecies_OlsenAstral_NoEntrapment_SplineTest.poin/TestSplineLibrary.json"
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

    qc_plot_folder, rt_alignment_folder, mass_err_estimation_folder, results_folder, temp_folder = makeOutputDirectories(
        joinpath(params_[:benchmark_params]["results_folder"], "RESULTS"),
        params
    )

    first_search_psms_folder = joinpath(temp_folder, "first_search_psms")
    if !isdir(first_search_psms_folder)
        mkpath(first_search_psms_folder)
    end

    irt_indices_folder = joinpath(temp_folder, "irt_indices_folder")
    if !isdir(irt_indices_folder)
        mkpath(irt_indices_folder)
    end

    quant_psms_folder = joinpath(temp_folder, "quant_psms_folder")
    if !isdir(quant_psms_folder )
        mkpath(quant_psms_folder )
    end

    passing_psms_folder = joinpath(temp_folder, "passing_psms")
    if !isdir( passing_psms_folder )
    mkpath( passing_psms_folder )
    end

    passing_proteins_folder = joinpath(temp_folder, "passing_proteins")
    if !isdir( passing_proteins_folder )
       mkpath( passing_proteins_folder )
    end

    second_quant_folder = joinpath(temp_folder, "second_quant")
    if !isdir( second_quant_folder)
        mkpath( second_quant_folder)
    end
    ###########
    #Load Spectral Libraries
    ###########
    spec_lib = loadSpectralLibrary(SPEC_LIB_DIR);
    precursors = spec_lib["precursors"];
    unique_proteins = unique(precursors[:accession_numbers]);
    accession_number_to_id = Dict(zip(unique_proteins, range(one(UInt32), UInt32(length(unique_proteins)))));
    ###########
    #Set CV Folds 
    ###########
    pid_to_cv_fold = getCVFolds(
        collect(range(UInt32(1), UInt32(length(spec_lib["precursors"][:sequence])))),#precursor id's, 
        spec_lib["precursors"][:accession_numbers]
        );

    N = Threads.nthreads()
    M = 250000
    n_precursors = length(spec_lib["precursors"][:mz])
    ionMatches = [[FragmentMatch{Float32}() for _ in range(1, M)] for _ in range(1, N)];
    ionMisses = [[FragmentMatch{Float32}() for _ in range(1, M)] for _ in range(1, N)];
    all_fmatches = [[FragmentMatch{Float32}() for _ in range(1, M)] for _ in range(1, N)];
    IDtoCOL = [ArrayDict(UInt32, UInt16, n_precursors ) for _ in range(1, N)];
    ionTemplates = [[DetailedFrag{Float32}() for _ in range(1, M)] for _ in range(1, N)];
    iso_splines = parseIsoXML(joinpath(@__DIR__,"../../data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml"));
    #iso_splines = parseIsoXML(joinpath(@__DIR__,"../data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml"));
    scored_PSMs = [Vector{SimpleScoredPSM{Float32, Float16}}(undef, 5000) for _ in range(1, N)];
    unscored_PSMs = [[SimpleUnscoredPSM{Float32}() for _ in range(1, 5000)] for _ in range(1, N)];
    spectral_scores = [Vector{SpectralScoresSimple{Float16}}(undef, 5000) for _ in range(1, N)];
    precursor_weights = [zeros(Float32, n_precursors ) for _ in range(1, N)];
    precs = [Counter(UInt32, UInt8,n_precursors ) for _ in range(1, N)];
    chromatograms = [Vector{ChromObject}(undef, 5000) for _ in range(1, N)];
    complex_scored_PSMs = [Vector{ComplexScoredPSM{Float32, Float16}}(undef, 5000) for _ in range(1, N)];
    complex_unscored_PSMs = [[ComplexUnscoredPSM{Float32}() for _ in range(1, 5000)] for _ in range(1, N)];
    complex_spectral_scores = [Vector{SpectralScoresComplex{Float16}}(undef, 5000) for _ in range(1, N)];

    file_id_to_parsed_name, parsed_fnames,file_path_to_parsed_name = parseFileNames(MS_TABLE_PATHS)



   # using Profile 
   # using PProf
   # @profile begin
        RT_to_iRT_map_dict, frag_err_dist_dict, irt_errs = parameterTuningSearch(rt_alignment_folder, #ms_file_idx_to_remove, failed_ms_file_idxs
        mass_err_estimation_folder,
        MS_TABLE_PATHS,
        params_,
        spec_lib,
        ionMatches,
        ionMisses,
        all_fmatches,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        scored_PSMs,
        unscored_PSMs,
        spectral_scores,
        precs);
    end
   # pprof()

   nce_model_dict = nceTuningSearch(
        RT_to_iRT_map_dict,
        frag_err_dist_dict,
        irt_errs,
        MS_TABLE_PATHS,
        params_,
        spec_lib,
        LinRange(21.0f0, 40.0f0, 15),
        ionMatches,
        ionMisses,
        all_fmatches,
        IDtoCOL,
        ionTemplates,
        iso_splines,
        scored_PSMs,
        unscored_PSMs,
        spectral_scores,
        precs);

    p = plot()
    for i in range(1, 18)
        model = nce_model_dict[i]
        plot!(p, LinRange(300, 1000, 100), model.(LinRange(300, 1000, 100), 3), label = string(i));
        #plot!(p, LinRange(300, 1000, 100), model.(LinRange(300, 1000, 100), 3), label = string(key));
        #break
    end
    p

    p = plot()
    for (key, model) in pairs(nce_model_dict)
        plot!(p, LinRange(300, 1000, 100), model.(LinRange(300, 1000, 100), 3), label = string(key));
    end
    p
    spsms = DataFrame(map(pairs(groupby(test_out,[:precursor_idx,:scan_idx]))) do (key, psms)
        max_arg = argmax(psms[!,:scribe])
        #=
        return (precursor_idx = psms[max_arg,:precursor_idx], 
                target = psms[max_arg,:target], 
                scribe = psms[max_arg,:scribe], 
                spectral_contrast = psms[max_arg,:spectral_contrast])
        =#
        return psms[max_arg,:]
    end)
    
    scorePresearch!(spsms)
    getQvalues!(spsms[!,:prob], spsms[!,:target], spsms[!,:q_value])
    filter!(x->(x.target)&(x.q_value<=0.01), spsms)
    passing_precs = Set(spsms[!,:precursor_idx])
    filter!(x->x.precursor_idx∈passing_precs, test_out)
    test_out_old = copy(test_out)
    #test_out = copy(test_out_old)
    
    n = 1
    sort!(test_out,:nce)
    gto = groupby(test_out,[:precursor_idx,:scan_idx]);
    #gto[(precursor_idx =  796635, scan_idx =  131676)]
    #gto[n][!,[:precursor_idx,:charge,:q_value,:nce,:spectral_contrast,:b_count,:y_count,:matched_ratio,:scribe,:scan_idx]]
    lib_path = "/Users/n.t.wamsley/Desktop/SplineLibTest/spline_examples"
    for n in range(1, 100)
        seq = precursors[:sequence][gto[n][1,:precursor_idx]]
        charge = precursors[:prec_charge][gto[n][1,:precursor_idx]]
        mz = precursors[:mz][gto[n][1,:precursor_idx]]
        p = plot(gto[n][!,:nce], gto[n][!,:spectral_contrast],
        xlabel = "Altimeter NCE",
        ylabel="Spectral Contrast",
        label="Spectral Contrast",
        legend = :bottomleft,
        title = "$seq^+$charge\n mz: $mz"
        )
        # Add second plot with right y-axis
        plot!(twinx(), gto[n][!,:nce], gto[n][!,:scribe],
            ylabel="Scribe",
            label="Scribe",
            legend = :bottomright,
            color=:red)  # Different color to distinguish the lines
        savefig(p, joinpath(lib_path, string(n)*".pdf"))
    end
    merge_pdfs([x for x in readdir(lib_path, join=true) if endswith(x, ".pdf")],
    joinpath(lib_path, "spline_examples.pdf"), cleanup=true)
    n += 1
    #precursor_idx =  0x001c0298
    test_out[!,:best_psms] .= false
    test_out[!,:prec_mz] = [precursors[:mz][pid] for pid in test_out[!,:precursor_idx]]
    gpsms = groupby(test_out,:precursor_idx)
    for (key, psms) in pairs(gpsms)
        #if maximum(psms[!,:spectral_contrast]) > 0.95
            psms[argmax(psms[!,:scribe]),:best_psms] = true
        #end
    end
    filter!(x->x.best_psms, test_out)

    p = plot(size = (300, 300))

    pwl = fit_nce_model(PiecewiseNceModel(0.0f0),
    test_out[!,:prec_mz], test_out[!,:nce], test_out[!,:charge], NCE_MODEL_BREAKPOINT)
    plot!(p, LinRange(300, 1000, 100), pwl.(LinRange(300, 1000, 100), 2));
    plot!(p, LinRange(300, 1000, 100), pwl.(LinRange(300, 1000, 100), 3));
    plot!(p, LinRange(300, 1000, 100), pwl.(LinRange(300, 1000, 100), 4));
    plot!(p, test_out2[!,:prec_mz], test_out2[!,:nce], seriestype=:scatter, alpha = 0.1, xlim = (390,1000), bins = 100);
  
    #test_out2[!,:spl_err] = abs.(test_out2[!,:nce] .- pwl_2.(test_out2[!,:prec_mz]))
    #spl_err = 2*mad(test_out2[!,:spl_err])
    #filter!(x->x.spl_err<spl_err, test_out2)
    #pwl_2 = fit_piecewise_linear(test_out2[!,:prec_mz], test_out2[!,:nce], 500.0f0);
    #plot!(p, LinRange(300, 1000, 100), pwl_2.(LinRange(300, 1000, 100)));
    plot!(p, test_out2[!,:prec_mz], test_out2[!,:nce], seriestype=:scatter, alpha = 0.1, xlim = (390,1000), bins = 100);
    p
    
    p = plot(size = (300, 300))
    test_out2 = test_out[test_out[!,:charge].==3,:];
    #plot!(p, test_out2[!,:prec_mz], test_out2[!,:nce], seriestype=:scatter, alpha = 0.1, xlim = (390,1000), bins = 100);
    pwl_2 = fit_piecewise_linear(test_out2[!,:prec_mz], test_out2[!,:nce], 500.0f0);
    plot!(p, LinRange(300, 1000, 100), pwl_2.(LinRange(300, 1000, 100)));
    test_out2[!,:spl_err] = abs.(test_out2[!,:nce] .- pwl_2.(test_out2[!,:prec_mz]))
    plot!(p, test_out2[!,:prec_mz], test_out2[!,:nce], seriestype=:scatter, alpha = 0.1, xlim = (390,1000), bins = 100);
    p



    test_out2 = test_out[test_out[!,:charge].==4,:]
    plot!(p, test_out2[!,:prec_mz], test_out2[!,:nce], seriestype=:scatter, alpha = 0.05, xlim = (390,1000), bins = 100)
    pwl_4 = fit_piecewise_linear(test_out2[!,:prec_mz], test_out2[!,:nce], 500.0f0)
    plot!(p, LinRange(300, 1000, 100), pwl_4.(LinRange(300, 1000, 100)))

    pwl_fits = (pwl_2, pwl_3,pwl_4)

    test_out2 = test_out[test_out[!,:charge].==4,:]
    #plot!(p, test_out2[!,:prec_mz], test_out2[!,:nce], seriestype=:scatter, alpha = 0.2, xlim = (390,1000))
    spl_test = UniformSpline(test_out2[!,:nce], test_out2[!,:prec_mz], 3, 3)
    #plot!(p, LinRange(300, 1000, 100), spl_test.(LinRange(300, 1000, 100)))
    test_out2[!,:spl_err] = abs.(test_out2[!,:nce] .- spl_test.(test_out2[!,:prec_mz]))
    spl_err = mad(test_out2[!,:spl_err])
    filter!(x->x.spl_err<spl_err, test_out2)
    #p = plot(size = (300, 300))
    plot!(p, test_out2[!,:prec_mz], test_out2[!,:nce], seriestype=:scatter, alpha = 0.2, xlim = (390,1000))
    spl_test = UniformSpline(test_out2[!,:nce], test_out2[!,:prec_mz], 3,5)
    plot!(p, LinRange(300, 1000, 100), spl_test.(LinRange(300, 1000, 100)))



    struct NceModel{N, T<:AbstractFloat}
        mz_bin::NTuple{N, T}
        slope::NTuple{N, T}
        intercept::NTuple{N, T}
    end

    struct NceModel{N,M,T<:AbstractFloat}
        mz_bin::NTuple{N, T}
        nce_val::NTuple{N, NTuple{M, T}}
    end

    
    function (nce::NceModel{N, T})(mz::T, charge::UInt8)where {N,T<:AbstractFloat}
        for i in range(1, N)
            a, b, c = nce.slope[1], nce.slope[1], nce.mz_bin[1]
            if mz < c
                return T(charge*a) + b
            end
        end
        return T(charge*a[end]) + b[end]
    end
    function getOptimalNCE()
        prec_charge_range = range(2, 4)
        maximum_prc
        bin_mzs = (500.0f0, 700.0f0, typemax(Float32))
        sort!(test_out,:charge)
        mz_bins = Float32[]
        slopes = Float32[]
        intercepts = Float32[]
        for i in range(1, length(bin_mzs))
            subdf = test_out[test_out[!,:prec_mz].<bin_mzs[i],:]
            charge_states = groupby(subdf, :charge)
            c2 = median(charge_states[(charge = 2,)][!,:nce])
            c3 = median(charge_states[(charge = 3,)][!,:nce])
            ((c3 - c2)/(3 - 2))  + c2
        end
    end


    p = plot()
    test_out[!,:prec_mz] = [precursors[:mz][pid] for pid in test_out[!,:precursor_idx]]
    for (key, psms) in pairs(groupby(test_out,:charge))
        plot!(psms[!,:nce], psms[!,:prec_mz], show = true, alpha = 0.1, seriestype=:scatter, label = string(key))
    end

    p = plot(size = (300, 300))
    test_out2 = test_out[test_out[!,:charge].==2,:]
    plot!(p, test_out2[!,:prec_mz], test_out2[!,:nce], seriestype=:scatter, alpha = 0.2, xlim = (390,1000))
    spl_test = UniformSpline(test_out2[!,:nce], test_out2[!,:prec_mz], 3, 3)
    plot!(p, LinRange(300, 1000, 100), spl_test.(LinRange(300, 1000, 100)))
# Create the linear model
fit = lm(@formula(nce ~ prec_mz), test_out2_float[test_out2[!,:prec_mz].<=550.0,:])
x_fit = Float64.(collect(390:550))
y_fit = GLM.predict(fit, DataFrame(prec_mz = x_fit))
plot!(x_fit, y_fit, color=:red, linewidth=2)
# Create the linear model
fit = lm(@formula(nce ~ prec_mz), test_out2_float[test_out2[!,:prec_mz].>550.0,:])
x_fit = Float64.(collect(550:900))
y_fit = GLM.predict(fit, DataFrame(prec_mz = x_fit))
plot!(x_fit, y_fit, color=:red, linewidth=2)



    # Option 1: Using histogram2d
plot(test_out2[!,:prec_mz], test_out2[!,:nce], 
seriestype=:histogram2d,
xlim=(390,1000),
bins=50,  # Adjust this to control resolution
colorbar=true)


# Option 2: Using density plot (smoother)
plot(test_out2[!,:prec_mz], test_out2[!,:nce], 
    seriestype=:density2d,
    xlim=(390,1000),
    colorbar=true)



# Option 3: hexbin plot for a hexagonal heatmap
plot(test_out2[!,:prec_mz], test_out2[!,:nce], 
    seriestype=:hexbin,
    xlim=(390,1000),
    bins=50,  # Adjust this to control resolution
    colorbar=true)
    
    using DataFrames
    using StatsPlots
    using Statistics

    # Create the boxplot
    p = @df test_out boxplot(
    :nce,                    # x-axis: NCE values
    :prec_mz,               # y-axis: precision m/z values
    group=:nce,             # group by NCE
    legend=false,           # no legend needed
    xlabel="NCE",
    ylabel="Precursor m/z",
    title="Precursor m/z Distribution by NCE"
    )

    # Add sample size annotations
    for (i, nce_val) in enumerate(sort(unique(test_out.nce)))
    # Count samples for this NCE value
    n = count(==(nce_val), test_out.nce)
    # Add annotation above the boxplot
    annotate!(i, maximum(filter(r -> r.nce == nce_val, test_out).prec_mz) + 50,
                text("n=$n", :center, 8))
    end

    display(p)


    nce_ = Dict{String, Any}()
    nce_median = Dict{String, Float32}()
    nce_psms_01fdr = Dict{String, Int64}()
    _nce_ = Float32[]
    _charge_ = UInt8[]
    _θ_ = Float32[]
    _n_ = Int64[]
    p = plot()
    filter!(x->!isnan(x.spectral_contrast), test_out)
    for (key, v) in pairs(groupby(copy(test_out),[:charge,:nce]))  
        push!(_nce_, key[:nce])
        push!(_charge_, key[:charge])
        push!(_θ_, mean(Float32.(v[!,:spectral_contrast])))
        push!(_n_, size(v, 1))
    end

    results = DataFrame(
        nce = _nce_,
        charge = _charge_,
        spectral_contrast = _θ_,
        n = _n_
    )
    filter!(x->x.n>200, results)
    sort!(results,:spectral_contrast)

    p = plot()
    for (k, v) in pairs(groupby(results, [:charge]))
        sort!(v, :nce)
        plot!(p, v[!,:nce], v[!,:spectral_contrast], label = string(k), show = true)
    end

    p = plot()
    for (i, (k, v)) in enumerate(pairs(groupby(results, [:charge])))
        sort!(v, :nce)
        plot!(p, v[!,:nce], v[!,:spectral_contrast], label = string(k), show = true, color = i)
        println(v[argmax(v[!,:spectral_contrast]), :nce])
        vline!(p, [v[argmax(v[!,:spectral_contrast]), :nce]], color = i)
    end

    function convertAltimeterFragment(f::SplineDetailedFrag)
        getIntensity(frag, getKnots)
    end

    using Profile
using PProf

# Collect a profile
Profile.clear()
    @profile begin
        ms_table_path_to_psms_path = quantSearch(
            frag_err_dist_dict,
            pid_to_cv_fold,
            prec_to_irt,
            quant_psms_folder,
            rt_index_paths,
            bin_rt_size,
            rt_irt,
            irt_errs,
            quad_model_dict,
            _ISOTOPE_TRACE_TYPE_,
            chromatograms,
            file_path_to_parsed_name,
            [MS_TABLE_PATHS[1]],
            params_,
            spec_lib,
            ionMatches,
            ionMisses,
            IDtoCOL,
            ionTemplates,
            iso_splines,
            complex_scored_PSMs,
            complex_unscored_PSMs,
            complex_spectral_scores,
            precursor_weights
            );
    
    end
    pprof()




end


using DataFrames
using StatsPlots
using Statistics
using Plots

# Function to create bins and assign bin centers
function assign_bins(x, bin_width)
    bins = minimum(x):bin_width:maximum(x)
    bin_centers = bins[1:end-1] .+ bin_width/2
    bin_indices = searchsortedlast.(Ref(bins), x)
    return bin_indices, bin_centers, bins
end

# Create separate plots for each charge state
bin_width = 75  # adjust this value based on your data distribution
plots = Plots.Plot[]  # Initialize as Plot array with proper namespace

for charge_val in sort(unique(test_out.charge))
    # Filter data for this charge
    charge_data = filter(row -> row.charge == charge_val, test_out)
    
    # Create bins for prec_mz
    bin_indices, bin_centers, bins = assign_bins(charge_data.prec_mz, bin_width)
    
    # Add bin information to the dataframe
    charge_data_with_bins = DataFrame(
        nce = charge_data.nce,
        prec_mz = charge_data.prec_mz,
        bin_index = bin_indices
    )
    
    # Create the boxplot
    p = @df charge_data_with_bins boxplot(
        string.(bin_indices),  # x-axis: bin indices
        :nce,                 # y-axis: NCE values
        fillcolor=:lightblue,
        color=:blue,
        xlabel="Precursor m/z",
        ylabel="NCE",
        title="NCE Distribution by Precursor m/z Bins (Charge $charge_val)",
        legend=false,
        xrotation=45  # Rotate x-axis labels as part of the initial plot
    )
    
    # Add sample size annotations
    for bin_idx in unique(bin_indices)
        subset = charge_data_with_bins[charge_data_with_bins.bin_index .== bin_idx, :]
        n = nrow(subset)
        mnce = median(subset[!,:nce])
        if n > 0
            y_max = maximum(subset.nce)
            annotate!(bin_idx, y_max + 2,  # adjust the +2 offset as needed
                     text("n=$n \n nce=$mnce", :center, 8))
        end
    end
    
    # Customize x-axis ticks to show actual m/z values
    xticks!(1:length(bin_centers), string.(round.(Int, bin_centers)))
    
    push!(plots, p)
end

# Display all plots in a grid
plot(plots..., layout=(length(plots), 1), size=(800, 300*length(plots)))