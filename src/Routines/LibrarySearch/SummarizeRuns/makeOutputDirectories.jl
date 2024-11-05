function makeOutputDirectories(
    OUT_DIR::String,
    params::Dict{String, Any}
)

    if !isabspath(OUT_DIR)
        OUT_DIR = joinpath(@__DIR__, "../../../../", OUT_DIR)
    end
    if !isdir(OUT_DIR)
        mkpath(OUT_DIR)
    end

    search_folder = joinpath(OUT_DIR)
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
    results_folder = joinpath(OUT_DIR, "RESULTS")
    if !isdir(results_folder)
        mkpath(results_folder)
    end
    params_folder = joinpath(results_folder, "PARAMS")
    if !isdir(params_folder )
        mkpath(params_folder)
    end
    #Write params to folder 
    open(joinpath(params_folder, "config.json"),"w") do f 
        JSON.print(f, params)
    end

    temp_folder = joinpath(search_folder, "temp")
    if !isdir(temp_folder)
        mkpath(temp_folder)
    end

    return qc_plot_folder, rt_alignment_folder, mass_err_estimation_folder, results_folder, temp_folder
end