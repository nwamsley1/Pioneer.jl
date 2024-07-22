function makeOutputDirectories(
    OUT_DIR::String
)
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
