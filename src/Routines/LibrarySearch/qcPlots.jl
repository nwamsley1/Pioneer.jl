#Group psms by file of origin
grouped_precursors = groupby(best_psms_passing, :parsed_fname)
#Number of files to parse
n_files = length(grouped_precursors)
n_files_per_plot = Int64(params_[:qc_plot_params]["n_files_per_plot"])
#Number of QC plots to build 
n_qc_plots = Int64(round(n_files/n_files_per_plot))
###############
#Plot precursor abundance distribution
function plotAbundanceECDF(
    gpsms::GroupedDataFrame,
    parsed_fnames::Vector{<:Any};
    column::Symbol = peak_area,
    title::String = "precursor_abundance_qc",
    f_out::String = "./test.pdf"
    )
    p = Plots.plot(title = title,
    legend=:outertopright, layout = (1, 1), show = true)

    for parsed_fname in parsed_fnames
        psms = gpsms[(parsed_fname = parsed_fname,)]

        sorted_precursor_abundances = sort(psms[!,column], 
                                            rev = true)

        sampled_range = 1:20:length(sorted_precursor_abundances)

        Plots.plot!(p, sampled_range,
                [log2.(sorted_precursor_abundances[x]) for x in sampled_range],
                ylim = (0, log2(maximum(sorted_precursor_abundances))),
                subplot = 1,
                label = parsed_fname
                )
    end
    savefig(p, f_out)
end

for n in 1:n_qc_plots
    start = (n - 1)*n_files_per_plot + 1
    stop = min(n*n_files_per_plot, length(parsed_fnames))
    plotAbundanceECDF(
        grouped_precursors,
        parsed_fnames[start:stop],
        column = :peak_area,
        title = "Precursor Abundance by Rank",
        f_out = joinpath(qc_plot_folder, "precursor_abundance_qc_"*string(n)*".pdf")
    )
end

###############
#Plot precursor IDs
function plotPrecursorIDBarPlot(
    gpsms::GroupedDataFrame,
    parsed_fnames::Vector{<:Any};
    title::String = "precursor_abundance_qc",
    f_out::String = "./test.pdf"
    )

    p = Plots.plot(title = title,
                    legend=:none, layout = (1, 1), show = true)

    ids = Vector{Int64}(undef, length(parsed_fnames))
    for (i, parsed_fname) in enumerate(parsed_fnames)
        ids[i] = size(gpsms[(parsed_fname = parsed_fname,)], 1)
    end

    Plots.bar!(p, 
    parsed_fnames,
    ids,
    subplot = 1,
    texts = [text(string(x), valign = :vcenter, halign = :right, rotation = 90) for x in ids],
    xrotation = 45,
    )

    savefig(p, f_out)
end


for n in 1:n_qc_plots
    start = (n - 1)*n_files_per_plot + 1
    stop = min(n*n_files_per_plot, length(parsed_fnames))

    plotPrecursorIDBarPlot(
        grouped_precursors,
        parsed_fnames[start:stop],
        title = "Precursor ID's per File",
        f_out = joinpath(qc_plot_folder, "precursor_ids_barplot_"*string(n)*".pdf")
    )

end

###############
#Plot missed cleavage rate
function plotMissedCleavageRate(
    gpsms::GroupedDataFrame,
    parsed_fnames::Vector{<:Any};
    abundance_col::Symbol = :weight,
    title::String = "precursor_abundance_qc",
    f_out::String = "./test.pdf"
    )

    p = Plots.plot(title = title,
                    legend=:none, layout = (1, 1), show = true)

    missed_cleavage_rates = Vector{Float64}(undef, length(parsed_fnames))

    for (i, parsed_fname) in enumerate(parsed_fnames)
        psms = gpsms[(parsed_fname = parsed_fname,)]
    
        missed_abundance = sum(psms[!,:missed_cleavage].*psms[!,abundance_col])

        missed_cleavage_rates[i] = round(missed_abundance/sum(psms[!,abundance_col]), sigdigits = 3)

    end

    Plots.bar!(p, 
    parsed_fnames,
    missed_cleavage_rates,
    subplot = 1,
    texts = [text(string(x)[1:4], valign = :vcenter, halign = :right, rotation = 90) for x in  missed_cleavage_rates],
    xrotation = 45,
    )

    savefig(p, f_out)
end

for n in 1:n_qc_plots
    start = (n - 1)*n_files_per_plot + 1
    stop = min(n*n_files_per_plot, length(parsed_fnames))

    plotMissedCleavageRate(
        grouped_precursors,
        parsed_fnames[start:stop],
        title = "Missed Cleavage Rate",
        f_out = joinpath(qc_plot_folder, "missed_cleavage_plot_"*string(n)*".pdf")
    )

end

###############
#Plot TIC
function plotTIC(
    ms_table_paths::Vector{String},
    file_id_to_parsed_name::Dict{Int64, String};
    title::String = "precursor_abundance_qc",
    f_out::String = "./test.pdf"
    )

    p = Plots.plot(title = title,
    legend=:outertopright, layout = (1, 1), show = true)

    for (id, ms_table_path) in enumerate(ms_table_paths)
        parsed_fname = file_id_to_parsed_name[id]
        ms_table = Arrow.Table(ms_table_path)
        ms1_scans = ms_table[:msOrder].==1


        Plots.plot!(p, 
                ms_table[:retentionTime][ms1_scans],
                ms_table[:TIC][ms1_scans],
                subplot = 1,
                xlabel = "retention time (min)",
                label = parsed_fname,
                alpha = 0.3
                )
    end
    savefig(p, f_out)
end

for n in 1:n_qc_plots
    start = (n - 1)*n_files_per_plot + 1
    stop = min(n*n_files_per_plot, length(parsed_fnames))
    plotTIC(
        MS_TABLE_PATHS[start:stop],
        file_id_to_parsed_name,
        title = "TIC plot",
        f_out = joinpath(qc_plot_folder, "tic_plot_"*string(n)*".pdf")
    )
end


###############
#Plot RT_align
function plotRTAlign(
    iRT_RT::Dictionary{String, Any},
    keys::Vector{String};
    title::String = "rt_align_plot",
    f_out::String = "./test.pdf"
    )

    p = Plots.plot(title = title,
    legend=:outertopright, show = true)

    for key in keys

        irt_range = LinRange(0, 35, 250)
        irt_rt = iRT_RT[key]
        Plots.plot!(p, 
                irt_range,
                irt_rt.(irt_range),
                xlabel = "iRT",
                ylabel = "RT (min)",
                label = basename(key),
                alpha = 0.5,
                lw = 3
                )
    end
    savefig(p, f_out)
end

for n in 1:n_qc_plots
    start = (n - 1)*n_files_per_plot + 1
    stop = min(n*n_files_per_plot, length(parsed_fnames))
    plotRTAlign(
        iRT_RT,
        collect(keys(iRT_RT))[start:stop],
        title = "RT Alignment Plot",
        f_out = joinpath(qc_plot_folder, "rt_align_plot_"*string(n)*".pdf")
    )
end


merge_pdfs([joinpath(qc_plot_folder, x) for x in readdir(qc_plot_folder) if endswith(x, ".pdf")], 
    joinpath(qc_plot_folder, "QC_PLOTS.pdf"), cleanup=true)