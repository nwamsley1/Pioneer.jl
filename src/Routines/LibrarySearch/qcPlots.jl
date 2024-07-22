#Group psms by file of origin
function qcPlots(
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
    grouped_precursors = groupby(best_psms, :file_name)
    grouped_protein_groups = groupby(protein_quant, :file_name)
    #Number of files to parse
    n_files = length(grouped_precursors)
    n_files_per_plot = Int64(params_[:qc_plot_params]["n_files_per_plot"])
    #Number of QC plots to build 
    n_qc_plots = Int64(round(n_files/n_files_per_plot, RoundUp))
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
            psms = gpsms[(file_name = parsed_fname,)]

            sorted_precursor_abundances = sort(psms[!,column], 
                                                rev = true)

            sampled_range = 1:20:length(sorted_precursor_abundances)

            Plots.plot!(p, sampled_range,
                    [log2.(max(sorted_precursor_abundances[x], 0.0)) for x in sampled_range],
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
            ids[i] = size(gpsms[(file_name = parsed_fname,)], 1)
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
    #Plot precursor IDs
    function plotProteinIDBarPlot(
        gpsms::GroupedDataFrame,
        parsed_fnames::Vector{<:Any};
        title::String = "precursor_abundance_qc",
        f_out::String = "./test.pdf"
        )

        p = Plots.plot(title = title,
                        legend=:none, layout = (1, 1), show = true)

        ids = Vector{Int64}(undef, length(parsed_fnames))
        for (i, parsed_fname) in enumerate(parsed_fnames)
            ids[i] = size(gpsms[(file_name = parsed_fname,)], 1)
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

        plotProteinIDBarPlot(
            grouped_protein_groups,
            parsed_fnames[start:stop],
            title = "Precursor ID's per File",
            f_out = joinpath(qc_plot_folder, "protein_ids_barplot_"*string(n)*".pdf")
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
            psms = gpsms[(file_name = parsed_fname,)]
        
            missed_abundance = sum(psms[!,:missed_cleavage].*psms[!,abundance_col])

            missed_cleavage_rates[i] = round(missed_abundance/sum(psms[!,abundance_col]), sigdigits = 3)

        end

        Plots.bar!(p, 
        parsed_fnames,
        missed_cleavage_rates,
        subplot = 1,
        texts = [text(string(x)[1:min(4, length(string(x)))], valign = :vcenter, halign = :right, rotation = 90) for x in  missed_cleavage_rates],
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


    ###############
    #Plot precursor IDs
    function plotMS2MassErrors(
        gpsms::GroupedDataFrame,
        parsed_fnames::Vector{<:Any};
        title::String = "precursor_abundance_qc",
        f_out::String = "./test.pdf"
        )

        p = Plots.plot(title = title,
                        legend=:none, layout = (1, 1), show = true)

        for parsed_fname in parsed_fnames
            Plots.boxplot!(p, 
            [parsed_fname],
            gpsms[(file_name = parsed_fname,)][!,:error_norm],
            subplot = 1,
            outliers = false,
            fill = :white,
            color = :black,
            lw = 3,
            #texts = [text(string(x), valign = :vcenter, halign = :right, rotation = 90) for x in ids],
            xrotation = 45
            )
        end

        savefig(p, f_out)
    end

    #=
    for n in 1:n_qc_plots
        start = (n - 1)*n_files_per_plot + 1
        stop = min(n*n_files_per_plot, length(parsed_fnames))

        plotMS2MassErrors(
            grouped_precursors,
            parsed_fnames[start:stop],
            title = "MS2 Mass Errors Distribution",
            f_out = joinpath(qc_plot_folder, "mass_err_boxplots_"*string(n)*".pdf")
        )

    end
    =#

    ###############
    #Plot precursor IDs
    function plotMS2MassErrorCorrection(
        file_id_to_parsed_name::Dict{Int64, String},
        frag_err_dist_dict::Dict{Int64, MassErrorModel},
        file_ids::Vector{Int64};
        title::String = "precursor_abundance_qc",
        f_out::String = "./test.pdf"
        )

        p = Plots.plot(title = title,
                        legend=:none, layout = (1, 1), show = true)

        parsed_fnames = [file_id_to_parsed_name[file_id] for file_id in file_ids]
        mass_corrections = [getMassCorrection(frag_err_dist_dict[file_id]) for file_id in file_ids]

        Plots.bar!(p, 
        parsed_fnames,
        mass_corrections,
        subplot = 1,
        yflip = true,
        texts = [text(string(x)[1:min(5, length(string(x)))], valign = :vcenter, halign = :right, rotation = 90) for x in mass_corrections],
        xrotation = 45,
        )

        savefig(p, f_out)
    end


    for n in 1:n_qc_plots
        start = (n - 1)*n_files_per_plot + 1
        stop = min(n*n_files_per_plot, length(parsed_fnames))

        plotMS2MassErrorCorrection(
            file_id_to_parsed_name,
            frag_err_dist_dict,
            collect(start:stop),
            title = "MS2 Mass Error Correction",
            f_out = joinpath(qc_plot_folder, "mass_err_correction_"*string(n)*".pdf")
        )

    end


    ###############
    #Plot precursor IDs
    function plotFWHM(
        gpsms::GroupedDataFrame,
        parsed_fnames::Vector{<:Any};
        title::String = "fwhm_qc",
        f_out::String = "./test.pdf"
        )

        p = Plots.plot(title = title,
                        legend=:none, layout = (1, 1), show = true)


        for parsed_fname in parsed_fnames
            psms =  gpsms[(file_name = parsed_fname,)]
            Plots.boxplot!(p, 
            [parsed_fname],
            psms[psms[!,:points_above_FWHM_01].>1,:FWHM],
            subplot = 1,
            outliers = false,
            fill = :white,
            color = :black,
            ylabel = "FWHM (min)",
            lw = 3,
            #texts = [text(string(x), valign = :vcenter, halign = :right, rotation = 90) for x in ids],
            xrotation = 45
            )
        end

        savefig(p, f_out)
    end

    #=
    for n in 1:n_qc_plots
        start = (n - 1)*n_files_per_plot + 1
        stop = min(n*n_files_per_plot, length(parsed_fnames))

        plotFWHM(
            grouped_precursors,
            parsed_fnames[start:stop],
            title = "FWHM Distribution",
            f_out = joinpath(qc_plot_folder, "FWHM_"*string(n)*".pdf")
        )

    end
    =#
    println("qc_plot_folder ", qc_plot_folder)
    plots_to_merge = [joinpath(qc_plot_folder, x) for x in readdir(qc_plot_folder) if endswith(x, ".pdf")]
    println("plots_to_merge ", plots_to_merge)

    if length(plots_to_merge)>1
        merge_pdfs(plots_to_merge, 
                    joinpath(qc_plot_folder, "QC_PLOTS.pdf"), cleanup=true)
    end
end