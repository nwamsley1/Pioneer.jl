function shortenFileNames(parsed_fnames,
    max_length = 20)

    result = Vector{String}(undef, length(parsed_fnames))

    for (i, s) in pairs(parsed_fnames)
        if length(s) <= max_length
            # String is already within length, keep it as is
            result[i] = s
        else
            # Calculate how many characters remain for the sides after inserting "..."
            side_space = max_length - 3  # e.g. if max_length=10, side_space=7
            left_len = div(side_space, 2)  # integer division
            right_len = side_space - left_len

            # left chunk from start, right chunk from end
            left_part = s[1:left_len]
            right_part = s[end - right_len + 1:end]
            result[i] = left_part * "..." * right_part
        end
    end

    return result
end

#Group psms by file of origin
function qcPlots(
    precursors_wide_path,
    protein_groups_wide_path,
    params_,
    precursors,
    parsed_fnames,
    qc_plot_folder,
    MS_TABLE_PATHS,
    irt_rt,
    frag_err_dist_dict
)
    short_fnames = shortenFileNames(parsed_fnames)
    #grouped_precursors = groupby(best_psms, :file_name)
    #grouped_protein_groups = groupby(protein_quant, :file_name)
    #Number of files to parse
    precursors_wide = Arrow.Table(precursors_wide_path)
    protein_groups_wide = Arrow.Table(protein_groups_wide_path)
    n_files = length(parsed_fnames)
    n_files_per_plot = Int64(params_.output[:plots_per_page])
    #Number of QC plots to build 
    n_qc_plots = Int64(round(n_files/n_files_per_plot, RoundUp))
    ###############
    #Plot precursor abundance distribution
    function plotAbundanceECDF(
        precursors_wide::Arrow.Table,
        parsed_fnames::Any,
        short_fnames::Any;
        title::String = "precursor_abundance_qc",
        f_out::String = "./test.pdf"
        )

        p = Plots.plot(title = title,
        xlabel = "Log2(precursor rank)",
        ylabel = "Log10(precursor abundance)",
        legend=:outertopright, layout = (1, 1))
        function getColumnECDF(
            abundance::AbstractVector{Union{Missing, Float32}})
            non_missing_count = 0
            for i in range(1, length(abundance))
                if !ismissing(abundance[i])
                    non_missing_count += 1
                end
            end
            sorted_precursor_abundances = zeros(eltype(abundance), non_missing_count)
            non_missing_count = 0
            for i in range(1, length(abundance))
                if !ismissing(abundance[i])
                    non_missing_count += 1
                    sorted_precursor_abundances[non_missing_count] = abundance[i]
                end
            end
            sort!(sorted_precursor_abundances, rev=true)
            sampled_range = 1:20:length(sorted_precursor_abundances)
            return log2.(sampled_range), [log10(x) for x in sorted_precursor_abundances[sampled_range]]
        end

        for (parsed_fname, short_fname) in zip(parsed_fnames, short_fnames)
            try
            sampled_range, sorted_precursor_abundances = getColumnECDF(precursors_wide[Symbol(parsed_fname)])
            Plots.plot!(p, 
                    collect(sampled_range),
                    sorted_precursor_abundances,
                    #ylim = (0, log2(maximum(sorted_precursor_abundances))),
                    subplot = 1,
                    label = short_fname
                    )
            catch
                continue
            end
        end
        savefig(p, f_out)
    end

    for n in 1:n_qc_plots
        start = (n - 1)*n_files_per_plot + 1
        stop = min(n*n_files_per_plot, length(parsed_fnames))

        plotAbundanceECDF(
            precursors_wide,
            [x for x in parsed_fnames[start:stop]],
            [x for x in short_fnames[start:stop]],
            title = "Precursor Abundance by Rank",
            f_out = joinpath(qc_plot_folder, "precursor_abundance_qc_"*string(n)*".pdf")
        )
    end

    ###############
    #Plot precursor IDs
    function plotPrecursorIDBarPlot(
        precursors_wide::Arrow.Table,
        parsed_fnames::Any,
        short_fnames::Any;
        title::String = "precursor_abundance_qc",
        f_out::String = "./test.pdf"
        )

        function getColumnIDs(
            abundance::AbstractVector{Union{Missing, Float32}})
            non_missing_count = 0
            for i in range(1, length(abundance))
                if !ismissing(abundance[i])
                    non_missing_count += 1
                end
            end
            return non_missing_count
        end
        ids = zeros(Int64, length(parsed_fnames))
        for (i, fname) in enumerate(parsed_fnames)
            try
            ids[i] = getColumnIDs(precursors_wide[Symbol(fname)])
            catch
                continue
            end
        end
        p = Plots.plot(title = title,
                        legend=:none, layout = (1, 1))

        Plots.bar!(p, 
        short_fnames,
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
            precursors_wide,
            [x for x in parsed_fnames[start:stop]],
            [x for x in short_fnames[start:stop]],
            title = "Precursor ID's per File",
            f_out = joinpath(qc_plot_folder, "precursor_ids_barplot_"*string(n)*".pdf")
        )

    end
    ###############
    #Plot Protein Group IDs
    
    function plotProteinGroupIDBarPlot(
        protein_groups_wide::Arrow.Table,
        parsed_fnames::Any,
        short_fnames::Any;
        title::String = "protein_abundance_qc",
        f_out::String = "./test.pdf"
        )
        function getColumnIDs(
            abundance::AbstractVector{Union{Missing, Float32}})
            non_missing_count = 0
            for i in range(1, length(abundance))
                if !ismissing(abundance[i])
                    non_missing_count += 1
                end
            end
            return non_missing_count
        end
        ids = zeros(Int64, length(parsed_fnames))
        for (i, fname) in enumerate(parsed_fnames)
            try
                ids[i] = getColumnIDs(protein_groups_wide[Symbol(fname)])
            catch
                continue
            end
        end
        p = Plots.plot(title = title,
                        legend=:none, layout = (1, 1))

        Plots.bar!(p, 
        short_fnames,
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

        plotProteinGroupIDBarPlot(
            protein_groups_wide,
            [x for x in parsed_fnames[start:stop]],
            [x for x in short_fnames[start:stop]],
            title = "Protein Group ID's per File",
            f_out = joinpath(qc_plot_folder, "protein_group_ids_barplot_"*string(n)*".pdf")
        )

    end
    
    ###############
    #Plot missed cleavage rate
    function plotMissedCleavageRate(
        precursors_wide::Arrow.Table,
        parsed_fnames::Any,
        short_fnames::Any;
        title::String = "precursor_abundance_qc",
        f_out::String = "./test.pdf"
        )
        function getMissedCleavageRate(
            abundance::AbstractVector{Union{Missing, Float32}},
            precursor_idx::AbstractVector{Union{UInt32}},
            missed_cleavages::AbstractVector{UInt8})

            non_missing_count = 0
            missed_cleavage_count = 0
            for i in range(1, length(abundance))
                if !ismissing(abundance[i])
                    non_missing_count += 1
                    missed_cleavage_count += missed_cleavages[precursor_idx[i]]
                end
            end
            return missed_cleavage_count/non_missing_count
        end
        ids = zeros(Float32, length(parsed_fnames))
        for (i, fname) in enumerate(parsed_fnames)
            try
                ids[i] = 100*getMissedCleavageRate(precursors_wide[Symbol(fname)],
                precursors_wide[:precursor_idx],
                getMissedCleavages(precursors))#[:missed_cleavages])
            catch
                continue
            end
        end
        p = Plots.plot(title = title,
                        legend=:none, layout = (1, 1))

        Plots.bar!(p, 
        short_fnames,
        ids,
        subplot = 1,
        texts = [text(string(x)[1:min(6,length(string(x)))], valign = :vcenter, halign = :right, rotation = 90) for x in ids],
        xrotation = 45,
        )

        savefig(p, f_out)
    end

    for n in 1:n_qc_plots
        start = (n - 1)*n_files_per_plot + 1
        stop = min(n*n_files_per_plot, length(parsed_fnames))

        plotMissedCleavageRate(
            precursors_wide,
            [x for x in parsed_fnames[start:stop]],
            [x for x in short_fnames[start:stop]],
            title = "Missed Cleavage Percentage",
            f_out = joinpath(qc_plot_folder, "missed_cleavage_plot_"*string(n)*".pdf")
        )

    end

    ###############
    #Plot TIC
    function plotTIC(
        ms_table_paths::Vector{String},
        parsed_fnames::Any,
        short_fnames::Any;
        title::String = "precursor_abundance_qc",
        f_out::String = "./test.pdf"
        )

        p = Plots.plot(title = title,
        legend=:outertopright, layout = (1, 1))

        for (id, ms_table_path) in enumerate(ms_table_paths)
            short_fname = short_fnames[id]
            ms_table = BasicMassSpecData(ms_table_path)
            ms1_scans = getMsOrders(ms_table).==1

            
            Plots.plot!(p, 
            [getRetentionTime(ms_table, scan_idx) for scan_idx in range(1, length(ms_table)) if getMsOrder(ms_table, scan_idx)==1],
            [getTIC(ms_table, scan_idx) for scan_idx in range(1, length(ms_table)) if getMsOrder(ms_table, scan_idx)==1],
                    subplot = 1,
                    xlabel = "retention time (min)",
                    label = short_fname,
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
            parsed_fnames,
            short_fnames,
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
        legend=:outertopright)

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
#=
    for n in 1:n_qc_plots
        start = (n - 1)*n_files_per_plot + 1
        stop = min(n*n_files_per_plot, length(parsed_fnames))
        plotRTAlign(
            irt_rt,
            collect(keys(irt_rt))[start:stop],
            title = "RT Alignment Plot",
            f_out = joinpath(qc_plot_folder, "rt_align_plot_"*string(n)*".pdf")
        )
    end
=#
    ###############
    #Plot precursor IDs
    function plotMS2MassErrors(
        gpsms::GroupedDataFrame,
        parsed_fnames::Any,
        short_fnames::Any;
        title::String = "precursor_abundance_qc",
        f_out::String = "./test.pdf"
        )

        p = Plots.plot(title = title,
                        legend=:none, layout = (1, 1))

        for (parsed_fname, short_fname) in zip(parsed_fnames, short_fnames)
            Plots.boxplot!(p, 
            [short_fname],
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
        parsed_fnames::Any,
        short_fnames::Any,
        frag_err_dist_dict::Dict{Int64, MassErrorModel},
        file_ids::Vector{Int64};
        title::String = "precursor_abundance_qc",
        f_out::String = "./test.pdf"
        )

        p = Plots.plot(title = title,
                        legend=:none, layout = (1, 1))

        parsed_fnames = [parsed_fnames[file_id] for file_id in file_ids]
        short_fnames = [short_fnames[file_id] for file_id in file_ids]
        mass_corrections = [getMassCorrection(frag_err_dist_dict[file_id]) for file_id in file_ids]

        Plots.bar!(p, 
        short_fnames,
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
        stop = min(n*n_files_per_plot, length(short_fnames))

        plotMS2MassErrorCorrection(
            parsed_fnames,
            short_fnames,
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
        parsed_fnames::Any,
        short_fnames::Any;
        title::String = "fwhm_qc",
        f_out::String = "./test.pdf"
        )

        p = Plots.plot(title = title,
                        legend=:none, layout = (1, 1))


        for (parsed_fname, short_fname) in zip(parsed_fnames, short_fnames)
            psms =  gpsms[(file_name = parsed_fname,)]
            Plots.boxplot!(p, 
            [short_fname],
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
    output_path = joinpath(qc_plot_folder, "QC_PLOTS.pdf")
    try
        if isfile(output_path)
            rm(output_path)
        end
    catch e
        @warn "Could not clear existing file: $e"
    end
    plots_to_merge = [joinpath(qc_plot_folder, x) for x in readdir(qc_plot_folder) if endswith(x, ".pdf")]
    if length(plots_to_merge)>1
        merge_pdfs(plots_to_merge, 
                    output_path, cleanup=true)
    end
end