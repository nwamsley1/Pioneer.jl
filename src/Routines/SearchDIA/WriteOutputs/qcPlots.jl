# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

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

"""
    _bar_value_label(x) -> String

Format the value annotation drawn inside a QC bar.

Replaces `string(x)[1:min(n, end)]`, which truncated the decimal representation by *character
count* and so silently corrupted small values: an MS2 mass correction of -2.11e-4 rendered as
`"-0.00"` (indistinguishable from zero) and 3.06e-5 as `"3.06e"` (malformed). `%.3g` keeps three
significant figures and switches to exponent notation on its own, so every magnitude stays
readable and no value is ever misreported.
"""
_bar_value_label(x::Integer) = string(x)
function _bar_value_label(x::Real)
    xf = float(x)
    isfinite(xf) || return string(x)
    return @sprintf("%.3g", xf)
end
_bar_value_label(x) = string(x)

"""
    _bar_value_texts(values)

Build the rotated in-bar value annotations. The anchor sits at the bar's tip and `halign = :right`
sends the text back down into the bar, so the label stays inside the plot as long as the bar is
taller than the label is long. Where that does not hold -- the MS2 mass-error plots, whose values
are ~1e-4 and whose bars can be a few percent of the axis -- put the value in the tick label via
`_bar_tick_labels` instead, rather than annotating inside the bar.
"""
function _bar_value_texts(values)
    return [text(_bar_value_label(v), valign = :vcenter, halign = :right, rotation = 90) for v in values]
end

"""
    _bar_tick_labels(names, values)

`"name\nvalue"` x tick labels, for bar plots where the bars are too short to hold the value
annotation. Guarantees the value is legible regardless of bar height.
"""
function _bar_tick_labels(names, values)
    return [string(names[i], "\n", _bar_value_label(values[i])) for i in eachindex(names, values)]
end

"""
    _bar_bottom_margin(labels)

Bottom margin for a bar plot whose x tick labels are rotated 45 degrees.

Plots.jl does not grow the plot area to fit rotated tick labels, so with the default margin a name
like "MA5117_Rep1" ran off the bottom edge with its leading characters clipped. A 45-degree label
occupies `len * cos(45deg)` character widths vertically; at roughly 2mm per character that is
`1.4 * len` mm, plus a few mm of padding. Capped so a pathological name cannot squeeze the axes
down to nothing.
"""
function _bar_bottom_margin(labels)
    isempty(labels) && return 3Plots.mm
    maxline = 0
    maxlines = 1
    for l in labels
        parts = split(string(l), '\n')
        maxlines = max(maxlines, length(parts))
        for part in parts
            maxline = max(maxline, length(part))
        end
    end
    pad = 3 + ceil(Int, 1.4 * maxline) + 4 * (maxlines - 1)
    return min(pad, 45) * Plots.mm
end

#Group psms by file of origin
function qcPlots(
    precursors_wide_path,
    precursors_long_path,
    protein_groups_wide_path,
    params_,
    precursors,
    parsed_fnames,
    qc_plot_folder,
    MS_TABLE_PATHS,
    irt_rt,
    frag_err_dist_dict,
    valid_file_indices
)
    # Get q-value column name
    qval_col = :qval

    short_fnames = shortenFileNames(parsed_fnames)
    #grouped_precursors = groupby(best_psms, :file_name)
    #grouped_protein_groups = groupby(protein_quant, :file_name)
    #Number of files to parse
    precursors_wide = Arrow.Table(precursors_wide_path)
    precursors_long = Arrow.Table(precursors_long_path)
    # Only `getIdCounts` reads this frame, and it touches exactly two columns: the file column and
    # the q-value column. Materialising the whole table (DataFrame(Tables.columntable(...)) copies
    # every column, and expands the dict-encoded string columns into millions of String objects)
    # pulled the entire final output into RAM just to count rows per file -- the one unbounded
    # allocation in the QC path, growing with the size of the experiment.
    # copycols = false wraps the Arrow columns in place, so even these two are not copied.
    precursors_long_df = DataFrame(
        :file_name => precursors_long[:file_name],
        qval_col   => precursors_long[qval_col];
        copycols = false,
    )
    protein_groups_wide = Arrow.Table(protein_groups_wide_path)
    n_files = length(parsed_fnames)
    # Number of files per QC plot page — hardcoded since no user override
    # was ever shipped (was global.output.plots_per_page).
    n_files_per_plot = 12
    #Number of QC plots to build
    n_qc_plots = Int64(round(n_files/n_files_per_plot, RoundUp))
    qc_plots = Plots.Plot[]
    ###############
    #Plot precursor abundance distribution

    # Get counts of precursor IDs per file meeting q-value threshold
    function getIdCounts(
        precursors_long::DataFrame;
        file_column::Symbol = :file_name,
        q_value_column::Symbol = :qval,
        q_value_threshold::Real = 0.01f0
        )

        # Group + count rows passing the q-value threshold in one pass.
        # Folding the threshold into `combine` skips an intermediate
        # `filter`-allocated DataFrame copy of the surviving subset.
        file_counts = combine(
            groupby(precursors_long, file_column),
            q_value_column => (qvals -> count(<=(q_value_threshold), qvals)) => :precursor_count,
        )

        return Dict(
            file_counts[!, file_column] .=> file_counts.precursor_count
        )
    end

    # Create abundance ECDF plots using the same approach as other plots
    function create_abundance_ecdf_plots(
        precursors_wide::Arrow.Table,
        precursors_long::DataFrame,
        parsed_fnames::Vector{String};
        n_files_per_plot::Int = 20,
        file_column::Symbol = :file_name,
        q_value_column::Symbol = :qval,
        q_value_threshold::Real = 0.01f0,
        qc_plot_folder::String = "./qc_plots/"
    )
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

        # Get file names from precursors data to ensure same ordering as other plots
        fname_to_id = getIdCounts(
            precursors_long,
            file_column = file_column,
            q_value_column = q_value_column,
            q_value_threshold = q_value_threshold
        )
        # Use original file order, filtering to only include files that have data
        fnames = [fname for fname in parsed_fnames if haskey(fname_to_id, fname)]

        # Calculate number of plots needed
        n_qc_plots = ceil(Int, length(fnames) / n_files_per_plot)
        plots = Plots.Plot[]

        # Create plots in chunks
        for n in 1:n_qc_plots
            start = (n - 1) * n_files_per_plot + 1
            stop = min(n * n_files_per_plot, length(fnames))

            chunk_fnames = fnames[start:stop]
            chunk_short_names = shortenFileNames(chunk_fnames)

            p = Plots.plot(title = "Precursor Abundance by Rank",
                xlabel = "Log2(precursor rank)",
                ylabel = "Log10(precursor abundance)",
                legend = :outertopright, layout = (1, 1))

            for (parsed_fname, short_fname) in zip(chunk_fnames, chunk_short_names)
                try
                    sampled_range, sorted_precursor_abundances = getColumnECDF(precursors_wide[Symbol(parsed_fname)])
                    plot!(p,
                        collect(sampled_range),
                        sorted_precursor_abundances,
                        subplot = 1,
                        label = short_fname
                    )
                catch
                    continue
                end
            end

            push!(plots, p)
        end

        return plots
    end

    # Create abundance ECDF plots using the same file ordering as other plots
    abundance_ecdf_plots = create_abundance_ecdf_plots(
        precursors_wide,
        precursors_long_df,
        parsed_fnames,
        n_files_per_plot = n_files_per_plot,
        file_column = :file_name,
        q_value_column = qval_col,
        q_value_threshold = _resolve_q_value_threshold(params_.global_settings),
        qc_plot_folder = qc_plot_folder
    )
    append!(qc_plots, abundance_ecdf_plots)

    ###############
    # Plot precursor IDs using long-format data
    ###############

    # Plot precursor IDs for a subset of files
    function plotPrecursorIDBarPlot(
        fname_to_id::Dict,
        parsed_fnames::Vector;
        title::String = "Precursor ID's per File",
        f_out::String = "./test.pdf"
        )

        # Filter to only files that have data (successful files)
        successful_files = String[]
        successful_ids = Int[]
        for fname in parsed_fnames
            if haskey(fname_to_id, fname)
                push!(successful_files, fname)
                push!(successful_ids, fname_to_id[fname])
            end
        end

        # If no successful files, create empty plot with message
        if isempty(successful_files)
            p = Plots.plot(
                title = title,
                legend = :none,
                layout = (1, 1),
                annotations = [(0.5, 0.5, text("No successful files", :center))]
            )
            return p
        end

        p = Plots.plot(
            title = title,
            legend = :none,
            layout = (1, 1)
        )

        Plots.bar!(
            p,
            successful_files,
            successful_ids,
            subplot = 1,
            texts = _bar_value_texts(successful_ids),
            xrotation = 45,
            bottom_margin = _bar_bottom_margin(successful_files),
            left_margin = 5Plots.mm,
        )

        return p
    end

    # Create multiple plots with chunking
    function create_precursor_id_plots(
        precursors_long::DataFrame,
        parsed_fnames::Vector{String};
        n_files_per_plot::Int = 20,
        file_column::Symbol = :ms_file_idx,
        q_value_column::Symbol = :qval,
        q_value_threshold::Real = 0.01f0,
        qc_plot_folder::String = "./qc_plots/"
    )
        # Get counts of IDs per file once
        fname_to_id = getIdCounts(
            precursors_long,
            file_column = file_column,
            q_value_column = q_value_column,
            q_value_threshold = q_value_threshold
        )
        # Use original file order for all successful files; show 0 when missing
        fnames = parsed_fnames
        # Calculate number of plots needed
        n_qc_plots = ceil(Int, length(fnames) / n_files_per_plot)
        plots = Plots.Plot[]
        
        # Create plots in chunks
        for n in 1:n_qc_plots
            start = (n - 1) * n_files_per_plot + 1
            stop = min(n * n_files_per_plot, length(fnames))
            
            push!(plots, plotPrecursorIDBarPlot(
                fname_to_id,
                fnames[start:stop],
                title = "Precursor ID's per File",
                f_out = joinpath(qc_plot_folder, "precursor_ids_barplot_" * string(n) * ".pdf")
            ))
        end
        return plots
    end

    # Create precursor ID plots using the same file ordering approach
    id_plots = create_precursor_id_plots(precursors_long_df,
                                parsed_fnames,
                                n_files_per_plot=n_files_per_plot,
                                file_column = :file_name,
                                q_value_column = qval_col,
                                q_value_threshold = _resolve_q_value_threshold(params_.global_settings),
                                qc_plot_folder = qc_plot_folder)
    append!(qc_plots, id_plots)
    ###############
    #Plot Protein Group IDs

    # Create protein group ID plots using the same approach as precursor plots
    function create_protein_group_id_plots(
        protein_groups_wide::Arrow.Table,
        precursors_long::DataFrame,
        parsed_fnames::Vector{String};
        n_files_per_plot::Int = 20,
        file_column::Symbol = :file_name,
        q_value_column::Symbol = :qval,
        q_value_threshold::Real = 0.01f0,
        qc_plot_folder::String = "./qc_plots/"
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

        # Create protein group counts only for successful files
        protein_fname_to_id = Dict{String, Int64}()
        for fname in parsed_fnames
            try
                count = getColumnIDs(protein_groups_wide[Symbol(fname)])
                if count > 0  # Only add files with actual protein groups
                    protein_fname_to_id[fname] = count
                end
            catch
                # Skip files that don't have data
            end
        end

        # Filter to only successful files
        successful_fnames = [fname for fname in parsed_fnames if haskey(protein_fname_to_id, fname)]

        # If no successful files, return empty array
        if isempty(successful_fnames)
            return Plots.Plot[]
        end

        # Calculate number of plots needed
        n_qc_plots = ceil(Int, length(successful_fnames) / n_files_per_plot)
        plots = Plots.Plot[]

        # Create plots in chunks
        for n in 1:n_qc_plots
            start = (n - 1) * n_files_per_plot + 1
            stop = min(n * n_files_per_plot, length(successful_fnames))

            chunk_fnames = successful_fnames[start:stop]
            chunk_ids = [protein_fname_to_id[fname] for fname in chunk_fnames]
            chunk_short_names = shortenFileNames(chunk_fnames)

            p = Plots.plot(title = "Protein Group ID's per File", legend = :none, layout = (1, 1))

            Plots.bar!(p,
                chunk_short_names,
                chunk_ids,
                subplot = 1,
                texts = _bar_value_texts(chunk_ids),
                xrotation = 45,
                bottom_margin = _bar_bottom_margin(chunk_fnames),
                left_margin = 5Plots.mm,
            )

            push!(plots, p)
        end

        return plots
    end

    # Create protein group ID plots using the same file ordering as precursor plots
    protein_id_plots = create_protein_group_id_plots(
        protein_groups_wide,
        precursors_long_df,
        parsed_fnames,
        n_files_per_plot = n_files_per_plot,
        file_column = :file_name,
        q_value_column = qval_col,
        q_value_threshold = _resolve_q_value_threshold(params_.global_settings),
        qc_plot_folder = qc_plot_folder
    )
    append!(qc_plots, protein_id_plots)
    
    ###############
    #Plot missed cleavage rate

    # Create missed cleavage plots using the same approach as precursor plots
    function create_missed_cleavage_plots(
        precursors_wide::Arrow.Table,
        precursors_long::DataFrame,
        precursors,
        parsed_fnames::Vector{String};
        n_files_per_plot::Int = 20,
        file_column::Symbol = :file_name,
        q_value_column::Symbol = :qval,
        q_value_threshold::Real = 0.01f0,
        qc_plot_folder::String = "./qc_plots/"
    )
        function getMissedCleavageRate(
            abundance,
            precursor_idx,
            missed_cleavages)

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

        # Get file names from precursors data to ensure same ordering as other plots
        fname_to_id = getIdCounts(
            precursors_long,
            file_column = file_column,
            q_value_column = q_value_column,
            q_value_threshold = q_value_threshold
        )
        # Use original file order, filtering to only include files that have data
        fnames = [fname for fname in parsed_fnames if haskey(fname_to_id, fname)]

        # Create missed cleavage rates using the same file order
        fname_to_cleavage_rate = Dict{String, Float32}()
        for fname in fnames
            try
                fname_to_cleavage_rate[fname] = 100 * getMissedCleavageRate(
                    precursors_wide[Symbol(fname)],
                    precursors_wide[:precursor_idx],
                    getMissedCleavages(precursors)
                )
            catch e
                @user_warn "Failed to compute missed cleavage rate for $(fname): $e"
                fname_to_cleavage_rate[fname] = 0.0f0
            end
        end

        # Calculate number of plots needed
        n_qc_plots = ceil(Int, length(fnames) / n_files_per_plot)
        plots = Plots.Plot[]

        # Create plots in chunks
        for n in 1:n_qc_plots
            start = (n - 1) * n_files_per_plot + 1
            stop = min(n * n_files_per_plot, length(fnames))

            chunk_fnames = fnames[start:stop]
            chunk_rates = [fname_to_cleavage_rate[fname] for fname in chunk_fnames]
            chunk_short_names = shortenFileNames(chunk_fnames)

            p = Plots.plot(title = "Missed Cleavage Percentage", legend = :none, layout = (1, 1))

            Plots.bar!(p,
                chunk_short_names,
                chunk_rates,
                subplot = 1,
                texts = _bar_value_texts(chunk_rates),
                xrotation = 45,
                bottom_margin = _bar_bottom_margin(chunk_fnames),
                left_margin = 5Plots.mm,
            )

            push!(plots, p)
        end

        return plots
    end

    # Create missed cleavage plots using the same file ordering as other plots
    missed_cleavage_plots = create_missed_cleavage_plots(
        precursors_wide,
        precursors_long_df,
        precursors,
        parsed_fnames,
        n_files_per_plot = n_files_per_plot,
        file_column = :file_name,
        q_value_column = qval_col,
        q_value_threshold = _resolve_q_value_threshold(params_.global_settings),
        qc_plot_folder = qc_plot_folder
    )
    append!(qc_plots, missed_cleavage_plots)

    ###############
    #Plot TIC
    function plotTIC(
        ms_table_paths::Vector{String},
        short_fnames::Vector{String};  # Now expects only the chunk of file names
        title::String = "precursor_abundance_qc",
        f_out::String = "./test.pdf"
        )

        p = Plots.plot(title = title,
        legend=:outertopright, layout = (1, 1))

        for (id, ms_table_path) in enumerate(ms_table_paths)
            short_fname = short_fnames[id]  # Now correctly indexed
            ms_table = BasicMassSpecData(ms_table_path)
            ms1_indices = findall(getMsOrders(ms_table) .== 1)
            rts  = [getRetentionTime(ms_table, i) for i in ms1_indices]
            tics = [getTIC(ms_table, i) for i in ms1_indices]

            plot!(p, rts, tics,
                    subplot = 1,
                    xlabel = "retention time (min)",
                    label = short_fname,
                    alpha = 0.3
                    )
        end
        return p
    end

    for n in 1:n_qc_plots
        start = (n - 1)*n_files_per_plot + 1
        stop = min(n*n_files_per_plot, length(parsed_fnames))
        p = plotTIC(
            MS_TABLE_PATHS[start:stop],
            short_fnames[start:stop],  # Pass only the chunk of file names
            title = "TIC plot",
            f_out = joinpath(qc_plot_folder, "tic_plot_"*string(n)*".pdf")
        )
        push!(qc_plots, p)
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
            plot!(p, 
                    irt_range,
                    irt_rt.(irt_range),
                    xlabel = "iRT",
                    ylabel = "RT (min)",
                    label = basename(key),
                    alpha = 0.5,
                    lw = 3
                    )
        end
        return p
    end
    ###############
    #Plot precursor IDs
    function plotMS2MassErrors(
        gpsms::GroupedDataFrame,
        parsed_fnames::AbstractVector{<:AbstractString},
        short_fnames::AbstractVector{<:AbstractString};
        title::String = "precursor_abundance_qc",
        f_out::String = "./test.pdf"
        )

        p = Plots.plot(title = title,
                        legend=:none, layout = (1, 1))

        for (parsed_fname, short_fname) in zip(parsed_fnames, short_fnames)
            boxplot!(p, 
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

        return p
    end


    ###############
    #Plot precursor IDs
    function plotMS2MassErrorCorrection(
        parsed_fnames::AbstractVector{<:AbstractString},
        short_fnames::AbstractVector{<:AbstractString},
        frag_err_dist_dict::Dict{Int64, AbstractMassErrorModel},
        file_ids::Vector{Int64};
        title::String = "precursor_abundance_qc",
        f_out::String = "./test.pdf"
        )

        p = Plots.plot(title = title,
                        legend=:none, layout = (1, 1))

        parsed_fnames = [parsed_fnames[file_id] for file_id in file_ids]
        short_fnames = [short_fnames[file_id] for file_id in file_ids]
        mass_corrections = [getMassCorrection(frag_err_dist_dict[file_id]) for file_id in file_ids]

        _mass_corr_labels = _bar_tick_labels(short_fnames, mass_corrections)
        Plots.bar!(p,
        _mass_corr_labels,
        mass_corrections,
        subplot = 1,
        yflip = true,
        xrotation = 45,
        bottom_margin = _bar_bottom_margin(_mass_corr_labels),
        left_margin = 5Plots.mm,
        )

        return p
    end

    function plotMS2MassErrorCorrectionFixed(
        parsed_fnames::Vector{String},
        short_fnames::Vector{String},
        frag_err_dist_dict::Dict{Int64, AbstractMassErrorModel},
        actual_file_indices::Vector{Int64};
        title::String = "precursor_abundance_qc",
        f_out::String = "./test.pdf"
        )

        p = Plots.plot(title = title,
                        legend=:none, layout = (1, 1))

        # Use pre-sliced file names and actual file indices for mass error lookup
        mass_corrections = [getMassCorrection(frag_err_dist_dict[file_id]) for file_id in actual_file_indices]

        _mass_corr_labels = _bar_tick_labels(short_fnames, mass_corrections)
        Plots.bar!(p,
        _mass_corr_labels,
        mass_corrections,
        subplot = 1,
        yflip = true,
        xrotation = 45,
        bottom_margin = _bar_bottom_margin(_mass_corr_labels),
        left_margin = 5Plots.mm,
        )

        return p
    end


    for n in 1:n_qc_plots
        start = (n - 1)*n_files_per_plot + 1
        stop = min(n*n_files_per_plot, length(short_fnames))

        # Get actual file indices for this range
        actual_file_indices = valid_file_indices[start:stop]
        # Pass both sequential indices and actual file indices
        p = plotMS2MassErrorCorrectionFixed(
            parsed_fnames[start:stop],
            short_fnames[start:stop],
            frag_err_dist_dict,
            actual_file_indices,
            title = "MS2 Mass Error Correction",
            f_out = joinpath(qc_plot_folder, "mass_err_correction_"*string(n)*".pdf")
        )
        push!(qc_plots, p)

    end


    ###############
    #Plot precursor IDs
    function plotFWHM(
        gpsms::GroupedDataFrame,
        parsed_fnames::AbstractVector{<:AbstractString},
        short_fnames::AbstractVector{<:AbstractString};
        title::String = "fwhm_qc",
        f_out::String = "./test.pdf"
        )

        p = plot(title = title,
                        legend=:none, layout = (1, 1))


        for (parsed_fname, short_fname) in zip(parsed_fnames, short_fnames)
            psms =  gpsms[(file_name = parsed_fname,)]
            boxplot!(p, 
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

        return p
    end

    output_path = joinpath(qc_plot_folder, "QC_PLOTS.pdf")
    try
        if isfile(output_path)
            safeRm(output_path, nothing)
        end
    catch e
        @user_warn "Could not clear existing file: $e"
    end
    save_multipage_pdf(qc_plots, output_path)
end
