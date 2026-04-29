module PDFGenerator
import ..Plots
const GR = Plots.GR

function create_multipage_pdf(plots::Vector, dest::String)
    # Force a non-interactive GR workstation before backend initialization so
    # headless environments do not attempt to open gksqt/Qt windows.
    withenv("GKSwstype" => "100", "GKS_WSTYPE" => "100") do
        Plots.gr()  # This will only initialize once
        GR.beginprint(dest)
        try
            for p in plots
                # Suppress all display output
                redirect_stdout(devnull) do
                    redirect_stderr(devnull) do
                        display(p)
                    end
                end
            end
        finally
            GR.endprint()
        end
    end
end
end

function render_plot(spec::MultiSeriesPlotSpec)
    plot_kwargs = Dict{Symbol, Any}(
        :title => spec.title,
        :xlabel => spec.xlabel,
        :ylabel => spec.ylabel,
        :grid => spec.grid,
    )
    if spec.legend !== nothing
        plot_kwargs[:legend] = spec.legend
    end
    if spec.right_margin_px > 0
        plot_kwargs[:right_margin] = spec.right_margin_px * Plots.px
    end

    plot_obj = Plots.plot(; plot_kwargs...)

    for series in spec.series
        series_kwargs = Dict{Symbol, Any}(
            :seriestype => series.seriestype,
            :alpha => series.alpha,
            :linewidth => series.linewidth,
            :linestyle => series.linestyle,
            :markersize => series.markersize,
        )
        if series.label !== nothing
            series_kwargs[:label] = series.label
        end
        if series.color !== nothing
            series_kwargs[:color] = series.color
        end
        Plots.plot!(plot_obj, series.x, series.y; series_kwargs...)
    end

    if spec.x_min !== nothing && spec.x_max !== nothing
        Plots.xlims!(plot_obj, spec.x_min, spec.x_max)
    end
    if spec.y_min !== nothing && spec.y_max !== nothing
        Plots.ylims!(plot_obj, spec.y_min, spec.y_max)
    end

    for annotation in spec.annotations
        Plots.annotate!(
            plot_obj,
            annotation.x,
            annotation.y,
            Plots.text(annotation.text, annotation.halign, annotation.size),
        )
    end

    return plot_obj
end

function render_plot(spec::RtAlignmentPlotSpec)
    plot_obj = Plots.plot(
        spec.rt,
        spec.irt;
        seriestype = :scatter,
        title = spec.title,
        xlabel = "Retention Time RT (min)",
        ylabel = "Indexed Retention Time iRT (min)",
        label = nothing,
        alpha = 0.1,
        size = 100 * [13.3, 7.5],
        grid = isempty(spec.rt),
    )

    if !isempty(spec.curve_rt) && !isempty(spec.curve_irt)
        Plots.plot!(
            plot_obj,
            spec.curve_rt,
            spec.curve_irt;
            lw = spec.curve_width,
            color = spec.curve_color,
            ls = spec.curve_style,
            label = spec.curve_label,
        )
    end

    if spec.annotation_text !== nothing &&
       spec.annotation_x !== nothing &&
       spec.annotation_y !== nothing
        Plots.annotate!(
            plot_obj,
            spec.annotation_x,
            spec.annotation_y,
            Plots.text(spec.annotation_text, :center, 10),
        )
    end

    return plot_obj
end

function render_plot(spec::MassErrorPlotSpec)
    spread = max(spec.left_spread, spec.right_spread)

    plot_obj = if spec.orientation == :horizontal
        if isempty(spec.errs)
            Plots.plot(
                title = spec.title,
                xlabel = spec.x_label,
                ylabel = spec.y_label,
                size = 100 * [13.3, 7.5],
                grid = true,
                yflip = true,
                orientation = :h,
                ylim = (spec.center - spread - 10.0f0, spec.center + spread + 10.0f0),
            )
        else
            bins = LinRange(spec.center - 2.0f0 * spec.left_spread, spec.center + 2.0f0 * spec.right_spread, 50)
            Plots.histogram(
                spec.errs;
                orientation = :h,
                yflip = true,
                title = spec.title,
                xlabel = spec.x_label,
                ylabel = spec.y_label,
                label = nothing,
                bins = bins,
                ylim = (spec.center - 2.0f0 * spec.left_spread, spec.center + 2.0f0 * spec.right_spread),
            )
        end
    elseif isempty(spec.errs)
        Plots.plot(
            Float32[],
            Float32[];
            xlabel = spec.x_label,
            ylabel = spec.y_label,
            title = spec.title,
            legend = :topright,
            grid = true,
        )
    else
        Plots.histogram(
            spec.errs;
            bins = 50,
            xlabel = spec.x_label,
            ylabel = spec.y_label,
            title = spec.title,
            label = "Mass Errors",
            alpha = 0.7,
            color = :blue,
        )
    end

    lower_bound = spec.center - spec.left_spread
    upper_bound = spec.center + spec.right_spread

    if spec.orientation == :horizontal
        Plots.hline!(plot_obj, [spec.center]; color = :black, lw = 2, label = spec.center_label)
        Plots.hline!(plot_obj, [lower_bound]; color = :red, lw = 1, ls = :dash, label = nothing)
        Plots.hline!(plot_obj, [upper_bound]; color = :red, lw = 1, ls = :dash, label = spec.bound_label)
    else
        Plots.vline!(plot_obj, [spec.center]; color = :red, linewidth = 2, label = spec.center_label)
        Plots.vline!(plot_obj, [lower_bound, upper_bound]; color = :green, linestyle = :dash, linewidth = 1.5, label = spec.bound_label)
    end

    if spec.annotation_text !== nothing &&
       spec.annotation_x !== nothing &&
       spec.annotation_y !== nothing
        Plots.annotate!(
            plot_obj,
            spec.annotation_x,
            spec.annotation_y,
            Plots.text(spec.annotation_text, :left, 8),
        )
    end

    return plot_obj
end

function save_multipage_pdf(plots::AbstractVector{<:Plots.Plot}, dest::String)
    ensure_directory_exists(dest)
    PDFGenerator.create_multipage_pdf(plots, dest)
    return dest
end

function save_multipage_pdf(plot_specs::AbstractVector{<:AbstractPioneerPlotSpec}, dest::String)
    ensure_directory_exists(dest)
    PDFGenerator.create_multipage_pdf(map(render_plot, plot_specs), dest)
    return dest
end
