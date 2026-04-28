function plotRTAlign(
    RT::Vector{T},
    iRT::Vector{T},
    rt_map::Any;
    out_fdir::String = "./",
    out_fname::String = "rt_align_plot",
) where {T<:AbstractFloat}
    plot_title = ""
    chunk_len = 0
    for idx in eachindex(out_fname)
        chunk_len += 1
        if chunk_len > 24
            chunk_len = 1
            plot_title *= "\n"
        end
        plot_title *= out_fname[idx]
    end

    n_points = length(RT)
    plot_obj = Plots.plot(
        RT,
        iRT;
        seriestype = :scatter,
        title = plot_title * "\n n = $n_points",
        xlabel = "Retention Time RT (min)",
        ylabel = "Indexed Retention Time iRT (min)",
        label = nothing,
        size = 100 * [13.3, 7.5],
        fontsize = 24,
        titlefontsize = 24,
        legendfontsize = 24,
        tickfontsize = 24,
        guidefontsize = 24,
        margin = 10Plots.mm,
        alpha = 0.1,
        dpi = 300,
    )
    Plots.savefig(plot_obj, joinpath(out_fdir, out_fname) * ".pdf")
    return nothing
end
