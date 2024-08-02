function plotRTAlign(RT::Vector{T}, 
                    iRT::Vector{T}, 
                    rt_map::Any; 
                    out_fdir::String = "./",
                    out_fname::String = "rt_align_plot") where {T<:AbstractFloat}
    n = length(RT)

    plot_title = ""
    n = 0
    for i in range(1, length(out_fname))
        n += 1
        if n > 24
            n = 1
            plot_title *= "\n"
        end
        plot_title *= out_fname[i]
    end
    n = length(RT)
    p = Plots.plot(RT, iRT, seriestype=:scatter,
                        title = plot_title*"\n n = $n",
                        xlabel = "Retention Time RT (min)",
                        ylabel ="Indexed Retention Time iRT (min)",
                        label = nothing,
                        size = 100*[13.3, 7.5]
            ,
            fontsize = 24,
            titlefontsize = 24,
            legendfontsize = 24,
            tickfontsize = 24,
            guidefontsize = 24,
            margin = 10Plots.mm,
            dpi = 300)

    Plots.plot!(p, (LinRange(minimum(RT), maximum(RT), 100)), 
            rt_map.(LinRange(minimum(RT), maximum(RT), 100)),
            lw = 6.0,
            label = "RT Spline")

    savefig(p, joinpath(out_fdir, out_fname)*".pdf")


end