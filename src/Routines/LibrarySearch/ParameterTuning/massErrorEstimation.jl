function ModelMassErrs(ppm_errs::Vector{Float32};
                       frag_err_quantile::Float32 = 0.01f0,
                       out_fdir::String = "./",
                       out_fname::String = "mass_err_estimate")

    bins = LinRange(minimum(ppm_errs), maximum(ppm_errs), 100)
    mass_err = median(ppm_errs)
    ppm_errs = ppm_errs .- mass_err
    l_bound = quantile(ppm_errs, frag_err_quantile)
    r_bound = quantile(ppm_errs, 1 - frag_err_quantile)
    errs = ppm_errs .+ mass_err
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
    n = length(errs)
    p = Plots.histogram(errs,
                    orientation = :h, 
                    yflip = true,
                    #seriestype=:scatter,
                    title = plot_title*"\n n = $n",
                    xlabel = "Count",
                    ylabel = "Mass Error (ppm)",
                    label = nothing,
                    bins = bins,
                    ylim = (minimum(ppm_errs)-2, maximum(ppm_errs)+2),
                    topmargin =15mm,
                    #bottommargin = 10mm,
                    )

    Plots.hline!([l_bound + mass_err, mass_err, r_bound + mass_err], label = nothing, color = :black, lw = 2)
    l_err = l_bound + mass_err
    Plots.annotate!(last(xlims(p)), l_err, text("$l_err", :black, :right, :bottom, 12))
    r_err = r_bound + mass_err
    Plots.annotate!(last(xlims(p)), r_err, text("$r_err", :black, :right, :bottom, 12))
    Plots.annotate!(last(xlims(p)), mass_err, text("$mass_err", :black, :right, :bottom, 12))
    savefig(p, joinpath(out_fdir, out_fname)*".pdf")

    MassErrorModel(
                    Float32(mass_err),
                    (Float32(abs(l_bound)), Float32(abs(r_bound)))
                    )
end

