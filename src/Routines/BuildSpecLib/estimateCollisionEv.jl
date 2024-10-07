function getMzToEvInterp(
    ms_table_path::String,
    plot_dir::String,
)
    ms_table= DataFrame(Arrow.Table(ms_table_path))
    ms_table = dropmissing(unique(ms_table[!,[:centerMz,:collisionEnergyEvField]]))
    filter!(x->x.collisionEnergyEvField>0, ms_table)
    sort!(ms_table, :centerMz)

    if size(ms_table, 1) > 3
        mz_to_ev_interp = linear_interpolation(
            ms_table[!,:centerMz], 
            ms_table[!,:collisionEnergyEvField],
            extrapolation_bc=Line()) 

        tbins = LinRange(
                    minimum(ms_table[!,:centerMz]),
                    maximum(ms_table[!,:centerMz]), 100)
        p = plot(
            ms_table[!,:centerMz], 
            ms_table[!,:collisionEnergyEvField],
            seriestype=:scatter,
            title = "Precursor m/z vs. HCD eV",
            xlabel = "Precursor m/z",
            ylabel = "HCD eV",
            label = "empirical"
            )
        plot!(p, tbins, mz_to_ev_interp.(tbins), label = "interpolation")
        savefig(joinpath(plot_dir, "mz_to_ev_plot.pdf"))
        return mz_to_ev_interp
    else
        @warn "Collision eV data missing from raw file. Resorting to defaults"
        return missing
    end
end