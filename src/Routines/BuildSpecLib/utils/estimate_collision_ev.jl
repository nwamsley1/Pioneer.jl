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

function get_mz_to_ev_interp(
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
        @user_warn "Collision eV data missing from raw file. Resorting to defaults"
        return missing
    end
end