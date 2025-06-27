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
            alpha = 0.1,
            dpi = 300)
    #=
    Plots.plot!(p, (LinRange(minimum(RT), maximum(RT), 100)), 
            rt_map.(LinRange(minimum(RT), maximum(RT), 100)),
            lw = 6.0,
            label = "RT Spline")
    =#
    savefig(p, joinpath(out_fdir, out_fname)*".pdf")


end