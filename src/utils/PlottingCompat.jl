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

"""
    PlottingCompat.jl
    
Compatibility layer for lazy-loaded plotting functionality.
Ensures backward compatibility while preventing Qt conflicts.
"""
module PlottingCompat

using ..Pioneer: ensure_plotting_loaded

"""
    @with_plots expr

Macro to ensure plotting libraries are loaded before executing an expression.
Use this when you need to access Plots or StatsPlots directly.

# Example
```julia
@with_plots begin
    p = plot(x, y)
    savefig(p, "output.png")
end
```
"""
macro with_plots(expr)
    quote
        ensure_plotting_loaded()
        $(esc(expr))
    end
end

export @with_plots

end # module