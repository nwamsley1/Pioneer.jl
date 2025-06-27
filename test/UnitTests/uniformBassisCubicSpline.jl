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


#Example data. Approximate a sine curve 
N = 200
t = collect(LinRange(0.0, 4*π, N))
u = sin.(t) 
#u .+= randn(N)./50
plot(t, u, seriestype=:scatter)
test_spline = UniformSpline(u, t, 3, 20)
#plot!(LinRange(-1, 4*π+1, 500), test_spline.(LinRange(0-1, 4*π+1, 500)), linewidth = 3, alpha = 0.5)
#using DataInterpolations
#test_old = BSplineApprox(u, t, 4, 20, :Uniform, :Uniform, extrapolate = true)
#plot!(LinRange(-1, 4*π+1, 500), test_old.(LinRange(0-1, 4*π+1, 500)), linewidth = 3, alpha = 0.5)
UniformSpline(u, t, 3, 3)

@test maximum(test_spline.(t) .- u) .< 1e-3