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

#############
#Simmulate data 

iso_splines = parseIsoXML(joinpath(dirname(dirname(@__DIR__)),"assets", "IsotopeSplines_10kDa_21isotopes.xml"));

SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/ASMS_2024/ASTRAL_THREE_PROTEOME/unispec_chronologer_1mc_1var_by_052724/spec_lib/pioneer_lib/"
 
test_rqm = RazoQuadParams(1.0, 2.0, 2.0, 10.0)
Random.seed!(42)
test_data = simmulateQuad(test_rqm, (-3.0, 3.0), (400.0, 1000.0), 1000)
#plot(plot_bins, test_rqm.(plot_bins), lw = 2.0, alpha = 0.5, xlabel = "m/z offset", ylabel = "probability", label = nothing)#"True model")
#plot(test_data[!,:x0], test_data[!,:yt], seriestype=:scatter, ylabel = "log2(f(z0)/f(z1))", xlabel = "m/z offset", label = nothing, alpha = 0.3)


#Estimate parameters of model that generated the data 
test_data[!,:prec_charge] .= 2
fitted_rqm = fitRazoQuadModel(
    2.0,
    test_data[!,:x0], test_data[!,:x1], test_data[!,:yt],
            λ0 = 1e1, #Initial L-M damping coeficient 
        ϵ1 = 1e-6, #Convergence in the gradient
        ϵ2 = 1e-5, #Convergence in the coeficients
        ϵ3 = 1e-6 #Conergence in squared error
    )
#Do estimated parameters agree with original parameters 
@test maximum(abs(fitted_rqm - test_rqm)) < (5e-3)

y = test_rqm.(test_data[!,:x0])
ŷ = fitted_rqm.(test_data[!,:x0])
@test maximum(abs(fitted_rqm - test_rqm)) < (5e-3)
@test (maximum(abs.(y - ŷ)) < 2e-3)
#Fit spline to data?
#=
quad_transmission_spline = splineQuadTransmissionModel(
    test_data[!,:yt],
    test_data[!,:x0],
    test_data[!,:x1], 3, 200
)

#Examine prediction error 
plot_bins = LinRange(-3, 3, 100)
y = test_rqm.(plot_bins)
ŷ_razo = fitted_rqm.(plot_bins)
norm_fac = maximum(quad_transmission_spline.(plot_bins))
ŷ_spline = exp.(quad_transmission_spline.(plot_bins) .- norm_fac)
@test maximum(abs.(ŷ_spline .- y)) < 1e-2
@test maximum(abs.(ŷ_razo .- y)) < 1e-2
=#
#=
test_rqm = RazoQuadParams(1.0, 2.0, 2.0, 10.0)
plot(plot_bins, test_rqm.(plot_bins), lw = 2.0, alpha = 0.5, label = "True model")
plot!(plot_bins, ŷ_razo, lw = 2.0, alpha = 0.5, label = "Fitted model")
plot!(plot_bins, ŷ_spline, lw = 2.0, alpha = 0.5, label = "Spline")
=#


#Now bin the data, and try to recover the original parameters from the binned data 
#=
Random.seed!(41)
sim_x0 = rand(Uniform(-3, 3), 1000);
test_bined = MergeBins(test_data, (-3.0, 3.0), min_bin_size = 50, min_bin_width = 0.1)
fitted_rqm = fitRazoQuadModel(
    3.5,
    test_bined[!,:median_x0],  test_bined[!,:median_x1],  test_bined[!,:median_yt],
    λ0 = 1e1, #Initial L-M damping coeficient 
    ϵ1 = 1e-6, #Convergence in the gradient
    ϵ2 = 1e-5, #Convergence in the coeficients
    ϵ3 = 1e-6 #Conergence in squared error
    )
=#
#=
plot_bins = LinRange(-3, 3, 100)
plot(plot_bins, test_rqm.(plot_bins), lw = 2.0, alpha = 0.5, label = "True model")
plot!(plot_bins, fitted_rqm.(plot_bins), lw = 2.0, alpha = 0.5, label = "Fitted model")
y = test_rqm.(plot_bins)
ŷ = fitted_rqm.(plot_bins)
@test maximum(abs(fitted_rqm - test_rqm)) < (5e-3)
@test (maximum(abs.(y - ŷ)) < 2e-3)
=#


##########
#Fit to actual data 
test_data = DataFrame(Tables.columntable(Arrow.Table(joinpath(@__DIR__, "../../data/quadTransmissionTests/test_to_bin.arrow"))))
#test_data = DataFrame(Tables.columntable(Arrow.Table(joinpath(@__DIR__, "./data/quadTransmissionTests/test_to_bin.arrow"))))
test_bined = MergeBins(test_data,(-2.0, 2.0), min_bin_size = 50, min_bin_width = 0.1)
@test size(test_bined, 1) > 10
@test size(test_bined, 1) < 40
fitted_rqm = fitRazoQuadModel(
    2.0,
    test_bined[!,:median_x0],  test_bined[!,:median_x1],  test_bined[!,:median_yt],
    λ0 = 1e-1, #Initial L-M damping coeficient 
    ϵ1 = 1e-5, #Convergence in the gradient
    ϵ2 = 1e-4, #Convergence in the coeficients
    ϵ3 = 1e-5 #Conergence in squared error
    )
plot_bins = LinRange(-3, 3, 100)
@test fitted_rqm.al < fitted_rqm.ar
@test fitted_rqm.bl < fitted_rqm.br
plot!(plot_bins, fitted_rqm.(plot_bins), lw = 2.0, alpha = 0.5, label = "Fitted model")

#=
fitted_rqm = fitRazoQuadModel(
    1.0,
    test_out[!,:x0], test_out[!,:x1], test_out[!,:yt],
    λ0 = 1e-1, #Initial L-M damping coeficient 
    ϵ1 = 1e-5, #Convergence in the gradient
    ϵ2 = 1e-4, #Convergence in the coeficients
    ϵ3 = 1e-5 #Conergence in squared error
    )


continuous_rqm = fitRazoQuadModel(
    14.0,
    test_out[!,:x0], test_out[!,:x1], test_out[!,:yt],
    λ0 = 1e-1, #Initial L-M damping coeficient 
    ϵ1 = 1e-5, #Convergence in the gradient
    ϵ2 = 1e-4, #Convergence in the coeficients
    ϵ3 = 1e-5 #Conergence in squared error
    )

binned_rqm = fitRazoQuadModel(
    14.0,
    test_bined[!,:median_x0],  test_bined[!,:median_x1],  test_bined[!,:median_yt],
    λ0 = 1e-1, #Initial L-M damping coeficient 
    ϵ1 = 1e-5, #Convergence in the gradient
    ϵ2 = 1e-4, #Convergence in the coeficients
    ϵ3 = 1e-5 #Conergence in squared error
    )

    plot_bins = LinRange(-9, 9, 200)
plot(plot_bins, continuous_rqm.(plot_bins), lw = 2.0, alpha = 0.5, label = "Fitted continuous model")
plot!(plot_bins, binned_rqm.(plot_bins), lw = 2.0, alpha = 0.5, label = "Fitted binned model")
=#