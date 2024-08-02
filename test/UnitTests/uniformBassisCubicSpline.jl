
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