using KernelDensity
kernel_dist(::Type{Laplace}, w::Real) = Laplace(0.0, w)

function KDEmapping(X::Vector{T}, Y::Vector{T}; n::Int = 200, bandwidth::AbstractFloat = 1.0, w = 11) where {T<:AbstractFloat}
    x_grid = LinRange(minimum(X), maximum(X), n)
    y_grid = LinRange(minimum(Y), maximum(Y), n)
    ys = zeros(T, n)
    z = zeros(T, (n, n))
    min_x, max_x = minimum(X), maximum(X)
    bound_x = (max_x - min_x)/2
    min_y, max_y = minimum(Y), maximum(Y)
    min_x, min_y = zero(T), zero(T)
    bound_y = (max_y - min_y)/2
    B = kde((X, Y), bandwidth = (bandwidth, bandwidth), boundary = ((min_x - bound_x, max_x + bound_x), (min_y - bound_y, max_y + bound_y)), kernel = Laplace) #Uses Silverman's rule by default
    ik = KernelDensity.InterpKDE(B)
    #Get KDE
    for i in eachindex(x_grid), j in eachindex(y_grid)
        z[i, j] = Distributions.pdf(ik, x_grid[i], y_grid[j])
    end

    min_z = sum(z)/length(z)
    for i in eachindex(z)
        if z[i] < min_z
            z[i] = zero(eltype(z))
        end
    end
    #Monotonic increasing walk along ridge
    max_j = 1
    for i in eachindex(x_grid)
        j = argmax(@view(z[i,:]))
        if y_grid[j] > y_grid[max_j]
            max_j = j
        end
        ys[i] = y_grid[max_j]
    end
    w = isodd(n÷5) ? n÷5 : n÷5 + 1
    return LinearInterpolation(x_grid, 
                                    ys,
                                    #savitzky_golay(ys, w, 3).y, 
                                    extrapolation_bc = Line())
    #return LinearInterpolation(x_grid, ys, extrapolation_bc = Line())
end

function plotRTAlign(RT::Vector{T}, iRT::Vector{T}, rt_map::Any; f_out::String = "./") where {T<:AbstractFloat}
    p = Plots.plot(RT, iRT, seriestype=:scatter,
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
    savefig(p, f_out*"_rt_alignment.pdf")


end