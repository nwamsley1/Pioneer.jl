x = collect(LinRange(0, 4*π, 100))
y = sin.(x)
test_spline = BSplineApprox(y, x, 4, 5, :Uniform, :Uniform, extrapolate = true)
plot_bins = LinRange(0, 5*pi, 1000)
plot(plot_bins, sin.(plot_bins))
plot!(plot_bins, test_spline.(plot_bins))


test_spline = BSplineApprox(y, x, 4, 10, :Uniform, :Uniform, extrapolate = true)
plot_bins = LinRange(0, 5*pi, 1000)
plot(plot_bins, sin.(plot_bins))
plot!(plot_bins, test_spline.(plot_bins))


function _interpolate(A::BSplineApprox{<:AbstractVector{<:Number}}, t::Number)
    println("INFUNC!")
    t < A.t[1] && return A.u[1], 1
    t > A.t[end] && return A.u[end], lastindex(t)
    # change t into param [0 1]
    idx = searchsortedlast(A.t, t)
    println("idx $idx")
    idx == length(A.t) ? idx -= 1 : nothing
    t = A.p[idx] + (t - A.t[idx]) / (A.t[idx + 1] - A.t[idx]) * (A.p[idx + 1] - A.p[idx])
    #For uniformly spaced, just proportion of distance of query from min to max. 
    println("t $t")
    n = length(A.t)
    N = spline_coefficients(A.h, A.d, A.k, t)
    ucum = zero(eltype(A.u))
    for i in 1:(A.h)
        ucum += N[i] * A.c[i]
    end
    ucum, idx
end

function spline_coefficients(n, d, k, u::Number)
    N = zeros(eltype(u), n)
    M = zeros(eltype(u), (n, d))
    if u == k[1]
        N[1] = one(u)
    elseif u == k[end]
        N[end] = one(u)
    else
        i = findfirst(x -> x > u, k) - 1
        N[i] = one(u)
        for deg in 1:d
            N[i - deg] = (k[i + 1] - u) / (k[i + 1] - k[i - deg + 1]) * N[i - deg + 1]
            for j in (i - deg + 1):(i - 1)
                N[j] = (u - k[j]) / (k[j + deg] - k[j]) * N[j] +
                       (k[j + deg + 1] - u) / (k[j + deg + 1] - k[j + 1]) * N[j + 1]
            end
            N[i] = (u - k[i]) / (k[i + deg] - k[i]) * N[i]
            M[:,d] = copy(N)
        end
    end
    N
end

@test abs(first(_interpolate(test_spline, 10.0)) - test_spline(10.0)) < 1e-10


_interpolate(test_spline, 4*π-0.3)
