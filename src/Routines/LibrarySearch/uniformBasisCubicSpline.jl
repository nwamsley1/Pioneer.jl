#Uniform Cubic Smoothing Spline
#Given number of evenly spaced control points and data (time, value)
#solve coefficients for B-spline basis. Then build a speedy implementation. 


#Example data
N = 200
t = collect(LinRange(0.0, 4*π, N))
u = sin.(t) 
u .+= randn(N)./50
plot(t, u, seriestype=:scatter)


using Polynomials, StaticArrays, Plots

CubicBSplineBasis = NTuple{4, Polynomial}([
    Polynomial([0, 0, 0, 1])/6, #b0
    Polynomial([1, 3, 3, -3])/6, #b1
    Polynomial([4, 0, -6, 3])/6, #b2
    Polynomial([1, -3, 3, -1])/6, #b3
])

d = 4 #For cubic spline 
n = 20 #Number of evenly spaced knots (includes endpoints)
N = n + 3 #Number of basis functions 
t0, tn = 0.0, 4*π
knots = collect(LinRange(t0, tn, n))
bin_width = (tn - t0)/(n - 1)

function buildDesignMat(t::Vector{T}, 
                        knots::Vector{T},
                        bin_width::T,
                        spline_basis::NTuple{4, Polynomial}) where {T<:AbstractFloat}
    
    function fillDesignMatRow!(X::Matrix{T}, 
                                row::Int,
                                knot_idx::Int,
                                u::T,
                                spline_basis::NTuple{4, Polynomial}) where {T<:AbstractFloat}
        i = 4
        #println(" t - knot_val: ", t - knot_val)
        for col in range(knot_idx, knot_idx + 3)
            X[row, col] = spline_basis[i](u)
            i -= 1
        end
    end

    X = zeros(T, (length(t), length(knots) + 3))
    for (i, t) in enumerate(t)
        knot_idx = min(Int64((t - first(knots))÷bin_width)+1, length(knots))
        fillDesignMatRow!(
            X,
            i,
            knot_idx,
            (t-knots[knot_idx])/bin_width,
            spline_basis
        )
    end

    return X
end
function buildPieceWise(
                        knots::Vector{T},
                        bin_width::T,
                        spline_basis::NTuple{4, Polynomial}) where {T<:AbstractFloat}
    
    function fillDesignMatRow!(X::Matrix{Polynomial}, 
                                row::Int,
                                knot_idx::Int,
                                spline_basis::NTuple{4, Polynomial}) where {T<:AbstractFloat}
        i = 4
        #println(" t - knot_val: ", t - knot_val)
        for col in range(knot_idx, knot_idx + 3)
            X[row, col] = spline_basis[i]
            i -= 1
        end
    end

    X = zeros(Polynomial, (length(knots), length(knots) + 3))
    for (i, t) in enumerate(knots)
        t = t + bin_width/2
        knot_idx = min(Int64((t - first(knots))÷bin_width)+1, length(knots))
        fillDesignMatRow!(
            X,
            i,
            knot_idx,
            spline_basis
        )
    end

    return X
end

X = buildDesignMat(t, collect(knots), bin_width, CubicBSplineBasis)
c = X\u 
plot(t, u, seriestype=:scatter)
plot!(t, X*c)
X = buildPieceWise(knots, bin_width, CubicBSplineBasis)
test_poly = X*c
#plot!(LinRange(0, 0.5, 50), test_poly[1].(LinRange(0, 1.0, 50)))
#plot(t, X[:,1])
n_coeffs = n*d
spline_coeffs = SVector{n_coeffs}(vcat([polynomial.coeffs for polynomial in test_poly]...))

struct CustomSpline{N, T<:AbstractFloat} 
    coeffs::SVector{N, T}
    degree::Int64
    first::T
    last::T
    bin_width::T
end

test_spline = CustomSpline(
    spline_coeffs,
    3,
    0.0, 
    4*π,
    (4*π - 0.0)/19
)

function (s::CustomSpline)(t::U) where {U<:AbstractFloat}
    t = max(min(t, s.last), s.first)
    idx = floor(Int32, 
                    (t - s.first)/s.bin_width
                )
    u = (t - (s.first + s.bin_width*(idx)))/s.bin_width
    x = zero(U)
    coeff = idx*(s.degree + 1) + 1
    c = one(U)
    x += s.coeffs[coeff]*c
    c *= u
    coeff += 1
    x += s.coeffs[coeff]*c
    c *= u
    coeff += 1
    x += s.coeffs[coeff]*c
    c *= u
    coeff += 1
    x += s.coeffs[coeff]*c
    c *= 1
    #for i in range(first_coeff, first_coeff + s.degree)
    #    x += s.coeffs[i]*c
    #    c *= u
    #end
    return x
end
plot(t, u, seriestype=:scatter)
plot!(LinRange(-1, 4*π+1, 500), test_spline.(LinRange(0-1, 4*π+1, 500)), linewidth = 3, alpha = 0.5)
#using DataInterpolations
#test_old = BSplineApprox(u, t, 4, 20, :Uniform, :Uniform, extrapolate = true)
#plot!(LinRange(-1, 4*π+1, 500), test_old.(LinRange(0-1, 4*π+1, 500)), linewidth = 3, alpha = 0.5)