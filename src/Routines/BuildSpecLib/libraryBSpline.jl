function _B_(x::T, k::Int, i::Int, t::NTuple{N, T}) where {N, T<:AbstractFloat}
    # Base case for k=0 (constant basis function)
    if k == 0
        return T(t[i] <= x < t[i+1])
    end
    
    # Compute first recursive term
    c1 = if t[i+k] == t[i]
        zero(T)
    else
        ((x - t[i]) / (t[i+k] - t[i])) * B(x, k-1, i, t)
    end
    
    # Compute second recursive term
    c2 = if t[i+k+1] == t[i+1]
        zero(T)
    else
        ((t[i+k+1] - x) / (t[i+k+1] - t[i+1])) * B(x, k-1, i+1, t)
    end
    
    return c1 + c2
end

function speval_naieve(x::T, t::NTuple{N, T}, c::NTuple{M, T}, k::Int) where {M,N,T<:AbstractFloat}
    # Compute number of basis functions
    n = length(t) - k - 1
    # Validate inputs (Julia's equivalent to Python's assert)
    @assert n >= k+1 && length(c) >= n "Invalid input sizes"
    
    # Compute B-spline value
    v = zero(T)
    for i in range(1, n)
        v += c[i] * B(x, k, i, t)
    end
    return v
end
knots = Tuple(x for x in SVector{length(knots), Float32}(knots))
df[1,:coefficients]
EvalBSpline(25.0f0, knots, df[1,:coefficients],3) 
p = plot()
for k in range(1, 30)
v = [splev_(Float64.(collect(knots)),
                Float64.(collect(df[k,:coefficients])), 
                3, x) for x in tbins]
plot!(tbins, v)
end

v, _ = splev!(
    Float64.(collect(knots)),
    8,
    Float64.(collect(df[1,:coefficients])),
    3,
    collect(LinRange(20, 40, 100)),
)

knots = Float64[6.0, 13.0, 20.0, 27.0, 34.0, 41.0, 48.0, 55.0]
coeffs = Float64[1.6181915e-6, 7.382022e-6, 7.887343e-5, 0.00023642876]

# Evaluate the B-spline at x = 25.0, with degree k = 3
result = splev(knots,
                Float64.(collect(df[1,:coefficients])), 
                3, 25.0)


tbins = LinRange(20, 40, 100)
v, _ = splev!(
    Float64.(collect(knots)),
    8,
    Float64.(test_coeff),
    3,
    collect(tbins),
)
plot(tbins, v)
tbins = LinRange(20, 40, 100)
v, _ = splev!(
    Float64.(collect(knots)),
    8,
    Float64.(collect(df[1,:coefficients])),
    3,
    collect(tbins),
)
plot!(tbins, v)

batch_path = ordered_altimeter_json_paths[5000]
df, frags_per_prec, knot_vec = parseBatchToTable(JSON.parse(read(batch_path, String))["outputs"], SplineCoefficientModel("test"))
#filterEachPrecursor!(df, model_type, intensity_threshold = intensity_threshold)

batch_id = 5000
precs_per_batch = 1000
pid = 1




pbins = LinRange(0.0f0, 60.0f0, 100)
pid = 1
color = 0

p = plot(legend = :outertopleft)
pid = 2
for k in range((pid-1)*50+1, pid*50)
    if startswith(df[k,:annotation], "I")
        color = 1
        continue
    elseif startswith(df[k,:annotation], "p")
        color = 2
    elseif startswith(df[k,:annotation], "y")
        color = 3
    elseif startswith(df[k,:annotation], "b")
        color = 4
    else
        continue
    end
    println("$pid $k ", df[k,:annotation])
    #plot!(p, pbins ,
    #    [EvalBSpline(x, knots, df[k,:coefficients], 3) for x in pbins], color = color,
    #    label = df[k,:annotation], show = true, xlim = (30, 31)
    #    )
    coefs = Float64.(collect( df[k,:coefficients]))
    plot!(p, pbins ,
    [splev(knots, coefs, 3, Float64(x)) for x in pbins], color = color,
    label = df[k,:annotation], show = true, xlim = (20, 41)
    )

end 

ttable[precs_per_batch*(batch_id - 1) + pid,:]
pid += 1

AAAAAAAAAAAAAAAAGATCLER 

#=
knots = Float32.(JSON.parse(read(batch_path, String))["outputs"][1]["data"]);
tus = UniformSpline(
     SVector(df[1,:coefficients]),
     3,
     6.0f0,
     55.0f0,
     7.0f0
)
     =
     =#

