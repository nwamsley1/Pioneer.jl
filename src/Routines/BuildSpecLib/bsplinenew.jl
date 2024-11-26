

Arrow.write("/Users/n.t.wamsley/Projects/Pioneer.jl/data/example_library_spline_coef.arrow", df[1:10,:])
df = DataFrame(Arrow.Table("/Users/n.t.wamsley/Projects/Pioneer.jl/data/library_spline_tests/example_library_spline_coef.arrow"))
knots = [6.0, 13.0, 20.0, 27.0, 34.0, 41.0, 48.0, 55.0]
degree = 3
xbins = LinRange(20.0, 40.0, 100)
spline_eval_dict = Dict{String, Vector{Float64}}()
for i in range(1, size(df, 1))
    coefs = Float64.(collect(df[i,:coefficients]))
    spl = BSpline(knots, coefs, 3)
    spline_eval_dict[df[i,:annotation]] = spl.(xbins)
end
scipy_spline_eval = Arrow.Table("/Users/n.t.wamsley/Projects/Pioneer.jl/data/library_spline_tests/evaluated_splines_scipy.arrow")
scipy_eval_dict = Dict{String, Vector{Float64}}()
for i in range(1, length(scipy_spline_eval[:spline_values]))
    scipy_eval_dict[scipy_spline_eval[:annotation][i]] = coalesce.(scipy_spline_eval[:spline_values][i], 0.0)
end

plot(xbins, spline_eval_dict["y3^1"])
plot!(xbins, scipy_eval_dict["y3^1"])
@testset "scipy_splines_compare" begin
    for frag_name in keys(spline_eval_dict)
        @test maximum(
            abs.(
                (spline_eval_dict[frag_name] .- scipy_eval_dict[frag_name])./spline_eval_dict[frag_name]
                )
                ) < 1e-6
    end
end

p = plot()
for frag_name in keys(spline_eval_dict)
    plot!(p, xbins, spline_eval_dict[frag_name], label = nothing, show = true)
    plot!(p, xbins, scipy_eval_dict[frag_name], label = frag_name, show = true)
end


# Run the example
knots = [6.0, 13.0, 20.0, 27.0, 34.0, 41.0, 48.0, 55.0]
coefs = [1.6181915043489425e-6, 7.382021976809483e-6, 7.887343235779554e-5, 0.00023642876476515085]
k = 3
tbins = LinRange(20.0, 40.0, 100)
coefs = Float64.(collect(df[k,:coefficients]))

knots = [6.0, 13.0, 20.0, 27.0, 34.0, 41.0, 48.0, 55.0]
coefs = [ 1.4280883533501765e-6, 5.625079666060628e-6, 5.846780913998373e-5, 0.00019516375323291868]
splev(knots, coefs, 3, 37.0)
speval_naive_(37.0, knots, coefs, 3)

coefs = Float64.(collect(df[k,:coefficients]))

tbins = LinRange(20.0, 40.0, 100)
knots = [6.0, 13.0, 20.0, 27.0, 34.0, 41.0, 48.0, 55.0]
coefs = [ 1.4280883533501765e-6, 5.625079666060628e-6, 5.846780913998373e-5, 0.00019516375323291868]
spl = BSpline(knots, coefs, 3)
plot(tbins, spl.(tbins))
pid = 500
p = plot()
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
    coeffs = Float64.(collect(df[k,:coefficients]))
    degree = 3
    spl = BSpline(knots, coeffs, degree)
    v = [spl(x) for x in tbins]
    plot!(tbins, v, color = color,
    label = df[k,:annotation], show = true, xlim = (20, 41)
    )

end
k += 1

knots = [6.0, 13.0, 20.0, 27.0, 34.0, 41.0, 48.0, 55.0]
coeffs = [1.4280883533501765e-6, 5.625079666060628e-6, 5.846780913998373e-5, 0.00019516375323291868]
degree = 3

# Create BSpline
spl = BSpline(knots, coeffs, degree)

# Evaluate at x = 37.0
result = spl(37.0)

v = [speval_naive_(x, knots, coefs, 3) for x in tbins]
plot!(tbins, v)
v = [speval_naieve(x, Tuple(knots), Tuple(coefs), 3) for x in tbins]
plot!(tbins, v)
k += 1

test_bspline(
    knots,
    coefs,
    collect(LinRange(20.0, 40.0, 100))
)
