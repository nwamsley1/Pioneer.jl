isoxml = readxml("./data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml")
isosplines = root(isoxml)

for node in eachelement(isosplines)
    if hasnode(node)
        if node.name == "model"
        for subnode in eachelement(node)
            println(node.name)
            println(node.content)
        end
        break
        end
    end
end
collect(eachelement(isosplines))[100]

xdoc = parse_file("./data/IsotopeSplines/IsotopeSplines_10kDa_21isotopes-1.xml")

i = 0
knots = nothing
coefficients = nothing
for c in root(xdoc)["model"]
    #println(attributes_dict(c))
    #println(attributes_dict(c["knots"][1]))
    #println(DecodeCoefficients(content(c["knots"][1])))
    knots = (DecodeCoefficients(content(c["knots"][1])))
    ##println(attributes_dict(c["coefficients"][1]))
    #println(DecodeCoefficients(content(c["coefficients"][1])))
    coefficients = DecodeCoefficients(content(c["coefficients"][1]))

end

function DecodeCoefficients(encoded::String)
    return reinterpret(Float64, Base64.base64decode(encoded))
end

struct PolynomialSpline{T<:Real}
    polynomials::Vector{Polynomial{T, :x}}
    knots::Vector{T}
end

function buildPolynomials(coefficients::Vector{T}, order::I) where {T<:Real, I<:Integer}
    order += 1
    n_knots = length(coefficients)÷order
    polynomials = Vector{Polynomial{T, :x}}()
    for n in range(1, n_knots)
        start = ((n - 1)*order + 1)
        stop = (n)*order
        push!(polynomials, Polynomial(coefficients[start:stop], :x))
    end
    return polynomials
end

cubic_spline = PolynomialSpline(
    buildPolynomials(collect(coefficients), 3),
    collect(knots)
)