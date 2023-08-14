using Compose
using Color

function venn(xs::Union{AbstractArray,Set}...; proportional::Bool = true, labels=Union{Bool,Vector{String}}, colors=Union{Bool,Vector{ColorValue},Vector{AlphaColorValue}})
    n = length(xs)

    if eltype(xs) != Set
        cols = map(Set, xs)
    else
        cols = xs
    end

    if proportional
        error("Not implemented yet!")
        #sizes = 0
    end

    return p
end
