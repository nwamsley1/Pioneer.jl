mutable struct ArrayDict{I<:Unsigned,C<:Real}
    keys::Vector{I}
    vals::Vector{C}
    size::Int64
    function ArrayDict(I::DataType, C::DataType, size::Int) #where {I,C<:Unsigned}
        new{I, C}(zeros(I, size), zeros(C, size), 0)
    end
end

function update!(c::ArrayDict{I,C}, key::I, val::C) where {I,C<:Unsigned}
    c.size += 1
    c.vals[key] = val
    c.keys[c.size] = key
end

function reset!(c::ArrayDict{I,C}) where {I,C<:Unsigned}
    @turbo for i in range(1, c.size)
        c.vals[c.keys[i]] = zero(I)
        c.keys[i] = zero(C)
    end
    c.size = 0
end

function Base.getindex(c::ArrayDict{I,C}, i::Ti) where {I,C<:Unsigned, Ti<:Integer}
    c.vals[i]
end