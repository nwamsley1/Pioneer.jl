abstract type Ion{T<:AbstractFloat} end

getMZ(I::Ion{T}) where {T<:AbstractFloat} = I.mz


import Base.<
import Base.>

<(y::Ion{T}, x::T) where {T<:AbstractFloat} = getMZ(y) < x
<(x::T, y::Ion{T}) where {T<:AbstractFloat} = <(y, x)
>(y::Ion{T}, x::T) where {T<:AbstractFloat} = getMZ(y) > x
>(x::T, y::Ion{T}) where {T<:AbstractFloat} = >(y, x)