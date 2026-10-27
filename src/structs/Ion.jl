# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

abstract type Ion{T<:AbstractFloat} end

getMZ(I::Ion{T}) where {T<:AbstractFloat} = I.mz


import Base.<
import Base.>

<(y::Ion{T}, x::T) where {T<:AbstractFloat} = getMZ(y) < x
<(x::T, y::Ion{T}) where {T<:AbstractFloat} = <(y, x)
>(y::Ion{T}, x::T) where {T<:AbstractFloat} = getMZ(y) > x
>(x::T, y::Ion{T}) where {T<:AbstractFloat} = >(y, x)