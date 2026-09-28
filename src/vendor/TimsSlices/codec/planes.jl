# Copyright (C) 2026 Nathan Wamsley
#
# This file is part of TimsSlices.jl
#
# TimsSlices.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# Byte transposition: n UInt32 words <-> four byte planes of length n (plane p holds byte p of every word).
# Bruker's layout, and ours. Both loops are plain SIMD loops (widen/shift/or, narrowing stores).
#
# Why: the words are mostly small numbers (bin deltas, intensities), so their high bytes are nearly all zero.
# Stored word by word, those zeros are interleaved with the busy low bytes; stored plane by plane, the high planes
# become long runs of zeros and near-zeros that zstd compresses to almost nothing.

"""
    untranspose!(words, planes, n)

`words[i] = planes[i] | planes[n+i] << 8 | planes[2n+i] << 16 | planes[3n+i] << 24` for i in 1:n.
`planes` must hold at least 4n bytes; `words` is resized to n.
"""
function untranspose!(words::Vector{UInt32}, planes::AbstractVector{UInt8}, n::Integer)
    length(planes) >= 4n || throw(ArgumentError("planes has $(length(planes)) bytes, need $(4n)"))
    length(words) < n && resize!(words, n)
    n == 0 && return words
    GC.@preserve planes words begin
        p = pointer(planes); w = pointer(words)
        @inbounds @simd for i in 0:n-1
            b0 = unsafe_load(p, i + 1); b1 = unsafe_load(p, n + i + 1)
            b2 = unsafe_load(p, 2n + i + 1); b3 = unsafe_load(p, 3n + i + 1)
            unsafe_store!(w, UInt32(b0) | (UInt32(b1) << 8) | (UInt32(b2) << 16) | (UInt32(b3) << 24), i + 1)
        end
    end
    words
end

"""
    transpose!(planes, words, n)

Inverse of `untranspose!`: writes the four byte planes of `words[1:n]` into `planes` (resized to 4n).
"""
function transpose!(planes::Vector{UInt8}, words::AbstractVector{UInt32}, n::Integer)
    length(words) >= n || throw(ArgumentError("words has $(length(words)) entries, need $n"))
    length(planes) < 4n && resize!(planes, 4n)
    n == 0 && return planes
    GC.@preserve planes words begin
        p = pointer(planes); w = pointer(words)
        for pl in 0:3
            sh = 8pl; base = pl * n
            @inbounds @simd for i in 1:n
                unsafe_store!(p, (unsafe_load(w, i) >> sh) % UInt8, base + i)
            end
        end
    end
    planes
end
