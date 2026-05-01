# Retention time index for chromatogram extraction.
#
# Maps RT ranges to precursors, enabling efficient lookup of which precursors
# to extract chromatograms for in a given RT window. Used by
# IntegrateChromatogramsSearch to avoid scanning all precursors for every scan.
#
# Structure: Vector of RT bins, each containing a sorted list of
# (precursor_id, precursor_mz) tuples covering that RT range.

"""
    rtIndexBin{T, U}

A single RT bin containing precursors whose iRT falls within [lb, ub].

# Fields
- `lb::T`: Lower bound of the RT bin (iRT of first precursor)
- `ub::T`: Upper bound of the RT bin (iRT of last precursor)
- `prec::Vector{Tuple{UInt32, U}}`: Precursors in this bin as (precursor_id, precursor_mz),
   sorted by m/z for efficient isolation window matching
"""
struct rtIndexBin{T, U<:AbstractFloat}
    lb::T
    ub::T
    prec::Vector{Tuple{UInt32, U}}
end

getLow(r::rtIndexBin) = r.lb
getHigh(r::rtIndexBin) = r.ub

"""
    retentionTimeIndex{T, U}

Collection of RT bins covering the full gradient. Bins are contiguous and
sorted by RT, with precursors sorted by m/z within each bin.

Used by IntegrateChromatogramsSearch to find precursors near a given scan's
RT, then filter by m/z against the quadrupole isolation window.
"""
struct retentionTimeIndex{T, U<:AbstractFloat}
    rt_bins::Vector{rtIndexBin{T, U}}
end

retentionTimeIndex(T::DataType, U::DataType) = retentionTimeIndex(Vector{rtIndexBin{T, U}}())
