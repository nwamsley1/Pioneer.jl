# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
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

"""
    setMatch!(matches, transition, mass, intensity, peak_ind, i)

Adds a FragmentMatch to `matches` at index `i+1`.
Grows the vector by `block_size` if needed.
"""
function setMatch!(matches::Vector{M},
                    transition::DetailedFrag{Float32},
                    mass::Float32,
                    intensity::Float32,
                    peak_ind::Int64,
                    i::Int64;
                    block_size = 10000) where {M<:MatchIon{Float32}}
    i += 1
    #Grow pre-allocated placeholder array if needed
    if i > length(matches)
        append!(matches, [FragmentMatch{Float32}() for _ in range(1, block_size)])
    end

    ion_type = zero(UInt8)
    if transition.is_y
        ion_type = UInt8(2)
    elseif transition.is_p
        ion_type = UInt8(3)
    elseif transition.is_b
        ion_type = UInt8(1)
    end
    matches[i] = FragmentMatch(
                                Float32(getIntensity(transition)),
                                 intensity,
                                 getMZ(transition),
                                 mass,
                                 UInt32(peak_ind),
                                 getPrecID(transition),
                                 getIonPosition(transition),
                                 getFragCharge(transition),
                                 ion_type,
                                 transition.is_isotope,
                                 getRank(transition)
                                 )
    return i
end

function setMatch!(matches::Vector{M},
                    precursor_isotope::Isotope{Float32},
                    mass::Float32,
                    intensity::Float32,
                    peak_ind::Int64,
                    i::Int64;
                    block_size = 10000) where {M<:MatchIon{Float32}}
    i += 1
    #Grow pre-allocated placeholder array if needed
    if i > length(matches)
        append!(matches, [PrecursorMatch{Float32}() for _ in range(1, block_size)])
    end
    matches[i] = PrecursorMatch(
                                Float32(getIntensity(precursor_isotope)),
                                 intensity,
                                 getMZ(precursor_isotope),
                                 mass,
                                 getIsoIdx(precursor_isotope),
                                 peak_ind,
                                 getPrecID(precursor_isotope)
                                 )
    return i
end

"""
    setUnmatched!(unmatched, ion, i)

Writes a lightweight UnmatchedIon to `unmatched` at index `i+1`.
Only stores the three fields needed by buildDesignMatrix!.
"""
function setUnmatched!(unmatched::Vector{UnmatchedIon},
                       ion::DetailedFrag{Float32},
                       i::Int64;
                       block_size = 10000)
    i += 1
    if i > length(unmatched)
        append!(unmatched, [UnmatchedIon() for _ in range(1, block_size)])
    end
    unmatched[i] = UnmatchedIon(getPrecID(ion), Float32(getIntensity(ion)), ion.is_isotope)
    return i
end

function setUnmatched!(unmatched::Vector{UnmatchedIon},
                       ion::Isotope{Float32},
                       i::Int64;
                       block_size = 10000)
    i += 1
    if i > length(unmatched)
        append!(unmatched, [UnmatchedIon() for _ in range(1, block_size)])
    end
    unmatched[i] = UnmatchedIon(getPrecID(ion), Float32(getIntensity(ion)), getIsoIdx(ion) > 0)
    return i
end

"""
    setNearest!(matches, unmatched, masses, intensities, Ion, ...)

Finds the nearest peak within tolerance and writes a match (or miss).
"""
function setNearest!(matches::Vector{M},
                    unmatched::Vector{UnmatchedIon},
                    masses::AbstractArray{Union{Missing, Float32}},
                    intensities::AbstractArray{Union{Missing, Float32}},
                    Ion::DetailedFrag{Float32},
                    mass_err_model::MassErrorModel,
                    low_theoretical_mz::Float32,
                    high_theoretical_mz::Float32,
                    corrected_empirical_mz::Float32,
                    peak_idx::Int64,
                    matched_idx::Int64,
                    unmatched_idx::Int64
                    )::Tuple{Int64,Int64} where {M<:MatchIon{Float32}}

    smallest_diff = typemax(Float32)
    #Have already verified that the current peak is within the mass tolerance
    best_mz, best_peak, i = corrected_empirical_mz, peak_idx, one(Int64)
    theoretical_mz = getMZ(Ion)

    #Iterate through peaks in  `masses` until a peak is encountered that is
    #greater in m/z than the upper bound of the theoretical ion mass tolerance
    #or until there are no more masses to check. Keep track of the best peak/transition match.

        if peak_idx + 1 <= length(masses)
            @inbounds @fastmath begin
                corrected_empirical_mz = getCorrectedMz(mass_err_model, masses[peak_idx + i])
                while (corrected_empirical_mz <= high_theoretical_mz)
                    if (corrected_empirical_mz >= low_theoretical_mz) #
                        mz_diff = abs(corrected_empirical_mz-theoretical_mz)
                        if mz_diff < smallest_diff
                            smallest_diff = mz_diff
                            best_peak = peak_idx + i
                            best_mz = corrected_empirical_mz
                        end
                    end
                    i+=1
                    if (peak_idx + 1 + i > length(masses)) break end
                end
            end
        end

        if best_peak > 0
            matched_idx = setMatch!(matches,
                                    Ion,
                                    best_mz,
                                    intensities[best_peak],
                                    best_peak,
                                    matched_idx);
        else
            unmatched_idx = setUnmatched!(unmatched,
                        Ion,
                        unmatched_idx);
        end
    return matched_idx, unmatched_idx
end
function setNearest!(matches::Vector{M},
                    unmatched::Vector{UnmatchedIon},
                    masses::AbstractArray{Union{Missing, Float32}},
                    intensities::AbstractArray{Union{Missing, Float32}},
                    Ion::Isotope{Float32},
                    mass_err_model::MassErrorModel,
                    low_theoretical_mz::Float32,
                    high_theoretical_mz::Float32,
                    corrected_empirical_mz::Float32,
                    peak_idx::Int64,
                    matched_idx::Int64,
                    unmatched_idx::Int64
                    )::Tuple{Int64,Int64} where {M<:MatchIon{Float32}}

    smallest_diff = typemax(Float32)
    #Have already verified that the current peak is within the mass tolerance
    best_mz, best_peak, i = corrected_empirical_mz, peak_idx, one(Int64)
    theoretical_mz = getMZ(Ion)

    #Iterate through peaks in  `masses` until a peak is encountered that is
    #greater in m/z than the upper bound of the theoretical ion mass tolerance
    #or until there are no more masses to check. Keep track of the best peak/transition match.

        if peak_idx + 1 <= length(masses)
            @inbounds @fastmath begin
                corrected_empirical_mz = getCorrectedMz(mass_err_model, masses[peak_idx + i])
                while (corrected_empirical_mz <= high_theoretical_mz)
                    if (corrected_empirical_mz >= low_theoretical_mz) #
                        mz_diff = abs(corrected_empirical_mz-theoretical_mz)
                        if mz_diff < smallest_diff
                            smallest_diff = mz_diff
                            best_peak = peak_idx + i
                            best_mz = corrected_empirical_mz
                        end
                    end
                    i+=1
                    if (peak_idx + 1 + i > length(masses)) break end
                end
            end
        end

        if best_peak > 0
            matched_idx = setMatch!(matches,
                                    Ion,
                                    best_mz,
                                    intensities[best_peak],
                                    best_peak,
                                    matched_idx);
        else
            unmatched_idx = setUnmatched!(unmatched,
                        Ion,
                        unmatched_idx);
        end
    return matched_idx, unmatched_idx
end
"""
    matchPeaks!(matches, unmatched, mz_index, templates, max_ions_idx, masses, intensities, mass_err_model, high_mass)

Sidecar variant: iterates mz-sorted MzSortEntry array for comparisons (8B sequential reads),
resolves full DetailedFrag template only on match/miss writes. Same O(T+P) algorithm.
"""
function matchPeaks!(matches::Vector{M},
                    unmatched::Vector{UnmatchedIon},
                    mz_index::Vector{MzSortEntry},
                    templates::Vector{DetailedFrag{Float32}},
                    max_ions_idx::Int64,
                    masses::AbstractArray{Union{Missing, Float32}},
                    intensities::AbstractArray{Union{Missing, Float32}},
                    mass_err_model::MassErrorModel,
                    high_mass::Float32
                    ) where {M<:MatchIon{Float32}}

    peak_idx, ion_idx, matched_idx, unmatched_idx = 1, 1, 0, 0

    if (max_ions_idx<1) | (iszero(length(masses)))
        return matched_idx, unmatched_idx
    end

    corrected_empirical_mz = getCorrectedMz(mass_err_model, masses[peak_idx])
    low, high = getMzBounds(mass_err_model, mz_index[ion_idx].mz)

    @inbounds @fastmath while (peak_idx <= length(masses)) & (ion_idx <= max_ions_idx)
        if (corrected_empirical_mz>=low)
            if (corrected_empirical_mz<=high)
                matched_idx, unmatched_idx = setNearest!(matches,
                                            unmatched,
                                            masses,
                                            intensities,
                                            templates[mz_index[ion_idx].idx],
                                            mass_err_model,
                                            low,
                                            high,
                                            corrected_empirical_mz,
                                            peak_idx,
                                            matched_idx,
                                            unmatched_idx)

                ion_idx += 1
                if ion_idx > max_ions_idx
                    return matched_idx, unmatched_idx
                end

                low, high = getMzBounds(mass_err_model, mz_index[ion_idx].mz)
                continue
            end
            unmatched_idx = setUnmatched!(unmatched,
                                        templates[mz_index[ion_idx].idx],
                                        unmatched_idx);
            ion_idx += 1
            low, high = getMzBounds(mass_err_model, mz_index[ion_idx].mz)

            if ion_idx > max_ions_idx
                return matched_idx, unmatched_idx
            end
            continue
        end

        peak_idx += 1
        if peak_idx <= length(masses)
            corrected_empirical_mz = getCorrectedMz(mass_err_model, masses[peak_idx])
        end
    end

    @inbounds @fastmath while ion_idx <= max_ions_idx
        if mz_index[ion_idx].mz > high_mass
            break
        end
        unmatched_idx = setUnmatched!(unmatched,
                                    templates[mz_index[ion_idx].idx],
                                    unmatched_idx);
        ion_idx += 1
    end

    return matched_idx, unmatched_idx
end

"""
    matchPeaks!(matches, unmatched, ions, max_ions_idx, masses, intensities, mass_err_model, high_mass)

Single-pass O(T+P) algorithm matching theoretical ions to empirical peaks.
"""
function matchPeaks!(matches::Vector{M}, #Pre-allocated container for Matched Ions
                    unmatched::Vector{UnmatchedIon}, #Pre-allocated container for Unmatched Ions
                    ions::Vector{I}, #Library Ion Templates to match to the spectrum
                    max_ions_idx::Int64, #Maximum element in Ions to consider
                    masses::AbstractArray{Union{Missing, Float32}},
                    intensities::AbstractArray{Union{Missing, Float32}},
                    mass_err_model::MassErrorModel,
                    high_mass::Float32
                    ) where {I<:LibraryIon{Float32},M<:MatchIon{Float32}}

    #matched_idx is a running count of the number of transitions that have been matched to a peak
    #unmatched_idx is a running count of the number of transitions that did not match any peak
    #peak_idx is the index of the peak in the spectrum (masses[peak_idx])
    #ion_idx is the index of the theoretical ion being matched to the spectrum (ions[ion])
    peak_idx, ion_idx, matched_idx, unmatched_idx = 1, 1, 0, 0

    if (max_ions_idx<1) | (iszero(length(masses)))
        return matched_idx, unmatched_idx
    end

    #Corrected m/z of the empirical peak
    corrected_empirical_mz = getCorrectedMz(mass_err_model, masses[peak_idx])
    #Mass tolerance of the theoretical ion
    low, high = getMzBounds(mass_err_model, getMZ(ions[ion_idx]))

    @inbounds @fastmath while (peak_idx <= length(masses)) & (ion_idx <= max_ions_idx)
        if (corrected_empirical_mz>=low) #Is the peak within the tolerance of the theoretical ion m/z?
            if (corrected_empirical_mz<=high) #Empirical peak is within mass tolerance
                #Set nearest matching peak if there are any.
                matched_idx, unmatched_idx = setNearest!(matches,
                                            unmatched,
                                            masses,
                                            intensities,
                                            ions[ion_idx],
                                            mass_err_model,
                                            low,
                                            high,
                                            corrected_empirical_mz,
                                            peak_idx,
                                            matched_idx,
                                            unmatched_idx)

                #Progress to the next theoretical ion
                #Do not advance the `peak_idx` counter
                #because multiple theoretical ions can match
                #to a single empirical peak, but not vice-versa
                ion_idx += 1
                if ion_idx > max_ions_idx
                    return matched_idx, unmatched_idx
                end

                #Mass tolerance of the new theoretical ion
                low, high = getMzBounds(mass_err_model, getMZ(ions[ion_idx]))
                continue
            end
            #Progress to next theoretical ion
            unmatched_idx = setUnmatched!(unmatched,
                                        ions[ion_idx],
                                        unmatched_idx);
            ion_idx += 1
            low, high = getMzBounds(mass_err_model, getMZ(ions[ion_idx]))

            if ion_idx > max_ions_idx
                return matched_idx, unmatched_idx
            end
            continue
        end

        peak_idx += 1 #Progress to next peak
        if peak_idx <= length(masses)
            corrected_empirical_mz = getCorrectedMz(mass_err_model, masses[peak_idx])
        end
    end

    #Remaining templates with higher m/z than the highest
    #m/z peak in the spectrum are written to 'unmatched'
    @inbounds @fastmath while ion_idx <= max_ions_idx
        if getMZ(ions[ion_idx]) > high_mass
            break
        end
        unmatched_idx = setUnmatched!(unmatched,
                                    ions[ion_idx],
                                    unmatched_idx);
        ion_idx += 1
    end

    return matched_idx, unmatched_idx
end
