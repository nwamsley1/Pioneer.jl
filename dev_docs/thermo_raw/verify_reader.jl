# End-to-end check of ThermoRaw (offsets from the file's own headers) against a RawFileReader-converted .arrow.
include(joinpath(@__DIR__, "ThermoRaw.jl")); using .ThermoRaw, Arrow, Tables
raw = length(ARGS) >= 1 ? ARGS[1] : expanduser("~/BrukerTims/rawformat/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.raw")
arw = length(ARGS) >= 2 ? ARGS[2] : expanduser("~/BrukerTims/tuning_ab/data/OlsenExploris500ng/20220909_EXPL8_Evo5_ZY_MixedSpecies_500ng_E5H50Y45_30SPD_DIA_1.arrow")
@time f = ThermoRaw.open_raw(raw)
t = Arrow.Table(arw)
println("version $(f.version); scans $(f.first_scan)-$(f.last_scan) ($(length(f.scans))) vs arrow $(length(t.msOrder))")
bad = Dict(:rt => 0, :tic => 0, :bp => 0, :range => 0, :ptype => 0, :order => 0, :center => 0, :width => 0, :ce => 0, :npk => 0, :peaks => 0)
@time for (i, s) in enumerate(f.scans)
    Float32(s.rt) == t.retentionTime[i] || (bad[:rt] += 1)
    Float32(s.tic) == t.TIC[i] || (bad[:tic] += 1)
    (Float32(s.base_peak_mz) == t.basePeakMz[i] && Float32(s.base_peak_intensity) == t.basePeakIntensity[i]) || (bad[:bp] += 1)
    (Float32(s.low_mz) == t.lowMz[i] && Float32(s.high_mz) == t.highMz[i]) || (bad[:range] += 1)
    s.packet_type == t.packetType[i] || (bad[:ptype] += 1)
    s.ms_order == t.msOrder[i] || (bad[:order] += 1)
    if s.ms_order > 1
        Float32(s.center_mz) == t.centerMz[i] || (bad[:center] += 1)
        Float32(s.isolation_width) == t.isolationWidthMz[i] || (bad[:width] += 1)
        Float32(s.collision_energy) == t.collisionEnergyField[i] || (bad[:ce] += 1)
    end
    mz, it = ThermoRaw.centroids(f, i)
    amz = collect(skipmissing(t.mz_array[i])); ait = collect(skipmissing(t.intensity_array[i]))
    length(mz) == length(amz) || (bad[:npk] += 1; continue)
    (mz == amz && it == ait) || (bad[:peaks] += 1)
end
println("mismatches: ", bad)
