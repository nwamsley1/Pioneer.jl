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

# analysis.tdf: the tables the converter needs, read once into plain vectors.
using SQLite

struct FrameTable
    id::Vector{Int32}
    time::Vector{Float64}          # seconds
    msms_type::Vector{Int32}       # 0 = MS1, 9 = diaPASEF MS2
    num_scans::Vector{Int32}
    num_peaks::Vector{Int64}
    tims_id::Vector{Int64}         # byte offset of the block in analysis.tdf_bin
    mz_calibration::Vector{Int32}
    ramp_time::Vector{Float64}     # ms
end
Base.length(f::FrameTable) = length(f.id)

struct DiaWindow
    scan_begin::Int32   # inclusive, 0-based
    scan_end::Int32     # exclusive
    center::Float32
    width::Float32
    ce::Float32
end

struct DiaScheme
    groups::Vector{Vector{DiaWindow}}     # windows of window group g at groups[g] (sorted by scan_begin)
    frame_group::Vector{Int32}            # per frame row: window group (0 for non-DIA frames)
end

function read_global_metadata(db::SQLite.DB)
    meta = Dict{String, String}()
    for r in DBInterface.execute(db, "SELECT Key, Value FROM GlobalMetadata")
        meta[string(r[:Key])] = string(r[:Value])
    end
    meta
end

function read_frames(db::SQLite.DB)
    q = DBInterface.execute(db, "SELECT Id, Time, MsMsType, NumScans, NumPeaks, TimsId, MzCalibration, RampTime FROM Frames ORDER BY Id")
    f = FrameTable(Int32[], Float64[], Int32[], Int32[], Int64[], Int64[], Int32[], Float64[])
    for r in q
        push!(f.id, r[:Id]); push!(f.time, r[:Time]); push!(f.msms_type, r[:MsMsType]); push!(f.num_scans, r[:NumScans])
        push!(f.num_peaks, r[:NumPeaks]); push!(f.tims_id, r[:TimsId]); push!(f.mz_calibration, r[:MzCalibration])
        push!(f.ramp_time, r[:RampTime])
    end
    f
end

function read_dia_scheme(db::SQLite.DB, frames::FrameTable)
    groups = Dict{Int, Vector{DiaWindow}}()
    for r in DBInterface.execute(db, "SELECT WindowGroup, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth, CollisionEnergy FROM DiaFrameMsMsWindows ORDER BY WindowGroup, ScanNumBegin")
        push!(get!(groups, Int(r[:WindowGroup]), DiaWindow[]),
              DiaWindow(Int32(r[:ScanNumBegin]), Int32(r[:ScanNumEnd]), Float32(r[:IsolationMz]), Float32(r[:IsolationWidth]), Float32(r[:CollisionEnergy])))
    end
    ng = isempty(groups) ? 0 : maximum(keys(groups))
    gv = [get(groups, g, DiaWindow[]) for g in 1:ng]
    row_of = Dict{Int32, Int}(frames.id[i] => i for i in 1:length(frames))
    fg = zeros(Int32, length(frames))
    for r in DBInterface.execute(db, "SELECT Frame, WindowGroup FROM DiaFrameMsMsInfo")
        i = get(row_of, Int32(r[:Frame]), 0)
        i > 0 && (fg[i] = Int32(r[:WindowGroup]))
    end
    DiaScheme(gv, fg)
end

# --- calibration ------------------------------------------------------------------------------------
"√(m/z) = intercept + slope · bin"
struct LinearMzCal
    intercept::Float64
    slope::Float64
end
@inline bin_to_mz(c::LinearMzCal, bin) = (c.intercept + c.slope * bin)^2
@inline mz_to_bin(c::LinearMzCal, mz) = (sqrt(mz) - c.intercept) / c.slope

"""
Least squares on the CalibrationInfo reference peaks (exact on timsControl files).

In a TOF analyser the flight time grows with the square root of m/z, so `sqrt(m/z)` is linear in the TOF bin.
timsControl stores the calibrant peaks it fitted: their masses (`ReferencePeakMasses`) and measured flight times
(`MeasuredTimesOfFlight`, in ns). A flight time becomes a bin as `(t - DigitizerDelay) / DigitizerTimebase`; a
straight-line fit of `sqrt(mass)` on bin then gives the file's calibration. It reproduces Bruker's own m/z to
within 0-2 ppm on the files checked; the model built from the acquisition range alone (`MzAcqRange` and the
digitizer sample count) is off by -20..+7 ppm, which is why it is not used. Returns the line and each reference
peak's residual in ppm.
"""
function regressed_mz_cal(db::SQLite.DB, mzcal_id::Integer)
    function blob(key)
        r = first(DBInterface.execute(db, "SELECT Value FROM CalibrationInfo WHERE KeyPolarity='+' AND KeyName='$key'"))
        v = r[:Value]
        collect(reinterpret(Float64, v isa Vector{UInt8} ? v : Vector{UInt8}(v)))
    end
    masses = blob("ReferencePeakMasses"); t = blob("MeasuredTimesOfFlight")
    r = first(DBInterface.execute(db, "SELECT DigitizerTimebase, DigitizerDelay FROM MzCalibration WHERE Id=$mzcal_id"))
    x = (t .- r[:DigitizerDelay]) ./ r[:DigitizerTimebase]
    y = sqrt.(masses)
    xm = sum(x) / length(x); ym = sum(y) / length(y)
    b = sum((x .- xm) .* (y .- ym)) / sum((x .- xm) .^ 2)
    a = ym - b * xm
    resid_ppm = [((a + b * xi)^2 - m) / m * 1e6 for (xi, m) in zip(x, masses)]
    LinearMzCal(a, b), resid_ppm
end

"""
1/K0 = intercept + slope · scan (boundary model; the search refits its own line).

The TIMS ramp runs from high to low mobility, so scan 0 is `OneOverK0AcqRangeUpper` and the last scan
`OneOverK0AcqRangeLower`, linearly in between. This is the nominal ramp, not a fitted calibration (it was ~4% off
in slope on one timsTOF Pro file); Pioneer only uses the slope, to turn 1/K0 widths into scan counts.
"""
struct LinearImCal
    intercept::Float64
    slope::Float64
end
@inline scan_to_im(c::LinearImCal, scan) = c.intercept + c.slope * scan
function boundary_im_cal(meta::Dict{String, String}, frames::FrameTable)
    lo = parse(Float64, meta["OneOverK0AcqRangeLower"]); hi = parse(Float64, meta["OneOverK0AcqRangeUpper"])
    LinearImCal(hi, (lo - hi) / maximum(frames.num_scans))
end

"""
Collision energy ramp eV = intercept + slope · scan, fit to the window table's mid-scan values.

diaPASEF ramps the collision energy with mobility. The window table stores one energy per window, the ramp's
value at the window's middle scan; a line through those points recovers the ramp, so each slice can be given the
energy at its own scan rather than its window's average.
"""
struct CeRamp
    intercept::Float64
    slope::Float64
end
@inline ce_at(c::CeRamp, scan) = c.intercept + c.slope * scan
function fit_ce_ramp(scheme::DiaScheme)
    xs = Float64[]; ys = Float64[]
    for ws in scheme.groups, w in ws
        push!(xs, (w.scan_begin + w.scan_end) / 2); push!(ys, w.ce)
    end
    length(xs) < 2 && return CeRamp(isempty(ys) ? 0.0 : ys[1], 0.0)
    xm = sum(xs) / length(xs); ym = sum(ys) / length(ys)
    den = sum((xs .- xm) .^ 2)
    sl = den > 0 ? sum((xs .- xm) .* (ys .- ym)) / den : 0.0
    CeRamp(ym - sl * xm, sl)
end
