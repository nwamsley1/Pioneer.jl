# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

"""
    SpectrumDecoder

Pure-Julia reader for Orbitrap-family acquisition files: scan headers, centroided peak lists, profile points and
auxiliary records, and a writer for Pioneer's spectrum `.arrow` input.

    f = open_acquisition("run.raw")
    s = scans(f)                 # Vector{ScanInfo}
    mz, intensity = peaks(f, 1)  # centroided peak list of scan 1
    write_pioneer_arrow(f, "run.arrow")
"""
module SpectrumDecoder

using Mmap
using Arrow

export open_acquisition, scans, peaks, profile, features, traces, status_logs, scan_extras, instrument_status, read_auxiliary, method_text, method_windows, acquisition_mode, centroid!, Centroid, ScanInfo, Acquisition, write_pioneer_arrow, check_consistency, write_synthetic_acquisition

include("layout.jl")
include("reader.jl")
include("centroid.jl")
include("writer.jl")
include("checks.jl")
include("devices.jl")
include("method.jl")
include("acquisition.jl")
include("synthetic.jl")

end # module
