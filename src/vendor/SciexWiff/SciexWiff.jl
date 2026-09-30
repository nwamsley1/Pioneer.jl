"""
    SciexWiff

Vendored into Pioneer.jl from SciexWiff.jl (github.com/nwamsley1/SciexWiff.jl, commit e4d9097, v0.1.1).

Pure-Julia reader for SCIEX `.wiff` + `.wiff.scan` DIA (stepped SWATH) runs,
writing Pioneer's Arrow input format and the `.scxs` container. See `notes/format.md` and
`notes/scxs_format.md`.
"""
module SciexWiff

include("CFB.jl")
using .CFB

include("Idx.jl")
include("ScanFile.jl")
include("Method.jl")
include("Run.jl")
include("Centroid.jl")
include("Codec.jl")
include("Container.jl")
include("ArrowOut.jl")
include("Convert.jl")

export WiffRun, ScanBuffer, CentroidBuffer, CentroidParams, read_scan!, centroid!,
    ms_order, window, cycle_experiment, retention_time_min, bin_to_mz,
    ConvertParams, convert_run, open_scxs

end # module
