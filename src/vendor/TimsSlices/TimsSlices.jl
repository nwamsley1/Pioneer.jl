# ©2026 Washington University
#
# This file is part of TimsSlices.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

module TimsSlices

include("codec/zstd.jl")
include("codec/planes.jl")
include("codec/words.jl")
include("tdf/sqlite.jl")
include("tdf/block.jl")
include("tdf/reader.jl")
include("smooth/params.jl")
include("smooth/kernels.jl")
include("smooth/im.jl")
include("smooth/centroid.jl")
include("smooth/mz.jl")
include("smooth/window.jl")
include("codec/tdfs.jl")
include("arrow.jl")
include("convert.jl")
include("expand.jl")
include("cli.jl")

export open_tdf, FrameBuffer, read_frame!, valid_frames, windows, n_frames,
       FrameSlices, SliceBlock, BlockCodec, quantize!, encode_block!, decode_block!,
       ConvertParams, LevelParams, LevelSetup, SmoothScratch, level_params, smooth_frame!, smooth_window!, gauss_kernel,
       convert_run, expand, open_tdfs, read_frame_block!, read_slice!, SliceBuffer, encode_slice!, decode_slice!, output_name, TdfsFile, n_slices

end
