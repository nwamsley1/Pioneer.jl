# Copyright (C) 2026 Nathan Wamsley
# Vendored into Pioneer.jl from TimsSlices.jl (github.com/nwamsley1/TimsSlices.jl, commit b3e6d12).
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
