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

# Command line: julia -t N bin/tims_convert.jl <run.d> <out_dir> [--key value ...]
#               julia bin/tims_convert.jl expand <name.tdfs> <out.arrow>

const CLI_HELP = """
tims_convert.jl <run.d> <out_dir> [options]      convert a .d bundle
tims_convert.jl expand <name.tdfs> <out.arrow>   expand a tdfs into Pioneer's slice Arrow

Options (defaults in brackets):
  --stride-k0 X [0.0065]  --im-sigma-k0 X [0.004325]   slice spacing / IM kernel sigma in 1/K0, converted to scans
                             per run (stride rounded up)
  --stride K  --ms1-stride K  --im-sigma S  --ms1-im-sigma S   override in scans (default: derived as above)
  --kernel-extent E [3]  --no-sum-scale
  --mz-sigma-ns X [0.3125]    m/z kernel sigma in ns of flight time, converted to TOF bins per run
  --mz-sigma S                override in bins (default: derived as above; 2.5 on timsTOF Ultra, 1.56 on Pro)
  --centroid wmean|gauss|none [wmean]  --max-half H [max(4, 4*mz-sigma)]
  --min-scans N [1]
  --max-peaks N [1500, 0 = off]  --ms1-max-peaks N [0]   keep the N most intense centroids per slice
  --bin-scale K [1]  --int-scale X [1]  --zstd-level L [3]  --format tdfs|arrow|both [tdfs]
  --batch-frames N [16*threads, results in flight]  --frames a:b (frame rows)  --name NAME
"""

function parse_cli(args::Vector{String})
    length(args) >= 2 || (println(CLI_HELP); error("not enough arguments"))
    kw = Dict{Symbol, Any}(); name = nothing
    i = 3
    flags = Dict("--stride-k0" => (:stride_k0, Float64), "--im-sigma-k0" => (:im_sigma_k0, Float64),
                 "--im-sigma" => (:im_sigma, Float64), "--ms1-im-sigma" => (:ms1_im_sigma, Float64), "--kernel-extent" => (:kernel_extent, Float64),
                 "--stride" => (:stride, Int), "--ms1-stride" => (:ms1_stride, Int), "--mz-sigma" => (:mz_sigma, Float64), "--mz-sigma-ns" => (:mz_sigma_ns, Float64),
                 "--centroid" => (:centroid, Symbol), "--max-half" => (:max_half, Int),
                 "--min-scans" => (:min_scans, Int),
                 "--max-peaks" => (:max_peaks, Int), "--ms1-max-peaks" => (:ms1_max_peaks, Int),
                 "--bin-scale" => (:bin_scale, Int), "--int-scale" => (:int_scale, Float64), "--zstd-level" => (:zstd_level, Int),
                 "--format" => (:format, Symbol), "--batch-frames" => (:batch_frames, Int))
    while i <= length(args)
        a = args[i]
        if a == "--no-sum-scale"; kw[:sum_scale] = false; i += 1
        elseif a == "--frames"
            lo, hi = split(args[i + 1], ':'); kw[:frames] = collect(parse(Int, lo):parse(Int, hi)); i += 2
        elseif a == "--name"; name = args[i + 1]; i += 2
        elseif haskey(flags, a)
            k, T = flags[a]; v = args[i + 1]
            kw[k] = T === Symbol ? Symbol(v) : T === Int ? parse(Int, v) : parse(Float64, v); i += 2
        else
            println(CLI_HELP); error("unknown option $a")
        end
    end
    ConvertParams(; kw...), name
end

function main(args::Vector{String} = ARGS)
    if !isempty(args) && args[1] == "expand"
        length(args) == 3 || (println(CLI_HELP); error("expand takes <name.tdfs> <out.arrow>"))
        expand(args[2], args[3])
        return
    end
    p, name = parse_cli(args)
    # `--frames` filter needs the file to validate; convert does that
    if name === nothing
        convert_run(args[1], args[2]; params = p)
    else
        convert_run(args[1], args[2]; params = p, name = name)
    end
    nothing
end
