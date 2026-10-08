# Synthetic ZT Scan DIA runs, so the release build precompiles the scanning-quad (ZT) search path
# without shipping or downloading real ZT data. Included by snoop.jl.
#
# A run is written as `<name>.zt.scxs`, exactly as `convertSciex(...; zt_scan = true)` writes a real
# one: the same container, and the `acquisition_type = zt_scan_dia` metadata that puts the search in
# ZT mode. Its contents are simulated from the committed E. coli test library:
#
#   - target precursors elute as Gaussians in RT, placed by their library iRT;
#   - each cycle is one MS1 scan, then MS2 scans sweeping Q1 across the library's precursor range in
#     contiguous bins (the recorded centerMz / isolationWidthMz of a ZT run);
#   - a precursor's library fragments (spline intensities at the search NCE) appear in every bin
#     within the quadrupole's transmission, scaled by a triangle profile in Q1 offset, as on the
#     instrument; MS1 carries its first isotopes;
#   - plus low-level random noise peaks in every scan.
#
# Two runs (different abundances, a small RT shift, some precursors missing from each) so the
# cross-run stages (MBR, normalization, protein quantification) also run.

using Arrow, Random

const ZT_SYNTH_BIN_STEP = 1.02          # Q1 bin width = step (Da), as on the ZenoTOF 8600
const ZT_SYNTH_HALF_BASE = 6.0          # transmission triangle half-base (Da); measured 6.2-6.5
const ZT_SYNTH_CYCLE_S = 0.5            # one MS1 + one Q1 sweep (s)
const ZT_SYNTH_GRADIENT_MIN = 3.0       # run length (min)
const ZT_SYNTH_SIGMA_MIN = 0.05         # elution peak sigma (min)
const ZT_SYNTH_NOISE_PEAKS = 25
const ZT_SYNTH_CAL_A = 4.8e-4           # m/z = (cal_a * (bin / bin_scale / 40 - cal_b))^2, ~0.3 ppm bins
const ZT_SYNTH_BIN_SCALE = 4
const ZT_SYNTH_INT_SCALE = 10.0

"Library precursors to simulate: m/z, charge, iRT and their fragments (intensities at `nce`)."
function _zt_synth_precursors(lib_dir::AbstractString, n::Int, nce::Real, rng::AbstractRNG)
    t = Arrow.Table(joinpath(lib_dir, "precursors_table.arrow"))
    frags, ranges = Pioneer.load_detailed_frags_and_ranges(lib_dir)
    knots = Tuple(Float32.(Pioneer.deserialize_from_jls(joinpath(lib_dir, "spline_knots.jls"))))
    spline = Pioneer.prepare_spline_fractions(Float32(nce), knots)
    targets = shuffle!(rng, findall(!, t.is_decoy))[1:min(n, count(!, t.is_decoy))]
    precs = map(targets) do pid
        fr = frags[ranges[pid]:ranges[pid+1]-1]
        # Every library fragment, floored at a few percent of the base peak: the fragment index scores
        # the precursor's top fragments by library rank, which need not be the most intense at `nce`,
        # and tuning requires all of them to match.
        ints = max.(Float64[Pioneer.splevl_prepared(f.intensity, spline) for f in fr], 0.0)
        ints = max.(ints ./ max(maximum(ints), eps()), 0.05)
        (mz = Float64(t.mz[pid]), charge = Int(t.prec_charge[pid]), irt = Float64(t.irt[pid]),
         frag_mz = Float64[f.mz for f in fr], frag_int = ints)
    end
    return sort!(filter(p -> !isempty(p.frag_mz), precs); by = p -> p.mz)
end

_zt_synth_triangle(delta) = max(0.0, 1.0 - abs(delta) / ZT_SYNTH_HALF_BASE)

# Stored bin of an m/z under the synthetic calibration (cal_b = 0).
_zt_synth_bin(mz) = round(UInt32, ZT_SYNTH_BIN_SCALE * 40 * sqrt(mz) / ZT_SYNTH_CAL_A)

"Quantize peaks to stored bins/intensities (sorted, equal bins merged, zero intensities dropped)."
function _zt_synth_quantize(mz::Vector{Float64}, it::Vector{Float64})
    order = sortperm(mz)
    bins = UInt32[]; ints = UInt32[]
    for k in order
        b = _zt_synth_bin(mz[k]); v = round(UInt64, ZT_SYNTH_INT_SCALE * it[k])
        if !isempty(bins) && bins[end] == b
            ints[end] = UInt32(min(UInt64(ints[end]) + v, typemax(UInt32)))
        elseif v > 0
            push!(bins, b); push!(ints, UInt32(min(v, typemax(UInt32))))
        end
    end
    return bins, ints
end

"""
    write_synthetic_zt_run(path, precs; rng, abundance, rt_shift_min, present)

One synthetic ZT run as a `.scxs` directory at `path`.
"""
function write_synthetic_zt_run(path::AbstractString, precs; rng::AbstractRNG,
                                abundance::Vector{Float64}, rt_shift_min::Float64, present::BitVector)
    irt_lo, irt_hi = extrema(p.irt for p in precs)
    rt_of(p) = 0.25 + (p.irt - irt_lo) / (irt_hi - irt_lo) * (ZT_SYNTH_GRADIENT_MIN - 0.5) + rt_shift_min
    rts = [rt_of(p) for p in precs]
    mzs = [p.mz for p in precs]
    lo = floor(minimum(mzs)) - 2; hi = ceil(maximum(mzs)) + 2
    centers = collect(lo + ZT_SYNTH_BIN_STEP / 2 : ZT_SYNTH_BIN_STEP : hi)
    n_cycles = floor(Int, ZT_SYNTH_GRADIENT_MIN * 60 / ZT_SYNTH_CYCLE_S)
    dwell_min = ZT_SYNTH_CYCLE_S / (length(centers) + 1) / 60

    meta = Dict{String, Any}(
        "source" => basename(path), "converter" => "Pioneer synthetic ZT (src/build/synthetic_zt.jl)",
        "bin_scale" => ZT_SYNTH_BIN_SCALE, "int_scale" => ZT_SYNTH_INT_SCALE,
        "acquisition_method" => "synthetic ZT Scan DIA", "acquisition_type_source" => "user",
        "acquisition_type" => "zt_scan_dia",
        "q1_bin_width_mz" => string(ZT_SYNTH_BIN_STEP), "q1_bin_step_mz" => string(ZT_SYNTH_BIN_STEP),
        "q1_bins_per_cycle" => string(length(centers)),
        "q1_first_bin_center_mz" => string(centers[1]), "q1_last_bin_center_mz" => string(centers[end]),
        "q1_bin_dwell_ms" => string(dwell_min * 60_000),
        "q1_scan_rate_mz_per_s" => string(ZT_SYNTH_BIN_STEP / (dwell_min * 60)))
    rm(path; force = true, recursive = true)
    writer = Pioneer.SciexWiff.ScxsWriter(path, meta)
    codec = Pioneer.SciexWiff.BlockCodec()
    record = 0
    function emit!(cycle, experiment, ms_order, rt, center, mz, it)
        bins, ints = _zt_synth_quantize(mz, it)
        n = length(bins)
        z = codec.zbuf[1:Pioneer.SciexWiff.encode_block!(codec, bins, ints, n, 3)]
        k = n == 0 ? 0 : argmax(ints)
        record += 1
        row = (record = Int32(record), cycle = Int32(cycle), experiment = Int16(experiment),
               ms_order = UInt8(ms_order), retention_time = Float32(rt), low_mz = 100f0, high_mz = 1800f0,
               tic = Float32(sum(ints; init = UInt32(0)) / ZT_SYNTH_INT_SCALE),
               base_peak_mz = n == 0 ? NaN32 : Float32((ZT_SYNTH_CAL_A * bins[k] / ZT_SYNTH_BIN_SCALE / 40)^2),
               base_peak_intensity = n == 0 ? 0f0 : Float32(ints[k] / ZT_SYNTH_INT_SCALE),
               center_mz = ms_order == 1 ? NaN32 : Float32(center),
               isolation_width = ms_order == 1 ? NaN32 : Float32(ZT_SYNTH_BIN_STEP),
               cal_a = ZT_SYNTH_CAL_A, cal_b = 0.0, n_peaks = Int32(n))
        Pioneer.SciexWiff.write_scan!(writer, row, z)
    end

    reach = 4 * ZT_SYNTH_SIGMA_MIN
    mz = Float64[]; it = Float64[]
    for c in 1:n_cycles
        t0 = (c - 1) * ZT_SYNTH_CYCLE_S / 60
        eluting = [i for i in eachindex(precs) if present[i] && abs(rts[i] - t0) < reach]
        # MS1: first three isotopes of every eluting precursor
        empty!(mz); empty!(it)
        for i in eluting
            p = precs[i]; e = abundance[i] * exp(-((t0 - rts[i]) / ZT_SYNTH_SIGMA_MIN)^2 / 2)
            for (iso, r) in ((0, 1.0), (1, 0.6), (2, 0.25))
                push!(mz, p.mz + iso * 1.00336 / p.charge); push!(it, e * r)
            end
        end
        _zt_synth_noise!(mz, it, rng)
        emit!(c, 1, 1, t0, NaN, mz, it)
        # MS2: the Q1 sweep
        for (b, center) in enumerate(centers)
            t = t0 + b * dwell_min
            empty!(mz); empty!(it)
            for i in eluting
                p = precs[i]
                w = _zt_synth_triangle(p.mz - center)
                w > 0 || continue
                e = abundance[i] * w * exp(-((t - rts[i]) / ZT_SYNTH_SIGMA_MIN)^2 / 2)
                for (fm, fi) in zip(p.frag_mz, p.frag_int)
                    push!(mz, fm * (1 + 2e-6 * randn(rng))); push!(it, e * fi)
                end
            end
            _zt_synth_noise!(mz, it, rng)
            emit!(c, b + 1, 2, t, center, mz, it)
        end
    end
    close(writer)
    return path
end

function _zt_synth_noise!(mz, it, rng)
    for _ in 1:ZT_SYNTH_NOISE_PEAKS
        push!(mz, 150 + 1450 * rand(rng)); push!(it, 50 + 200 * rand(rng))
    end
end

"""
    generate_synthetic_zt_runs(lib_dir, out_dir; n_precursors = 1000, n_runs = 2, nce = 26, seed = 1)
        -> Vector{String}

Write `n_runs` synthetic ZT `.scxs` runs into `out_dir` (replacing it) and return their paths.
"""
function generate_synthetic_zt_runs(lib_dir::AbstractString, out_dir::AbstractString;
                                    n_precursors::Int = 1000, n_runs::Int = 2, nce::Real = 26, seed::Int = 1)
    rng = MersenneTwister(seed)
    precs = _zt_synth_precursors(lib_dir, n_precursors, nce, rng)
    base = [exp(log(2e4) + 1.0 * randn(rng)) for _ in precs]
    rm(out_dir; force = true, recursive = true); mkpath(out_dir)
    map(1:n_runs) do r
        abundance = base .* exp.(0.2 .* randn(rng, length(precs)))
        present = BitVector(rand(rng, length(precs)) .> 0.1)
        write_synthetic_zt_run(joinpath(out_dir, "synthetic_zt_$(r).zt.scxs"), precs; rng = rng,
                               abundance = abundance, rt_shift_min = 0.02 * (r - 1), present = present)
    end
end

"""
    check_zt_search_ran(results_dir)

Fail unless the search's debug log shows every ZT stage taking its ZT path: scanning-quad mode on,
a successful transmission-triangle fit in quad tuning, the meta-scan main search and the meta-scan
chromatogram collapse. A search that finishes after silently falling back would otherwise leave
those code paths uncompiled.
"""
function check_zt_search_ran(results_dir::AbstractString)
    log = read(joinpath(results_dir, "pioneer_search_debug.log"), String)
    occursin("triangle fit failed", log) &&
        error("ZT quad tuning fell back (triangle fit failed); the fit was not exercised")
    for (what, pattern) in (("ZT mode", r"Scanning-quad \(ZT\) \[[^\]]+\]: ON"),
                            ("quad tuning triangle fit", r"ZT quad tuning \[file \d+\]: h="),
                            ("ZT main search", r"ZT main search: "),
                            ("ZT chromatogram collapse", r"ZT chromatogram collapse "))
        occursin(pattern, log) || error("ZT stage did not run on the synthetic data: $what")
    end
    return nothing
end
