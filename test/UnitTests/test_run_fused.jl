# Unit tests for run_fused! orchestration.
#
# The five extracted helpers (passes_irt_filter, quad_window_with_iso_bounds,
# passes_prec_mz_filter, iso_mz_for, in_frag_mz_window, conservative_half_width,
# match_window, compute_ppm_err, push_match!, push_miss!) have their own tests
# in test_fused_prec_filters.jl. This file targets the **orchestration**:
# loop control flow, branch coverage, column lifecycle, scratch accounting,
# Hs construction, and id_to_col writes.
#
# Run: julia --project=. test/UnitTests/test_run_fused.jl

using Test
using Pioneer
using Pioneer: parseIsoXML

const FusedStandard          = Pioneer.FusedStandard
const FusedRTIndexed         = Pioneer.FusedRTIndexed
const FusedQuadEst           = Pioneer.FusedQuadEst
const FullPrecCapture        = Pioneer.FullPrecCapture
const PartialPrecCapture     = Pioneer.PartialPrecCapture
const SparseArrayFused       = Pioneer.SparseArrayFused
const FusedScratch           = Pioneer.FusedScratch
const MainUnscoredPSM     = Pioneer.MainUnscoredPSM
const TuningUnscoredPSM      = Pioneer.TuningUnscoredPSM
const DensePrecMap           = Pioneer.DensePrecMap
const StandardFragmentLookup = Pioneer.StandardFragmentLookup
const CompactFrag            = Pioneer.CompactFrag
const PiecewiseNceModel      = Pioneer.PiecewiseNceModel
const SquareQuadFunction     = Pioneer.SquareQuadFunction
const SimpleMassErrorModel   = Pioneer.SimpleMassErrorModel
const run_fused!             = Pioneer.run_fused!

# Load the real isotope splines once (small XML, fast parse).
const ISO_XML = joinpath(dirname(dirname(@__DIR__)), "assets",
                         "IsotopeSplines_10kDa_21isotopes.xml")
const ISO_SPLINES = parseIsoXML(ISO_XML)

# -----------------------------------------------------------------------------
# Fixture builder.
# -----------------------------------------------------------------------------
# Constructs all 25+ args run_fused! needs. Returns a NamedTuple so individual
# tests can override fields and re-call run_fused!.
#
# Defaults: 1 precursor (idx=1), 2 fragments at known m/z, 1 scan with peaks
# placed at the fragment m/z values exactly so they should both match.
#
# Keyword args let each test override what it needs.
function make_fused_fixture(;
        # Precursor properties
        n_precursors::Int = 1,
        prec_mzs::Vector{Float32} = Float32[500.0],
        prec_charges::Vector{UInt8} = UInt8[2],
        prec_sulfur_counts::Vector{UInt8} = UInt8[0],
        prec_irts::Vector{Float32} = Float32[50.0],
        # Fragment library: vector of (prec_id, mz, intensity, rank)
        frags::Vector = [(UInt32(1), 200.0f0, Float32(1000), UInt8(0)),
                         (UInt32(1), 300.0f0, Float32( 800), UInt8(0))],
        # Per-precursor [start_idx, end_idx) into frags
        prec_frag_ranges::Vector{UInt64} = UInt64[1, 3],
        # Scan peaks: Vector{Float32} of peak m/z (sorted ascending)
        peak_mz::Vector{Float32} = Float32[200.0, 300.0],
        peak_int::AbstractVector = Float32[1500.0, 1200.0],
        # Quad window (defaults wide enough to admit prec_mz=500)
        quad_min::Float32 = 480.0f0,
        quad_max::Float32 = 520.0f0,
        # Mass error tolerance (ppm, symmetric). Either a scalar (applied to
        # every peak) or a per-peak Vector{Float32} the same length as peak_mz
        # — useful for stressing intensity-dependent tolerance, where high-
        # intensity peaks have tight windows and low-intensity peaks have wide
        # ones (mimics IntensityMassErrorModel's getCorrectedMzAndBounds).
        ppm_tol::Union{Float32, Vector{Float32}} = 20.0f0,
        # Scan / iRT
        scan_irt::Float32 = 50.0f0,
        irt_tol::Float32 = 5.0f0,
        # Other knobs
        n_frag_isotopes::Int64 = 1,
        isotope_err_bounds::Tuple{Int,Int} = (0, 0),
        kind = FusedStandard(FullPrecCapture(), UInt8(7)),
        prec_range::UnitRange{Int64} = 1:1,
    )

    # Hs: design matrix (constructor takes a single capacity Integer)
    Hs = SparseArrayFused(UInt32(64))

    # PSM accumulator (one slot per possible column = up to n_precursors)
    unscored_psms = [MainUnscoredPSM{Float32}() for _ in 1:max(8, n_precursors)]

    # id_to_col: map prec_id -> column. FusedQuadEst uses (prec_id-1)*3 + iso_pass,
    # so size is 3*n_precursors + 1; for Standard it's just n_precursors.
    id_to_col_size = max(3 * n_precursors + 1, 16)
    id_to_col = DensePrecMap{UInt16}(id_to_col_size)

    # Scratch: small initial cap so growth path is exercisable.
    scratch = FusedScratch(8)

    # Convert peak m/z to corrected_mz / obs_low / obs_high. Per-peak tol
    # supports stressing intensity-dependent mass error (peaks of varying
    # quality get different acceptance widths). The mem instance keeps a
    # conservative ppm so conservative_half_width(mem, frag_mz) doesn't
    # gate matches that the per-peak window admits.
    n_peaks = length(peak_mz)
    tol_vec = ppm_tol isa Vector{Float32} ?
        ppm_tol :
        fill(Float32(ppm_tol), n_peaks)
    @assert length(tol_vec) == n_peaks "ppm_tol vector must match length(peak_mz)"
    conservative_ppm = max(maximum(tol_vec; init=0.0f0), 20.0f0)
    mem = SimpleMassErrorModel(0.0f0, (conservative_ppm, conservative_ppm))
    scan_corrected_mz = copy(peak_mz)            # already "corrected" (no offset)
    scan_obs_low  = Float32[peak_mz[i] - peak_mz[i] * tol_vec[i] * 1f-6 for i in 1:n_peaks]
    scan_obs_high = Float32[peak_mz[i] + peak_mz[i] * tol_vec[i] * 1f-6 for i in 1:n_peaks]
    scan_int = Vector{Union{Missing, Float32}}(undef, n_peaks)
    for i in 1:n_peaks
        scan_int[i] = peak_int[i] isa Missing ? missing : Float32(peak_int[i])
    end

    # Working buffers used inside the iso loop. `getFragAbundance!` (called
    # transitively from `getFragIsotopes!`) iterates `f ∈ 0:length(prec_isotopes)-1`
    # and writes `frag_isotopes[f+1]` — so `isotopes_buf` must be ≥ the
    # length of `prec_trans_buf` or `--check-bounds=yes` runs surface a
    # `BoundsError` that's silenced by `@inbounds` otherwise.
    prec_trans_buf = zeros(Float32, 5)   # 5 precursor isotopes is the standard size
    isotopes_buf   = zeros(Float32, max(n_frag_isotopes, length(prec_trans_buf)))

    # Library fragments: pack into CompactFrag. Tuple may be:
    #   (pid, mz, intensity, rank)             — defaults to y-ion, ion_pos=1
    #   (pid, mz, intensity, rank, ion_type, ion_pos)
    #     where ion_type ∈ (:y, :b, :p), ion_pos is the ladder position.
    compact_frags = CompactFrag{Float32}[]
    for f in frags
        pid, mz, int, rank = f[1], f[2], f[3], f[4]
        ion_type = length(f) >= 5 ? f[5] : :y
        ion_pos  = length(f) >= 6 ? f[6] : UInt8(1)
        is_y = ion_type === :y
        is_b = ion_type === :b
        is_p = ion_type === :p
        push!(compact_frags, CompactFrag(
            UInt32(pid), Float32(mz), Float16(int),
            is_y, is_b, is_p, false,           # flags
            UInt8(1),                          # frag_charge
            UInt8(ion_pos),                    # ion_position
            UInt8(prec_charges[pid]),          # prec_charge
            UInt8(rank),                       # rank
            UInt8(0),                          # sulfur_count
        ))
    end
    ion_list = StandardFragmentLookup{Float32}(compact_frags, prec_frag_ranges)

    # NCE model — only used when SplineFragmentLookup is in play; CompactFrag
    # path returns ConstantType regardless. A no-op model is fine.
    nce_model = PiecewiseNceModel(0.0f0, 0.0f0, 30.0f0, 30.0f0, 0.0f0)

    # Quad transmission: flat passthrough between [quad_min, quad_max].
    qfunc = SquareQuadFunction(quad_min, quad_max, (quad_min + quad_max) / 2)

    # precursors_passed: identity ordering by default
    precursors_passed = UInt32.(1:n_precursors)

    return (
        kind = kind,
        Hs = Hs,
        unscored_psms = unscored_psms,
        id_to_col = id_to_col,
        scratch = scratch,
        scan_corrected_mz = scan_corrected_mz,
        scan_obs_low = scan_obs_low,
        scan_obs_high = scan_obs_high,
        peak_mz_len = n_peaks,
        isotopes_buf = isotopes_buf,
        prec_trans_buf = prec_trans_buf,
        ion_list = ion_list,
        nce_model = nce_model,
        precursors_passed = precursors_passed,
        prec_range = prec_range,
        prec_mzs = prec_mzs,
        prec_charges = prec_charges,
        prec_sulfur_counts = prec_sulfur_counts,
        prec_irts = prec_irts,
        iso_splines = ISO_SPLINES,
        qfunc = qfunc,
        mem = mem,
        scan_int = scan_int,
        scan_irt = scan_irt,
        irt_tol = irt_tol,
        frag_mz_bounds = (100.0f0, 2000.0f0),
        n_frag_isotopes = n_frag_isotopes,
        isotope_err_bounds = isotope_err_bounds,
    )
end

# Thin wrapper: invoke run_fused! with all the fixture's args.
function call_run_fused!(fx; kwargs...)
    f = merge(fx, values(kwargs))
    return run_fused!(
        f.kind, f.Hs, f.unscored_psms, f.id_to_col, f.scratch,
        f.scan_corrected_mz, f.scan_obs_low, f.scan_obs_high, f.peak_mz_len,
        f.isotopes_buf, f.prec_trans_buf, f.ion_list, f.nce_model,
        f.precursors_passed, f.prec_range,
        f.prec_mzs, f.prec_charges, f.prec_sulfur_counts, f.prec_irts,
        f.iso_splines, f.qfunc, f.mem, f.scan_int,
        f.scan_irt, f.irt_tol, f.frag_mz_bounds, f.n_frag_isotopes,
        f.isotope_err_bounds)
end

# =============================================================================
@testset "run_fused! orchestration" begin

    # ------------------------------------------------------------
    @testset "empty prec_range returns (0,0), Hs untouched" begin
        fx = make_fused_fixture(prec_range = 1:0)
        n_match, n_miss = call_run_fused!(fx)
        @test n_match == 0
        @test n_miss == 0
        @test fx.Hs.n == 0
        @test fx.Hs.n_vals == 0
        @test fx.scratch.n == 0
        @test fx.scratch.miss_n == 0
    end

    # ------------------------------------------------------------
    @testset "single precursor, all fragments match" begin
        fx = make_fused_fixture()  # 2 fragments at 200 and 300; peaks placed there
        n_match, n_miss = call_run_fused!(fx)
        @test n_match == 2          # both fragments match their peak
        @test n_miss  == 0
        @test fx.Hs.n == 1          # one column allocated
        @test fx.Hs.n_vals == 2     # two entries in that column
        # id_to_col must record the column for prec_idx=1 (Standard maps to prec_idx).
        @test fx.id_to_col[1] == UInt16(1)
        # MainUnscoredPSM should have been updated for col=1.
        @test fx.unscored_psms[1].precursor_idx == UInt32(1)
        # Both peaks landed in the matched buffer (rows 1 and 2).
        rows_in_col = sort([fx.Hs.rowval[i] for i in 1:fx.Hs.n_vals])
        @test rows_in_col == UInt32[1, 2]
    end

    # ------------------------------------------------------------
    @testset "grows tuning PSM accumulator before writing new columns" begin
        n_precs = 10
        frags = [(UInt32(pid), Float32(200 + pid), 1000.0f0, UInt8(0))
                 for pid in 1:n_precs]
        fx = make_fused_fixture(
            n_precursors = n_precs,
            prec_mzs = fill(500.0f0, n_precs),
            prec_charges = fill(UInt8(2), n_precs),
            prec_sulfur_counts = fill(UInt8(0), n_precs),
            prec_irts = fill(50.0f0, n_precs),
            frags = frags,
            prec_frag_ranges = UInt64.(1:(n_precs + 1)),
            peak_mz = Float32[200 + pid for pid in 1:n_precs],
            peak_int = fill(1000.0f0, n_precs),
            prec_range = 1:n_precs,
        )
        tuning_psms = [TuningUnscoredPSM{Float32}() for _ in 1:4]

        n_match, n_miss = call_run_fused!(fx, unscored_psms = tuning_psms)

        @test n_match == n_precs
        @test n_miss == 0
        @test fx.Hs.n == n_precs
        @test length(tuning_psms) >= n_precs
        @test [tuning_psms[i].precursor_idx for i in 1:n_precs] == UInt32.(1:n_precs)
    end

    # ------------------------------------------------------------
    @testset "single precursor, no peak in window → all misses, no column" begin
        # Peaks far away from fragment m/z (200, 300) — no match possible.
        fx = make_fused_fixture(peak_mz = Float32[1000.0, 1100.0],
                                peak_int = Float32[500.0, 500.0])
        n_match, n_miss = call_run_fused!(fx)
        @test n_match == 0
        @test n_miss  == 2          # both fragments recorded as misses
        @test fx.Hs.n == 0          # no column allocated (col_started never true)
        @test fx.Hs.n_vals == 0
        # id_to_col stays untouched for prec_idx=1.
        @test fx.id_to_col[1] == UInt16(0)
    end

    # ------------------------------------------------------------
    @testset "do_prec_check: precursor outside iRT window skipped" begin
        # scan_irt=50, irt_tol=5 — set prec_irt=100 so |100-50|=50 > 5.
        fx = make_fused_fixture(prec_irts = Float32[100.0])
        n_match, n_miss = call_run_fused!(fx)
        @test n_match == 0
        @test n_miss  == 0          # filter triggers `continue`, no fragments seen
        @test fx.Hs.n == 0
        @test fx.id_to_col[1] == UInt16(0)
    end

    # ------------------------------------------------------------
    @testset "do_prec_check: precursor outside m/z window skipped" begin
        # quad_min=480, quad_max=520, isotope_err_bounds=(0,0) → window=[480,520].
        # Set prec_mz=600 to fall outside.
        fx = make_fused_fixture(prec_mzs = Float32[600.0])
        n_match, n_miss = call_run_fused!(fx)
        @test n_match == 0
        @test n_miss  == 0
        @test fx.Hs.n == 0
        @test fx.id_to_col[1] == UInt16(0)
    end

    # ------------------------------------------------------------
    @testset "FusedRTIndexed: iRT filter bypassed" begin
        # iRT 100 with scan_irt=50, irt_tol=5 → FusedStandard skips this
        # precursor; FusedRTIndexed processes it (check_prec_filters=false).
        # Keep prec_mz inside the quad window so transmission is non-zero.
        rt_kind = FusedRTIndexed(FullPrecCapture(), UInt8(7))
        fx_skip = make_fused_fixture(prec_irts = Float32[100.0])
        n1, _ = call_run_fused!(fx_skip)
        @test n1 == 0   # FusedStandard filters it out

        fx_keep = make_fused_fixture(prec_irts = Float32[100.0], kind = rt_kind)
        n2, _ = call_run_fused!(fx_keep)
        @test n2 == 2   # FusedRTIndexed bypasses iRT check, both fragments match
    end

    # ------------------------------------------------------------
    @testset "rank filter: high-rank fragments skipped" begin
        # Make the second fragment rank=15 (above FusedStandard's threshold=7).
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0)),
                 (UInt32(1), 300.0f0,  800.0f0, UInt8(15))]
        fx = make_fused_fixture(frags = frags)
        n_match, n_miss = call_run_fused!(fx)
        @test n_match == 1          # only the rank-0 fragment matches
        @test n_miss  == 0
        @test fx.Hs.n == 1
    end

    # ------------------------------------------------------------
    @testset "frag_mz_bounds: out-of-window iso m/z dropped" begin
        # Restrict the scan window so fragment at 300 falls outside.
        fx = make_fused_fixture()  # default n_frag_isotopes=1 means iso_idx=0 only
        # Override frag_mz_bounds to exclude m/z=300.
        n_match, n_miss = call_run_fused!(fx, frag_mz_bounds = (100.0f0, 250.0f0))
        @test n_match == 1          # only m/z=200 fragment matches
        @test n_miss  == 0          # 300 was filtered before match attempt
    end

    # ------------------------------------------------------------
    @testset "mixed match + miss: scratch records both" begin
        # Two fragments at 200, 300; only one peak at 200.
        fx = make_fused_fixture(peak_mz = Float32[200.0],
                                peak_int = Float32[1500.0])
        n_match, n_miss = call_run_fused!(fx)
        @test n_match == 1
        @test n_miss  == 1
        @test fx.Hs.n == 1          # column allocated because at least one match
        @test fx.Hs.n_vals == 2     # 1 match + 1 miss flushed into the column
    end

    # ------------------------------------------------------------
    @testset "multi-iso: n_frag_isotopes=2 → 2 entries per fragment" begin
        # Charge-1 fragment at m/z 200 → iso 0 at 200, iso 1 at ~201.0033.
        # Place peaks at both m/z values.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0))]
        peak_mz  = Float32[200.0, 201.0033]
        peak_int = Float32[1500.0, 800.0]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 2],
            peak_mz = peak_mz,
            peak_int = peak_int,
            n_frag_isotopes = 2,
        )
        n_match, n_miss = call_run_fused!(fx)
        @test n_match == 2          # iso 0 and iso 1 both match
        @test n_miss  == 0
        @test fx.Hs.n == 1          # one precursor → one column
        @test fx.Hs.n_vals == 2     # two entries in that column
        # Both isotope flags present (0 and 1).
        isos = sort([Pioneer.isotope_at(fx.Hs, i) for i in 1:fx.Hs.n_vals])
        @test isos == UInt8[0, 1]
    end

    # ------------------------------------------------------------
    @testset "iso_anchor advances within fragment (iso 1 search starts ≥ iso 0 anchor)" begin
        # Two fragments with overlapping iso m/z ranges. With proper iso_anchor
        # advancement, iso 1's bsearch lo ≥ iso 0's start_idx — testing here
        # means the higher iso of fragment A doesn't accidentally swallow
        # fragment B's monoisotopic peak.
        # Fragment A at 200 charge 1 → iso 0 at 200, iso 1 at 201.0033.
        # Fragment B at 201 charge 1 → iso 0 at 201, iso 1 at 202.0033.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0)),
                 (UInt32(1), 201.0f0,  900.0f0, UInt8(0))]
        peak_mz  = Float32[200.0, 201.0033, 201.0, 202.0033]
        peak_int = Float32[1000, 800, 900, 700]
        sort_order = sortperm(peak_mz)
        peak_mz  = peak_mz[sort_order]
        peak_int = peak_int[sort_order]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 3],
            peak_mz = peak_mz,
            peak_int = peak_int,
            n_frag_isotopes = 2,
        )
        n_match, _ = call_run_fused!(fx)
        # 2 fragments × 2 isos = up to 4 matches. With 20-ppm tolerance the
        # peaks at 201.0 and 201.0033 are very close (~16 ppm apart at m/z 201);
        # depending on bsearch exact behavior, the sort+nearest scan may match
        # both or fewer. Assert at least 3 (fragment A iso 0 + fragment B both
        # isos) to confirm iso_anchor doesn't *over*-skip.
        @test n_match >= 3
    end

    # ------------------------------------------------------------
    @testset "FusedQuadEst: 3 iso passes → 3 distinct columns per precursor" begin
        # FusedQuadEst always runs 3 iso passes per precursor and writes a
        # separate column per pass (match_column_id = (prec_idx-1)*3 + iso_pass).
        # Each pass uses one-hot precursor transmission for iso (pass-1).
        # With prec_mz=500 charge=2, prec isotopes 0,1,2 sit at 500, 500.5, 501.
        # Default quad window [480, 520] passes all three.
        kind = FusedQuadEst()
        # Single fragment at m/z 200 charge 1 → fragment isotopes at
        # 200, 201.0033, 202.0066, 203.01 (FusedQuadEst hardcodes max_iso=3).
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0))]
        peak_mz  = Float32[200.0, 201.0033, 202.0066, 203.01]
        peak_int = fill(1000.0f0, 4)
        fx = make_fused_fixture(
            kind = kind,
            frags = frags,
            prec_frag_ranges = UInt64[1, 2],
            peak_mz = peak_mz,
            peak_int = peak_int,
            n_frag_isotopes = 4,    # ignored by FusedQuadEst (uses 0:3)
        )
        n_match, _ = call_run_fused!(fx)
        @test n_match >= 4          # at least the monoisotopic match per pass

        # 3 columns allocated, one per iso pass.
        @test fx.Hs.n == 3
        # match_column_id for FusedQuadEst: key = (prec-1)*3 + pass
        # prec=1 → keys 1, 2, 3 each map to a column
        @test fx.id_to_col[1] == UInt16(1)
        @test fx.id_to_col[2] == UInt16(2)
        @test fx.id_to_col[3] == UInt16(3)
    end

    # ------------------------------------------------------------
    @testset "FusedQuadEst: record_match! is a no-op (unscored_psms untouched)" begin
        kind = FusedQuadEst()
        fx = make_fused_fixture(kind = kind, n_frag_isotopes = 4)
        # Snapshot psm fields before.
        before_pidx = fx.unscored_psms[1].precursor_idx
        before_brank = fx.unscored_psms[1].best_rank
        n_match, _ = call_run_fused!(fx)
        @test n_match >= 1
        @test fx.unscored_psms[1].precursor_idx == before_pidx
        @test fx.unscored_psms[1].best_rank     == before_brank
    end

    # ------------------------------------------------------------
    @testset "Hs.colptr correctly delimits each column" begin
        # Two precursors, each with one matching fragment. Resulting Hs has
        # 2 columns, colptr should be [1, 2, 3] (1 entry each, 1-indexed CSC).
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0)),
                 (UInt32(2), 300.0f0,  800.0f0, UInt8(0))]
        fx = make_fused_fixture(
            n_precursors = 2,
            prec_mzs = Float32[500.0, 510.0],
            prec_charges = UInt8[2, 2],
            prec_sulfur_counts = UInt8[0, 0],
            prec_irts = Float32[50.0, 50.0],
            frags = frags,
            prec_frag_ranges = UInt64[1, 2, 3],   # prec1 → frag 1; prec2 → frag 2
            peak_mz  = Float32[200.0, 300.0],
            peak_int = Float32[1500.0, 1200.0],
            prec_range = 1:2,
        )
        n_match, _ = call_run_fused!(fx)
        @test n_match == 2
        @test fx.Hs.n == 2
        @test fx.Hs.n_vals == 2
        # colptr[c]:colptr[c+1]-1 gives the entry range for column c.
        # With one entry per column, expect colptr[1]=1, colptr[2]=2, colptr[3]=3.
        @test fx.Hs.colptr[1] == UInt32(1)
        @test fx.Hs.colptr[2] == UInt32(2)
        @test fx.Hs.colptr[3] == UInt32(3)
    end

    # ------------------------------------------------------------
    @testset "intensity-dependent tol: tight peak rejects, wide peak admits" begin
        # Single fragment at m/z 200; single peak shifted 100 ppm to 200.02.
        # Tight peak (20 ppm tol → ±0.004 obs window centered at 200.02) does
        # NOT contain iso_mz=200 → no match. Wide peak (500 ppm → ±0.1) does
        # contain it → matches. Demonstrates the per-peak obs window gates
        # matches independently of the conservative_half_width fragment window.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0))]

        fx_tight = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 2],
            peak_mz = Float32[200.02],
            peak_int = Float32[1500.0],
            ppm_tol = 20.0f0,
        )
        n_t, _ = call_run_fused!(fx_tight)
        @test n_t == 0                       # tight obs window excludes iso_mz=200

        fx_wide = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 2],
            peak_mz = Float32[200.02],
            peak_int = Float32[1500.0],
            ppm_tol = 500.0f0,
        )
        n_w, _ = call_run_fused!(fx_wide)
        @test n_w == 1                       # wide obs window admits it
    end

    # ------------------------------------------------------------
    @testset "intensity-dependent tol: heterogeneous per-peak widths" begin
        # Same fragment, two peaks: one tight (rejects), one wide (admits).
        # The wide peak is intentionally placed *farther* from iso_mz so
        # that nearest-by-abs would prefer the tight one — but the tight
        # one's obs window excludes iso_mz, forcing the algorithm to fall
        # through to the wide one.
        # frag at 200, iso_mz=200. Peak 1 at 200.005 tight (±0.004 → [200.001, 200.009])
        # → excludes iso_mz=200. Peak 2 at 200.05 wide (±0.1 → [199.95, 200.15])
        # → includes iso_mz. Wide peak should match.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0))]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 2],
            peak_mz  = Float32[200.005, 200.05],
            peak_int = Float32[1500.0, 1500.0],
            ppm_tol  = Float32[20.0, 500.0],   # tight, wide
        )
        n_match, _ = call_run_fused!(fx)
        @test n_match == 1
        @test fx.Hs.rowval[1] == UInt32(2)    # the wide peak (index 2) matched
    end

    # ------------------------------------------------------------
    @testset "intensity-dependent tol: nearest-by-abs wins when both admit" begin
        # Two peaks both within wide tolerance, fragment iso_mz between them.
        # Algorithm should pick the closer one (nearest-by-abs).
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0))]
        peak_mz  = Float32[199.96, 200.02]   # iso 200 closer to peak 2 (Δ=0.02 vs 0.04)
        peak_int = Float32[1500.0, 1500.0]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 2],
            peak_mz = peak_mz,
            peak_int = peak_int,
            ppm_tol = 500.0f0,                # both peaks well within tol
        )
        n_match, _ = call_run_fused!(fx)
        @test n_match == 1
        @test fx.Hs.rowval[1] == UInt32(2)   # peak 2 (closer) wins
    end

    # ------------------------------------------------------------
    @testset "M+1 of frag A overlaps with M+0 of frag B (no cross-talk)" begin
        # Charge-1 fragments: A at 200.0, B at 200.5.
        # A's M+1 = 201.0033 (well above B's M+0 = 200.5), B's M+1 = 201.5066.
        # Peak layout (sorted ascending): 200.0, 200.5, 201.0033, 201.5066.
        # All four isos should match their own peak — no cross-talk.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0)),
                 (UInt32(1), 200.5f0,  900.0f0, UInt8(0))]
        peak_mz  = Float32[200.0, 200.5, 201.0033, 201.5066]
        peak_int = fill(1500.0f0, 4)
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 3],
            peak_mz = peak_mz,
            peak_int = peak_int,
            n_frag_isotopes = 2,
        )
        n_match, _ = call_run_fused!(fx)
        @test n_match == 4                   # all 4 isos match
        @test fx.Hs.n_vals == 4
        # Matched rows must include each peak index exactly once.
        rows = sort([fx.Hs.rowval[i] for i in 1:fx.Hs.n_vals])
        @test rows == UInt32[1, 2, 3, 4]
    end

    # ------------------------------------------------------------
    @testset "two fragments hit same peak: finalize_column! dedupes with sum" begin
        # Two fragments at the same m/z (200) both match the single peak at 200.
        # `finalize_column!` deduplicates entries that share (column, row) by
        # summing their predicted intensities — the count stays at 1 entry per
        # row, not 2. (The two scratch entries are collapsed into one Hs row.)
        # Stresses dedup logic in finalize_column! and confirms iso_anchor /
        # lower advancement does NOT block fragment B from re-finding the peak.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(1)),
                 (UInt32(1), 200.0f0,  900.0f0, UInt8(2))]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 3],
            peak_mz = Float32[200.0],
            peak_int = Float32[1500.0],
        )
        n_match, _ = call_run_fused!(fx)
        @test n_match == 2                   # both fragments matched (counter)
        @test fx.Hs.n == 1                   # one column
        @test fx.Hs.n_vals == 1              # but only ONE entry after dedup
        @test fx.Hs.rowval[1] == UInt32(1)
        # Predicted intensities summed: pred_int_A + pred_int_B.
        @test fx.Hs.nzval[1] > Float32(1000) # > either fragment alone
        @test Pioneer.rank_at(fx.Hs, 1) == UInt8(2)
    end

    # ------------------------------------------------------------
    @testset "all fragment isos outside frag_mz_bounds: 0 match, 0 miss" begin
        # Single fragment at m/z 500; restrict frag_mz_bounds to [100, 250].
        # iso 0 = 500 is outside → continue (no match, no miss recorded).
        frags = [(UInt32(1), 500.0f0, 1000.0f0, UInt8(0))]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 2],
        )
        n_match, n_miss = call_run_fused!(fx, frag_mz_bounds = (100.0f0, 250.0f0))
        @test n_match == 0
        @test n_miss  == 0                   # filtered before match attempt
        @test fx.Hs.n == 0
    end

    # ------------------------------------------------------------
    @testset "Missing peak intensity: match recorded with int_obs=0" begin
        # Peak m/z is fine, intensity is missing → ismissing-branch fires,
        # int_obs falls back to 0. Match is still recorded.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0))]
        peak_int = Vector{Any}([missing, Float32(1200.0)])
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 2],
            peak_mz = Float32[200.0, 300.0],
            peak_int = peak_int,
        )
        n_match, n_miss = call_run_fused!(fx)
        @test n_match == 1
        @test fx.Hs.n_vals == 1
        # Found row 1 (m/z 200), and x (observed intensity) is 0 since missing.
        @test fx.Hs.rowval[1] == UInt32(1)
        @test fx.Hs.x[1] === 0.0f0
    end

    # ------------------------------------------------------------
    @testset "empty fragment list for a precursor: passes filters but 0 matches" begin
        # Two precursors; precursor 1 has no fragments (empty range), 2 does.
        # prec_frag_ranges = [1, 1, 2] → prec1 covers [1,1) (empty), prec2 covers [1,2).
        frags = [(UInt32(2), 200.0f0, 1000.0f0, UInt8(0))]
        fx = make_fused_fixture(
            n_precursors = 2,
            prec_mzs = Float32[500.0, 510.0],
            prec_charges = UInt8[2, 2],
            prec_sulfur_counts = UInt8[0, 0],
            prec_irts = Float32[50.0, 50.0],
            frags = frags,
            prec_frag_ranges = UInt64[1, 1, 2],
            peak_mz = Float32[200.0],
            peak_int = Float32[1500.0],
            prec_range = 1:2,
        )
        n_match, n_miss = call_run_fused!(fx)
        @test n_match == 1                   # only precursor 2 contributes
        @test n_miss  == 0
        @test fx.Hs.n == 1                   # one column allocated (for precursor 2)
        @test fx.id_to_col[1] == UInt16(0)   # precursor 1 never started a column
        @test fx.id_to_col[2] == UInt16(1)
    end

    # ------------------------------------------------------------
    @testset "fragment all-miss but precursor matches another fragment" begin
        # Two fragments: A matches a peak, B has no peak in window.
        # Column allocated by A; B's miss flushed alongside A's match.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0)),
                 (UInt32(1), 999.0f0,  500.0f0, UInt8(0))]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 3],
            peak_mz  = Float32[200.0, 300.0],   # 300 is far from 999
            peak_int = Float32[1500.0, 800.0],
        )
        n_match, n_miss = call_run_fused!(fx)
        @test n_match == 1
        @test n_miss  == 1
        @test fx.Hs.n == 1                   # column allocated
        @test fx.Hs.n_vals == 2              # 1 match + 1 miss in this column
        # The miss has matched=false in the packed meta.
        flags = sort([Pioneer.matched_at(fx.Hs, i) for i in 1:fx.Hs.n_vals])
        @test flags == Bool[false, true]
    end

    # ------------------------------------------------------------
    @testset "unscored_psms field values: y-ion match populates y_count, longest_y" begin
        # Single y-ion at ion_pos=5, m/z 200, rank 0. Match against a peak.
        # Expect: y_count=1, y_int=peak intensity, longest_y=5, b_count=0,
        # best_rank=0, isotope_count=0.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0), :y, UInt8(5))]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 2],
            peak_mz = Float32[200.0],
            peak_int = Float32[1500.0],
        )
        n_match, _ = call_run_fused!(fx)
        @test n_match == 1
        psm = fx.unscored_psms[1]
        @test psm.y_count   == UInt8(1)
        @test psm.y_int     === Float32(1500.0)
        @test psm.longest_y == UInt8(5)
        @test psm.b_count   == UInt8(0)
        @test psm.best_rank == UInt8(0)
        @test psm.isotope_count == UInt8(0)     # iso 0 only
        @test psm.precursor_idx == UInt32(1)
    end

    # ------------------------------------------------------------
    @testset "unscored_psms: b-ion match populates b_count, longest_b" begin
        # Two b-ions at different ion_pos. longest_b should be the larger.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0), :b, UInt8(3)),
                 (UInt32(1), 300.0f0, 1000.0f0, UInt8(0), :b, UInt8(7))]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 3],
        )
        n_match, _ = call_run_fused!(fx)
        @test n_match == 2
        psm = fx.unscored_psms[1]
        @test psm.b_count   == UInt8(2)
        @test psm.longest_b == UInt8(7)         # max of (3, 7)
        @test psm.b_int     === Float32(1500.0 + 1200.0)
        @test psm.y_count   == UInt8(0)
        @test psm.longest_y == UInt8(0)
    end

    # ------------------------------------------------------------
    @testset "unscored_psms: M+1 isotope updates isotope_count + y_count_iso" begin
        # y-ion charge 1, n_frag_isotopes=2. Both M+0 and M+1 match.
        # M+0 → y_count=1, longest_y. M+1 → isotope_count=1, y_count_iso=1,
        # longest_y_iso = ion_pos.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0), :y, UInt8(4))]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 2],
            peak_mz  = Float32[200.0, 201.0033],
            peak_int = Float32[1500.0, 800.0],
            n_frag_isotopes = 2,
        )
        n_match, _ = call_run_fused!(fx)
        @test n_match == 2
        psm = fx.unscored_psms[1]
        @test psm.y_count       == UInt8(1)     # only mono counted in y_count
        @test psm.y_count_iso   == UInt8(1)     # M+1 counted in y_count_iso
        @test psm.isotope_count == UInt8(1)
        @test psm.longest_y     == UInt8(4)
        @test psm.longest_y_iso == UInt8(4)
    end

    # ------------------------------------------------------------
    @testset "unscored_psms: error accumulates |ppm_err|" begin
        # Place peak slightly off-center: peak at 200.005 vs frag at 200.0
        # → ppm_err = (200.0 - 200.005) / (200.0 * 1e-6) ≈ -25 ppm.
        # error field accumulates abs(ppm_err) ≈ 25.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0))]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 2],
            peak_mz = Float32[200.005],
            peak_int = Float32[1500.0],
            ppm_tol = 50.0f0,                    # ensure peak is admitted
        )
        n_match, _ = call_run_fused!(fx)
        @test n_match == 1
        @test fx.unscored_psms[1].error ≈ 25.0f0  atol=1.0
    end

    # ------------------------------------------------------------
    @testset ">2 precursors interleaved with mixed match/no-match" begin
        # 4 precursors: 1 matches, 2 has no peak in window (all-miss → no col),
        # 3 has empty fragments, 4 matches. Expect 2 columns total.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0)),
                 (UInt32(2), 999.0f0,  500.0f0, UInt8(0)),  # peak doesn't reach
                 (UInt32(4), 300.0f0, 1000.0f0, UInt8(0))]
        # Per-precursor ranges: 1→[1,2), 2→[2,3), 3→[3,3) (empty), 4→[3,4).
        fx = make_fused_fixture(
            n_precursors = 4,
            prec_mzs = Float32[500.0, 510.0, 505.0, 515.0],
            prec_charges = UInt8[2, 2, 2, 2],
            prec_sulfur_counts = UInt8[0, 0, 0, 0],
            prec_irts = Float32[50.0, 50.0, 50.0, 50.0],
            frags = frags,
            prec_frag_ranges = UInt64[1, 2, 3, 3, 4],
            peak_mz  = Float32[200.0, 300.0],
            peak_int = Float32[1500.0, 1200.0],
            prec_range = 1:4,
        )
        n_match, n_miss = call_run_fused!(fx)
        @test n_match == 2                       # precursor 1 + 4
        @test n_miss  == 1                       # precursor 2's miss (col not allocated → discarded)
        @test fx.Hs.n == 2                       # only 2 columns
        @test fx.id_to_col[1] == UInt16(1)       # precursor 1 got column 1
        @test fx.id_to_col[2] == UInt16(0)       # precursor 2 had no match
        @test fx.id_to_col[3] == UInt16(0)       # precursor 3 had no fragments
        @test fx.id_to_col[4] == UInt16(2)       # precursor 4 got column 2
    end

    # ------------------------------------------------------------
    @testset "FusedStandard + PartialPrecCapture: same outcome as FullPrecCapture" begin
        # Default fixture uses FullPrecCapture; swap to PartialPrecCapture and
        # confirm match counts are identical for a simple case where partial-
        # vs-full precursor capture produces the same fragment intensities.
        # (Full charge-state coverage with sulfur=0 keeps both branches equivalent.)
        kind_full    = FusedStandard(FullPrecCapture(),    UInt8(7))
        kind_partial = FusedStandard(PartialPrecCapture(), UInt8(7))
        fx_full = make_fused_fixture(kind = kind_full)
        n_full, _ = call_run_fused!(fx_full)

        fx_partial = make_fused_fixture(kind = kind_partial)
        n_partial, _ = call_run_fused!(fx_partial)

        @test n_full == n_partial
        @test fx_full.Hs.n == fx_partial.Hs.n
        @test fx_full.Hs.n_vals == fx_partial.Hs.n_vals
    end

    # ------------------------------------------------------------
    @testset "unscored_psms: p-ion match increments p_count, leaves y/b at 0" begin
        # P-ion (precursor ion) match. apply_main_scoring! routes via the
        # isP branch on the mono path, incrementing p_count. b_count, y_count,
        # longest_b, longest_y stay at 0.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0), :p, UInt8(2))]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 2],
            peak_mz = Float32[200.0],
            peak_int = Float32[1500.0],
        )
        n_match, _ = call_run_fused!(fx)
        @test n_match == 1
        psm = fx.unscored_psms[1]
        @test psm.p_count             == UInt8(1)
        @test psm.b_count             == UInt8(0)
        @test psm.y_count             == UInt8(0)
        @test psm.longest_b           == UInt8(0)
        @test psm.longest_y           == UInt8(0)
        @test psm.non_cannonical_count == UInt8(0)
        @test psm.b_int               === Float32(0)
        @test psm.y_int               === Float32(0)
    end

    # ------------------------------------------------------------
    @testset "unscored_psms: non-canonical ion (no flag) increments fallback counter" begin
        # apply_main_scoring! has an `else` branch that catches frags
        # which are neither y, b, nor p — increments non_cannonical_count.
        # Build a frag with all flags false (we hijack the constructor by
        # passing :other ion_type → none of is_y/is_b/is_p set).
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(0), :other, UInt8(1))]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 2],
            peak_mz = Float32[200.0],
            peak_int = Float32[1500.0],
        )
        n_match, _ = call_run_fused!(fx)
        @test n_match == 1
        psm = fx.unscored_psms[1]
        @test psm.non_cannonical_count == UInt8(1)
        @test psm.b_count              == UInt8(0)
        @test psm.y_count              == UInt8(0)
        @test psm.p_count              == UInt8(0)
    end

    # ------------------------------------------------------------
    @testset "unscored_psms iso branch: isotope_count + y_count_iso updated" begin
        # Two fragments with different ranks; both M+1 isotopes match.
        # iso branch increments isotope_count and y_count_iso / longest_y_iso.
        # Frag A rank=2, Frag B rank=5. Both M+1 should match.
        # Expect: isotope_count=2, y_count_iso=2, longest_y_iso=6.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(2), :y, UInt8(3)),
                 (UInt32(1), 300.0f0,  900.0f0, UInt8(5), :y, UInt8(6))]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 3],
            peak_mz  = Float32[200.0, 201.0033, 300.0, 301.0033],
            peak_int = fill(1500.0f0, 4),
            n_frag_isotopes = 2,
        )
        n_match, _ = call_run_fused!(fx)
        @test n_match == 4                           # 2 fragments × 2 isos
        psm = fx.unscored_psms[1]
        @test psm.isotope_count == UInt8(2)          # 2 M+1 matches
        @test psm.y_count_iso   == UInt8(2)          # both Y-ions on iso pass
        @test psm.longest_y_iso == UInt8(6)          # max ion_pos on iso Y
    end

    # ------------------------------------------------------------
    @testset "unscored_psms: best_rank tracks minimum across mono matches" begin
        # Three Y-ions: ranks 4, 1, 7. After all match, best_rank should be 1.
        frags = [(UInt32(1), 200.0f0, 1000.0f0, UInt8(4), :y, UInt8(1)),
                 (UInt32(1), 300.0f0,  900.0f0, UInt8(1), :y, UInt8(2)),
                 (UInt32(1), 400.0f0,  800.0f0, UInt8(7), :y, UInt8(3))]
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, 4],
            peak_mz  = Float32[200.0, 300.0, 400.0],
            peak_int = fill(1500.0f0, 3),
        )
        n_match, _ = call_run_fused!(fx)
        @test n_match == 3
        psm = fx.unscored_psms[1]
        @test psm.best_rank == UInt8(1)
        @test psm.y_count   == UInt8(3)
    end

    # ------------------------------------------------------------
    @testset "scratch growth path (precursor exceeds initial capacity)" begin
        # Initial scratch capacity is 8. Build 12 fragments at distinct m/z
        # with peaks at each — forces grow_fused_scratch! mid-loop.
        n = 12
        frags = [(UInt32(1), Float32(200 + 10*i), 1000.0f0, UInt8(0)) for i in 0:(n-1)]
        peak_mz  = Float32[200 + 10*i for i in 0:(n-1)]
        peak_int = fill(1500.0f0, n)
        fx = make_fused_fixture(
            frags = frags,
            prec_frag_ranges = UInt64[1, n + 1],
            peak_mz = peak_mz,
            peak_int = peak_int,
        )
        @test length(fx.scratch.row) == 8   # confirm initial cap
        n_match, _ = call_run_fused!(fx)
        @test n_match == n
        @test length(fx.scratch.row) >= n   # grew to fit
        @test fx.Hs.n_vals == n
    end

end
