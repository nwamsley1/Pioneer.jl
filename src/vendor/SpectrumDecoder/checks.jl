# ©2026 Washington University
#
# This file is part of Pioneer.jl and is licensed under the
# Washington University non-commercial research license.
# See the LICENSE file in the project root for full terms.

# Internal-consistency suite. Each check runs over every scan of a file and reports pass/total.

struct CheckResult
    group::String
    name::String
    pass::Int
    total::Int
    note::String
end
ok(r::CheckResult) = r.pass == r.total

const ISOTOPE_SPACING = 1.00335  # 13C - 12C, Da

"""
    check_consistency(f::Acquisition) -> Vector{CheckResult}
"""
median_width(wl) = (w = sort(last.(wl)); w[(length(w) + 1) ÷ 2])

function check_consistency(f::Acquisition)
    b = f.bytes; S = scans(f); n = length(S)
    R = CheckResult[]
    add!(g, name, p, t, note = "") = push!(R, CheckResult(g, name, p, t, note))

    # ---------------- structure closes ----------------
    _, dirs = device_directory(b)
    add!("structure", "device directory offsets inside file", count(d -> 0 < d[3] < length(b), dirs), length(dirs))
    add!("structure", "scan count == last - first + 1", n == f.last_scan - f.first_scan + 1, 1)
    add!("structure", "scan table index / scan number sequential",
         count(j -> ru32(b, f.scan_table + (j-1)*f.scan_record + ST_INDEX) == j - 1 &&
                    S[j].scan_number == f.first_scan + j - 1, 1:n), n)
    # an acquisition can hold no scans at all (e.g. an aborted run); nothing else applies
    if n == 0
        add!("info", "acquisition holds no scans", 1, 1, "first/last scan $(f.first_scan)/$(f.last_scan)")
        return R
    end
    tile = 0; prev_end = f.packet_start
    for s in S
        tile += s.packet_offset == prev_end; prev_end = s.packet_offset + s.packet_size
    end
    add!("structure", "packets tile in scan order (no gap/overlap)", tile, n)
    add!("structure", "packets end before scan table", prev_end <= f.scan_table, 1,
         "last packet ends $(f.scan_table - prev_end) bytes before the scan table")
    closes = 0; pkw = 0
    for s in S
        p = s.packet_offset
        W = Int(ru32(b, p + PKT_W_PROFILE)) + Int(ru32(b, p + PKT_W_PEAKS)) + sum(Int(ru32(b, p + k)) for k in PKT_W_REST)
        closes += PKT_HEADER + 4W == s.packet_size
        _, np = peak_list_location(b, p)
        w2 = ru32(b, p + PKT_W_PEAKS)
        pkw += w2 == 1 + f.layout.peak_words * np || (w2 == 0 && np == 0)   # W2 == 0: no peak section
    end
    add!("structure", "packet section word counts sum to packet size", closes, n)
    add!("structure", "peak-list word count == 1 + words-per-peak × n", pkw, n)

    # ---------------- physics / acquisition logic ----------------
    sorted = 0; positive = 0; inrange = 0; edge = 0; inversions = 0; maxinv = 0.0; onepeak = 0
    mz = Float32[]; it = Float32[]
    tic_dev = 0.0; tic_ok = 0; tic_exact = 0; bp_ok = 0; hdr_ok = 0
    # histogram of MS1 neighbour differences in 0.5 mDa bins over [1.00335/z - 0.05, + 0.05], z = 1..4
    nb = 200; hist = zeros(Int, 4, nb); nms1 = 0
    for (i, s) in enumerate(S)
        peaks!(mz, it, f, i)
        # A stored list may hold a stray peak out of sequence (e.g. a near-duplicate centroid of an ion stored one or
        # two places late). Each descent must be one stray peak: dropping either mz[j-1] or mz[j] restores the local
        # order; at most 5 per scan. The writer sorts these scans. The largest step back is reported.
        desc = 0; single = true; np_ = length(mz)
        for j in 2:np_
            if mz[j] < mz[j-1]
                desc += 1; inversions += 1; maxinv = max(maxinv, (mz[j-1] - mz[j]) / mz[j] * 1e6)
                prev_stray = j == 2 || mz[j] >= mz[j-2]          # mz[j-1] is the stray
                this_stray = j == np_ || mz[j+1] >= mz[j-1]      # mz[j] is the stray
                single &= prev_stray || this_stray
            end
        end
        sorted += desc <= 5 && single; onepeak += desc >= 1
        positive += all(x -> isfinite(x) && x > 0, it)
        # Edge-truncated centroids may sit just outside the stated range. Allowed: high-res up to 3 peaks within
        # 50 ppm of the edge; ion trap 1 peak within 1.0 m/z (unit-resolution peaks are ~0.5-1 wide). Counted over
        # all peaks.
        nlo = 0; nhi = 0; dlo = 0.0; dhi = 0.0
        for m in mz
            if m < s.low_mz - 1e-3; nlo += 1; dlo = max(dlo, s.low_mz - m); end
            if m > s.high_mz + 1e-3; nhi += 1; dhi = max(dhi, m - s.high_mz); end
        end
        it_ = s.analyser == 0x00
        lo_good = nlo == 0 || (nlo <= (it_ ? 1 : 3) && dlo <= (it_ ? 1.0 : 50e-6 * s.low_mz))
        hi_good = nhi == 0 || (nhi <= (it_ ? 1 : 3) && dhi <= (it_ ? 1.0 : 50e-6 * s.high_mz))
        inrange += lo_good && hi_good
        edge += (nlo + nhi > 0) && lo_good && hi_good
        # redundant fields
        if s.packet_type == 0x00000013   # profile-only ion-trap scan: compare with its profile points
            pmz, pit = profile(f, i); tot = sum(Float64, pit); rel = abs(s.tic / tot - 1)
            tic_ok += rel < 1e-5 || abs(s.tic - tot) < 0.5; tic_exact += rel < 1e-5 || abs(s.tic - tot) < 0.5
            # same profile point: m/z within 0.001 (points are 0.125 apart). The stored point spacing and base-peak
            # m/z may be single-precision values.
            m = argmax(pit)
            # an all-zero profile is an empty scan: stored base peak 0 / 0
            bp_ok += pit[m] == 0 ? (s.base_peak_intensity == 0) : (pit[m] == Float32(s.base_peak_intensity) && abs(s.base_peak_mz - pmz[m]) <= 1e-3)
        elseif !isempty(it)
            tot = sum(Float64, it); rel = abs(s.tic / tot - 1); tic_dev = max(tic_dev, rel)
            # the stored TIC may be coarsely rounded (e.g. multiples of 1/8): visible only in near-empty scans
            exact = rel < 1e-5 || abs(s.tic - tot) < 0.5
            # the stored TIC may differ from the sum by a zero-mean amount that shrinks as 1/sqrt(peaks); accepted
            # within 1e-4/sqrt(n). The strict count is reported in the note.
            tic_exact += exact
            tic_ok += exact || rel <= 1e-4 / sqrt(length(it))
            # table holds f64, peak list f32. Intensity may be stored rounded to 6 significant digits (rel. 1e-5);
            # the table m/z can be a separate estimate (within 2.5 ppm)
            m = argmax(it); bp_ok += abs(s.base_peak_intensity - it[m]) <= 1e-5 * it[m] && abs(s.base_peak_mz - mz[m]) <= max(eps(mz[m]), 2.5e-6 * mz[m])
        else
            tic_ok += s.tic == 0; bp_ok += 1
        end
        p = s.packet_offset
        hdr_ok += abs(rf32(b, p + PKT_LOW_MZ) - s.low_mz) < 1e-3 && abs(rf32(b, p + PKT_HIGH_MZ) - s.high_mz) < 1e-2
        # MS1 isotope spacing: count neighbour differences near 1.00335/z vs. a shifted control window
        if s.ms_order == 0x01 && s.analyser != 0x00 && !isempty(it)   # ion-trap MS1: isotopes not resolved
            nms1 += 1
            # strongest 10 % of peaks (dense low-intensity peaks otherwise dominate), but at least 50 peaks
            keep = max(cld(length(it), 10), min(50, length(it)))
            thr = partialsort(it, length(it) - keep + 1); sm = mz[it .>= thr]
            for j in 1:length(sm)-1
                for k in j+1:length(sm)   # separate loop: `break` must only end the inner scan
                    d = sm[k] - sm[j]; d > 1.06 && break
                    for z in 1:4
                        x = d - (ISOTOPE_SPACING / z - 0.05)
                        0 <= x < 0.1 && (hist[z, 1 + floor(Int, x / 5e-4)] += 1)
                    end
                end
            end
        end
    end
    add!("physics", "peaks sorted by m/z (isolated stray peaks allowed, <= 5 per scan)", sorted, n,
         inversions == 0 ? "" : "$onepeak scans with stray peaks ($inversions in all), largest step back $(round(maxinv, digits = 2)) ppm (written sorted)")
    add!("physics", "intensities positive and finite", positive, n)
    add!("physics", "peaks inside scan m/z range", inrange, n,
         "$edge scans have edge-truncated centroids just outside the range (high-res <= 3 within 50 ppm; ion trap 1 within 1.0)")
    rts = [s.retention_time for s in S]
    add!("physics", "retention time strictly increasing", count(>(0), diff(rts)) + 1, n,
         "RT $(round(rts[1], digits=3))-$(round(rts[end], digits=3)) min")
    hres = Tuple{Int,Bool,String,Int}[]
    for z in 1:(nms1 == 0 ? 0 : 4)
        h = hist[z, :]; pk = argmax(h); bg = sort(h)[nb ÷ 2]
        mode = ISOTOPE_SPACING / z - 0.05 + (pk - 0.5) * 5e-4
        # the peak bin must sit on 1.00335/z and stand > 10 Poisson sigma above the median bin
        # 1 mDa: strong 2+ isotope peaks partly coalesce on low-resolution analysers
        good = abs(mode - ISOTOPE_SPACING / z) <= 1e-3 && (h[pk] - bg) > 10 * sqrt(max(bg, 1))
        note = "mode $(round(mode, digits=4)) (expect $(round(ISOTOPE_SPACING/z, digits=4))), peak bin $(h[pk]) vs median bin $bg ($(round((h[pk]-bg)/sqrt(max(bg,1)), digits=1)) sigma), $nms1 MS1 scans"
        push!(hres, (z, good, note, h[pk]))   # reported after the per-peak charge test below
    end
    ms1 = findall(s -> s.ms_order == 0x01, S); ms2 = findall(s -> s.ms_order == 0x02, S)
    add!("physics", "MS order is 1 or 2", length(ms1) + length(ms2), n)
    cls = [ru32(b, f.scan_table + (j-1)*f.scan_record + ST_MSCLASS) for j in 1:n]
    if all(==(0xffffffff), cls)
        add!("physics", "scan-table +8 is 0xffffffff on every scan (no class stored)", 1, 1)
    else
        # +8 = the scan's position (event number) in its acquisition cycle. MS1 may span several events (split MS1
        # ranges), so: +8 rises by one scan to the next or restarts at 0, and within a cycle every MS1 event precedes
        # every MS2 event.
        okstep = count(j -> j == 1 || cls[j] == 0 || cls[j] == cls[j-1] + 1 || (S[j].ms_order == 0x02 && cls[j] == cls[j-1]), 1:n)
        add!("physics", "scan-table +8 (cycle event) rises within a cycle, restarts at 0", okstep, n)
        maxms1 = maximum((cls[j] for j in ms1); init = 0); minms2 = minimum((cls[j] for j in ms2); init = typemax(UInt32))
        add!("physics", "MS1 events precede MS2 events", maxms1 < minms2, 1, "MS1 events 0-$(maxms1), MS2 events from $(minms2)")
    end
    nprec = [ru32(b, c + f.region_c_fixed) for c in f.region_c_starts]
    add!("physics", "every MS2 has a precursor entry in region C", count(j -> nprec[j] >= 1, ms2), length(ms2))
    blocks = diff([ms1; n + 1]) .- 1
    modal = isempty(blocks) ? 0 : argmax(x -> count(==(x), blocks), unique(blocks))
    # DIA = the MS2 windows repeat in complete passes (Pioneer's cycle rule), not necessarily a fixed MS2 count per MS1
    # (MS1 scans may be inserted by time)
    cyc = pioneer_cycle_idxs([s.ms_order for s in S], [s.center_mz for s in S], [s.isolation_width for s in S])
    windows = sort(unique(round.([S[j].center_mz for j in ms2], digits=3)))
    passes = Dict{Int32,Vector{Float64}}()
    for j in ms2; push!(get!(passes, cyc[j], Float64[]), round(S[j].center_mz, digits=3)); end
    lastpass = isempty(passes) ? Int32(0) : maximum(keys(passes))
    full = [sort(v) == windows for (k, v) in passes if k != lastpass]
    regular = length(full) >= 3 && count(full) >= 0.9 * length(full)
    an = sort(unique([(s.ms_order, s.analyser) for s in S]))
    add!("acquisition", "analyser code is a known value (0, 4, 7)", count(s -> haskey(ANALYSER_NAMES, s.analyser), S), n,
         join(["MS$(o): $(get(ANALYSER_NAMES, a, string(a)))" for (o, a) in an], ", "))
    add!("acquisition", "acquisition type", 1, 1, regular ? "DIA: $(length(windows)) windows repeating in $(length(full)) complete passes" :
         "windows do not repeat in complete passes (targeted/DDA): DIA window checks skipped")
    if !isempty(ms2)
        lo1 = minimum(S[j].low_mz for j in ms1; init=Inf); hi1 = maximum(S[j].high_mz for j in ms1; init=-Inf)
        add!("acquisition", "MS2 isolation centres inside MS1 range", count(j -> lo1 <= S[j].center_mz <= hi1, ms2), length(ms2))
        known = filter(j -> !isnan(S[j].isolation_width), ms2)
        add!("acquisition", "MS2 isolation width positive and finite (where known)", count(j -> isfinite(S[j].isolation_width) && S[j].isolation_width > 0, known), length(known))
        tw = method_windows(f)
        if !isempty(tw)
            # every MS2 (centre, width) should be a window listed in the embedded method's table. Width: the table's
            # (within 0.1 %: tables after window placement optimization list slightly wider edges) or the method's
            # stated window + overlap
            sw_m = method_single_width(f)
            inlist(j) = any(t -> abs(t[1] - S[j].center_mz) < 1e-3 * max(1, t[1] / 1000) + 5e-4 &&
                                 (abs(t[2] - S[j].isolation_width) <= 1e-3 * t[2] + 1e-3 ||
                                  (sw_m !== nothing && abs(sw_m - S[j].isolation_width) < 0.01)), tw)
            add!("redundant", "MS2 isolation window is listed in the method's window table", count(inlist, known), length(known),
                 "$(length(tw)) windows in the method text")
        end
        mw = method_text_widths(f)
        if length(mw) == 1 && !isempty(known)
            add!("redundant", "MS2 isolation width == method text ($(mw[1]))", count(j -> abs(S[j].isolation_width - mw[1]) < 1e-3, known), length(known))
        end
        length(known) < length(ms2) && add!("info", "isolation width unknown for $(length(ms2) - length(known)) MS2 (pre-66 file without an extras width field)", 1, 1)
    end
    if !isempty(ms2) && !regular
        # DDA / targeted: each MS2 precursor centre should be a peak (10 ppm) in its MS1; control: centre + 0.3.
        # The MS1 is the triggering scan when the per-scan extras record it ("Master Scan Number"; with parallel
        # MS2 the latest MS1 is not always the trigger), otherwise the latest MS1 at the same FAIMS voltage whose
        # m/z range contains the precursor.
        lastms1 = Dict{Tuple{Float64,Float64,Float64},Vector{Float32}}(); nm = 0; hit = 0; ctl = 0; by_master = 0
        has(v, t) = (j = searchsortedfirst(v, t * (1 - 10e-6)); j <= length(v) && v[j] <= t * (1 + 10e-6))
        cvkey(x) = isnan(x) ? -1e9 : x
        ex = try scan_extras(f) catch; nothing end
        mcol = (ex === nothing || !haskey(ex.columns, "Master Scan Number")) ? nothing : ex.columns["Master Scan Number"]
        sidx = Dict(s.scan_number => i for (i, s) in enumerate(S)); ms1peaks = Dict{Int,Vector{Float32}}()
        for (i, s) in enumerate(S)
            if s.ms_order == 0x01
                lastms1[(cvkey(s.faims_cv), s.low_mz, s.high_mz)] = copy(peaks!(mz, it, f, i)[1])
            elseif isfinite(s.center_mz)
                mj = mcol === nothing ? 0 : (x = mcol[i]; x isa Real && isfinite(x) ? get(sidx, Int(x), 0) : 0)
                if mj > 0 && S[mj].ms_order == 0x01
                    v = get!(() -> copy(peaks!(mz, it, f, mj)[1]), ms1peaks, mj); by_master += 1
                else
                    ks = [k for k in keys(lastms1) if k[1] == cvkey(s.faims_cv) && k[2] <= s.center_mz <= k[3]]
                    isempty(ks) && continue
                    v = lastms1[ks[1]]
                end
                nm += 1; hit += has(v, s.center_mz); ctl += has(v, s.center_mz + 0.3)
            end
        end
        # data-dependent runs pick a new precursor nearly every MS2; targeted (PRM) runs revisit a fixed list
        distinct = length(unique(round.([S[j].center_mz for j in ms2], digits=3)))
        # a targeted method (tSIM / tMS2 / tMSn experiment or a mass list, and no data-dependent section) chooses its
        # precursors from the list, not from MS1 peaks
        mt = join(method_text(f), "\n")
        targeted_method = !isempty(mt) && !occursin("Data Dependent", mt) &&
                          occursin(r"Experiment Name = t(SIM|MS2|MSn)|Mass List Table|Targeted"i, mt)
        if nm > 0 && targeted_method
            add!("info", "targeted method (no data-dependent section): DDA precursor check skipped", 1, 1,
                 "$(round(100hit/nm, digits=1)) % of $nm MS2 precursors are MS1 peaks (control $(round(100ctl/nm, digits=1)) %)")
        elseif nm > 0 && distinct < 0.3 * length(ms2)
            add!("info", "targeted acquisition: $distinct distinct precursors over $(length(ms2)) MS2; DDA precursor check skipped", 1, 1)
        elseif nm > 0
            add!("acquisition", "DDA: MS2 precursor is a peak in its master (else latest same-voltage) MS1 (>= 80 %, >= 5x control)",
                 hit >= 0.8nm && hit >= 5ctl, 1, "$(round(100hit/nm, digits=1)) % of $nm MS2 (control +0.3: $(round(100ctl/nm, digits=1)) %; $by_master matched to their recorded master MS1)")
        end
    end
    if !isempty(ms2) && regular
        # DIA windows: a full pass (Pioneer cycle) contains each distinct window once; windows are evenly spaced,
        # equal width, overlap-or-touch (no gaps) and lie inside the MS1 range.
        add!("acquisition", "each window pass contains every window once", count(full), length(full),
             "$(length(windows)) windows; final pass ($(length(passes[lastpass])) windows) excluded as cut off by run end")
        add!("acquisition", "final pass is a subset of the windows", all(in(windows), passes[lastpass]), 1)
        # windows may be fixed- or variable-width: test with each window's own edges
        wl = sort(unique([(round(S[j].center_mz, digits=3), round(S[j].isolation_width, digits=3)) for j in ms2]))
        lows = [c - w / 2 for (c, w) in wl]; highs = [c + w / 2 for (c, w) in wl]
        ov = [highs[k] - lows[k+1] for k in 1:length(wl)-1]                 # overlap with the next window (>= 0: no gap)
        widths = [S[j].isolation_width for j in ms2]
        known_ms2 = filter(j -> !isnan(S[j].isolation_width), ms2)
        # Tiling checks compare centres/widths with a regular layout. When the method's own window table lists every
        # MS2 window, the layout is whatever the method defines (gaps/overlaps by design) and is reported as info.
        # A repeating short list with gaps wider than the windows is a targeted list, not a tiling.
        tw = method_windows(f)
        sw_t = method_single_width(f)
        by_table = !isempty(tw) && all(j -> any(t -> abs(t[1] - S[j].center_mz) < 1e-3 * max(1, t[1] / 1000) + 5e-4 &&
                                                    (abs(t[2] - S[j].isolation_width) <= 1e-3 * t[2] + 1e-3 ||
                                                     (sw_t !== nothing && abs(sw_t - S[j].isolation_width) < 0.01)), tw), known_ms2)
        sparse = !isempty(ov) && count(o -> -o > median_width(wl), ov) > length(ov) / 2
        gapstr = isempty(ov) ? "" : "neighbour overlap $(round(minimum(ov), digits=3)) to $(round(maximum(ov), digits=3)) (negative = gap)"
        touch_ok = count(>=(-5e-3), ov); const_ok = count(o -> abs(o - ov[1]) < 0.01, ov)
        if by_table && (touch_ok < length(ov) || const_ok < length(ov))
            add!("info", "window layout as listed in the method table (gaps/overlaps by design)", 1, 1, gapstr)
        elseif sparse
            add!("info", "repeating targeted list (gaps wider than the windows): tiling checks skipped", 1, 1, gapstr)
        else
            add!("acquisition", "adjacent windows overlap or touch (no gaps)", count(>=(-5e-3), ov), length(ov))   # 3-decimal rounding
            add!("acquisition", "overlap between neighbouring windows constant (±0.01)", count(o -> abs(o - ov[1]) < 0.01, ov), length(ov),
                 isempty(ov) ? "" : "overlap $(round(ov[1], digits=3))")
        end
        fixed = all(x -> abs(x[2] - wl[1][2]) < 1e-3, wl)
        add!("info", "window widths: $(fixed ? "fixed $(wl[1][2])" : "variable $(minimum(last.(wl)))-$(maximum(last.(wl)))")", 1, 1)
        sp = diff(first.(wl))
        ew = method_experiment_widths(f)
        if length(ew) == length(wl) && length(ew) > 1
            add!("redundant", "window widths (by centre) == method experiment widths, in order", count(k -> abs(ew[k] - wl[k][2]) < 1e-3, eachindex(ew)), length(ew))
        elseif isempty(ew) && method_single_width(f) !== nothing && isempty(method_windows(f))
            sw = method_single_width(f)
            add!("redundant", "window widths == method isolation window (+ stated overlap) ($(sw), ±0.01)", count(x -> abs(x[2] - sw) < 0.01, wl), length(wl))
        elseif length(unique(ew)) == 1
            add!("redundant", "window widths == the method's single isolation window ($(ew[1]))", count(x -> abs(x[2] - ew[1]) < 1e-3, wl), length(wl))
        end
        add!("acquisition", "window coverage", 1, 1,
             "windows cover $(round(minimum(lows), digits=2))-$(round(maximum(highs), digits=2)); MS1 $(lo1)-$(hi1)")
        cm = count(==(modal), blocks[1:end-1])
        add!("info", "MS2 count between MS1 scans ($(cm == length(blocks) - 1 ? "constant" : "varies"))", 1, 1,
             "$modal MS2 per MS1 in $cm/$(length(blocks) - 1) blocks; final block $(blocks[end])")
    end

    # ---------------- redundant fields agree ----------------
    add!("redundant", "scan-table TIC == sum of peak intensities (within 1e-4/sqrt(peaks))", tic_ok, n,
         "exact (rel < 1e-5 or abs < 0.5) in $tic_exact/$n; max relative deviation $(round(tic_dev, sigdigits=3))")
    add!("redundant", "scan-table base peak == highest peak", bp_ok, n)
    add!("redundant", "packet header m/z range == scan-table range", hdr_ok, n)
    # per-peak word (section 3): index bits == peak position; charge byte agrees with isotope spacing
    with3 = 0; idx_ok = 0; zt = zeros(Int, 6, 3)
    for (i, s) in enumerate(S)
        p = s.packet_offset; pk, np = peak_list_location(b, p)
        (np > 0 && ru32(b, p + PKT_W_REST[1]) == np) || continue
        with3 += 1; s3 = pk + 4 * Int64(ru32(b, p + PKT_W_PEAKS))
        u = [ru32(b, s3 + 4k) for k in 0:np-1]
        idx_ok += all(k -> (u[k+1] & PW_INDEX_MASK) == k, 0:np-1)
        s.ms_order == 0x01 || continue
        peaks!(mz, it, f, i)
        for k in 1:np
            z = Int(u[k] >> PW_CHARGE_SHIFT); 2 <= z <= 6 || continue
            near(sp) = any(j -> abs(abs(mz[j] - mz[k]) - sp) < 3e-3, max(1, k-40):min(np, k+40))
            zt[z, 1] += 1; zt[z, 2] += near(ISOTOPE_SPACING / z); zt[z, 3] += near(ISOTOPE_SPACING / (z + 1))
        end
    end
    if with3 > 0
        add!("redundant", "per-peak word: index bits == peak position", idx_ok, with3,
             "$with3 of $n scans carry the per-peak section")
        for z in 2:3
            zt[z, 1] == 0 && continue
            r = zt[z, 2] / zt[z, 1]; w = zt[z, 3] / zt[z, 1]
            note = "$(round(100r, digits=1))% at 1.00335/$z vs $(round(100w, digits=1))% at 1.00335/$(z+1); $(zt[z, 1]) peaks"
            if zt[z, 1] < 100   # near-empty MS1: too few charged peaks to test
                add!("info", "MS1 peaks with charge byte $z: too few to test", 1, 1, note)
            else
                add!("redundant", "MS1 peaks with charge byte $z have an isotope neighbour at 1.00335/$z", r >= 5w && r > 0.25, 1, note)
            end
        end
    end
    # MS1 isotope histograms. z = 2, 3 are the peptide test, except when (a) there are too few MS1 peaks to fill it
    # (peak bin < 50) or (b) a low-complexity sample lets another recurring spacing dominate while the direct per-peak
    # test for the same charge passes on >= 1,000 peaks. z = 1 can be dominated by singly charged background (~1.0006);
    # z = 4 can be rare: always info.
    for (z, good, note, pkc) in hres
        direct = with3 > 0 && zt[min(z, 6), 1] >= 1000 && zt[z, 2] >= 5 * zt[z, 3] && zt[z, 2] > 0.25 * zt[z, 1]
        if z in (2, 3) && (good || (pkc >= 50 && !direct))
            add!("physics", "MS1 neighbour-difference histogram peaks at 1.00335/$z", good, 1, note)
        elseif z in (2, 3)
            add!("info", "MS1 neighbour-difference histogram at 1.00335/$z: $(pkc < 50 ? "too few MS1 peaks" : "other spacing dominates; per-peak charge test passes")", 1, 1, note)
        else
            add!("info", "MS1 neighbour-difference histogram at 1.00335/$z ($(good ? "peak present" : "no clear peak"))", 1, 1, note)
        end
    end
    rng(c) = c + f.region_c_fixed + 4 + f.region_c_prec * Int64(ru32(b, c + f.region_c_fixed)) + 4   # first range entry
    add!("redundant", "region C first m/z range == scan-table range",
         count(j -> begin r = rng(f.region_c_starts[j]); rf64(b, r) == S[j].low_mz && rf64(b, r + 8) == S[j].high_mz end, 1:n), n)
    ms1p = filter(j -> nprec[j] >= 1, ms1)
    isempty(ms1p) || add!("redundant", "MS1 precursor entry centre ± width/2 == scan range",
         count(j -> begin pe = f.region_c_starts[j] + f.region_c_fixed + 4; ce = rf64(b, pe + PE_CENTER); w = rf64(b, pe + PE_WIDTH)
                    abs(ce - w/2 - S[j].low_mz) < 1e-6 && abs(ce + w/2 - S[j].high_mz) < 1e-6 end, ms1p), length(ms1p))
    add!("redundant", "device block first/last RT == first/last scan RT",
         (rf64(b, f.device_block + DEV_FIRST_RT) == rts[1]) + (rf64(b, f.device_block + DEV_LAST_RT) == rts[end]), 2)
    R
end

function report(io::IO, f::Acquisition, R::Vector{CheckResult})
    println(io, "# Consistency report: ", basename(f.path))
    println(io, "\n", length(f), " scans; ", count(ok, R), "/", length(R), " checks fully pass\n")
    println(io, "| group | check | pass / total | | note |")
    println(io, "|---|---|---|---|---|")
    for r in R
        println(io, "| ", r.group, " | ", r.name, " | ", r.pass, " / ", r.total, " | ", ok(r) ? "PASS" : "FAIL", " | ", r.note, " |")
    end
end
report(f::Acquisition, R) = report(stdout, f, R)
