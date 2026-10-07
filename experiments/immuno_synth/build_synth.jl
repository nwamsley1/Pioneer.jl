# EXPERIMENT (immunopeptidomics feasibility): build a synthetic fragment index = a real library's index
# fragments (real predictions, real ranks) + "fake" non-enzymatic human 8-12mer precursors (charges 1-3, <=1 Mox,
# shuffled decoys) whose fragment ranks are synthetic. Used by the PIONEER_SYNTH_DIR harness in LibrarySearch.jl
# to measure fragment-index candidate volumes and bitvec calibration at immunopeptidome library size.
#
# usage: julia --project=PIONEER_DIR build_synth.jl CONFIG.json
# CONFIG keys: out_dir, real_lib, template_lib, fasta, fraction (0..1 of unique fake peptides, nested by hash),
#   frag_model ("template" | "uniform"), widths ([5.0, ...]), id_type ("UInt16" | "UInt32" | "auto"),
#   min_len, max_len, min_charge, max_charge, prec_mz_min, prec_mz_max, frag_mz_min, frag_mz_max, seed,
#   check_index (optional path: compare a fraction-0 build against this serialized index)

using Pioneer, Arrow, DataFrames, JSON, Random, LinearAlgebra, Statistics, Printf
const P = Pioneer
const SimpleFrag = P.SimpleFrag

const C_CARB = 57.021464
const MOX = 15.99491
const MOD_RX = r"\((\d+),(\w),([^)]+)\)"

logmsg(s) = (println(Dates_now(), "  ", s); flush(stdout))
Dates_now() = string(round(time() - T0, digits = 1), "s")
const T0 = time()

# ── library loading: fragments back in rank order (as add_fragment_indexes! does) ─────────────────────────────
function load_rank_ordered_lib(lib_path::AbstractString)
    frags = P.deserialize_from_jls(joinpath(lib_path, "detailed_fragments.jls"))
    pid_to_fid = P.deserialize_from_jls(joinpath(lib_path, "precursor_to_fragment_indices.jls"))
    Threads.@threads :static for k in 1:length(pid_to_fid)-1
        lo, hi = Int(pid_to_fid[k]), Int(pid_to_fid[k+1]) - 1
        lo < hi && sort!(view(frags, lo:hi), by = P.getRank, alg = InsertionSort)
    end
    precursors = P.SetPrecursors(Arrow.Table(joinpath(lib_path, "precursors_table.arrow")))
    proteins = P.SetProteins(Arrow.Table(joinpath(lib_path, "proteins_table.arrow")))
    empty_pfi = P.LocalPartitionedFragmentIndex{Float32}(P.LocalPartition{Float32}[], Tuple{Float32,Float32}[], 0)
    if eltype(frags) <: P.SplineCompactFrag || eltype(frags) <: P.SplineDetailedFrag
        knots = P.deserialize_from_jls(joinpath(lib_path, "spline_knots.jls"))
        return P.SplineFragmentIndexLibrary(empty_pfi, empty_pfi, precursors, proteins,
            P.SplineFragmentLookup(frags, pid_to_fid, Tuple(knots)), P.OutputSchemaPolicy())
    end
    return P.FragmentIndexLibrary(empty_pfi, empty_pfi, precursors, proteins,
        P.StandardFragmentLookup(frags, pid_to_fid), P.OutputSchemaPolicy())
end

# ── residue masses with Pioneer's fixed C carbamidomethyl; mox_pos = 0 for none ──────────────────────────────
@inline res_mass(c::Char) = P.AA_to_mass[c] + (c == 'C' ? C_CARB : 0.0)
function residue_masses!(buf::Vector{Float64}, seq::AbstractString, mox_pos::Int)
    resize!(buf, length(seq))
    @inbounds for (i, c) in enumerate(seq)
        buf[i] = res_mass(c) + (i == mox_pos ? MOX : 0.0)
    end
    return buf
end
@inline function ion_mz(rm::Vector{Float64}, isy::Bool, pos::Int, z::Int)
    L = length(rm)
    s = 0.0
    if isy
        @inbounds for i in (L - pos + 1):L; s += rm[i]; end
        s += P.H2O
    else
        @inbounds for i in 1:pos; s += rm[i]; end
    end
    return (s + z * P.PROTON) / z
end
"Position of the j-th M (1-based) in seq; 0 when j == 0."
function nth_m(seq::AbstractString, j::Int)
    j == 0 && return 0
    n = 0
    for (i, c) in enumerate(seq)
        c == 'M' && (n += 1) == j && return i
    end
    error("no M #$j in $seq")
end
function parse_mox(mods)::Tuple{Bool, Int}   # (only C-carb / Mox mods?, Mox position or 0)
    ismissing(mods) && return (true, 0)
    pos = 0
    for m in eachmatch(MOD_RX, mods)
        name = m.captures[3]
        if name == "Unimod:35"
            pos != 0 && return (false, 0)
            pos = parse(Int, m.captures[1])
        elseif name != "Unimod:4"
            return (false, 0)
        end
    end
    return (true, pos)
end

# ── 1. fake peptides: unique (I/L-folded, first seen kept) valid 8-12mers, nested subsample by hash ────────────
function digest_fakes(fasta::AbstractString, min_len::Int, max_len::Int, fraction::Float64)
    prots = P.parse_fasta(fasta, "HUMAN")
    seen = Set{String}()
    peps = String[]
    n_windows = 0
    thresh = fraction >= 1.0 ? typemax(UInt64) : UInt64(floor(fraction * typemax(UInt64)))
    for e in prots
        s = P.get_sequence(e)
        n = length(s)
        for L in min_len:max_len, i in 1:(n - L + 1)
            sub = s[i:i+L-1]
            n_windows += 1
            all(c -> c in P.VALID_AAS, sub) || continue
            key = replace(sub, 'I' => 'L')
            key in seen && continue
            push!(seen, key)
            (fraction >= 1.0 || hash(key, UInt64(0x5eed)) <= thresh) && push!(peps, sub)
        end
    end
    return peps, seen, n_windows, length(seen)
end

# ── 2. decoys: Pioneer's "shuffle" (last residue fixed, others permuted), avoiding all target sequences ───────
function make_decoys(peps::Vector{String}, target_keys::Set{String}, seed::Int)
    decoys = Vector{String}(undef, length(peps))
    ok = trues(length(peps))
    used = Set{String}()
    rng = Xoshiro(seed)
    for (k, s) in enumerate(peps)
        L = length(s)
        chars = collect(s)
        found = false
        for _ in 1:20
            perm = randperm(rng, L - 1)
            d = String(vcat(chars[perm], chars[L]))
            key = replace(d, 'I' => 'L')
            if !(key in target_keys) && !(key in used)
                push!(used, key); decoys[k] = d; found = true
                break
            end
        end
        if !found
            ok[k] = false; decoys[k] = ""
        end
    end
    return decoys, ok
end

# ── 3. iRT: linear composition model fitted on the real library's targets (+ Gaussian residual noise) ─────────
const AAS = collect("ACDEFGHIKLMNPQRSTVWY")
const AA_IDX = Dict(c => i for (i, c) in enumerate(AAS))
function comp_features!(x::AbstractVector{Float64}, seq::AbstractString, has_mox::Bool)
    fill!(x, 0.0)
    for c in seq; x[AA_IDX[c]] += 1; end
    x[21] = has_mox ? 1.0 : 0.0
    x[22] = 1.0
    return x
end
function fit_irt_model(real_prec)
    seqs = P.getSequence(real_prec); mods = P.getStructuralMods(real_prec)
    irts = P.getIrt(real_prec); dec = P.getIsDecoy(real_prec)
    idx = [i for i in eachindex(irts) if !dec[i] && all(c -> haskey(AA_IDX, c), seqs[i])]
    X = Matrix{Float64}(undef, length(idx), 22); y = Vector{Float64}(undef, length(idx))
    for (r, i) in enumerate(idx)
        _, mp = parse_mox(mods[i])
        comp_features!(view(X, r, :), seqs[i], mp != 0)
        y[r] = irts[i]
    end
    β = X \ y
    resid = y .- X * β
    σ = std(resid)
    r2 = 1 - sum(abs2, resid) / sum(abs2, y .- mean(y))
    short = [r for (r, i) in enumerate(idx) if 8 <= length(seqs[i]) <= 12]
    σ_short = std(resid[short])
    return β, σ, σ_short, r2
end

# ── 4. fragment templates from a real library: index fragments of targets by (length, charge) ──────────────────
struct Templates
    isy::Vector{Bool}; pos::Vector{UInt8}; fz::Vector{UInt8}; delta::Vector{Float32}
    offsets::Vector{Int}                       # template t's ions: offsets[t]:offsets[t+1]-1
    by_lz::Dict{Tuple{Int,Int}, Vector{Int}}   # (length, prec charge) -> template ids
end
function build_templates(lib, min_len, max_len, charges)
    prec = P.getPrecursors(lib); lookup = P.getFragmentLookupTable(lib); dfr = P.getFragments(lookup)
    seqs = P.getSequence(prec); mods = P.getStructuralMods(prec); zs = P.getCharge(prec); dec = P.getIsDecoy(prec)
    isy = Bool[]; pos = UInt8[]; fz = UInt8[]; delta = Float32[]; offsets = [1]
    by_lz = Dict{Tuple{Int,Int}, Vector{Int}}()
    rm = Float64[]; n_exact = 0; n_all = 0; n_skipped = 0
    for pid in UInt32(1):UInt32(length(seqs))
        dec[pid] && continue
        s = seqs[pid]; L = length(s); z = Int(zs[pid])
        (min_len <= L <= max_len && z in charges) || continue
        ok, mp = parse_mox(mods[pid]); ok || continue
        all(c -> haskey(AA_IDX, c), s) || continue
        residue_masses!(rm, s, mp)
        nb = length(isy)
        bad = false
        P._visit_index_frags(lookup, dfr, pid, (UInt8(4), UInt8(3), false)) do rank, f
            (P.isY(f) || P.isB(f)) || (bad = true; return)
            p = Int(P.getIonPosition(f)); c = Int(P.getFragCharge(f))
            (1 <= p < L && c >= 1) || (bad = true; return)
            d = Float32(P.getMz(f) - ion_mz(rm, P.isY(f), p, c))
            push!(isy, P.isY(f)); push!(pos, UInt8(p)); push!(fz, UInt8(c)); push!(delta, d)
            n_all += 1; abs(d) < 0.01f0 && (n_exact += 1)
        end
        if bad || length(isy) == nb
            resize!(isy, nb); resize!(pos, nb); resize!(fz, nb); resize!(delta, nb); n_skipped += 1
            continue
        end
        push!(offsets, length(isy) + 1)
        push!(get!(by_lz, (L, z), Int[]), length(offsets) - 1)
    end
    logmsg(@sprintf("templates: %d precursors (%d skipped), %d ions, %.1f%% with |delta|<0.01 (plain b/y)",
        length(offsets) - 1, n_skipped, n_all, 100 * n_exact / max(n_all, 1)))
    for k in sort(collect(keys(by_lz))); logmsg("  templates L=$(k[1]) z=$(k[2]): $(length(by_lz[k]))"); end
    return Templates(isy, pos, fz, delta, offsets, by_lz)
end

# ── 5. one fake precursor's index fragments (up to 8, rank order) ─────────────────────────────────────────────
function fake_frags!(out::Vector{SimpleFrag{Float32}}, rng, rm::Vector{Float64}, z::Int, gpid::UInt32,
                     pmz::Float32, pirt::Float32, model::Symbol, T::Templates, fmin, fmax,
                     cand::Vector{Tuple{Bool,Int,Int}})
    empty!(out); L = length(rm)
    tids = model === :template ? get(T.by_lz, (L, z), nothing) : nothing
    if tids !== nothing && !isempty(tids)
        t = tids[rand(rng, 1:length(tids))]
        for i in T.offsets[t]:(T.offsets[t+1]-1)
            length(out) == 8 && break
            mz = Float32(ion_mz(rm, T.isy[i], Int(T.pos[i]), Int(T.fz[i])) + T.delta[i])
            (fmin <= mz <= fmax) || continue
            push!(out, SimpleFrag{Float32}(mz, gpid, pmz, pirt, UInt8(0), UInt8(1) << UInt8(length(out))))
        end
    else   # uniform: random eligible b/y ions (y>=4, b>=3, frag charge <= min(z,2)), random order
        empty!(cand)
        for c in 1:min(z, 2), p in 4:(L-1); push!(cand, (true, p, c)); end
        for c in 1:min(z, 2), p in 3:(L-1); push!(cand, (false, p, c)); end
        shuffle!(rng, cand)
        for (isy, p, c) in cand
            length(out) == 8 && break
            mz = Float32(ion_mz(rm, isy, p, c))
            (fmin <= mz <= fmax) || continue
            push!(out, SimpleFrag{Float32}(mz, gpid, pmz, pirt, UInt8(0), UInt8(1) << UInt8(length(out))))
        end
    end
    return out
end

function main(cfg_path)
    cfg = JSON.parsefile(cfg_path)
    out = cfg["out_dir"]; mkpath(out)
    fraction = Float64(cfg["fraction"]); model = Symbol(cfg["frag_model"]); seed = Int(get(cfg, "seed", 1844))
    min_len, max_len = Int(cfg["min_len"]), Int(cfg["max_len"])
    charges = Int(cfg["min_charge"]):Int(cfg["max_charge"])
    pmin, pmax = Float32(cfg["prec_mz_min"]), Float32(cfg["prec_mz_max"])
    fmin, fmax = Float32(cfg["frag_mz_min"]), Float32(cfg["frag_mz_max"])
    meta = Dict{String, Any}("config" => cfg, "threads" => Threads.nthreads())

    logmsg("loading real library $(cfg["real_lib"])")
    real = load_rank_ordered_lib(cfg["real_lib"])
    rprec = P.getPrecursors(real)
    sel_real = P.select_index_fragments(real)
    Nr = length(sel_real.prec_mzs)
    logmsg("real: $Nr precursors, $(length(sel_real.frags)) index fragments")

    if haskey(cfg, "check_index") && fraction == 0
        for w in cfg["widths"]
            idx = P.build_partitioned_index_from_selection(sel_real; partition_width = Float32(w),
                frag_bin_tol_ppm = 0.0f0, frag_bin_tol_mda = 2.0f0, rt_bin_tol = 3.0f0, id_type = UInt16)
            ref = P.deserialize_from_jls(cfg["check_index"])
            same = idx.n_partitions == ref.n_partitions && all(1:idx.n_partitions) do k
                a, b = idx.partitions[k], ref.partitions[k]
                a.fragments == b.fragments && a.local_to_global == b.local_to_global &&
                a.fragment_bins.lows == b.fragment_bins.lows && a.fragment_bins.first_bins == b.fragment_bins.first_bins &&
                [(r.lb, r.ub, r.first_bin, r.last_bin) for r in a.rt_bins] == [(r.lb, r.ub, r.first_bin, r.last_bin) for r in b.rt_bins]
            end
            logmsg("CHECK width=$w: rebuilt index identical to $(cfg["check_index"]): $same")
            meta["check_identical_w$(w)"] = same
        end
    end

    # real-library keys: target sequences (I/L folded) and (folded seq, mox pos, charge) of real precursors
    rseq = P.getSequence(rprec); rmods = P.getStructuralMods(rprec); rz = P.getCharge(rprec); rdec = P.getIsDecoy(rprec)
    real_target_keys = Set{String}(); real_prec_keys = Set{Tuple{String, Int, UInt8}}()
    for i in eachindex(rseq)
        rdec[i] && continue
        key = replace(rseq[i], 'I' => 'L'); push!(real_target_keys, key)
        if min_len <= length(key) <= max_len
            ok, mp = parse_mox(rmods[i]); ok && push!(real_prec_keys, (key, mp, rz[i]))
        end
    end

    logmsg("digesting $(cfg["fasta"])")
    peps, fake_keys, n_windows, n_unique = digest_fakes(cfg["fasta"], min_len, max_len, fraction)
    logmsg("fake peptides: $n_windows windows, $n_unique unique (I/L), $(length(peps)) kept at fraction=$fraction")
    union!(fake_keys, real_target_keys)
    decoys, dec_ok = make_decoys(peps, fake_keys, seed)
    fake_keys = nothing
    logmsg("decoys: $(count(dec_ok)) / $(length(peps)) (exhausted $(count(!, dec_ok)))")

    β, σ, σ_short, r2 = fit_irt_model(rprec)
    logmsg(@sprintf("iRT model: R2=%.3f residual sd=%.3f (8-12mers %.3f)", r2, σ, σ_short))
    meta["irt_model"] = Dict("beta" => β, "sd" => σ, "sd_short" => σ_short, "r2" => r2)

    # enumerate fake precursors: (peptide k, decoy?, mox variant j, charge z) with m/z in range
    fk = UInt32[]; fd = Bool[]; fv = UInt8[]; fz = UInt8[]; fmz = Float32[]; n_dup = 0
    rm = Float64[]
    for (k, s) in enumerate(peps)
        nM = count(==('M'), s)
        for isdec in (false, true)
            isdec && !dec_ok[k] && continue
            seq = isdec ? decoys[k] : s
            for j in 0:nM
                mp = nth_m(seq, j)
                residue_masses!(rm, seq, mp)
                mass = sum(rm) + P.H2O
                for z in charges
                    mz = Float32((mass + z * P.PROTON) / z)
                    (pmin <= mz <= pmax) || continue
                    if !isdec && (replace(s, 'I' => 'L'), mp, UInt8(z)) in real_prec_keys
                        n_dup += 1; continue
                    end
                    push!(fk, k); push!(fd, isdec); push!(fv, j); push!(fz, z); push!(fmz, mz)
                end
            end
        end
    end
    Nf = length(fk)
    logmsg("fake precursors: $Nf ($(count(!, fd)) targets, $(count(fd)) decoys; $n_dup dropped as real duplicates)")

    templates = model === :template ?
        build_templates(load_rank_ordered_lib(cfg["template_lib"]), min_len, max_len, charges) :
        Templates(Bool[], UInt8[], UInt8[], Float32[], [1], Dict{Tuple{Int,Int}, Vector{Int}}())
    GC.gc()

    # fake iRT + fragments (threaded; per-precursor RNG -> independent of thread count)
    firt = Vector{Float32}(undef, Nf)
    slot = Vector{SimpleFrag{Float32}}(undef, 8 * Nf); cnt = zeros(UInt8, Nf)
    n_fallback = Threads.Atomic{Int}(0)
    Threads.@threads :static for chunk in collect(Iterators.partition(1:Nf, max(1, cld(Nf, 4 * Threads.nthreads()))))
        rmb = Float64[]; buf = SimpleFrag{Float32}[]; x = zeros(22); cand = Tuple{Bool,Int,Int}[]
        for i in chunk
            rng = Xoshiro(hash((seed, i)))
            seq = fd[i] ? decoys[fk[i]] : peps[fk[i]]
            mp = nth_m(seq, Int(fv[i]))
            comp_features!(x, seq, mp != 0)
            firt[i] = Float32(dot(x, β) + σ_short * randn(rng))
            residue_masses!(rmb, seq, mp)
            model === :template && !haskey(templates.by_lz, (length(seq), Int(fz[i]))) && Threads.atomic_add!(n_fallback, 1)
            fake_frags!(buf, rng, rmb, Int(fz[i]), UInt32(Nr + i), fmz[i], firt[i], model, templates, fmin, fmax, cand)
            cnt[i] = length(buf)
            for (r, f) in enumerate(buf); slot[8 * (i - 1) + r] = f; end
        end
    end
    logmsg(@sprintf("fake fragments: mean %.2f per precursor, %d precursors fell back to uniform", mean(cnt), n_fallback[]))

    # combined selection: real (pids 1..Nr) then fakes (Nr+1..Nr+Nf)
    nfr = length(sel_real.frags) + sum(Int, cnt; init = 0)
    frags = Vector{SimpleFrag{Float32}}(undef, nfr)
    copyto!(frags, sel_real.frags)
    offsets = Vector{Int}(undef, Nr + Nf + 1); offsets[1:Nr+1] .= sel_real.offsets
    o = length(sel_real.frags)
    for i in 1:Nf
        for r in 1:cnt[i]; frags[o += 1] = slot[8 * (i - 1) + r]; end
        offsets[Nr + i + 1] = o + 1
    end
    slot = nothing; GC.gc()
    sel = P.IndexFragSelection(frags, offsets, vcat(sel_real.prec_mzs, fmz))
    logmsg("combined selection: $(Nr + Nf) precursors, $nfr fragments")

    # precursor arrays for the harness
    Arrow.write(joinpath(out, "precursors.arrow"), DataFrame(
        mz = sel.prec_mzs,
        irt = vcat(Vector{Float32}(P.getIrt(rprec)), firt),
        charge = vcat(Vector{UInt8}(rz), fz),
        is_decoy = vcat(Vector{Bool}(rdec), fd),
        is_fake = vcat(falses(Nr), trues(Nf)),
        fake_pep = vcat(zeros(UInt32, Nr), fk),
        fake_mox = vcat(zeros(UInt8, Nr), fv)))
    Arrow.write(joinpath(out, "fake_peptides.arrow"), DataFrame(target = peps, decoy = decoys))
    meta["n_real"] = Nr; meta["n_fake"] = Nf; meta["n_fake_peptides"] = length(peps)
    meta["n_fake_targets"] = count(!, fd); meta["n_fake_decoys"] = count(fd); meta["n_dup_dropped"] = n_dup
    meta["n_fragments"] = nfr; meta["fallback_uniform"] = n_fallback[]
    meta["indexes"] = Dict{String, Any}[]

    for w in cfg["widths"]
        idt = cfg["id_type"] == "auto" ? first(P.resolve_local_id_type("auto", sel.prec_mzs, Float32(w))) :
              (cfg["id_type"] == "UInt32" ? UInt32 : UInt16)
        t = @timed P.build_partitioned_index_from_selection(sel; partition_width = Float32(w),
            frag_bin_tol_ppm = 0.0f0, frag_bin_tol_mda = 2.0f0, rt_bin_tol = 3.0f0, id_type = idt)
        idx = t.value
        file = "index_w$(w)_$(idt).jls"
        P.serialize_to_jls(joinpath(out, file), idx)
        maxl = maximum(p -> Int(p.n_local_precs), idx.partitions)
        nbins = sum(p -> length(p.fragment_bins), idx.partitions)
        logmsg(@sprintf("index w=%.1f %s: %d partitions, max %d precursors/partition, %d frag bins, %.1f GB in memory, built in %.0f s",
            w, idt, idx.n_partitions, maxl, nbins, Base.summarysize(idx) / 1e9, t.time))
        push!(meta["indexes"], Dict("file" => file, "width" => w, "id_type" => string(idt),
            "n_partitions" => idx.n_partitions, "max_local" => maxl, "frag_bins" => nbins,
            "bytes" => Base.summarysize(idx), "build_s" => t.time))
        idx = nothing; GC.gc()
    end
    meta["pieced_indexes"] = Dict{String, Any}[]
    for w in get(cfg, "piece_widths", Any[])
        dir = joinpath(out, "pieces_w$(w)")
        t = @timed P.build_index_pieces(sel, dir; partition_width = Float32(w), frag_bin_tol_ppm = 0.0f0,
            frag_bin_tol_mda = 2.0f0, rt_bin_tol = 3.0f0, id_type_request = String(cfg["id_type"]),
            max_piece_bytes = round(Int, Float64(cfg["piece_max_gb"]) * 1e9))
        pfi = t.value
        logmsg(@sprintf("pieced index w=%.1f: %d pieces, sizes %s GB, built in %.0f s", w, length(pfi.pieces),
            string([round(p.bytes / 1e9, digits = 2) for p in pfi.pieces]), t.time))
        push!(meta["pieced_indexes"], Dict("dir" => dir, "width" => w, "n_pieces" => length(pfi.pieces),
            "bytes" => [p.bytes for p in pfi.pieces], "build_s" => t.time))
    end
    write(joinpath(out, "synth_meta.json"), JSON.json(meta, 2))
    logmsg("done")
end

if abspath(PROGRAM_FILE) == @__FILE__()
    main(ARGS[1])
end
