# Library-free Q1 deconvolution of a scanning-quad file into a conventional 1-Da DIA file.
# Per cycle: bin every MS2 on an m/z grid, smooth in m/z, deconvolve each grid row along Q1 with
# a triangle of half-base H bins (NNLS), emit one MS2 spectrum per 1-Da slice. Peak m/z is the
# intensity-weighted raw m/z of the row's peaks in the bins nearest the slice (instrument
# precision, not grid precision). MS1 scans pass through. Writes Pioneer's Arrow schema.
# usage: julia zt_q1_deconvolve.jl <in.arrow> <out.arrow> <rt_lo> <rt_hi> [also write raw slice to <out_raw.arrow>]
using Arrow, Statistics, Printf, Base.Threads
const IN, OUT = ARGS[1], ARGS[2]; const RTLO, RTHI = parse(Float64, ARGS[3]), parse(Float64, ARGS[4])
const OUTRAW = length(ARGS) >= 5 ? ARGS[5] : ""
const GRID = parse(Float64, get(ENV, "GRID_MDA", "5")) / 1000; const SIGMA = parse(Float64, get(ENV, "SIGMA_BINS", "1"))
const H = parse(Float32, get(ENV, "H_BINS", "6.49")); const K = parse(Int, get(ENV, "K_BINS", "6"))
const MINX = parse(Float32, get(ENV, "MIN_INT", "20"))
const MODE = get(ENV, "MODE", "nnls")
const OMEGA = parse(Float32, get(ENV, "CD_OMEGA", "1.6"))   # over-relaxation for coordinate descent
const POST_Q1 = parse(Int, get(ENV, "POST_Q1_TRI", "0"))   # after NNLS: triangle smoothing along Q1, half-base in slices (0 = off)
const POST_H = parse(Float32, get(ENV, "POST_Q1_H", "0"))   # alt: 3-tap kernel [1-1/H, 1, 1-1/H] (e.g. H=6 -> 5/6,1,5/6); overrides POST_Q1_TRI                            # nnls | centroid
const SIGQ1 = parse(Float64, get(ENV, "SIGMA_Q1", "1.5"))         # centroid mode: Gaussian σ (bins) along Q1
const TOL_PPM = parse(Float64, get(ENV, "TRACE_PPM", "10"))     # trace mode: bin-to-bin link tolerance
const TRACE_MINLEN = parse(Int, get(ENV, "TRACE_MINLEN", "3"))
const kq1 = Float32[exp(-0.5 * (d / SIGQ1)^2) for d in -4:4] ./ sum(exp(-0.5 * (d / SIGQ1)^2) for d in -4:4)          # drop deconvolved peaks below this
const t = Arrow.Table(IN); const ms = t.msOrder; const cyc = t.cycle_idx; const rt = t.retentionTime; const cm = t.centerMz
const FLO = 140.0; const FHI = 1750.0
const nrows = round(Int, (FHI - FLO) / GRID) + 1
row_of(mz) = round(Int, (mz - FLO) / GRID) + 1
const kern = Float32[exp(-0.5 * (d / SIGMA)^2) for d in -3:3] ./ sum(exp(-0.5 * (d / SIGMA)^2) for d in -3:3)
tri(d) = max(0f0, 1f0 - abs(Float32(d)) / H)
tri(d, h) = max(0f0, 1f0 - abs(Float32(d)) / Float32(h))

function nnls_row!(x, r, y, nb; iters = 100)
    # Solve only over Q1 segments where the row has signal: bin j is reached by slices s in j..j+2K,
    # so a run of non-zero bins [a,b] needs slices a..b+2K. Runs closer than 2K+1 bins are merged.
    fill!(x, 0f0); copyto!(r, 1, y, 1, nb)
    j = 1
    while j <= nb
        if y[j] == 0; j += 1; continue; end
        a = j; b = j
        while true                                    # extend run, bridging gaps < 2K+1
            nxt = b + 1
            while nxt <= nb && y[nxt] == 0; nxt += 1; end
            (nxt > nb || nxt - b > 2K + 1) && break
            b = nxt
        end
        _nnls_segment!(x, r, a, b, nb, iters)
        j = b + 1
    end
end
function _nnls_segment!(x, r, a, b, nb, iters)
    scale = max(maximum(@view r[a:b]), 1f-6)         # stop test relative to the segment's signal
    for it in 1:iters
        maxd = 0f0
        @inbounds for s in a:(b + 2K)
            num = 0f0; den = 0f0
            for j in max(a, s - 2K):min(b, s)
                w = tri(j - (s - K)); w == 0 && continue
                num += w * r[j]; den += w * w
            end
            den == 0 && continue
            xn = max(0f0, x[s] + OMEGA * num / den); d = xn - x[s]; d == 0 && continue
            for j in max(a, s - 2K):min(b, s); r[j] -= tri(j - (s - K)) * d; end
            x[s] = xn; maxd = max(maxd, abs(d) / scale)
        end
        maxd < 1f-3 && break
    end
end


# ---- trace mode: link peaks across adjacent Q1 bins (robust to correlated m/z drift), then deconvolve each trace
struct Trace; bins::Vector{Int}; mz::Vector{Float32}; int::Vector{Float32}; end
function build_traces(ms2)
    nb = length(ms2); done = Trace[]
    open_tr = Trace[]                                  # traces extendable at this bin (sorted by last m/z)
    for j in 1:nb
        s = ms2[j]; mzv = t.mz_array[s]; itv = t.intensity_array[s]
        keys = Float32[tr.mz[end] for tr in open_tr]   # open_tr sorted by last mz
        used = falses(length(open_tr)); nxt = Trace[]
        for (m, x) in zip(mzv, itv)
            (ismissing(m) || ismissing(x) || x <= 0) && continue
            m32 = Float32(m); tol = m32 * TOL_PPM * 1f-6
            lo = searchsortedfirst(keys, m32 - tol); hi = searchsortedlast(keys, m32 + tol)
            best = 0; bd = Inf32
            for i in lo:hi
                used[i] && continue
                d = abs(keys[i] - m32); d < bd && (bd = d; best = i)
            end
            if best > 0
                tr = open_tr[best]; used[best] = true
                push!(tr.bins, j); push!(tr.mz, m32); push!(tr.int, Float32(x)); push!(nxt, tr)
            else
                push!(nxt, Trace([j], [m32], [Float32(x)]))
            end
        end
        for (i, tr) in enumerate(open_tr)               # traces not extended: allow one missed bin, else close
            used[i] && continue
            if j - tr.bins[end] <= 1; push!(nxt, tr) else; length(tr.bins) >= TRACE_MINLEN && push!(done, tr); end
        end
        sort!(nxt; by = tr -> tr.mz[end]); open_tr = nxt
    end
    for tr in open_tr; length(tr.bins) >= TRACE_MINLEN && push!(done, tr); end
    done
end

# output columns
out = (mz = Vector{Vector{Union{Missing,Float32}}}(), int = Vector{Vector{Union{Missing,Float32}}}(), scanHeader = String[], scanNumber = Int32[],
       packetType = Int32[], retentionTime = Float32[], lowMz = Float32[], highMz = Float32[], TIC = Float32[],
       centerMz = Union{Missing,Float32}[], isolationWidthMz = Union{Missing,Float32}[], collisionEnergyField = Union{Missing,Float32}[],
       collisionEnergyEvField = Float32[], msOrder = UInt8[], cycle_idx = Int32[])
function push_scan!(o, mz, int, s; center = missing, width = missing)
    push!(o.mz, mz); push!(o.int, int); push!(o.scanHeader, ""); push!(o.scanNumber, Int32(length(o.mz)))
    push!(o.packetType, 0); push!(o.retentionTime, Float32(rt[s])); push!(o.lowMz, Float32(t.lowMz[s])); push!(o.highMz, Float32(t.highMz[s]))
    push!(o.TIC, Float32(sum(int; init = 0f0))); push!(o.centerMz, center); push!(o.isolationWidthMz, width)
    push!(o.collisionEnergyField, t.collisionEnergyField[s]); push!(o.collisionEnergyEvField, Float32(t.collisionEnergyEvField[s]))
    push!(o.msOrder, UInt8(ms[s])); push!(o.cycle_idx, Int32(cyc[s]))
end
raw = OUTRAW == "" ? nothing : deepcopy(out)
mutable struct GridBuf; Y::Matrix{Float32}; SM::Matrix{Float32}; SI::Matrix{Float32}; X::Matrix{Float32}; active::Vector{Int}; end
const BUF = GridBuf(zeros(Float32, 0, 0), zeros(Float32, 0, 0), zeros(Float32, 0, 0), zeros(Float32, 0, 0), Int[])

sel = [i for i in 1:length(ms) if RTLO <= rt[i] <= RTHI]
cycles = sort(unique(cyc[sel]))
println("[$MODE] deconvolving $(length(cycles)) cycles, RT $RTLO-$RTHI, grid $(GRID*1000) mDa σ=$SIGMA, H=$H K=$K, $(nthreads()) threads")
t_all = time()
for (ci, c) in enumerate(cycles)
    scans = [i for i in sel if cyc[i] == c]
    ms1 = [i for i in scans if ms[i] == 1]; ms2 = sort([i for i in scans if ms[i] == 2 && !ismissing(cm[i])])
    isempty(ms2) && continue
    nb = length(ms2); centers = Float32.(cm[ms2]); step = nb > 1 ? (centers[end] - centers[1]) / (nb - 1) : 1.02f0
    for i in ms1
        v = Float32.(coalesce.(t.mz_array[i], NaN32)); w = Float32.(coalesce.(t.intensity_array[i], 0f0))
        push_scan!(out, v, w, i); raw === nothing || push_scan!(raw, v, w, i)
    end
    if MODE in ("trace", "tracecen")
        if raw !== nothing
            for s in ms2
                push_scan!(raw, Float32.(coalesce.(t.mz_array[s], NaN32)), Float32.(coalesce.(t.intensity_array[s], 0f0)), s; center = cm[s], width = t.isolationWidthMz[s])
            end
        end
        traces = build_traces(ms2)
        ntr = length(traces)
        bufmz = [[Float32[] for _ in 1:nb] for _ in 1:Threads.maxthreadid()]; bufint = [[Float32[] for _ in 1:nb] for _ in 1:Threads.maxthreadid()]
        tt2 = sum(tri(d)^2 for d in -K:K)
        @threads for ti in 1:ntr
            local tr = traces[ti]; local tid = threadid()
            local a = tr.bins[1]; local b = tr.bins[end]; local len = b - a + 1
            local y = zeros(Float32, len); local ymz = zeros(Float32, len)
            for (bb, m, x) in zip(tr.bins, tr.mz, tr.int); y[bb - a + 1] = x; ymz[bb - a + 1] = m; end
            local mzat = function(jj)                   # m/z of the trace nearest bin jj (intensity-weighted over ±1)
                local num = 0f0; local den = 0f0
                for d in -1:1; local q = jj + d; 1 <= q <= len && y[q] > 0 && (num += y[q] * ymz[q]; den += y[q]); end
                den > 0 ? num / den : ymz[argmax(y)]
            end
            if MODE == "trace"
                local x = zeros(Float32, len + 2K); local r = zeros(Float32, len)
                nnls_row!(x, r, y, len)
                for sidx in 1:(len + 2K)
                    local xv = x[sidx]; xv >= MINX || continue
                    local j = a + sidx - K - 1; 1 <= j <= nb || continue          # slice apex bin
                    push!(bufmz[tid][j], mzat(j - a + 1)); push!(bufint[tid][j], xv)
                end
            else                                        # local maximum along the (lightly smoothed) trace
                local z = similar(y)
                for q in 1:len; local acc = 0f0; for (d, kk) in zip(-4:4, kq1); local qq = q + d; 1 <= qq <= len && (acc += y[qq] * kk); end; z[q] = acc; end
                for q in 1:len
                    z[q] > 0 || continue
                    (q > 1 && z[q-1] >= z[q]) && continue; (q < len && z[q+1] > z[q]) && continue
                    local num = 0f0
                    for d in -K:K; local qq = q + d; 1 <= qq <= len && (num += tri(d) * y[qq]); end
                    local xv = num / tt2; xv >= MINX || continue
                    push!(bufmz[tid][a + q - 1], mzat(q)); push!(bufint[tid][a + q - 1], xv)
                end
            end
        end
        for j in 1:nb
            mzs = reduce(vcat, (bufmz[tid][j] for tid in 1:Threads.maxthreadid())); ints = reduce(vcat, (bufint[tid][j] for tid in 1:Threads.maxthreadid()))
            ord = sortperm(mzs)
            push_scan!(out, mzs[ord], ints[ord], ms2[j]; center = centers[j], width = step)
        end
        ci % 10 == 0 && (@printf("  %d/%d cycles, %d traces, %.0f s\n", ci, length(cycles), ntr, time() - t_all); flush(stdout))
        continue
    end
    # ---- grid path. Buffers are (Q1 × row) so a fragment row's Q1 profile is contiguous; they are
    # allocated once and only the rows touched in the previous cycle are zeroed.
    NS = nb + 2K
    if size(BUF.Y, 1) != nb
        BUF.Y = zeros(Float32, nb, nrows); BUF.SM = zeros(Float32, nb, nrows); BUF.SI = zeros(Float32, nb, nrows)
        BUF.X = zeros(Float32, NS, nrows); empty!(BUF.active)
    else
        @threads for r in BUF.active
            fill!(@view(BUF.Y[:, r]), 0f0); fill!(@view(BUF.SM[:, r]), 0f0); fill!(@view(BUF.SI[:, r]), 0f0); fill!(@view(BUF.X[:, r]), 0f0)
        end
    end
    Y, SM, SI, X = BUF.Y, BUF.SM, BUF.SI, BUF.X
    tg = time()
    touched = [falses(nrows) for _ in 1:Threads.maxthreadid()]
    @threads for j in 1:nb
        local s = ms2[j]; local tid = threadid(); local tb = touched[tid]
        local mzv = t.mz_array[s]; local itv = t.intensity_array[s]
        @inbounds for (m, x) in zip(mzv, itv)
            (ismissing(m) || ismissing(x) || m < FLO || m > FHI) && continue
            local r = row_of(m); local x32 = Float32(x)
            SM[j, r] += x32 * Float32(m); SI[j, r] += x32
            for (d, kk) in zip(-3:3, kern); local rr = r + d; 1 <= rr <= nrows || continue; Y[j, rr] += x32 * kk; tb[rr] = true; end
        end
    end
    if raw !== nothing                        # serial: push! is not thread-safe
        for s in ms2
            push_scan!(raw, Float32.(coalesce.(t.mz_array[s], NaN32)), Float32.(coalesce.(t.intensity_array[s], 0f0)), s; center = cm[s], width = t.isolationWidthMz[s])
        end
    end
    anyt = touched[1]; for tb in touched[2:end]; anyt .|= tb; end
    active = findall(anyt); BUF.active = active
    tgrid = time() - tg; tg = time()
    if MODE == "nnls"
        @threads for r in active
            local x = @view X[:, r]; local rr = zeros(Float32, nb); local y = @view Y[:, r]
            nnls_row!(x, rr, y, nb)
        end
    else
        tt2 = sum(tri(d)^2 for d in -K:K)
        @threads for r in active
            local y = @view Y[:, r]; local z = zeros(Float32, nb)
            @inbounds for j in 1:nb
                local acc = 0f0
                for (d, kk) in zip(-4:4, kq1); local jj = j + d; 1 <= jj <= nb && (acc += y[jj] * kk); end
                z[j] = acc
            end
            @inbounds for j in 1:nb
                z[j] > 0 || continue
                (j > 1 && z[j-1] >= z[j]) && continue; (j < nb && z[j+1] > z[j]) && continue
                local num = 0f0
                for d in -K:K; local jj = j + d; 1 <= jj <= nb && (num += tri(d) * y[jj]); end
                X[j + K, r] = num / tt2
            end
        end
    end
    tsolve = time() - tg; tg = time()
    if POST_Q1 > 0 || POST_H > 0                  # absorb ±1-slice placement errors
        local PW = POST_H > 0 ? 1 : POST_Q1
        kp = POST_H > 0 ? Float32[1f0 - 1f0 / POST_H, 1f0, 1f0 - 1f0 / POST_H] : Float32[tri(d, POST_Q1) for d in -POST_Q1:POST_Q1]; kp ./= sum(kp)
        @threads for r in active
            local xr = copy(@view X[:, r]); any(>(0f0), xr) || continue
            @inbounds for sidx in 1:NS
                local acc = 0f0
                for (d, kk) in zip(-PW:PW, kp); local q = sidx + d; 1 <= q <= NS && (acc += xr[q] * kk); end
                X[sidx, r] = acc
            end
        end
    end
    # emit: sparse over active rows, per-thread per-slice buffers, then merge + sort by m/z
    bufmz = [[Float32[] for _ in 1:nb] for _ in 1:Threads.maxthreadid()]; bufint = [[Float32[] for _ in 1:nb] for _ in 1:Threads.maxthreadid()]
    @threads for r in active
        local tid = threadid()
        @inbounds for sidx in (K + 1):(K + nb)
            local xv = X[sidx, r]; xv >= MINX || continue
            local j = sidx - K; local num = 0f0; local den = 0f0
            for jj in max(1, j - 1):min(nb, j + 1); num += SM[jj, r]; den += SI[jj, r]; end
            push!(bufmz[tid][j], den > 0 ? num / den : Float32(FLO + (r - 1) * GRID)); push!(bufint[tid][j], xv)
        end
    end
    for j in 1:nb
        mzs = reduce(vcat, (bufmz[tid][j] for tid in 1:Threads.maxthreadid())); ints = reduce(vcat, (bufint[tid][j] for tid in 1:Threads.maxthreadid()))
        ord = sortperm(mzs)
        push_scan!(out, mzs[ord], ints[ord], ms2[j]; center = centers[j], width = step)
    end
    temit = time() - tg
    ci % 10 == 0 && (@printf("  %d/%d cycles, %.0f s  (last cycle: grid %.1f solve %.1f emit %.1f s, %d active rows)\n", ci, length(cycles), time() - t_all, tgrid, tsolve, temit, length(active)); flush(stdout))

end
function write(o, path)
    Arrow.write(path, (mz_array = o.mz, intensity_array = o.int, scanHeader = o.scanHeader, scanNumber = o.scanNumber, packetType = o.packetType,
        retentionTime = o.retentionTime, lowMz = o.lowMz, highMz = o.highMz, TIC = o.TIC, centerMz = o.centerMz, isolationWidthMz = o.isolationWidthMz,
        collisionEnergyField = o.collisionEnergyField, collisionEnergyEvField = o.collisionEnergyEvField, msOrder = o.msOrder, cycle_idx = o.cycle_idx))
end
write(out, OUT); raw === nothing || write(raw, OUTRAW)
n2 = count(==(2), out.msOrder); npk = mean(length(out.mz[i]) for i in 1:length(out.mz) if out.msOrder[i] == 2)
@printf("wrote %s: %d scans (%d MS2), mean %.0f peaks per deconvolved MS2, %.0f s\n", OUT, length(out.mz), n2, npk, time() - t_all)
