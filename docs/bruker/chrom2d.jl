# Offline 2D chromatogram integration on the PIONEER_CHROM_DUMP_DIR weight dumps (cycle x IM slice per precursor).
# Steps mirror Pioneer's 1D path: weighted Whittaker-Henderson smoothing (2nd differences, both axes), baseline
# plane through the window's border, region grown from the apex, 2D trapezoid.
# Usage: julia --project=<Pioneer> chrom2d.jl <dump.arrow> <passing_psms.arrow> <out.arrow> [lambda_rt lambda_im]
using Arrow, DataFrames, LinearAlgebra, Statistics, Printf

const STRIDE = 8

"""
Second divided-difference penalty DᵀD for n points at coordinates x (the 1D smoother's convention: the window's
axis is normalised to 0-1 and the penalty is on divided differences, so λ is comparable across axes and windows).
"""
function d2tD2(x::AbstractVector{Float64})
    n = length(x)
    n < 3 && return zeros(n, n)
    D = zeros(n - 2, n)
    for i in 1:n-2
        h1 = x[i+1] - x[i]; h2 = x[i+2] - x[i+1]
        # second divided difference f[x_i, x_{i+1}, x_{i+2}] scaled like a second derivative
        D[i, i] = 2 / (h1 * (h1 + h2)); D[i, i+1] = -2 / (h1 * h2); D[i, i+2] = 2 / (h2 * (h1 + h2))
    end
    D' * D
end

"""
    wh2d(Y, W, λr, λi) -> Z
Weighted 2D Whittaker-Henderson on the unit square: rows = cycles (RT normalised to 0-1 over the window),
cols = IM slices (normalised to 0-1); penalties on second divided differences along each axis.
"""
function wh2d(Y::Matrix{Float64}, W::Matrix{Float64}, λr, λi)
    nr, nc = size(Y)
    xr = nr > 1 ? collect(range(0.0, 1.0, length = nr)) : [0.0]; xc = nc > 1 ? collect(range(0.0, 1.0, length = nc)) : [0.0]
    A = Diagonal(vec(W)) + λr * kron(I(nc), d2tD2(xr)) + λi * kron(d2tD2(xc), I(nr))
    reshape(Symmetric(Matrix(A)) \ (vec(W) .* vec(Y)), nr, nc)
end

"Least-squares plane a + b·r + c·c through the given cells; returns the plane evaluated on the grid."
function border_plane(Z::Matrix{Float64}, ring::Int = 1)
    nr, nc = size(Z)
    rows = Float64[]; cols = Float64[]; vals = Float64[]
    for r in 1:nr, c in 1:nc
        (r <= ring || r > nr - ring || c <= ring || c > nc - ring) || continue
        push!(rows, r); push!(cols, c); push!(vals, Z[r, c])
    end
    X = hcat(ones(length(rows)), rows, cols)
    β = X \ vals
    [β[1] + β[2] * r + β[3] * c for r in 1:nr, c in 1:nc]
end

"Region grown from the apex over cells above frac·apex (4-connected), bounded by the grid."
function grow_region(Z::Matrix{Float64}, apex::CartesianIndex, frac::Float64)
    nr, nc = size(Z); thr = frac * Z[apex]
    inreg = falses(nr, nc); stack = [apex]; inreg[apex] = true
    while !isempty(stack)
        p = pop!(stack)
        for d in (CartesianIndex(1, 0), CartesianIndex(-1, 0), CartesianIndex(0, 1), CartesianIndex(0, -1))
            q = p + d
            (1 <= q[1] <= nr && 1 <= q[2] <= nc) || continue
            inreg[q] && continue
            if Z[q] > thr && Z[q] <= Z[p] * 1.15      # keep descending (allow 15% noise), as the 1D running-min rule
                inreg[q] = true; push!(stack, q)
            end
        end
    end
    inreg
end

"2D trapezoid over the region: cells outside the region count as 0; Δrt in minutes per cycle, Δim in slices."
function trapezoid2d(Z::Matrix{Float64}, inreg::BitMatrix, drt::Float64, dim::Float64)
    nr, nc = size(Z); V = ifelse.(inreg, Z, 0.0)
    total = 0.0
    for r in 1:nr-1, c in 1:nc-1
        total += (V[r, c] + V[r+1, c] + V[r, c+1] + V[r+1, c+1]) / 4
    end
    total * drt * dim
end

function integrate_all(dumpf, psmsf; λr = 1e-6, λi = 1e-6, frac = 0.03, verbose = true)
    d = DataFrame(Arrow.Table(dumpf))
    p = DataFrame(Arrow.Table(psmsf)); p = p[(p.target) .& (p.qval .<= 0.01), :]
    best = Dict(r.precursor_idx => (scan = Int(r.scan_idx), rt = Float64(r.rt)) for r in eachrow(p))
    # cycle -> RT (minutes) map, from the dump itself
    cyc_rt = Dict{UInt32, Float64}()
    for r in eachrow(d); cyc_rt[r.cycle_idx] = Float64(r.rt); end
    out = DataFrame(precursor_idx = UInt32[], raw_sum = Float64[], area2d = Float64[], area_imsum1d = Float64[], apex_smoothed = Float64[],
                    n_cycles = Int[], n_im = Int[], region_cells = Int[], apex_on_edge = Bool[], baseline_frac = Float64[])
    for sub in groupby(d, :precursor_idx)
        pid = sub.precursor_idx[1]
        cyc = Int.(sub.cycle_idx); ims = Int.(sub.im_scan) .÷ STRIDE
        c0, c1 = extrema(cyc); i0, i1 = extrema(ims)
        nr = c1 - c0 + 1; nc = i1 - i0 + 1
        (nr >= 3 && nc >= 3) || continue
        Y = zeros(nr, nc); W = zeros(nr, nc)
        for k in eachindex(cyc); Y[cyc[k]-c0+1, ims[k]-i0+1] += Float64(sub.weight[k]); W[cyc[k]-c0+1, ims[k]-i0+1] = 1.0; end
        raw_sum = sum(Y)
        raw_sum > 0 || continue
        scale = maximum(Y); Yn = Y ./ scale
        Z = wh2d(Yn, W, λr, λi)
        B = border_plane(Z)
        Zb = max.(Z .- B, 0.0)
        apex = argmax(Zb)
        edge = apex[1] == 1 || apex[1] == nr || apex[2] == 1 || apex[2] == nc
        reg = grow_region(Zb, apex, frac)
        rts = sort(unique(cyc)); drt = length(rts) > 1 ? median(diff([cyc_rt[UInt32(c)] for c in rts])) : 1.0
        area = trapezoid2d(Zb, reg, drt, 1.0) * scale
        # IM-summed 1D: sum the region's IM columns per cycle, trapezoid along RT
        trace = [sum(ifelse.(reg[r, :], Zb[r, :], 0.0)) for r in 1:nr]
        area1d = (sum(trace) - (trace[1] + trace[end]) / 2) * drt * scale
        push!(out, (pid, raw_sum, area, area1d, Zb[apex] * scale, nr, nc, count(reg), edge, sum(max.(B, 0.0)) / max(sum(Z), 1e-9)))
    end
    verbose && @printf("%d precursors integrated; apex on window edge %.1f%%; median region cells %d; median baseline share %.2f\n",
        nrow(out), 100mean(out.apex_on_edge), median(out.region_cells), median(out.baseline_frac))
    out
end

if abspath(PROGRAM_FILE) == @__FILE__
    λr = length(ARGS) >= 5 ? parse(Float64, ARGS[4]) : 1e-6; λi = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 1e-6
    res = integrate_all(ARGS[1], ARGS[2]; λr = λr, λi = λi)
    Arrow.write(ARGS[3], res)
end
