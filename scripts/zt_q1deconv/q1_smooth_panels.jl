# Raw / smoothed / deconvolved 2D grids (fragment m/z × Q1 centre) for a sub-region of one cycle.
using Arrow, Statistics, Plots
const RAW = "/Users/nathanwamsley/SciexZT_5min/arrow_a1/LFQ_ZenoTOF8600_ZTScanDIA_5Da_5min_50ng_Condition_A_REP1.arrow"
const OUT = "/Users/nathanwamsley/SciexZT_5min/plots"
const t = Arrow.Table(RAW); const ms = t.msOrder; const cyc = t.cycle_idx; const rt = t.retentionTime; const cm = t.centerMz
const FLO, FHI, GRID = parse(Float64, get(ENV, "FLO", "608.0")), parse(Float64, get(ENV, "FHI", "609.0")), 0.005
const Q1LO, Q1HI, SHOW_LO, SHOW_HI = 585.0, 655.0, 600.0, 640.0     # compute wide, display narrow (edge effects)
const H = 6.49f0; const K = 6; const QSM = parse(Int, get(ENV, "Q1_TRI", "3"))   # Q1 smoothing triangle half-base (bins)
i7 = findfirst(i -> ms[i] == 2 && rt[i] >= 7.0, 1:length(ms)); const c = cyc[i7]
scans = [i for i in 1:length(ms) if ms[i] == 2 && cyc[i] == c && !ismissing(cm[i]) && Q1LO <= cm[i] <= Q1HI]
centers = Float64.(cm[scans]); nb = length(scans)
nr = round(Int, (FHI - FLO) / GRID) + 1
tri(d, h) = max(0f0, 1f0 - abs(Float32(d)) / Float32(h))

function rawgrid()
    M = zeros(Float32, nr, nb)
    for (j, s) in enumerate(scans), (m, x) in zip(t.mz_array[s], t.intensity_array[s])
        (ismissing(m) || ismissing(x) || m < FLO || m > FHI) && continue
        M[round(Int, (m - FLO) / GRID) + 1, j] += Float32(x)
    end
    M
end
function smooth_mz(M, sigma)
    kern = [exp(-0.5 * (d / sigma)^2) for d in -3:3]; kern ./= sum(kern)
    W = zeros(Float32, size(M))
    for j in 1:nb, r in 1:nr
        M[r, j] == 0 && continue
        for (d, kk) in zip(-3:3, kern); rr = r + d; 1 <= rr <= nr && (W[rr, j] += M[r, j] * kk); end
    end
    W
end
function smooth_q1(M, hb)               # triangular kernel of half-base hb bins, normalised
    ks = [tri(d, hb) for d in -hb:hb]; ks ./= sum(ks)
    W = zeros(Float32, size(M))
    for r in 1:nr, j in 1:nb
        M[r, j] == 0 && continue
        for (d, kk) in zip(-hb:hb, ks); jj = j + d; 1 <= jj <= nb && (W[r, jj] += M[r, j] * kk); end
    end
    W
end
# generic kernel: column s of A places kern[] centred at bin s-K (kern is symmetric, length 2L+1)
function nnls_row(y, kern, L)
    ns = nb + 2K; x = zeros(Float32, ns); r = copy(y); scale = max(maximum(y), 1f-6)
    for it in 1:200
        maxd = 0f0
        for s in 1:ns
            c0 = s - K; num = 0f0; den = 0f0
            for d in -L:L; j = c0 + d; 1 <= j <= nb || continue; a = kern[d + L + 1]; num += a * r[j]; den += a * a; end
            den == 0 && continue
            xn = max(0f0, x[s] + num / den); dd = xn - x[s]; dd == 0 && continue
            for d in -L:L; j = c0 + d; 1 <= j <= nb || continue; r[j] -= kern[d + L + 1] * dd; end
            x[s] = xn; maxd = max(maxd, abs(dd) / scale)
        end
        maxd < 1f-3 && break
    end
    x[K+1:K+nb]                          # interior slices only, aligned with bins
end
function deconv(M, kern, L)
    X = zeros(Float32, size(M))
    Threads.@threads for r in 1:nr
        maximum(@view M[r, :]) > 0 || continue
        X[r, :] .= nnls_row(Vector{Float32}(@view M[r, :]), kern, L)
    end
    X
end
# kernels
Kt = Float32[tri(d, H) for d in -K:K]                                       # measured triangle
qs = Float32[tri(d, QSM) for d in -QSM:QSM]; qs ./= sum(qs)
Lc = K + QSM; Kc = zeros(Float32, 2Lc + 1)                                   # composite = triangle ⊛ Q1-smoother
for (i, a) in enumerate(-K:K), (j, b) in enumerate(-QSM:QSM); Kc[a + b + Lc + 1] += Kt[i] * qs[j]; end

R = rawgrid(); S1 = smooth_mz(R, 1.0); S2 = smooth_q1(S1, QSM)
D1 = deconv(S1, Kt, K); D2 = deconv(S2, Kc, Lc); D3 = deconv(S2, Kt, K)
sel = [j for j in 1:nb if SHOW_LO <= centers[j] <= SHOW_HI]
xs = centers[sel]; ys = [FLO + (r - 1) * GRID for r in 1:nr]
lg(M) = log10.(max.(M[:, sel], 1f0))
panels = [(R, "1. raw grid (5 mDa)"), (S1, "2. + m/z Gaussian σ=1 row"), (S2, "3. + Q1 triangle half-base $QSM bins"),
          (D1, "4. NNLS on (2), measured triangle kernel"), (D2, "5. NNLS on (3), composite kernel"), (D3, "6. NNLS on (3), original kernel (mismatched)")]
ps = [heatmap(xs, ys, lg(M); title = ttl, titlefontsize = 9, xlabel = "Q1 centre m/z", ylabel = "fragment m/z", clims = (0, 5), c = :viridis, colorbar = false)
      for (M, ttl) in panels]
mkpath(OUT)
plot(ps...; layout = (2, 3), size = (1800, 1000), margin = 5Plots.mm)
savefig(joinpath(OUT, "q1_smooth_panels_cycle$(c)_tri$(QSM)_mz$(FLO)-$(FHI).png"))
nz(M) = count(>(20f0), M[:, sel])
println("cycle $c, Q1 $(SHOW_LO)-$(SHOW_HI), m/z $FLO-$FHI: nonzero(>20) cells: raw=$(nz(R)) mz-smooth=$(nz(S1)) mz+q1-smooth=$(nz(S2)) | deconv: (2)→$(nz(D1)) (3,composite)→$(nz(D2)) (3,mismatched)→$(nz(D3))")
# strongest rows: show Q1 profiles before/after for a few
top = sortperm([maximum(R[r, sel]) for r in 1:nr]; rev = true)[1:3]
for r in top
    println("row m/z $(round(ys[r], digits = 3)):")
    println("  raw     ", join(round.(Int, R[r, sel]), " "))
    println("  smooth  ", join(round.(Int, S2[r, sel]), " "))
    println("  dec(2)  ", join(round.(Int, D1[r, sel]), " "))
    println("  dec(3c) ", join(round.(Int, D2[r, sel]), " "))
end
