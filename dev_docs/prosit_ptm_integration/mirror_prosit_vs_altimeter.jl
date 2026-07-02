# Prosit vs Altimeter mirror-plot sanity check (Phase 0 / Tier-1 validation).
#
# For each peptide+charge, overlay predicted fragment spectra on a shared m/z axis:
#   - Prosit  (InstrumentAgnosticModel)  as UP stems   (blue)
#   - Altimeter (SplineCoefficientModel) as DOWN stems (red), spline evaluated at NCE
# swept over an NCE grid. Each spectrum is self-normalized to its base peak so the
# SHAPES are comparable regardless of absolute intensity scale. This is a visual
# sanity check that the re-added Prosit path predicts sensible spectra relative to
# the model Pioneer ships (Altimeter).
#
# Usage:  julia --project=. dev_docs/prosit_ptm_integration/mirror_prosit_vs_altimeter.jl
# Output: dev_docs/prosit_ptm_integration/mirror_prosit_vs_altimeter.pdf

using Pioneer, DataFrames, Plots
gr()

const PROSIT = Pioneer.InstrumentAgnosticModel("prosit_2020_hcd")
const ALT    = Pioneer.SplineCoefficientModel("altimeter")

# Representative, well-behaved tryptic peptides (unmodified; Phase 0 is non-PTM).
const PEPTIDES = ["LGGNEQVTR", "GISNEGQNASIK", "TASEFDSAIAQDK", "LFLQFGAQGSPFLK"]
const CHARGE   = Int32(2)
const NCE_GRID = Float32[22, 26, 30]
const MIN_INT  = 0.01f0   # drop tiny peaks for legibility

# --- Prosit: one request per NCE (collision energy is a per-precursor input) ---
function prosit_spectra(nce::Float32)
    df = DataFrame(koina_sequence = PEPTIDES,
                   precursor_charge = fill(CHARGE, length(PEPTIDES)),
                   collision_energy = fill(nce, length(PEPTIDES)))
    batch = Pioneer.prepare_koina_batch(PROSIT, df; batch_size = 1000)[1]
    resp  = Pioneer.make_koina_request(batch, Pioneer.KOINA_URLS["prosit_2020_hcd"])
    res   = Pioneer.parse_koina_batch(PROSIT, resp)
    F = res.frags_per_precursor
    d = res.fragments
    # slice per precursor
    return [(mz = d.mz[(i-1)*F+1:i*F], int = d.intensities[(i-1)*F+1:i*F])
            for i in 1:length(PEPTIDES)]
end

# --- Altimeter: one request; evaluate the spline at each NCE ---
function altimeter_coeffs()
    df = DataFrame(koina_sequence = PEPTIDES,
                   precursor_charge = fill(CHARGE, length(PEPTIDES)))
    batch = Pioneer.prepare_koina_batch(ALT, df; batch_size = 1000)[1]
    resp  = Pioneer.make_koina_request(batch, Pioneer.KOINA_URLS["altimeter"])
    res   = Pioneer.parse_koina_batch(ALT, resp)
    return res  # KoinaBatchResult{Vector{Float32}} (extra_data = knots)
end

function altimeter_spectrum(res, pep_idx::Int, nce::Float32)
    F = res.frags_per_precursor
    d = res.fragments
    knots = Tuple(res.extra_data)
    rng = (pep_idx-1)*F+1 : pep_idx*F
    mz  = d.mz[rng]
    coefs = d.coefficients[rng]
    inten = Float32[max(0f0, Pioneer.splevl(nce, knots, c, 3)) for c in coefs]
    return (mz = mz, int = inten)
end

norm_int(v) = (m = maximum(v; init = 0f0); m > 0 ? v ./ m : v)

function mirror_panel(pep, nce, ps, as)
    p = plot(legend = false, title = "$pep  z$CHARGE  NCE $(Int(nce))",
             titlefontsize = 7, xlabel = "m/z", ylabel = "rel. int",
             xguidefontsize = 6, yguidefontsize = 6, tickfontsize = 5,
             ylims = (-1.15, 1.15), grid = false)
    # Prosit up
    pm, pi = ps.mz, norm_int(ps.int)
    keep = pi .> MIN_INT
    sticks!(p, pm[keep], pi[keep], color = :steelblue, linewidth = 1)
    # Altimeter down
    am, ai = as.mz, norm_int(as.int)
    keepa = ai .> MIN_INT
    sticks!(p, am[keepa], -ai[keepa], color = :firebrick, linewidth = 1)
    hline!(p, [0.0], color = :black, linewidth = 0.5)
    return p
end

println("Requesting Altimeter coefficients...")
alt = altimeter_coeffs()
println("Requesting Prosit spectra across NCE grid $NCE_GRID ...")
prosit_by_nce = Dict(nce => prosit_spectra(nce) for nce in NCE_GRID)

panels = Plots.Plot[]
for (pi_, pep) in enumerate(PEPTIDES)
    for nce in NCE_GRID
        ps = prosit_by_nce[nce][pi_]
        as = altimeter_spectrum(alt, pi_, nce)
        push!(panels, mirror_panel(pep, nce, ps, as))
    end
end

fig = plot(panels...; layout = (length(PEPTIDES), length(NCE_GRID)),
           size = (300 * length(NCE_GRID), 200 * length(PEPTIDES)),
           plot_title = "Prosit (up, blue) vs Altimeter (down, red) — self-normalized",
           plot_titlefontsize = 10)
out = joinpath(@__DIR__, "mirror_prosit_vs_altimeter.pdf")
savefig(fig, out)
println("Wrote $out")
