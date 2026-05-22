module PDFGenerator
import ..Plots
end

# QC plots: rasterized PNG output.
#
# We previously wrote multipage vector PDFs via GR.beginprint / GR.endprint.
# Vector-PDF size scales with the number of graphical primitives; histogram-
# heavy diagnostic plots routinely produced 1+ MB per file.
#
# For QC purposes (visual inspection, not publication-grade figures) PNG is
# fine. PNG at the default Plots.jl size (900×600 here, or whatever the plot
# specified) is ~30-150 KB and renders identically at the diagnostic level.
#
# The save_multipage_pdf API is kept for back-compat at call sites: if the
# destination ends in ".pdf" we strip the extension and write per-page PNGs
# next to it. A single-page input produces "{base}.png"; multipage produces
# "{base}_p1.png", "{base}_p2.png", ...

function save_multipage_pdf(plots::AbstractVector{<:Plots.Plot}, dest::String)
    ensure_directory_exists(dest)
    base = endswith(dest, ".pdf") ? dest[1:end-4] : dest
    if length(plots) == 0
        return dest
    elseif length(plots) == 1
        Plots.savefig(plots[1], base * ".png")
    else
        for (i, p) in enumerate(plots)
            Plots.savefig(p, base * "_p$(i).png")
        end
    end
    return dest
end
