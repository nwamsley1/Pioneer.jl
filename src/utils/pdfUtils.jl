module PDFGenerator
import ..Plots
const GR = Plots.GR

function create_multipage_pdf(plots::Vector, dest::String)
    # Force a non-interactive GR workstation before backend initialization so
    # headless environments do not attempt to open gksqt/Qt windows.
    withenv("GKSwstype" => "100", "GKS_WSTYPE" => "100") do
        Plots.gr()  # This will only initialize once
        GR.beginprint(dest)
        try
            for p in plots
                # Suppress all display output
                redirect_stdout(devnull) do
                    redirect_stderr(devnull) do
                        display(p)
                    end
                end
            end
        finally
            GR.endprint()
        end
    end
end
end

function save_multipage_pdf(plots::AbstractVector{<:Plots.Plot}, dest::String)
    ensure_directory_exists(dest)
    PDFGenerator.create_multipage_pdf(plots, dest)
    return dest
end
