# Module for PDF generation - loads lazily with plotting
module PDFGenerator

function create_multipage_pdf(plots::Vector, dest::String)
    # Ensure plotting is loaded first
    if !isdefined(Main.Pioneer, :PLOTTING_LOADED) || !Main.Pioneer.PLOTTING_LOADED[]
        Main.Pioneer.ensure_plotting_loaded()
    end
    
    # Now we can use Plots and GR
    Plots = Main.Pioneer.Plots
    GR = Plots.GR
    
    # Ensure GR backend is active and set wstype for PDF output
    Plots.gr()  # This will only initialize once
    
    withenv("GKSwstype" => "100") do
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

function save_multipage_pdf(plots::Vector, dest::String)
    ensure_directory_exists(dest)
    PDFGenerator.create_multipage_pdf(plots, dest)
    return dest
end