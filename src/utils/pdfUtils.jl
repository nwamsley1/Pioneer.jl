module PDFGenerator
    import GR

    function create_multipage_pdf(plots::Vector, dest::String)
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

function save_multipage_pdf(plots::Vector{Plots.Plot}, dest::String)
    ensure_directory_exists(dest)
    PDFGenerator.create_multipage_pdf(plots, dest)
    return dest
end
