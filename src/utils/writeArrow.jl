function writeArrow(fpath::String, df::AbstractDataFrame)
    fpath = normpath(fpath)
    if Sys.iswindows()
        tpath = tempname()
        Arrow.write(tpath, df)
        if isfile(fpath)
            run(`cmd /c del /f "$fpath"`)
        end
        mv(tpath, fpath, force = true)
    else     #If Linux/MacOS easy
        Arrow.write(fpath, df)
    end
    return nothing
end