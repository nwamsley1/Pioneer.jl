function main_pioneer()::Cint
    if isempty(ARGS)
        println("Usage: pioneer <subcommand> [args...]")
        println("Available subcommands: searchdia, buildspeclib, parsespeclib, getsearchparams, getbuildlibparams, convertmzml")
        return 1
    end
    cmd = lowercase(ARGS[1])
    args = ARGS[2:end]
    try
        if cmd == "searchdia"
            length(args) == 1 || error("SearchDIA requires one argument")
            SearchDIA(args[1])
        elseif cmd == "buildspeclib"
            length(args) == 1 || error("BuildSpecLib requires one argument")
            BuildSpecLib(args[1])
        elseif cmd == "parsespeclib"
            length(args) == 1 || error("ParseSpecLib requires one argument")
            ParseSpecLib(args[1])
        elseif cmd == "getsearchparams"
            length(args) >= 3 || error("GetSearchParams requires at least three arguments")
            params_path = length(args) >= 4 ? args[4] : missing
            GetSearchParams(args[1], args[2], args[3]; params_path=params_path)
        elseif cmd == "getbuildlibparams"
            length(args) >= 3 || error("GetBuildLibParams requires at least three arguments")
            params_path = length(args) >= 4 ? args[4] : missing
            GetBuildLibParams(args[1], args[2], args[3]; params_path=params_path)
        elseif cmd == "convertmzml"
            length(args) >= 1 || error("convertMzML requires at least one argument")
            skip_scan_header = length(args) >= 2 ? parse(Bool, args[2]) : true
            convertMzML(args[1]; skip_scan_header=skip_scan_header)
        else
            println("Unknown subcommand: $cmd")
            return 1
        end
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end
