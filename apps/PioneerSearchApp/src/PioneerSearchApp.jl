module PioneerSearchApp

using PioneerSearch

const SEARCH_APP_NAME = PioneerSearch.SEARCH_APP_NAME

GetSearchParams(args...; kwargs...) = PioneerSearch.GetSearchParams(args...; kwargs...)
SearchDIA(args...; kwargs...) = PioneerSearch.SearchDIA(args...; kwargs...)
main_GetSearchParams(argv=ARGS)::Cint = PioneerSearch.main_GetSearchParams(argv)
main_SearchDIA(argv=ARGS)::Cint = PioneerSearch.main_SearchDIA(argv)
julia_main()::Cint = main_SearchDIA(ARGS)

export SEARCH_APP_NAME
export GetSearchParams
export SearchDIA
export julia_main
export main_GetSearchParams
export main_SearchDIA

end
