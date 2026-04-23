module PioneerParamsApp

using PioneerPredict
using PioneerSearch

main_GetBuildLibParams(argv=ARGS)::Cint = PioneerPredict.main_GetBuildLibParams(argv)
main_GetParseSpecLibParams(argv=ARGS)::Cint = PioneerSearch.main_GetParseSpecLibParams(argv)
main_GetSearchParams(argv=ARGS)::Cint = PioneerSearch.main_GetSearchParams(argv)
julia_main()::Cint = main_GetSearchParams(ARGS)

export julia_main
export main_GetBuildLibParams
export main_GetParseSpecLibParams
export main_GetSearchParams

end
