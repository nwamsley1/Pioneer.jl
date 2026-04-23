module PioneerParamsApp

using PioneerParams

main_GetBuildLibParams(argv=ARGS)::Cint = PioneerParams.main_GetBuildLibParams(argv)
main_GetParseSpecLibParams(argv=ARGS)::Cint = PioneerParams.main_GetParseSpecLibParams(argv)
main_GetSearchParams(argv=ARGS)::Cint = PioneerParams.main_GetSearchParams(argv)
julia_main()::Cint = main_GetSearchParams(ARGS)

export julia_main
export main_GetBuildLibParams
export main_GetParseSpecLibParams
export main_GetSearchParams

end
