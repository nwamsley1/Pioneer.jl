module PioneerPredictApp

using PioneerPredict

const PREDICT_APP_NAME = PioneerPredict.PREDICT_APP_NAME

BuildSpecLib(args...; kwargs...) = PioneerPredict.BuildSpecLib(args...; kwargs...)
GetBuildLibParams(args...; kwargs...) = PioneerPredict.GetBuildLibParams(args...; kwargs...)
main_BuildSpecLib(argv=ARGS)::Cint = PioneerPredict.main_BuildSpecLib(argv)
main_GetBuildLibParams(argv=ARGS)::Cint = PioneerPredict.main_GetBuildLibParams(argv)
julia_main()::Cint = main_BuildSpecLib(ARGS)

export PREDICT_APP_NAME
export BuildSpecLib
export GetBuildLibParams
export julia_main
export main_BuildSpecLib
export main_GetBuildLibParams

end
