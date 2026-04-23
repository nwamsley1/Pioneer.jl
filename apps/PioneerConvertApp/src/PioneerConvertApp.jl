module PioneerConvertApp

using PioneerConvert

const CONVERT_MZML_APP_NAME = PioneerConvert.CONVERT_MZML_APP_NAME

convertMzML(args...; kwargs...) = PioneerConvert.convertMzML(args...; kwargs...)
main_convertMzML(argv=ARGS)::Cint = PioneerConvert.main_convertMzML(argv)
julia_main()::Cint = main_convertMzML(ARGS)

export CONVERT_MZML_APP_NAME
export convertMzML
export julia_main
export main_convertMzML

end
