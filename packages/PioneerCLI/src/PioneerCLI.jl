__precompile__(false)

module PioneerCLI

import PioneerConvert
import PioneerParams
import PioneerPredict
import PioneerSearch

main_GetSearchParams(argv=ARGS)::Cint = PioneerParams.main_GetSearchParams(argv)
main_GetBuildLibParams(argv=ARGS)::Cint = PioneerParams.main_GetBuildLibParams(argv)
main_GetParseSpecLibParams(argv=ARGS)::Cint = PioneerParams.main_GetParseSpecLibParams(argv)
main_BuildSpecLib(argv=ARGS)::Cint = PioneerPredict.main_BuildSpecLib(argv)
main_SearchDIA(argv=ARGS)::Cint = PioneerSearch.main_SearchDIA(argv)
main_convertMzML(argv=ARGS)::Cint = PioneerConvert.main_convertMzML(argv)

export main_GetSearchParams
export main_GetBuildLibParams
export main_GetParseSpecLibParams
export main_BuildSpecLib
export main_SearchDIA
export main_convertMzML

end
