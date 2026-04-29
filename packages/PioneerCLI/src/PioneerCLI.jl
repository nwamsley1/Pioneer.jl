module PioneerCLI

const PIONEER_CONVERT_PKGID = Base.PkgId(Base.UUID("d0bd7c82-dff2-47fc-80b8-8504b1185c05"), "PioneerConvert")
const PIONEER_PARAMS_PKGID = Base.PkgId(Base.UUID("fc87d0e0-10c0-40c1-9f6d-0d6a27f22e8d"), "PioneerParams")
const PIONEER_PREDICT_PKGID = Base.PkgId(Base.UUID("4cddfed8-3fb1-4302-8bc9-cb29ad99671b"), "PioneerPredict")
const PIONEER_SEARCH_PKGID = Base.PkgId(Base.UUID("151a6990-cb6c-4706-a64b-c15e62217f67"), "PioneerSearch")

function _require_module(pkgid::Base.PkgId)
    Base.require(pkgid)
    mod = get(Base.loaded_modules, pkgid, nothing)
    mod === nothing && error("Failed to load $(pkgid.name)")
    return mod
end

_pioneer_convert() = _require_module(PIONEER_CONVERT_PKGID)
_pioneer_params() = _require_module(PIONEER_PARAMS_PKGID)
_pioneer_predict() = _require_module(PIONEER_PREDICT_PKGID)
_pioneer_search() = _require_module(PIONEER_SEARCH_PKGID)

function _call_entrypoint(mod::Module, name::Symbol, argv)::Cint
    return Base.invokelatest(getfield(mod, name), argv)::Cint
end

main_GetSearchParams(argv=ARGS)::Cint = _call_entrypoint(_pioneer_params(), :main_GetSearchParams, argv)
main_GetBuildLibParams(argv=ARGS)::Cint = _call_entrypoint(_pioneer_params(), :main_GetBuildLibParams, argv)
main_GetParseSpecLibParams(argv=ARGS)::Cint = _call_entrypoint(_pioneer_params(), :main_GetParseSpecLibParams, argv)
main_BuildSpecLib(argv=ARGS)::Cint = _call_entrypoint(_pioneer_predict(), :main_BuildSpecLib, argv)
main_SearchDIA(argv=ARGS)::Cint = _call_entrypoint(_pioneer_search(), :main_SearchDIA, argv)
main_convertMzML(argv=ARGS)::Cint = _call_entrypoint(_pioneer_convert(), :main_convertMzML, argv)

export main_GetSearchParams
export main_GetBuildLibParams
export main_GetParseSpecLibParams
export main_BuildSpecLib
export main_SearchDIA
export main_convertMzML

end
