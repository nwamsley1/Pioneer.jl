@echo off
setlocal enabledelayedexpansion

if "%JULIA_NUM_THREADS%"=="" set JULIA_NUM_THREADS=auto

set SCRIPT_DIR=%~dp0
if not defined PIONEER_VERSION (
    if exist "%SCRIPT_DIR%VERSION" (
        set /p PIONEER_VERSION=<"%SCRIPT_DIR%VERSION"
    ) else (
        set "PROJECT_TOML=%SCRIPT_DIR%..\..\..\Project.toml"
        if exist "!PROJECT_TOML!" (
            for /f "tokens=2 delims=\" %%V in ('findstr /R /C:"^version *= *\".*\"" "!PROJECT_TOML!"') do (
                set "PIONEER_VERSION=%%V"
                goto version_ready
            )
        )
    )
)
if not defined PIONEER_VERSION set "PIONEER_VERSION=unknown"

:version_ready

set SUBCOMMAND=
set SUBCOMMAND_ARGS=
set VALID_COMMANDS=search predict params-search params-predict convert-raw convert-mzml
rem empirical params-empirical 

:parse_args
if "%~1"=="" goto check_subcommand
if "%~1"=="--threads" (
    if "%~2"=="" (
        echo Error: --threads requires a value
        exit /b 1
    )
    set JULIA_NUM_THREADS=%~2
    shift
    shift
    goto parse_args
)
for %%H in (--help -h -help /?) do if /I "%~1"=="%%H" goto handle_help
for %%V in (--version -V) do if /I "%~1"=="%%V" goto handle_version

rem Check for --threads=value format
echo %~1 | findstr /C:"--threads=" >nul
if !errorlevel! equ 0 (
    set THREADS_ARG=%~1
    set JULIA_NUM_THREADS=!THREADS_ARG:--threads=!
    shift
    goto parse_args
)

rem Check for unknown options before subcommand; otherwise pass through
echo %~1 | findstr /R "^-" >nul
if !errorlevel! equ 0 (
    if "%SUBCOMMAND%"=="" (
        echo Error: Unknown option %~1
        echo Use --help for usage information
        exit /b 1
    ) else (
        if "%SUBCOMMAND_ARGS%"=="" (
            set SUBCOMMAND_ARGS=%~1
        ) else (
            set SUBCOMMAND_ARGS=%SUBCOMMAND_ARGS% %~1
        )
        shift
        goto parse_args
    )
)

rem First non-option argument should be subcommand
if "%SUBCOMMAND%"=="" (
    rem Check if it's a valid subcommand
    echo %VALID_COMMANDS% | findstr /C:"%~1" >nul
    if !errorlevel! equ 0 (
        set SUBCOMMAND=%~1
    ) else (
        echo Error: Unknown subcommand '%~1'
        echo Valid subcommands: %VALID_COMMANDS%
        echo Use --help for usage information
        exit /b 1
    )
) else (
    rem All remaining arguments go to the subcommand
    if "%SUBCOMMAND_ARGS%"=="" (
        set SUBCOMMAND_ARGS=%~1
    ) else (
        set SUBCOMMAND_ARGS=%SUBCOMMAND_ARGS% %~1
    )
)
shift
goto parse_args

:handle_help
if "%SUBCOMMAND%"=="" (
    goto show_help
) else (
    if "%SUBCOMMAND_ARGS%"=="" (
        set SUBCOMMAND_ARGS=%~1
    ) else (
        set SUBCOMMAND_ARGS=%SUBCOMMAND_ARGS% %~1
    )
    shift
    goto parse_args
)


:show_help
echo Pioneer - Mass Spectrometry Data Analysis
echo Version: %PIONEER_VERSION%
echo.
echo Usage: pioneer [options] ^<subcommand^> [subcommand-args...]
echo.
echo Options:
echo   --threads N        Set number of Julia threads (default: auto)
echo   --threads=N        Alternative syntax for setting threads
echo   --help, -h         Show this help message
echo   --version, -V      Show Pioneer version
echo.
echo Subcommands:
echo   search ^<params.json^>                       Perform DIA search analysis
echo   predict ^<params.json^>                      Predict spectral library
rem echo   empirical ^<params.json^>                    Parse spectral library
echo   params-search ^<library_path^> ^<ms_data_path^> ^<results_path^> [--params-path ^<params_out_path^>]
echo                                                 Generate search parameter template
echo                                                 Default params output path: ./search_parameters.json
echo   params-predict ^<library_outpath^> ^<fasta_path^> [--params-path ^<params_out_path^>]
echo                                                 Generate library build parameter template
echo                                                 Default params output path: ./buildspeclib_params.json
rem echo   params-empirical ^<empirical_lib_path^> ^<library_outpath^> [--params-path ^<params_out_path^>]
rem echo                                                 Generate parse parameter template
rem echo                                                 Default params output path: ./parsespeclib_params.json
echo   convert-raw ^<data_path^> [--output-dir ^<dir^>] [--skip-existing] [--concurrent-files N] [--threads-per-file N]
echo                                                 Convert Thermo RAW files via PioneerConverter
echo                                                 For additional converter options, run: pioneer convert-raw --help
echo   convert-mzml ^<data_path^> [options]
echo                                                 Convert mzML files
echo.
echo Examples:
echo   pioneer params-predict yeast.poin fasta/ --params-path predict_params.json
echo   pioneer predict predict_params.json
echo   pioneer params-search yeast.poin data/ results/ --params-path search_params.json
echo   pioneer search search_params.json                    # Use auto threading
echo   pioneer --threads=8 search search_params.json        # Use 8 threads
echo.
echo For subcommand-specific help:
echo   pioneer ^<subcommand^> --help
exit /b 0

:handle_version
if "%SUBCOMMAND%"=="" (
    goto show_version
) else (
    if "%SUBCOMMAND_ARGS%"=="" (
        set SUBCOMMAND_ARGS=%~1
    ) else (
        set SUBCOMMAND_ARGS=%SUBCOMMAND_ARGS% %~1
    )
    shift
    goto parse_args
)

:show_version
echo Pioneer version: %PIONEER_VERSION%
exit /b 0

:check_subcommand
if "%SUBCOMMAND%"=="" (
    echo Error: Subcommand required
    echo Valid subcommands: %VALID_COMMANDS%
    echo Use --help for usage information
    exit /b 1
)

rem Auto-set GC threads to half of user threads (mark threads, 1 sweep thread)
if "%JULIA_NUM_THREADS%"=="auto" (
    set /a _THREAD_COUNT=%NUMBER_OF_PROCESSORS%
) else (
    set /a _THREAD_COUNT=%JULIA_NUM_THREADS%
)
set /a _GC_MARK_THREADS=(_THREAD_COUNT + 1) / 2
if %_GC_MARK_THREADS% LSS 1 set _GC_MARK_THREADS=1
set JULIA_NUM_GC_THREADS=%_GC_MARK_THREADS%,1

rem Map aliases to canonical executable names
if /I "%SUBCOMMAND%"=="search" set SUBCOMMAND=SearchDIA
if /I "%SUBCOMMAND%"=="predict" set SUBCOMMAND=BuildSpecLib
rem if /I "%SUBCOMMAND%"=="empirical" set SUBCOMMAND=ParseSpecLib
rem if /I "%SUBCOMMAND%"=="params-empirical" set SUBCOMMAND=GetParseSpecLibParams
if /I "%SUBCOMMAND%"=="params-search" set SUBCOMMAND=GetSearchParams
if /I "%SUBCOMMAND%"=="params-predict" set SUBCOMMAND=GetBuildLibParams
if /I "%SUBCOMMAND%"=="convert-mzml" set SUBCOMMAND=convertMzML
if /I "%SUBCOMMAND%"=="convert-raw" set SUBCOMMAND=PioneerConverter


:run_pioneer
set "EXEC=%SCRIPT_DIR%bin\%SUBCOMMAND%.exe"
if /I "%SUBCOMMAND%"=="GetSearchParams" if exist "%SCRIPT_DIR%apps\params\bin\%SUBCOMMAND%.exe" set "EXEC=%SCRIPT_DIR%apps\params\bin\%SUBCOMMAND%.exe"
if /I "%SUBCOMMAND%"=="GetBuildLibParams" if exist "%SCRIPT_DIR%apps\params\bin\%SUBCOMMAND%.exe" set "EXEC=%SCRIPT_DIR%apps\params\bin\%SUBCOMMAND%.exe"
if /I "%SUBCOMMAND%"=="GetParseSpecLibParams" if exist "%SCRIPT_DIR%apps\params\bin\%SUBCOMMAND%.exe" set "EXEC=%SCRIPT_DIR%apps\params\bin\%SUBCOMMAND%.exe"
if /I "%SUBCOMMAND%"=="BuildSpecLib" if exist "%SCRIPT_DIR%apps\predict\bin\%SUBCOMMAND%.exe" set "EXEC=%SCRIPT_DIR%apps\predict\bin\%SUBCOMMAND%.exe"
if /I "%SUBCOMMAND%"=="SearchDIA" if exist "%SCRIPT_DIR%apps\search\bin\%SUBCOMMAND%.exe" set "EXEC=%SCRIPT_DIR%apps\search\bin\%SUBCOMMAND%.exe"
if /I "%SUBCOMMAND%"=="convertMzML" if exist "%SCRIPT_DIR%apps\convert\bin\%SUBCOMMAND%.exe" set "EXEC=%SCRIPT_DIR%apps\convert\bin\%SUBCOMMAND%.exe"

set "SUBAPP_ROOT="
if /I "%SUBCOMMAND%"=="GetSearchParams" set "SUBAPP_ROOT=%SCRIPT_DIR%apps\params"
if /I "%SUBCOMMAND%"=="GetBuildLibParams" set "SUBAPP_ROOT=%SCRIPT_DIR%apps\params"
if /I "%SUBCOMMAND%"=="GetParseSpecLibParams" set "SUBAPP_ROOT=%SCRIPT_DIR%apps\params"
if /I "%SUBCOMMAND%"=="BuildSpecLib" set "SUBAPP_ROOT=%SCRIPT_DIR%apps\predict"
if /I "%SUBCOMMAND%"=="SearchDIA" set "SUBAPP_ROOT=%SCRIPT_DIR%apps\search"
if /I "%SUBCOMMAND%"=="convertMzML" set "SUBAPP_ROOT=%SCRIPT_DIR%apps\convert"

if defined SUBAPP_ROOT (
    if exist "%SUBAPP_ROOT%\share\julia" (
        set "BUNDLE_SHARE=%SUBAPP_ROOT%\share\julia"
        set "BUNDLE_DEV=%BUNDLE_SHARE%\dev"
        set "RUNTIME_DEPOT=%LOCALAPPDATA%\Pioneer\julia\%SUBCOMMAND%"
        if not defined LOCALAPPDATA set "RUNTIME_DEPOT=%USERPROFILE%\AppData\Local\Pioneer\julia\%SUBCOMMAND%"
        if not defined USERPROFILE set "RUNTIME_DEPOT=%TEMP%\Pioneer\julia\%SUBCOMMAND%"
        if not exist "%RUNTIME_DEPOT%" mkdir "%RUNTIME_DEPOT%" >nul 2>&1
        if defined JULIA_LOAD_PATH (
            set "JULIA_LOAD_PATH=%BUNDLE_SHARE%;%BUNDLE_DEV%;@stdlib;%JULIA_LOAD_PATH%"
        ) else (
            set "JULIA_LOAD_PATH=%BUNDLE_SHARE%;%BUNDLE_DEV%;@stdlib"
        )
        if defined JULIA_DEPOT_PATH (
            set "JULIA_DEPOT_PATH=%RUNTIME_DEPOT%;%BUNDLE_SHARE%;%JULIA_DEPOT_PATH%"
        ) else (
            set "JULIA_DEPOT_PATH=%RUNTIME_DEPOT%;%BUNDLE_SHARE%"
        )
    )
)

if /I "%SUBCOMMAND%"=="SearchDIA" call :precompile_search_plot_runtime
if %ERRORLEVEL% NEQ 0 exit /b %ERRORLEVEL%

set "JULIA_ARGS_EXTRA="
if /I "%SUBCOMMAND%"=="BuildSpecLib" set "JULIA_ARGS_EXTRA=--pkgimages=no"
if /I "%SUBCOMMAND%"=="SearchDIA" set "JULIA_ARGS_EXTRA=--pkgimages=no"
if /I "%SUBCOMMAND%"=="BuildSpecLib" set JULIA_PKG_PRECOMPILE_AUTO=0
if /I "%SUBCOMMAND%"=="SearchDIA" set JULIA_PKG_PRECOMPILE_AUTO=0

if "%SUBCOMMAND_ARGS%"=="" (
    if defined JULIA_ARGS_EXTRA (
        "%EXEC%" --julia-args %JULIA_ARGS_EXTRA%
    ) else (
        "%EXEC%"
    )
) else (
    if defined JULIA_ARGS_EXTRA (
        "%EXEC%" %SUBCOMMAND_ARGS% --julia-args %JULIA_ARGS_EXTRA%
    ) else (
        "%EXEC%" %SUBCOMMAND_ARGS%
    )
)
exit /b %ERRORLEVEL%

:precompile_search_plot_runtime
if not defined BUNDLE_SHARE exit /b 0
set "PLOTS_MARKER=%RUNTIME_DEPOT%\pioneer-plots-runtime-smoke-%PIONEER_VERSION%-v1"
if exist "%PLOTS_MARKER%" exit /b 0
set "JULIA_BIN=%SUBAPP_ROOT%\bin\julia.exe"
if not exist "%JULIA_BIN%" (
    echo Error: bundled Julia runtime not found at %JULIA_BIN%
    exit /b 1
)
set "GKSwstype=100"
set "GKS_WSTYPE=100"
set "JULIA_PKG_PRECOMPILE_AUTO=0"
"%JULIA_BIN%" --startup-file=no --pkgimages=no -C generic --project="%BUNDLE_SHARE%" -e "pkgid=Base.PkgId(Base.UUID(0x9adec8dfeb604bac9b0a4137287d3bbd), String(:PioneerPlots)); Base.require(pkgid); m=Base.root_module(pkgid); p=m.Plots.plot([1.0,2.0], [1.0,2.0]; label=nothing, title=String(:PioneerPlots_runtime_smoke)); m.save_multipage_pdf([p], joinpath(mktempdir(), string(:pioneer_plots_runtime_smoke, :., :pdf)))"
if %ERRORLEVEL% NEQ 0 exit /b %ERRORLEVEL%
type nul > "%PLOTS_MARKER%"
exit /b 0
