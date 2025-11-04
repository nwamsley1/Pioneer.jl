@echo off
setlocal enabledelayedexpansion

if "%JULIA_NUM_THREADS%"=="" set JULIA_NUM_THREADS=auto

set SCRIPT_DIR=%~dp0

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
echo.
echo Usage: pioneer [options] ^<subcommand^> [subcommand-args...]
echo.
echo Options:
echo   --threads N        Set number of Julia threads (default: auto)
echo   --threads=N        Alternative syntax for setting threads
echo   --help, -h         Show this help message
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
echo   convert-raw ^<data_path^> [options]
echo                                                 Convert Thermo RAW files
echo   convert-mzml ^<data_path^> [skip_header]
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


:check_subcommand
if "%SUBCOMMAND%"=="" (
    echo Error: Subcommand required
    echo Valid subcommands: %VALID_COMMANDS%
    echo Use --help for usage information
    exit /b 1
)

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
rem The executables are in the bin\ subdirectory
set "EXEC=%SCRIPT_DIR%bin\%SUBCOMMAND%.exe"
if "%SUBCOMMAND_ARGS%"=="" (
    "%EXEC%"
) else (
    "%EXEC%" %SUBCOMMAND_ARGS%
)
