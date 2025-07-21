@echo off
setlocal enabledelayedexpansion

if "%JULIA_NUM_THREADS%"=="" set JULIA_NUM_THREADS=auto

set SUBCOMMAND=
set SUBCOMMAND_ARGS=
set VALID_COMMANDS=SearchDIA BuildSpecLib GetSearchParams GetBuildLibParams convertMzML ParseSpecLib

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
if "%~1"=="--help" goto show_help
if "%~1"=="-h" goto show_help
if "%~1"=="-help" goto show_help
if "%~1"=="/?" goto show_help

rem Check for --threads=value format
echo %~1 | findstr /C:"--threads=" >nul
if !errorlevel! equ 0 (
    set THREADS_ARG=%~1
    set JULIA_NUM_THREADS=!THREADS_ARG:--threads=!
    shift
    goto parse_args
)

rem Check for unknown options
echo %~1 | findstr /R "^-" >nul
if !errorlevel! equ 0 (
    echo Error: Unknown option %~1
    echo Use --help for usage information
    exit /b 1
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
echo   SearchDIA          Perform DIA search analysis
echo   BuildSpecLib       Build spectral library
echo   GetSearchParams    Generate search parameter template
echo   GetBuildLibParams  Generate library build parameter template
echo   convertMzML        Convert mzML files
echo   ParseSpecLib       Parse spectral library
echo.
echo Examples:
echo   pioneer SearchDIA config.json                    # Use auto threading
echo   pioneer --threads=8 SearchDIA config.json        # Use 8 threads
echo   pioneer --threads 4 BuildSpecLib config.json     # Use 4 threads
echo.
echo For subcommand-specific help:
echo   pioneer ^<subcommand^> --help
exit /b 0

:check_subcommand
if "%SUBCOMMAND%"=="" (
    echo Error: Subcommand required
    echo Valid subcommands: %VALID_COMMANDS%
    echo Use --help for usage information"
    exit /b 1
)

:run_pioneer
rem The executables are in the bin\ subdirectory
if "%SUBCOMMAND_ARGS%"=="" (
    "%~dp0bin\%SUBCOMMAND%.exe"
) else (
    "%~dp0bin\%SUBCOMMAND%.exe" %SUBCOMMAND_ARGS%
)