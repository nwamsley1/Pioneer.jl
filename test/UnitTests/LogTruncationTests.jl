using Test
using Pioneer
import Pioneer: truncate_for_log, MAX_LOG_MSG_BYTES
import Pioneer: DEBUG_CONSOLE_LEVEL, ESSENTIAL_FILE, CONSOLE_FILE, DEBUG_FILE, WARNINGS_FILE

function _with_temp_logs_trunc(f::Function)
    tmp = mktempdir()
    ess_path = joinpath(tmp, "essential.log")
    con_path = joinpath(tmp, "console.log")
    dbg_path = joinpath(tmp, "debug.log")
    wrn_path = joinpath(tmp, "warnings.log")

    ess = open(ess_path, "w")
    con = open(con_path, "w")
    dbg = open(dbg_path, "w")
    wrn = open(wrn_path, "w")

    # Save originals
    old_level = DEBUG_CONSOLE_LEVEL[]
    old_ess = ESSENTIAL_FILE[]
    old_con = CONSOLE_FILE[]
    old_dbg = DEBUG_FILE[]
    old_wrn = WARNINGS_FILE[]
    try
        ESSENTIAL_FILE[] = ess
        CONSOLE_FILE[] = con
        DEBUG_FILE[] = dbg
        WARNINGS_FILE[] = wrn
        f((ess_path, con_path, dbg_path, wrn_path))
    finally
        DEBUG_CONSOLE_LEVEL[] = old_level
        ESSENTIAL_FILE[] = old_ess
        CONSOLE_FILE[] = old_con
        DEBUG_FILE[] = old_dbg
        WARNINGS_FILE[] = old_wrn
        close(ess); close(con); close(dbg); close(wrn)
    end
end

@testset "Log message truncation" begin
    # Keep original setting
    original_max = MAX_LOG_MSG_BYTES[]
    try
        # Short ASCII stays intact
        MAX_LOG_MSG_BYTES[] = 32
        short = "hello world"
        @test truncate_for_log(short) == short

        # Long ASCII truncates and appends suffix
        long_ascii = repeat("x", 1000)
        t = truncate_for_log(long_ascii)
        @test occursin("[truncated", t)
        # Prefix should be exactly max bytes of 'x'
        @test startswith(t, repeat("x", MAX_LOG_MSG_BYTES[]))

        # UTF-8 safe truncation (😀 is 4 bytes)
        MAX_LOG_MSG_BYTES[] = 10
        utf = repeat("😀", 10)
        tu = truncate_for_log(utf)
        @test occursin("[truncated", tu)
        # Should contain two whole emojis (8 bytes) but not split
        @test startswith(tu, repeat("😀", 2))
        @test isvalid(tu)

        # Integration: user_error writes truncated text to files
        _with_temp_logs_trunc() do (ess_path, con_path, dbg_path, wrn_path)
            MAX_LOG_MSG_BYTES[] = 16
            DEBUG_CONSOLE_LEVEL[] = 1
            msg = repeat("A", 200)
            @user_error msg
            ess = read(ess_path, String)
            con = read(con_path, String)
            dbg = read(dbg_path, String)
            wrn = read(wrn_path, String)  # errors not written to warnings
            @test occursin("[truncated", ess)
            @test occursin("[ERROR]", con)
            @test occursin("[error]", dbg)
            @test !isempty(wrn) == false
        end
    finally
        MAX_LOG_MSG_BYTES[] = original_max
    end
end

