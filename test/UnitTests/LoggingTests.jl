using Test
using Dates

# Import macros and internal logging state
using Pioneer: @user_info, @user_warn, @user_error, @user_print,
               @debug_l1, @debug_l2, @debug_l3, @trace
import Pioneer: DEBUG_CONSOLE_LEVEL, ESSENTIAL_FILE, CONSOLE_FILE, DEBUG_FILE, WARNINGS_FILE

function _with_temp_logs(f::Function)
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
        # Restore refs and level
        DEBUG_CONSOLE_LEVEL[] = old_level
        ESSENTIAL_FILE[] = old_ess
        CONSOLE_FILE[] = old_con
        DEBUG_FILE[] = old_dbg
        WARNINGS_FILE[] = old_wrn
        # Ensure streams are closed
        close(ess); close(con); close(dbg); close(wrn)
    end
end

@testset "Logging macros" begin
    # user_info with essential pattern
    _with_temp_logs() do (ess_path, con_path, dbg_path, wrn_path)
        DEBUG_CONSOLE_LEVEL[] = 0
        @user_info "Starting search at: dummy"
        # Flush by reopening
        ess = read(ess_path, String)
        con = read(con_path, String)
        dbg = read(dbg_path, String)
        @test occursin("Starting search at:", ess)
        @test occursin("[INFO]", con)
        @test occursin("[info]", dbg)
    end

    # warn shows a two-line console entry when level > 0
    _with_temp_logs() do (ess_path, con_path, dbg_path, wrn_path)
        DEBUG_CONSOLE_LEVEL[] = 1
        @user_warn "Something happened"
        con = read(con_path, String)
        wrn = read(wrn_path, String)
        dbg = read(dbg_path, String)
        @test occursin("[WARN] Something happened", con)
        # Should include a second line with source location marker
        @test occursin("└ @", con)
        @test occursin("Something happened", wrn)
        @test occursin("[warn] Something happened", dbg)
    end

    # debug_l1 gated by level
    _with_temp_logs() do (ess_path, con_path, dbg_path, wrn_path)
        DEBUG_CONSOLE_LEVEL[] = 0
        @debug_l1 "no-console"
        dbg = read(dbg_path, String)
        @test !occursin("DEBUG1", dbg)
    end
    _with_temp_logs() do (ess_path, con_path, dbg_path, wrn_path)
        DEBUG_CONSOLE_LEVEL[] = 1
        @debug_l1 "hello"
        dbg = read(dbg_path, String)
        @test occursin("[DEBUG1] hello", dbg)
    end

    # debug_l2 includes source location in debug file when present
    _with_temp_logs() do (ess_path, con_path, dbg_path, wrn_path)
        DEBUG_CONSOLE_LEVEL[] = 2
        @debug_l2 "step"
        dbg = read(dbg_path, String)
        @test occursin("[DEBUG2] step", dbg)
        # The macro passes file/line; debug file entry includes " @ "
        @test occursin(" @ ", dbg)
    end

    # user_print mirrors console and essential/debug
    _with_temp_logs() do (ess_path, con_path, dbg_path, wrn_path)
        DEBUG_CONSOLE_LEVEL[] = 0
        @user_print "plain-line"
        ess = read(ess_path, String)
        con = read(con_path, String)
        dbg = read(dbg_path, String)
        @test occursin("plain-line", ess)
        @test occursin("plain-line", con)
        @test occursin("plain-line", dbg)
    end
end

