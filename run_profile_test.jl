#!/usr/bin/env julia

using Pkg
Pkg.activate(".")

using Pioneer
using Profile, PProf

function run_profiled_secondpass()
    println("🎯 Running SearchDIA with SecondPassSearch profiling...")
    
    params_path = "./data/ecoli_test/ecoli_test_params.json"
    
    if !isfile(params_path)
        error("E. coli test parameters not found: $params_path")
    end
    
    println("📋 Starting SearchDIA with profiling enabled...")
    
    # Clear any existing profile data
    Profile.clear()
    
    # Enable profiling around the SearchDIA call
    @profile SearchDIA(params_path)
    
    println("✅ SearchDIA completed. Analyzing profile data...")
    
    # Save profile data
    PProf.pprof(; out="profile_secondpass.pb.gz", web=false)
    
    println("📊 Profile saved to: profile_secondpass.pb.gz")
    println("🌐 To view flame graph: PProf.pprof(\"profile_secondpass.pb.gz\"; web=true)")
    
    # Print top functions by time
    println("\n🔥 Top functions by sample count:")
    Profile.print(; maxdepth=15, noisefloor=0.01)
    
    return "profile_secondpass.pb.gz"
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        profile_file = run_profiled_secondpass()
        println("\n✅ Profiling completed!")
        println("📈 Open flame graph with: PProf.pprof(\"$profile_file\"; web=true)")
    catch e
        println("❌ Error during profiling: $e")
        rethrow(e)
    end
end