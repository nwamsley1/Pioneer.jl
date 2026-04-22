import HostCPUFeatures

root = joinpath(@__DIR__)
data_dir = joinpath(root, "..", "..", "data")

cmd = get(ENV, "PIONEER_CMD", nothing)

function maybe_mask_precompile_cpu_features()
    target = get(ENV, "PIONEER_PRECOMPILE_CPU_TARGET", nothing)
    if !isnothing(target) && ((Sys.ARCH === :x86_64) || (Sys.ARCH === :i686))
        @info "Overriding HostCPUFeatures during precompile for target $target"
        # PackageCompiler records precompile signatures while tracing on the host CPU.
        # On AVX-512 runners that can capture VectorizationBase methods which LLVM
        # cannot legalize for the portable app CPU target later in the build.
        @eval HostCPUFeatures begin
            has_feature(::Val{:x86_64_avx512vl}) = False()
            has_feature(::Val{:x86_64_avx512f}) = False()
            has_feature(::Val{:x86_64_avx2}) = False()
            has_feature(::Val{:x86_64_3}) = False()
            has_feature(::Val{:x86_64_avx}) = False()
        end
    end
end

maybe_mask_precompile_cpu_features()

using Pioneer

function maybe_run(f, name)
    if cmd === nothing || cmd == name
        start_time = time()
        @user_info "Starting precompile target $name"
        try
            f()
            elapsed = round(time() - start_time; digits=1)
            @user_info "Finished precompile target $name in $(elapsed)s"
        catch e
            bt = catch_backtrace()
            target_cmd = cmd === nothing ? "<all>" : cmd
            elapsed = round(time() - start_time; digits=1)
            @user_warn "Error executing $name during precompile of $target_cmd: $(sprint(showerror, e))"
            @warn "Precompile exception details for $name after $(elapsed)s" exception=(e, bt)
        end
    end
end

##########################################
# Generate params
##########################################

# search
maybe_run("GetSearchParams") do
    Pioneer.GetSearchParams(
        joinpath(data_dir, "ecoli_test", "altimeter_ecoli.poin"),
        joinpath(data_dir, "ecoli_test", "raw"),
        mktempdir(),
    )
end

# predict
maybe_run("GetBuildLibParams") do
    Pioneer.GetBuildLibParams(mktempdir(), "test_lib", joinpath(data_dir, "fasta"))
end

# empirical
maybe_run("GetParseSpecLibParams") do
    Pioneer.GetParseSpecLibParams("test_lib", mktempdir())
end


##########################################
# Empirical libraries
##########################################
# ParseSpecLib is currently disabled in importScripts.jl due EmpiricalLibrary loading issues.
# Keep this precompile target disabled until ParseSpecLib is re-enabled in the module.
# maybe_run("ParseSpecLib") do
#     Pioneer.ParseSpecLib(joinpath(data_dir, "precompile", "build_empirical.json"))
# end



##########################################
# Predict libraries
##########################################

# Build a tiny Prosit library and search it with low memory thresholds
#maybe_run("BuildSpecLib") do
#    Pioneer.BuildSpecLib(joinpath(data_dir, "precompile", "build_ecoli_prosit.json"))
#end
# Build a tiny Altimeter library
maybe_run("BuildSpecLib") do
    Pioneer.BuildSpecLib(joinpath(data_dir, "precompile", "build_ecoli_altimeter.json"))
end


##########################################
# Search
##########################################

maybe_run("SearchDIA") do
#    Pioneer.SearchDIA(joinpath(data_dir, "precompile", "search_ecoli_prosit.json"))         # prosit
    Pioneer.SearchDIA(joinpath(data_dir, "precompile", "search_yeast_altimeter.json"))      # altimeter + MBR
    Pioneer.SearchDIA(joinpath(data_dir, "precompile", "search_yeast_altimeter_OOM.json"))  # altimeter + MBR + OOM
end


##########################################
# ConvertMzML
##########################################
maybe_run("convertMzML") do
    Pioneer.convertMzML(joinpath(data_dir, "precompile", "convert_example.mzML"))
end
