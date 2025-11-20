using Pioneer

param_dir = "C:/pioneer-runner/params" 

for path in readdir(param_dir; join=true)
    endswith(path, ".json") || continue
    @info "Running SearchDIA" path
    SearchDIA(path)
end