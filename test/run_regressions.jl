using Pioneer

param_dir = get(
    ENV,
    "PIONEER_PARAMS_DIR",
    joinpath(@__DIR__, "..", "pioneer-regression-configs", "params"),
)

for path in readdir(param_dir; join=true)
    endswith(path, ".json") || continue
    @info "Running SearchDIA" path
    SearchDIA(path)
end