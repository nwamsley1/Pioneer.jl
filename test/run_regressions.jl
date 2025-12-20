using Pioneer

param_dir = get(
    ENV,
    "PIONEER_PARAMS_DIR",
    joinpath(@__DIR__, "..", "pioneer-regression-configs", "params"),
)

for path in readdir(param_dir; join=true)
    endswith(path, ".json") || continue
    params_name, _ = splitext(basename(path))
    data_path = joinpath(raw"C:\pioneer-runner\data", params_name)
    if !isdir(data_path)
        @info "Skipping SearchDIA due to missing data path" path data_path
        continue
    end
    @info "Running SearchDIA" path
    SearchDIA(path)
end
