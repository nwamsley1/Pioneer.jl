using Pioneer, JSON, Statistics, Serialization
# REPLAY_ROOT contains replay_inputs/metadata.jls and replays/matrices.jsonl.
# These are the saved pre-quantification inputs from the MTAC replay experiment.
2 <= length(ARGS) <= 4 || error("Usage: julia --project scripts/evaluate_sparse_maxlfq.jl REPLAY_ROOT OUTPUT_DIR [PARTNERS_CSV] [SEEDS]")
const PARTNER_COUNTS = length(ARGS) >= 3 ? parse.(Int, split(ARGS[3], ',')) : collect(0:5)
const SEEDS = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 20
all(>=(0), PARTNER_COUNTS) && SEEDS > 0 || error("Invalid comparison budgets or seed count")
const EROOT = abspath(ARGS[1])
const OUT = abspath(ARGS[2])
mkpath(OUT)
const SP = Pioneer
function evaluate_sparse_lfq()
    # Restore the production run priorities for exact component-selection parity.
    meta = deserialize(joinpath(EROOT,"replay_inputs","metadata.jls"))
    prior = Dict{Tuple{String,Int,String},Float64}()
    for path in meta.paths
        table = SP.Arrow.Table(path)
        for i in 1:length(table.precursor_idx)
            table.target[i] || continue
            ismissing(table.inferred_protein_group[i]) && continue
            coalesce(table.pg_qval[i] <= meta.q, false) || continue
            coalesce(table.global_pg_qval[i] <= meta.q, false) || continue
            table.use_for_protein_quant[i] || continue
            key = (String(table.inferred_protein_group[i]), Int(table.entrapment_group_id[i]), meta.names[table.ms_file_idx[i]])
            get!(prior, key, Float64(table.pg_score[i]))
        end
    end
    # Per-run profiles and graph diagnostics allow paired accuracy evaluation.
    nproteins=0
    open(joinpath(OUT,"profiles.jsonl"),"w") do out
        for line in eachline(joinpath(EROOT,"replays","matrices.jsonl"))
            item = JSON.parse(line)
            nr = length(item["runs"])
            np = length(item["matrix"])
            X = Matrix{Union{Missing,Float64}}(missing,np,nr)
            for p in 1:np,r in 1:nr
                v = item["matrix"][p][r]
                v === nothing || (X[p,r]=v)
            end
            priorities = Union{Missing,Float64}[prior[(item["protein"],item["entrap_id"],r)] for r in item["runs"]]
            full, _ = SP.solve_maxlfq(X, priorities)
            records = Any[]
            for partners in PARTNER_COUNTS, seed in 1:SEEDS
                result = SP.solve_sparse_maxlfq(X, priorities; partners, seed)
                estimates = result.estimates
                # Match production's single quantified-run fallback.
                quantified = [r for r in 1:nr if any(!ismissing,view(X,:,r))]
                if length(quantified)==1
                    r=only(quantified)
                    estimates[r]=log2(sum(exp2,skipmissing(view(X,:,r))))
                    full[r]=estimates[r]
                end
                push!(records,(;partners,seed,edges=result.edge_count,iterations=result.iterations,
                    values=[ismissing(v) ? nothing : Float64(v) for v in estimates]))
            end
            item["full"]=[ismissing(v) ? nothing : Float64(v) for v in full]
            item["sparse"]=records
            delete!(item,"matrix")
            println(out,JSON.json(item))
            nproteins+=1
            nproteins%2000==0 && (println("Evaluated $nproteins proteins"); flush(stdout))
        end
    end
    println("Finished $nproteins proteins")
end
evaluate_sparse_lfq()
