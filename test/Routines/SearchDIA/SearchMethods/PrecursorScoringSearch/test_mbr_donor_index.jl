using Test, Random, Arrow, DataFrames
using Pioneer

@testset "Indexed integrated MBR donors preserve selection" begin
    function donor(pid, run, score, weight)
        Pioneer._MBRDonorEntry(score, UInt32(pid), weight, -1f0, 0f0, 10f0,
            5f0, (1f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0), 0x03, UInt16(2), UInt32(run))
    end
    rng = MersenneTwister(812)
    entries = Dict{UInt32, Vector{Pioneer._MBRDonorEntry}}()
    for pid in 1:32
        runs = shuffle(rng, collect(1:8))[1:mod(pid, 9)]
        entries[UInt32(pid)] = [donor(pid, run, rand(rng, Float32[0.5, 0.9, 0.9, 1]),
            rand(rng, Float32[-Inf, 0, 1, 1, Inf, NaN])) for run in runs]
    end
    index = Pioneer._MBRDonorIndex(entries)
    @test length(index) == length(entries)
    @test sum(length, values(index)) == sum(length, values(entries))
    @test all(index[pid] === donors for (pid, donors) in entries)
    shared = Dict((UInt32(i), UInt32(j)) => Float32(mod(i+j, 3) / 3)
        for i in 1:8 for j in i+1:8)
    atlas = Pioneer.RunSimilarityAtlas(shared, Dict{Pioneer.RunPair, Float32}(),
        Dict(UInt32(i) => 1f0 for i in 1:8), Dict{UInt32, Float32}(), BitSet(1:8),
        Dict{UInt32, BitSet}(), Dict{UInt32, Float32}())
    for receiver in UInt32.(0:9), current_atlas in (nothing, atlas)
        context = Pioneer._mbr_receiver_donors(index, receiver, current_atlas)
        for pid in UInt32.(1:33)
            expected = Pioneer._mbr_select_donor(entries, pid, receiver, current_atlas)
            @test isequal(Pioneer._mbr_select_donor(context, pid, receiver, current_atlas), expected)
            @test isequal(Pioneer._mbr_top_scoring_donor(context, pid, receiver),
                Pioneer._mbr_top_scoring_donor(entries, pid, receiver))
            @test isequal(Pioneer._mbr_donor_in_file(context, pid, receiver),
                Pioneer._mbr_donor_in_file(entries, pid, receiver))
            for best in get(entries, pid, Pioneer._MBRDonorEntry[])
                @test isequal(Pioneer._mbr_worst_alternate_donor(context, pid, receiver, best),
                    Pioneer._mbr_worst_alternate_donor(entries, pid, receiver, best))
            end
        end
    end
    tied = Dict(UInt32(1) => [donor(1, 7, 0.9f0, 1f0), donor(1, 2, 0.9f0, 1f0),
        donor(1, 5, 0.9f0, 1f0)])
    tied_index = Pioneer._MBRDonorIndex(tied)
    context = Pioneer._mbr_receiver_donors(tied_index, UInt32(8), nothing)
    @test Pioneer._mbr_select_donor(context, UInt32(1), UInt32(8), nothing).ms_file_idx == 7
    @test Pioneer._mbr_worst_alternate_donor(context, UInt32(1), UInt32(8), tied[UInt32(1)][1]).ms_file_idx == 2
    sorted_index = Pioneer._MBRDonorIndex(Dict(UInt32(1) => sort(tied[UInt32(1)]; by=d -> d.ms_file_idx)))
    @test isempty(sorted_index.lookups[UInt32(1)].file_order)
    @test !isempty(tied_index.lookups[UInt32(1)].file_order)
end

@testset "Integrated donor construction deduplicates without changing order" begin
    mktempdir() do directory
        paths = String[]
        for (part, rows) in enumerate((
            [(1, 5, 0.8f0, 10f0), (1, 5, 0.9f0, 20f0), (1, 2, 0.9f0, 30f0), (2, 2, 0.4f0, 40f0)],
            [(1, 5, 0.9f0, 50f0), (1, 2, 0.95f0, 60f0), (2, 3, 0.9f0, 70f0)],
        ))
            n = length(rows)
            table = DataFrame(precursor_idx=UInt32[row[1] for row in rows],
                ms_file_idx=UInt32[row[2] for row in rows],
                trace_prob_prepass=Float32[row[3] for row in rows],
                qval=fill(0.001f0, n), global_qval=fill(0.001f0, n), irt_pred=fill(11f0, n))
            table[!, Pioneer.MBR_INTEGRATED_WEIGHT_COLUMN] = Float32[row[4] for row in rows]
            table[!, Pioneer.MBR_INTEGRATED_LOG2_INTENSITY_EXPLAINED_COLUMN] = fill(-1f0, n)
            table[!, Pioneer.MBR_INTEGRATED_APEX_IRT_COLUMN] = fill(10f0, n)
            table[!, Pioneer.MBR_INTEGRATED_FRAG_CORR_BITVEC_COLUMN] = fill(0x03, n)
            table[!, Pioneer.MBR_INTEGRATED_N_CORRELATED_FRAGMENTS_BITVEC_RANK_COLUMN] = fill(UInt16(2), n)
            table[!, Pioneer.MBR_INTEGRATED_N_SCANS_COLUMN] = fill(5f0, n)
            for (rank, name) in enumerate(Pioneer.MBR_INTEGRATED_FRAGMENT_SQRT_COLUMNS)
                table[!, name] = fill(rank == 1 ? 1f0 : 0f0, n)
            end
            path = joinpath(directory, "part_$part.arrow")
            open(Arrow.Writer, path; file=true) do writer
                Arrow.write(writer, table[1:1, :])
                Arrow.write(writer, table[2:end, :])
            end
            push!(paths, path)
        end
        index = Pioneer.build_mbr_integrated_donor_dict(paths, 0.5f0; q_value_threshold=0.01f0)
        @test getproperty.(index[UInt32(1)], :ms_file_idx) == UInt32[5, 2]
        @test getproperty.(index[UInt32(1)], :weight) == Float32[20, 60]
        @test getproperty.(index[UInt32(1)], :trace_prob) == Float32[0.9, 0.95]
        @test only(index[UInt32(2)]).ms_file_idx == 3
        @test isempty(Pioneer.build_mbr_integrated_donor_dict(String[], 0.5f0; q_value_threshold=0.01f0).entries)
    end
end
