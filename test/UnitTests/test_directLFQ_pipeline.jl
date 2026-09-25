using Test, DataFrames, Arrow

@testset "directLFQ chunked protein output" begin
    mktempdir() do dir
        function observations(protein, nprec)
            DataFrame(inferred_protein_group=fill(protein, 3nprec),
                precursor_idx=repeat(UInt32.(1:nprec); inner=3),
                ms_file_idx=repeat(UInt32[1,2,3]; outer=nprec),
                peak_area=repeat(Float32[10,20,40]; outer=nprec),
                peak_area_normalized=repeat(Float32[100,200,400]; outer=nprec),
                use_for_protein_quant=fill(true, 3nprec),
                target=fill(true, 3nprec), entrapment_group_id=zeros(UInt8, 3nprec),
                global_pg_qval=fill(0.001f0, 3nprec), pg_qval=fill(0.001f0, 3nprec),
                pg_pep=fill(0.002f0, 3nprec), pg_score=fill(4.0f0, 3nprec),
                global_pg_score=fill(5.0f0, 3nprec))
        end
        a = observations("A", 101)
        b = observations("B", 1)[1:2, :]
        rejected = observations("C", 1)
        rejected.pg_qval .= 0.1f0
        refs = Pioneer.PSMFileReference[]
        for (i, table) in enumerate((a, vcat(b, rejected)))
            path = joinpath(dir, "chunk_$i.arrow")
            Arrow.write(path, table)
            push!(refs, Pioneer.PSMFileReference(path))
        end
        seqs = ["PEPTIDE$i" for i in 1:101]
        for (column, factor) in ((:peak_area, 1), (:peak_area_normalized, 10))
            path = joinpath(dir, "$column.arrow")
            Pioneer.LFQ_chunked(refs, path, column, ["r1","r2","r3"], seqs,
                fill(missing, 101), fill(missing, 101), 0.01f0, Dict("A"=>"human", "B"=>"yeast");
                batch_size=2, quantification_method=:directlfq)
            out = DataFrame(Arrow.Table(path))
            @test nrow(out) == 5
            @test out.protein == ["A","A","A","B","B"]
            @test out.file_name == ["r1","r2","r3","r1","r2"]
            @test out.species == ["human","human","human","yeast","yeast"]
            @test out.n_precursors == [101,101,101,1,1]
            @test out.n_precursors_quantified == [100,100,100,1,1]
            @test out.abundance ≈ factor .* [1000,2000,4000,10,20] rtol=2e-6
            @test out.total_peak_area ≈ factor .* [1010,2020,4040,10,20]
            @test out.qval == fill(0.001f0, 5)
        end
    end
end

@testset "Protein quantification parameter selection" begin
    defaults = JSON.parsefile(joinpath(@__DIR__, "..", "..", "assets", "example_config", "defaultSearchParams.json"))
    @test defaults["maxLFQ"]["quantification_method"] == "sparsemaxlfq"
    mktempdir() do dir
        for method in ("sparsemaxlfq", "directlfq", "maxlfq")
            defaults["maxLFQ"]["quantification_method"] = method
            path = joinpath(dir, "params.json")
            write(path, JSON.json(defaults))
            params = Pioneer.parse_pioneer_parameters(path)
            selected = Pioneer.ProteinQuantificationSearchParameters(params)
            @test selected.quantification_method == Symbol(method)
        end
        delete!(defaults["maxLFQ"], "quantification_method")
        path = joinpath(dir, "omitted.json")
        write(path, JSON.json(defaults))
        selected = Pioneer.ProteinQuantificationSearchParameters(Pioneer.parse_pioneer_parameters(path))
        @test selected.quantification_method == :sparsemaxlfq
        defaults["maxLFQ"]["quantification_method"] = "typo"
        path = joinpath(dir, "bad.json")
        write(path, JSON.json(defaults))
        @test_throws Pioneer.InvalidParametersError Pioneer.checkParams(path)
    end
end

@testset "Sparse MaxLFQ production routing" begin
    mktempdir() do dir
        nr, np = 37, 5
        areas = Float32[exp2(5 + 0.03r + sin(p*r)) for p in 1:np, r in 1:nr]
        table = DataFrame(inferred_protein_group=fill("A", nr*np),
            precursor_idx=repeat(UInt32.(1:np); outer=nr),
            ms_file_idx=repeat(UInt32.(1:nr); inner=np), peak_area=vec(areas),
            use_for_protein_quant=fill(true, nr*np), target=fill(true,nr*np),
            entrapment_group_id=zeros(UInt8,nr*np),
            global_pg_qval=fill(0.001f0,nr*np), pg_qval=fill(0.001f0,nr*np),
            pg_pep=fill(0.002f0,nr*np), pg_score=fill(4f0,nr*np),
            global_pg_score=fill(5f0,nr*np))
        singleton = copy(table[1:1,:])
        singleton.inferred_protein_group .= "B"
        append!(table,singleton)
        input = joinpath(dir,"input.arrow")
        Arrow.write(input,table)
        output = joinpath(dir,"output.arrow")
        Pioneer.LFQ_chunked([Pioneer.PSMFileReference(input)], output, :peak_area,
            ["r$i" for i in 1:nr], ["PEPTIDE$i" for i in 1:np],
            fill(missing,np),fill(missing,np),0.01f0,Dict("A"=>"human","B"=>"human");
            quantification_method=:sparsemaxlfq)
        out=DataFrame(Arrow.Table(output))
        X=Union{Missing,Float64}[log2(Float64(v)) for v in areas]
        expected=Pioneer.solve_sparse_maxlfq(X,fill(4f0,nr))
        @test Pioneer.SPARSE_MAXLFQ_PARTNERS == 16
        @test Pioneer.SPARSE_MAXLFQ_SEED == 0
        @test out[out.protein .== "A",:abundance] ≈ exp2.(expected.estimates) rtol=2e-6
        full,_=Pioneer.solve_maxlfq(X,Union{Missing,Float64}[4.0 for _ in 1:nr])
        @test maximum(abs.(full .- expected.estimates)) > 1e-4
        @test only(out[out.protein .== "B",:abundance]) ≈ areas[1,1]
        @test all(==(np),out[out.protein .== "A",:n_precursors_quantified])
    end
end
