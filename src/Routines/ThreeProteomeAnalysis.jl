function ThreeProteomeAnalysis(
    results_folder::String, #/Users/n.t.wamsley/Desktop/astral_test_out/RESULTS/RESULTS/
    key_file_path::String,  #"/Users/n.t.wamsley/Desktop/astral_test_key_080324.txt"
    out_dir::String #"/Users/n.t.wamsley/Desktop/astral_test_out"
)
    #=
    ThreeProteomeAnalysis(
    "/Users/n.t.wamsley/Desktop/ALTIMETER_PIONEER_111824/AstralCombineTraces/RESULTS/RESULTS",
    "/Users/n.t.wamsley/Desktop/astral_test_key_080324.txt",
    "/Users/n.t.wamsley/Desktop/CombinedTraces111824"
    )


    ThreeProteomeAnalysis(
    "/Users/n.t.wamsley/Desktop/ALTIMETER_PIONEER_111824/AstralSeperateTraces/RESULTS/RESULTS",
    "/Users/n.t.wamsley/Desktop/astral_test_key_080324.txt",
    "/Users/n.t.wamsley/Desktop/SeperateTraces111824"
    )



    ThreeProteomeAnalysis(
    "/Users/n.t.wamsley/Desktop/ALTIMETER_PIONEER_111824/AstralSeperateTracesNoMax/RESULTS/RESULTS",
    "/Users/n.t.wamsley/Desktop/astral_test_key_080324.txt",
    "/Users/n.t.wamsley/Desktop/SeperateTracesNoMax111824"
    )

    =#
    #=
    Example key_file

    ```
    #condition keys
    E10H50Y40,E30H50Y20,E40H50Y10,E45H50Y5
    #condition pairs
    E10H50Y40,E30H50Y20,ECOLI;3.0|HUMAN;1.0|YEAST;0.5
    E10H50Y40,E40H50Y10,ECOLI;4.0|HUMAN;1.0|YEAST;0.25
    E10H50Y40,E45H50Y5,ECOLI;4.5|HUMAN;1.0|YEAST;0.1111111
    ```
    
    =#
    function summarizeCondition(
        condition_columns::Matrix{Union{Missing,  Float32}}
        )
        function summarizeRows(x::AbstractVector{Union{Missing, Float32}})
            x = collect(skipmissing(x))
            if length(x) > 2
                return (CV = Float32(100*(1 + 1/(4*length(x)))*std(x)/mean(x)), Mean = mean(x), N = length(x))
            else
                return (CV = missing, Mean = missing, N = length(x))
            end
        end
        return map(x->summarizeRows(x), eachcol(condition_columns))
    end
    function getModsFilter(
        structural_mods::AbstractVector{Union{Missing,String}})
        has_Mox = zeros(Bool, length(structural_mods))
        for (i, mods) in enumerate(structural_mods)
            has_Mox[i] = !occursin("M,Unimod:35", mods)
        end
        has_Mox
    end
    function getSpeciesFilter(
        species::AbstractVector{String})
        multiple_species = zeros(Bool, length(species))
        for (i, spec) in enumerate(species)
            multiple_species[i] = length(unique(split(spec,';')))==1
        end
        return multiple_species
    end
    function getSpeciesFilter(
        species::AbstractVector{Union{Missing, String}})
        multiple_species = zeros(Bool, length(species))
        for (i, spec) in enumerate(species)
            multiple_species[i] = length(unique(split(spec,';')))==1
        end
        return multiple_species
    end
    function compareConditions(
        condition_a::String,
        condition_b::String,
        base_cols::Vector{Symbol},
        precursors_wide::DataFrame,
        species_to_log2fc::Dict{String, Float64},
        out_dir::String,
        log_scale::Bool,
        table_type::String;
        max_cv::Float32 = 20.0f0,
        min_non_missing::Int64 = 3)
        text_results_path = joinpath(out_dir, condition_a*"_vs_"*condition_b*"_"*string(table_type)*"_results.txt")
        condition_a_cols = Symbol.([condition_a*"_CV", condition_a*"_Mean", condition_a*"_N"])
        condition_b_cols = Symbol.([condition_b*"_CV", condition_b*"_Mean", condition_b*"_N"])
        sub_precursors_wide = precursors_wide[!,vcat([base_cols, condition_a_cols, condition_b_cols]...)]
        row_filter = (sub_precursors_wide[!,condition_a_cols[1]].<=max_cv
        ).&(sub_precursors_wide[!,condition_b_cols[1]].<=max_cv
        ).&(sub_precursors_wide[!,condition_a_cols[end]].>=min_non_missing
        ).&(sub_precursors_wide[!,condition_b_cols[end]].>=min_non_missing)
        sub_precursors_wide[!,:unique_species] = join.(unique.(split.(sub_precursors_wide[!,:species],';')),';')
        sub_precursors_wide = sub_precursors_wide[row_filter,:]
        total_passing_precs = size(sub_precursors_wide,1)
        println("total passing precs $total_passing_precs")
        if !log_scale
            sub_precursors_wide[!,:log2_mean] = (log2.(sub_precursors_wide[!,condition_a_cols[2]]) .+ log2.(sub_precursors_wide[!,condition_b_cols[2]]))./2
            sub_precursors_wide[!,:log2_diff] = (log2.(sub_precursors_wide[!,condition_b_cols[2]]) .- log2.(sub_precursors_wide[!,condition_a_cols[2]]))
        else
            sub_precursors_wide[!,:log2_mean] = (sub_precursors_wide[!,condition_a_cols[2]] .+ sub_precursors_wide[!,condition_b_cols[2]])./2
            sub_precursors_wide[!,:log2_diff] = (sub_precursors_wide[!,condition_b_cols[2]] .- sub_precursors_wide[!,condition_a_cols[2]])
        end
        gdf = groupby(sub_precursors_wide,:unique_species)
        for (species, precursors) in pairs(gdf)
            println("$species has ", size(precursors, 1), " precursors")
            println("$species has ", std(precursors[!,:log2_diff]))
        end

        open(text_results_path,"w") do io1
            println(io1, "total passing precs $total_passing_precs")
            for (species, precursors) in pairs(gdf)
                println(io1, "$species has ", size(precursors, 1), " precursors")
            end
        end
        nt = NamedTuple.(keys(gdf))
        p = plot(legenZd=:outertopright, show = true, title = "Olsen Astral Benchmark Pioneer \n 3of3, CV<20%, 1% FDR \n $total_passing_precs Precursors",
        topmargin=5mm, dpi = 300)
        i = 1
        #xlim = (quantile(sub_precursors_wide[!,:log2_diff], 0.001), quantile(sub_precursors_wide[!,:log2_diff], 0.999))
        xlim = (-3.5, 3.5)
        ylim = (0, 2.5)
        #SPECIES_TO_LOG2FC = Dict("HUMAN" => 0.0,
        #                     "YEAST" => -1.0,
        #                     "ECOLI" => log2(3.0))
        for k in ["HUMAN","YEAST","ECOLI"]
            density!(p, gdf[(unique_species = k,)][:,:log2_diff], color = i, label=nothing, bins = LinRange(-3, 3, 100), xlim = xlim, ylim = ylim, show = true, normalize = :probability)
            vline!([species_to_log2fc[k]], color = i, legend = true,
                label="$(nt[i][:unique_species])"
            )
        
            i += 1
        end

        savefig(joinpath(out_dir, condition_a*"_vs_"*condition_b*"_"*string(table_type)*"_densityplot.pdf"))
        Arrow.write(
            joinpath(proteome_results_dir, condition_a*"_vs_"*condition_b*"_"*string(table_type)*"_table.arrow"),
            sub_precursors_wide)
        return 
    end

    pg_wide = DataFrame(Tables.columntable(Arrow.Table(joinpath(results_folder, "protein_groups_wide.arrow"))))
    proteome_results_dir = joinpath(out_dir,"three_proteome_results")
    if !isdir(proteome_results_dir)
        mkdir(proteome_results_dir)
    end

    content_ = read(key_file_path, String)
    lines = split(content_, '\n')
    condition_keys = Vector{String}(split(lines[2],','))
    comparisons = lines[4:end]
    condition_to_columns = Dictionary{String, Vector{Symbol}}()
    for condition_key in condition_keys
        insert!(
            condition_to_columns,
            condition_key,
            [Symbol(colname) for colname in names(pg_wide) if occursin(condition_key, colname)]
        )
    end
    #pg_wide = pg_wide[getModsFilter(pg_wide[!,:structural_mods]),:]
    pg_wide = pg_wide[getSpeciesFilter(pg_wide[!,:species]),:]
    for (condition_key, condition_columns) in pairs(condition_to_columns)
        row_major = collect(Matrix(pg_wide[!,condition_columns])')
        condition_newcols = DataFrame(summarizeCondition(row_major))
        pg_wide[!,Symbol(condition_key*"_CV")] = condition_newcols[!,:CV]
        pg_wide[!,Symbol(condition_key*"_Mean")] = condition_newcols[!,:Mean]
        pg_wide[!,Symbol(condition_key*"_N")] = condition_newcols[!,:N]
    end
    #base_cols = [Symbol(x) for x in ["species","accession_numbers","sequence","structural_mods","isotopic_mods","precursor_idx","target"]]
    base_cols = [Symbol(x) for x in ["species","protein","target"]]
    for comparison in comparisons
        if length(comparison)==0
            continue
        end
        split_comparison = split(comparison,',')
        condition_a = String(split_comparison[1])
        condition_b = String(split_comparison[2])
        species_to_log2fc = Dict{String, Float64}()
        for _species_ in split(split_comparison[3],'|')
            spec, fc = split(_species_,';')
            species_to_log2fc[spec] = log2(parse(Float64, fc))
        end
        println("species_to_log2fc $species_to_log2fc")
        compareConditions(
            condition_a,
            condition_b,
            base_cols,
            pg_wide,
            species_to_log2fc,
            proteome_results_dir,
            false,
            "protein_groups",
            max_cv = 20.0f0#typemax(Float32)
        )
    end

    precursors_wide = DataFrame(Tables.columntable(Arrow.Table(joinpath(results_folder, "precursors_wide.arrow"))))
    precursors_results_dir = joinpath(out_dir,"three_proteome_results")
    if !isdir(proteome_results_dir)
        mkdir(proteome_results_dir)
    end

    content = read(key_file_path, String)
    lines = split(content, '\n')
    condition_keys = Vector{String}(split(lines[2],','))
    comparisons = lines[4:end]
    condition_to_columns = Dictionary{String, Vector{Symbol}}()
    for condition_key in condition_keys
        insert!(
            condition_to_columns,
            condition_key,
            [Symbol(colname) for colname in names(precursors_wide) if occursin(condition_key, colname)]
        )
    end
    #precursors_wide = precursors_wide[getModsFilter(precursors_wide[!,:structural_mods]),:]
    precursors_wide = precursors_wide[getSpeciesFilter(precursors_wide[!,:species]),:]
    for (condition_key, condition_columns) in pairs(condition_to_columns)
        row_major = collect(Matrix(precursors_wide[!,condition_columns])')
        condition_newcols = DataFrame(summarizeCondition(row_major))
        precursors_wide[!,Symbol(condition_key*"_CV")] = condition_newcols[!,:CV]
        precursors_wide[!,Symbol(condition_key*"_Mean")] = condition_newcols[!,:Mean]
        precursors_wide[!,Symbol(condition_key*"_N")] = condition_newcols[!,:N]
    end
    base_cols = [Symbol(x) for x in ["species","accession_numbers","sequence","structural_mods","isotopic_mods","precursor_idx","target"]]
    #base_cols = [Symbol(x) for x in ["species","protein","target"]]
    for comparison in comparisons
        if length(comparison)==0
            continue
        end
        split_comparison = split(comparison,',')
        condition_a = String(split_comparison[1])
        condition_b = String(split_comparison[2])
        species_to_log2fc = Dict{String, Float64}()
        for _species_ in split(split_comparison[3],'|')
            spec, fc = split(_species_,';')
            species_to_log2fc[spec] = log2(parse(Float64, fc))
        end
        println("species_to_log2fc $species_to_log2fc")
        compareConditions(
            condition_a,
            condition_b,
            base_cols,
            precursors_wide,
            species_to_log2fc,
            precursors_results_dir,
            false,
            "precursors",
            max_cv = 20.0f0,#typemax(Float32)
        )
    end

end
function ThreeProteomeAnalysis(
    params_path::String,
    key_file_path::String
)
    params = JSON.parse(read(params_path, String));
    out_dir = params["benchmark_params"]["results_folder"]
    results_folder = joinpath(out_dir, "RESULTS","RESULTS")

    ThreeProteomeAnalysis(results_folder,
                        key_file_path,
                        out_dir
                        )
end
#condition_a = [Symbol("E10H50Y40_1"), Symbol("E10H50Y40_2"), Symbol("E10H50Y40_3"), Symbol("E10H50Y40_4")]
#condition_b = [Symbol("E30H50Y20_1"), Symbol("E30H50Y20_2"), Symbol("E30H50Y20_3"), Symbol("E30H50Y20_4")]
#N = size(precursors_wide,1)
#precursors_wide[!,:mean_a], precursors_wide[!,:mean_b] = zeros(Float32, N), zeros(Float32, N)