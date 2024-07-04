function getCondition(row, file_col::Symbol)
    file_name = split(split(row[file_col], '/')[end], '.')[1]
    file_name_split = split(file_name, '_')
    condition, biological, technical = file_name_split[end - 1], "", file_name_split[end]
    return condition, biological, technical
end
function getLog2Diff(prec_group::SubDataFrame,
                     group_keys::Tuple{String, String})
    #println(size(prec_group))
    if size(prec_group)[1].!=2
        return (log2_diff = zeros(Float64, 1),
                 log2_mean = zeros(Float64, 1),
                 log2_a = zeros(Float64, 1),
                 log2_b = zeros(Float64, 1),
                 CV_a = zeros(Float64, 1),
                 CV_b = zeros(Float64, 1),
                 nonMissing_a = zeros(Int64, 1),
                 nonMissing_b = zeros(Int64, 1))
    end
    group_A = occursin.(first(group_keys), prec_group[:,:condition])
    group_B = occursin.(last(group_keys), prec_group[:,:condition])
    #if any(group_A) == false
    log2_diff = prec_group[:,:log2_mean][group_B] - prec_group[:,:log2_mean][group_A]
    out = (
     log2_diff = log2_diff, 
     log2_mean = (prec_group[:,:log2_mean][group_A] + prec_group[:,:log2_mean][group_B])/2,
     log2_a = prec_group[:,:log2_mean][group_A],
     log2_b = prec_group[:,:log2_mean][group_B],
     CV_a = prec_group[:,:CV][group_A],
     CV_b = prec_group[:,:CV][group_B],
     nonMissing_a = prec_group[:,:non_missing][group_A],
     nonMissing_b = prec_group[:,:non_missing][group_B]
        )
    return out
end
function groupPSMS(
    psms::DataFrame,
    quant_col::Symbol)
    psms = psms[(psms[:,:decoy].==false
                ).&(psms[:,:Mox].==0
                #).&(psms[:,:points_above_FWHM_01].>1
                ),:]
    psms[!,quant_col] = max.(psms[!,quant_col], 0.0)
    gpsms = groupby(psms, 
                            [:species,:accession_numbers,:modified_sequence,:charge,:isotopes_captured,:condition]
                        );

    pioneer_groupmeans = combine(prec_group -> (log2_mean = log2(mean(prec_group[!,quant_col])), 
                                        CV = 100*(1 + (1/(4*size(prec_group, 1))))*std(prec_group[!,quant_col])/mean(prec_group[!,quant_col]),
                                        non_missing = length(prec_group[!,quant_col]),
                                        min_q_value = minimum(prec_group[!,:q_value])), 
    gpsms);
    filter!(x->!isnan(x.CV), pioneer_groupmeans)
    return pioneer_groupmeans#combine(prec_group -> getLog2Diff(prec_group), groupby(pioneer_groupmeans, [:species, :accession_numbers, :modified_sequence, :charge,:isotopes_captured]));
end
function filterResults(
        psms::DataFrame,
        min_non_missing::Int,
        max_cv::AbstractFloat
    )
    cols_to_keep = (psms[!,:nonMissing_a].>=min_non_missing).&(
     psms[!,:nonMissing_b].>=min_non_missing).&(
     psms[!,:CV_a].<=max_cv).&(   
     psms[!,:CV_b].<=max_cv)
    return psms[cols_to_keep,:]
end
function parseMods(mods_string::AbstractString)::Base.RegexMatchIterator
    #Example: "1(7,M,Oxidation)(13,K,AnExampleMod)"
    mods_regex = r"(?<=\().*?(?=\))"
    return eachmatch(mods_regex, mods_string)
end
function getModIndex(mod_string::AbstractString)::UInt8
    parse(UInt8, match(r"^[0-9]+(?=,)", mod_string).match)
end
function getModName(mod_string::AbstractString)::String
    match(r"[^,]+(?=$)", mod_string).match
end
function insert_at_indices(original::String, insertions::Vector{Tuple{String, UInt8}})
    # Convert the original string into an array of characters for easier manipulation
    char_array = collect(original)

    # Sort the insertions by index in ascending order
    sorted_insertions = sort(insertions, by = x -> x[2])

    # Adjust the index for each insertion
    offset = 0
    for (substr, idx) in sorted_insertions
        # Adjust the index with the current offset
        insertion_index = idx + offset
        # Insert each character of the substring at the specified index
        for (i, char) in enumerate(substr)
            insert!(char_array, insertion_index + i, char)
        end
        # Update the offset by the length of the inserted substring
        offset += length(substr)
    end

    # Join the array of characters back into a single string
    return join(char_array)
end
function getModifiedSequence(
    sequence::String,
    isotope_mods::String,
    structural_mods::String)

    mods = structural_mods*isotope_mods
    mods = [("("*getModName(mod.match)*")", getModIndex(mod.match)) for mod in parseMods(mods)]
    return insert_at_indices(sequence, mods)
end

using FASTX, CodecZlib
function parseFasta(fasta_path::String, parse_identifier::Function = x -> split(x,"|")[2])

    function getReader(fasta_path::String)
        if endswith(fasta_path, ".fasta.gz")
            return FASTA.Reader(GzipDecompressorStream(open(fasta_path)))
        elseif endswith(fasta_path, ".fasta")
            return FASTA.Reader(open(fasta_path))
        else
            throw(ErrorException("fasta_path \"$fasta_path\" did not end with `.fasta` or `.fasta.gz`"))
        end
    end

    #I/O for Fasta
    reader = getReader(fasta_path)

    #In memory representation of FastaFile
    #fasta = Vector{FastaEntry}()
    fasta = Vector{Tuple{String, String}}()
    @time begin
        for record in reader
                push!(fasta, 
                        (parse_identifier(FASTA.identifier(record)),
                         split(split(split(FASTA.description(record), ' ')[1], '|')[end], '_')[end],
                                #FASTA.sequence(record),
                                #false
                        )
                )
        end
    end

    return fasta
end
total_time = @timed begin
    println("A")
end
#Create results folder 
benchmark_results_folder = results_folder#params_[:benchmark_params]["results_folder"]
if !isdir(benchmark_results_folder )
    mkpath(benchmark_results_folder )
end

open(joinpath(benchmark_results_folder , "config.json"),"w") do f
    JSON.print(f, params)
end

rpathtext = joinpath(benchmark_results_folder, "results.txt")
rm(rpathtext, force = true)
touch(rpathtext)

open(rpathtext, "w") do file
    write(file, "")
end

open(joinpath(benchmark_results_folder, "results.txt"), "a") do io
    original_stdout = stdout
    redirect_stdout(io) do 
        #println("total time $total_time")
        for (key, value) in pairs(PSMs_Dict)
            hits_at_01 = sum(value[!,:q_value].<=0.01)
            hits_at_10 = sum(value[!,:q_value].<=0.1)
            println("$hits_at_01 precursors at 1% FDR in first search for $key")
            println("$hits_at_10 precursors at 10% FDR in first search for $key")
        end
       println
       println(IDs_PER_FILE)
    end
    redirect_stdout(original_stdout)
    #[println("Precursors w/ CV under $x ", sum(gpsms[!,:CV].<=x)) for x in [5.0, 10.0, 20.0, 30.0]]
end

best_psms_passing[!,:sequence]
best_psms_passing[!,:structural_mods] = [precursors[:structural_mods][pid] for pid in best_psms_passing[!,:precursor_idx]]
best_psms_passing[!,:isotopic_mods] = [precursors[:isotopic_mods][pid] for pid in best_psms_passing[!,:precursor_idx]]
best_psms_passing[!,:modified_sequence] .= ""
for i in range(1, size(best_psms_passing, 1))
    best_psms_passing[i,:modified_sequence] = getModifiedSequence(
        best_psms_passing[i,:sequence],
        best_psms_passing[i,:structural_mods],
        ""
    )
    best_psms_passing[i,:modified_sequence] = "_"*best_psms_passing[i,:modified_sequence]*"_."*string(best_psms_passing[i,:charge])
end
jldsave(joinpath(benchmark_results_folder, "pioneer_all_passing.jld2"); best_psms_passing)

transform!(best_psms_passing, AsTable(:) => ByRow(precursor -> getCondition(precursor, :file_name)) => [:condition, :biological, :technical]);
best_psms_passing[!,:species] .= "HUMAN"
println("size(best_psms_passing) ", size(best_psms_passing))
conditions = unique(best_psms_passing[!,:condition])
conditions = ["E10H50Y40","E30H50Y20"]
quant_name = :trapezoid_area_normalized
cv_param = 20.0
#best_psms_passing_old = copy(best_psms_passing)
ACC_TO_SPEC = Dict(vcat([parseFasta("/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/PROTEOMES/UP000000625_83333_Escherichia_coli.fasta.gz"),
                   parseFasta("/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/PROTEOMES/UP000002311_559292_Saccharomyces_cerevisiae.fasta.gz"),
                   parseFasta("/Users/n.t.wamsley/RIS_temp/ASMS_2024/UNISPEC_LIBS/PROTEOMES/UP000005640_9606_human.fasta.gz")]...));
best_psms_passing[!,:species_names] = [Set([ACC_TO_SPEC[id] for id in ids]) for ids in split.(best_psms_passing[!,:accession_numbers],';')];
best_psms_passing = best_psms_passing[length.(best_psms_passing[!,:species_names]).==1,:];
best_psms_passing[!,:species] = first.(best_psms_passing[!,:species_names]);
best_psms_passing = best_psms_passing[best_psms_passing[!,:Mox].==0,:]

to_keep = occursin.(first(conditions), best_psms_passing[!,:condition]) .| occursin.(last(conditions), best_psms_passing[!,:condition])

best_psms_passing = best_psms_passing[to_keep,:]
for quant_name in [:peak_area, :peak_area_normalized]
    #Get grouped psms 
    gpsms = groupPSMS(best_psms_passing, quant_name)
    #Filter out psms not common to all three replicates 
    filter!(x->x.non_missing>2,gpsms)
    #Get three proteome results
    three_proteome_results = combine(prec_group -> getLog2Diff(prec_group, 
    (string(first(conditions)), string(last(conditions)))#("E5H50Y45","E20H50Y30")
    ), 
    groupby(gpsms, [:species, :accession_numbers, :modified_sequence, :charge, :isotopes_captured]));

    jldsave(joinpath(benchmark_results_folder, string(quant_name)*"_pioneer_combined.jld2"); three_proteome_results)

    three_proteome_results = filterResults(three_proteome_results, 3, 20.0)
    prec_counts = size(three_proteome_results, 1)
    gdf = groupby(three_proteome_results, [:species])
    nt = NamedTuple.(keys(gdf))
    p = plot(legenZd=:outertopright, show = true, title = "Olsen Astral Benchmark Pioneer \n 3of3, CV<20%, 1% FDR \n $prec_counts Precursors",
    topmargin=5mm, dpi = 300)
    i = 1
    xlim = (quantile(gdf[(species = "HUMAN", )][:,:log2_mean], 0.001), quantile(gdf[(species = "HUMAN", )][:,:log2_mean], 0.999))
    ylim = (-3, 3)
    for k in ["HUMAN","YEAST","ECOLI"]
        plot!(p, gdf[(species = k,)][:,:log2_mean], gdf[(species = k,)][:,:log2_diff], color = i, label=nothing, xlim = xlim,
            ylim = ylim,
            alpha = 0.025, seriestype=:scatter)
        hline!([0.0], color = i, legend = false,
            #label = label="$(nt[i][:species])"
        )
    
        i += 1
    end

    savefig(p, joinpath(benchmark_results_folder, "PIONEER_SCATTER_"*string(quant_name)*".pdf"))

    p = plot(legenZd=:outertopright, show = true, title = "Olsen Astral Benchmark Pioneer \n 3of3, CV<20%, 1% FDR \n $prec_counts Precursors",
    topmargin=5mm, dpi = 300)
    i = 1
    xlim = (quantile(gdf[(species = "HUMAN", )][:,:log2_mean], 0.001), quantile(gdf[(species = "HUMAN", )][:,:log2_mean], 0.999))
    ylim = (-3, 3)
    SPECIES_TO_LOG2FC = Dict("HUMAN" => 0.0,
                         "YEAST" => -1.0,
                         "ECOLI" => log2(3.0))
    for k in ["HUMAN","YEAST","ECOLI"]
        density!(p, gdf[(species = k,)][:,:log2_diff], color = i, label=nothing, bins = LinRange(-3, 3, 100), xlim = (-3, 3.5), show = true, normalize = :probability)
        vline!([0.0], color = i, legend = false,
            label = label="$(nt[i][:species])"
        )
    
        i += 1
    end

    savefig(p, joinpath(benchmark_results_folder, "PIONEER_HIST_"*string(quant_name)*".pdf"))

    p = plot(legenZd=:outertopright, show = true, title = "Olsen Astral Benchmark Pioneer \n 3of3, CV<20%, 1% FDR \n $prec_counts Precursors",
    topmargin=5mm, dpi = 300)
    i = 1
    xlim = (quantile(gdf[(species = "HUMAN", )][:,:log2_mean], 0.001), quantile(gdf[(species = "HUMAN", )][:,:log2_mean], 0.999))
    ylim = (-3, 3)
    SPECIES_TO_LOG2FC = Dict("HUMAN" => 0.0,
                         "YEAST" => -1.0,
                         "ECOLI" => log2(3.0))
    for k in ["HUMAN","YEAST","ECOLI"]
        violin!(p, gdf[(species = k,)][:,:log2_diff], color = i, label=nothing, ylim = (-3, 3.5), show = true)
    
        i += 1
    end
    i = 1
    for k in ["HUMAN","YEAST","ECOLI"]
      hline!([SPECIES_TO_LOG2FC[nt[i][:species]]], color = i, label = label="$(nt[i][:species])", xlabel="Log2(A/B)")
        i += 1
    end
    savefig(p, joinpath(benchmark_results_folder, "PIONEER_VIOLIN_"*string(quant_name)*".pdf"))


    open(joinpath(benchmark_results_folder, "results.txt"), "a") do io
        original_stdout = stdout
        redirect_stdout(io) do 
            print("Summary Stats for ", string(quant_name), " \n ")
            println("CV :")
            describe(gpsms[!,:CV])
            println("log2_diff ")
            describe(gdf[1][:,:log2_diff])
            [println("Precursors w/ CV under $x ", sum(gpsms[!,:CV].<=x)) for x in [5.0, 10.0, 20.0, 30.0]]
        end
        redirect_stdout(original_stdout)
        #[println("Precursors w/ CV under $x ", sum(gpsms[!,:CV].<=x)) for x in [5.0, 10.0, 20.0, 30.0]]
    end
end



<configuration>
    <startup> 
        <supportedRuntime version="v4.0" sku=".NETFramework,Version=v4.7"/>
        <GenerateRuntimeConfigurationFiles>true</GenerateRuntimeConfigurationFiles>
    </startup>
</configuration>

