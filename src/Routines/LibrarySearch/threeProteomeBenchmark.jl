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
                ).&(psms[:,:points_above_FWHM_01].>1
                ),:]
    psms[!,quant_col] = max.(psms[!,quant_col], 0.0)
    gpsms = groupby(psms, 
                            [:species,:accession_numbers,:sequence,:charge,:isotopes_captured,:condition]
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

test_time = @timed begin
    println("A")
end
#Create results folder 
benchmark_results_folder = params_[:benchmark_params]["results_folder"]
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
        println("total time $total_time")
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

transform!(best_psms_passing, AsTable(:) => ByRow(precursor -> getCondition(precursor, :file_name)) => [:condition, :biological, :technical]);
best_psms_passing[!,:species] .= "HUMAN"

conditions = unique(best_psms_passing[!,:condition])

quant_name = :trapezoid_area_normalized
cv_param = 20.0
for quant_name in [:peak_area_normalized, :trapezoid_area_normalized, :weight_normalized]
    #Get grouped psms 
    gpsms = groupPSMS(best_psms_passing, quant_name)
    #Filter out psms not common to all three replicates 
    filter!(x->x.non_missing>2,gpsms)
    #Get three proteome results
    three_proteome_results = combine(prec_group -> getLog2Diff(prec_group, 
    (string(first(conditions)), string(last(conditions)))#("E5H50Y45","E20H50Y30")
    ), 
    groupby(gpsms, [:species, :accession_numbers, :sequence, :charge, :isotopes_captured]));


    three_proteome_results = filterResults(three_proteome_results, 3, 20.0)
    prec_counts = size(three_proteome_results, 1)
    gdf = groupby(three_proteome_results, [:species])
    nt = NamedTuple.(keys(gdf))
    p = plot(legenZd=:outertopright, show = true, title = "Puyvelde et al. 2023 w/ PIONEER \n 6of9, CV<$cv_param%, 1% FDR \n $prec_counts Precursors",
    topmargin=5mm, dpi = 300)
    i = 1
    xlim = (quantile(gdf[(species = "HUMAN", )][:,:log2_mean], 0.001), quantile(gdf[(species = "HUMAN", )][:,:log2_mean], 0.999))
    ylim = (-3, 3)
    for (k,v) in pairs(gdf)
        plot!(p, gdf[k][:,:log2_mean], gdf[k][:,:log2_diff], color = i, label=nothing, xlim = xlim,
            ylim = ylim,
            alpha = 0.01, seriestype=:scatter)
        hline!([0.0], color = i, legend = false,
            #label = label="$(nt[i][:species])"
        )
    
        i += 1
    end

    savefig(p, joinpath(benchmark_results_folder, "PIONEER_SCATTER_"*string(quant_name)*".pdf"))

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