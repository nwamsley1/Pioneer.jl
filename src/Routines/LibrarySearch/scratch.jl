fixed_mods = [(p=r"C", r="C[Carb]")]
mods_dict = Dict("Carb" => Float64(57.021464),
                 "Ox" => Float64(15.994915)
                 )
f_simp, f_det, f_precs = parsePrositLib("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/Prosit_HumanYeastEcoli_NCE33_fixed_091623.csv", fixed_mods, mods_dict);

TEST_CSV = DataFrame(CSV.File("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/Prosit_HumanYeastEcoli_NCE33_corrected_091623.csv"))
value_counts(df, col) = combine(groupby(df, col), nrow)



grouped_df = groupby(TEST_CSV, :modified_sequence);
TEST_CSV[:,:accession_numbers_comb] = (combine(grouped_df) do sub_df
    accession_numbers = join(sort(unique(sub_df[:,:accession_numbers])),';')
    repeat([accession_numbers], size(sub_df)[1])
end)[:,:x1]

TEST_CSV = combine(sdf -> sdf[1,:], groupby(TEST_CSV, [:modified_sequence,:Charge]));
TEST_CSV[:,:accession_numbers] = TEST_CSV[:,:accession_numbers_comb]
select!(TEST_CSV, Not(:accession_numbers_comb))
TEST_CSV[TEST_CSV[:,:modified_sequence].=="YPKEGTHIK",:]
maximum(value_counts(TEST_CSV,:modified_sequence)[:,:nrow])

CSV.write("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/Prosit_HumanYeastEcoli_NCE33_corrected_092823.csv", TEST_CSV)



fixed_mods = [(p=r"C", r="C[Carb]")]
mods_dict = Dict("Carb" => Float64(57.021464),
                 "Ox" => Float64(15.994915)
                 )
f_simp, f_det, precursors = parsePrositLib("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/Prosit_HumanYeastEcoli_NCE33_corrected_092823.csv", fixed_mods, mods_dict, start_ion = 1, 
                                            low_frag_mz = 80.0, high_frag_mz = 2000.0,
                                            max_rank_index = 3,
                                            index_start_ion = 3);
println("SIZE f_smp ", length(f_simp))
@time begin
    @save "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/PrositRAW_HumanYeastEcoli_NCE33_corrected_100723_nOf3_indexStart3_2ratios_allStart1.jld2" f_simp f_det precursors
end

@time begin
    prosit_lib  = load("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/PrositRAW_HumanYeastEcoli_NCE33_corrected_100723_nOf3_indexStart3_2ratios_allStart1.jld2")
end

f_index = buildFragmentIndex!(prosit_lib["f_simp"], Float32(5.0), Float32(20.2))
f_det = prosit_lib["f_det"]
precursors = prosit_lib["precursors"]
@time begin
    @save  "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/PrositINDEX_HumanYeastEcoli_NCE33_corrected_100723_nOf3_indexStart3_2ratios_allStart1.jld2" f_index f_det precursors
end







println("Loaded spectral libraries in ", spec_load_time.time, " seconds")
@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_testing_082423/precursors_mouse_detailed_33NCEcorrected_start1.jld2" precursors_mouse_detailed_33NCEcorrected_start1
@load "/Users/n.t.wamsley/Projects/PROSIT/mouse_testing_082423/precursors_mouse_detailed_33NCEcorrected_start1.jld2" precursors_mouse_detailed_33NCEcorrected_start1

@time @load "/Users/n.t.wamsley/Projects/PROSIT/mouse_testing_082423/frags_mouse_detailed_33NCEcorrected_start1.jld2" frags_mouse_detailed_33NCEcorrected_start1
@time @load "/Users/n.t.wamsley/Projects/PROSIT/mouse_testing_082423/precursors_mouse_detailed_33NCEcorrected_start1.jld2" precursors_mouse_detailed_33NCEcorrected_start1.jld2
@time @load "/Users/n.t.wamsley/Projects/PROSIT/mouse_testing_082423/prosit_mouse_33NCEcorrected_start1_5ppm_15irt.jld2" prosit_mouse_33NCEcorrected_start1_5ppm_15irt

frags_mouse_detailed_33NCEcorrected_chronologer...

test_psms = best_psms[ismissing.(best_psms[:,:ρ]).==false,:]
sort(test_psms[:,[:sequence,:peak_area,:peak_area_ms1,:ρ,:precursor_idx]], :ρ)

integratePrecursor(ms2_chroms,UInt32(best_psms[N,:precursor_idx]), (0.1f0, 0.15f0, 0.15f0, Float32(66.2004), Float32(1e4)), isplot = true)
ms2_chroms[(precursor_idx=UInt32(best_psms[N,:precursor_idx]),)][:,:]

best_psms_passing = DataFrame(CSV.File("/Users/n.t.wamsley/TEST_DATA/best_psms_100123_mod.csv"));



passing_idx = Set(best_psms_passing[(best_psms_passing[:,:q_value].<=0.01).&(best_psms_passing[:,:decoy].==false),:precursor_idx])
best_psms_pass = best_psms[[id ∈ passing_idx for id in best_psms[:,:precursor_idx]],:]

include("src/Routines/LibrarySearch/integrateChroms.jl")

integratePrecursor(ms2_chroms, UInt32(best_psms_pass[N,:precursor_idx]), (0.1f0, 0.15f0, 0.15f0, Float32(best_psms_pass[N,:RT]), Float32(best_psms_pass[N,:weight])), isplot = true)
huber_loss = ms2_chroms[(precursor_idx=UInt32(best_psms_pass[N,:precursor_idx]),)][:,:]
N += 1

integratePrecursor(ms2_chroms_square, UInt32(best_psms_pass[N,:precursor_idx]), (0.1f0, 0.15f0, 0.15f0, Float32(best_psms_pass[N,:RT]), Float32(best_psms_pass[N,:weight])), isplot = true)
ms2_chroms_square[(precursor_idx=UInt32(best_psms_pass[N,:precursor_idx]),)][:,:]
#best_psms[best_psms[:,:precursor_idx].==443817,:]


CSV.write("/Users/n.t.wamsley/TEST_DATA/best_psms_pass_test.csv", best_psms_pass)
@save "/Users/n.t.wamsley/TEST_DATA/ms2_chroms.jld2" ms2_chroms


best_psms_pass = DataFrame(CSV.File("/Users/n.t.wamsley/TEST_DATA/best_psms_pass_test.csv"))
@load "/Users/n.t.wamsley/TEST_DATA/ms2_chroms.jld2" ms2_chroms

wα(α::T, x::T) where {T<:AbstractFloat} = exp(-α*x^2) + exp(-α*(x + 2)^2) + exp(-α*(x - 2)^2) - 2*exp(-α) - exp(-9*α)

function MSF(α::T, m::Int64, n::Int64) where {T<:AbstractFloat}
    out = Vector{T}(undef, 2*m + 1)
    j = 1
    for i in range(-m, m)
        x = i/(m + 1)
        if i ==0
            out[j] = 1
            j += 1
            continue
        end
        out[j] = wα(α,x)*sin( (((n + 4)/2)*π*x) )/(((n + 4)/2)*π*x)
        j += 1
    end
    return out./sum(out)
end


CSV.write("/Users/n.t.wamsley/TEST_DATA/TEST_PSMS.csv", PSMs)
model_fit = glm(@formula(target ~ entropy_sim + poisson + hyperscore +
scribe_score + topn + spectral_contrast + 
n_obs + RT_error + missed_cleavage + Mox + intensity_explained + err_log2 + total_ions), PSMs, 
Binomial(), 
ProbitLink())


PSMS_SUB = PSMs[PSMs[:,:n_obs].>2,:]

model_fit = glm(@formula(target ~ entropy_sim + poisson + hyperscore +
scribe_score + weight_log2 + topn + spectral_contrast + 
n_obs + RT_error + missed_cleavage + Mox + intensity_explained + err_log2 + total_ions), PSMS_SUB, 
Binomial(), 
ProbitLink())
Y′ = GLM.predict(model_fit, PSMS_SUB);
getQvalues!(PSMS_SUB, allowmissing(Y′),  allowmissing(PSMS_SUB[:,:decoy]));
println("Target PSMs at 25% FDR: ", sum((PSMS_SUB[:,:q_value].<=0.25).&(PSMS_SUB[:,:decoy].==false)))
PSMS_SUB[:,:prob] = allowmissing(Y′)

PSMS_SUB = PSMS_SUB[(PSMS_SUB[:,:q_value].<=0.25),:]
sort!(PSMS_SUB,:RT);
test_chroms = groupby(PSMS_SUB[:,[:precursor_idx,:q_value,:prob,:decoy,:scan_idx,:topn,:total_ions,:scribe_score,:spectral_contrast,:entropy_sim,:matched_ratio,:hyperscore,:weight,:RT,:RT_pred]],:precursor_idx);


N = 20000
test_chroms[N]

plot(test_chroms[(precursor_idx =  4190469,)][:,:RT], test_chroms[(precursor_idx =  4190469,)][:,:weight], seriestype=:scatter)
prec_id = test_chroms[N][:,:precursor_idx][1]
huber_loss = ms2_chroms[(precursor_idx=prec_id,)][:,:]
plot!(huber_loss[:,:rt], huber_loss[:,:weight], seriestype=:scatter)

N = 20000
include("src/Routines/LibrarySearch/scratch_newchroms.jl")


N = 5000


best_psms_passing= best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false),:]
prec_id = best_psms_passing[N,:precursor_idx]
integratePrecursorMS2(test_chroms,UInt32(prec_id), (0.1f0, 0.15f0, 0.15f0, Float32(66.2004), Float32(1e4)), isplot = true)
#test_chroms[(precursor_idx = prec_id,)]
N += 1

PSMs[PSMs[:,:precursor_idx].==prec_id,[:precursor_idx,:decoy,:scan_idx,:total_ions,:RT,:matched_ratio,:spectral_contrast,:weight]]
#plot(test[:,:RT], test[:,:weight], seriestype=:scatter)
best_psms[:,[:precursor_idx,:sequence,:total_ions,:entropy_sim,:matched_ratio,:spectral_contrast,:n_obs]]
test_chroms[(precursor_idx=2348785  ,)]


best_psms = combine(sdf -> getBestPSM(sdf), groupby(PSMs, [:precursor_idx]))
@profview 
@time transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursorMS2(test_chroms, 
                                                UInt32(psm[:precursor_idx]), 
                                                (0.1f0, #Fraction peak height α
                                                0.15f0, #Distance from vertical line containing the peak maximum to the leading edge at α fraction of peak height
                                                0.15f0, #Distance from vertical line containing the peak maximum to the trailing edge at α fraction of peak height
                                                Float32(psm[:RT]), psm[:weight]), isplot = false)) => [:peak_area,:GOF,:FWHM,:FWHM_01,:asymmetry,:points_above_FWHM,:points_above_FWHM_01,:σ,:tᵣ,:τ,:H,:weight_sum,:cosine_product,:entropy_sum,:scribe_sum,:ions_sum,:data_points,:ratio_sum,:base_width]);

features = [:hyperscore,:total_ions,:intensity_explained,
            :poisson,:spectral_contrast,:entropy_sim,
            :RT_error,:scribe_score,:RT,:charge,
            :city_block,:matched_ratio,:weight,:missed_cleavage,:Mox,:best_rank,:topn,:err_norm_log2,:error]
append!(features, [:peak_area,:GOF,:points_above_FWHM_01,:points_above_FWHM,:H,:ρ,:FWHM,:FWHM_01, :y_ladder,:b_ladder,:peak_area_ms1,:prec_mz,:sequence_length,:spectrum_peaks,
:log_sum_of_weights,:mean_log_spectral_contrast,:mean_log_entropy,:mean_log_probability,:mean_scribe_score,:ions_sum,:mean_matched_ratio,:data_points,:base_width_min,:ms1_ms2_diff]);
for column in names(best_psms)
    if eltype(best_psms[!,column]) == String
        continue
    end
    if  any(isnan.(skipmissing(best_psms[!,column])))
        println("$column has NaNs? ", any(isnan.(skipmissing(best_psms[!,column]))))
    end
    if any(isinf.(skipmissing(best_psms[!,column])))
        println("$column has Inf? ", any(isinf.(skipmissing(best_psms[!,column]))))
    end
end
#best_psms[:,:prob_over_data] = best_psms[:,:scribe_sum]./best_psms[:,:data_points]
#best_psms[:,:hyperscore_sum] = best_psms[:,:hyperscore_sum]./best_psms[:,:data_points]
#best_psms[:,:entropy_sum] = best_psms[:,:entropy_sum].*best_psms[:,:data_points]
#best_psms[:,:mean_ratio] = best_psms[:,:ratio_sum]./best_psms[:,:data_points]
#Train Model 
#best_psms_old = best_psms[:,:]
#best_psms = best_psms[(ismissing.(best_psms[:,:peak_area]).==false),:]
best_psms[:,:q_value] .= 0.0
for i in range(1, size(best_psms)[1])
    if ismissing(best_psms[i,:ρ])
        continue
    end
    if isnan(best_psms[i,:ρ])
best_psms[i,:ρ] = missing
    end
end

for i in range(1, size(best_psms)[1])
    if isinf(best_psms[i,:mean_log_entropy])
        best_psms[i,:mean_log_entropy] = Float16(-3)
    end
end

#best_psms[isnan.(best_psms[:,:entropy_sum]).&(ismissing.(best_psms[:,:entropy_sum]).==false),:entropy_sum] .= 0.0;
xgboost_time = @timed bst = rankPSMs!(best_psms, 
                        features,
                        colsample_bytree = 1.0, 
                        min_child_weight = 5, 
                        gamma = 1, 
                        subsample = 0.5, 
                        n_folds = 2, 
                        num_round = 200, 
                        max_depth = 10, 
                        eta = 0.05, 
                        #eta = 0.0175,
                        train_fraction = 9.0/9.0,
                        n_iters = 2);

getQvalues!(best_psms, allowmissing(best_psms[:,:prob]), allowmissing(best_psms[:,:decoy]));
best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false).&(ismissing.(best_psms[:,:peak_area]).==false),:]

best_psms[best_psms[:,:precursor_idx].==1049334,[:precursor_idx,:prob,:q_value,:weight,:sequence,:n_obs]]

best_psms[best_psms[:,:precursor_idx].==prec_id,[:precursor_idx,:prob,:q_value,:weight,:sequence,:n_obs]]


main_search_params[:min_spectral_contrast] = 0.5f0
main_search_params[:min_matched_ratio] = 0.5f0
main_search_params[:min_frag_count] = 2
main_search_params[:topN] = 1000000
sub_search_time = @timed 
PSMs = mainLibrarySearch(
    MS_TABLE,
    prosit_lib["f_index"],
    prosit_lib["f_det"],
    RT_to_iRT_map_dict[ms_file_idx], #RT to iRT map'
    UInt32(ms_file_idx), #MS_FILE_IDX
    frag_err_dist_dict[ms_file_idx],
    main_search_params,
    scan_range = (201389, 204389),
    
    #scan_range = (1, length(MS_TABLE[:masses]))
);

PSMS_TEST01 = X, Hs, IDtoROW, last_matched_col, ionMatches, ionMisses, filtered_nmatches, filtered_nmisses, nmatches, nmisses = mainLibrarySearch(
    MS_TABLE,
    prosit_lib["f_index"],
    prosit_lib["f_det"],
    RT_to_iRT_map_dict[ms_file_idx], #RT to iRT map'
    UInt32(ms_file_idx), #MS_FILE_IDX
    frag_err_dist_dict[ms_file_idx],
    main_search_params,
    scan_range = (201389, 204389),
    
    #scan_range = (0, length(MS_TABLE[:masses]))
);

println("before filter ", size(PSMS_TEST01))
filter!(x->x.best_rank==1, PSMS_TEST01);
println("after filter", size(PSMS_TEST01))
filter!(x->x.matched_ratio>0.45, PSMS_TEST01);
println("after filter", size(PSMS_TEST01))
println("Search time for $ms_file_idx ", sub_search_time.time)


main_search_params[:min_spectral_contrast] = 0.5f0
main_search_params[:min_matched_ratio] = 0.45f0
main_search_params[:min_frag_count] = 2
sub_search_time = @timed PSMs = mainLibrarySearch(
    MS_TABLE,
    prosit_lib["f_index"],
    prosit_lib["f_det"],
    RT_to_iRT_map_dict[ms_file_idx], #RT to iRT map'
    UInt32(ms_file_idx), #MS_FILE_IDX
    frag_err_dist_dict[ms_file_idx],
    main_search_params,
    scan_range = (201389, 221389),
    #scan_range = (208389, 221389),
    #scan_range = (0, length(MS_TABLE[:masses]))
);
println("before filter ", size(PSMs))
filter!(x->x.topn>1, PSMs);
println("after filter", size(PSMs))
filter!(x->x.matched_ratio>0.5, PSMs);
println("after filter", size(PSMs))
println("Search time for $ms_file_idx ", sub_search_time.time)


@profview PSMs = mainLibrarySearch(
    MS_TABLE,
    prosit_lib["f_index"],
    prosit_lib["f_det"],
    RT_to_iRT_map_dict[ms_file_idx], #RT to iRT map'
    UInt32(ms_file_idx), #MS_FILE_IDX
    frag_err_dist_dict[ms_file_idx],
    main_search_params,
    #scan_range = (201389, 204389),
    
    scan_range = (1, length(MS_TABLE[:masses]))
);

test_diff = setdiff(Set(PSMS_TEST01[:,:precursor_idx]), Set(PSMs[:,:precursor_idx]))

setdiff(Set(PSMs[:,:precursor_idx]), Set(PSMS_TEST01[:,:precursor_idx]))

PSMS_TEST01[[x ∈ test_diff for x in PSMS_TEST01[:,:precursor_idx]],:]



sum(Hs[1:nmatches,:], dims = 1)

cols_to_keep = sum(Hs[1:nmatches,:], dims = 1)./ sum(Hs[(nmatches+1):end,:], dims = 1)

IDtoROW2 =  UnorderedDictionary{UInt32, Tuple{UInt32, UInt8}}()
for (id, row) in pairs(IDtoROW)
    if cols_to_keep[first(row)] > 0.5
        insert!(IDtoROW2, id, IDtoROW[id])
    end
end
Hs = Hs[:, [true for x in cols_to_keep if x > 0.5]]
X2 = X[sum(Hs, dims = 2)[:,1].>0.0]
Hs = Hs[sum(Hs, dims = 2)[:,1].>0.0,:]

we


Hs_new = Hs[:,(scores[:matched_ratio].>=0.5).&(scores[:spectral_contrast].>=0.5)]

w = zeros(eltype(Hs_new), Hs_new.n)
#for (id, row) in pairs(IDtoROW_weights)
#    weights[row] = precursor_weights[id]# = precursor_weights[id]
#end


i = 1
for (id, row) in pairs(IDtoROW)
    println(sum(Hs_new[:,first(row)] .- Hs[:,first(IDtoROW_weights[id])]))
    i += 1
end
rtPSMs, all_matches  = firstSearch(
    MS_TABLE,
    prosit_lib["f_index"],
    prosit_lib["f_det"],
    x->x, #RT to iRT map'
    UInt32(ms_file_idx), #MS_FILE_IDX
    Laplace(zero(Float64), first_search_params[:fragment_tolerance]),
    first_search_params,
    scan_range = (200000, 210000)
    #scan_range = (1, length(MS_TABLE[:masses]))
    );
first(rtPSMs, 5)
include("src/Routines/LibrarySearch/spectralDistanceMetrics.jl")
rtPSMs, all_matches  = firstSearch(
    MS_TABLE,
    prosit_lib["f_index"],
    prosit_lib["f_det"],
    x->x, #RT to iRT map'
    UInt32(ms_file_idx), #MS_FILE_IDX
    Laplace(zero(Float64), first_search_params[:fragment_tolerance]),
    first_search_params,
    scan_range = (200000, 210000)
    #scan_range = (1, length(MS_TABLE[:masses]))
    );
first(rtPSMs, 5)

X = LinRange(-1.0, 1.0, 100)

Plots.plot(X,  
GAUSS(collect(X), Tuple((3e5, 0.0, 0.1))), 
fillrange = [0.0 for x in 1:length(X)], 
alpha = 0.25, color = :blue, show = true
); 

gx, gw = gausslegendre(100)
test_time = @timed begin for N in range(5000, 40000)
    prec_id = best_psms_passing[N,:precursor_idx]
    integratePrecursorMS2(test_chroms,
    gx, gw,
                            UInt32(prec_id), 
                            (0.1f0, 0.15f0, 0.15f0, Float32(66.2004), Float32(1e4)), 
                            isplot = false, 
                            integration_points = 60)
    end
end
N += 1

test_chroms_keys =  Set([key.precursor_idx for key in keys(test_chroms)])
test_time = @timed transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursorMS2(test_chroms, 
                                                test_chroms_keys,
                                                gw,
                                                gx,
                                                UInt32(psm[:precursor_idx]), 
                                                isplot = false)) => [:peak_area,
                                                                                :GOF,
                                                                                :FWHM,
                                                                                :FWHM_01,
                                                                                :asymmetry,
                                                                                :points_above_FWHM,
                                                                                :points_above_FWHM_01,
                                                                                :σ,
                                                                                :tᵣ,
                                                                                :τ,
                                                                                :H,
                                                                                :sum_of_weights,
                                                                                :mean_spectral_contrast,
                                                                                :entropy_sum,
                                                                                :mean_log_probability,
                                                                                :ions_sum,
                                                                                :data_points,
                                                                                :mean_matched_ratio,
                                                                                :base_width_min]);

best_psms = combine(sdf -> getBestPSM(sdf), groupby(PSMS_SUB, [:precursor_idx]))

best_psms_dict = Dict{Symbol, AbstractVector}()
for col_name in names(PSMS_SUB)
    col_type = eltype(PSMS_SUB[!,col_name])
    best_psms_dict[Symbol(col_name)] = Vector{col_type}(undef, size(test_chroms)[1])
end

best_psms_test = DataFrame(best_psms_dict)
new_cols = [(:peak_area,                   Union{Float32, Missing})
    (:GOF,                      Union{Float16, Missing})
    (:FWHM,                     Union{Float16, Missing})
    (:FWHM_01,                  Union{Float16, Missing})
    (:points_above_FWHM,        Union{UInt16, Missing})
    (:points_above_FWHM_01,     Union{UInt16, Missing})
    (:σ,                        Union{Float32, Missing})
    (:tᵣ,                       Union{Float16, Missing})
    (:τ,                        Union{Float32, Missing})
    (:H,                        Union{Float32, Missing})
    (:sum_of_weights,           Union{Float32, Missing})
    (:mean_spectral_contrast,   Union{Float32, Missing})
    (:entropy_sum,              Union{Float32, Missing})
    (:mean_log_probability,     Union{Float16, Missing})
    (:ions_sum,                 Union{UInt32, Missing})
    (:data_points,              Union{UInt32, Missing})
    (:mean_matched_ratio,       Union{Float32, Missing})
    (:base_width_min,           Union{Float16, Missing})
    (:best_scan, Union{Bool, Missing})]


best_psms_dict = Dict{Symbol, AbstractVector}()
for col_name in names(PSMS_SUB)
    col_type = eltype(PSMS_SUB[!,col_name])
    best_psms_dict[Symbol(col_name)] = Vector{col_type}(undef, size(test_chroms)[1])
end

time_test = @timed integratePrecursors(test_chroms, best_psms_test, names(test_chroms[1]))

[(Symbol(col_name), eltype(PSMs[!, col_name])) for col_name in names(PSMs)]
#combine(test_chroms, AsTable(:) => ByRow 
result_df = combine(groupby(test_chroms, :precursor_idx)) do precursor_group
    #idx_max_c = argmax(group.c)
    #row_max_c = group[idx_max_c, :]
    #sum_of_other_entries = sum(row_max_c[!, [:a, :b]])
    
    integratePrecursorMS2(precursor_group, 
                            test_chroms_keys,
                            gw,
                            gx
                        )


    return DataFrame(
        a = row_max_c.a,
        b = row_max_c.b,
        c = row_max_c.c,
        sum_of_other_entries = sum_of_other_entries
    )
end

Profile.Allocs.Clear()
using PProf

gx, gw = gausslegendre(100)
ms1_chrom_keys = keys(ms1_chroms)
best_psms_sub = copy(first(best_psms, 100000))
time_test = @timed transform!(best_psms_sub, AsTable(:) => ByRow(psm -> integratePrecursor(ms1_chroms, 
                                    ms1_chrom_keys,
                                    gx,
                                    gw,
                                    UInt32(psm[:precursor_idx]), 
                                    [Float32(psm[:σ]),
                                     Float32(psm[:tᵣ]),
                                     Float32(psm[:τ]), 
                                     Float32(1e4)],
                                    isplot = false)) => [:peak_area_ms1,
                                    :GOF_ms1,
                                    :FWHM_ms1,
                                    :FWHM_01_ms1,
                                    :asymmetry_ms1,
                                    :points_above_FWHM_ms1,
                                    :points_above_FWHM_01_ms1,
                                    :σ_ms1,:tᵣ_ms1,:τ_ms1,:H_ms1]);
test_psm = best_psms[10000,:]
integratePrecursor(ms1_chroms, 
                                    ms1_chrom_keys,
                                    gx,
                                    gw,
                                    UInt32(test_psm[:precursor_idx]), 
                                    [Float32(test_psm[:σ]),
                                     Float32(test_psm[:tᵣ]),
                                     Float32(test_psm[:τ]), 
                                     Float32(1e4)],
                                    isplot = false)
best_psms_passing = best_psms[(best_psms[:,:q_value].<=0.01).&(best_psms[:,:decoy].==false).&(ismissing.(best_psms[:,:peak_area]).==false),:]
N = 10000
test_psm = best_psms_passing[N,:]
integratePrecursor(ms1_chroms, 
                                           ms1_chrom_keys,
                                           gx,
                                           gw,
                                           UInt32(test_psm[:precursor_idx]), 
                                           [Float32(test_psm[:σ]),
                                            Float32(test_psm[:tᵣ]),
                                            Float32(test_psm[:τ]), 
                                            Float32(1e4)],
                                           isplot = true)
                                           N += 1
#Hs_new[:,1]
 #=
            IDtoMatchedRatio = UnorderedDictionary{UInt32, Float32}()
            for i in range(1, nmatches)
                m = ionMatches[i]
                if haskey(IDtoMatchedRatio, getPrecID(m))
                    IDtoMatchedRatio[getPrecID(m)] += getPredictedIntenisty(m)
                else
                    insert!(IDtoMatchedRatio, getPrecID(m), getPredictedIntenisty(m))
                end
            end

            for i in range(1, nmisses)
                m = ionMisses[i]
                if !haskey(IDtoMatchedRatio, getPrecID(m))
                    insert!(IDtoMatchedRatio, getPrecID(m), zero(Float32))
                end
            end


            #println("nmatches before ", nmatches)
            #println("nmisses before ", nmisses)

            last_open = 1
            for i in range(1, nmatches)
                m = ionMatches[i]
                if IDtoMatchedRatio[getPrecID(m)] >= 0.5
                    if getPredictedIntenisty(ionMatches[i]) > 0.0
                        filtered_nmatches += 1
                        ionMatches[last_open] = ionMatches[i]
                        last_open = filtered_nmatches + 1
                    end
                end
            end

            last_open = 1
            for i in range(1, nmisses)
                m = ionMisses[i]
                if IDtoMatchedRatio[getPrecID(m)] >= 0.5
                    if getPredictedIntenisty(ionMisses[i]) > 0.0
                        filtered_nmisses += 1
                        ionMisses[last_open] = ionMisses[i]
                        last_open = filtered_nmisses + 1
                    end
                end
            end
            =#


             #Initial guess from non-negative least squares. May need to reconsider if this is beneficial
                #weights = sparseNMF(Hs, X, λ, γ, regularize, max_iter=max_iter, tol=nmf_tol)[:]
                #weights = sparseNMF(Hs, X, λ, γ, regularize, max_iter=max_iter, tol=Hs.n)[:]
                #println("i $i")
                #if i == 201394  
                #    return X, Hs, IDtoROW, last_matched_col, ionMatches, ionMisses, filtered_nmatches, filtered_nmisses, nmatches, nmisses
                #end
                #=
                build_matrix_time += @elapsed begin 
                cols_to_keep = sum(Hs[1:nmatches,:], dims = 1)./sum(Hs[(nmatches+1):end,:], dims = 1)
                
                IDtoROW2 =  UnorderedDictionary{UInt32, Tuple{UInt32, UInt8}}()
                col = UInt32(1)
                for (id, row) in pairs(IDtoROW)
                    if cols_to_keep[first(row)] > 0.6
                        insert!(IDtoROW2, id, (col , last(IDtoROW[id])))
                        col += one(UInt32)
                    end
                end
                IDtoROW = IDtoROW2
                select_ions_time += @elapsed begin
                    Hs = Hs[:, [true for x in cols_to_keep if x > 0.6]]
                    #X = X[sum(Hs, dims = 2)[:,1].>0.0]
                    #select_ions_time += @elapsed Hs = Hs[sum(Hs, dims = 2)[:,1].>0.0,:]
                    weights = zeros(eltype(Hs), Hs.n)
                end
                last_matched_col = findfirst(x->iszero(x), X)
                if last_matched_col === nothing
                    last_matched_col = length(X)
                else
                    last_matched_col = last_matched_col - 1
                end

                filtered_nmatches = last_matched_col
                filtered_nmisses = length(X) - last_matched_col
                for (id, row) in pairs(IDtoROW)
                    weights[first(row)] = precursor_weights[id]
                end
                
                end
                =(#