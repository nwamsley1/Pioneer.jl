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
f_simp, f_det, precursors = parsePrositLib("/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/Prosit_HumanYeastEcoli_NCE33_corrected_092823.csv", fixed_mods, mods_dict, start_ion = 1, low_frag_mz = 80.0, high_frag_mz = 2000.0);
@time begin
    @save "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/PrositRAW_HumanYeastEcoli_NCE33_corrected_100423_01.jld2" f_simp f_det precursors
end

@time begin
    prosit_lib  = load( "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/PrositRAW_HumanYeastEcoli_NCE33_corrected_100423_01.jld2")
end

f_index = buildFragmentIndex!(prosit_lib["f_simp"], Float32(5.0), Float32(20.2))
f_det = prosit_lib["f_det"]
precursors = prosit_lib["precursors"]
@time begin
    @save "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/PrositINDEXED_HumanYeastEcoli_NCE33_corrected_100423_5ppm_20irt.jld2" f_index f_det precursors
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
plot(test_chroms[N][:,:RT], test_chroms[N][:,:weight], seriestype=:scatter)
prec_id = test_chroms[N][:,:precursor_idx][1]
huber_loss = ms2_chroms[(precursor_idx=prec_id,)][:,:]
plot!(huber_loss[:,:rt], huber_loss[:,:weight], seriestype=:scatter)

N = 20000
include("src/Routines/LibrarySearch/scratch_newchroms.jl")


N = 20000
prec_id = test_chroms[N][:,:precursor_idx][1]
test_chroms[N]
integratePrecursor(test_chroms,UInt32(prec_id), (0.1f0, 0.15f0, 0.15f0, Float32(66.2004), Float32(1e4)), isplot = true)
N += 1
#plot(test[:,:RT], test[:,:weight], seriestype=:scatter)

best_psms = combine(sdf -> getBestPSM(sdf), groupby(PSMs, [:precursor_idx]))
sub_search_time = @timed transform!(best_psms, AsTable(:) => ByRow(psm -> integratePrecursor(test_chroms, 
                                                UInt32(psm[:precursor_idx]), 
                                                (0.1f0, #Fraction peak height α
                                                0.15f0, #Distance from vertical line containing the peak maximum to the leading edge at α fraction of peak height
                                                0.15f0, #Distance from vertical line containing the peak maximum to the trailing edge at α fraction of peak height
                                                Float32(psm[:RT]), psm[:weight]), isplot = false)) => [:peak_area,:GOF,:FWHM,:FWHM_01,:asymmetry,:points_above_FWHM,:points_above_FWHM_01,:σ,:tᵣ,:τ,:H,:weight_sum,:hyperscore_sum,:entropy_sum,:scribe_sum,:ions_sum,:data_points,:ratio_sum,:base_width]);

                                                features = [:hyperscore,:total_ions,:intensity_explained,
:poisson,:spectral_contrast,:entropy_sim,
:RT_error,:scribe_score,:RT,:charge,
:city_block,:matched_ratio,:weight,:missed_cleavage,:Mox,:best_rank,:topn,:err_norm,:error]
append!(features, [:peak_area,:GOF,:asymmetry,:points_above_FWHM_01,:H,:weight_sum,:hyperscore_sum,:entropy_sum,:scribe_sum,:ions_sum,:data_points,:ratio_sum,:base_width]);
#Train Model 
xgboost_time = @timed bst = rankPSMs!(best_psms, 
            features,
            colsample_bytree = 1.0, 
            min_child_weight = 10, 
            gamma = 10, 
            subsample = 0.5, 
            n_folds = 2, 
            num_round = 200, 
            max_depth = 10, 
            eta = 0.0375, 
            max_train_size = size(best_psms)[1]);