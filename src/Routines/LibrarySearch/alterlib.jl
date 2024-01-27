using CSV, JLD2
include("src\\Routines\\LibrarySearch\\test.jl")

params = JSON.parse(read("./data/example_config/LibrarySearch.json", String));
SPEC_LIB_DIR = pwd()*"\\..\\data\\nOf3_y4b3_102123_sulfur"#"/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123_sulfur"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf3_y4b3_102123"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf5_ally3b2/"
#SPEC_LIB_DIR = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/nOf1_indy6b5_ally3b2/"
MS_DATA_DIR = pwd()*"\\..\\data\\RAW"#"/Users/n.t.wamsley/TEST_DATA/PXD028735/"
#MS_DATA_DIR = "/Users/n.t.wamsley/TEST_DATA/mzXML/"
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, file) for file in filter(file -> isfile(joinpath(MS_DATA_DIR, file)) && match(r"\.arrow$", file) != nothing, readdir(MS_DATA_DIR))];
EXPERIMENT_NAME = "TEST_y4b3_nOf5"

f_det_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"f_det\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
f_index_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"f_index\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
precursors_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"precursors\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
#prosit_lib_path = "/Users/n.t.wamsley/RIS_temp/BUILD_PROSIT_LIBS/PrositINDEX_HumanYeastEcoli_NCE33_corrected_100723_nOf3_indexStart3_2ratios_allStart2.jld2"
println("Loading spectral libraries into main memory...")
prosit_lib = Dict{String, Any}()
spec_load_time = @timed begin
    const f_index = load(f_index_path)["f_index"];
    prosit_lib["f_index"] = f_index#["f_index"]
    const f_det = load(f_det_path)["f_det"]
    prosit_lib["f_det"] = f_det#["f_det"];
    const precursors = load(precursors_path)["precursors"]
    prosit_lib["precursors"] = precursors#["precursors"];
end


fsimp_path = [joinpath(SPEC_LIB_DIR, file) for file in filter(file -> isfile(joinpath(SPEC_LIB_DIR, file)) && match(r"f_simp\.jld2$", file) != nothing, readdir(SPEC_LIB_DIR))][1];
const f_simp = load(fsimp_path)["f_simp"]

simple_frags = Vector{SimpleFrag{Float32}}()
for i in ProgressBar(range(1, length(f_simp)))
    score = zero(UInt8)
    if f_simp[i].prec_intensity > 0.6f0
        score = UInt8(2)
    else
        score = one(UInt8)
    end
    push!(
        simple_frags,
        SimpleFrag(
                    f_simp[i].frag_mz,
                    f_simp[i].prec_id,
                    f_simp[i].prec_mz,
                    f_simp[i].prec_rt,
                    f_simp[i].prec_charge,
                    score
        )
    )
end

@save "C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\nOf3_y4b3_102123_sulfur\\HumanYeastEcoli_NCE33COR_101723_nOf3_indy4b3_ally3b2_simple_frags.jld2" simple_frags

detailed_frags = Vector{Vector{DetailedFrag{Float32}}}()

for prec_frags in ProgressBar(f_det)
    push!(detailed_frags,
          Vector{DetailedFrag{Float32}}())
    for frag in prec_frags
        push!(detailed_frags[end],
                DetailedFrag(
                    frag.prec_id,

                    frag.frag_mz,
                    Float16(frag.intensity),

                    frag.is_y_ion,
                    frag.is_isotope,

                    frag.frag_charge,
                    frag.ion_position,
                    frag.prec_charge,
                    frag.rank,
                    frag.sulfur_count
                ))
    end
end

@save "C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\nOf3_y4b3_102123_sulfur\\HumanYeastEcoli_NCE33COR_101723_nOf3_indy4b3_ally3b2_detailed_frags.jld2" detailed_frags


library_precursors = Vector{LibraryPrecursorIon{Float32}}()
for prec in ProgressBar(precursors)
    push!(
        library_precursors,
        LibraryPrecursorIon(
                    prec.iRT,
                    prec.mz,

                    prec.isDecoy,

                    prec.accession_numbers,
                    prec.sequence,

                    prec.charge,
                    prec.missed_cleavages,
                    prec.variable_mods,
                    prec.length,
                    prec.sulfur_count
        )
    )
end

@save "C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\nOf3_y4b3_102123_sulfur\\HumanYeastEcoli_NCE33COR_101723_nOf3_indy4b3_ally3b2_library_precursors.jld2" library_precursors


simple_frags = load("C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\nOf3_y4b3_102123_sulfur\\HumanYeastEcoli_NCE33COR_101723_nOf3_indy4b3_ally3b2_simple_frags.jld2")["simple_frags"]

function buildFragmentIndex!(frag_ions::Vector{SimpleFrag{T}}, bin_ppm::AbstractFloat, rt_size::AbstractFloat; 
    low_frag_mz::AbstractFloat = 80.0, high_frag_mz::AbstractFloat = 3000.0, low_prec_mz::AbstractFloat = 300.0, high_prec_mz::AbstractFloat = 1100.0) where {T<:AbstractFloat}
    println("sorting...")
    sort!(frag_ions, by = x->getMZ(x))
    println("sorted")
    #The fragment ions are divided into bins of roughtly equal m/z width.
    #That should correspond to roughly half the fragment mass accuracy of the detector?
    frag_index = FragmentIndex(T) 

    function fillPrecursorBin!(frag_index::FragmentIndex{<:AbstractFloat}, frag_ions::Vector{SimpleFrag{T}}, frag_bin_idx::UInt32, start::Int, stop::Int, low_prec_mz::AbstractFloat, high_prec_mz::AbstractFloat)
        for ion_index in range(start, stop)
            prec_id = getPrecID(frag_ions[ion_index])
                prec_mz = getPrecMZ(frag_ions[ion_index])

                if (prec_mz < low_prec_mz) | (prec_mz > high_prec_mz) #Precursor m/z outside the bounds
                    continue
                end
                #Add precursor corresponding to the charge state
                addPrecursorBinFragment!(frag_index,
                                         frag_bin_idx,
                                         PrecursorBinFragment(
                                                                prec_id, 
                                                                prec_mz, 
                                                                getScore(frag_ions[ion_index]), 
                                                                getPrecCharge(frag_ions[ion_index])
                                                            )
                                    )
        end
    end

    rt_bin_idx = UInt32(1) #Current fragment bin index
    prec_bin_idx = one(UInt32)
    start = 1 #Fragment index of the first fragment in the current bin
    stop = 2

    diff = getPPM(getMZ(frag_ions[start]), bin_ppm) #ppm tolerance of the current fragment bin
    #Build bins 
    while stop < (length(frag_ions)) + 1

        if (stop % 1_000_000) == 0
            println(stop/1_000_000)
        end
        #Haven't reached minimum fragment m/z yet
        if getMZ(frag_ions[start]) < low_frag_mz
            start += 1
            stop += 1
            continue
        end

        #Does the difference between the highest and lowest m/z in the bin 
        #enought exceed 10 ppm of the lowest m/z?
        if ((getMZ(frag_ions[min(stop + 1, length(frag_ions))]) - getMZ(frag_ions[start])) > diff) | ((stop + 1) >= length(frag_ions))
            
            #Add a new fragment bin
            addFragmentBin!(frag_index, 
                            FragBin(getMZ(frag_ions[start]),
                                    getMZ(frag_ions[stop]),
                                    rt_bin_idx
                                    )
                            )

            #Holds the RT bins for these fragments. 
            addRTBin!(frag_index)
            
            #Sort fragments in the frag_bin by retention time
            frag_ions[start:stop] = sort(frag_ions[start:stop], by = frag -> frag.prec_irt)
            start_rt_idx = start
            #Make new precursor and retention time bins. 
            for i in start:stop
                #The first and last fragment differ in retention time by a threshold
                if ((frag_ions[i].prec_irt + (-1)*frag_ions[start_rt_idx].prec_irt) > rt_size) | (i == stop)

                    #Create a retention time bin to hold the start_idx:i-1 fragments 
                    push!(frag_index.rt_bins[rt_bin_idx] , RTBin(frag_ions[start_rt_idx].prec_irt,
                                                            frag_ions[max(i-1, start_rt_idx)].prec_irt,
                                                            prec_bin_idx)
                            )

                    #Create a precursor bin to hold the start_idx:i-1fragments. 
                    addPrecursorBin!(frag_index, 
                                        #Preallocate an empty precursor bin of the correct length 
                                        PrecursorBin{T}(Vector{PrecursorBinFragment{T}}())#undef, (last_frag_in_bin - start + 1)*length(charges)))
                                        )
                    if i == stop
                        fillPrecursorBin!(frag_index, frag_ions, prec_bin_idx, start_rt_idx, i, low_prec_mz, high_prec_mz)
                    else
                        fillPrecursorBin!(frag_index, frag_ions, prec_bin_idx, start_rt_idx, i-1, low_prec_mz, high_prec_mz)
                    end

                    start_rt_idx = i
                    prec_bin_idx += one(UInt32)
                end
            end
            #Update counters and ppm tolerance 
            rt_bin_idx += UInt32(1)
            start = stop + 1
            diff = getPPM(getMZ(frag_ions[min(start, length(frag_ions))]), bin_ppm)
            stop = start
            #Maximum fragment m/z reached. Stop adding bins.
            if getMZ(frag_ions[min(start, length(frag_ions))]) > high_frag_mz
                break
            end
            
        else
            stop += 1
        end
    end

    function sortPrecursorBins!(frag_index::FragmentIndex{<:AbstractFloat})
        for i in 1:length(frag_index.precursor_bins)
            sort!(frag_index.precursor_bins[i].precs, by = x->getPrecMZ(x));
        end
        return nothing
    end
    sortPrecursorBins!(frag_index)
    return frag_index
end

library_fragment_index = buildFragmentIndex!(
    simple_frags,
    5.0f0,
    15.0f0
)

@save "C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\nOf3_y4b3_102123_sulfur\\HumanYeastEcoli_NCE33COR_101723_nOf3_indy4b3_ally3b2_5ppm_15irt_library_fragment_index.jld2" library_fragment_index


library_fragment_index = buildFragmentIndex!(
    simple_frags,
    5.0f0,
    20.0f0
)


@save "C:\\Users\\n.t.wamsley\\Pioneer.jl-1\\..\\data\\nOf3_y4b3_102123_sulfur\\HumanYeastEcoli_NCE33COR_101723_nOf3_indy4b3_ally3b2_5ppm_20irt_library_fragment_index.jld2" library_fragment_index
